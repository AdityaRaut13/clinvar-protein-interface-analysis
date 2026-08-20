#!/usr/bin/env python

import hashlib
import json
import os
import re
import time
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm


ROOT = Path(os.environ["IF_PROJECT"])

PDB_FILE = (
    ROOT
    / "data/processed/structure_discovery"
    / "variant_observed_pdb_manifest.tsv.gz"
)
SEGMENT_FILE = (
    ROOT
    / "data/processed/structure_discovery"
    / "variant_observed_pdb_segments.parquet"
)
TARGET_FILE = (
    ROOT
    / "data/processed/uniprot_ready"
    / "unique_uniprot_structural_targets.tsv.gz"
)

CACHE_DIR = ROOT / "cache/rcsb/graphql_assemblies"
OUTPUT_DIR = ROOT / "data/processed/assembly_discovery"

LIMIT = int(os.environ.get("RCSB_LIMIT", "0"))
BATCH_SIZE = int(os.environ.get("RCSB_BATCH_SIZE", "100"))
SUFFIX = f"_test{LIMIT}" if LIMIT else ""

ENTRY_OUTPUT = OUTPUT_DIR / f"rcsb_entry_metadata{SUFFIX}.parquet"
ASSEMBLY_OUTPUT = OUTPUT_DIR / f"rcsb_assembly_metadata{SUFFIX}.parquet"
MANIFEST_OUTPUT = (
    OUTPUT_DIR / f"target_complex_assembly_manifest{SUFFIX}.tsv.gz"
)
SUMMARY_OUTPUT = (
    OUTPUT_DIR / f"rcsb_assembly_summary{SUFFIX}.json"
)

GRAPHQL_URL = "https://data.rcsb.org/graphql"

QUERY = """
query Entries($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id

    rcsb_entry_info {
      experimental_method
      resolution_combined
      assembly_count
    }

    assemblies {
      rcsb_id

      rcsb_assembly_container_identifiers {
        assembly_id
      }

      rcsb_assembly_info {
        polymer_entity_instance_count
        polymer_entity_instance_count_protein
        polymer_entity_count_protein
      }

      pdbx_struct_assembly {
        details
        method_details
        oligomeric_details
        oligomeric_count
      }

      polymer_entity_instances {
        rcsb_polymer_entity_instance_container_identifiers {
          entity_id
          asym_id
          auth_asym_id
        }

        polymer_entity {
          entity_poly {
            rcsb_entity_polymer_type
          }

          rcsb_polymer_entity_container_identifiers {
            reference_sequence_identifiers {
              database_name
              database_accession
            }
          }
        }
      }
    }
  }
}
"""

CACHE_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def join_unique(values):
    return ";".join(sorted({
        str(value)
        for value in values
        if value is not None and str(value).strip()
    }))


def as_number(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def batch_cache_path(ids):
    digest = hashlib.sha256(
        ",".join(ids).encode()
    ).hexdigest()[:20]

    return CACHE_DIR / f"entries_{digest}.json"


def fetch_batch(session, ids):
    cache_path = batch_cache_path(ids)

    if cache_path.exists():
        with open(cache_path) as handle:
            return json.load(handle), "cached"

    last_error = None

    for attempt in range(6):
        try:
            response = session.post(
                GRAPHQL_URL,
                json={
                    "query": QUERY,
                    "variables": {"ids": ids},
                },
                timeout=300,
            )
            response.raise_for_status()
            payload = response.json()

            if payload.get("errors"):
                raise RuntimeError(payload["errors"])

            temporary = cache_path.with_suffix(".json.tmp")

            with open(temporary, "w") as handle:
                json.dump(payload, handle)

            temporary.replace(cache_path)
            return payload, "downloaded"

        except Exception as error:
            last_error = str(error)
            time.sleep(2 ** attempt)

    raise RuntimeError(
        f"RCSB batch failed after retries: {last_error}"
    )


pdb_manifest = pd.read_csv(
    PDB_FILE,
    sep="\t",
    dtype={"pdb_id": str},
)

pdb_ids = (
    pdb_manifest["pdb_id"]
    .dropna()
    .astype(str)
    .str.upper()
    .drop_duplicates()
    .sort_values()
    .tolist()
)

if LIMIT:
    pdb_ids = pdb_ids[:LIMIT]

segments = pd.read_parquet(
    SEGMENT_FILE,
    columns=[
        "pdb_id",
        "auth_chain_id",
        "uniprot_canonical",
    ],
)

segments["pdb_id"] = (
    segments["pdb_id"].astype(str).str.upper()
)
segments["auth_chain_id"] = (
    segments["auth_chain_id"].astype(str)
)

chain_targets = {
    key: set(group["uniprot_canonical"].dropna().astype(str))
    for key, group in segments.groupby(
        ["pdb_id", "auth_chain_id"],
        sort=False,
    )
}

target_table = pd.read_csv(TARGET_FILE, sep="\t")

target_variant_counts = dict(
    zip(
        target_table["uniprot_canonical"],
        target_table["n_variants"],
    )
)

batches = [
    pdb_ids[index:index + BATCH_SIZE]
    for index in range(0, len(pdb_ids), BATCH_SIZE)
]

session = requests.Session()
session.headers.update({
    "User-Agent": "IITK-ClinVar-interface/1.0"
})

entry_rows = []
assembly_rows = []
downloaded_batches = 0
cached_batches = 0

for batch in tqdm(batches, desc="RCSB GraphQL batches"):
    payload, source = fetch_batch(session, batch)

    if source == "cached":
        cached_batches += 1
    else:
        downloaded_batches += 1

    entries = (
        payload.get("data", {}).get("entries")
        or []
    )

    for entry in entries:
        if not entry:
            continue

        pdb_id = str(entry["rcsb_id"]).upper()
        entry_info = entry.get("rcsb_entry_info") or {}

        resolutions = entry_info.get(
            "resolution_combined"
        ) or []

        valid_resolutions = [
            float(value)
            for value in resolutions
            if value is not None
        ]

        best_resolution = (
            min(valid_resolutions)
            if valid_resolutions
            else None
        )

        assemblies = entry.get("assemblies") or []

        entry_rows.append({
            "pdb_id": pdb_id.lower(),
            "experimental_method":
                entry_info.get("experimental_method"),
            "best_resolution": best_resolution,
            "reported_assembly_count":
                as_number(entry_info.get("assembly_count")),
            "returned_assembly_count": len(assemblies),
        })

        for assembly in assemblies:
            identifiers = (
                assembly.get(
                    "rcsb_assembly_container_identifiers"
                )
                or {}
            )
            assembly_info = (
                assembly.get("rcsb_assembly_info")
                or {}
            )
            description = (
                assembly.get("pdbx_struct_assembly")
                or {}
            )

            assembly_id = str(
                identifiers.get("assembly_id") or ""
            )

            protein_instances = as_number(
                assembly_info.get(
                    "polymer_entity_instance_count_protein"
                )
            )

            protein_entities = as_number(
                assembly_info.get(
                    "polymer_entity_count_protein"
                )
            )

            all_protein_label_chains = set()
            all_protein_auth_chains = set()
            target_label_chains = set()
            target_auth_chains = set()
            target_accessions = set()

            for instance in (
                assembly.get("polymer_entity_instances")
                or []
            ):
                container = (
                    instance.get(
                        "rcsb_polymer_entity_instance_container_identifiers"
                    )
                    or {}
                )

                polymer_entity = (
                    instance.get("polymer_entity")
                    or {}
                )
                entity_poly = (
                    polymer_entity.get("entity_poly")
                    or {}
                )

                polymer_type = str(
                    entity_poly.get(
                        "rcsb_entity_polymer_type"
                    )
                    or ""
                ).lower()

                if "protein" not in polymer_type:
                    continue

                label_chain = str(
                    container.get("asym_id") or ""
                )
                auth_chain = str(
                    container.get("auth_asym_id") or ""
                )

                if label_chain:
                    all_protein_label_chains.add(label_chain)
                if auth_chain:
                    all_protein_auth_chains.add(auth_chain)

                accessions = chain_targets.get(
                    (pdb_id, auth_chain),
                    set(),
                )

                if accessions:
                    target_accessions.update(accessions)

                    if label_chain:
                        target_label_chains.add(label_chain)
                    if auth_chain:
                        target_auth_chains.add(auth_chain)

            is_numeric_assembly = bool(
                re.fullmatch(r"\d+", assembly_id)
            )
            is_protein_complex = (
                protein_instances is not None
                and protein_instances >= 2
            )
            target_present = bool(target_accessions)

            is_candidate = (
                is_numeric_assembly
                and is_protein_complex
                and target_present
            )

            assembly_rows.append({
                "pdb_id": pdb_id.lower(),
                "assembly_id": assembly_id,
                "experimental_method":
                    entry_info.get("experimental_method"),
                "best_resolution": best_resolution,
                "protein_instance_count": protein_instances,
                "protein_entity_count": protein_entities,
                "all_protein_label_chains":
                    join_unique(all_protein_label_chains),
                "all_protein_auth_chains":
                    join_unique(all_protein_auth_chains),
                "target_label_chains":
                    join_unique(target_label_chains),
                "target_auth_chains":
                    join_unique(target_auth_chains),
                "target_accessions":
                    join_unique(target_accessions),
                "n_target_accessions":
                    len(target_accessions),
                "target_present": target_present,
                "is_protein_complex": is_protein_complex,
                "is_numeric_assembly": is_numeric_assembly,
                "is_candidate": is_candidate,
                "assembly_details":
                    description.get("details"),
                "assembly_method_details":
                    description.get("method_details"),
                "oligomeric_details":
                    description.get("oligomeric_details"),
                "oligomeric_count":
                    as_number(
                        description.get("oligomeric_count")
                    ),
            })

entries_df = pd.DataFrame(entry_rows)
assemblies_df = pd.DataFrame(assembly_rows)

entries_df.to_parquet(
    ENTRY_OUTPUT,
    index=False,
    compression="zstd",
)

assemblies_df.to_parquet(
    ASSEMBLY_OUTPUT,
    index=False,
    compression="zstd",
)

candidates = assemblies_df[
    assemblies_df["is_candidate"]
].copy()

candidates["filename"] = (
    candidates["pdb_id"]
    + "-assembly"
    + candidates["assembly_id"]
    + ".cif.gz"
)

candidates["download_url"] = (
    "https://files.rcsb.org/pub/pdb/data/assemblies/"
    "mmCIF/all/"
    + candidates["filename"]
)

candidates.to_csv(
    MANIFEST_OUTPUT,
    sep="\t",
    index=False,
    compression="gzip",
)

candidate_accessions = {
    accession
    for value in candidates["target_accessions"]
    for accession in str(value).split(";")
    if accession
}

summary = {
    "requested_pdb_entries": len(pdb_ids),
    "returned_pdb_entries": int(
        entries_df["pdb_id"].nunique()
    ),
    "downloaded_batches": downloaded_batches,
    "cached_batches": cached_batches,
    "assemblies_examined": int(len(assemblies_df)),
    "protein_complex_assemblies": int(
        assemblies_df["is_protein_complex"].sum()
    ),
    "target_complex_assemblies": int(len(candidates)),
    "unique_target_complex_pdb_entries": int(
        candidates["pdb_id"].nunique()
    ),
    "target_accessions_in_complexes": len(
        candidate_accessions
    ),
    "variants_on_targets_in_complexes": int(sum(
        target_variant_counts.get(accession, 0)
        for accession in candidate_accessions
    )),
}

with open(SUMMARY_OUTPUT, "w") as handle:
    json.dump(summary, handle, indent=2)

print()
print(json.dumps(summary, indent=2))
print()
print("Assembly manifest:", MANIFEST_OUTPUT)