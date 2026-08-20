#!/usr/bin/env python

import gzip
import json
import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from tqdm import tqdm


ROOT = Path(os.environ["IF_PROJECT"])

MANIFEST_FILE = (
    ROOT
    / "data/processed/assembly_discovery"
    / "target_complex_assembly_manifest.tsv.gz"
)
FULL_SIFTS_FILE = (
    ROOT
    / "data/input/sifts"
    / "pdb_chain_uniprot.tsv.gz"
)
TARGET_SEGMENT_FILE = (
    ROOT
    / "data/processed/structure_discovery"
    / "variant_observed_pdb_segments.parquet"
)
ASSEMBLY_DIR = ROOT / "structures/assemblies"
OUTPUT_DIR = ROOT / "data/processed/chain_pairs"

LIMIT = int(os.environ.get("PARSE_LIMIT", "0"))
WORKERS = int(os.environ.get("PARSE_WORKERS", "6"))
CHUNK_SIZE = int(os.environ.get("PARSE_CHUNK_SIZE", "250"))
PINDER_CUTOFF = 10.0
HEAVY_CUTOFF = 5.0
BACKBONE_ATOMS = {"N", "CA", "C", "O"}

SUFFIX = (
    f"_union10A5A_test{LIMIT}"
    if LIMIT
    else "_union10A5A"
)

CHAIN_PARTS = OUTPUT_DIR / f"protein_chains{SUFFIX}"
PAIR_PARTS = OUTPUT_DIR / f"contacting_chain_pairs{SUFFIX}"
STATUS_PARTS = OUTPUT_DIR / f"assembly_parse_status{SUFFIX}"
SUMMARY_OUTPUT = OUTPUT_DIR / f"chain_pair_summary{SUFFIX}.json"

CHAIN_UNIPROT = {}
TARGET_CHAIN_UNIPROT = {}

CHAIN_COLUMNS = [
    "pdb_id", "assembly_id", "label_asym_id", "auth_asym_id",
    "original_label_asym_id", "original_auth_asym_id",
    "assembly_entity_id", "original_entity_id",
    "n_residues", "n_heavy_atoms","n_backbone_atoms", "chain_uniprot_accessions",
    "target_uniprot_accessions", "is_target_chain",
]

PAIR_COLUMNS = [
    "pdb_id", "assembly_id", "target_uniprot_canonical",
    "pair_type", "target_label_asym_id",
    "target_original_auth_asym_id", "partner_label_asym_id",
    "partner_original_auth_asym_id",
    "partner_uniprot_accessions", "target_entity_id",
    "partner_entity_id", "n_protein_chains",
    "initial_min_backbone_distance",
    "initial_min_heavy_atom_distance",
    "pinder_pair_10A",
    "heavy_contact_pair_5A",
    "best_resolution",
    "experimental_method",
]

STATUS_COLUMNS = [
    "pdb_id", "assembly_id", "status", "n_protein_chains",
    "n_target_chains", "n_contacting_pairs", "error",
]


def values(block, tag):
    return [str(value) for value in block.find_values(tag)]


def optional(items, index, default=""):
    if index >= len(items):
        return default

    value = str(items[index])

    if value in {"", ".", "?"}:
        return default

    return value


def clean_id(value):
    value = str(value)
    return value[:-2] if value.endswith(".0") else value


def model_sort_key(value):
    try:
        return 0, int(value)
    except (TypeError, ValueError):
        return 1, str(value)


def mapped_accessions(pdb_id, auth_value, mapping):
    accessions = set()

    for auth_id in str(auth_value).split(";"):
        auth_id = auth_id.strip()

        if auth_id:
            accessions.update(
                mapping.get((pdb_id, auth_id), set())
            )

    return accessions


def bounding_box_distance(first, second):
    first_min = first.min(axis=0)
    first_max = first.max(axis=0)
    second_min = second.min(axis=0)
    second_max = second.max(axis=0)

    gap = np.maximum(
        np.maximum(first_min - second_max, second_min - first_max),
        0,
    )

    return float(np.linalg.norm(gap))

def minimum_distance_within(tree, coordinates, cutoff):
    distances, _ = tree.query(
        coordinates,
        k=1,
        distance_upper_bound=cutoff,
    )

    finite = distances[np.isfinite(distances)]

    if not len(finite):
        return None

    return float(finite.min())


def parse_assembly(path):
    with gzip.open(
        path,
        "rt",
        encoding="utf-8",
        errors="replace",
    ) as handle:
        document = gemmi.cif.read_string(handle.read())

    block = document.sole_block()

    remapped_entities = values(
        block, "_pdbx_entity_remapping.entity_id"
    )
    original_entities = values(
        block, "_pdbx_entity_remapping.orig_entity_id"
    )

    entity_to_original = {
        clean_id(entity): clean_id(original)
        for entity, original in zip(
            remapped_entities,
            original_entities,
        )
    }

    entity_ids = values(block, "_entity_poly.entity_id")
    entity_types = values(block, "_entity_poly.type")

    protein_entities = {
        clean_id(entity)
        for entity, entity_type in zip(entity_ids, entity_types)
        if "polypeptide" in entity_type.lower()
    }

    label_asym = values(block, "_atom_site.label_asym_id")
    auth_asym = values(block, "_atom_site.auth_asym_id")
    atom_entity = values(block, "_atom_site.label_entity_id")
    label_seq = values(block, "_atom_site.label_seq_id")
    elements = values(block, "_atom_site.type_symbol")
    atom_names = values(block, "_atom_site.label_atom_id")
    x_values = values(block, "_atom_site.Cartn_x")
    y_values = values(block, "_atom_site.Cartn_y")
    z_values = values(block, "_atom_site.Cartn_z")
    models = values(block, "_atom_site.pdbx_PDB_model_num")

    if not label_asym:
        raise ValueError("Missing atom_site records")

    if not models:
        models = ["1"] * len(label_asym)

    model_used = sorted(
        set(models),
        key=model_sort_key,
    )[0]

    chain_atoms = defaultdict(
        lambda: {
            "entity_id": "",
            "auth_ids": set(),
            "residues": set(),
            "coordinates": [],
            "backbone_coordinates": [],
        }
    )

    for index, label_id in enumerate(label_asym):
        if optional(models, index, "1") != model_used:
            continue

        element = optional(elements, index).upper()

        if element in {"H", "D"}:
            continue

        try:
            coordinate = (
                float(x_values[index]),
                float(y_values[index]),
                float(z_values[index]),
            )
        except (ValueError, IndexError):
            continue

        label_id = clean_id(label_id)
        record = chain_atoms[label_id]

        record["entity_id"] = clean_id(
            optional(atom_entity, index)
        )

        auth_id = optional(auth_asym, index)

        if auth_id:
            record["auth_ids"].add(auth_id)

        seq_id = optional(label_seq, index)

        if seq_id:
            record["residues"].add(seq_id)

        record["coordinates"].append(coordinate)
        
        atom_name = optional(atom_names, index).upper()

        if atom_name in BACKBONE_ATOMS:
            record["backbone_coordinates"].append(coordinate)

    remap_entities = values(
        block, "_pdbx_chain_remapping.entity_id"
    )
    remap_labels = values(
        block, "_pdbx_chain_remapping.label_asym_id"
    )
    remap_auth = values(
        block, "_pdbx_chain_remapping.auth_asym_id"
    )
    remap_original_labels = values(
        block, "_pdbx_chain_remapping.orig_label_asym_id"
    )
    remap_original_auth = values(
        block, "_pdbx_chain_remapping.orig_auth_asym_id"
    )

    chains = []

    if remap_labels:
        for index, label_id in enumerate(remap_labels):
            label_id = clean_id(label_id)
            atom_record = chain_atoms.get(label_id)

            if not atom_record or not atom_record["coordinates"]:
                continue

            entity_id = clean_id(
                optional(
                    remap_entities,
                    index,
                    atom_record["entity_id"],
                )
            )

            if entity_id not in protein_entities:
                continue

            auth_id = optional(
                remap_auth,
                index,
                ";".join(sorted(atom_record["auth_ids"])),
            )

            original_auth = optional(
                remap_original_auth,
                index,
                auth_id,
            )

            chains.append({
                "label_asym_id": label_id,
                "auth_asym_id": auth_id,
                "original_label_asym_id": optional(
                    remap_original_labels,
                    index,
                    label_id,
                ),
                "original_auth_asym_id": original_auth,
                "assembly_entity_id": entity_id,
                "original_entity_id": entity_to_original.get(
                    entity_id,
                    entity_id,
                ),
                "n_residues": len(atom_record["residues"]),
                "coordinates": np.asarray(
                    atom_record["coordinates"],
                    dtype=np.float32,
                ),
                "backbone_coordinates": np.asarray(
                    atom_record["backbone_coordinates"],
                    dtype=np.float32
                ).reshape(-1, 3),
            })

    else:
        for label_id, atom_record in chain_atoms.items():
            entity_id = atom_record["entity_id"]

            if (
                entity_id not in protein_entities
                or not atom_record["coordinates"]
            ):
                continue

            auth_id = ";".join(sorted(atom_record["auth_ids"]))

            chains.append({
                "label_asym_id": label_id,
                "auth_asym_id": auth_id,
                "original_label_asym_id": label_id,
                "original_auth_asym_id": auth_id,
                "assembly_entity_id": entity_id,
                "original_entity_id": entity_to_original.get(
                    entity_id,
                    entity_id,
                ),
                "n_residues": len(atom_record["residues"]),
                "coordinates": np.asarray(
                    atom_record["coordinates"],
                    dtype=np.float32,
                ),
                "backbone_coordinates": np.asarray(
                    atom_record["backbone_coordinates"],
                    dtype=np.float32,
                ).reshape(-1, 3),
            })

    return chains


def process_assembly(record):
    pdb_id = str(record["pdb_id"]).lower()
    assembly_id = str(record["assembly_id"])
    filename = f"{pdb_id}-assembly{assembly_id}.cif.gz"
    path = ASSEMBLY_DIR / filename

    if not path.exists():
        return [], [], {
            "pdb_id": pdb_id,
            "assembly_id": assembly_id,
            "status": "missing_file",
            "n_protein_chains": 0,
            "n_target_chains": 0,
            "n_contacting_pairs": 0,
            "error": "",
        }

    try:
        chains = parse_assembly(path)

        for chain in chains:
            original_auth = chain["original_auth_asym_id"]

            chain["chain_accessions"] = mapped_accessions(
                pdb_id,
                original_auth,
                CHAIN_UNIPROT,
            )

            chain["target_accessions"] = mapped_accessions(
                pdb_id,
                original_auth,
                TARGET_CHAIN_UNIPROT,
            )

        target_chains = [
            chain
            for chain in chains
            if chain["target_accessions"]
        ]

        chain_rows = []

        for chain in chains:
            chain_rows.append({
                "pdb_id": pdb_id,
                "assembly_id": assembly_id,
                "label_asym_id": chain["label_asym_id"],
                "auth_asym_id": chain["auth_asym_id"],
                "original_label_asym_id":
                    chain["original_label_asym_id"],
                "original_auth_asym_id":
                    chain["original_auth_asym_id"],
                "assembly_entity_id":
                    chain["assembly_entity_id"],
                "original_entity_id":
                    chain["original_entity_id"],
                "n_residues": chain["n_residues"],
                "n_heavy_atoms": len(chain["coordinates"]),
                "n_backbone_atoms": len(chain["backbone_coordinates"]),
                "chain_uniprot_accessions":
                    ";".join(sorted(chain["chain_accessions"])),
                "target_uniprot_accessions":
                    ";".join(sorted(chain["target_accessions"])),
                "is_target_chain":
                    bool(chain["target_accessions"]),
            })

        pair_rows = []

        for target in target_chains:
            target_heavy_tree = cKDTree(
                target["coordinates"]
            )

            target_backbone_tree = (
                cKDTree(target["backbone_coordinates"])
                if len(target["backbone_coordinates"])
                else None
            )

            for partner in chains:
                if (
                    partner["label_asym_id"]
                    == target["label_asym_id"]
                ):
                    continue

                backbone_distance = None
                heavy_atom_distance = None

                if (
                    target_backbone_tree is not None
                    and len(partner["backbone_coordinates"])
                    and bounding_box_distance(
                        target["backbone_coordinates"],
                        partner["backbone_coordinates"],
                    ) <= PINDER_CUTOFF
                ):
                    backbone_distance = minimum_distance_within(
                        target_backbone_tree,
                        partner["backbone_coordinates"],
                        PINDER_CUTOFF,
                    )

                if (
                    bounding_box_distance(
                        target["coordinates"],
                        partner["coordinates"],
                    ) <= HEAVY_CUTOFF
                ):
                    heavy_atom_distance = minimum_distance_within(
                        target_heavy_tree,
                        partner["coordinates"],
                        HEAVY_CUTOFF,
                    )

                pinder_pair = backbone_distance is not None
                heavy_contact_pair = heavy_atom_distance is not None

                if not (pinder_pair or heavy_contact_pair):
                    continue

                pair_type = (
                    "homomeric"
                    if target["original_entity_id"]
                    == partner["original_entity_id"]
                    else "heteromeric"
                )

                for accession in sorted(
                    target["target_accessions"]
                ):
                    pair_rows.append({
                        "pdb_id": pdb_id,
                        "assembly_id": assembly_id,
                        "target_uniprot_canonical": accession,
                        "pair_type": pair_type,
                        "target_label_asym_id":
                            target["label_asym_id"],
                        "target_original_auth_asym_id":
                            target["original_auth_asym_id"],
                        "partner_label_asym_id":
                            partner["label_asym_id"],
                        "partner_original_auth_asym_id":
                            partner["original_auth_asym_id"],
                        "partner_uniprot_accessions":
                            ";".join(
                                sorted(
                                    partner["chain_accessions"]
                                )
                            ),
                        "target_entity_id":
                            target["original_entity_id"],
                        "partner_entity_id":
                            partner["original_entity_id"],
                        "n_protein_chains": len(chains),
                        "initial_min_backbone_distance":
                            backbone_distance,
                        "initial_min_heavy_atom_distance":
                            heavy_atom_distance,
                        "pinder_pair_10A":
                            pinder_pair,
                        "heavy_contact_pair_5A":
                            heavy_contact_pair,
                        "best_resolution":
                            record.get("best_resolution"),
                        "experimental_method":
                            record.get("experimental_method"),
                    })

        if not target_chains:
            state = "target_chain_not_mapped"
        elif not pair_rows:
            state = "no_contacting_partner"
        else:
            state = "ok"

        return chain_rows, pair_rows, {
            "pdb_id": pdb_id,
            "assembly_id": assembly_id,
            "status": state,
            "n_protein_chains": len(chains),
            "n_target_chains": len(target_chains),
            "n_contacting_pairs": len(pair_rows),
            "error": "",
        }

    except Exception as error:
        return [], [], {
            "pdb_id": pdb_id,
            "assembly_id": assembly_id,
            "status": "parse_error",
            "n_protein_chains": 0,
            "n_target_chains": 0,
            "n_contacting_pairs": 0,
            "error": str(error),
        }


def build_mapping(frame, pdb_column, chain_column, accession_column):
    mapping = defaultdict(set)

    for record in frame.itertuples(index=False):
        mapping[
            (
                str(getattr(record, pdb_column)).lower(),
                str(getattr(record, chain_column)),
            )
        ].add(str(getattr(record, accession_column)))

    return dict(mapping)


def main():
    global CHAIN_UNIPROT, TARGET_CHAIN_UNIPROT

    for directory in [
        CHAIN_PARTS,
        PAIR_PARTS,
        STATUS_PARTS,
    ]:
        directory.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(
        MANIFEST_FILE,
        sep="\t",
        dtype={
            "pdb_id": str,
            "assembly_id": str,
        },
    ).drop_duplicates(["pdb_id", "assembly_id"])

    manifest = manifest.sort_values(
        ["pdb_id", "assembly_id"]
    )

    if LIMIT:
        manifest = manifest.head(LIMIT)

    selected_pdbs = set(
        manifest["pdb_id"].astype(str).str.lower()
    )

    full_sifts = pd.read_csv(
        FULL_SIFTS_FILE,
        sep="\t",
        comment="#",
        dtype=str,
        usecols=["PDB", "CHAIN", "SP_PRIMARY"],
    )

    full_sifts["PDB"] = full_sifts["PDB"].str.lower()
    full_sifts = full_sifts[
        full_sifts["PDB"].isin(selected_pdbs)
    ]

    target_segments = pd.read_parquet(
        TARGET_SEGMENT_FILE,
        columns=[
            "pdb_id",
            "auth_chain_id",
            "uniprot_canonical",
        ],
    )

    target_segments["pdb_id"] = (
        target_segments["pdb_id"].str.lower()
    )

    CHAIN_UNIPROT = build_mapping(
        full_sifts,
        "PDB",
        "CHAIN",
        "SP_PRIMARY",
    )

    TARGET_CHAIN_UNIPROT = build_mapping(
        target_segments,
        "pdb_id",
        "auth_chain_id",
        "uniprot_canonical",
    )

    records = manifest.to_dict("records")
    chunks = [
        records[index:index + CHUNK_SIZE]
        for index in range(0, len(records), CHUNK_SIZE)
    ]

    with ProcessPoolExecutor(
        max_workers=WORKERS
    ) as executor:
        for chunk_index, chunk in enumerate(
            tqdm(chunks, desc="Assembly chunks")
        ):
            chain_path = (
                CHAIN_PARTS / f"part_{chunk_index:05d}.parquet"
            )
            pair_path = (
                PAIR_PARTS / f"part_{chunk_index:05d}.parquet"
            )
            status_path = (
                STATUS_PARTS / f"part_{chunk_index:05d}.parquet"
            )

            if (
                chain_path.exists()
                and pair_path.exists()
                and status_path.exists()
            ):
                continue

            futures = [
                executor.submit(process_assembly, record)
                for record in chunk
            ]

            chain_rows = []
            pair_rows = []
            status_rows = []

            for future in as_completed(futures):
                chains, pairs, status = future.result()
                chain_rows.extend(chains)
                pair_rows.extend(pairs)
                status_rows.append(status)

            pd.DataFrame(
                chain_rows,
                columns=CHAIN_COLUMNS,
            ).to_parquet(
                chain_path,
                index=False,
                compression="zstd",
            )

            pd.DataFrame(
                pair_rows,
                columns=PAIR_COLUMNS,
            ).to_parquet(
                pair_path,
                index=False,
                compression="zstd",
            )

            pd.DataFrame(
                status_rows,
                columns=STATUS_COLUMNS,
            ).to_parquet(
                status_path,
                index=False,
                compression="zstd",
            )

    status_files = sorted(STATUS_PARTS.glob("part_*.parquet"))
    statuses = pd.concat(
        [pd.read_parquet(path) for path in status_files],
        ignore_index=True,
    )

    summary = {
        "assemblies_requested": len(manifest),
        "assemblies_processed": len(statuses),
        "status_counts": {
            str(key): int(value)
            for key, value in statuses[
                "status"
            ].value_counts().items()
        },
        "protein_chains": sum(
            len(pd.read_parquet(path))
            for path in CHAIN_PARTS.glob("part_*.parquet")
        ),
        "contacting_target_partner_pairs": sum(
            len(pd.read_parquet(path))
            for path in PAIR_PARTS.glob("part_*.parquet")
        ),
    }

    with open(SUMMARY_OUTPUT, "w") as handle:
        json.dump(summary, handle, indent=2)

    print()
    print(json.dumps(summary, indent=2))
    print()
    print("Chain dataset:", CHAIN_PARTS)
    print("Pair dataset:", PAIR_PARTS)


if __name__ == "__main__":
    main()