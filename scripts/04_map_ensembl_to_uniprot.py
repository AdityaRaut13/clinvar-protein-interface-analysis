#!/usr/bin/env python

import io
import json
import time
from pathlib import Path

import pandas as pd
import requests


ROOT = Path(os.environ["IF_PROJECT"])

STATUS_FILE = (
    ROOT
    / "data/processed/uniprot_mapping"
    / "protein_uniprot_exact_status.tsv.gz"
)
VARIANT_FILE = (
    ROOT
    / "data/processed/source_validation"
    / "clinvar_source_valid.parquet"
)
ANNOTATED_CANDIDATES = (
    ROOT
    / "data/processed/uniprot_mapping"
    / "protein_uniprot_candidates.parquet"
)

OUTPUT_DIR = ROOT / "data/processed/uniprot_mapping"
RAW_OUTPUT = OUTPUT_DIR / "ensembl_uniprot_raw.tsv"
JOB_OUTPUT = OUTPUT_DIR / "ensembl_uniprot_job.json"
MAPPING_OUTPUT = OUTPUT_DIR / "ensembl_uniprot_mapping.parquet"
COMBINED_OUTPUT = OUTPUT_DIR / "protein_uniprot_candidates_combined.parquet"
SUMMARY_OUTPUT = OUTPUT_DIR / "ensembl_uniprot_mapping_summary.json"

BASE_URL = "https://rest.uniprot.org"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def canonical(accession):
    return str(accession).split("-")[0]


status = pd.read_csv(
    STATUS_FILE,
    sep="\t",
    dtype={"source_sequence_key": str},
)

unresolved = status[
    ~status["mapping_status"].eq("unique_exact_canonical")
].copy()

ensembl_ids = sorted(
    unresolved["source_sequence_key"]
    .dropna()
    .astype(str)
    .unique()
)

invalid_ids = [
    value
    for value in ensembl_ids
    if not value.startswith("ENSP")
]

if invalid_ids:
    raise ValueError(
        f"Non-Ensembl source identifiers found: {invalid_ids[:10]}"
    )

print(f"Unresolved Ensembl proteins: {len(ensembl_ids):,}")

session = requests.Session()
session.headers.update({
    "User-Agent": "IITK-ClinVar-interface/1.0"
})

if not RAW_OUTPUT.exists():
    if JOB_OUTPUT.exists():
        with open(JOB_OUTPUT) as handle:
            job = json.load(handle)

        job_id = job["jobId"]
        print("Resuming UniProt job:", job_id)

    else:
        response = session.post(
            f"{BASE_URL}/idmapping/run",
            data={
                "from": "Ensembl_Protein",
                "to": "UniProtKB",
                "ids": ",".join(ensembl_ids),
            },
            timeout=180,
        )
        response.raise_for_status()

        job = response.json()
        job_id = job["jobId"]

        with open(JOB_OUTPUT, "w") as handle:
            json.dump(job, handle, indent=2)

        print("Submitted UniProt job:", job_id)

    while True:
        response = session.get(
            f"{BASE_URL}/idmapping/status/{job_id}",
            timeout=180,
        )
        response.raise_for_status()

        job_status = response.json()
        state = job_status.get("jobStatus")

        if state in {"RUNNING","NEW"}:
            print("Mapping job still running")
            time.sleep(5)
            continue

        if state and state not in {"FINISHED"}:
            raise RuntimeError(
                f"UniProt mapping failed: {job_status}"
            )

        if (
            "results" in job_status
            or "failedIds" in job_status
            or state == "FINISHED"
        ):
            break

        time.sleep(5)

    response = session.get(
        f"{BASE_URL}/idmapping/stream/{job_id}",
        params={
            "format": "tsv",
        },
        timeout=600,
    )
    response.raise_for_status()
    RAW_OUTPUT.write_text(response.text)

else:
    print("Using cached raw mapping:", RAW_OUTPUT)

raw = pd.read_csv(
    RAW_OUTPUT,
    sep="\t",
    dtype=str,
)

print("Raw columns:", list(raw.columns))

if len(raw.columns) < 2:
    raise ValueError(
        f"Unexpected UniProt mapping output: {list(raw.columns)}"
    )

raw = raw.iloc[:, :2].copy()
raw.columns = [
    "source_sequence_key",
    "uniprot_accession",
]

raw = (
    raw.dropna()
    .drop_duplicates()
)

raw["uniprot_canonical"] = (
    raw["uniprot_accession"].map(canonical)
)
raw["is_isoform_accession"] = (
    raw["uniprot_accession"].str.contains("-", regex=False)
)
raw["annotation_sources"] = "ensembl_id_mapping"
raw["annotation_priority"] = 4
raw["protein_ids"] = raw["source_sequence_key"]

variants = pd.read_parquet(
    VARIANT_FILE,
    columns=["source_sequence_key", "gene_symbol"],
)

gene_map = (
    variants.groupby("source_sequence_key")["gene_symbol"]
    .agg(
        lambda values: ";".join(
            sorted(set(values.dropna().astype(str)))
        )
    )
    .to_dict()
)

raw["gene_symbols"] = (
    raw["source_sequence_key"]
    .map(gene_map)
    .fillna("")
)

raw.to_parquet(
    MAPPING_OUTPUT,
    index=False,
    compression="zstd",
)

annotated = pd.read_parquet(ANNOTATED_CANDIDATES)

combined = pd.concat(
    [annotated, raw],
    ignore_index=True,
)

combined = (
    combined.sort_values(
        [
            "source_sequence_key",
            "annotation_priority",
            "uniprot_accession",
        ]
    )
    .drop_duplicates(
        [
            "source_sequence_key",
            "uniprot_accession",
        ],
        keep="first",
    )
)

combined.to_parquet(
    COMBINED_OUTPUT,
    index=False,
    compression="zstd",
)

mapped_ids = set(raw["source_sequence_key"])

summary = {
    "submitted_ensembl_proteins": int(len(ensembl_ids)),
    "mapped_ensembl_proteins": int(len(mapped_ids)),
    "unmapped_ensembl_proteins": int(
        len(set(ensembl_ids) - mapped_ids)
    ),
    "ensembl_uniprot_mapping_rows": int(len(raw)),
    "unique_mapped_uniprot_accessions": int(
        raw["uniprot_accession"].nunique()
    ),
    "original_candidate_rows": int(len(annotated)),
    "combined_candidate_rows": int(len(combined)),
}

with open(SUMMARY_OUTPUT, "w") as handle:
    json.dump(summary, handle, indent=2)

print()
print(json.dumps(summary, indent=2))
print()
print("Raw mapping:", RAW_OUTPUT)
print("Mapping:", MAPPING_OUTPUT)
print("Combined candidates:", COMBINED_OUTPUT)