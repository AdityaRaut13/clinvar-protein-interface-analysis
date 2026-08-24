#!/usr/bin/env python

import json
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd
import os


ROOT = Path(os.environ["IF_PROJECT"])

VARIANT_FILE = (
    ROOT
    / "data/processed/source_validation"
    / "clinvar_source_valid.parquet"
)
PROTEIN_FILE = (
    ROOT
    / "data/processed/source_validation"
    / "unique_source_proteins.parquet"
)
OUTPUT_DIR = ROOT / "data/processed/uniprot_mapping"

CANDIDATE_OUTPUT = OUTPUT_DIR / "protein_uniprot_candidates.parquet"
COVERAGE_OUTPUT = OUTPUT_DIR / "protein_uniprot_candidate_coverage.tsv.gz"
INVALID_OUTPUT = OUTPUT_DIR / "invalid_uniprot_annotations.tsv.gz"
SUMMARY_OUTPUT = OUTPUT_DIR / "uniprot_candidate_summary.json"

UNIPROT_RE = re.compile(
    r"^(?:"
    r"[OPQ][0-9][A-Z0-9]{3}[0-9]|"
    r"[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}"
    r")(?:-\d+)?$"
)

MAPPING_FIELDS = {
    "uniprot_isoform": 1,
    "swissprot": 2,
    "trembl": 3,
}

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def split_accessions(value):
    if pd.isna(value):
        return []
    accessions = set()

    for item in re.split(r"[\s,;&|]+", str(value)):
        item = item.strip().upper()

        if not item:
            continue

        item = re.sub(r"\.\d+$", "", item)
        accessions.add(item)

    return sorted(accessions)


variants = pd.read_parquet(VARIANT_FILE)
proteins = pd.read_parquet(PROTEIN_FILE)

missing = set(MAPPING_FIELDS) - set(variants.columns)

if missing:
    raise ValueError(f"Missing mapping columns: {sorted(missing)}")

annotation_columns = [
    "source_sequence_key",
    "protein_id",
    "gene_symbol",
    "sequence_hash",
    *MAPPING_FIELDS,
]

annotations = variants[annotation_columns].drop_duplicates()

candidate_data = defaultdict(
    lambda: {
        "sources": set(),
        "priorities": [],
        "protein_ids": set(),
        "gene_symbols": set(),
    }
)

invalid_data = set()

for record in annotations.itertuples(index=False):
    source_key = str(record.source_sequence_key)

    for field, priority in MAPPING_FIELDS.items():
        for accession in split_accessions(getattr(record, field)):
            if not UNIPROT_RE.fullmatch(accession):
                invalid_data.add((source_key, field, accession))
                continue

            key = (source_key, accession)
            candidate_data[key]["sources"].add(field)
            candidate_data[key]["priorities"].append(priority)
            candidate_data[key]["protein_ids"].add(
                str(record.protein_id)
            )
            candidate_data[key]["gene_symbols"].add(
                str(record.gene_symbol)
            )

candidate_rows = []

for (source_key, accession), data in candidate_data.items():
    candidate_rows.append({
        "source_sequence_key": source_key,
        "uniprot_accession": accession,
        "uniprot_canonical": accession.split("-")[0],
        "is_isoform_accession": "-" in accession,
        "annotation_sources": ";".join(sorted(data["sources"])),
        "annotation_priority": min(data["priorities"]),
        "protein_ids": ";".join(sorted(data["protein_ids"])),
        "gene_symbols": ";".join(sorted(data["gene_symbols"])),
    })

candidates = pd.DataFrame(candidate_rows)

if len(candidates):
    candidates = candidates.sort_values(
        [
            "source_sequence_key",
            "annotation_priority",
            "uniprot_accession",
        ]
    )

    candidates.to_parquet(
        CANDIDATE_OUTPUT,
        index=False,
        compression="zstd",
    )
else:
    raise ValueError("No valid UniProt candidates were found")

invalid = pd.DataFrame(
    sorted(invalid_data),
    columns=[
        "source_sequence_key",
        "annotation_field",
        "annotation_value",
    ],
)

invalid.to_csv(
    INVALID_OUTPUT,
    sep="\t",
    index=False,
    compression="gzip",
)

candidate_counts = (
    candidates.groupby("source_sequence_key", as_index=False)
    .agg(
        n_uniprot_candidates=("uniprot_accession", "nunique"),
        n_canonical_candidates=("uniprot_canonical", "nunique"),
        has_isoform_candidate=("is_isoform_accession", "max"),
        best_annotation_priority=("annotation_priority", "min"),
    )
)

coverage = proteins[
    [
        "source_sequence_key",
        "sequence_hash",
        "sequence_length",
        "n_variants",
        "n_benign",
        "n_pathogenic",
    ]
].merge(
    candidate_counts,
    on="source_sequence_key",
    how="left",
)

coverage["n_uniprot_candidates"] = (
    coverage["n_uniprot_candidates"].fillna(0).astype(int)
)
coverage["n_canonical_candidates"] = (
    coverage["n_canonical_candidates"].fillna(0).astype(int)
)
coverage["has_isoform_candidate"] = (
    coverage["has_isoform_candidate"].fillna(False).astype(bool)
)

coverage["candidate_status"] = "has_candidate"
coverage.loc[
    coverage["n_uniprot_candidates"].eq(0),
    "candidate_status",
] = "no_annotated_candidate"

coverage.to_csv(
    COVERAGE_OUTPUT,
    sep="\t",
    index=False,
    compression="gzip",
)

mapped_keys = set(candidates["source_sequence_key"])

variant_candidate_coverage = variants[
    "source_sequence_key"
].isin(mapped_keys)

summary = {
    "unique_source_proteins": int(len(proteins)),
    "proteins_with_candidate": int(
        coverage["candidate_status"].eq("has_candidate").sum()
    ),
    "proteins_without_candidate": int(
        coverage["candidate_status"]
        .eq("no_annotated_candidate")
        .sum()
    ),
    "candidate_rows": int(len(candidates)),
    "unique_uniprot_accessions": int(
        candidates["uniprot_accession"].nunique()
    ),
    "unique_canonical_accessions": int(
        candidates["uniprot_canonical"].nunique()
    ),
    "proteins_with_multiple_candidates": int(
        coverage["n_uniprot_candidates"].gt(1).sum()
    ),
    "invalid_annotation_values": int(len(invalid)),
    "variants_with_candidate": int(
        variant_candidate_coverage.sum()
    ),
    "variants_without_candidate": int(
        (~variant_candidate_coverage).sum()
    ),
}

with open(SUMMARY_OUTPUT, "w") as handle:
    json.dump(summary, handle, indent=2)

print(json.dumps(summary, indent=2))
print()
print("Candidates:", CANDIDATE_OUTPUT)
print("Coverage:", COVERAGE_OUTPUT)
print("Invalid annotations:", INVALID_OUTPUT)
print("Summary:", SUMMARY_OUTPUT)
