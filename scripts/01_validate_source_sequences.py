#!/usr/bin/env python

import gzip
import hashlib
import json
from pathlib import Path

import pandas as pd
from Bio import SeqIO
import os


ROOT = Path(os.environ["IF_PROJECT"])

VARIANT_FILE = (
    ROOT / "data/input/clinvar_missense_protein_ready.csv.gz"
)
FASTA_FILE = (
    ROOT / "data/input/protein_reference_sequences.fasta.gz"
)
OUTPUT_DIR = ROOT / "data/processed/source_validation"

VALID_OUTPUT = OUTPUT_DIR / "clinvar_source_valid.parquet"
FAILURE_OUTPUT = OUTPUT_DIR / "clinvar_source_failures.tsv.gz"
PROTEIN_OUTPUT = OUTPUT_DIR / "unique_source_proteins.parquet"
SUMMARY_OUTPUT = OUTPUT_DIR / "source_validation_summary.json"

VALID_AA = set("ACDEFGHIKLMNPQRSTVWYOU")

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def clean_text(series):
    return series.astype("string").str.strip()


def sequence_hash(sequence):
    return hashlib.sha256(sequence.encode()).hexdigest()


print("Loading FASTA")

sequences = {}

with gzip.open(FASTA_FILE, "rt") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        record_id = str(record.id).strip()
        sequence = str(record.seq).upper()

        if record_id in sequences:
            raise ValueError(f"Duplicate FASTA identifier: {record_id}")

        sequences[record_id] = sequence

if not sequences:
    raise ValueError("No protein sequences found in FASTA")

fasta_ids = set(sequences)

print(f"FASTA sequences: {len(sequences):,}")
print("Loading ClinVar variants")

variants = pd.read_csv(
    VARIANT_FILE,
    dtype={
        "variant_id": str,
        "chrom": str,
        "protein_id": str,
        "protein_sequence_id": str,
        "protein_position": str,
        "aa_ref": str,
        "aa_alt": str,
    },
)

for column in variants.select_dtypes(include=["object"]).columns:
    variants[column] = variants[column].astype("string")

for column in [
    "variant_id",
    "protein_id",
    "protein_sequence_id",
    "aa_ref",
    "aa_alt",
]:
    variants[column] = clean_text(variants[column])

variants["aa_ref"] = variants["aa_ref"].str.upper()
variants["aa_alt"] = variants["aa_alt"].str.upper()

candidate_columns = ["protein_sequence_id", "protein_id"]
match_counts = {
    column: variants[column].isin(fasta_ids).sum()
    for column in candidate_columns
}

fasta_key_column = max(match_counts, key=match_counts.get)

if match_counts[fasta_key_column] == 0:
    raise ValueError(
        "Neither protein_sequence_id nor protein_id matches FASTA IDs"
    )

print("FASTA ID matches:", match_counts)
print("Selected FASTA key:", fasta_key_column)

variants["source_sequence_key"] = variants[fasta_key_column]

variants["protein_position_numeric"] = pd.to_numeric(
    variants["protein_position"],
    errors="coerce",
)

variants["validation_status"] = "valid"

variants.loc[
    variants["source_sequence_key"].isna(),
    "validation_status",
] = "missing_sequence_identifier"

variants.loc[
    variants["validation_status"].eq("valid")
    & ~variants["source_sequence_key"].isin(fasta_ids),
    "validation_status",
] = "sequence_not_found"

variants.loc[
    variants["validation_status"].eq("valid")
    & variants["protein_position_numeric"].isna(),
    "validation_status",
] = "invalid_protein_position"

variants.loc[
    variants["validation_status"].eq("valid")
    & (variants["protein_position_numeric"] % 1 != 0),
    "validation_status",
] = "non_integer_protein_position"

variants.loc[
    variants["validation_status"].eq("valid")
    & ~variants["aa_ref"].isin(VALID_AA),
    "validation_status",
] = "invalid_reference_amino_acid"

variants.loc[
    variants["validation_status"].eq("valid")
    & ~variants["aa_alt"].isin(VALID_AA),
    "validation_status",
] = "invalid_alternate_amino_acid"

variants.loc[
    variants["validation_status"].eq("valid")
    & variants["aa_ref"].eq(variants["aa_alt"]),
    "validation_status",
] = "identical_reference_and_alternate"

variants["source_sequence_length"] = variants[
    "source_sequence_key"
].map(
    lambda key: len(sequences[key]) if key in sequences else pd.NA
)

variants.loc[
    variants["validation_status"].eq("valid")
    & (
        variants["protein_position_numeric"].lt(1)
        | variants["protein_position_numeric"].gt(
            variants["source_sequence_length"]
        )
    ),
    "validation_status",
] = "position_out_of_range"

reported_length = pd.to_numeric(
    variants["protein_length"],
    errors="coerce",
)

variants.loc[
    variants["validation_status"].eq("valid")
    & reported_length.notna()
    & reported_length.ne(variants["source_sequence_length"]),
    "validation_status",
] = "protein_length_mismatch"


def residue_at_position(row):
    if row["validation_status"] != "valid":
        return pd.NA

    sequence = sequences[row["source_sequence_key"]]
    position = int(row["protein_position_numeric"])
    return sequence[position - 1]


variants["sequence_reference_aa"] = variants.apply(
    residue_at_position,
    axis=1,
)

variants.loc[
    variants["validation_status"].eq("valid")
    & variants["sequence_reference_aa"].ne(variants["aa_ref"]),
    "validation_status",
] = "reference_residue_mismatch"

valid = variants[
    variants["validation_status"].eq("valid")
].copy()

valid["protein_position"] = (
    valid["protein_position_numeric"].astype("int32")
)

valid["source_sequence_length"] = (
    valid["source_sequence_length"].astype("int32")
)

valid["sequence_hash"] = valid["source_sequence_key"].map(
    lambda key: sequence_hash(sequences[key])
)

valid["variant_key"] = (
    valid["variant_id"].astype(str)
    + "|"
    + valid["protein_id"].astype(str)
    + "|"
    + valid["protein_position"].astype(str)
    + "|"
    + valid["aa_ref"]
    + ">"
    + valid["aa_alt"]
)

if valid["variant_key"].duplicated().any():
    duplicated = valid.loc[
        valid["variant_key"].duplicated(keep=False),
        [
            "variant_key",
            "variant_id",
            "protein_id",
            "protein_position",
            "aa_ref",
            "aa_alt",
        ],
    ]

    duplicate_path = (
        OUTPUT_DIR / "duplicate_variant_keys.tsv.gz"
    )

    duplicated.to_csv(
        duplicate_path,
        sep="\t",
        index=False,
        compression="gzip",
    )

    raise ValueError(
        f"Duplicate variant keys detected; inspect {duplicate_path}"
    )

drop_columns = ["protein_position_numeric"]
valid = valid.drop(columns=drop_columns)

failures = variants[
    ~variants["validation_status"].eq("valid")
].copy()

valid.to_parquet(
    VALID_OUTPUT,
    index=False,
    compression="zstd",
)

failures.to_csv(
    FAILURE_OUTPUT,
    sep="\t",
    index=False,
    compression="gzip",
)

protein_rows = []

for key in sorted(valid["source_sequence_key"].unique()):
    sequence = sequences[key]

    protein_rows.append({
        "source_sequence_key": key,
        "sequence_hash": sequence_hash(sequence),
        "protein_sequence": sequence,
        "sequence_length": len(sequence),
    })

proteins = pd.DataFrame(protein_rows)

protein_metadata = (
    valid.groupby("source_sequence_key", as_index=False)
    .agg(
        n_variants=("variant_key", "nunique"),
        n_protein_ids=("protein_id", "nunique"),
        n_genes=("gene_symbol", "nunique"),
        n_benign=("label", lambda values: values.eq(0).sum()),
        n_pathogenic=("label", lambda values: values.eq(1).sum()),
    )
)

proteins = proteins.merge(
    protein_metadata,
    on="source_sequence_key",
    how="left",
)

proteins.to_parquet(
    PROTEIN_OUTPUT,
    index=False,
    compression="zstd",
)

status_counts = (
    variants["validation_status"]
    .value_counts(dropna=False)
    .to_dict()
)

summary = {
    "input_variants": int(len(variants)),
    "valid_variants": int(len(valid)),
    "failed_variants": int(len(failures)),
    "unique_valid_source_sequences": int(len(proteins)),
    "unique_valid_protein_ids": int(valid["protein_id"].nunique()),
    "fasta_records": int(len(sequences)),
    "selected_fasta_key_column": fasta_key_column,
    "fasta_match_counts": {
        key: int(value)
        for key, value in match_counts.items()
    },
    "validation_status_counts": {
        str(key): int(value)
        for key, value in status_counts.items()
    },
}

with open(SUMMARY_OUTPUT, "w") as handle:
    json.dump(summary, handle, indent=2)

print()
print(json.dumps(summary, indent=2))
print()
print("Valid variants:", VALID_OUTPUT)
print("Failures:", FAILURE_OUTPUT)
print("Unique proteins:", PROTEIN_OUTPUT)
print("Summary:", SUMMARY_OUTPUT)
