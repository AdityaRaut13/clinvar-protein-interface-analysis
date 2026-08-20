#!/usr/bin/env python

import json
from pathlib import Path

import pandas as pd


ROOT = Path(os.environ["IF_PROJECT"])

TARGET_FILE = (
    ROOT
    / "data/processed/uniprot_ready"
    / "unique_uniprot_structural_targets.tsv.gz"
)
SIFTS_FILE = (
    ROOT
    / "data/input/sifts"
    / "pdb_chain_uniprot.tsv.gz"
)

OUTPUT_DIR = ROOT / "data/processed/structure_discovery"

CHAIN_OUTPUT = OUTPUT_DIR / "target_pdb_chains.parquet"
PDB_OUTPUT = OUTPUT_DIR / "unique_target_pdb_entries.tsv.gz"
COVERAGE_OUTPUT = OUTPUT_DIR / "target_pdb_coverage.tsv.gz"
SUMMARY_OUTPUT = OUTPUT_DIR / "sifts_pdb_manifest_summary.json"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def join_unique(values):
    return ";".join(sorted({
        str(value)
        for value in values.dropna()
        if str(value).strip()
    }))


targets = pd.read_csv(
    TARGET_FILE,
    sep="\t",
    dtype={"uniprot_canonical": str},
)

sifts = pd.read_csv(
    SIFTS_FILE,
    sep="\t",
    comment="#",
    dtype=str,
)

sifts.columns = [
    str(column).strip().upper()
    for column in sifts.columns
]

required = {
    "PDB",
    "CHAIN",
    "SP_PRIMARY",
    "RES_BEG",
    "RES_END",
    "PDB_BEG",
    "PDB_END",
    "SP_BEG",
    "SP_END",
}

missing = required - set(sifts.columns)

if missing:
    raise ValueError(
        f"Missing SIFTS columns: {sorted(missing)}; "
        f"found: {list(sifts.columns)}"
    )

target_accessions = set(
    targets["uniprot_canonical"].dropna()
)

overlap = sifts[
    sifts["SP_PRIMARY"].isin(target_accessions)
].copy()

overlap = overlap.rename(
    columns={
        "PDB": "pdb_id",
        "CHAIN": "auth_chain_id",
        "SP_PRIMARY": "uniprot_canonical",
        "RES_BEG": "pdb_residue_begin",
        "RES_END": "pdb_residue_end",
        "PDB_BEG": "pdb_sequence_begin",
        "PDB_END": "pdb_sequence_end",
        "SP_BEG": "uniprot_begin",
        "SP_END": "uniprot_end",
    }
)

overlap["pdb_id"] = overlap["pdb_id"].str.lower()

overlap = overlap.drop_duplicates(
    [
        "pdb_id",
        "auth_chain_id",
        "uniprot_canonical",
        "uniprot_begin",
        "uniprot_end",
    ]
)

overlap.to_parquet(
    CHAIN_OUTPUT,
    index=False,
    compression="zstd",
)

pdb_manifest = (
    overlap.groupby("pdb_id", as_index=False)
    .agg(
        target_accessions=(
            "uniprot_canonical",
            join_unique,
        ),
        n_target_accessions=(
            "uniprot_canonical",
            "nunique",
        ),
        n_target_auth_chains=(
            "auth_chain_id",
            "nunique",
        ),
    )
    .sort_values("pdb_id")
)

pdb_manifest.to_csv(
    PDB_OUTPUT,
    sep="\t",
    index=False,
    compression="gzip",
)

accession_coverage = (
    overlap.groupby("uniprot_canonical", as_index=False)
    .agg(
        n_pdb_entries=("pdb_id", "nunique"),
        n_sifts_chains=("auth_chain_id", "nunique"),
    )
)

coverage = targets.merge(
    accession_coverage,
    on="uniprot_canonical",
    how="left",
)

coverage["n_pdb_entries"] = (
    coverage["n_pdb_entries"].fillna(0).astype(int)
)
coverage["n_sifts_chains"] = (
    coverage["n_sifts_chains"].fillna(0).astype(int)
)
coverage["has_experimental_pdb"] = (
    coverage["n_pdb_entries"].gt(0)
)

coverage.to_csv(
    COVERAGE_OUTPUT,
    sep="\t",
    index=False,
    compression="gzip",
)

has_pdb = coverage["has_experimental_pdb"]

summary = {
    "structural_targets": int(len(targets)),
    "targets_with_experimental_pdb": int(has_pdb.sum()),
    "targets_without_experimental_pdb": int((~has_pdb).sum()),
    "unique_target_pdb_entries": int(len(pdb_manifest)),
    "target_pdb_chain_mapping_rows": int(len(overlap)),
    "ready_variants": int(targets["n_variants"].sum()),
    "variants_on_targets_with_pdb": int(
        coverage.loc[has_pdb, "n_variants"].sum()
    ),
    "variants_on_targets_without_pdb": int(
        coverage.loc[~has_pdb, "n_variants"].sum()
    ),
}

with open(SUMMARY_OUTPUT, "w") as handle:
    json.dump(summary, handle, indent=2)

print(json.dumps(summary, indent=2))
print()
print("Target chains:", CHAIN_OUTPUT)
print("PDB manifest:", PDB_OUTPUT)
print("Coverage:", COVERAGE_OUTPUT)