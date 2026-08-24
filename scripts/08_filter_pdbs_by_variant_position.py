#!/usr/bin/env python

import json
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm
import os


ROOT = Path(os.environ["IF_PROJECT"])

VARIANT_FILE = (
    ROOT
    / "data/processed/uniprot_ready"
    / "clinvar_variants_uniprot_ready.parquet"
)
TARGET_FILE = (
    ROOT
    / "data/processed/uniprot_ready"
    / "unique_uniprot_structural_targets.tsv.gz"
)
SIFTS_FILE = (
    ROOT
    / "data/input/sifts"
    / "uniprot_segments_observed.tsv.gz"
)

OUTPUT_DIR = ROOT / "data/processed/structure_discovery"

SEGMENT_OUTPUT = OUTPUT_DIR / "variant_observed_pdb_segments.parquet"
VARIANT_OUTPUT = OUTPUT_DIR / "clinvar_variants_with_observed_pdb.parquet"
PDB_OUTPUT = OUTPUT_DIR / "variant_observed_pdb_manifest.tsv.gz"
COVERAGE_OUTPUT = OUTPUT_DIR / "variant_observed_target_coverage.tsv.gz"
SUMMARY_OUTPUT = OUTPUT_DIR / "variant_observed_pdb_summary.json"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def join_unique(values):
    return ";".join(sorted({
        str(value)
        for value in values.dropna()
        if str(value).strip()
    }))


variants = pd.read_parquet(VARIANT_FILE)
targets = pd.read_csv(TARGET_FILE, sep="\t")

target_accessions = set(
    targets["uniprot_canonical"].dropna().astype(str)
)

print("Reading observed SIFTS segments")

segments = pd.read_csv(
    SIFTS_FILE,
    sep="\t",
    comment="#",
    dtype=str,
)

segments.columns = [
    str(column).strip().upper()
    for column in segments.columns
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

missing = required - set(segments.columns)

if missing:
    raise ValueError(f"Missing SIFTS columns: {sorted(missing)}")

segments = segments[
    segments["SP_PRIMARY"].isin(target_accessions)
].copy()

segments["SP_BEG"] = pd.to_numeric(
    segments["SP_BEG"],
    errors="coerce",
)
segments["SP_END"] = pd.to_numeric(
    segments["SP_END"],
    errors="coerce",
)

segments = segments.dropna(
    subset=["SP_BEG", "SP_END"]
)

segments["SP_BEG"] = segments["SP_BEG"].astype(int)
segments["SP_END"] = segments["SP_END"].astype(int)

segments["uniprot_start"] = segments[
    ["SP_BEG", "SP_END"]
].min(axis=1)

segments["uniprot_end"] = segments[
    ["SP_BEG", "SP_END"]
].max(axis=1)

positions_by_accession = {
    accession: np.sort(
        group["uniprot_position"]
        .dropna()
        .astype(int)
        .unique()
    )
    for accession, group in variants.groupby(
        "uniprot_canonical",
        sort=False,
    )
}

retained_groups = []
covered_position_pairs = set()

print(
    "Filtering observed segments:",
    f"{len(segments):,} candidate segments",
)

for accession, group in tqdm(
    segments.groupby("SP_PRIMARY", sort=False),
    total=segments["SP_PRIMARY"].nunique(),
):
    positions = positions_by_accession.get(accession)

    if positions is None or len(positions) == 0:
        continue

    starts = group["uniprot_start"].to_numpy()
    ends = group["uniprot_end"].to_numpy()

    left = np.searchsorted(
        positions,
        starts,
        side="left",
    )
    right = np.searchsorted(
        positions,
        ends,
        side="right",
    )

    keep = left < right

    if not keep.any():
        continue

    retained = group.loc[keep].copy()
    retained_left = left[keep]
    retained_right = right[keep]

    retained["n_clinvar_positions_in_segment"] = (
        retained_right - retained_left
    )

    for begin, end in zip(
        retained_left,
        retained_right,
    ):
        for position in positions[begin:end]:
            covered_position_pairs.add(
                (accession, int(position))
            )

    retained_groups.append(retained)

if retained_groups:
    retained = pd.concat(
        retained_groups,
        ignore_index=True,
    )
else:
    raise ValueError("No observed segments overlap ClinVar positions")

retained = retained.rename(
    columns={
        "PDB": "pdb_id",
        "CHAIN": "auth_chain_id",
        "SP_PRIMARY": "uniprot_canonical",
        "RES_BEG": "observed_residue_begin",
        "RES_END": "observed_residue_end",
        "PDB_BEG": "auth_residue_begin",
        "PDB_END": "auth_residue_end",
        "SP_BEG": "uniprot_segment_begin",
        "SP_END": "uniprot_segment_end",
    }
)

retained["pdb_id"] = retained["pdb_id"].str.lower()

retained = retained.drop_duplicates(
    [
        "pdb_id",
        "auth_chain_id",
        "uniprot_canonical",
        "uniprot_segment_begin",
        "uniprot_segment_end",
    ]
)

retained.to_parquet(
    SEGMENT_OUTPUT,
    index=False,
    compression="zstd",
)

variant_has_observed_structure = [
    (
        str(accession),
        int(position),
    ) in covered_position_pairs
    for accession, position in zip(
        variants["uniprot_canonical"],
        variants["uniprot_position"],
    )
]

variant_has_observed_structure = np.asarray(
    variant_has_observed_structure,
    dtype=bool,
)

covered_variants = variants[
    variant_has_observed_structure
].copy()

covered_variants.to_parquet(
    VARIANT_OUTPUT,
    index=False,
    compression="zstd",
)

retained["chain_key"] = (
    retained["pdb_id"]
    + ":"
    + retained["auth_chain_id"]
)

pdb_manifest = (
    retained.groupby("pdb_id", as_index=False)
    .agg(
        target_accessions=(
            "uniprot_canonical",
            join_unique,
        ),
        n_target_accessions=(
            "uniprot_canonical",
            "nunique",
        ),
        n_target_chains=(
            "chain_key",
            "nunique",
        ),
        n_observed_segments=(
            "uniprot_canonical",
            "size",
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
    retained.groupby("uniprot_canonical", as_index=False)
    .agg(
        n_position_relevant_pdb_entries=(
            "pdb_id",
            "nunique",
        ),
        n_position_relevant_chains=(
            "chain_key",
            "nunique",
        ),
    )
)

coverage = targets.merge(
    accession_coverage,
    on="uniprot_canonical",
    how="left",
)

for column in [
    "n_position_relevant_pdb_entries",
    "n_position_relevant_chains",
]:
    coverage[column] = (
        coverage[column].fillna(0).astype(int)
    )

coverage["has_observed_variant_position"] = (
    coverage["n_position_relevant_pdb_entries"].gt(0)
)

coverage.to_csv(
    COVERAGE_OUTPUT,
    sep="\t",
    index=False,
    compression="gzip",
)

summary = {
    "uniprot_ready_variants": int(len(variants)),
    "variants_with_observed_pdb_position": int(
        len(covered_variants)
    ),
    "variants_without_observed_pdb_position": int(
        len(variants) - len(covered_variants)
    ),
    "structural_targets": int(len(targets)),
    "targets_with_observed_variant_position": int(
        coverage["has_observed_variant_position"].sum()
    ),
    "targets_without_observed_variant_position": int(
        (~coverage["has_observed_variant_position"]).sum()
    ),
    "retained_sifts_segments": int(len(retained)),
    "position_relevant_pdb_entries": int(len(pdb_manifest)),
}

with open(SUMMARY_OUTPUT, "w") as handle:
    json.dump(summary, handle, indent=2)

print()
print(json.dumps(summary, indent=2))
print()
print("Observed segments:", SEGMENT_OUTPUT)
print("Covered variants:", VARIANT_OUTPUT)
print("PDB manifest:", PDB_OUTPUT)
print("Coverage:", COVERAGE_OUTPUT)