#!/usr/bin/env python

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from tqdm import tqdm


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
STATUS_FILE = (
    ROOT
    / "data/processed/uniprot_mapping"
    / "protein_uniprot_final_status.tsv.gz"
)
UNIPROT_FASTA = (
    ROOT
    / "data/input/uniprot"
    / "human_uniprot_UP000005640_isoforms.fasta.gz"
)

OUTPUT_DIR = ROOT / "data/processed/uniprot_ready"

READY_OUTPUT = OUTPUT_DIR / "clinvar_variants_uniprot_ready.parquet"
FAILURE_OUTPUT = OUTPUT_DIR / "clinvar_variants_uniprot_failures.tsv.gz"
POSITION_OUTPUT = OUTPUT_DIR / "source_to_canonical_positions.parquet"
TARGET_OUTPUT = OUTPUT_DIR / "unique_uniprot_structural_targets.tsv.gz"
SUMMARY_OUTPUT = OUTPUT_DIR / "uniprot_ready_summary.json"

ACCEPTED = {
    "unique_exact_canonical",
    "global_unique_exact_canonical",
}

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def join_unique(values):
    return ";".join(sorted({
        str(value)
        for value in values.dropna()
        if str(value).strip()
    }))


def build_position_map(source_sequence, canonical_sequence):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    alignment = aligner.align(
        source_sequence,
        canonical_sequence,
    )[0]

    coordinates = alignment.coordinates
    position_map = {}

    for index in range(coordinates.shape[1] - 1):
        source_start = int(coordinates[0, index])
        source_end = int(coordinates[0, index + 1])
        canonical_start = int(coordinates[1, index])
        canonical_end = int(coordinates[1, index + 1])

        source_advance = source_end - source_start
        canonical_advance = canonical_end - canonical_start

        if (
            source_advance > 0
            and source_advance == canonical_advance
        ):
            for offset in range(source_advance):
                position_map[source_start + offset + 1] = (
                    canonical_start + offset + 1
                )

    return position_map


print("Loading input tables")

variants = pd.read_parquet(VARIANT_FILE)
proteins = pd.read_parquet(PROTEIN_FILE)

status = pd.read_csv(
    STATUS_FILE,
    sep="\t",
    dtype={"source_sequence_key": str},
    keep_default_na=False,
)

print("Loading UniProt sequences")

uniprot_sequences = {}

with gzip.open(UNIPROT_FASTA, "rt") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        parts = record.id.split("|")
        accession = parts[1] if len(parts) >= 3 else record.id
        uniprot_sequences[accession] = str(record.seq).upper()

protein_sequences = proteins.set_index(
    "source_sequence_key"
)["protein_sequence"].to_dict()

status_map = status.set_index("source_sequence_key")

accepted_status = status[
    status["mapping_status"].isin(ACCEPTED)
].copy()

positions_by_protein = (
    variants.groupby("source_sequence_key")
    ["protein_position"]
    .apply(lambda values: sorted(set(map(int, values))))
    .to_dict()
)

position_rows = []

print(
    "Mapping source positions to canonical UniProt coordinates:"
    f" {len(accepted_status):,} proteins"
)

for record in tqdm(
    accepted_status.itertuples(index=False),
    total=len(accepted_status),
):
    source_key = record.source_sequence_key
    selected = record.selected_uniprot_accession
    canonical = record.selected_uniprot_canonical

    source_sequence = protein_sequences[source_key]
    canonical_sequence = uniprot_sequences.get(canonical)
    positions = positions_by_protein.get(source_key, [])

    if canonical_sequence is None:
        for position in positions:
            position_rows.append({
                "source_sequence_key": source_key,
                "source_position": position,
                "selected_uniprot_accession": selected,
                "uniprot_canonical": canonical,
                "canonical_position": pd.NA,
                "coordinate_mapping_status":
                    "canonical_sequence_missing",
            })
        continue

    if selected == canonical:
        position_map = {
            position: position
            for position in positions
        }
        mapping_type = "exact_canonical_direct"

    else:
        position_map = build_position_map(
            source_sequence,
            canonical_sequence,
        )
        mapping_type = "isoform_to_canonical_alignment"

    for position in positions:
        canonical_position = position_map.get(position)

        if canonical_position is None:
            mapping_status = "position_unmapped_in_canonical"

        elif not (
            1 <= canonical_position <= len(canonical_sequence)
        ):
            canonical_position = None
            mapping_status = "canonical_position_out_of_range"

        elif (
            source_sequence[position - 1]
            != canonical_sequence[canonical_position - 1]
        ):
            canonical_position = None
            mapping_status = "canonical_reference_mismatch"

        else:
            mapping_status = mapping_type

        position_rows.append({
            "source_sequence_key": source_key,
            "source_position": position,
            "selected_uniprot_accession": selected,
            "uniprot_canonical": canonical,
            "canonical_position": canonical_position,
            "coordinate_mapping_status": mapping_status,
        })

position_mapping = pd.DataFrame(position_rows)

position_mapping["canonical_position"] = pd.to_numeric(
    position_mapping["canonical_position"],
    errors="coerce",
).astype("Int32")

position_mapping.to_parquet(
    POSITION_OUTPUT,
    index=False,
    compression="zstd",
)

mapping_columns = status[
    [
        "source_sequence_key",
        "selected_uniprot_accession",
        "selected_uniprot_canonical",
        "mapping_status",
        "mapping_source",
    ]
].rename(
    columns={
        "mapping_status": "protein_mapping_status",
        "mapping_source": "protein_mapping_source",
    }
)

result = variants.merge(
    mapping_columns,
    on="source_sequence_key",
    how="left",
)

result = result.merge(
    position_mapping[
        [
            "source_sequence_key",
            "source_position",
            "canonical_position",
            "coordinate_mapping_status",
        ]
    ],
    left_on=[
        "source_sequence_key",
        "protein_position",
    ],
    right_on=[
        "source_sequence_key",
        "source_position",
    ],
    how="left",
)

result["uniprot_accession"] = (
    result["selected_uniprot_accession"]
)
result["uniprot_canonical"] = (
    result["selected_uniprot_canonical"]
)
result["uniprot_position"] = (
    result["canonical_position"]
)

protein_mapping_accepted = result[
    "protein_mapping_status"
].isin(ACCEPTED)

coordinate_mapping_accepted = result[
    "coordinate_mapping_status"
].isin({
    "exact_canonical_direct",
    "isoform_to_canonical_alignment",
})

ready_mask = (
    protein_mapping_accepted
    & coordinate_mapping_accepted
    & result["uniprot_position"].notna()
)

ready = result[ready_mask].copy()
failures = result[~ready_mask].copy()

ready["uniprot_position"] = (
    ready["uniprot_position"].astype("int32")
)

ready = ready.drop(
    columns=[
        "source_position",
        "canonical_position",
    ],
    errors="ignore",
)

ready.to_parquet(
    READY_OUTPUT,
    index=False,
    compression="zstd",
)

failures.to_csv(
    FAILURE_OUTPUT,
    sep="\t",
    index=False,
    compression="gzip",
)

targets = (
    ready.groupby("uniprot_canonical", as_index=False)
    .agg(
        n_variants=("variant_key", "nunique"),
        n_source_proteins=("source_sequence_key", "nunique"),
        n_selected_accessions=("uniprot_accession", "nunique"),
        selected_accessions=("uniprot_accession", join_unique),
        gene_symbols=("gene_symbol", join_unique),
        n_benign=("label", lambda values: values.eq(0).sum()),
        n_pathogenic=("label", lambda values: values.eq(1).sum()),
    )
    .sort_values(
        ["n_variants", "uniprot_canonical"],
        ascending=[False, True],
    )
)

targets.to_csv(
    TARGET_OUTPUT,
    sep="\t",
    index=False,
    compression="gzip",
)

summary = {
    "input_variants": int(len(variants)),
    "uniprot_ready_variants": int(len(ready)),
    "uniprot_failed_variants": int(len(failures)),
    "unique_structural_targets": int(len(targets)),
    "direct_canonical_variants": int(
        ready["coordinate_mapping_status"]
        .eq("exact_canonical_direct")
        .sum()
    ),
    "isoform_lifted_variants": int(
        ready["coordinate_mapping_status"]
        .eq("isoform_to_canonical_alignment")
        .sum()
    ),
    "failure_status_counts": {
        str(key): int(value)
        for key, value in failures[
            "protein_mapping_status"
        ].fillna("missing").value_counts().items()
    },
    "coordinate_status_counts": {
        str(key): int(value)
        for key, value in result[
            "coordinate_mapping_status"
        ].fillna("not_attempted").value_counts().items()
    },
}

with open(SUMMARY_OUTPUT, "w") as handle:
    json.dump(summary, handle, indent=2)

print()
print(json.dumps(summary, indent=2))
print()
print("Ready variants:", READY_OUTPUT)
print("Failures:", FAILURE_OUTPUT)
print("Position mapping:", POSITION_OUTPUT)
print("Structural targets:", TARGET_OUTPUT)