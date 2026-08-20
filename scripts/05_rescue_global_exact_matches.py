#!/usr/bin/env python

import json
from pathlib import Path

import pandas as pd


ROOT = Path(os.environ["IF_PROJECT"])

STATUS_FILE = (
    ROOT
    / "data/processed/uniprot_mapping"
    / "protein_uniprot_exact_status.tsv.gz"
)
CATALOG_FILE = (
    ROOT
    / "data/processed/uniprot_mapping"
    / "human_uniprot_sequence_catalog.parquet"
)

OUTPUT_DIR = ROOT / "data/processed/uniprot_mapping"
FINAL_OUTPUT = OUTPUT_DIR / "protein_uniprot_final_status.tsv.gz"
SUMMARY_OUTPUT = OUTPUT_DIR / "uniprot_global_exact_summary.json"

RESCUE_STATUSES = {
    "no_annotated_candidate",
    "no_exact_candidate_family_match",
}


status = pd.read_csv(
    STATUS_FILE,
    sep="\t",
    dtype={
        "source_sequence_key": str,
        "sequence_hash": str,
    },
)

catalog = pd.read_parquet(CATALOG_FILE)

catalog_by_hash = {
    sequence_hash: group
    for sequence_hash, group in catalog.groupby(
        "sequence_hash",
        sort=False,
    )
}

status["mapping_source"] = "annotated_candidate_exact"
status.loc[
    ~status["mapping_status"].eq("unique_exact_canonical"),
    "mapping_source",
] = ""

rescued = 0
ambiguous = 0

for index, row in status.iterrows():
    if row["mapping_status"] not in RESCUE_STATUSES:
        continue

    matches = catalog_by_hash.get(row["sequence_hash"])

    if matches is None:
        continue

    accessions = sorted(
        matches["uniprot_accession"].unique()
    )
    canonicals = sorted(
        matches["uniprot_canonical"].unique()
    )

    status.at[index, "n_exact_accessions"] = len(accessions)
    status.at[index, "n_exact_canonical_accessions"] = (
        len(canonicals)
    )
    status.at[index, "exact_accessions"] = ";".join(accessions)
    status.at[index, "exact_canonical_accessions"] = (
        ";".join(canonicals)
    )

    if len(canonicals) == 1:
        selected_canonical = canonicals[0]

        selected_accession = (
            selected_canonical
            if selected_canonical in accessions
            else sorted(
                accessions,
                key=lambda accession: (
                    len(accession),
                    accession,
                ),
            )[0]
        )

        status.at[
            index, "selected_uniprot_accession"
        ] = selected_accession
        status.at[
            index, "selected_uniprot_canonical"
        ] = selected_canonical
        status.at[
            index, "mapping_status"
        ] = "global_unique_exact_canonical"
        status.at[
            index, "mapping_source"
        ] = "global_sequence_exact"

        rescued += 1

    elif len(canonicals) > 1:
        status.at[
            index, "mapping_status"
        ] = "global_ambiguous_multiple_canonicals"
        status.at[
            index, "mapping_source"
        ] = "global_sequence_exact_ambiguous"

        ambiguous += 1

status.to_csv(
    FINAL_OUTPUT,
    sep="\t",
    index=False,
    compression="gzip",
)

accepted_statuses = {
    "unique_exact_canonical",
    "global_unique_exact_canonical",
}

accepted = status["mapping_status"].isin(accepted_statuses)

summary = {
    "source_proteins": int(len(status)),
    "proteins_rescued_global_exact": int(rescued),
    "proteins_global_exact_ambiguous": int(ambiguous),
    "total_exact_mapped_proteins": int(accepted.sum()),
    "remaining_unmapped_proteins": int((~accepted).sum()),
    "total_exact_mapped_variants": int(
        status.loc[accepted, "n_variants"].sum()
    ),
    "remaining_unmapped_variants": int(
        status.loc[~accepted, "n_variants"].sum()
    ),
    "status_counts": {
        str(key): int(value)
        for key, value in status[
            "mapping_status"
        ].value_counts().items()
    },
}

with open(SUMMARY_OUTPUT, "w") as handle:
    json.dump(summary, handle, indent=2)

print(json.dumps(summary, indent=2))
print()
print("Final status:", FINAL_OUTPUT)