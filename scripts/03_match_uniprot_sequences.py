#!/usr/bin/env python

import gzip
import hashlib
import json
import os
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO


ROOT = Path(os.environ["IF_PROJECT"])

SOURCE_PROTEINS = (
    ROOT
    / "data/processed/source_validation"
    / "unique_source_proteins.parquet"
)
SOURCE_VARIANTS = (
    ROOT
    / "data/processed/source_validation"
    / "clinvar_source_valid.parquet"
)
INITIAL_CANDIDATES = (
    ROOT
    / "data/processed/uniprot_mapping"
    / "protein_uniprot_candidates.parquet"
)
COMBINED_CANDIDATES = (
    ROOT
    / "data/processed/uniprot_mapping"
    / "protein_uniprot_candidates_combined.parquet"
)
CANDIDATE_FILE = (
    COMBINED_CANDIDATES
    if COMBINED_CANDIDATES.exists()
    else INITIAL_CANDIDATES
)
UNIPROT_FASTA = (
    ROOT
    / "data/input/uniprot"
    / "human_uniprot_UP000005640_isoforms.fasta.gz"
)

OUTPUT_DIR = ROOT / "data/processed/uniprot_mapping"

CATALOG_OUTPUT = OUTPUT_DIR / "human_uniprot_sequence_catalog.parquet"
MATCH_OUTPUT = OUTPUT_DIR / "protein_uniprot_exact_matches.parquet"
STATUS_OUTPUT = OUTPUT_DIR / "protein_uniprot_exact_status.tsv.gz"
SUMMARY_OUTPUT = OUTPUT_DIR / "uniprot_exact_match_summary.json"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def sha256(sequence):
    return hashlib.sha256(sequence.encode()).hexdigest()


def canonical_accession(accession):
    return accession.split("-")[0]


def lookup_accession(accession):
    # UniProt canonical isoform "-1" is normally stored without "-1".
    if accession.endswith("-1"):
        return canonical_accession(accession)
    return accession


print("Reading source proteins")

source_proteins = pd.read_parquet(SOURCE_PROTEINS)
variants = pd.read_parquet(
    SOURCE_VARIANTS,
    columns=["source_sequence_key", "variant_key"],
)
candidates = pd.read_parquet(CANDIDATE_FILE)

print("Reading UniProt FASTA")

catalog_rows = []
family_hash_index = defaultdict(list)

with gzip.open(UNIPROT_FASTA, "rt") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        parts = record.id.split("|")
        accession = parts[1] if len(parts) >= 3 else record.id
        accession = accession.strip()

        sequence = str(record.seq).upper()
        sequence_hash = sha256(sequence)
        canonical = canonical_accession(accession)

        row = {
            "uniprot_accession": accession,
            "uniprot_canonical": canonical,
            "is_isoform": "-" in accession,
            "sequence_length": len(sequence),
            "sequence_hash": sequence_hash,
        }

        catalog_rows.append(row)
        family_hash_index[(canonical, sequence_hash)].append(row)

catalog = pd.DataFrame(catalog_rows)

if catalog["uniprot_accession"].duplicated().any():
    raise ValueError("Duplicate accessions found in UniProt FASTA")

catalog.to_parquet(
    CATALOG_OUTPUT,
    index=False,
    compression="zstd",
)

print(f"UniProt sequences indexed: {len(catalog):,}")

candidate_groups = {
    key: group
    for key, group in candidates.groupby(
        "source_sequence_key",
        sort=False,
    )
}

variant_counts = (
    variants.groupby("source_sequence_key")
    ["variant_key"]
    .nunique()
    .to_dict()
)

match_rows = []
status_rows = []

for protein in source_proteins.itertuples(index=False):
    source_key = str(protein.source_sequence_key)
    source_hash = str(protein.sequence_hash)
    candidate_group = candidate_groups.get(source_key)

    if candidate_group is None:
        status_rows.append({
            "source_sequence_key": source_key,
            "sequence_hash": source_hash,
            "sequence_length": int(protein.sequence_length),
            "n_variants": int(variant_counts.get(source_key, 0)),
            "n_candidates": 0,
            "n_exact_accessions": 0,
            "n_exact_canonical_accessions": 0,
            "exact_accessions": "",
            "exact_canonical_accessions": "",
            "selected_uniprot_accession": "",
            "selected_uniprot_canonical": "",
            "mapping_status": "no_annotated_candidate",
        })
        continue

    annotated_priority = {}

    for candidate in candidate_group.itertuples(index=False):
        accession = str(candidate.uniprot_accession)
        lookup = lookup_accession(accession)
        priority = int(candidate.annotation_priority)

        current = annotated_priority.get(lookup)

        if current is None or priority < current:
            annotated_priority[lookup] = priority

    exact_by_accession = {}

    for canonical in sorted(
        candidate_group["uniprot_canonical"].unique()
    ):
        for match in family_hash_index.get(
            (canonical, source_hash),
            [],
        ):
            exact_by_accession[
                match["uniprot_accession"]
            ] = match

    exact_accessions = sorted(exact_by_accession)
    exact_canonicals = sorted({
        exact_by_accession[accession]["uniprot_canonical"]
        for accession in exact_accessions
    })

    for accession in exact_accessions:
        record = exact_by_accession[accession]

        match_rows.append({
            "source_sequence_key": source_key,
            "source_sequence_hash": source_hash,
            "uniprot_accession": accession,
            "uniprot_canonical": record["uniprot_canonical"],
            "uniprot_sequence_length": record["sequence_length"],
            "is_isoform": record["is_isoform"],
            "was_directly_annotated":
                accession in annotated_priority,
            "annotation_priority":
                annotated_priority.get(accession),
            "match_type": (
                "exact_annotated_accession"
                if accession in annotated_priority
                else "exact_candidate_family_isoform"
            ),
        })

    selected_accession = ""
    selected_canonical = ""

    if not exact_accessions:
        mapping_status = "no_exact_candidate_family_match"

    elif len(exact_canonicals) > 1:
        mapping_status = "ambiguous_multiple_exact_canonicals"

    else:
        selected_canonical = exact_canonicals[0]

        def rank(accession):
            directly_annotated = accession in annotated_priority
            annotation_priority = annotated_priority.get(
                accession,
                999,
            )
            is_canonical = accession == selected_canonical

            return (
                0 if directly_annotated else 1,
                annotation_priority,
                0 if is_canonical else 1,
                accession,
            )

        selected_accession = min(exact_accessions, key=rank)
        mapping_status = "unique_exact_canonical"

    status_rows.append({
        "source_sequence_key": source_key,
        "sequence_hash": source_hash,
        "sequence_length": int(protein.sequence_length),
        "n_variants": int(variant_counts.get(source_key, 0)),
        "n_candidates": int(
            candidate_group["uniprot_accession"].nunique()
        ),
        "n_exact_accessions": len(exact_accessions),
        "n_exact_canonical_accessions": len(exact_canonicals),
        "exact_accessions": ";".join(exact_accessions),
        "exact_canonical_accessions":
            ";".join(exact_canonicals),
        "selected_uniprot_accession": selected_accession,
        "selected_uniprot_canonical": selected_canonical,
        "mapping_status": mapping_status,
    })

matches = pd.DataFrame(match_rows)
status = pd.DataFrame(status_rows)

matches.to_parquet(
    MATCH_OUTPUT,
    index=False,
    compression="zstd",
)

status.to_csv(
    STATUS_OUTPUT,
    sep="\t",
    index=False,
    compression="gzip",
)

selected_keys = set(
    status.loc[
        status["mapping_status"].eq("unique_exact_canonical"),
        "source_sequence_key",
    ]
)

variant_exact_coverage = variants[
    "source_sequence_key"
].isin(selected_keys)

summary = {
    "source_proteins": int(len(source_proteins)),
    "uniprot_sequences": int(len(catalog)),
    "proteins_unique_exact": int(
        status["mapping_status"]
        .eq("unique_exact_canonical")
        .sum()
    ),
    "proteins_ambiguous_exact": int(
        status["mapping_status"]
        .eq("ambiguous_multiple_exact_canonicals")
        .sum()
    ),
    "proteins_no_exact_match": int(
        status["mapping_status"]
        .eq("no_exact_candidate_family_match")
        .sum()
    ),
    "proteins_without_candidate": int(
        status["mapping_status"]
        .eq("no_annotated_candidate")
        .sum()
    ),
    "variants_unique_exact": int(
        variant_exact_coverage.sum()
    ),
    "variants_not_unique_exact": int(
        (~variant_exact_coverage).sum()
    ),
}

with open(SUMMARY_OUTPUT, "w") as handle:
    json.dump(summary, handle, indent=2)

print()
print(json.dumps(summary, indent=2))
print()
print("Catalog:", CATALOG_OUTPUT)
print("Matches:", MATCH_OUTPUT)
print("Status:", STATUS_OUTPUT)
print("Summary:", SUMMARY_OUTPUT)
