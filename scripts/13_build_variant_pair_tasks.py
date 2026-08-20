#!/usr/bin/env python3

import json
import os
from pathlib import Path

import duckdb

PROJECT = Path(os.environ["IF_PROJECT"])

VARIANTS = (
    PROJECT / "data/processed/variant_structure_mapping/"
    "variant_sifts_residue_mapping.parquet"
)

PAIRS = (
    PROJECT / "data/processed/chain_pairs/"
    "contacting_chain_pairs_union10A5A/**/*.parquet"
)

OUTDIR = PROJECT / "data/processed/variant_pair_tasks"
TASKS = OUTDIR / "variant_chain_pair_tasks.parquet"
COVERAGE = OUTDIR / "variant_pair_coverage.parquet"
SUMMARY = OUTDIR / "variant_pair_task_summary.json"

TEMP = Path(os.environ["IF_TEMP"])

OUTDIR.mkdir(parents=True, exist_ok=True)
TEMP.mkdir(parents=True, exist_ok=True)

TASKS.unlink(missing_ok=True)
COVERAGE.unlink(missing_ok=True)

con = duckdb.connect()
con.execute("SET threads = 8")
con.execute("SET memory_limit = '24GB'")
con.execute(f"SET temp_directory = '{TEMP}'")
con.execute("SET preserve_insertion_order = false")

print("Building variant-chain-pair tasks")

con.execute(
    f"""
    COPY (
        SELECT DISTINCT
            v.variant_key,
            v.variant_id,
            v.gene_symbol,
            v.protein_id,
            v.selected_uniprot_accession,
            v.uniprot_canonical,
            v.uniprot_position,
            v.aa_ref,
            v.aa_alt,
            v.label,
            v.review_stars,

            lower(v.pdb_id) AS pdb_id,
            p.assembly_id,
            v.pdb_label_seq_id,

            p.target_label_asym_id,
            p.target_original_auth_asym_id,
            p.partner_label_asym_id,
            p.partner_original_auth_asym_id,
            p.partner_uniprot_accessions,
            p.pair_type,

            p.initial_min_backbone_distance,
            p.initial_min_heavy_atom_distance,
            p.pinder_pair_10A,
            p.heavy_contact_pair_5A,
            p.best_resolution,
            p.experimental_method

        FROM read_parquet('{VARIANTS}') AS v

        INNER JOIN read_parquet('{PAIRS}') AS p
          ON lower(v.pdb_id) = lower(p.pdb_id)
         AND v.uniprot_canonical = p.target_uniprot_canonical
         AND trim(v.original_auth_asym_id)
             = trim(p.target_original_auth_asym_id)
    )
    TO '{TASKS}'
    (
        FORMAT PARQUET,
        COMPRESSION ZSTD,
        ROW_GROUP_SIZE 100000
    )
    """
)

print("Building per-variant coverage table")

con.execute(
    f"""
    COPY (
        WITH mapped_variants AS (
            SELECT DISTINCT
                variant_key,
                variant_id,
                gene_symbol,
                protein_id,
                uniprot_canonical,
                uniprot_position,
                aa_ref,
                aa_alt,
                label,
                review_stars
            FROM read_parquet('{VARIANTS}')
        ),

        task_counts AS (
            SELECT
                variant_key,
                count(*) AS n_candidate_chain_pairs,
                count(
                    DISTINCT pdb_id || ':' || assembly_id
                ) AS n_candidate_assemblies,
                count(DISTINCT pdb_id) AS n_candidate_pdb_entries
            FROM read_parquet('{TASKS}')
            GROUP BY variant_key
        )

        SELECT
            v.*,
            coalesce(t.n_candidate_chain_pairs, 0)
                AS n_candidate_chain_pairs,
            coalesce(t.n_candidate_assemblies, 0)
                AS n_candidate_assemblies,
            coalesce(t.n_candidate_pdb_entries, 0)
                AS n_candidate_pdb_entries,
            t.variant_key IS NOT NULL
                AS has_contacting_chain_pair

        FROM mapped_variants AS v
        LEFT JOIN task_counts AS t USING (variant_key)
    )
    TO '{COVERAGE}'
    (
        FORMAT PARQUET,
        COMPRESSION ZSTD
    )
    """
)

def scalar(query):
    return int(con.execute(query).fetchone()[0])

label_rows = con.execute(
    f"""
    SELECT
        label,
        count(*) AS sifts_mapped_variants,
        count(*) FILTER (
            WHERE has_contacting_chain_pair
        ) AS variants_with_contacting_pair
    FROM read_parquet('{COVERAGE}')
    GROUP BY label
    ORDER BY label
    """
).fetchall()

summary = {
    "sifts_mapping_rows": scalar(
        f"SELECT count(*) FROM read_parquet('{VARIANTS}')"
    ),
    "sifts_mapped_variants": scalar(
        f"SELECT count(DISTINCT variant_key) FROM read_parquet('{VARIANTS}')"
    ),
    "retained_chain_pairs": scalar(
        f"SELECT count(*) FROM read_parquet('{PAIRS}')"
    ),
    "variant_chain_pair_tasks": scalar(
        f"SELECT count(*) FROM read_parquet('{TASKS}')"
    ),
    "variants_with_contacting_pair": scalar(
        f"""
        SELECT count(*)
        FROM read_parquet('{COVERAGE}')
        WHERE has_contacting_chain_pair
        """
    ),
    "unique_task_assemblies": scalar(
        f"""
        SELECT count(DISTINCT pdb_id || ':' || assembly_id)
        FROM read_parquet('{TASKS}')
        """
    ),
    "coverage_by_label": [
        {
            "label": int(label),
            "sifts_mapped_variants": int(total),
            "variants_with_contacting_pair": int(covered),
            "coverage_percent": round(100 * covered / total, 2),
        }
        for label, total, covered in label_rows
    ],
}

with open(SUMMARY, "w") as handle:
    json.dump(summary, handle, indent=2)

print(json.dumps(summary, indent=2))
print("\nTasks:", TASKS)
print("Coverage:", COVERAGE)
print("Summary:", SUMMARY)

con.close()