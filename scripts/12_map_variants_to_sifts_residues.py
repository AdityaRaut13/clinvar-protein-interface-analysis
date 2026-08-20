#!/usr/bin/env python

import json
from pathlib import Path

import duckdb


ROOT = Path(os.environ["IF_PROJECT"])

VARIANT_FILE = (
    ROOT
    / "data/processed/uniprot_ready"
    / "clinvar_variants_uniprot_ready.parquet"
)
SEGMENT_FILE = (
    ROOT
    / "data/processed/structure_discovery"
    / "variant_observed_pdb_segments.parquet"
)

OUTPUT_DIR = ROOT / "data/processed/variant_structure_mapping"

MAPPING_OUTPUT = OUTPUT_DIR / "variant_sifts_residue_mapping.parquet"
INVALID_SEGMENT_OUTPUT = OUTPUT_DIR / "invalid_sifts_segments.parquet"
SUMMARY_OUTPUT = OUTPUT_DIR / "variant_sifts_mapping_summary.json"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

connection = duckdb.connect()

variant_path = str(VARIANT_FILE)
segment_path = str(SEGMENT_FILE)
mapping_path = str(MAPPING_OUTPUT)
invalid_path = str(INVALID_SEGMENT_OUTPUT)

connection.execute(f"""
CREATE TEMP TABLE segments AS
SELECT
    lower(CAST(pdb_id AS VARCHAR)) AS pdb_id,
    CAST(auth_chain_id AS VARCHAR) AS auth_chain_id,
    CAST(uniprot_canonical AS VARCHAR) AS uniprot_canonical,

    TRY_CAST(observed_residue_begin AS INTEGER)
        AS pdb_label_begin,

    TRY_CAST(observed_residue_end AS INTEGER)
        AS pdb_label_end,

    TRY_CAST(uniprot_segment_begin AS INTEGER)
        AS uniprot_begin,

    TRY_CAST(uniprot_segment_end AS INTEGER)
        AS uniprot_end

FROM read_parquet('{segment_path}')
""")

connection.execute("""
CREATE TEMP TABLE valid_segments AS
SELECT *
FROM segments
WHERE
    pdb_label_begin IS NOT NULL
    AND pdb_label_end IS NOT NULL
    AND uniprot_begin IS NOT NULL
    AND uniprot_end IS NOT NULL
    AND abs(pdb_label_end - pdb_label_begin)
        = abs(uniprot_end - uniprot_begin)
""")

connection.execute(f"""
COPY (
    SELECT *
    FROM segments
    WHERE
        pdb_label_begin IS NULL
        OR pdb_label_end IS NULL
        OR uniprot_begin IS NULL
        OR uniprot_end IS NULL
        OR abs(pdb_label_end - pdb_label_begin)
            != abs(uniprot_end - uniprot_begin)
)
TO '{invalid_path}'
(
    FORMAT PARQUET,
    COMPRESSION ZSTD
)
""")

connection.execute(f"""
COPY (
    SELECT DISTINCT
        CAST(v.variant_key AS VARCHAR) AS variant_key,
        CAST(v.variant_id AS VARCHAR) AS variant_id,
        CAST(v.gene_symbol AS VARCHAR) AS gene_symbol,
        CAST(v.protein_id AS VARCHAR) AS protein_id,
        CAST(v.source_sequence_key AS VARCHAR)
            AS source_sequence_key,

        CAST(v.uniprot_accession AS VARCHAR)
            AS selected_uniprot_accession,

        CAST(v.uniprot_canonical AS VARCHAR)
            AS uniprot_canonical,

        CAST(v.uniprot_position AS INTEGER)
            AS uniprot_position,

        CAST(v.aa_ref AS VARCHAR) AS aa_ref,
        CAST(v.aa_alt AS VARCHAR) AS aa_alt,
        CAST(v.label AS INTEGER) AS label,
        CAST(v.review_stars AS INTEGER) AS review_stars,

        s.pdb_id,
        s.auth_chain_id AS original_auth_asym_id,

        CAST(
            CASE
                WHEN s.uniprot_end = s.uniprot_begin
                THEN s.pdb_label_begin

                ELSE
                    s.pdb_label_begin
                    + (
                        CAST(v.uniprot_position AS INTEGER)
                        - s.uniprot_begin
                    )
                    * CASE
                        WHEN
                            (
                                s.pdb_label_end
                                - s.pdb_label_begin
                            )
                            * (
                                s.uniprot_end
                                - s.uniprot_begin
                            ) >= 0
                        THEN 1
                        ELSE -1
                    END
            END
            AS INTEGER
        ) AS pdb_label_seq_id,

        s.pdb_label_begin,
        s.pdb_label_end,
        s.uniprot_begin,
        s.uniprot_end

    FROM read_parquet('{variant_path}') AS v

    INNER JOIN valid_segments AS s
        ON CAST(v.uniprot_canonical AS VARCHAR)
            = s.uniprot_canonical

        AND CAST(v.uniprot_position AS INTEGER)
            BETWEEN
                least(s.uniprot_begin, s.uniprot_end)
            AND greatest(s.uniprot_begin, s.uniprot_end)
)
TO '{mapping_path}'
(
    FORMAT PARQUET,
    COMPRESSION ZSTD
)
""")

input_variants = connection.execute(f"""
SELECT count(DISTINCT variant_key)
FROM read_parquet('{variant_path}')
""").fetchone()[0]

mapped_rows = connection.execute(f"""
SELECT count(*)
FROM read_parquet('{mapping_path}')
""").fetchone()[0]

mapped_variants = connection.execute(f"""
SELECT count(DISTINCT variant_key)
FROM read_parquet('{mapping_path}')
""").fetchone()[0]

mapped_pdb_entries = connection.execute(f"""
SELECT count(DISTINCT pdb_id)
FROM read_parquet('{mapping_path}')
""").fetchone()[0]

mapped_chains = connection.execute(f"""
SELECT count(
    DISTINCT pdb_id || ':' || original_auth_asym_id
)
FROM read_parquet('{mapping_path}')
""").fetchone()[0]

valid_segments = connection.execute("""
SELECT count(*) FROM valid_segments
""").fetchone()[0]

invalid_segments = connection.execute(f"""
SELECT count(*)
FROM read_parquet('{invalid_path}')
""").fetchone()[0]

summary = {
    "input_uniprot_ready_variants": int(input_variants),
    "sifts_mapped_variants": int(mapped_variants),
    "sifts_unmapped_variants": int(
        input_variants - mapped_variants
    ),
    "variant_pdb_chain_mapping_rows": int(mapped_rows),
    "mapped_pdb_entries": int(mapped_pdb_entries),
    "mapped_pdb_chains": int(mapped_chains),
    "valid_observed_segments": int(valid_segments),
    "invalid_observed_segments": int(invalid_segments),
}

with open(SUMMARY_OUTPUT, "w") as handle:
    json.dump(summary, handle, indent=2)

print(json.dumps(summary, indent=2))
print()
print("Variant residue mapping:", MAPPING_OUTPUT)
print("Invalid segments:", INVALID_SEGMENT_OUTPUT)