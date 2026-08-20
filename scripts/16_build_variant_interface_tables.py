#!/usr/bin/env python3

import json
import os
from pathlib import Path

import duckdb
import gemmi
import pyarrow as pa


PROJECT = Path(os.environ["IF_PROJECT"])

TASKS = (
    PROJECT / "data/processed/variant_pair_tasks/"
    "variant_chain_pair_tasks.parquet"
)

GEOMETRY = (
    PROJECT / "data/processed/variant_geometry/"
    "residue_partner_geometry_union10A5A/*.parquet"
)

OUTDIR = PROJECT / "data/processed/interface_mapping"

OBSERVATIONS = (
    OUTDIR / "variant_structure_interface_observations.parquet"
)

VARIANT_SUMMARY = (
    OUTDIR / "variant_interface_summary_all_structures.parquet"
)

RESIDUE_CODES = OUTDIR / "pdb_residue_code_mapping.tsv"
SUMMARY = OUTDIR / "variant_interface_mapping_summary.json"

TEMP = Path(os.environ["IF_TEMP"])

OUTDIR.mkdir(parents=True, exist_ok=True)
TEMP.mkdir(parents=True, exist_ok=True)

OBSERVATIONS.unlink(missing_ok=True)
VARIANT_SUMMARY.unlink(missing_ok=True)

con = duckdb.connect()
con.execute("SET threads = 8")
con.execute("SET memory_limit = '24GB'")
con.execute(f"SET temp_directory = '{TEMP}'")
con.execute("SET preserve_insertion_order = false")
con.execute("PRAGMA enable_progress_bar")

residue_names = [
    row[0]
    for row in con.execute(
        f"""
        SELECT DISTINCT pdb_residue_name
        FROM read_parquet('{GEOMETRY}')
        WHERE pdb_residue_name IS NOT NULL
        """
    ).fetchall()
]

code_rows = []

for name in residue_names:
    try:
        code = gemmi.find_tabulated_residue(
            name
        ).one_letter_code
        code = str(code).upper()

        if len(code) != 1 or code == "?":
            code = "X"
    except Exception:
        code = "X"

    code_rows.append({
        "pdb_residue_name": name,
        "pdb_residue_aa": code,
    })

code_table = pa.Table.from_pylist(code_rows)
con.register("residue_code_map", code_table)

con.execute(
    f"""
    COPY residue_code_map
    TO '{RESIDUE_CODES}'
    (
        FORMAT CSV,
        DELIMITER '\t',
        HEADER
    )
    """
)

print("Joining variant observations to residue geometry", flush=True)

con.execute(
    f"""
    COPY (
        SELECT
            t.*,

            g.pdb_residue_name,
            m.pdb_residue_aa,
            g.n_target_heavy_atoms,
            g.n_target_backbone_atoms,
            g.min_backbone_distance_to_partner,
            g.min_heavy_atom_distance_to_partner,
            g.residue_pinder_interface_10A,
            g.residue_heavy_interface_5A,

            CASE
                WHEN g.residue_pinder_interface_10A IS NULL
                 AND g.residue_heavy_interface_5A IS NULL
                    THEN NULL
                ELSE
                    coalesce(
                        g.residue_pinder_interface_10A,
                        false
                    )
                    OR
                    coalesce(
                        g.residue_heavy_interface_5A,
                        false
                    )
            END AS residue_union_interface,

            g.min_backbone_distance_to_partner IS NOT NULL
                AS pinder_geometry_available,

            g.min_heavy_atom_distance_to_partner IS NOT NULL
                AS heavy_geometry_available,

            CASE
                WHEN g.geometry_status IS NULL
                    THEN 'geometry_row_missing'
                WHEN g.pdb_residue_name IS NULL
                    THEN 'pdb_residue_unavailable'
                WHEN m.pdb_residue_aa IS NULL
                  OR m.pdb_residue_aa = 'X'
                    THEN 'unknown_pdb_residue'
                WHEN upper(t.aa_ref) = m.pdb_residue_aa
                    THEN 'reference_match'
                WHEN upper(t.aa_alt) = m.pdb_residue_aa
                    THEN 'alternate_match'
                ELSE 'other_mismatch'
            END AS residue_match_status,

            g.geometry_status

        FROM read_parquet('{TASKS}') AS t

        LEFT JOIN read_parquet('{GEOMETRY}') AS g
          ON t.pdb_id = g.pdb_id
         AND t.assembly_id = g.assembly_id
         AND t.target_label_asym_id
             = g.target_label_asym_id
         AND t.partner_label_asym_id
             = g.partner_label_asym_id
         AND t.pdb_label_seq_id
             = g.pdb_label_seq_id

        LEFT JOIN residue_code_map AS m
          ON g.pdb_residue_name = m.pdb_residue_name
    )
    TO '{OBSERVATIONS}'
    (
        FORMAT PARQUET,
        COMPRESSION ZSTD,
        ROW_GROUP_SIZE 100000
    )
    """
)

print("Collapsing observations to ClinVar variants", flush=True)

con.execute(
    f"""
    COPY (
        SELECT
            variant_key,
            any_value(variant_id) AS variant_id,
            any_value(gene_symbol) AS gene_symbol,
            any_value(protein_id) AS protein_id,
            any_value(selected_uniprot_accession)
                AS selected_uniprot_accession,
            any_value(uniprot_canonical)
                AS uniprot_canonical,
            any_value(uniprot_position)
                AS uniprot_position,
            any_value(aa_ref) AS aa_ref,
            any_value(aa_alt) AS aa_alt,
            any_value(label) AS label,
            any_value(review_stars) AS review_stars,

            count(*) AS n_structure_pair_observations,
            count(DISTINCT pdb_id) AS n_pdb_entries,

            count(
                DISTINCT pdb_id || ':' || assembly_id
            ) AS n_assemblies,

            count(DISTINCT
                pdb_id || ':' ||
                assembly_id || ':' ||
                target_label_asym_id || ':' ||
                partner_label_asym_id
            ) AS n_chain_pairs,

            count(*) FILTER (
                WHERE pinder_geometry_available
            ) AS n_pinder_observations,

            count(*) FILTER (
                WHERE heavy_geometry_available
            ) AS n_heavy_observations,

            count(*) FILTER (
                WHERE residue_pinder_interface_10A
            ) AS n_pinder_interface_observations,

            count(*) FILTER (
                WHERE residue_heavy_interface_5A
            ) AS n_heavy_interface_observations,

            count(*) FILTER (
                WHERE residue_union_interface
            ) AS n_union_interface_observations,

            count(*) FILTER (
                WHERE residue_match_status = 'reference_match'
            ) AS n_reference_match_observations,

            count(*) FILTER (
                WHERE residue_match_status = 'alternate_match'
            ) AS n_alternate_match_observations,

            count(*) FILTER (
                WHERE residue_match_status = 'other_mismatch'
            ) AS n_other_mismatch_observations,

            bool_or(
                coalesce(residue_pinder_interface_10A, false)
            ) AS ever_pinder_interface_10A,

            bool_or(
                coalesce(residue_heavy_interface_5A, false)
            ) AS ever_heavy_interface_5A,

            bool_or(
                coalesce(residue_union_interface, false)
            ) AS ever_union_interface,

            avg(
                CASE
                    WHEN pinder_geometry_available
                    THEN cast(
                        residue_pinder_interface_10A AS DOUBLE
                    )
                END
            ) AS pinder_interface_fraction,

            avg(
                CASE
                    WHEN heavy_geometry_available
                    THEN cast(
                        residue_heavy_interface_5A AS DOUBLE
                    )
                END
            ) AS heavy_interface_fraction,

            avg(
                CASE
                    WHEN pinder_geometry_available
                      OR heavy_geometry_available
                    THEN cast(
                        residue_union_interface AS DOUBLE
                    )
                END
            ) AS union_interface_fraction,

            bool_or(
                pair_type = 'homomeric'
                AND coalesce(residue_union_interface, false)
            ) AS ever_homomeric_union_interface,

            bool_or(
                pair_type = 'heteromeric'
                AND coalesce(residue_union_interface, false)
            ) AS ever_heteromeric_union_interface

        FROM read_parquet('{OBSERVATIONS}')
        GROUP BY variant_key
    )
    TO '{VARIANT_SUMMARY}'
    (
        FORMAT PARQUET,
        COMPRESSION ZSTD
    )
    """
)

observation_counts = dict(
    con.execute(
        f"""
        SELECT residue_match_status, count(*)
        FROM read_parquet('{OBSERVATIONS}')
        GROUP BY residue_match_status
        """
    ).fetchall()
)

variant_match_row = con.execute(
    f"""
    SELECT
        count(*) FILTER (
            WHERE n_reference_match_observations > 0
        ),
        count(*) FILTER (
            WHERE n_alternate_match_observations > 0
        ),
        count(*) FILTER (
            WHERE n_other_mismatch_observations > 0
        )
    FROM read_parquet('{VARIANT_SUMMARY}')
    """
).fetchone()

label_rows = con.execute(
    f"""
    SELECT
        label,
        count(*) AS variants,
        count(*) FILTER (
            WHERE ever_pinder_interface_10A
        ) AS pinder_interface_variants,
        count(*) FILTER (
            WHERE ever_heavy_interface_5A
        ) AS heavy_interface_variants,
        count(*) FILTER (
            WHERE ever_union_interface
        ) AS union_interface_variants
    FROM read_parquet('{VARIANT_SUMMARY}')
    GROUP BY label
    ORDER BY label
    """
).fetchall()

summary = {
    "variant_structure_observations": int(
        con.execute(
            f"SELECT count(*) FROM read_parquet('{OBSERVATIONS}')"
        ).fetchone()[0]
    ),
    "variants_with_contacting_chain_pair": int(
        con.execute(
            f"SELECT count(*) FROM read_parquet('{VARIANT_SUMMARY}')"
        ).fetchone()[0]
    ),
    "observation_residue_match_counts": {
        str(key): int(value)
        for key, value in observation_counts.items()
    },
    "variants_with_any_reference_match": int(variant_match_row[0]),
    "variants_with_any_alternate_match": int(variant_match_row[1]),
    "variants_with_any_other_mismatch": int(variant_match_row[2]),
    "interface_by_label": [
        {
            "label": int(label),
            "variants": int(total),
            "ever_pinder_interface_10A": int(pinder),
            "ever_heavy_interface_5A": int(heavy),
            "ever_union_interface": int(union),
            "union_percent": round(100 * union / total, 2),
        }
        for label, total, pinder, heavy, union in label_rows
    ],
}

with open(SUMMARY, "w") as handle:
    json.dump(summary, handle, indent=2)

print(json.dumps(summary, indent=2))
print("\nObservations:", OBSERVATIONS)
print("Variant summary:", VARIANT_SUMMARY)
print("Residue mapping:", RESIDUE_CODES)
print("Summary:", SUMMARY)

con.close()