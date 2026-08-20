#!/usr/bin/env python3

import json
import os
from pathlib import Path

import duckdb


PROJECT = Path(os.environ["IF_PROJECT"])

OBSERVATIONS = (
    PROJECT / "data/processed/interface_mapping/"
    "variant_structure_interface_observations.parquet"
)

OUTDIR = PROJECT / "data/processed/interface_catalog"

PDB_EVIDENCE = OUTDIR / "protein_residue_pdb_evidence.parquet"
PARTNER_EVIDENCE = OUTDIR / "protein_residue_partner_evidence.parquet"
CATALOG = OUTDIR / "unique_protein_residue_interface_catalog.parquet"
MAPPED_VARIANTS = OUTDIR / "clinvar_variants_interface_mapped.parquet"
SUMMARY = OUTDIR / "interface_catalog_summary.json"

TEMP = Path(os.environ["IF_TEMP"])

OUTDIR.mkdir(parents=True, exist_ok=True)
TEMP.mkdir(parents=True, exist_ok=True)

for path in [
    PDB_EVIDENCE,
    PARTNER_EVIDENCE,
    CATALOG,
    MAPPED_VARIANTS,
]:
    path.unlink(missing_ok=True)

con = duckdb.connect()
con.execute("SET threads = 8")
con.execute("SET memory_limit = '24GB'")
con.execute(f"SET temp_directory = '{TEMP}'")
con.execute("SET preserve_insertion_order = false")
con.execute("PRAGMA enable_progress_bar")

print("Collapsing assembly copies within each PDB entry", flush=True)

con.execute(f"""
COPY (
    WITH unique_observations AS (
        SELECT DISTINCT
            uniprot_canonical,
            uniprot_position,
            aa_ref,
            pdb_id,
            assembly_id,
            target_label_asym_id,
            partner_label_asym_id,
            partner_original_auth_asym_id,
            partner_uniprot_accessions,
            pair_type,
            best_resolution,
            experimental_method,
            min_backbone_distance_to_partner,
            min_heavy_atom_distance_to_partner,
            residue_pinder_interface_10A,
            residue_heavy_interface_5A,
            residue_union_interface,

            CASE
                WHEN partner_uniprot_accessions IS NOT NULL
                 AND trim(partner_uniprot_accessions) <> ''
                THEN 'UNIPROT:' || partner_uniprot_accessions
                ELSE
                    'PDB:' || pdb_id || ':AUTH:' ||
                    coalesce(
                        partner_original_auth_asym_id,
                        partner_label_asym_id
                    )
            END AS partner_key

        FROM read_parquet('{OBSERVATIONS}')
        WHERE residue_match_status = 'reference_match'
    )

    SELECT
        uniprot_canonical,
        uniprot_position,
        any_value(aa_ref) AS aa_ref,
        pdb_id,
        partner_key,
        any_value(partner_uniprot_accessions)
            AS partner_uniprot_accessions,
        pair_type,

        count(*) AS n_assembly_chain_copies,
        count(DISTINCT assembly_id) AS n_assemblies,

        min(best_resolution) AS best_resolution,
        any_value(experimental_method) AS experimental_method,

        min(min_backbone_distance_to_partner)
            AS min_backbone_distance_to_partner,
        min(min_heavy_atom_distance_to_partner)
            AS min_heavy_atom_distance_to_partner,

        bool_or(
            coalesce(residue_pinder_interface_10A, false)
        ) AS pinder_interface_10A,

        bool_or(
            coalesce(residue_heavy_interface_5A, false)
        ) AS heavy_interface_5A,

        bool_or(
            coalesce(residue_union_interface, false)
        ) AS union_interface

    FROM unique_observations
    GROUP BY
        uniprot_canonical,
        uniprot_position,
        pdb_id,
        partner_key,
        pair_type
)
TO '{PDB_EVIDENCE}'
(
    FORMAT PARQUET,
    COMPRESSION ZSTD,
    ROW_GROUP_SIZE 100000
)
""")

print("Collapsing repeated PDB entries by interaction partner", flush=True)

con.execute(f"""
COPY (
    SELECT
        uniprot_canonical,
        uniprot_position,
        any_value(aa_ref) AS aa_ref,
        partner_key,
        any_value(partner_uniprot_accessions)
            AS partner_uniprot_accessions,
        pair_type,

        count(*) AS n_pdb_entries,

        min(min_backbone_distance_to_partner)
            AS min_backbone_distance_to_partner,
        min(min_heavy_atom_distance_to_partner)
            AS min_heavy_atom_distance_to_partner,

        count(*) FILTER (
            WHERE pinder_interface_10A
        ) AS n_pinder_pdb_entries,

        count(*) FILTER (
            WHERE heavy_interface_5A
        ) AS n_heavy_pdb_entries,

        count(*) FILTER (
            WHERE union_interface
        ) AS n_union_pdb_entries,

        avg(cast(pinder_interface_10A AS DOUBLE))
            AS pinder_pdb_fraction,

        avg(cast(heavy_interface_5A AS DOUBLE))
            AS heavy_pdb_fraction,

        avg(cast(union_interface AS DOUBLE))
            AS union_pdb_fraction,

        bool_or(pinder_interface_10A)
            AS ever_pinder_interface_10A,

        bool_or(heavy_interface_5A)
            AS ever_heavy_interface_5A,

        bool_or(union_interface)
            AS ever_union_interface

    FROM read_parquet('{PDB_EVIDENCE}')
    GROUP BY
        uniprot_canonical,
        uniprot_position,
        partner_key,
        pair_type
)
TO '{PARTNER_EVIDENCE}'
(
    FORMAT PARQUET,
    COMPRESSION ZSTD,
    ROW_GROUP_SIZE 100000
)
""")

print("Building unique protein-residue catalog", flush=True)

con.execute(f"""
COPY (
    WITH pdb_statistics AS (
        SELECT
            uniprot_canonical,
            uniprot_position,
            any_value(aa_ref) AS aa_ref,

            count(DISTINCT pdb_id) AS n_pdb_entries,
            count(*) AS n_pdb_partner_evidence,
            count(DISTINCT partner_key) AS n_partner_contexts,

            count(*) FILTER (
                WHERE pinder_interface_10A
            ) AS n_pinder_pdb_evidence,

            count(*) FILTER (
                WHERE heavy_interface_5A
            ) AS n_heavy_pdb_evidence,

            count(*) FILTER (
                WHERE union_interface
            ) AS n_union_pdb_evidence,

            avg(cast(pinder_interface_10A AS DOUBLE))
                AS pinder_pdb_evidence_fraction,

            avg(cast(heavy_interface_5A AS DOUBLE))
                AS heavy_pdb_evidence_fraction,

            avg(cast(union_interface AS DOUBLE))
                AS union_pdb_evidence_fraction,

            bool_or(pinder_interface_10A)
                AS ever_pinder_interface_10A,

            bool_or(heavy_interface_5A)
                AS ever_heavy_interface_5A,

            bool_or(union_interface)
                AS ever_union_interface,

            bool_or(
                pair_type = 'homomeric'
                AND union_interface
            ) AS ever_homomeric_union_interface,

            bool_or(
                pair_type = 'heteromeric'
                AND union_interface
            ) AS ever_heteromeric_union_interface

        FROM read_parquet('{PDB_EVIDENCE}')
        GROUP BY
            uniprot_canonical,
            uniprot_position
    ),

    partner_statistics AS (
        SELECT
            uniprot_canonical,
            uniprot_position,

            count(*) FILTER (
                WHERE ever_union_interface
            ) AS n_union_interface_partners,

            avg(cast(
                ever_pinder_interface_10A AS DOUBLE
            )) AS pinder_partner_fraction,

            avg(cast(
                ever_heavy_interface_5A AS DOUBLE
            )) AS heavy_partner_fraction,

            avg(cast(
                ever_union_interface AS DOUBLE
            )) AS union_partner_fraction

        FROM read_parquet('{PARTNER_EVIDENCE}')
        GROUP BY
            uniprot_canonical,
            uniprot_position
    )

    SELECT
        p.*,
        s.n_union_interface_partners,
        s.pinder_partner_fraction,
        s.heavy_partner_fraction,
        s.union_partner_fraction

    FROM pdb_statistics AS p
    LEFT JOIN partner_statistics AS s
      USING (uniprot_canonical, uniprot_position)
)
TO '{CATALOG}'
(
    FORMAT PARQUET,
    COMPRESSION ZSTD
)
""")

print("Mapping the catalog back to ClinVar variants", flush=True)

con.execute(f"""
COPY (
    WITH variants AS (
        SELECT DISTINCT
            variant_key,
            variant_id,
            gene_symbol,
            protein_id,
            selected_uniprot_accession,
            uniprot_canonical,
            uniprot_position,
            aa_ref,
            aa_alt,
            label,
            review_stars
        FROM read_parquet('{OBSERVATIONS}')
    )

    SELECT
        v.*,

        c.uniprot_canonical IS NOT NULL
            AS has_reference_structure_coverage,

        c.n_pdb_entries,
        c.n_pdb_partner_evidence,
        c.n_partner_contexts,
        c.n_union_interface_partners,

        c.ever_pinder_interface_10A,
        c.ever_heavy_interface_5A,
        c.ever_union_interface,

        c.ever_homomeric_union_interface,
        c.ever_heteromeric_union_interface,

        c.pinder_pdb_evidence_fraction,
        c.heavy_pdb_evidence_fraction,
        c.union_pdb_evidence_fraction,

        c.pinder_partner_fraction,
        c.heavy_partner_fraction,
        c.union_partner_fraction

    FROM variants AS v
    LEFT JOIN read_parquet('{CATALOG}') AS c
      ON v.uniprot_canonical = c.uniprot_canonical
     AND v.uniprot_position = c.uniprot_position
)
TO '{MAPPED_VARIANTS}'
(
    FORMAT PARQUET,
    COMPRESSION ZSTD
)
""")

counts = con.execute(f"""
SELECT
    count(*) AS mapped_variants,
    count(*) FILTER (
        WHERE has_reference_structure_coverage
    ) AS reference_covered_variants
FROM read_parquet('{MAPPED_VARIANTS}')
""").fetchone()

label_rows = con.execute(f"""
SELECT
    label,
    count(*) AS covered_variants,
    count(*) FILTER (
        WHERE ever_pinder_interface_10A
    ) AS pinder_variants,
    count(*) FILTER (
        WHERE ever_heavy_interface_5A
    ) AS heavy_variants,
    count(*) FILTER (
        WHERE ever_union_interface
    ) AS union_variants
FROM read_parquet('{MAPPED_VARIANTS}')
WHERE has_reference_structure_coverage
GROUP BY label
ORDER BY label
""").fetchall()

summary = {
    "protein_residue_pdb_evidence_rows": int(
        con.execute(
            f"SELECT count(*) FROM read_parquet('{PDB_EVIDENCE}')"
        ).fetchone()[0]
    ),
    "protein_residue_partner_evidence_rows": int(
        con.execute(
            f"SELECT count(*) FROM read_parquet('{PARTNER_EVIDENCE}')"
        ).fetchone()[0]
    ),
    "unique_protein_residues": int(
        con.execute(
            f"SELECT count(*) FROM read_parquet('{CATALOG}')"
        ).fetchone()[0]
    ),
    "variants_with_contacting_chain_pair": int(counts[0]),
    "variants_with_reference_structure_coverage": int(counts[1]),
    "interface_by_label": [
        {
            "label": int(label),
            "covered_variants": int(total),
            "pinder_interface_variants": int(pinder),
            "heavy_interface_variants": int(heavy),
            "union_interface_variants": int(union),
            "union_percent": round(100 * union / total, 2),
        }
        for label, total, pinder, heavy, union in label_rows
    ],
}

with open(SUMMARY, "w") as handle:
    json.dump(summary, handle, indent=2)

print(json.dumps(summary, indent=2))
print("\nPDB evidence:", PDB_EVIDENCE)
print("Partner evidence:", PARTNER_EVIDENCE)
print("Residue catalog:", CATALOG)
print("Mapped variants:", MAPPED_VARIANTS)

con.close()