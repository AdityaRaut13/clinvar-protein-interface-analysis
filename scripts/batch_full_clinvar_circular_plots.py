#!/usr/bin/env python3
"""Select informative proteins and create full-ClinVar circular interface plots."""

from __future__ import annotations

import argparse
import gzip
import os
from pathlib import Path
import subprocess
import sys

import duckdb
import pandas as pd


def fasta_lengths(path: Path) -> dict[str, int]:
    lengths = {}
    accession = None
    length = 0

    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if accession and "-" not in accession:
                    lengths[accession] = length
                token = line[1:].split()[0]
                accession = token.split("|")[1] if "|" in token else token
                length = 0
            else:
                length += len(line.strip())

    if accession and "-" not in accession:
        lengths[accession] = length
    return lengths


def select_proteins(args):
    con = duckdb.connect()
    definition = args.definition
    ever = {
        "union": "ever_union_interface",
        "pinder": "ever_pinder_interface_10A",
        "heavy": "ever_heavy_interface_5A",
    }[definition]

    frame = con.execute(f"""
    WITH variant_stats AS (
        SELECT
            uniprot_canonical,
            any_value(gene_symbol) AS gene_symbol,
            count(*) AS n_variants,
            count(*) FILTER (WHERE label = 0) AS n_benign,
            count(*) FILTER (WHERE label = 1) AS n_pathogenic,
            count(*) FILTER (WHERE {ever}) AS n_interface_variants
        FROM read_parquet('{args.variants}')
        WHERE has_reference_structure_coverage
        GROUP BY uniprot_canonical
    ),
    partner_stats AS (
        SELECT
            uniprot_canonical,
            count(DISTINCT partner_key) FILTER (
                WHERE {ever}
                  AND partner_uniprot_accessions IS NOT NULL
                  AND trim(partner_uniprot_accessions) <> ''
            ) AS n_mapped_interface_partners,
            count(DISTINCT uniprot_position) FILTER (
                WHERE {ever}
            ) AS n_interface_positions
        FROM read_parquet('{args.partner_evidence}')
        GROUP BY uniprot_canonical
    )
    SELECT v.*, p.n_mapped_interface_partners, p.n_interface_positions
    FROM variant_stats v
    JOIN partner_stats p USING (uniprot_canonical)
    """).fetchdf()
    con.close()

    lengths = fasta_lengths(args.fasta)
    frame["protein_length"] = frame["uniprot_canonical"].map(lengths)
    frame = frame.dropna(subset=["protein_length"]).copy()
    frame["protein_length"] = frame["protein_length"].astype(int)

    frame = frame.loc[
        (frame["protein_length"] <= args.max_length)
        & (frame["n_benign"] >= args.min_each_label)
        & (frame["n_pathogenic"] >= args.min_each_label)
        & (frame["n_mapped_interface_partners"] >= args.min_partners)
    ].copy()

    frame["selection_score"] = (
        (frame["n_benign"] * frame["n_pathogenic"]) ** 0.5
        * (1 + frame["n_mapped_interface_partners"]).map(lambda value: value ** 0.5)
    )
    return frame.sort_values(
        ["selection_score", "n_interface_positions"],
        ascending=False,
    ).head(args.top_n)


def main():
    args = parse_args()
    selected = select_proteins(args)
    if selected.empty:
        raise SystemExit("No proteins pass the requested selection filters")

    args.plot_dir.mkdir(parents=True, exist_ok=True)
    args.table_dir.mkdir(parents=True, exist_ok=True)
    selected.to_csv(args.table_dir / "selected_proteins.tsv", sep="\t", index=False)

    summaries = []
    for rank, row in enumerate(selected.itertuples(index=False), start=1):
        accession = row.uniprot_canonical
        stem = f"{rank:02d}_{accession}_{row.gene_symbol}_{row.protein_length}aa"
        output = args.plot_dir / f"{stem}_radial.png"
        support = args.table_dir / f"{stem}_support.tsv.gz"

        command = [
            sys.executable,
            str(args.plotter),
            "--protein", accession,
            "--length", str(row.protein_length),
            "--partner-evidence", str(args.partner_evidence),
            "--variants", str(args.variants),
            "--definition", args.definition,
            "--min-pdb-entries", str(args.min_pdb_entries),
            "--max-partners", str(args.max_partners),
            "--support-output", str(support),
            "--output", str(output),
        ]

        print(f"[{rank}/{len(selected)}] {accession} {row.gene_symbol}", flush=True)
        result = subprocess.run(command, text=True, capture_output=True)
        status = "ok" if result.returncode == 0 else "error"
        detail = (result.stdout if status == "ok" else result.stderr).strip()
        print(detail)
        summaries.append({
            "rank": rank,
            "uniprot_canonical": accession,
            "gene_symbol": row.gene_symbol,
            "protein_length": row.protein_length,
            "status": status,
            "plot": str(output) if status == "ok" else "",
            "detail": detail.replace("\n", " | "),
        })

    pd.DataFrame(summaries).to_csv(
        args.table_dir / "run_summary.tsv",
        sep="\t",
        index=False,
    )


def parse_args():
    project = Path(os.environ.get("IF_PROJECT", "."))
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--variants",
        type=Path,
        default=project / "data/processed/interface_catalog/clinvar_variants_interface_mapped.parquet",
    )
    parser.add_argument(
        "--partner-evidence",
        type=Path,
        default=project / "data/processed/interface_catalog/protein_residue_partner_evidence.parquet",
    )
    parser.add_argument(
        "--fasta",
        type=Path,
        default=project / "data/input/uniprot/human_uniprot_UP000005640_isoforms.fasta.gz",
    )
    parser.add_argument(
        "--plotter",
        type=Path,
        default=project / "scripts/full_clinvar_circular_interface_map.py",
    )
    parser.add_argument("--definition", choices=["union", "pinder", "heavy"], default="union")
    parser.add_argument("--top-n", type=int, default=10)
    parser.add_argument("--max-length", type=int, default=1000)
    parser.add_argument("--min-each-label", type=int, default=10)
    parser.add_argument("--min-partners", type=int, default=3)
    parser.add_argument("--min-pdb-entries", type=int, default=2)
    parser.add_argument("--max-partners", type=int, default=20)
    parser.add_argument(
        "--plot-dir",
        type=Path,
        default=project / "plots/full_clinvar_radial/top10",
    )
    parser.add_argument(
        "--table-dir",
        type=Path,
        default=project / "results/full_clinvar_radial/top10",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
