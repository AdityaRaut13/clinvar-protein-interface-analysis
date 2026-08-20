#!/usr/bin/env python3
"""Circular partner-interface map for the full ClinVar interface catalog."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import re

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import numpy as np
import pandas as pd


COLORS = {
    "homomeric": "#315F91",
    "heteromeric": "#CF873D",
    "pathogenic": "#B44335",
    "benign": "#4C78A8",
}

DEFINITIONS = {
    "union": {
        "count": "n_union_pdb_entries",
        "fraction": "union_pdb_fraction",
        "ever": "ever_union_interface",
        "title": "10 A backbone OR 5 A heavy atom",
    },
    "pinder": {
        "count": "n_pinder_pdb_entries",
        "fraction": "pinder_pdb_fraction",
        "ever": "ever_pinder_interface_10A",
        "title": "PINDER backbone <=10 A",
    },
    "heavy": {
        "count": "n_heavy_pdb_entries",
        "fraction": "heavy_pdb_fraction",
        "ever": "ever_heavy_interface_5A",
        "title": "heavy atom <=5 A",
    },
}


def angle(position, length: int, gap_degrees: float):
    usable = 2 * np.pi - np.deg2rad(gap_degrees)
    return np.deg2rad(gap_degrees / 2) + (np.asarray(position) - 1) * usable / length


def load_inputs(args):
    evidence = pd.read_parquet(
        args.partner_evidence,
        filters=[("uniprot_canonical", "==", args.protein)],
    )
    variants = pd.read_parquet(
        args.variants,
        filters=[("uniprot_canonical", "==", args.protein)],
    )

    if evidence.empty:
        raise SystemExit(f"No partner evidence for {args.protein}")
    if variants.empty:
        raise SystemExit(f"No ClinVar variants for {args.protein}")

    variants = variants.loc[variants["has_reference_structure_coverage"]].copy()
    return evidence, variants


def prepare_support(evidence, args):
    definition = DEFINITIONS[args.definition]

    if not args.include_unmapped_partners:
        evidence = evidence.loc[
            evidence["partner_uniprot_accessions"].notna()
            & evidence["partner_uniprot_accessions"].astype(str).str.strip().ne("")
        ].copy()

    evidence = evidence.loc[
        (evidence["n_pdb_entries"] >= args.min_pdb_entries)
        & (evidence[definition["count"]] > 0)
    ].copy()

    if evidence.empty:
        raise SystemExit("No partners pass the requested evidence filters")

    evidence["track"] = evidence["partner_key"].astype(str)
    evidence["partner_label"] = evidence["partner_uniprot_accessions"].fillna(
        evidence["partner_key"]
    )
    evidence["group"] = evidence["partner_uniprot_accessions"].apply(
        lambda value: (
            "homomeric"
            if {
                item.split("-")[0]
                for item in re.split(r"[;,]", str(value))
                if item and item.lower() != "nan"
            } == {args.protein}
            else "heteromeric"
        )
    )
    evidence["value"] = evidence[definition["fraction"]].clip(0, 1)

    # Merge duplicate coordinate/entity classifications for the same mapped
    # biological partner and residue. Maxima are conservative with respect to
    # overlapping PDB evidence and avoid counting one partner as two rings.
    evidence = (
        evidence.groupby(
            ["track", "partner_label", "group", "uniprot_position"],
            as_index=False,
        )
        .agg(
            n_pdb_entries=("n_pdb_entries", "max"),
            direct_count=(definition["count"], "max"),
            value=("value", "max"),
        )
    )

    ranking = (
        evidence.groupby(["track", "partner_label", "group"], as_index=False)
        .agg(
            supported_positions=("uniprot_position", "nunique"),
            direct_pdb_evidence=("direct_count", "sum"),
            observed_pdb_evidence=("n_pdb_entries", "sum"),
            maximum_pdb_support=("n_pdb_entries", "max"),
        )
        .sort_values(
            ["direct_pdb_evidence", "supported_positions", "track"],
            ascending=[False, False, True],
        )
        .head(args.max_partners)
    )

    support = evidence.loc[evidence["track"].isin(ranking["track"])].copy()
    support = support.merge(
        ranking,
        on=["track", "partner_label", "group"],
        how="left",
    )
    return support, ranking


def prepare_variant_overlay(variants, definition):
    direct_column = DEFINITIONS[definition]["ever"]
    variants["label_name"] = variants["label"].map(
        {0: "benign", 1: "pathogenic"}
    )
    variants["interface_class"] = np.where(
        variants[direct_column].fillna(False),
        "direct",
        "observed_noninterface",
    )
    return (
        variants.groupby(
            ["uniprot_position", "label_name", "interface_class"],
            as_index=False,
        )
        .agg(n_variants=("variant_key", "nunique"))
        .rename(columns={"uniprot_position": "position"})
    )


def plot(args):
    evidence, variants = load_inputs(args)
    support, ranking = prepare_support(evidence, args)
    overlay = prepare_variant_overlay(variants, args.definition)

    inferred_length = int(max(
        support["uniprot_position"].max(),
        variants["uniprot_position"].max(),
    ))
    length = args.length or inferred_length

    support = support.loc[support["uniprot_position"].between(1, length)].copy()
    overlay = overlay.loc[overlay["position"].between(1, length)].copy()

    ranking = ranking.loc[ranking["track"].isin(support["track"])].copy()
    ranking["group_order"] = ranking["group"].map(
        {"homomeric": 0, "heteromeric": 1}
    ).fillna(2)
    ranking = ranking.sort_values(
        ["group_order", "direct_pdb_evidence"],
        ascending=[True, False],
    ).reset_index(drop=True)

    tracks = ranking["track"].tolist()
    track_index = {track: index for index, track in enumerate(tracks)}
    group_by_track = dict(zip(ranking["track"], ranking["group"]))

    matrix = np.zeros((len(tracks), length), dtype=np.float32)
    for row in support.itertuples(index=False):
        index = track_index[row.track]
        position = int(row.uniprot_position) - 1
        matrix[index, position] = max(matrix[index, position], float(row.value))

    figure, axis = plt.subplots(
        figsize=(14, 10),
        subplot_kw={"projection": "polar"},
    )
    figure.subplots_adjust(left=0.03, right=0.72, top=0.97, bottom=0.03)
    axis.set_theta_zero_location("N")
    axis.set_theta_direction(-1)

    inner, outer = 2.6, 8.2
    radial_edges = np.linspace(inner, outer, len(tracks) + 1)
    theta_edges = angle(
        np.arange(1, length + 2) - 0.5,
        length,
        args.gap,
    )

    axis.pcolormesh(
        theta_edges,
        radial_edges,
        np.ones_like(matrix),
        cmap=ListedColormap(["#F4F6F3"]),
        shading="flat",
        edgecolors="#DDE2DC",
        linewidth=0.025,
        rasterized=True,
    )

    for group in ["homomeric", "heteromeric"]:
        rows = np.asarray([
            group_by_track[track] == group for track in tracks
        ])
        visible = np.ma.masked_where(
            ~rows[:, None] | (matrix <= 0),
            matrix,
        )
        axis.pcolormesh(
            theta_edges,
            radial_edges,
            visible,
            cmap=LinearSegmentedColormap.from_list(
                group,
                ["#F4F6F3", COLORS[group]],
            ),
            vmin=0,
            vmax=1,
            shading="flat",
            edgecolors="none",
            rasterized=True,
        )

    shell_radius = {
        "direct": outer + 0.35,
        "observed_noninterface": outer + 0.82,
    }
    marker = {"direct": "o", "observed_noninterface": "X"}

    for (label, interface_class), frame in overlay.groupby(
        ["label_name", "interface_class"],
        sort=False,
    ):
        jitter = 0.09 if label == "pathogenic" else -0.09
        theta = angle(frame["position"].to_numpy() + jitter, length, args.gap)
        sizes = 18 + 12 * np.sqrt(frame["n_variants"].to_numpy())
        axis.scatter(
            theta,
            np.full(len(frame), shell_radius[interface_class]),
            s=sizes,
            marker=marker[interface_class],
            c=COLORS[label],
            edgecolors="#202020",
            linewidths=0.35,
            zorder=8,
        )

    gene = variants["gene_symbol"].dropna().astype(str).head(1)
    center = f"{gene.iloc[0]}\n{args.protein}" if len(gene) else args.protein
    axis.text(0, 0, center, ha="center", va="center", fontsize=17)
    axis.text(
        0.5,
        0.46,
        DEFINITIONS[args.definition]["title"],
        transform=axis.transAxes,
        ha="center",
        va="center",
        fontsize=9,
        color="#555555",
    )

    axis.text(
        angle(1, length, args.gap), inner - 0.25, "1",
        ha="right", va="center", fontsize=9,
    )
    axis.text(
        angle(length + 1, length, args.gap),
        inner - 0.25,
        str(length),
        ha="left",
        va="center",
        fontsize=9,
    )

    partner_handles = [
        plt.Line2D([], [], marker="o", ls="", color=COLORS[group], label=group)
        for group in ["homomeric", "heteromeric"]
        if group in set(ranking["group"])
    ]
    first_legend = figure.legend(
        handles=partner_handles,
        title="Partner type\nIntensity: PDB support fraction",
        loc="upper left",
        bbox_to_anchor=(0.75, 0.97),
        frameon=False,
    )

    variant_handles = [
        plt.Line2D([], [], marker="o", ls="", color=COLORS["pathogenic"],
                   markeredgecolor="#202020", label="Pathogenic"),
        plt.Line2D([], [], marker="o", ls="", color=COLORS["benign"],
                   markeredgecolor="#202020", label="Benign"),
        plt.Line2D([], [], marker="o", ls="", color="#666666", label="Direct"),
        plt.Line2D([], [], marker="X", ls="", color="#666666",
                   label="Observed non-interface"),
    ]
    figure.legend(
        handles=variant_handles,
        title="ClinVar overlay",
        loc="lower left",
        bbox_to_anchor=(0.75, 0.10),
        frameon=False,
    )

    partner_lines = ["Partner rings (inner -> outer)"]
    for index, row in ranking.iterrows():
        partner_lines.append(
            f"{index + 1}. {row.partner_label} "
            f"(positions={int(row.supported_positions)}, "
            f"evidence={int(row.direct_pdb_evidence)})"
        )
    ring_font = 8 if len(ranking) <= 15 else 7
    figure.text(
        0.75, 0.72, "\n".join(partner_lines),
        va="top", fontsize=ring_font,
    )

    axis.set_ylim(0, outer + 1.35)
    axis.set_axis_off()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(
        args.output,
        dpi=args.dpi,
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(figure)

    if args.support_output:
        args.support_output.parent.mkdir(parents=True, exist_ok=True)
        support.to_csv(
            args.support_output,
            sep="\t",
            index=False,
            compression="gzip" if args.support_output.name.endswith(".gz") else None,
        )

    print("Protein:", args.protein)
    print("Length:", length)
    print("Partner rings:", len(ranking))
    print("ClinVar variants:", variants["variant_key"].nunique())
    print("Saved:", args.output)


def parse_args():
    project = Path(os.environ.get("IF_PROJECT", "."))
    parser = argparse.ArgumentParser()
    parser.add_argument("--protein", required=True)
    parser.add_argument("--length", type=int)
    parser.add_argument(
        "--partner-evidence",
        type=Path,
        default=project / "data/processed/interface_catalog/protein_residue_partner_evidence.parquet",
    )
    parser.add_argument(
        "--variants",
        type=Path,
        default=project / "data/processed/interface_catalog/clinvar_variants_interface_mapped.parquet",
    )
    parser.add_argument("--definition", choices=DEFINITIONS, default="union")
    parser.add_argument("--min-pdb-entries", type=int, default=2)
    parser.add_argument("--max-partners", type=int, default=20)
    parser.add_argument("--include-unmapped-partners", action="store_true")
    parser.add_argument("--support-output", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--gap", type=float, default=8.0)
    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args()


if __name__ == "__main__":
    plot(parse_args())
