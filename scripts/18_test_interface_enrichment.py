#!/usr/bin/env python3

import json
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.discrete.conditional_models import ConditionalLogit
from statsmodels.stats.contingency_tables import StratifiedTable


PROJECT = Path(os.environ["IF_PROJECT"])

INPUT = (
    PROJECT / "data/processed/interface_catalog/"
    "clinvar_variants_interface_mapped.parquet"
)

OUTDIR = PROJECT / "results/interface_enrichment"
TABLE = OUTDIR / "interface_enrichment_results.tsv"
SUMMARY = OUTDIR / "interface_enrichment_summary.json"

OUTDIR.mkdir(parents=True, exist_ok=True)

DEFINITIONS = {
    "pinder_10A": "ever_pinder_interface_10A",
    "heavy_5A": "ever_heavy_interface_5A",
    "union_10A_or_5A": "ever_union_interface",
}


def raw_test(frame, interface_column):
    interface = frame[interface_column].astype(bool)
    label = frame["label"].astype(int)

    a = int(((label == 1) & interface).sum())
    b = int(((label == 1) & ~interface).sum())
    c = int(((label == 0) & interface).sum())
    d = int(((label == 0) & ~interface).sum())

    odds_ratio, p_value = fisher_exact(
        [[a, b], [c, d]],
        alternative="two-sided",
    )

    standard_error = math.sqrt(
        1 / a + 1 / b + 1 / c + 1 / d
    )

    lower = math.exp(math.log(odds_ratio) - 1.96 * standard_error)
    upper = math.exp(math.log(odds_ratio) + 1.96 * standard_error)

    return {
        "pathogenic_interface": a,
        "pathogenic_noninterface": b,
        "benign_interface": c,
        "benign_noninterface": d,
        "pathogenic_interface_percent":
            round(100 * a / (a + b), 4),
        "benign_interface_percent":
            round(100 * c / (c + d), 4),
        "odds_ratio": float(odds_ratio),
        "ci_lower": float(lower),
        "ci_upper": float(upper),
        "p_value": float(p_value),
    }


def protein_stratified_test(frame, interface_column):
    tables = []

    for _, group in frame.groupby(
        "uniprot_canonical",
        sort=False,
    ):
        interface = group[interface_column].astype(bool)
        label = group["label"].astype(int)

        table = np.array([
            [
                ((label == 1) & interface).sum(),
                ((label == 1) & ~interface).sum(),
            ],
            [
                ((label == 0) & interface).sum(),
                ((label == 0) & ~interface).sum(),
            ],
        ], dtype=float)

        if (
            np.all(table.sum(axis=0) > 0)
            and np.all(table.sum(axis=1) > 0)
        ):
            tables.append(table)

    combined = np.stack(tables, axis=2)
    result = StratifiedTable(combined)
    lower, upper = result.oddsratio_pooled_confint()

    return {
        "contributing_proteins": len(tables),
        "cmh_odds_ratio": float(result.oddsratio_pooled),
        "cmh_ci_lower": float(lower),
        "cmh_ci_upper": float(upper),
        "cmh_p_value": float(
            result.test_null_odds().pvalue
        ),
    }


df = pd.read_parquet(INPUT)
df = df.loc[df["has_reference_structure_coverage"]].copy()

variant_results = []

for definition, column in DEFINITIONS.items():
    result = raw_test(df, column)
    result.update(
        protein_stratified_test(df, column)
    )
    result.update({
        "analysis_level": "variant",
        "definition": definition,
        "n_observations": len(df),
    })
    variant_results.append(result)

position_df = (
    df.groupby(
        ["uniprot_canonical", "uniprot_position"],
        as_index=False,
    )
    .agg(
        n_benign=("label", lambda x: int((x == 0).sum())),
        n_pathogenic=("label", lambda x: int((x == 1).sum())),
        ever_pinder_interface_10A=(
            "ever_pinder_interface_10A", "max"
        ),
        ever_heavy_interface_5A=(
            "ever_heavy_interface_5A", "max"
        ),
        ever_union_interface=(
            "ever_union_interface", "max"
        ),
    )
)

position_df["label"] = np.select(
    [
        (position_df["n_pathogenic"] > 0)
        & (position_df["n_benign"] == 0),

        (position_df["n_benign"] > 0)
        & (position_df["n_pathogenic"] == 0),
    ],
    [1, 0],
    default=-1,
)

mixed_positions = int((position_df["label"] == -1).sum())
pure_positions = position_df.loc[
    position_df["label"].isin([0, 1])
].copy()

position_results = []

for definition, column in DEFINITIONS.items():
    result = raw_test(pure_positions, column)
    result.update(
        protein_stratified_test(pure_positions, column)
    )
    result.update({
        "analysis_level": "unique_residue",
        "definition": definition,
        "n_observations": len(pure_positions),
    })
    position_results.append(result)

model_df = df[
    [
        "label",
        "uniprot_canonical",
        "ever_union_interface",
        "n_pdb_entries",
        "n_partner_contexts",
        "review_stars",
    ]
].dropna().copy()

valid_proteins = (
    model_df.groupby("uniprot_canonical")["label"]
    .nunique()
)
valid_proteins = valid_proteins[valid_proteins == 2].index

model_df = model_df.loc[
    model_df["uniprot_canonical"].isin(valid_proteins)
].copy()

model_df["interface"] = (
    model_df["ever_union_interface"].astype(float)
)
model_df["log_pdb"] = np.log1p(model_df["n_pdb_entries"])
model_df["log_partners"] = np.log1p(
    model_df["n_partner_contexts"]
)

for column in ["log_pdb", "log_partners", "review_stars"]:
    standard_deviation = model_df[column].std()

    if standard_deviation > 0:
        model_df[column] = (
            model_df[column] - model_df[column].mean()
        ) / standard_deviation

exog_columns = [
    "interface",
    "log_pdb",
    "log_partners",
    "review_stars",
]

conditional_model = ConditionalLogit(
    model_df["label"],
    model_df[exog_columns],
    groups=model_df["uniprot_canonical"],
)

conditional_result = conditional_model.fit(
    method="bfgs",
    maxiter=300,
    disp=False,
)

coefficient = float(conditional_result.params["interface"])
standard_error = float(conditional_result.bse["interface"])

conditional_summary = {
    "n_variants": int(len(model_df)),
    "n_proteins": int(model_df["uniprot_canonical"].nunique()),
    "interface_coefficient": coefficient,
    "adjusted_interface_odds_ratio": float(math.exp(coefficient)),
    "ci_lower": float(
        math.exp(coefficient - 1.96 * standard_error)
    ),
    "ci_upper": float(
        math.exp(coefficient + 1.96 * standard_error)
    ),
    "p_value": float(
        conditional_result.pvalues["interface"]
    ),
    "adjusted_for": [
        "protein",
        "log1p_n_pdb_entries",
        "log1p_n_partner_contexts",
        "review_stars",
    ],
}

results = variant_results + position_results
pd.DataFrame(results).to_csv(TABLE, sep="\t", index=False)

summary = {
    "reference_covered_variants": int(len(df)),
    "unique_protein_residues": int(len(position_df)),
    "pure_label_residues": int(len(pure_positions)),
    "mixed_label_residues_excluded": mixed_positions,
    "variant_level_results": variant_results,
    "unique_residue_results": position_results,
    "conditional_logistic_union_result": conditional_summary,
}

with open(SUMMARY, "w") as handle:
    json.dump(summary, handle, indent=2)

print(json.dumps(summary, indent=2))
print("\nResults table:", TABLE)
print("Summary:", SUMMARY)