#!/usr/bin/env python3
"""Local SLIT/ROBO2 sex comparisons on the exported limma expression matrix.

Usage: python3 Local_Sex_Welch.py input.xlsx output.csv
Two-sided Welch tests compare female minus male expression within region.
BH adjustment uses all 32 planned comparisons; non-estimable rows remain NA.
Output: one CSV with 32 rows, four target genes by eight regions, reporting group
sizes, means, SDs, the female-minus-male difference, Welch t/df, raw p, the BH
p over 32 tests, the estimability status and the contributing sample IDs.
Region labels follow the source workbook: SUB denotes substantia nigra, reported
as SN in Supplementary Data 2, and HCP/HPC denote hippocampus.
The 63 rows are tissue samples, not independent donors; the sample-to-donor
mapping of the source workbook is unresolved, so these are sample-level tests.
Dependencies: numpy, pandas, openpyxl, scipy (>=1.11).
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

GENES = ["SLIT1", "SLIT2", "SLIT3", "ROBO2"]
SITES = ["CER", "SUB", "HTL", "OCC", "TEMP", "HPC", "FRO", "PAR"]


def bh32(values):
    """Adjust finite P values with a planned family of 32, preserving missingness."""
    adjusted = np.full(len(values), np.nan)
    indices = np.flatnonzero(np.isfinite(values))
    order = indices[np.argsort(values[indices])]
    ranked = values[order] * 32 / np.arange(1, len(order) + 1)
    adjusted[order] = np.minimum(1, np.minimum.accumulate(ranked[::-1])[::-1])
    return adjusted


def analyse(input_path):
    data = pd.read_excel(input_path, sheet_name="Sheet1")
    required = {"#", "site", "Gender", "Unnamed: 0", "Unnamed: 1", "Unnamed: 2"}
    if not required.issubset(data.columns):
        raise ValueError(f"Missing workbook columns: {sorted(required - set(data.columns))}")
    sample_columns = [c for c in data.columns if isinstance(c, str) and c.startswith("X") and "_S" in c]
    # This is a sample-level comparison; the workbook's Case field is not a donor key.
    metadata = data[["#", "site", "Gender"]].dropna(subset=["#"]).copy()
    metadata["#"] = metadata["#"].astype(str)
    metadata["site"] = metadata["site"].replace({"HCP": "HPC"})
    if metadata["#"].duplicated().any() or set(sample_columns) != set(metadata["#"]):
        raise ValueError("Expression sample IDs and unique metadata IDs must match exactly.")
    if len(metadata) != 63 or set(metadata["site"]) != set(SITES):
        raise ValueError("Expected the original 63-sample, eight-region local cohort.")
    if not metadata["Gender"].isin(["F", "M"]).all():
        raise ValueError("Every local sample must have sex F or M.")
    expression = data[data["Unnamed: 1"].isin(GENES)].rename(columns={
        "Unnamed: 0": "entrez", "Unnamed: 1": "gene", "Unnamed: 2": "ensembl"
    })
    if len(expression) != 4 or set(expression["gene"]) != set(GENES):
        raise ValueError("Expected exactly one expression row for each target gene.")

    rows = []
    for gene in GENES:
        gene_row = expression.loc[expression["gene"] == gene].iloc[0]
        for site in SITES:
            samples, values = {}, {}
            for sex in ["F", "M"]:
                samples[sex] = metadata.loc[(metadata["site"] == site) & (metadata["Gender"] == sex), "#"].tolist()
                array = pd.to_numeric(gene_row[samples[sex]], errors="coerce").to_numpy(dtype=float)
                values[sex] = array[np.isfinite(array)]
            female, male = values["F"], values["M"]
            mean_f = female.mean() if len(female) else np.nan
            mean_m = male.mean() if len(male) else np.nan
            sd_f = female.std(ddof=1) if len(female) >= 2 else np.nan
            sd_m = male.std(ddof=1) if len(male) >= 2 else np.nan
            row = dict(gene=gene, entrez=gene_row["entrez"], ensembl=gene_row["ensembl"],
                       site=site, status="ok", n_F=len(female), n_M=len(male),
                       mean_F=mean_f, mean_M=mean_m, diff_F_minus_M=mean_f-mean_m,
                       sd_F=sd_f, sd_M=sd_m, t_stat=np.nan, df_welch=np.nan, p_value=np.nan,
                       samples_F=",".join(samples["F"]), samples_M=",".join(samples["M"]))
            if len(female) < 2 or len(male) < 2:
                row["status"] = "insufficient_group_size_for_welch"
            elif not np.isfinite(sd_f**2 / len(female) + sd_m**2 / len(male)) or sd_f**2 + sd_m**2 == 0:
                row["status"] = "zero_or_invalid_variance"
            else:
                test = ttest_ind(female, male, equal_var=False, alternative="two-sided")
                row.update(t_stat=float(test.statistic), df_welch=float(test.df), p_value=float(test.pvalue))
            rows.append(row)

    result = pd.DataFrame(rows)
    result["p_adj_BH_m32"] = bh32(result["p_value"].to_numpy())
    result["significant_BH_0_05"] = result["p_adj_BH_m32"] < 0.05
    columns = ["gene", "entrez", "ensembl", "site", "status", "n_F", "n_M", "mean_F", "mean_M",
               "diff_F_minus_M", "sd_F", "sd_M", "t_stat", "df_welch", "p_value", "p_adj_BH_m32",
               "significant_BH_0_05", "samples_F", "samples_M"]
    return result[columns]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Output already exists; choose a new output path.")
    result = analyse(args.input)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False, na_rep="NA")
    print(f"32 comparisons; {(result.status == 'ok').sum()} estimable; "
          f"{result.significant_BH_0_05.sum()} BH-significant.")


if __name__ == "__main__":
    main()
