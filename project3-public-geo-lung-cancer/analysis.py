"""
Probe-level differential-expression analysis for public GEO dataset GSE10072.

The script downloads public GEO files if needed, parses sample metadata, compares
lung tumor vs non-tumor samples, applies Welch t-tests and Benjamini-Hochberg
FDR correction, and writes result tables plus heatmap and volcano figures.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from fetch_data import (
    HERE,
    PROCESSED_DIR,
    build_sample_metadata,
    prepare_processed_files,
    read_series_matrix,
)


OUTPUT_DIR = HERE / "outputs"
FIGURE_DIR = HERE / "figures"


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    p_values = np.asarray(p_values, dtype=float)
    q_values = np.full(p_values.shape, np.nan, dtype=float)
    valid = np.isfinite(p_values)
    if not valid.any():
        return q_values

    p = p_values[valid]
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty_like(p)
    adjusted[order] = np.clip(ranked, 0, 1)
    q_values[valid] = adjusted
    return q_values


def _primary_symbol(value: object) -> str:
    text = "" if pd.isna(value) else str(value)
    return text if text and text != "nan" else ""


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    prepare_processed_files()
    metadata_lines, expression = read_series_matrix()
    sample_metadata = build_sample_metadata(metadata_lines)

    annotation_path = PROCESSED_DIR / "platform_annotation.csv"
    if annotation_path.exists():
        annotation = pd.read_csv(annotation_path, dtype=str).fillna("")
    else:
        annotation = pd.DataFrame(columns=["probe_id", "primary_gene_symbol", "gene_title"])

    return sample_metadata, expression, annotation


def run_analysis() -> pd.DataFrame:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    sample_metadata, expression, annotation = load_inputs()

    usable_metadata = sample_metadata[sample_metadata["group"].isin(["Tumor", "Non-tumor"])].copy()
    tumor_cols = usable_metadata.loc[usable_metadata["group"] == "Tumor", "sample_accession"].tolist()
    non_tumor_cols = usable_metadata.loc[
        usable_metadata["group"] == "Non-tumor", "sample_accession"
    ].tolist()

    if len(tumor_cols) < 2 or len(non_tumor_cols) < 2:
        raise ValueError("Could not identify enough Tumor and Non-tumor samples for comparison.")

    expression = expression[["probe_id", *tumor_cols, *non_tumor_cols]].copy()
    values = expression[tumor_cols + non_tumor_cols].to_numpy(dtype=float)

    transformed = False
    if np.nanpercentile(values, 99) > 100:
        values = np.log2(np.clip(values, a_min=0, a_max=None) + 1)
        transformed = True

    tumor_values = values[:, : len(tumor_cols)]
    non_tumor_values = values[:, len(tumor_cols) :]

    tumor_mean = np.nanmean(tumor_values, axis=1)
    non_tumor_mean = np.nanmean(non_tumor_values, axis=1)
    log2fc = tumor_mean - non_tumor_mean

    ttest = stats.ttest_ind(
        tumor_values,
        non_tumor_values,
        axis=1,
        equal_var=False,
        nan_policy="omit",
    )
    q_values = benjamini_hochberg(ttest.pvalue)

    results = pd.DataFrame(
        {
            "probe_id": expression["probe_id"].astype(str),
            "non_tumor_mean": non_tumor_mean,
            "tumor_mean": tumor_mean,
            "log2FC_tumor_vs_non_tumor": log2fc,
            "t_stat": ttest.statistic,
            "p_value": ttest.pvalue,
            "q_value": q_values,
            "significant_q_lt_0_05_abs_log2fc_ge_1": (q_values < 0.05) & (np.abs(log2fc) >= 1.0),
        }
    )

    if not annotation.empty:
        keep_cols = [c for c in ["probe_id", "primary_gene_symbol", "gene_symbol", "gene_title"] if c in annotation]
        results = results.merge(annotation[keep_cols], on="probe_id", how="left")
    else:
        results["primary_gene_symbol"] = ""
        results["gene_title"] = ""

    results["primary_gene_symbol"] = results.get("primary_gene_symbol", "").map(_primary_symbol)
    results["gene_label"] = np.where(
        results["primary_gene_symbol"].astype(str).str.len() > 0,
        results["primary_gene_symbol"].astype(str),
        results["probe_id"].astype(str),
    )

    results = results.sort_values(["q_value", "p_value"], na_position="last").reset_index(drop=True)
    results.to_csv(OUTPUT_DIR / "differential_expression_gse10072.csv", index=False)
    results.head(30).to_csv(OUTPUT_DIR / "top_30_differential_probes.csv", index=False)

    make_heatmap(expression, usable_metadata, values, results)
    make_volcano(results)
    write_summary(results, sample_metadata, transformed)

    return results


def make_heatmap(
    expression: pd.DataFrame,
    sample_metadata: pd.DataFrame,
    values: np.ndarray,
    results: pd.DataFrame,
) -> None:
    top = results.dropna(subset=["q_value"]).head(30)
    top_probe_ids = top["probe_id"].tolist()
    probe_positions = expression.index[expression["probe_id"].isin(top_probe_ids)].tolist()

    ordered_samples = sample_metadata.sort_values(["group", "sample_accession"])
    ordered_cols = ordered_samples["sample_accession"].tolist()
    all_cols = expression.columns[1:].tolist()
    col_positions = [all_cols.index(col) for col in ordered_cols if col in all_cols]

    heat_values = values[np.ix_(probe_positions, col_positions)]
    row_mean = np.nanmean(heat_values, axis=1, keepdims=True)
    row_std = np.nanstd(heat_values, axis=1, keepdims=True)
    row_std[row_std == 0] = 1
    heat_z = (heat_values - row_mean) / row_std
    heat_z = np.clip(heat_z, -2.5, 2.5)

    labels = (
        top.set_index("probe_id")
        .loc[expression.loc[probe_positions, "probe_id"], "gene_label"]
        .astype(str)
        .tolist()
    )

    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(heat_z, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xticks([])
    ax.set_xlabel("Samples ordered by group: Non-tumor then Tumor")
    ax.set_title("GSE10072 top differential probes - z-score expression")
    colorbar = plt.colorbar(im, ax=ax, shrink=0.72, pad=0.02)
    colorbar.set_label("z-score per probe")
    plt.tight_layout()
    plt.savefig(FIGURE_DIR / "gse10072_heatmap.png", dpi=220)
    plt.close()


def make_volcano(results: pd.DataFrame) -> None:
    plot_df = results.replace([np.inf, -np.inf], np.nan).dropna(
        subset=["log2FC_tumor_vs_non_tumor", "q_value"]
    )
    x = plot_df["log2FC_tumor_vs_non_tumor"].to_numpy()
    y = -np.log10(np.clip(plot_df["q_value"].to_numpy(dtype=float), 1e-300, 1.0))
    sig = (plot_df["q_value"] < 0.05) & (plot_df["log2FC_tumor_vs_non_tumor"].abs() >= 1.0)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(x[~sig], y[~sig], s=9, c="#8a8f98", alpha=0.55, linewidth=0)
    ax.scatter(x[sig], y[sig], s=16, c="#0f766e", alpha=0.75, linewidth=0)

    label_df = (
        pd.concat(
            [
                plot_df[sig].sort_values("q_value").head(4),
                plot_df[sig].sort_values("log2FC_tumor_vs_non_tumor").head(2),
                plot_df[sig].sort_values("log2FC_tumor_vs_non_tumor", ascending=False).head(2),
            ]
        )
        .drop_duplicates("probe_id")
        .head(8)
        .reset_index(drop=True)
    )
    for idx, row in label_df.iterrows():
        point_x = float(row["log2FC_tumor_vs_non_tumor"])
        if point_x < -3.2:
            dx, dy, ha = 8, 0, "left"
        elif point_x < 0:
            dx, dy, ha = -8, 8 if idx % 2 == 0 else -8, "right"
        elif point_x > 3.2:
            dx, dy, ha = -8, 0, "right"
        else:
            dx, dy, ha = 8, 8 if idx % 2 == 0 else -8, "left"
        ax.annotate(
            row["gene_label"],
            (row["log2FC_tumor_vs_non_tumor"], -np.log10(max(row["q_value"], 1e-300))),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=7,
            ha=ha,
            bbox={"boxstyle": "round,pad=0.15", "fc": "white", "ec": "none", "alpha": 0.72},
        )

    ax.axhline(-np.log10(0.05), color="#d0d5dd", lw=0.8, ls="--")
    ax.axvline(-1, color="#d0d5dd", lw=0.8, ls=":")
    ax.axvline(1, color="#d0d5dd", lw=0.8, ls=":")
    ax.axvline(0, color="#e5e7eb", lw=0.8)
    ax.set_xlabel("log2 fold change (Tumor vs Non-tumor)")
    ax.set_ylabel("-log10(q-value)")
    ax.set_title("GSE10072 volcano plot")
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    plt.savefig(FIGURE_DIR / "gse10072_volcano.png", dpi=220)
    plt.close()


def write_summary(results: pd.DataFrame, sample_metadata: pd.DataFrame, transformed: bool) -> None:
    significant = results[results["significant_q_lt_0_05_abs_log2fc_ge_1"]]
    top = results.head(10)
    group_counts = sample_metadata["group"].value_counts()

    with open(OUTPUT_DIR / "results_summary.txt", "w", encoding="utf-8") as handle:
        handle.write("Public GEO lung cancer expression analysis - GSE10072\n")
        handle.write("=" * 64 + "\n\n")
        handle.write("Source: NCBI GEO GSE10072 series matrix\n")
        handle.write("Comparison: Tumor vs Non-tumor lung tissue\n")
        handle.write(f"Samples: {group_counts.to_dict()}\n")
        handle.write(f"Expression values log2-transformed by script: {transformed}\n")
        handle.write(
            "Significance rule: q < 0.05 and absolute log2 fold change >= 1.0\n"
        )
        handle.write(f"Significant probes: {len(significant)}\n\n")
        handle.write("Top 10 probes by q-value:\n")
        handle.write(
            top[
                [
                    "probe_id",
                    "gene_label",
                    "log2FC_tumor_vs_non_tumor",
                    "p_value",
                    "q_value",
                ]
            ].to_string(index=False)
        )
        handle.write("\n")


if __name__ == "__main__":
    res = run_analysis()
    sig_count = int(res["significant_q_lt_0_05_abs_log2fc_ge_1"].sum())
    print(f"Saved outputs to: {OUTPUT_DIR}")
    print(f"Saved figures to: {FIGURE_DIR}")
    print(f"Significant probes: {sig_count}")
