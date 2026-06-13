"""
Differential gene-expression analysis: Healthy vs Disease.

Workflow:
  expression_data.csv  ->  per-gene t-test  ->  log2 fold change
  ->  Benjamini-Hochberg FDR correction  ->  heatmap + volcano plot
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path

HERE = Path(__file__).parent
df = pd.read_csv(HERE / "expression_data.csv")

GENES = [c for c in df.columns if c not in ("Sample_ID", "Group")]
healthy = df[df["Group"] == "Healthy"]
disease = df[df["Group"] == "Disease"]

# --- Per-gene statistics ----------------------------------------------------
results = []
for g in GENES:
    h = healthy[g].to_numpy()
    d = disease[g].to_numpy()
    t = stats.ttest_ind(d, h, equal_var=False)
    # log2 fold change: disease vs healthy.  Values are already log2-units,
    # so subtraction = log2(disease / healthy).
    log2fc = float(d.mean() - h.mean())
    results.append({
        "Gene": g,
        "Healthy_mean": round(float(h.mean()), 3),
        "Disease_mean": round(float(d.mean()), 3),
        "log2FC": round(log2fc, 3),
        "t_stat": round(float(t.statistic), 3),
        "p_value": float(t.pvalue),
    })
res = pd.DataFrame(results)

# Benjamini-Hochberg FDR correction
m = len(res)
res_sorted = res.sort_values("p_value").reset_index(drop=True)
res_sorted["rank"] = res_sorted.index + 1
res_sorted["q_value"] = (res_sorted["p_value"] * m / res_sorted["rank"]).clip(upper=1.0)
# enforce monotonicity
q = res_sorted["q_value"].to_numpy()
for i in range(len(q) - 2, -1, -1):
    q[i] = min(q[i], q[i + 1])
res_sorted["q_value"] = np.round(q, 5)
res = res_sorted.drop(columns="rank").sort_values("Gene").reset_index(drop=True)
res["p_value"] = res["p_value"].apply(lambda x: f"{x:.2e}")
res["significant_q<0.05"] = res["q_value"] < 0.05
res.to_csv(HERE / "differential_expression.csv", index=False)

print("=== Differential expression results ===")
print(res.to_string(index=False))

# --- Heatmap ----------------------------------------------------------------
mat = df.set_index("Sample_ID")[GENES].copy()
# z-score per gene so the colour scale shows relative expression
mat_z = (mat - mat.mean()) / mat.std()

fig, ax = plt.subplots(figsize=(7, 5.5))
im = ax.imshow(mat_z.values, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
ax.set_xticks(range(len(GENES))); ax.set_xticklabels(GENES, rotation=0)
ax.set_yticks(range(len(mat_z))); ax.set_yticklabels(mat_z.index, fontsize=7)
# colour the sample labels by group
for tick_label, group in zip(ax.get_yticklabels(), df["Group"]):
    tick_label.set_color("#c0392b" if group == "Disease" else "#2980b9")
ax.set_title("Gene-expression heatmap (z-score per gene)")
cb = plt.colorbar(im, ax=ax, shrink=0.7, pad=0.02)
cb.set_label("z-score")
plt.tight_layout()
plt.savefig(HERE / "heatmap.png", dpi=200)
plt.close()
print("Saved: heatmap.png")

# --- Volcano plot ------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 5))
log2fc = [r["log2FC"] for r in results]
qvals  = res.set_index("Gene").loc[GENES, "q_value"].to_numpy()
neglogq = -np.log10(np.clip(qvals, 1e-10, 1.0))
colors = ["#c0392b" if (abs(fc) > 1 and q < 0.05) else "#888888"
          for fc, q in zip(log2fc, qvals)]
ax.scatter(log2fc, neglogq, c=colors, s=140, edgecolor="black", linewidth=0.6)
for g, fc, nq in zip(GENES, log2fc, neglogq):
    ax.annotate(g, (fc, nq), xytext=(6, 4), textcoords="offset points", fontsize=9)
ax.axvline(0, color="#cccccc", lw=0.6)
ax.axhline(-np.log10(0.05), color="#cccccc", lw=0.6, ls="--")
ax.axvline(-1, color="#eeeeee", lw=0.6, ls=":")
ax.axvline(1,  color="#eeeeee", lw=0.6, ls=":")
ax.set_xlabel("log2 fold change (Disease / Healthy)")
ax.set_ylabel("-log10(q-value)")
ax.set_title("Volcano plot — differential expression")
ax.spines[["top", "right"]].set_visible(False)
plt.tight_layout()
plt.savefig(HERE / "volcano.png", dpi=200)
plt.close()
print("Saved: volcano.png")

# --- Plain-text summary ------------------------------------------------------
with open(HERE / "results.txt", "w", encoding="utf-8") as f:
    f.write("Differential gene-expression analysis  -  Disease vs Healthy\n")
    f.write("=" * 64 + "\n\n")
    f.write(res.to_string(index=False) + "\n")
print("Saved: results.txt, differential_expression.csv")
