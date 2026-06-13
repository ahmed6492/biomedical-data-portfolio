"""
Comparative Analysis of Inflammatory Markers Between Control and Treatment Groups
Author: Ahmed Tharwat
"""
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path

HERE = Path(__file__).parent
df = pd.read_csv(HERE / "dataset.csv")

summary = df.groupby("Group")[["IL-6", "TNF-alpha"]].agg(["mean", "std"]).round(2)
summary.to_csv(HERE / "summary_stats.csv")
print("=== Summary statistics ===")
print(summary, "\n")

ctrl = df[df["Group"] == "Control"]
trt = df[df["Group"] == "Treatment"]

il6_t  = stats.ttest_ind(ctrl["IL-6"],       trt["IL-6"],       equal_var=False)
tnf_t  = stats.ttest_ind(ctrl["TNF-alpha"],  trt["TNF-alpha"],  equal_var=False)

print(f"IL-6   Welch t-test: t = {il6_t.statistic:.3f}, p = {il6_t.pvalue:.5f}")
print(f"TNF-a  Welch t-test: t = {tnf_t.statistic:.3f}, p = {tnf_t.pvalue:.5f}")

means = df.groupby("Group")[["IL-6", "TNF-alpha"]].mean()
sds   = df.groupby("Group")[["IL-6", "TNF-alpha"]].std()

fig, ax = plt.subplots(figsize=(7, 4.5))
x = range(len(means.columns))
w = 0.35
ax.bar([i - w/2 for i in x], means.loc["Control"],   w,
       yerr=sds.loc["Control"],   capsize=4, label="Control",   color="#c0392b")
ax.bar([i + w/2 for i in x], means.loc["Treatment"], w,
       yerr=sds.loc["Treatment"], capsize=4, label="Treatment", color="#2980b9")
ax.set_xticks(list(x))
ax.set_xticklabels(["IL-6 (pg/mL)", "TNF-α (pg/mL)"])
ax.set_ylabel("Mean concentration ± SD")
ax.set_title("Inflammatory Markers: Control vs Treatment")
ax.legend()
ax.spines[["top", "right"]].set_visible(False)
plt.tight_layout()
plt.savefig(HERE / "inflammation_chart.png", dpi=200)
print("\nSaved: inflammation_chart.png")

with open(HERE / "results.txt", "w", encoding="utf-8") as f:
    f.write("Comparative Analysis of Inflammatory Markers\n")
    f.write("=" * 50 + "\n\n")
    f.write(str(summary) + "\n\n")
    f.write(f"IL-6   Welch t-test: t = {il6_t.statistic:.3f}, p = {il6_t.pvalue:.5f}\n")
    f.write(f"TNF-a  Welch t-test: t = {tnf_t.statistic:.3f}, p = {tnf_t.pvalue:.5f}\n")
print("Saved: results.txt, summary_stats.csv")
