"""
Generate a simulated gene-expression dataset for a portfolio project.

Design:
- 6 genes (inflammation / regenerative axes), 20 samples (10 Healthy, 10 Disease).
- Disease pattern: pro-inflammatory genes UP, anti-inflammatory genes DOWN,
  VEGFA UP (compensatory angiogenic response).
- Values are log2-normalised expression units (typical range 3–9).

Reproducible: fixed seed so re-running produces the same dataset.
"""
import numpy as np
import pandas as pd
from pathlib import Path

HERE = Path(__file__).parent
rng = np.random.default_rng(seed=42)

GENES = ["IL6", "TNF", "IL1B", "IL10", "TGFB1", "VEGFA"]

# (gene -> (healthy_mean, disease_mean, sd))
PROFILES = {
    "IL6":   (4.0, 7.2, 0.8),   # pro-inflammatory up
    "TNF":   (4.5, 7.0, 0.9),   # pro-inflammatory up
    "IL1B":  (4.2, 6.9, 0.8),   # pro-inflammatory up
    "IL10":  (5.5, 4.2, 0.7),   # anti-inflammatory down
    "TGFB1": (6.0, 5.2, 0.7),   # regenerative/anti-inflammatory mildly down
    "VEGFA": (5.0, 6.4, 0.9),   # angiogenic compensatory up
}

N_PER_GROUP = 10
rows = []
for i in range(1, N_PER_GROUP + 1):
    sample = {"Sample_ID": f"H{i:02d}", "Group": "Healthy"}
    for g in GENES:
        h_mean, d_mean, sd = PROFILES[g]
        sample[g] = round(float(rng.normal(h_mean, sd)), 3)
    rows.append(sample)
for i in range(1, N_PER_GROUP + 1):
    sample = {"Sample_ID": f"D{i:02d}", "Group": "Disease"}
    for g in GENES:
        h_mean, d_mean, sd = PROFILES[g]
        sample[g] = round(float(rng.normal(d_mean, sd)), 3)
    rows.append(sample)

df = pd.DataFrame(rows)
df.to_csv(HERE / "expression_data.csv", index=False)
df.to_excel(HERE / "expression_data.xlsx", index=False)
print(f"Saved: expression_data.csv  ({len(df)} samples x {len(GENES)} genes)")
print(df.head().to_string(index=False))
