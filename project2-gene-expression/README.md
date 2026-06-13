# Gene-Expression Panel Analysis - Healthy vs Disease

End-to-end differential analysis on a simulated six-gene inflammation/regeneration panel: IL6, TNF, IL1B, IL10, TGFB1, and VEGFA.

This compact portfolio project shows fold-change calculation, Welch t-tests, Benjamini-Hochberg FDR correction, heatmap visualisation, volcano plotting, and one-page scientific reporting. It is not presented as a full RNA-seq production pipeline.

## Workflow

`expression_data.csv` -> per-gene t-test -> log2 fold change -> Benjamini-Hochberg FDR -> heatmap + volcano plot -> one-page report.

## Files

| File | Purpose |
|---|---|
| `generate_data.py` | Builds the simulated 20-sample x 6-gene dataset |
| `expression_data.csv` / `expression_data.xlsx` | Input dataset |
| `analysis.py` | Stats and figure generation |
| `make_report.py` | Builds the one-page PDF |
| `heatmap.png`, `volcano.png` | Figures |
| `differential_expression.csv`, `results.txt` | Numerical outputs |
| `Gene_Expression_Report.pdf` | Final report |

## Key Result

All six genes are significant at q < 0.05 in the simulated dataset. The pro-inflammatory markers IL6, TNF, and IL1B are upregulated in Disease, while IL10 is downregulated.
