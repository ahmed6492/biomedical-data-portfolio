# Biomedical Data Portfolio

Two reproducible mini-projects demonstrating biomedical data analysis in Python.

The projects are intentionally small and transparent: dataset, analysis script, statistical output, figures, and one-page report. They show how I move from biological question to reproducible analysis and clear scientific communication.

## Project 1 - Comparative Analysis of Inflammatory Markers

**Question:** Does a simulated treatment reduce inflammatory cytokines compared with control?

**Dataset:** 20 simulated patient records, split into Control and Treatment groups, with IL-6 and TNF-alpha values.

**Methods:**

- Group means and standard deviations
- Welch two-sample t-test
- Bar chart with error bars
- One-page scientific report

**Key result:** IL-6 and TNF-alpha were significantly lower in the treatment group, consistent with an anti-inflammatory effect.

**Files:**

| File | Purpose |
|---|---|
| `dataset.csv` / `dataset.xlsx` | Input dataset |
| `analysis.py` | Statistical analysis and chart generation |
| `make_report.py` | Builds the one-page PDF report |
| `summary_stats.csv` | Group means and standard deviations |
| `results.txt` | Welch t-test output |
| `inflammation_chart.png` | Figure |
| `Portfolio_Report.pdf` | Final report |

## Project 2 - Differential Gene-Expression Panel Analysis

**Question:** Which genes in a small inflammation/regeneration panel differ between Healthy and Disease groups?

**Dataset:** Simulated 20-sample, six-gene panel: IL6, TNF, IL1B, IL10, TGFB1, VEGFA.

**Methods:**

- Per-gene mean expression
- log2 fold change
- Welch t-test
- Benjamini-Hochberg FDR correction
- z-score heatmap
- volcano plot
- One-page report

**Key result:** IL6, TNF, and IL1B were upregulated in Disease, while IL10 was downregulated. All six genes were significant at q < 0.05 in the simulated dataset.

**Files:**

| File | Purpose |
|---|---|
| `project2-gene-expression/generate_data.py` | Creates the simulated dataset |
| `project2-gene-expression/expression_data.csv` | Input dataset |
| `project2-gene-expression/analysis.py` | Differential analysis and figure generation |
| `project2-gene-expression/differential_expression.csv` | Statistics output |
| `project2-gene-expression/heatmap.png` | Heatmap |
| `project2-gene-expression/volcano.png` | Volcano plot |
| `project2-gene-expression/Gene_Expression_Report.pdf` | Final report |

## Reproduce Locally

Install dependencies:

```bash
pip install pandas matplotlib scipy reportlab openpyxl
```

Run Project 1:

```bash
python analysis.py
python make_report.py
```

Run Project 2:

```bash
cd project2-gene-expression
python generate_data.py
python analysis.py
python make_report.py
```

## Notes

- The datasets are simulated for portfolio demonstration.
- The purpose is to show workflow quality: clean data handling, appropriate basic statistics, figures, and reporting.
- This is a gene-expression panel analysis, not a full production RNA-seq pipeline.

## Author

Ahmed Tharwat Hassan  
M.Sc. Immunology and Regenerative Medicine  
Email: athyassen@gmail.com  
LinkedIn: https://www.linkedin.com/in/ahmed-tharwat-b905461a0/
