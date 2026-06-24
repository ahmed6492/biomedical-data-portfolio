# Biomedical Data Portfolio

Three reproducible mini-projects demonstrating biomedical data analysis and junior bioinformatics support work in Python.

The projects are intentionally small and transparent: dataset/source, analysis script, statistical output, figures, and one-page report. They show how I move from biological question to reproducible analysis and clear scientific communication.

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

## Project 3 - Public GEO Lung Cancer Expression Analysis

**Question:** Which expression probes differ most strongly between lung tumor and non-tumor lung tissue samples in public GEO series GSE10072?

**Dataset:** Official public NCBI GEO dataset GSE10072, using processed microarray expression values from 58 tumor and 49 non-tumor lung tissue samples.

**Methods:**

- Official GEO series matrix download and parsing
- Sample metadata extraction and group labeling
- GPL96 probe annotation parsing
- Probe-level Tumor vs Non-tumor mean comparison
- Welch t-test
- Benjamini-Hochberg FDR correction
- Heatmap and volcano plot
- One-page public-data report

**Key result:** The workflow identified 859 significant probes at q < 0.05 and |log2FC| >= 1. Top signals included lower endothelial/vascular markers in tumor tissue and strong SPP1 upregulation in tumor tissue. This project demonstrates public-data handling rather than only simulated data analysis.

**Files:**

| File | Purpose |
|---|---|
| `project3-public-geo-lung-cancer/fetch_data.py` | Downloads and parses official NCBI GEO source files |
| `project3-public-geo-lung-cancer/analysis.py` | Runs probe-level differential-expression analysis |
| `project3-public-geo-lung-cancer/make_report.py` | Builds the one-page PDF report |
| `project3-public-geo-lung-cancer/data/processed/sample_metadata.csv` | Parsed GEO sample metadata |
| `project3-public-geo-lung-cancer/outputs/differential_expression_gse10072.csv` | Full statistics output |
| `project3-public-geo-lung-cancer/figures/gse10072_heatmap.png` | Heatmap |
| `project3-public-geo-lung-cancer/figures/gse10072_volcano.png` | Volcano plot |
| `project3-public-geo-lung-cancer/outputs/GSE10072_Public_GEO_Lung_Cancer_Report.pdf` | Final public-data report |

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

Run Project 3:

```bash
cd project3-public-geo-lung-cancer
pip install -r requirements.txt
python fetch_data.py
python analysis.py
python make_report.py
```

## Notes

- Projects 1 and 2 use simulated datasets for portfolio demonstration.
- Project 3 uses official public NCBI GEO data and keeps raw GEO files out of the repository by default.
- The purpose is to show workflow quality: clean data handling, appropriate basic statistics, figures, and reporting.
- These are junior bioinformatics and biomedical data demonstrations, not a full production RNA-seq or clinical diagnostic pipeline.

## Author

Ahmed Tharwat Hassan  
M.Sc. Immunology and Regenerative Medicine  
Email: athyassen@gmail.com  
LinkedIn: https://www.linkedin.com/in/ahmed-tharwat-b905461a0/
