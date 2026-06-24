# Project 3 - Public GEO Lung Cancer Expression Analysis

This project analyzes an official public gene-expression dataset from NCBI GEO:

- **Dataset:** GSE10072
- **Source:** NCBI Gene Expression Omnibus
- **Design:** Lung tumor tissue vs non-tumor lung tissue
- **Platform:** GPL96 Affymetrix Human Genome U133A Array

The goal is to show a reproducible junior bioinformatics workflow using a real public dataset instead of simulated data.

## Question

Which expression probes differ most strongly between lung tumor and non-tumor lung tissue samples in public GEO series GSE10072?

## Workflow

1. Download the official GSE10072 series matrix from NCBI GEO.
2. Parse sample metadata and expression values.
3. Infer sample groups from GEO titles/source names.
4. Join probe IDs to GPL96 platform annotation when available.
5. Run Tumor vs Non-tumor Welch t-tests probe by probe.
6. Apply Benjamini-Hochberg FDR correction.
7. Generate a top-probe table, heatmap, volcano plot, and one-page PDF report.

## Files

| File | Purpose |
|---|---|
| `fetch_data.py` | Downloads and parses official NCBI GEO source files |
| `analysis.py` | Runs differential-expression analysis and generates figures |
| `make_report.py` | Builds the one-page PDF report |
| `requirements.txt` | Python dependencies |
| `data/processed/sample_metadata.csv` | Parsed GEO sample metadata |
| `data/processed/platform_annotation.csv` | GPL96 probe annotation parsed from NCBI GEO |
| `data/processed/expression_preview_top_variable_probes.csv` | Small preview of high-variance probes |
| `outputs/differential_expression_gse10072.csv` | Full probe-level statistics output |
| `outputs/top_30_differential_probes.csv` | Top differential probes by q-value |
| `outputs/results_summary.txt` | Plain-text analysis summary |
| `figures/gse10072_heatmap.png` | Heatmap of top differential probes |
| `figures/gse10072_volcano.png` | Volcano plot |
| `outputs/GSE10072_Public_GEO_Lung_Cancer_Report.pdf` | Final one-page report |

## Reproduce Locally

Install dependencies:

```bash
pip install -r requirements.txt
```

Run the project:

```bash
python fetch_data.py
python analysis.py
python make_report.py
```

## Notes

- Raw GEO files are downloaded into `data/raw/` and are ignored by git.
- The analysis uses GEO processed microarray values supplied in the series matrix.
- Results are reported at probe level. Gene symbols are included when available from GPL96 annotation.
- This is a portfolio demonstration of data handling, statistical analysis, visualization, and scientific reporting. It is not a clinical diagnostic analysis.

## Author

Ahmed Tharwat Hassan  
M.Sc. Immunology and Regenerative Medicine  
Email: athyassen@gmail.com  
LinkedIn: https://www.linkedin.com/in/ahmed-tharwat-b905461a0/  
GitHub: https://github.com/ahmed6492/biomedical-data-portfolio
