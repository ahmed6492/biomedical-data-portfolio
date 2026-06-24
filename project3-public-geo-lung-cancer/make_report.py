"""Build a one-page PDF report for the public GEO GSE10072 project."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.platypus import Image, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle


HERE = Path(__file__).resolve().parent
OUTPUT_DIR = HERE / "outputs"
FIGURE_DIR = HERE / "figures"
PDF_PATH = OUTPUT_DIR / "GSE10072_Public_GEO_Lung_Cancer_Report.pdf"


def _p_value(value: float) -> str:
    try:
        return f"{float(value):.2e}"
    except Exception:
        return ""


def build_report() -> None:
    results = pd.read_csv(OUTPUT_DIR / "differential_expression_gse10072.csv")
    sample_metadata = pd.read_csv(HERE / "data" / "processed" / "sample_metadata.csv")
    group_counts = sample_metadata["group"].value_counts()
    sig_count = int(results["significant_q_lt_0_05_abs_log2fc_ge_1"].sum())
    top = results.head(8).copy()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    doc = SimpleDocTemplate(
        str(PDF_PATH),
        pagesize=A4,
        leftMargin=1.25 * cm,
        rightMargin=1.25 * cm,
        topMargin=1.0 * cm,
        bottomMargin=1.0 * cm,
        title="GSE10072 Public GEO Lung Cancer Expression Report",
        author="Ahmed Tharwat Hassan",
    )

    base = getSampleStyleSheet()
    styles = {
        "Title": ParagraphStyle(
            "Title",
            parent=base["Title"],
            fontName="Helvetica-Bold",
            fontSize=14.5,
            leading=17,
            textColor=colors.HexColor("#111827"),
            spaceAfter=3,
        ),
        "Sub": ParagraphStyle(
            "Sub",
            parent=base["BodyText"],
            fontSize=8,
            leading=9.5,
            textColor=colors.HexColor("#4b5563"),
            spaceAfter=6,
        ),
        "Sec": ParagraphStyle(
            "Sec",
            parent=base["Heading2"],
            fontName="Helvetica-Bold",
            fontSize=9.8,
            leading=11,
            textColor=colors.HexColor("#0f766e"),
            spaceBefore=4,
            spaceAfter=2,
        ),
        "Body": ParagraphStyle(
            "Body",
            parent=base["BodyText"],
            fontSize=8.4,
            leading=10.4,
            textColor=colors.HexColor("#111827"),
            spaceAfter=3,
        ),
        "Small": ParagraphStyle(
            "Small",
            parent=base["BodyText"],
            fontSize=7.4,
            leading=8.7,
            textColor=colors.HexColor("#111827"),
            spaceAfter=1.5,
        ),
    }

    story = []
    story.append(Paragraph("Public GEO Lung Cancer Expression Analysis - GSE10072", styles["Title"]))
    story.append(
        Paragraph(
            "Ahmed Tharwat Hassan | Junior bioinformatics portfolio project | Source: NCBI Gene Expression Omnibus",
            styles["Sub"],
        )
    )

    story.append(Paragraph("Objective", styles["Sec"]))
    story.append(
        Paragraph(
            "Compare processed microarray expression values between lung tumor and non-tumor tissue samples "
            "from public GEO series GSE10072, then summarize the strongest probe-level differences using "
            "a reproducible Python workflow.",
            styles["Body"],
        )
    )

    story.append(Paragraph("Dataset and Methods", styles["Sec"]))
    story.append(
        Paragraph(
            f"Dataset: GSE10072 ({int(group_counts.get('Tumor', 0))} tumor samples, "
            f"{int(group_counts.get('Non-tumor', 0))} non-tumor samples). "
            "Analysis used pandas, SciPy, and matplotlib. For each probe, I computed group means, "
            "Tumor vs Non-tumor log2 fold change, Welch two-sample t-test p-values, and "
            "Benjamini-Hochberg FDR q-values. Significant probes were flagged at q < 0.05 and "
            "absolute log2FC >= 1.0.",
            styles["Body"],
        )
    )

    story.append(Paragraph("Results", styles["Sec"]))
    story.append(
        Paragraph(
            f"The workflow identified {sig_count} significant probes by the stated threshold. "
            "The table below lists the top probes by q-value; symbols are shown when available from "
            "the GPL96 platform annotation.",
            styles["Body"],
        )
    )

    data = [["Probe", "Gene", "Non-tumor", "Tumor", "log2FC", "p-value", "q-value"]]
    for _, row in top.iterrows():
        data.append(
            [
                str(row["probe_id"]),
                str(row.get("gene_label", ""))[:13],
                f"{row['non_tumor_mean']:.2f}",
                f"{row['tumor_mean']:.2f}",
                f"{row['log2FC_tumor_vs_non_tumor']:+.2f}",
                _p_value(row["p_value"]),
                _p_value(row["q_value"]),
            ]
        )

    table = Table(
        data,
        colWidths=[2.2 * cm, 2.2 * cm, 1.8 * cm, 1.8 * cm, 1.6 * cm, 2.1 * cm, 2.1 * cm],
    )
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#0f766e")),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ("FONTSIZE", (0, 0), (-1, -1), 6.8),
                ("ALIGN", (2, 0), (-1, -1), "CENTER"),
                ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor("#d0d5dd")),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f3f4f6")]),
            ]
        )
    )
    story.append(table)
    story.append(Spacer(1, 0.18 * cm))

    heatmap = FIGURE_DIR / "gse10072_heatmap.png"
    volcano = FIGURE_DIR / "gse10072_volcano.png"
    if heatmap.exists() and volcano.exists():
        story.append(
            Table(
                [[Image(str(heatmap), width=8.1 * cm, height=6.4 * cm), Image(str(volcano), width=7.8 * cm, height=5.85 * cm)]],
                colWidths=[8.5 * cm, 8.1 * cm],
                style=[("VALIGN", (0, 0), (-1, -1), "TOP")],
            )
        )

    story.append(Paragraph("Interpretation and Limitations", styles["Sec"]))
    story.append(
        Paragraph(
            "This is a compact public-data demonstration of expression analysis, not a clinical diagnostic claim. "
            "It shows data access, metadata parsing, differential-expression statistics, FDR correction, "
            "visualization, and written interpretation. Because the source is microarray data, results are "
            "reported at probe level and should be interpreted with platform annotation context.",
            styles["Small"],
        )
    )

    doc.build(story)
    print(f"Saved: {PDF_PATH}")


if __name__ == "__main__":
    build_report()
