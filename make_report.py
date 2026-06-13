"""
Builds the 1-page PDF report from analysis outputs.
Run analysis.py FIRST so that inflammation_chart.png and summary_stats.csv exist.
"""
import pandas as pd
from scipy import stats
from pathlib import Path
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle
)

HERE = Path(__file__).parent
df   = pd.read_csv(HERE / "dataset.csv")
ctrl = df[df["Group"] == "Control"]
trt  = df[df["Group"] == "Treatment"]

il6_t = stats.ttest_ind(ctrl["IL-6"],      trt["IL-6"],      equal_var=False)
tnf_t = stats.ttest_ind(ctrl["TNF-alpha"], trt["TNF-alpha"], equal_var=False)

styles = getSampleStyleSheet()
styles.add(ParagraphStyle(name="H", fontSize=15, leading=18, spaceAfter=4,
                          textColor=colors.HexColor("#1a2a4f"), fontName="Helvetica-Bold"))
styles.add(ParagraphStyle(name="Sub", fontSize=9, leading=11, spaceAfter=10,
                          textColor=colors.grey))
styles.add(ParagraphStyle(name="Sec", fontSize=11, leading=14, spaceBefore=6, spaceAfter=2,
                          textColor=colors.HexColor("#1a2a4f"), fontName="Helvetica-Bold"))
styles.add(ParagraphStyle(name="Body2", fontSize=9.5, leading=12.5, spaceAfter=2))

doc = SimpleDocTemplate(
    str(HERE / "Portfolio_Report.pdf"),
    pagesize=A4,
    leftMargin=1.8*cm, rightMargin=1.8*cm,
    topMargin=1.5*cm, bottomMargin=1.5*cm,
)

story = []
story.append(Paragraph("Comparative Analysis of Inflammatory Markers Between Control and Treatment Groups", styles["H"]))
story.append(Paragraph("Ahmed Tharwat &middot; M.Sc. Immunology and Regenerative Medicine &middot; Portfolio mini-project", styles["Sub"]))

story.append(Paragraph("Objective", styles["Sec"]))
story.append(Paragraph(
    "To evaluate whether a simulated treatment reduces the pro-inflammatory cytokines IL-6 and TNF-α "
    "compared to a control group, and to demonstrate a reproducible biomedical-data workflow "
    "(dataset &rarr; descriptive statistics &rarr; inferential test &rarr; visualisation &rarr; interpretation).",
    styles["Body2"]))

story.append(Paragraph("Methods", styles["Sec"]))
story.append(Paragraph(
    "A simulated dataset of 20 patients (10 Control, 10 Treatment) with measurements of IL-6 and TNF-α (pg/mL) "
    "was analysed in Python (pandas, SciPy, matplotlib). Group means and standard deviations were computed. "
    "Between-group differences were tested with Welch's two-sample t-test (unequal variances assumed). "
    "Significance threshold: p &lt; 0.05.",
    styles["Body2"]))

story.append(Paragraph("Results", styles["Sec"]))

mean_ctrl_il6 = ctrl["IL-6"].mean();      sd_ctrl_il6 = ctrl["IL-6"].std()
mean_trt_il6  = trt["IL-6"].mean();       sd_trt_il6  = trt["IL-6"].std()
mean_ctrl_tnf = ctrl["TNF-alpha"].mean(); sd_ctrl_tnf = ctrl["TNF-alpha"].std()
mean_trt_tnf  = trt["TNF-alpha"].mean();  sd_trt_tnf  = trt["TNF-alpha"].std()

data = [
    ["Group", "IL-6 (mean ± SD)", "TNF-α (mean ± SD)", "n"],
    ["Control",   f"{mean_ctrl_il6:.2f} ± {sd_ctrl_il6:.2f}", f"{mean_ctrl_tnf:.2f} ± {sd_ctrl_tnf:.2f}", "10"],
    ["Treatment", f"{mean_trt_il6:.2f} ± {sd_trt_il6:.2f}",   f"{mean_trt_tnf:.2f} ± {sd_trt_tnf:.2f}",   "10"],
]
tbl = Table(data, colWidths=[3.2*cm, 4.8*cm, 4.8*cm, 1.5*cm])
tbl.setStyle(TableStyle([
    ("BACKGROUND", (0,0), (-1,0), colors.HexColor("#1a2a4f")),
    ("TEXTCOLOR",  (0,0), (-1,0), colors.white),
    ("FONTNAME",   (0,0), (-1,0), "Helvetica-Bold"),
    ("FONTSIZE",   (0,0), (-1,-1), 9),
    ("ALIGN",      (1,0), (-1,-1), "CENTER"),
    ("GRID",       (0,0), (-1,-1), 0.4, colors.HexColor("#cccccc")),
    ("ROWBACKGROUNDS", (0,1), (-1,-1), [colors.white, colors.HexColor("#f4f6fa")]),
    ("BOTTOMPADDING", (0,0), (-1,-1), 4),
    ("TOPPADDING",    (0,0), (-1,-1), 4),
]))
story.append(tbl)
story.append(Spacer(1, 0.25*cm))

story.append(Paragraph(
    f"Welch's t-test &mdash; IL-6: t = {il6_t.statistic:.2f}, p = {il6_t.pvalue:.2e}. "
    f"TNF-α: t = {tnf_t.statistic:.2f}, p = {tnf_t.pvalue:.2e}. "
    "Both cytokines were significantly lower in the Treatment group.",
    styles["Body2"]))

img_path = HERE / "inflammation_chart.png"
if img_path.exists():
    story.append(Spacer(1, 0.2*cm))
    story.append(Image(str(img_path), width=13*cm, height=8*cm))

story.append(Paragraph("Conclusion", styles["Sec"]))
story.append(Paragraph(
    "The treatment group showed markedly lower mean IL-6 and TNF-α than the control group, consistent with "
    "an anti-inflammatory effect. The result is statistically significant for both markers, suggesting the "
    "treatment may modulate pro-inflammatory cytokine activity &mdash; a finding relevant to immunological "
    "and regenerative-medicine contexts. Tools used: Python (pandas, SciPy, matplotlib).",
    styles["Body2"]))

doc.build(story)
print("Saved: Portfolio_Report.pdf")
