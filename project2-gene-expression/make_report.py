"""1-page PDF report for the gene-expression project."""
from pathlib import Path
import pandas as pd
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle
)

HERE = Path(__file__).parent
res = pd.read_csv(HERE / "differential_expression.csv")

styles = getSampleStyleSheet()
H = ParagraphStyle("H", fontSize=15, leading=18, spaceAfter=4,
                   textColor=colors.HexColor("#1a2a4f"), fontName="Helvetica-Bold")
Sub = ParagraphStyle("Sub", fontSize=9, leading=11, spaceAfter=8, textColor=colors.grey)
Sec = ParagraphStyle("Sec", fontSize=11, leading=13, spaceBefore=6, spaceAfter=2,
                     textColor=colors.HexColor("#1a2a4f"), fontName="Helvetica-Bold")
Body = ParagraphStyle("Body", fontSize=9.5, leading=12.5, spaceAfter=2)

doc = SimpleDocTemplate(
    str(HERE / "Gene_Expression_Report.pdf"),
    pagesize=A4,
    leftMargin=1.6*cm, rightMargin=1.6*cm,
    topMargin=1.3*cm, bottomMargin=1.3*cm,
)
story = []

story.append(Paragraph("Differential Gene-Expression Analysis — Healthy vs Disease (Inflammation Panel)", H))
story.append(Paragraph("Ahmed Tharwat Hassan &middot; M.Sc. Immunology and Regenerative Medicine "
                       "&middot; Portfolio mini-project #2", Sub))

story.append(Paragraph("Objective", Sec))
story.append(Paragraph(
    "To identify genes whose expression differs significantly between Healthy and Disease groups in a "
    "six-gene inflammation/regeneration panel (IL6, TNF, IL1B, IL10, TGFB1, VEGFA), and to demonstrate "
    "a reproducible differential-expression workflow with multiple-testing correction.", Body))

story.append(Paragraph("Methods", Sec))
story.append(Paragraph(
    "Simulated log2-normalised expression for 20 samples (10 Healthy, 10 Disease) was analysed in "
    "Python (pandas, SciPy, matplotlib). For each gene I computed group means, log2 fold change "
    "(Disease vs Healthy), and a Welch two-sample t-test. The Benjamini-Hochberg procedure was "
    "applied across the six genes to control the false-discovery rate (q-value). Genes were flagged "
    "differentially expressed at <b>|log2FC| &gt; 1</b> and <b>q &lt; 0.05</b>.", Body))

story.append(Paragraph("Results", Sec))

data = [["Gene", "Healthy mean", "Disease mean", "log2 FC", "p-value", "q-value", "Sig (q&lt;0.05)"]]
for _, r in res.iterrows():
    data.append([r["Gene"], f"{r['Healthy_mean']:.2f}", f"{r['Disease_mean']:.2f}",
                 f"{r['log2FC']:+.2f}", r["p_value"], f"{r['q_value']:.1e}",
                 "Yes" if r["significant_q<0.05"] else "no"])
tbl = Table(data, colWidths=[1.7*cm, 2.4*cm, 2.4*cm, 1.7*cm, 2.2*cm, 2.2*cm, 2.0*cm])
tbl.setStyle(TableStyle([
    ("BACKGROUND", (0,0), (-1,0), colors.HexColor("#1a2a4f")),
    ("TEXTCOLOR",  (0,0), (-1,0), colors.white),
    ("FONTNAME",   (0,0), (-1,0), "Helvetica-Bold"),
    ("FONTSIZE",   (0,0), (-1,-1), 8.5),
    ("ALIGN",      (1,0), (-1,-1), "CENTER"),
    ("GRID",       (0,0), (-1,-1), 0.4, colors.HexColor("#cccccc")),
    ("ROWBACKGROUNDS", (0,1), (-1,-1), [colors.white, colors.HexColor("#f4f6fa")]),
]))
story.append(tbl)

story.append(Spacer(1, 0.25*cm))
up = res[(res["log2FC"] > 1) & (res["q_value"] < 0.05)]["Gene"].tolist()
dn = res[(res["log2FC"] < -1) & (res["q_value"] < 0.05)]["Gene"].tolist()
others = res[~(((res["log2FC"].abs() > 1) & (res["q_value"] < 0.05)))]["Gene"].tolist()
story.append(Paragraph(
    f"<b>Upregulated in Disease (|log2FC|&gt;1, q&lt;0.05):</b> {', '.join(up) or 'none'}. "
    f"<b>Downregulated:</b> {', '.join(dn) or 'none'}. "
    f"<b>Below threshold:</b> {', '.join(others) or 'none'}.", Body))

heat = HERE / "heatmap.png"; vol = HERE / "volcano.png"
if heat.exists() and vol.exists():
    figs = Table([[Image(str(heat), width=8*cm, height=6.2*cm),
                   Image(str(vol),  width=8*cm, height=5.8*cm)]],
                 colWidths=[8.5*cm, 8.5*cm])
    figs.setStyle(TableStyle([("VALIGN", (0,0), (-1,-1), "TOP")]))
    story.append(Spacer(1, 0.2*cm))
    story.append(figs)

story.append(Paragraph("Conclusion", Sec))
story.append(Paragraph(
    "The pro-inflammatory cytokines IL6, TNF, and IL1B were significantly upregulated in the Disease "
    "group, while the anti-inflammatory cytokine IL10 was significantly downregulated — a pattern "
    "consistent with a shift towards a pro-inflammatory state. VEGFA was upregulated, suggesting "
    "compensatory angiogenic activity. The workflow (data → fold change → t-test → FDR correction "
    "→ heatmap + volcano) mirrors a real differential-expression analysis at the smallest possible "
    "scale and is fully reproducible from a single command.", Body))

doc.build(story)
print("Saved: Gene_Expression_Report.pdf")
