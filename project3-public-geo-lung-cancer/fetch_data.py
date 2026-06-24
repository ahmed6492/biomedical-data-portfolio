"""
Download and parse the public GEO source files for project 3.

Dataset: GSE10072, lung tumor and non-tumor expression profiling by array.
Source: NCBI Gene Expression Omnibus.
"""
from __future__ import annotations

import csv
import gzip
import io
import re
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
RAW_DIR = HERE / "data" / "raw"
PROCESSED_DIR = HERE / "data" / "processed"

SERIES_ACCESSION = "GSE10072"
PLATFORM_ACCESSION = "GPL96"

SERIES_MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE10nnn/GSE10072/"
    "matrix/GSE10072_series_matrix.txt.gz"
)
PLATFORM_ANNOTATION_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/"
    "annot/GPL96.annot.gz"
)

SERIES_MATRIX_PATH = RAW_DIR / "GSE10072_series_matrix.txt.gz"
PLATFORM_ANNOTATION_PATH = RAW_DIR / "GPL96.annot.gz"


def _download(url: str, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and target.stat().st_size > 0:
        print(f"Already exists: {target}")
        return

    print(f"Downloading: {url}")
    with urllib.request.urlopen(url, timeout=120) as response:
        target.write_bytes(response.read())
    print(f"Saved: {target}")


def ensure_raw_files() -> None:
    """Download the public source files if they are not already present."""
    _download(SERIES_MATRIX_URL, SERIES_MATRIX_PATH)
    try:
        _download(PLATFORM_ANNOTATION_URL, PLATFORM_ANNOTATION_PATH)
    except Exception as exc:  # annotation is useful, but the project can run without it
        print(f"Warning: platform annotation download failed: {exc}")


def _clean_cell(value: str) -> str:
    return value.strip().strip('"')


def read_series_matrix(matrix_path: Path = SERIES_MATRIX_PATH) -> tuple[dict[str, list[str]], pd.DataFrame]:
    """
    Return sample metadata lines and the expression matrix from a GEO series matrix file.

    The expression table is probe rows x sample columns. Values are GEO processed
    expression values supplied by the series submitter.
    """
    metadata: dict[str, list[str]] = {}
    table_lines: list[str] = []
    in_table = False

    with gzip.open(matrix_path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line == "!series_matrix_table_begin":
                in_table = True
                continue
            if line == "!series_matrix_table_end":
                in_table = False
                continue
            if in_table:
                table_lines.append(line)
                continue
            if line.startswith("!Sample_"):
                parts = next(csv.reader([line], delimiter="\t"))
                metadata[parts[0].lstrip("!")] = [_clean_cell(v) for v in parts[1:]]

    if not table_lines:
        raise ValueError(f"No expression table found in {matrix_path}")

    expression = pd.read_csv(io.StringIO("\n".join(table_lines)), sep="\t")
    expression = expression.rename(columns={expression.columns[0]: "probe_id"})

    for column in expression.columns[1:]:
        expression[column] = pd.to_numeric(expression[column], errors="coerce")

    return metadata, expression


def _infer_group(title: str, source_name: str) -> str:
    text = f"{title} {source_name}".lower()
    if "normal" in text or "non-tumor" in text or "nontumor" in text:
        return "Non-tumor"
    if "tumor" in text or "cancer" in text or "adenocarcinoma" in text:
        return "Tumor"
    return "Unknown"


def _infer_subject_id(title: str) -> str:
    match = re.search(r"\bGT\d+\b", title)
    return match.group(0) if match else ""


def build_sample_metadata(metadata: dict[str, list[str]]) -> pd.DataFrame:
    accessions = metadata.get("Sample_geo_accession")
    if not accessions:
        raise ValueError("GEO sample accessions are missing from the series matrix metadata.")

    titles = metadata.get("Sample_title", [""] * len(accessions))
    sources = metadata.get("Sample_source_name_ch1", [""] * len(accessions))
    organisms = metadata.get("Sample_organism_ch1", [""] * len(accessions))

    rows = []
    for accession, title, source, organism in zip(accessions, titles, sources, organisms):
        rows.append(
            {
                "sample_accession": accession,
                "sample_title": title,
                "source_name": source,
                "organism": organism,
                "group": _infer_group(title, source),
                "subject_id": _infer_subject_id(title),
            }
        )

    return pd.DataFrame(rows)


def _read_platform_annotation(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(columns=["probe_id", "gene_symbol", "gene_title", "entrez_gene_id"])

    rows: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        in_table = False
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("!platform_table_begin"):
                in_table = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if line.startswith("#"):
                continue
            if in_table or line.startswith("ID\t"):
                rows.append(line)

    if not rows:
        return pd.DataFrame(columns=["probe_id", "gene_symbol", "gene_title", "entrez_gene_id"])

    annot = pd.read_csv(io.StringIO("\n".join(rows)), sep="\t", dtype=str)
    annot.columns = [str(c).strip() for c in annot.columns]

    def find_col(candidates: tuple[str, ...]) -> str | None:
        normalized = {c.lower().replace("_", " ").strip(): c for c in annot.columns}
        for candidate in candidates:
            if candidate.lower() in normalized:
                return normalized[candidate.lower()]
        return None

    id_col = find_col(("id", "probe id")) or annot.columns[0]
    symbol_col = find_col(("gene symbol", "gene symbols", "symbol"))
    title_col = find_col(("gene title", "gene title"))
    entrez_col = find_col(("entrez gene id", "entrez_gene_id", "gene id"))

    out = pd.DataFrame({"probe_id": annot[id_col].astype(str)})
    out["gene_symbol"] = annot[symbol_col].fillna("").astype(str) if symbol_col else ""
    out["gene_title"] = annot[title_col].fillna("").astype(str) if title_col else ""
    out["entrez_gene_id"] = annot[entrez_col].fillna("").astype(str) if entrez_col else ""

    for col in ("gene_symbol", "gene_title", "entrez_gene_id"):
        out[col] = out[col].replace({"nan": "", "---": ""})

    out["primary_gene_symbol"] = out["gene_symbol"].map(
        lambda value: re.split(r"\s*///\s*", value)[0] if value else ""
    )
    return out.drop_duplicates("probe_id")


def prepare_processed_files() -> None:
    ensure_raw_files()
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    metadata, expression = read_series_matrix()
    sample_metadata = build_sample_metadata(metadata)
    annotation = _read_platform_annotation(PLATFORM_ANNOTATION_PATH)

    sample_metadata.to_csv(PROCESSED_DIR / "sample_metadata.csv", index=False)
    annotation.to_csv(PROCESSED_DIR / "platform_annotation.csv", index=False)

    sample_cols = [c for c in expression.columns if c in set(sample_metadata["sample_accession"])]
    preview = expression[["probe_id", *sample_cols]].copy()
    preview["variance"] = np.nanvar(preview[sample_cols].to_numpy(dtype=float), axis=1)
    preview = preview.sort_values("variance", ascending=False).head(50).drop(columns="variance")
    preview.to_csv(PROCESSED_DIR / "expression_preview_top_variable_probes.csv", index=False)

    print("Processed files written:")
    print(f"- {PROCESSED_DIR / 'sample_metadata.csv'}")
    print(f"- {PROCESSED_DIR / 'platform_annotation.csv'}")
    print(f"- {PROCESSED_DIR / 'expression_preview_top_variable_probes.csv'}")
    print("Sample groups:")
    print(sample_metadata["group"].value_counts().to_string())


if __name__ == "__main__":
    prepare_processed_files()
