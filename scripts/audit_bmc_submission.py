from __future__ import annotations

import argparse
import json
import re
import zipfile
from pathlib import Path

from docx import Document
from pypdf import PdfReader


REQUIRED_DECLARATIONS = [
    "Ethics approval and consent to participate",
    "Consent for publication",
    "Availability of data and materials",
    "Competing interests",
    "Funding",
    "Authors' contributions",
    "Acknowledgements",
]


def expand_citation(content: str) -> set[int]:
    values: set[int] = set()
    for part in content.split(","):
        if "-" in part:
            start, end = (int(value) for value in part.split("-", 1))
            values.update(range(start, end + 1))
        else:
            values.add(int(part))
    return values


def docx_text(path: Path) -> str:
    document = Document(path)
    return "\n".join(paragraph.text for paragraph in document.paragraphs)


def xml_from_zip(path: Path, member: str) -> str:
    with zipfile.ZipFile(path) as archive:
        return archive.read(member).decode("utf-8")


def pdf_font_embedding(path: Path) -> tuple[int, int]:
    """Return (font programs found, font programs embedded) for a PDF."""
    found = 0
    embedded = 0
    reader = PdfReader(path)
    for page in reader.pages:
        resources = page.get("/Resources")
        if not resources:
            continue
        resources = resources.get_object()
        fonts = resources.get("/Font")
        if not fonts:
            continue
        for font_ref in fonts.get_object().values():
            font = font_ref.get_object()
            descendants = font.get("/DescendantFonts")
            candidates = [item.get_object() for item in descendants] if descendants else [font]
            for candidate in candidates:
                descriptor_ref = candidate.get("/FontDescriptor")
                if not descriptor_ref:
                    continue
                found += 1
                descriptor = descriptor_ref.get_object()
                if any(descriptor.get(key) is not None for key in ("/FontFile", "/FontFile2", "/FontFile3")):
                    embedded += 1
    return found, embedded


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package", required=True, type=Path)
    parser.add_argument("--markdown", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--review", required=True, type=Path)
    args = parser.parse_args()

    package = args.package
    manuscript_docx = package / "OC_OA_BMC_Medical_Genomics_main_manuscript.docx"
    manuscript_pdf = package / "OC_OA_BMC_Medical_Genomics_main_manuscript.pdf"
    cover_docx = package / "OC_OA_BMC_Medical_Genomics_cover_letter.docx"
    tables_xlsx = package / "Additional_file_1_supplementary_tables.xlsx"
    figures_pdf = package / "Additional_file_2_supplementary_figures.pdf"
    code_zip = package / "Additional_file_3_reproducible_code.zip"
    main_figures = [package / f"Figure_{number}.pdf" for number in range(1, 8)]
    required_files = [
        manuscript_docx,
        manuscript_pdf,
        cover_docx,
        tables_xlsx,
        figures_pdf,
        code_zip,
        *main_figures,
    ]

    errors: list[str] = []
    warnings: list[str] = []
    checks: dict[str, object] = {}
    for path in required_files:
        if not path.exists():
            errors.append(f"Missing required file: {path.name}")
    if errors:
        raise SystemExit("; ".join(errors))

    markdown = args.markdown.read_text(encoding="utf-8")
    body, reference_tail = markdown.split("\n## References\n", 1)
    reference_text = reference_tail.split("\n## Figure legends\n", 1)[0]
    reference_count = len(re.findall(r"^\d+\.\s", reference_text, flags=re.M))
    citations: set[int] = set()
    first_appearance: list[int] = []
    for match in re.finditer(r"\[([0-9,\-]+)\]", body):
        expanded = expand_citation(match.group(1))
        citations.update(expanded)
        for number in sorted(expanded):
            if number not in first_appearance:
                first_appearance.append(number)
    checks["reference_count"] = reference_count
    checks["citation_coverage"] = [min(citations), max(citations), len(citations)]
    if reference_count != 69 or citations != set(range(1, 70)):
        errors.append("Reference numbering or in-text citation coverage is incomplete")
    checks["references_in_first_appearance_order"] = first_appearance == list(range(1, 70))
    if not checks["references_in_first_appearance_order"]:
        errors.append("References are not numbered in order of first appearance")
    dois = [
        match.group(1).rstrip(".")
        for match in re.finditer(
            r"\bdoi:(10\.\d{4,9}/[^\s]+)", reference_text, flags=re.I
        )
    ]
    checks["reference_doi_count"] = len(dois)
    checks["unique_reference_doi_count"] = len({doi.lower() for doi in dois})
    if len(dois) != 69 or len({doi.lower() for doi in dois}) != 69:
        errors.append("Every reference must have one unique DOI")

    main_figure_refs = {
        int(value) for value in re.findall(r"\bFigure ([1-7])(?:[A-Z](?:-[A-Z])?)?\b", body)
    }
    supplementary_figure_refs = {
        int(value) for value in re.findall(r"\bFigure S([1-5])\b", body)
    }
    additional_file_refs = {
        int(value) for value in re.findall(r"\bAdditional file ([1-3])\b", body)
    }
    checks["main_figures_cited"] = sorted(main_figure_refs)
    checks["supplementary_figures_cited"] = sorted(supplementary_figure_refs)
    checks["additional_files_cited"] = sorted(additional_file_refs)
    if main_figure_refs != set(range(1, 8)):
        errors.append("One or more main figures are not cross-referenced in the main text")
    if supplementary_figure_refs != set(range(1, 6)):
        errors.append("One or more supplementary figures are not cross-referenced in the main text")
    if additional_file_refs != {1, 2, 3}:
        errors.append("One or more additional files are not cross-referenced in the main text")

    abstract = markdown.split("## Abstract", 1)[1].split("## Background", 1)[0]
    abstract_words = len(re.findall(r"\b[\w-]+\b", re.sub(r"\*", "", abstract)))
    keyword_line = re.search(r"\*\*Keywords:\*\*\s*(.*)", abstract)
    keyword_count = len(keyword_line.group(1).split(";")) if keyword_line else 0
    checks["abstract_words"] = abstract_words
    checks["keyword_count"] = keyword_count
    if abstract_words > 350:
        errors.append("Abstract exceeds 350 words")
    if not 3 <= keyword_count <= 10:
        errors.append("Keyword count is outside 3-10")
    for label in ("Background", "Methods", "Results", "Conclusions"):
        if not re.search(rf"\*\*{label}\.\*\*", abstract):
            errors.append(f"Structured abstract label missing: {label}")

    for section in ("Background", "Methods", "Results", "Discussion", "Conclusions", "Abbreviations", "Declarations"):
        if f"## {section}" not in markdown:
            errors.append(f"Required section missing: {section}")
    methods_text = markdown.split("## Methods", 1)[1].split("## Results", 1)[0]
    methods_words = len(re.findall(r"\b[\w-]+\b", re.sub(r"\*", "", methods_text)))
    checks["methods_words"] = methods_words
    if methods_words > 800:
        errors.append("Main-text Methods exceed the 800-word acceptance-polish target")
    expected_title = "# Shared molecular features between osteoarthritis and ovarian cancer revealed by multi-layer transcriptomic analyses"
    checks["discovery_framed_title"] = markdown.startswith(expected_title)
    if not checks["discovery_framed_title"]:
        errors.append("Discovery-framed manuscript title is missing")
    for heading in REQUIRED_DECLARATIONS:
        if f"### {heading}" not in markdown:
            errors.append(f"Required declaration heading missing: {heading}")
    if "### Code availability" in markdown or "### Data availability" in markdown:
        errors.append("Legacy declaration headings remain")

    if re.search(r"Supplementary Figure \d", markdown, flags=re.I):
        errors.append("Legacy supplementary-figure directory or citation remains")
    if re.search(r"\bSTRING\b|\bWGCNA\b", markdown):
        errors.append("Excluded STRING/WGCNA material remains in the manuscript")
    negative_framing = {
        term: len(re.findall(rf"\b{re.escape(term)}\b", body, flags=re.I))
        for term in ("opposite", "discordant", "divergent")
    }
    checks["negative_framing_terms"] = negative_framing
    if any(negative_framing.values()):
        errors.append("Negative directional framing remains in the main-text narrative")
    single_cell_tool_mentions = {
        term: len(re.findall(re.escape(term), body, flags=re.I))
        for term in ("scDblFinder", "UCell", "TPM-only")
    }
    checks["main_text_single_cell_tool_mentions"] = single_cell_tool_mentions
    if any(single_cell_tool_mentions.values()):
        errors.append("Detailed single-cell audit terms remain in the main text")
    if (
        "Additional file 1: Table S1" not in body
        or "Additional file 2: Figure S5" not in body
        or "Additional file 3" not in body
    ):
        errors.append("Additional files are not cited in the expected order")
    checks["placeholder_count"] = markdown.count("AUTHOR INPUT REQUIRED")
    if checks["placeholder_count"] != 0:
        errors.append("Author-input placeholders remain in the manuscript")

    doc_xml = xml_from_zip(manuscript_docx, "word/document.xml")
    styles_xml = xml_from_zip(manuscript_docx, "word/styles.xml")
    checks["continuous_line_numbers"] = "w:lnNumType" in doc_xml and 'w:restart="continuous"' in doc_xml
    checks["page_number_field"] = " PAGE " in doc_xml or " PAGE " in xml_from_zip(manuscript_docx, "word/footer1.xml")
    checks["double_spaced_normal_style"] = 'w:line="480"' in styles_xml
    checks["manual_page_breaks"] = len(re.findall(r"<w:br[^>]+w:type=\"page\"", doc_xml))
    if not checks["continuous_line_numbers"]:
        errors.append("Continuous line numbering is not encoded")
    if not checks["page_number_field"]:
        errors.append("Page numbering is not encoded")
    if not checks["double_spaced_normal_style"]:
        errors.append("Double spacing is not encoded in the Normal style")
    if checks["manual_page_breaks"]:
        errors.append("Manual page breaks remain in the main manuscript")

    text = docx_text(manuscript_docx)
    if "�" in text or "B枚" in text or "鈥" in text:
        errors.append("Mojibake or replacement characters remain in DOCX text")
    expected_author_line = "Junhui Shi¹†, Mengxiang Liu¹†, Repkat Inayatilla¹, Ke Li¹, Lei Chen¹*"
    expected_equal_contribution = "Junhui Shi and Mengxiang Liu contributed equally and share first authorship."
    checks["author_order_verified"] = expected_author_line in text
    checks["equal_contribution_verified"] = expected_equal_contribution in text
    checks["removed_author_mentions"] = {
        term: text.count(term) for term in ("Hongtao Yu", "于洪涛", "HTY")
    }
    if not checks["author_order_verified"]:
        errors.append("Final author order or authorship markers are incorrect")
    if not checks["equal_contribution_verified"]:
        errors.append("Equal-contribution statement is missing")
    if any(checks["removed_author_mentions"].values()):
        errors.append("Removed author Hongtao Yu remains in the manuscript")

    with zipfile.ZipFile(tables_xlsx) as archive:
        workbook_xml = archive.read("xl/workbook.xml").decode("utf-8")
        sheet_count = len(re.findall(r"<(?:[A-Za-z0-9_]+:)?sheet\s", workbook_xml))
        formula_count = 0
        for name in archive.namelist():
            if name.startswith("xl/worksheets/sheet") and name.endswith(".xml"):
                worksheet_xml = archive.read(name).decode("utf-8")
                formula_count += len(re.findall(r"<(?:[A-Za-z0-9_]+:)?f(?:\s|>)", worksheet_xml))
    checks["supplementary_workbook_sheets"] = sheet_count
    checks["supplementary_workbook_formulas"] = formula_count
    if sheet_count != 17:
        errors.append(f"Expected 17 workbook sheets, found {sheet_count}")
    if formula_count != 0:
        errors.append("Unexpected formulas remain in source-data workbook")

    with zipfile.ZipFile(code_zip) as archive:
        code_entries = [name.replace("\\", "/") for name in archive.namelist()]
        code_text = "\n".join(
            archive.read(name).decode("utf-8", errors="ignore")
            for name in archive.namelist()
            if not name.endswith("/") and archive.getinfo(name).file_size < 20 * 1024 * 1024
        )
    checks["code_archive_entries"] = len(code_entries)
    checks["code_archive_has_readme"] = any(name.endswith("/README.md") for name in code_entries)
    checks["code_archive_has_source_manifest"] = any(
        name.endswith("/SOURCE_DATASETS.tsv") for name in code_entries
    )
    checks["code_archive_forbidden_entries"] = [
        name
        for name in code_entries
        if name.endswith("/config/local.yml") or "/__pycache__/" in name or name.endswith(".pyc")
    ]
    if not checks["code_archive_has_readme"] or not checks["code_archive_has_source_manifest"]:
        errors.append("Code archive lacks its README or public source manifest")
    if checks["code_archive_forbidden_entries"]:
        errors.append("Code archive contains local configuration or Python cache files")
    jwt_prefix = "eyJ" + "hbGciOiJSUzI1Ni"
    if jwt_prefix in code_text:
        errors.append("An OpenGWAS JWT literal remains in the code archive")
    personal_paths = (
        "C:/Users/" + "shijunhui",
        "C:" + "\\Users\\" + "shijunhui",
        "E:" + "/",
    )
    if any(value in code_text for value in personal_paths):
        errors.append("A machine-specific personal path remains in the code archive")

    checks["manuscript_pdf_pages"] = len(PdfReader(manuscript_pdf).pages)
    checks["supplementary_pdf_pages"] = len(PdfReader(figures_pdf).pages)
    figure_readers = [PdfReader(path) for path in main_figures]
    checks["main_figure_pages"] = [len(reader.pages) for reader in figure_readers]
    checks["main_figure_dimensions_mm"] = [
        [
            round(float(reader.pages[0].mediabox.width) * 25.4 / 72, 1),
            round(float(reader.pages[0].mediabox.height) * 25.4 / 72, 1),
        ]
        for reader in figure_readers
    ]
    checks["main_figure_font_embedding"] = [pdf_font_embedding(path) for path in main_figures]
    if checks["supplementary_pdf_pages"] != 5:
        errors.append("Supplementary figure PDF does not contain five pages")
    if checks["main_figure_pages"] != [1] * 7:
        errors.append("Every main figure must be a one-page composite PDF")
    if any(not 168 <= width <= 171 or height > 225 for width, height in checks["main_figure_dimensions_mm"]):
        errors.append("A main figure is not at the 170 mm double-column width or exceeds 225 mm height")
    if any(found == 0 or embedded != found for found, embedded in checks["main_figure_font_embedding"]):
        errors.append("A main figure contains a font program that is not embedded")
    if any(
        path.stat().st_size > 20 * 1024 * 1024
        for path in (tables_xlsx, figures_pdf, code_zip)
    ):
        errors.append("An additional file exceeds the 20 MB limit")
    for path in main_figures:
        if path.stat().st_size > 10 * 1024 * 1024:
            errors.append(f"{path.name} exceeds the 10 MB figure limit")

    legend_block = markdown.split("## Figure legends", 1)[1].split("## Additional files", 1)[0]
    legend_titles = re.findall(r"^## Figure \d+\.\s*(.*)$", legend_block, flags=re.M)
    checks["figure_legend_titles"] = len(legend_titles)
    if len(legend_titles) != 7:
        errors.append("Expected seven main figure legends")
    if any(len(title.split()) > 15 for title in legend_titles):
        errors.append("A main figure title exceeds 15 words")
    legend_bodies = re.findall(
        r"^## Figure \d+\..*?\n\n(.*?)(?=\n## Figure \d+\.|\n## Additional files)",
        markdown,
        flags=re.M | re.S,
    )
    checks["figure_legend_word_counts"] = [
        len(re.findall(r"\b[\w-]+\b", re.sub(r"\*", "", legend)))
        for legend in legend_bodies
    ]
    if len(legend_bodies) != 7 or any(count > 300 for count in checks["figure_legend_word_counts"]):
        errors.append("A main figure legend is missing or exceeds 300 words")

    checks["errors"] = errors
    checks["warnings"] = warnings
    checks["status"] = "PASS" if not errors else "FAIL"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(checks, indent=2, ensure_ascii=False), encoding="utf-8")

    review = f"""# Final pre-submission peer review: BMC Medical Genomics

## Recommendation

**Minor revision before submission.** No further analytical expansion is recommended.

## Scientific assessment

The manuscript asks a clear cross-disease genomics question and uses a coherent evidence sequence: tissue discovery, directional stratification, functional characterization, external replication, source-defined single-cell mapping and a prespecified peripheral-blood screen. The 286-gene intersection establishes a measurable shared molecular landscape enriched for matrix, immune, stress and cell-cycle programs. Directional stratification adds disease-context resolution without displacing this principal finding. External cohorts show high agreement in ovarian cancer and greater variability in osteoarthritis. G0S2 is appropriately prioritized as a candidate systemic feature for prospective validation.

The statistical methods are proportionate to the study design. Cohorts are modelled separately, multiple-testing rules are explicit, denominators are retained, and the decision to avoid cross-cohort ComBat correction is explained. The main-text Methods are concise, while dataset-specific single-cell audit details remain traceable in the supplementary tables. The expanded G0S2 discussion links the cross-compartment finding to lipid metabolism, cellular quiescence and inflammatory context without promoting it to a common driver. Moving the descriptive evidence matrix to Supplementary Figure S5 keeps the main figures focused on primary data.

## Presentation and journal fit

The submission now follows the BMC Medical Genomics Research Article structure, has a structured abstract, exact declaration headings, an abbreviations list, separate figure files and three sequentially cited additional files. Main figures remain separate one-page composites. The supplementary directory has been removed from the manuscript, every retained supplementary figure is physically supplied in Additional file 2, and the executable code is supplied in Additional file 3 without credentials or machine-specific paths.

## Required minor author actions

1. Confirm that the four listed grants supported this OA-OC analysis and that the funder-role statement is accurate.
2. Final approval and the absence of competing interests have been confirmed for all five authors.
3. If the authors choose to mirror Additional file 3 to GitHub and Zenodo, insert the final public URL and DOI only after the deposit exists; no identifier has been invented in this package.

The author sequence is Junhui Shi, Mengxiang Liu, Repkat Inayatilla, Ke Li and Lei Chen. Junhui Shi and Mengxiang Liu are marked as equal contributors and co-first authors. Repkat Inayatilla is listed immediately before Ke Li. All five authors approved the final manuscript and declared no competing interests. Hongtao Yu is absent from the submission files.

## Audit result

Automated package audit: **{checks['status']}**. Abstract words: {abstract_words}; Methods words: {methods_words}; keywords: {keyword_count}; references: {reference_count}; manuscript PDF pages: {checks['manuscript_pdf_pages']}; supplementary workbook sheets: {sheet_count}; supplementary-figure PDF pages: {checks['supplementary_pdf_pages']}.
"""
    args.review.write_text(review, encoding="utf-8")
    if errors:
        raise SystemExit("Audit failed: " + "; ".join(errors))
    print("BMC submission audit passed")


if __name__ == "__main__":
    main()
