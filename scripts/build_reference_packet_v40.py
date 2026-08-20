from __future__ import annotations

import argparse
import html
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ReferenceSeed:
    category: str
    doi: str | None = None
    manual: str | None = None


SEEDS = [
    # Osteoarthritis background, tissue transcriptomics and single-cell biology.
    ReferenceSeed(
        "OA",
        manual=(
            "Glyn-Jones S, Palmer AJR, Agricola R, et al. Osteoarthritis. "
            "Lancet. 2015;386:376-387. doi:10.1016/S0140-6736(14)60802-3."
        ),
    ),
    ReferenceSeed("OA", "10.1016/S0140-6736(19)30417-9"),
    ReferenceSeed("OA", "10.1038/nrdp.2016.72"),
    ReferenceSeed("OA", "10.1002/art.34453"),
    ReferenceSeed("OA", "10.1016/S0140-6736(11)60243-2"),
    ReferenceSeed("OA", "10.1002/jcp.21258"),
    ReferenceSeed("OA", "10.1016/j.joca.2015.03.036"),
    ReferenceSeed("OA", "10.1016/j.joca.2018.07.012"),
    ReferenceSeed("OA", "10.1136/annrheumdis-2017-212863"),
    ReferenceSeed(
        "OA",
        manual=(
            "Fu W, Hettinghouse A, Chen Y, et al. 14-3-3 epsilon is an intracellular "
            "component of the TNFR2 receptor complex and its activation protects against "
            "osteoarthritis. Ann Rheum Dis. 2021;80:1615-1627. "
            "doi:10.1136/annrheumdis-2021-220000."
        ),
    ),
    ReferenceSeed("OA", "10.1136/ard-2023-224420"),
    ReferenceSeed("OA", "10.1016/j.gendis.2024.101241"),
    ReferenceSeed("OA", "10.1038/s41598-020-67730-y"),
    ReferenceSeed("OA", "10.1016/j.joca.2024.10.008"),
    ReferenceSeed("OA", "10.1136/annrheumdis-2013-203405"),
    # Ovarian cancer background, bulk molecular classes and cellular ecosystems.
    ReferenceSeed("OC", "10.1016/S0140-6736(18)32552-2"),
    ReferenceSeed("OC", "10.1038/nrc4019"),
    ReferenceSeed("OC", "10.1038/nrdp.2016.61"),
    ReferenceSeed(
        "OC",
        manual=(
            "Cancer Genome Atlas Research Network. Integrated genomic analyses of "
            "ovarian carcinoma. Nature. 2011;474:609-615. doi:10.1038/nature10166."
        ),
    ),
    ReferenceSeed("OC", "10.1038/nature14410"),
    ReferenceSeed("OC", "10.1158/1078-0432.CCR-08-0196"),
    ReferenceSeed("OC", "10.1093/jnci/dju249"),
    ReferenceSeed("OC", "10.1016/j.ccr.2009.10.018"),
    ReferenceSeed("OC", "10.1111/j.1349-7006.2009.01204.x"),
    ReferenceSeed("OC", "10.3390/ijms20040952"),
    ReferenceSeed("OC", "10.1038/s41591-020-0926-0"),
    ReferenceSeed("OC", "10.1158/1078-0432.CCR-22-0296"),
    ReferenceSeed("OC", "10.1016/j.ccell.2021.04.004"),
    ReferenceSeed("OC", "10.1158/0008-5472.CAN-20-0521"),
    ReferenceSeed("OC", "10.1038/s41586-022-05496-1"),
    ReferenceSeed("OC", "10.1186/1471-2407-13-178"),
    # Cross-disease and systems-level interpretation.
    ReferenceSeed(
        "Cross-disease",
        manual=(
            "Dudley JT, Tibshirani R, Deshpande T, Butte AJ. Disease signatures are "
            "robust across tissues and experiments. Mol Syst Biol. 2009;5:307. "
            "doi:10.1038/msb.2009.66."
        ),
    ),
    ReferenceSeed("Cross-disease", "10.1126/science.1257601"),
    ReferenceSeed("Cross-disease", "10.1038/nrg2918"),
    ReferenceSeed("Cross-disease", "10.1186/s13059-017-1215-1"),
    ReferenceSeed("Cross-disease", "10.1016/j.cell.2017.10.049"),
    ReferenceSeed(
        "Cross-disease",
        manual=(
            "Pinero J, Ramirez-Anguita JM, Sauch-Pitarch J, et al. The DisGeNET "
            "knowledge platform for disease genomics: 2019 update. Nucleic Acids Res. "
            "2020;48(D1):D845-D855. doi:10.1093/nar/gkz1021."
        ),
    ),
    ReferenceSeed("Cross-disease", "10.1371/journal.pone.0006536"),
    ReferenceSeed("Cross-disease", "10.1038/tp.2016.87"),
    # Single-cell analysis and context-aware interpretation.
    ReferenceSeed("Single-cell", "10.1038/s41576-019-0093-7"),
    ReferenceSeed("Single-cell", "10.15252/msb.20188746"),
    ReferenceSeed("Single-cell", "10.1038/s41592-019-0654-x"),
    ReferenceSeed("Single-cell", "10.1016/j.cell.2021.04.048"),
    ReferenceSeed("Single-cell", "10.1186/s13059-017-1382-0"),
    ReferenceSeed("Single-cell", "10.1038/nbt.4096"),
    ReferenceSeed("Single-cell", "10.1038/s41592-019-0619-0"),
    ReferenceSeed("Single-cell", "10.1186/s13059-019-1874-1"),
    ReferenceSeed("Single-cell", "10.1016/j.cels.2019.03.003"),
    ReferenceSeed("Single-cell", "10.1016/j.cels.2018.11.005"),
    ReferenceSeed("Single-cell", "10.1038/s41590-018-0276-y"),
    ReferenceSeed("Single-cell", "10.1126/science.aad0501"),
    ReferenceSeed("Single-cell", "10.12688/f1000research.73600.2"),
    ReferenceSeed("Single-cell", "10.1016/j.csbj.2021.06.043"),
    ReferenceSeed("Single-cell", "10.1002/ctm2.500"),
    ReferenceSeed("Single-cell", "10.1101/cshperspect.a041314"),
    # Peripheral blood interpretation and the G0S2 candidate signal.
    ReferenceSeed("Blood/G0S2", "10.1038/nature09247"),
    ReferenceSeed("Blood/G0S2", "10.1371/journal.pcbi.1002240"),
    ReferenceSeed("Blood/G0S2", "10.1093/jnci/93.14.1054"),
    ReferenceSeed("Blood/G0S2", "10.1038/s41591-020-0752-4"),
    ReferenceSeed("Blood/G0S2", "10.1073/pnas.252784499"),
    ReferenceSeed("Blood/G0S2", "10.1158/0008-5472.CAN-15-2265"),
    ReferenceSeed("Blood/G0S2", "10.1016/j.cmet.2010.02.003"),
    ReferenceSeed("Blood/G0S2", "10.1016/j.bbalip.2017.06.007"),
    ReferenceSeed("Blood/G0S2", "10.1158/1078-0432.CCR-18-2693"),
    # Statistical, enrichment, network and reproducibility methods.
    ReferenceSeed("Methods", "10.1093/nar/gkv007"),
    ReferenceSeed("Methods", "10.2202/1544-6115.1027"),
    ReferenceSeed(
        "Methods",
        manual=(
            "Benjamini Y, Hochberg Y. Controlling the false discovery rate: "
            "a practical and powerful approach to multiple testing. J R Stat Soc B. "
            "1995;57:289-300."
        ),
    ),
    ReferenceSeed("Methods", "10.1073/pnas.0506580102"),
    ReferenceSeed("Methods", "10.1016/j.cels.2015.12.004"),
    ReferenceSeed("Methods", "10.1186/1471-2105-9-559"),
    ReferenceSeed("Methods", "10.1093/nar/gkac1000"),
    ReferenceSeed("Methods", "10.1038/sdata.2016.18"),
    ReferenceSeed("Methods", "10.1038/nrg2825"),
    ReferenceSeed("Methods", "10.1089/omi.2011.0118"),
    ReferenceSeed("Methods", "10.1016/j.xinn.2021.100141"),
    ReferenceSeed("Methods", "10.1038/75556"),
    ReferenceSeed("Methods", "10.1093/nar/gkaa970"),
    ReferenceSeed("Methods", "10.1093/nar/30.1.207"),
    ReferenceSeed(
        "Methods",
        manual=(
            "Barrett T, Wilhite SE, Ledoux P, et al. NCBI GEO: archive for functional "
            "genomics data sets-update. Nucleic Acids Res. 2013;41(D1):D991-D995. "
            "doi:10.1093/nar/gks1193."
        ),
    ),
    ReferenceSeed("Methods", "10.1186/gb-2004-5-10-r80"),
]


def clean(text: str) -> str:
    text = html.unescape(text or "")
    text = re.sub(r"<[^>]+>", "", text)
    text = text.translate(
        str.maketrans(
            {
                "\u2010": "-",
                "\u2011": "-",
                "\u2012": "-",
                "\u2013": "-",
                "\u2014": "-",
                "\u2212": "-",
            }
        )
    )
    return re.sub(r"\s+", " ", text).strip()


def initials(given: str) -> str:
    parts = re.findall(r"[A-Za-zÀ-ÖØ-öø-ÿ]+", given or "")
    return "".join(part[0].upper() for part in parts if part)


def format_authors(authors: list[dict]) -> str:
    if not authors:
        return "Author information unavailable"
    formatted = []
    for author in authors[:6]:
        family = clean(author.get("family", ""))
        given = initials(author.get("given", ""))
        name = " ".join(piece for piece in (family, given) if piece)
        if name:
            formatted.append(name)
    suffix = ", et al." if len(authors) > 6 else "."
    return ", ".join(formatted) + suffix


def date_year(message: dict) -> str:
    for key in ("published-print", "published-online", "published", "issued"):
        try:
            return str(message[key]["date-parts"][0][0])
        except (KeyError, IndexError, TypeError):
            continue
    return "Year unavailable"


def get_crossref(doi: str) -> dict:
    encoded = urllib.parse.quote(doi, safe="")
    url = f"https://api.crossref.org/works/{encoded}?mailto=research@example.com"
    last_error: Exception | None = None
    for attempt in range(4):
        request = urllib.request.Request(
            url,
            headers={"User-Agent": "AcademicReferenceAudit/1.0 (mailto:research@example.com)"},
        )
        try:
            with urllib.request.urlopen(request, timeout=30) as response:
                payload = json.load(response)
            return payload["message"]
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as exc:
            last_error = exc
            if isinstance(exc, urllib.error.HTTPError) and exc.code == 404:
                raise
            time.sleep(1.0 + attempt * 1.5)
    raise last_error or RuntimeError(f"Crossref lookup failed: {doi}")


def format_reference(message: dict, doi: str) -> str:
    title = clean((message.get("title") or [""])[0])
    journal = clean((message.get("container-title") or [""])[0])
    volume = clean(str(message.get("volume", "")))
    issue = clean(str(message.get("issue", "")))
    pages = clean(str(message.get("page", "") or message.get("article-number", "")))
    year = date_year(message)
    author_text = format_authors(message.get("author") or [])
    journal_block = journal
    if volume:
        journal_block += f". {year};{volume}"
        if issue:
            journal_block += f"({issue})"
        if pages:
            journal_block += f":{pages}"
        journal_block += "."
    else:
        journal_block += f". {year}."
    return f"{author_text} {title}. {journal_block} doi:{doi}."


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    records: list[dict] = []
    lines: list[str] = []
    seen = set()
    indexed_seeds = list(enumerate(SEEDS, start=1))
    doi_tasks = []
    manual_by_index: dict[int, dict] = {}
    for index, seed in indexed_seeds:
        if seed.manual:
            manual_by_index[index] = {
                "number": index,
                "category": seed.category,
                "doi": None,
                "status": "manual_verified",
                "formatted": seed.manual,
            }
            continue
        doi = seed.doi.strip().lower()
        if doi in seen:
            raise RuntimeError(f"Duplicate DOI in reference seed list: {doi}")
        seen.add(doi)
        doi_tasks.append((index, seed.category, doi))

    fetched: dict[int, dict] = {}
    with ThreadPoolExecutor(max_workers=3) as pool:
        future_map = {
            pool.submit(get_crossref, doi): (index, category, doi)
            for index, category, doi in doi_tasks
        }
        for future in as_completed(future_map):
            index, category, doi = future_map[future]
            try:
                message = future.result()
                formatted = format_reference(message, doi)
                fetched[index] = {
                    "number": index,
                    "category": category,
                    "doi": doi,
                    "status": "crossref_verified",
                    "title": clean((message.get("title") or [""])[0]),
                    "journal": clean((message.get("container-title") or [""])[0]),
                    "year": date_year(message),
                    "formatted": formatted,
                }
            except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as exc:
                fetched[index] = {
                    "number": index,
                    "category": category,
                    "doi": doi,
                    "status": "needs_manual_check",
                    "error": str(exc),
                    "formatted": f"[REFERENCE METADATA NEEDS CHECK: doi:{doi}]",
                }

    for index, _seed in indexed_seeds:
        record = manual_by_index.get(index) or fetched[index]
        records.append(record)
        lines.append(f"{index}. {record['formatted']}")

    (args.output_dir / "reference_audit_v40.json").write_text(
        json.dumps(records, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    (args.output_dir / "references_v40.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    unresolved = [record for record in records if record["status"] == "needs_manual_check"]
    print(f"references={len(records)} unresolved={len(unresolved)}")
    if unresolved:
        for record in unresolved:
            print(f"UNRESOLVED {record['number']}: {record['doi']}")
        raise SystemExit(2)


if __name__ == "__main__":
    main()
