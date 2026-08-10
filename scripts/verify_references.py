from __future__ import annotations

import argparse
import json
import re
import time
import unicodedata
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from difflib import SequenceMatcher
from pathlib import Path


REFERENCE_RE = re.compile(r"^(\d+)\.\s+(.*)$")
CITATION_RE = re.compile(r"\[([0-9,\-]+)\]")
DOI_RE = re.compile(r"doi:\s*(10\.\d{4,9}/\S+)", re.I)


def normalize(value: str) -> str:
    value = unicodedata.normalize("NFKD", value)
    value = "".join(char for char in value if not unicodedata.combining(char))
    return re.sub(r"[^a-z0-9]+", " ", value.lower()).strip()


def expand(content: str) -> list[int]:
    result: list[int] = []
    for part in content.split(","):
        if "-" in part:
            start, end = (int(value) for value in part.split("-", 1))
            result.extend(range(start, end + 1))
        else:
            result.append(int(part))
    return result


def parse_manuscript(path: Path) -> tuple[str, dict[int, str]]:
    text = path.read_text(encoding="utf-8")
    body, tail = text.split("## References", 1)
    reference_text = tail.split("## Figure legends", 1)[0]
    references: dict[int, str] = {}
    for line in reference_text.splitlines():
        match = REFERENCE_RE.match(line.strip())
        if match:
            references[int(match.group(1))] = match.group(2)
    return body, references


def crossref(doi: str, retries: int = 3) -> dict:
    encoded = urllib.parse.quote(doi, safe="")
    request = urllib.request.Request(
        f"https://api.crossref.org/works/{encoded}",
        headers={"User-Agent": "OC-OA-reference-audit/1.0 (mailto:564386249@qq.com)"},
    )
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(request, timeout=30) as response:
                return json.loads(response.read().decode("utf-8"))["message"]
        except (urllib.error.URLError, TimeoutError):
            if attempt + 1 == retries:
                raise
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError("unreachable")


def datacite(doi: str, retries: int = 3) -> dict:
    encoded = urllib.parse.quote(doi, safe="")
    request = urllib.request.Request(
        f"https://api.datacite.org/dois/{encoded}",
        headers={"User-Agent": "OC-OA-reference-audit/1.0 (mailto:564386249@qq.com)"},
    )
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(request, timeout=30) as response:
                return json.loads(response.read().decode("utf-8"))["data"]["attributes"]
        except (urllib.error.URLError, TimeoutError):
            if attempt + 1 == retries:
                raise
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError("unreachable")


def verify_item(pair: tuple[int, str]) -> dict:
    number, record = pair
    doi_match = DOI_RE.search(record)
    doi = doi_match.group(1).rstrip(".") if doi_match else None
    item = {
        "number": number,
        "record": record,
        "doi": doi,
        "status": "unverifiable",
        "title_similarity": None,
        "crossref_title": None,
        "crossref_first_author": None,
        "crossref_year": None,
        "crossref_volume": None,
        "crossref_issue": None,
        "crossref_pages": None,
        "crossref_journal": None,
        "crossref_journal_short": None,
        "crossref_authors": [],
        "error": None,
    }
    if not doi:
        return item
    try:
        try:
            metadata = crossref(doi)
            source = "crossref"
        except urllib.error.HTTPError as exc:
            if exc.code != 404:
                raise
            metadata = datacite(doi)
            source = "datacite"
        if source == "datacite":
            title = ((metadata.get("titles") or [{}])[0].get("title") or "")
            journal = metadata.get("publisher") or ""
            creators = metadata.get("creators") or []
            first_author = (creators[0].get("familyName") or creators[0].get("name") or "") if creators else ""
            year = metadata.get("publicationYear")
            similarity = SequenceMatcher(None, normalize(title), normalize(record)).ratio()
            title_present = bool(normalize(title)) and normalize(title) in normalize(record)
            item.update(
                {
                    "status": "verified" if title_present or similarity >= 0.35 else "needs_fix",
                    "title_similarity": round(similarity, 3),
                    "crossref_title": title,
                    "crossref_first_author": first_author,
                    "crossref_year": year,
                    "crossref_journal": journal,
                    "crossref_journal_short": journal,
                    "crossref_authors": [
                        {
                            "given": creator.get("givenName", ""),
                            "family": creator.get("familyName", ""),
                            "name": creator.get("name", ""),
                        }
                        for creator in creators
                    ],
                }
            )
            return item
        title = (metadata.get("title") or [""])[0]
        journal = (metadata.get("container-title") or [""])[0]
        authors = metadata.get("author") or []
        first_author = authors[0].get("family", "") if authors else ""
        dates = metadata.get("published-print") or metadata.get("published") or metadata.get("published-online") or {}
        date_parts = dates.get("date-parts") or [[]]
        year = date_parts[0][0] if date_parts and date_parts[0] else None
        similarity = SequenceMatcher(None, normalize(title), normalize(record)).ratio()
        title_present = bool(normalize(title)) and normalize(title) in normalize(record)
        item.update(
            {
                "status": "verified" if title_present or similarity >= 0.35 else "needs_fix",
                "title_similarity": round(similarity, 3),
                "crossref_title": title,
                "crossref_first_author": first_author,
                "crossref_year": year,
                "crossref_volume": metadata.get("volume"),
                "crossref_issue": metadata.get("issue"),
                "crossref_pages": metadata.get("page") or metadata.get("article-number"),
                "crossref_journal": journal,
                "crossref_journal_short": (metadata.get("short-container-title") or [""])[0],
                "crossref_authors": [
                    {
                        "given": author.get("given", ""),
                        "family": author.get("family", ""),
                        "name": author.get("name", ""),
                    }
                    for author in authors
                ],
            }
        )
    except Exception as exc:  # network/HTTP details retained in the report
        item["error"] = f"{type(exc).__name__}: {exc}"
    return item


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manuscript", required=True, type=Path)
    parser.add_argument("--json", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()

    body, references = parse_manuscript(args.manuscript)
    first_order: list[int] = []
    all_citations: list[int] = []
    for match in CITATION_RE.finditer(body):
        for number in expand(match.group(1)):
            all_citations.append(number)
            if number not in first_order:
                first_order.append(number)

    with ThreadPoolExecutor(max_workers=6) as executor:
        results = list(executor.map(verify_item, sorted(references.items())))

    payload = {
        "reference_count": len(references),
        "cited_reference_count": len(set(all_citations)),
        "uncited_references": sorted(set(references) - set(all_citations)),
        "undefined_citations": sorted(set(all_citations) - set(references)),
        "first_appearance_order": first_order,
        "strict_first_appearance_order": first_order == list(range(1, len(references) + 1)),
        "results": results,
    }
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")

    counts = {key: sum(item["status"] == key for item in results) for key in ("verified", "needs_fix", "unverifiable")}
    lines = [
        "# Reference verification report",
        "",
        f"- References: {len(references)}",
        f"- DOI registry verified (Crossref/DataCite): {counts['verified']}",
        f"- Potential metadata mismatch: {counts['needs_fix']}",
        f"- Unverifiable/no DOI: {counts['unverifiable']}",
        f"- All references cited: {not payload['uncited_references']}",
        f"- All citations defined: {not payload['undefined_citations']}",
        f"- Numbered in first-appearance order: {payload['strict_first_appearance_order']}",
        "",
        "## Items requiring review",
        "",
        "| Current no. | Status | DOI | Crossref title / error |",
        "|---:|---|---|---|",
    ]
    for item in results:
        if item["status"] != "verified":
            detail = item["crossref_title"] or item["error"] or "No DOI in manuscript"
            lines.append(f"| {item['number']} | {item['status']} | {item['doi'] or ''} | {detail} |")
    lines.extend(["", "## First-appearance order", "", ", ".join(map(str, first_order)), ""])
    args.report.write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
