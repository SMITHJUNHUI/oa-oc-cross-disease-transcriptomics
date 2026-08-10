from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


# Low-contribution or redundant references removed during strategic compression.
# The remaining set retains disease background, primary datasets, core methods,
# single-cell standards, blood-context evidence and G0S2 biology.
REMOVED_V40_NUMBERS = {34, 36, 37, 38, 39, 44, 45, 46, 50, 51}


def parse_reference_lines(text: str) -> dict[int, str]:
    records: dict[int, str] = {}
    for line in text.splitlines():
        match = re.match(r"^(\d+)\.\s+(.+)$", line.strip())
        if match:
            records[int(match.group(1))] = match.group(2)
    return records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    source_md = args.source_dir / "references_v40.md"
    source_json = args.source_dir / "reference_audit_v40.json"
    references = parse_reference_lines(source_md.read_text(encoding="utf-8"))
    audit = json.loads(source_json.read_text(encoding="utf-8"))
    audit_by_number = {int(item["number"]): item for item in audit}

    expected = set(range(1, 81))
    if set(references) != expected or set(audit_by_number) != expected:
        raise RuntimeError("The V4.0 reference packet is incomplete.")

    retained = [number for number in range(1, 81) if number not in REMOVED_V40_NUMBERS]
    if len(retained) != 70:
        raise RuntimeError(f"Expected 70 retained references, found {len(retained)}")

    output_lines: list[str] = []
    output_audit: list[dict] = []
    for new_number, old_number in enumerate(retained, start=1):
        output_lines.append(f"{new_number}. {references[old_number]}")
        record = dict(audit_by_number[old_number])
        record["original_v40_number"] = old_number
        record["number"] = new_number
        output_audit.append(record)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "references_v41.md").write_text(
        "\n".join(output_lines) + "\n", encoding="utf-8"
    )
    (args.output_dir / "reference_audit_v41.json").write_text(
        json.dumps(output_audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    (args.output_dir / "reference_pruning_v41.json").write_text(
        json.dumps(
            {
                "source_count": 80,
                "retained_count": 70,
                "removed_v40_numbers": sorted(REMOVED_V40_NUMBERS),
                "principle": "Remove redundant network-medicine and single-cell method citations without weakening the evidence chain.",
            },
            indent=2,
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"Wrote 70-reference V4.1 packet to {args.output_dir}")


if __name__ == "__main__":
    main()
