from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


REMOVED_V41_NUMBERS = {60}


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

    references = parse_reference_lines((args.source_dir / "references_v41.md").read_text(encoding="utf-8"))
    audit = json.loads((args.source_dir / "reference_audit_v41.json").read_text(encoding="utf-8"))
    audit_by_number = {int(item["number"]): item for item in audit}
    expected = set(range(1, 71))
    if set(references) != expected or set(audit_by_number) != expected:
        raise RuntimeError("The V4.1 reference packet is incomplete.")

    retained = [number for number in range(1, 71) if number not in REMOVED_V41_NUMBERS]
    output_lines: list[str] = []
    output_audit: list[dict] = []
    for new_number, old_number in enumerate(retained, start=1):
        output_lines.append(f"{new_number}. {references[old_number]}")
        record = dict(audit_by_number[old_number])
        record["original_v41_number"] = old_number
        record["number"] = new_number
        output_audit.append(record)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "references_v42.md").write_text("\n".join(output_lines) + "\n", encoding="utf-8")
    (args.output_dir / "reference_audit_v42.json").write_text(
        json.dumps(output_audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    (args.output_dir / "reference_pruning_v42.json").write_text(
        json.dumps(
            {
                "source_count": 70,
                "retained_count": 69,
                "removed_v41_numbers": [60],
                "principle": "Remove the WGCNA citation after excluding WGCNA from the V4.2 submission scope.",
            },
            indent=2,
            ensure_ascii=False,
        ) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote 69-reference V4.2 packet to {args.output_dir}")


if __name__ == "__main__":
    main()

