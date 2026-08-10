from __future__ import annotations

import argparse
from pathlib import Path

import pypdfium2 as pdfium


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--dpi", type=int, default=130)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    document = pdfium.PdfDocument(str(args.input))
    scale = args.dpi / 72.0
    for zero_index in range(len(document)):
        page = document[zero_index]
        bitmap = page.render(scale=scale)
        image = bitmap.to_pil().convert("RGB")
        image.save(args.output_dir / f"page-{zero_index + 1}.png")
        bitmap.close()
        page.close()
    print(f"Rendered {len(document)} pages from {args.input}")


if __name__ == "__main__":
    main()
