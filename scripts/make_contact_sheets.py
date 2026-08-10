from __future__ import annotations

import argparse
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--pages-per-sheet", type=int, default=4)
    args = parser.parse_args()
    pages = sorted(args.input_dir.glob("page-*.png"), key=lambda p: int(p.stem.split("-")[-1]))
    if not pages:
        raise SystemExit("No rendered page images found")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    font = ImageFont.load_default()
    for offset in range(0, len(pages), args.pages_per_sheet):
        batch = pages[offset:offset + args.pages_per_sheet]
        thumbs: list[tuple[Path, Image.Image]] = []
        for path in batch:
            image = Image.open(path).convert("RGB")
            image.thumbnail((620, 880), Image.Resampling.LANCZOS)
            thumbs.append((path, image.copy()))
        sheet = Image.new("RGB", (1280, 1820), "#D0D5DD")
        draw = ImageDraw.Draw(sheet)
        for index, (path, image) in enumerate(thumbs):
            row, col = divmod(index, 2)
            x = 10 + col * 635
            y = 35 + row * 900
            sheet.paste(image, (x, y))
            page_number = int(path.stem.split("-")[-1])
            draw.text((x, 8 + row * 900), f"Page {page_number}", fill="black", font=font)
        first = int(batch[0].stem.split("-")[-1])
        last = int(batch[-1].stem.split("-")[-1])
        sheet.save(args.output_dir / f"contact_{first:02d}_{last:02d}.png")


if __name__ == "__main__":
    main()
