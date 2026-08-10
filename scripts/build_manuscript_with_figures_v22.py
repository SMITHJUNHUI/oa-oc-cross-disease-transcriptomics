from __future__ import annotations

import build_manuscript_with_figures as builder


builder.FIGURES = [
    *builder.FIGURES,
    (
        "Supplementary Figure 6",
        "SupplementaryFigure6_cell_type_functional_annotation.png",
    ),
]


if __name__ == "__main__":
    builder.main()
