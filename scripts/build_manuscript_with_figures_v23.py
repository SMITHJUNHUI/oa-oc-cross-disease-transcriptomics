from __future__ import annotations

import build_manuscript_with_figures as builder


builder.FIGURES = [
    *builder.FIGURES,
    (
        "Supplementary Figure 6",
        "SupplementaryFigure6_cell_type_functional_annotation.png",
    ),
    (
        "Supplementary Figure 7",
        "SupplementaryFigure7_pathway_direction.png",
    ),
    (
        "Supplementary Figure 8",
        "SupplementaryFigure8_detailed_feature_stability.png",
    ),
]


if __name__ == "__main__":
    builder.main()
