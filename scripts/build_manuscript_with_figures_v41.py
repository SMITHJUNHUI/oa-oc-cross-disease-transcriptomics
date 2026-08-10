from __future__ import annotations

import build_manuscript_with_figures as builder


builder.FIGURES = [
    ("Figure 1", "Figure1_study_design.png"),
    ("Figure 2", "Figure2_shared_tissue_transcriptomics.png"),
    ("Figure 3", "Figure3_shared_biological_themes.png"),
    ("Figure 4", "Figure4_external_tissue_and_illustrative_genes.png"),
    ("Figure 5", "Figure5_single_cell_localization.png"),
    ("Figure 6", "Figure6_peripheral_blood_validation.png"),
    ("Figure 7", "Figure7_integrated_interpretation.png"),
    ("Supplementary Figure 1", "SupplementaryFigure1_core_sensitivity.png"),
    ("Supplementary Figure 2", "SupplementaryFigure2_external_tissue_Hallmark.png"),
    ("Supplementary Figure 3", "SupplementaryFigure3_all_single_cell_UMAPs.png"),
    ("Supplementary Figure 4", "SupplementaryFigure4_bulk_PCA_QC.png"),
    ("Supplementary Figure 5", "SupplementaryFigure5_direction_aware_STRING_network.png"),
]


if __name__ == "__main__":
    builder.main()
