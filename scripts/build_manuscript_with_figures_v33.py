from __future__ import annotations

import build_manuscript_with_figures as builder


builder.FIGURES = [
    ("Figure 1", "Figure1_study_workflow_with_blood.png"),
    ("Figure 2", "Figure2_shared_transcriptomic_alterations.png"),
    ("Figure 3", "Figure3_functional_characterization.png"),
    ("Figure 4", "Figure4_transcriptional_divergence.png"),
    ("Figure 5", "Figure5_peripheral_blood_systemic_component.png"),
    ("Figure 6", "Figure6_single_cell_localization.png"),
    ("Figure 7", "Figure7_integrated_summary_model.png"),
    ("Supplementary Figure 1", "SupplementaryFigure1_core_sensitivity.png"),
    ("Supplementary Figure 2", "SupplementaryFigure2_complete_Hallmark_direction.png"),
    ("Supplementary Figure 3", "SupplementaryFigure3_all_single_cell_UMAPs.png"),
    ("Supplementary Figure 4", "SupplementaryFigure4_bulk_PCA_QC.png"),
    ("Supplementary Figure 5", "SupplementaryFigure5_direction_aware_STRING_network.png"),
]


if __name__ == "__main__":
    builder.main()
