from __future__ import annotations

import build_manuscript_with_figures as builder


builder.FIGURES = [
    ("Figure 1", "Figure1_question_driven_framework.png"),
    ("Figure 2", "Figure2_gene_and_direction_heterogeneity.png"),
    ("Figure 3", "Figure3_pathway_direction.png"),
    ("Figure 4", "Figure4_cellular_context.png"),
    ("Figure 5", "Figure5_molecular_separability.png"),
    ("Figure 6", "Figure6_integrated_context_model.png"),
    ("Supplementary Figure 1", "SupplementaryFigure1_sensitivity_details.png"),
    ("Supplementary Figure 2", "SupplementaryFigure2_negative_bidirectional_MR.png"),
    ("Supplementary Figure 3", "SupplementaryFigure3_single_cell_UMAPs.png"),
    ("Supplementary Figure 4", "SupplementaryFigure4_HPA_normal_tissue_context.png"),
    ("Supplementary Figure 5", "SupplementaryFigure5_TCGA_relative_context.png"),
    ("Supplementary Figure 6", "SupplementaryFigure6_cell_type_functional_annotation.png"),
    ("Supplementary Figure 7", "SupplementaryFigure7_complete_pathway_direction.png"),
    ("Supplementary Figure 8", "SupplementaryFigure8_candidate_evidence_stability.png"),
    ("Supplementary Figure 9", "SupplementaryFigure9_candidate_centered_Hallmark_context.png"),
    ("Supplementary Figure 10", "SupplementaryFigure10_design_and_reliability.png"),
    ("Supplementary Figure 11", "SupplementaryFigure11_upstream_regulatory_context.png"),
    ("Supplementary Figure 12", "SupplementaryFigure12_direction_aware_STRING_network.png"),
    ("Supplementary Figure 13", "SupplementaryFigure13_sample_consensus_CellChat.png"),
    ("Supplementary Figure 14", "SupplementaryFigure14_NicheNet_prior_overlay.png"),
    ("Supplementary Figure 15", "SupplementaryFigure15_panel_size_sensitivity.png"),
    ("Supplementary Figure 16", "SupplementaryFigure16_bulk_PCA_and_QC.png"),
]


if __name__ == "__main__":
    builder.main()
