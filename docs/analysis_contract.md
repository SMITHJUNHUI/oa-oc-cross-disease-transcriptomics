# Analysis contract and scientific boundaries

## Primary question

The project evaluates whether OA and OC show context-dependent transcriptional convergence and divergence. It does not test whether one disease clinically causes the other, and it does not generate treatment recommendations.

## Fixed primary settings

- Differential expression: BH FDR <0.05 and absolute log2 fold change >=1.0.
- Shared-gene directions are always retained; 146/286 are concordant and 140/286 discordant.
- WGCNA is performed separately by disease.
- Feature screening, scaling, and tuning remain inside resampling folds.
- The fixed ten-gene evidence set is SOX9, ELF3, JUNB, AKAP12, BNC1, CFI, DDIT3, DIRAS3, EFEMP1, and HK2.
- AUC is a secondary retrospective molecular-separability measure, not clinical performance.

## V3.0 network boundaries

- STRING uses all mapped shared DEGs with no added neighbors. Physical score >=0.700 is primary; the functional graph is sensitivity only.
- Mapped isolates remain in topology denominators.
- Direction-label permutation tests organization in a fixed database graph; it does not establish tissue-specific interaction or mechanism.

## V3.0 single-cell boundaries

- OA and OC atlases are never integrated into a common latent space.
- Exact source labels are preserved. `Fibroblast` is not relabeled as CAF, and `Ovarian.cancer.cell` is not relabeled as tumor epithelial.
- CCSS is a within-dataset detection-specificity score, not differential expression or pathway activity.
- UCell is summarized within biological sample first. Absolute scores are not compared across diseases.
- Pseudobulk inference requires biological replication and an interpretable contrast.

## Communication boundaries

- CellChat is fitted separately per biological sample, with deterministic balanced subsampling and explicit minimum cell counts.
- Consensus requires at least 3/5 OC samples, 2/3 OA-control donors, or 4/8 OA donors.
- CellChat probabilities are not compared numerically between OA and OC.
- NicheNet is restricted to an official v2 prior-consistency overlay. Ligand activity, differential regulation, mediation, spatial contact, signaling flux, and causality are not inferred.

## MR boundary

Bidirectional MR asks whether genetic liability to OA affects OC risk or genetic liability to OC affects OA risk under the selected GWAS datasets, instruments, and assumptions. It is not a test of shared heritability, common susceptibility loci, or genetic correlation.

## Completion criteria

- One-command V3.0 regression exits successfully.
- All 6 main and 14 supplementary figures have machine-readable source data and legends.
- Tables S1-S28, parameters, data sources, package versions, and cached external responses are retained.
- Claims are checked against the claim-evidence registry.
- No credential-like string is stored in the project or submission archive.
