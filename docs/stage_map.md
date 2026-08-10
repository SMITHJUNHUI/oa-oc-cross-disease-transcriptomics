# Pipeline stage map

| Stage | Inputs | Primary outputs | Cache |
|---|---|---|---|
| preflight | YAML, files, installed packages | manifests | always rerun |
| bulk_training | OA/OC training cohorts | harmonized training matrices, QC | yes |
| differential | OA/OC training matrices | full DEG tables, volcano plots | yes |
| shared | OA/OC DEG tables | shared gene evidence table | yes |
| enrichment | shared genes, GMT files | GO/KEGG/GSEA | yes |
| wgcna | OA/OC training matrices | module-trait tables and module genes | yes |
| machine_learning | shared/WGCNA candidates | LASSO, RF, final genes | yes |
| bulk_validation | four external cohorts and platform annotations | harmonized validation matrices, QC | yes |
| validation | four external cohorts | ROC/AUC tables and figures | yes |
| immune | OA/OC training matrices, signatures | scores, tests and figures | yes |
| tcga | TCGA-OV expression/clinical data | Cox, LASSO-Cox, KM, timeROC | yes |
| regulatory | final genes, miRTarBase/KnockTF | supplementary networks | yes |
| drug | final genes, DGIdb/CTD | supplementary candidate tables | yes |
| mr | pre-specified GWAS IDs and JWT | harmonized and MR estimates | yes |
| single_cell | five dataset-specific local inputs | adapter audits, per-batch QC, scDblFinder when supported, capability/blocker reports | yes |
| single_cell_downstream | five validated single-cell gates, published labels, hub genes | data-set-level embeddings/labels, composition, hub localization, replicate-aware pseudobulk | yes |
| report | all prior stage summaries | Markdown report and session info | always rerun |
