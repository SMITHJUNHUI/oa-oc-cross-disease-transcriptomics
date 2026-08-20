from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import build_manuscript_v24 as v24


TITLE = "# A systems-level analysis of context-dependent transcriptional convergence and divergence between osteoarthritis and ovarian cancer"


ABSTRACT = """## Abstract

### Background

Osteoarthritis (OA) and ovarian cancer (OC) arise in fundamentally different tissues. We tested whether their transcriptomes show context-dependent molecular convergence and divergence across gene, pathway, protein-association, cellular, intercellular-communication, and inherited-liability layers.

### Methods

Public human bulk and single-cell transcriptomic datasets were analyzed in a reproducible workflow. Differential expression was estimated separately in OA cartilage and ovarian tissue, followed by direction-aware overlap, disease-specific co-expression analysis, transparent ten-gene prioritization, strict nested resampling, and direction-fixed external molecular-separability assessment. STRING v12.0 supplied bounded protein-association context. Five single-cell atlases were evaluated with sample-aware UCell scores, a detection-based Candidate Cell Context Specificity Score (CCSS), and eligible pseudobulk contrasts. CellChat v2 was fitted separately to biological samples in two eligible atlases and aggregated by sample consensus. Consensus ligands were projected onto the official NicheNet v2 human ligand-target prior without ligand-activity inference. Paired Hallmark direction, HPA, TCGA-OV, focused regulatory resources, and supplementary bidirectional Mendelian randomization (MR) constrained interpretation.

### Results

We identified 2,008 OA and 2,310 OC differentially expressed genes, with 286 shared: 146 (51.0%) were concordant and 140 (49.0%) discordant. Among 10 Hallmark sets significant in both diseases, 6 had opposite directions. STRING mapped 275 shared genes; the high-confidence physical network contained 62 edges among 46 connected products, and concordant labels were more densely connected than expected under fixed-size label permutation. Across 1,025,361 quality-control-pass single cells, CCSS and sample-aware UCell localized the ten-gene evidence set to distinct dataset-specific cartilage, ovarian-cancer-cell, stromal, and immune contexts. Sample-resolved CellChat yielded 49,647 significant sample-level interactions and 199 consensus context-pathway records; 526 consensus interactions overlapped a shared DEG, but none contained a fixed candidate as a direct ligand or receptor. NicheNet supplied 10,240 prior ligand-target weights involving all ten candidates, but ligand activity was not estimable under the available design. External separability remained modest in OA and strong in retrospective OC tumor-reference contrasts. Bidirectional MR found no evidence of disease-to-disease genetic-liability effects under the selected datasets, instruments, and assumptions.

### Conclusions

OA and OC share a partial transcriptional landscape, but shared membership does not imply shared direction, network architecture, pathway state, cellular meaning, intercellular signaling, or inherited causality. The study defines an auditable systems-level resource and testable contextual hypotheses rather than a shared mechanism, diagnostic signature, therapeutic target set, or causal OA-OC relationship.

**Keywords:** osteoarthritis; ovarian cancer; cross-disease transcriptomics; single-cell RNA sequencing; STRING; UCell; CellChat; NicheNet; molecular separability
"""


STUDY_DESIGN = """### Study design and reproducibility

OA and OC were analyzed separately at every disease-specific stage and connected only by explicit gene-, pathway-, and evidence-level comparisons. The evidence chain comprised bulk discovery, gene direction, bounded STRING interaction context, paired pathway direction, disease-specific WGCNA, transparent candidate evidence, strict nested feature stability, cross-cohort molecular separability, five-atlas single-cell localization, sample-consensus CellChat in eligible atlases, a bounded NicheNet prior overlay, HPA/TCGA context audits, focused upstream resources, and supplementary bidirectional MR (Figure 1). Each layer answered a distinct question and carried an explicit inference boundary. Fixed configurations, manifests, seeds, cached external responses, tests, result tables, exact figure source data, and a claim-evidence registry are stored in the one-command project; randomized analyses used seed 20260726.
"""


PPI_METHOD = """### Direction-aware protein-association context

The 286 shared DEGs were submitted without additional neighbors to STRING v12.0 for *Homo sapiens* identifier mapping and interaction retrieval [36]. The primary graph retained physical STRING associations with score >=0.700; a high-confidence functional-association graph was retained as sensitivity. Mapped isolates remained in all denominators. Degree, normalized betweenness, density, component count, and largest-component size were descriptive, and network topology was not used to choose the ten genes. Concordant and discordant induced subgraphs were compared with 10,000 fixed-size permutations of direction labels on the fixed mapped graph. These analyses test organization within a database-derived association graph; they do not demonstrate physical binding, pathway activity, or tissue-specific interaction in OA or OC.
"""


SINGLE_CELL_METHOD = """### Single-cell localization, cell-context specificity, and UCell scoring

OA GSE104782, GSE169454, and GSE255460 and OC GSE154600 and GSE180661 were processed with dataset-specific adapters [17-21]. Released count layers and metadata were preserved; count, feature, and mitochondrial-fraction outliers were assessed within reliable sample partitions, and scDblFinder was required where count-level data and suitable partitions allowed it [22]. OA and OC were not integrated into a shared latent space, and labels remained dataset specific. Disease inference used pseudobulk only where biological replication and an interpretable contrast were available [23,24].

For each candidate and eligible exact source label, CCSS was calculated from within-dataset detection fractions. Labels named Unassigned, Ambiguous, or Other and strata below the cell threshold were excluded. If r was the detection fraction divided by the maximum eligible detection fraction and K the number of eligible labels, tau was sum(1-r)/(K-1), and CCSS was tau multiplied by r. Thresholds of 50, 100, and 200 cells were audited; 100 cells was primary. Medians were combined across datasets only for identical labels, without relabeling Fibroblast as cancer-associated fibroblast or Ovarian.cancer.cell as another epithelial state.

UCell v2.14.0 scored the unsigned ten-gene evidence set and a discovery-direction-compatible version from cell-level ranks with maxRank 1500 [37]. Cell scores were summarized first within biological sample and then within cell type; figures show within-atlas ranks rather than cross-disease absolute comparisons. The unsigned score represents joint rank abundance, not pathway activity, and the direction-compatible score is not a single-cell differential-expression effect.

### Sample-consensus CellChat and bounded NicheNet prior context

Communication inference was restricted to datasets with sufficient sample-resolved expression and usable cell labels: five OC samples from GSE154600 and three control plus eight OA donors from GSE255460. CellChat v2 was fitted independently to each biological sample with all human CellChatDB categories, deterministic cell-type-balanced subsampling of at most 150 cells per label, at least 20 cells per retained label, and 50 bootstrap iterations [38,39]. Sample-level significant ligand-receptor results were retained. Context consensus required support in at least 3/5 OC samples, 2/3 OA-control donors, or 4/8 OA donors. OA and OC probabilities were not compared because the tissues, label systems, and comparator designs differ; cells were never treated as independent patient replicates.

Consensus CellChat ligands were projected onto the official NicheNet v2 human ligand-target matrix from Zenodo record 7074291 [40,41]. Regulatory-potential weights linking those ligands to the fixed ten-gene set were retained as an external prior-consistency overlay. Full NicheNet ligand-activity inference was not performed because OC lacked a symmetric disease/reference receiver-cell contrast and the fixed target set was too small for an unbiased activity test. The overlay therefore does not establish ligand activity, disease-specific regulation, mediation, spatial contact, signaling flux, or mechanism.

For exploratory functional localization, the top 25 mean-expression contrast genes per analysis cluster in GSE104782 and GSE154600 were mapped to the majority dataset-specific cell label, pooled within cell type, and tested for Gene Ontology Biological Process over-representation. Candidate cell-function joins are descriptive and do not demonstrate conserved function.

### Focused upstream regulatory context

KnockTF 2.0 perturbational target sets and curated miRTarBase interactions were retained for SOX9, DDIT3, BNC1, and AKAP12 [34,35]. Target-set enrichment does not estimate transcription-factor activity, and miRNA abundance was not measured. These matrices provide testable upstream hypotheses rather than a disease-specific regulatory network.
"""


CELL_RESULTS = """### Single-cell localization was dataset specific rather than homologous across diseases

The five adapters audited 1,187,436 cells and retained 1,025,361 quality-control-pass cells (Figure 1; Table S9). The primary CCSS analysis retained exact source labels and a 100-cell minimum (Figure 4A; Tables S24a-b). In OA, recurrent high-scoring examples included CFI, EFEMP1, SOX9, and BNC1 in preHTC across three atlases. In OC, EFEMP1, DIRAS3, CFI, and AKAP12 scored highly in Fibroblast, whereas ELF3 and SOX9 scored highly in Ovarian.cancer.cell across the two OC atlases. These are within-atlas localization summaries and are not numerical OA-OC comparisons.

Sample-aware UCell placed the highest median unsigned ten-gene scores in HomC for GSE104782 and GSE255460, preHTC for GSE169454, Fibroblast for GSE154600, and Ovarian.cancer.cell for GSE180661 (Figure 4B; Table S24c). Eligible pseudobulk effects remained in Figure 4C and Table S10. The gene-cell-function matrix further linked SOX9 to an OA hypertrophic-chondrocyte mineralization context and an OC ovarian-cancer-cell translation context, and AKAP12 to OA preHTC extracellular-matrix organization and OC Fibroblast collagen-fibril organization (Table S19). Same gene therefore did not imply the same cellular meaning.

### Sample-consensus communication supplied hypotheses, not cross-disease signaling estimates

CellChat models completed independently for all 16 eligible biological samples: five GSE154600 OC samples, three GSE255460 controls, and eight GSE255460 OA donors (Figure S13; Tables S26a-b). Across samples, 49,647 significant ligand-receptor interactions were retained. Context-specific consensus rules yielded 199 context-pathway records (84 OC, 54 OA control, and 61 OA) and 526 consensus interactions whose ligand or receptor overlapped a shared DEG (Tables S26c-d). None of the fixed ten candidates occurred as a direct ligand or receptor in those consensus interactions. Accordingly, the communication layer does not provide a direct signaling explanation for the ten-gene set.

The NicheNet prior overlay contained 10,240 regulatory-potential weights connecting consensus ligands to all ten fixed candidate targets (Figure S14; Table S27). This result shows that indirect paths exist in the official prior network; it is not evidence that those ligands were active, that the candidates were regulated in vivo, or that signaling caused the observed expression directions. The predeclared feasibility table records why full ligand-activity inference was excluded (Table S28).
"""


FIGURE_LEGENDS_HEADER = "## Figure legends"


SUPPLEMENTARY_INDEX = """## Supplementary table index

- **Table S1:** Datasets, cohorts, and clinical-metadata availability.
- **Table S2:** Primary shared OA-OC DEGs.
- **Table S3a-b:** DEG threshold sensitivity.
- **Table S4:** WGCNA stability.
- **Table S5:** Strict nested machine-learning resampling and feature stability.
- **Table S6a-c:** Direction-fixed external validation and unsupervised PCA source data.
- **Table S7:** TCGA-OV model sensitivity.
- **Table S8:** Immune-signature comparisons.
- **Table S9:** Single-cell QC and downstream status.
- **Table S10:** Single-cell candidate localization and eligible pseudobulk results.
- **Table S11a-c:** Hallmark, Gene Ontology, and KEGG enrichment.
- **Table S12a-b:** Bidirectional MR estimates, provenance, and diagnostics.
- **Table S13:** Regulatory and compound interaction lookups retained as supplementary lookup data only.
- **Table S14:** HPA normal-tissue and normal-cell context.
- **Table S15a-b:** TCGA-OV relative stromal/immune context scores and correlations.
- **Table S16:** Transparent ten-gene evidence-summary matrix.
- **Table S17:** Exploratory cell-type marker GO annotation.
- **Table S18:** Complete paired Hallmark direction matrix.
- **Table S19:** Candidate gene-cell-function context matrix.
- **Table S20:** Candidate-centered Hallmark association context.
- **Table S21:** Cross-cohort molecular-separability context.
- **Table S22a-b:** External signed-score effect sizes and cross-fitted calibration sensitivity.
- **Table S23a-b:** Focused KnockTF and miRTarBase context.
- **Table S24a-c:** Dataset/context CCSS, exact-label consensus, and sample-aware UCell summaries.
- **Table S25a-e:** STRING mapping, edges, node/subgraph topology, and 10,000 label permutations.
- **Table S26a-d:** CellChat sample audit, sample interactions, consensus pathways, and shared-DEG-anchored interactions.
- **Table S27:** NicheNet v2 prior-consistency overlay.
- **Table S28:** Communication-analysis feasibility decisions and inference boundaries.
"""


REFERENCES = """36. Szklarczyk D, Kirsch R, Koutrouli M, et al. The STRING database in 2023: protein-protein association networks and functional enrichment analyses for any sequenced genome of interest. *Nucleic Acids Res*. 2023;51:D638-D646. doi:10.1093/nar/gkac1000.
37. Andreatta M, Carmona SJ. UCell: robust and scalable single-cell gene signature scoring. *Comput Struct Biotechnol J*. 2021;19:3796-3798. doi:10.1016/j.csbj.2021.06.043.
38. Jin S, Guerrero-Juarez CF, Zhang L, et al. Inference and analysis of cell-cell communication using CellChat. *Nat Commun*. 2021;12:1088. doi:10.1038/s41467-021-21246-9.
39. Jin S, Plikus MV, Nie Q. CellChat for systematic analysis of cell-cell communication from single-cell transcriptomics. *Nat Protoc*. 2025;20:180-219. doi:10.1038/s41596-024-01045-4.
40. Browaeys R, Saelens W, Saeys Y. NicheNet: modeling intercellular communication by linking ligands to target genes. *Nat Methods*. 2020;17:159-162. doi:10.1038/s41592-019-0667-5.
41. Browaeys R. NicheNet-v2: final networks and ligand-target matrices. Zenodo. 2022. doi:10.5281/zenodo.7074291.
"""


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def validate_v30(ppi_topology: Path, cellchat_audit: Path, cellchat_interactions: Path,
                 anchored: Path, nichenet: Path) -> None:
    topology = read_csv(ppi_topology)
    physical_all = [
        row for row in topology
        if row["network_type"] == "high-confidence physical"
        and row["subset"] == "all_mapped_shared_DEGs"
    ]
    if len(physical_all) != 1 or int(physical_all[0]["mapped_nodes"]) != 275:
        raise ValueError("Unexpected STRING mapping/topology result")
    if int(physical_all[0]["edges"]) != 62 or int(physical_all[0]["connected_nodes"]) != 46:
        raise ValueError("Unexpected high-confidence physical STRING graph")
    audit = [row for row in read_csv(cellchat_audit) if row["status"] == "completed"]
    if len(audit) != 16:
        raise ValueError(f"Expected 16 completed CellChat samples, found {len(audit)}")
    if len(read_csv(cellchat_interactions)) != 49647:
        raise ValueError("Unexpected number of sample-level CellChat interactions")
    if len(read_csv(anchored)) != 526:
        raise ValueError("Unexpected number of shared-DEG-anchored CellChat interactions")
    if len(read_csv(nichenet)) != 10240:
        raise ValueError("Unexpected NicheNet prior-overlay size")


def replace_section(text: str, start: str, end: str, replacement: str) -> str:
    pattern = re.compile(re.escape(start) + r".*?(?=" + re.escape(end) + r")", re.S)
    updated, count = pattern.subn(replacement.rstrip() + "\n\n", text, count=1)
    if count != 1:
        raise ValueError(f"Could not replace section beginning {start!r}")
    return updated


def build(args: argparse.Namespace) -> str:
    validate_v30(
        args.ppi_topology,
        args.cellchat_audit,
        args.cellchat_interactions,
        args.anchored,
        args.nichenet,
    )
    text = v24.build_manuscript(
        args.source,
        args.candidate_pathways,
        args.validation_context,
        args.calibration,
        args.tf_context,
        args.mirna_context,
    )
    text = re.sub(r"^# .*?$", TITLE, text, count=1, flags=re.M)
    text = text.replace(
        "**Running title:** Heterogeneous molecular overlap between OA and ovarian cancer",
        "**Running title:** Context-dependent convergence and divergence",
    )
    text = text.replace(
        "rather than shared genetic causality.",
        "rather than a simple disease-to-disease inherited-liability relationship.",
    )
    text = replace_section(text, "## Abstract", "## Introduction", ABSTRACT)
    text = text.replace(v24.STUDY_DESIGN_METHOD.rstrip(), STUDY_DESIGN.rstrip())
    text = text.replace(
        "### Functional, pathway-direction, and immune annotation",
        PPI_METHOD.rstrip() + "\n\n### Functional, pathway-direction, and immune annotation",
        1,
    )
    text = text.replace(v24.SINGLE_CELL_AND_REGULATORY_METHOD.rstrip(), SINGLE_CELL_METHOD.rstrip())

    ppi_results = """### Direction-aware STRING context showed unequal graph organization

STRING mapped 275/286 shared DEGs; 11 were retained as unmapped rather than silently discarded (Figure 2E; Tables S25a-c). The primary high-confidence physical graph contained 62 edges among 46 connected products, with 275 mapped nodes retained in density denominators. Concordant and discordant induced subgraphs contained 51 and 8 physical edges, respectively. Concordant edge count exceeded the fixed-size label-permutation null (mean 15.76; z=5.74; two-sided empirical P=9.999x10^-5), and the observed concordant-minus-discordant density difference was 0.00445 (z=3.64; P=9.999x10^-5). Discordant edge count did not differ from its null (P=0.241). The high-confidence functional graph produced a concordant sensitivity pattern (Figure S12; Tables S25d-e). These results describe organization in STRING and do not establish a common OA-OC protein mechanism.

"""
    text = text.replace(
        "### Pathway-level convergence was also directionally heterogeneous",
        ppi_results + "### Pathway-level convergence was also directionally heterogeneous",
        1,
    )
    text = replace_section(
        text,
        "### Single-cell data resolved distinct cellular and functional contexts",
        "### Focused upstream evidence defined hypotheses rather than a regulatory mechanism",
        CELL_RESULTS,
    )
    text = text.replace(
        "MR detected no causal effect in either direction.",
        "MR found no evidence of disease-to-disease genetic-liability effects under the selected datasets, instruments, and assumptions.",
    )
    text = text.replace(
        "MR was performed not to make OA or OC the presumed cause of the other disease, but to test whether transcriptomic convergence was supported by inherited genetic liability. It was not detected under the available instruments and assumptions.",
        "MR was performed not to make OA or OC the presumed cause of the other disease, but to test whether genetic liability to either disease affected risk of the other. No such evidence was detected under the available instruments and assumptions; this analysis does not test shared heritability, common susceptibility loci, or genetic correlation.",
    )

    discussion_old = """### Shared genes have different cellular and regulatory meanings

Single-cell localization showed candidate expression in cartilage cell states in OA but malignant, stromal, endothelial, mast-cell, and other immune-associated contexts in OC. The gene-cell-function matrix provided a transparent descriptive bridge, while the candidate-centered Hallmark analysis showed different residual association structures after disease-status adjustment. Same gene did not equal same cellular or pathway meaning.

The focused upstream analysis adds hypotheses rather than a completed mechanism. KnockTF target-set direction cannot be read as transcription-factor activity, and miRTarBase interaction presence cannot establish regulation in an unmeasured tissue context. The value of these matrices is to nominate experimentally testable axes while displaying exactly which evidence is absent. They complement rather than replace the stronger transcriptomic, cellular, and genetic boundary analyses.
"""
    discussion_new = """### Shared genes occupy different cellular and communication contexts

CCSS and sample-aware UCell moved the single-cell analysis beyond a presence map while retaining exact source labels and sample structure. The ten-gene set localized to HomC or preHTC contexts in OA atlases and to Fibroblast or Ovarian.cancer.cell contexts in OC atlases. These labels are not asserted to be homologous, and absolute scores were not compared across diseases. The gene-cell-function matrix and candidate-centered Hallmark analysis likewise showed that the same gene could occur in different cellular and pathway contexts.

Sample-consensus CellChat added an intercellular layer without pooling patients. Its most important result was a boundary: although 526 consensus interactions touched the broader shared-DEG set, none used a fixed candidate as a direct ligand or receptor. The NicheNet overlay supplied indirect prior paths to all ten candidates, but the available design did not support ligand-activity inference. Communication and prior-network results therefore nominate experimentally testable contexts; they do not explain the observed directions as a confirmed signaling mechanism.

STRING supplied a complementary intracellular association context. The concordant subset was unusually dense under fixed-size label permutation, driven in part by a connected cell-cycle-related component. This structure may reflect database coverage and established complexes rather than OA-OC tissue biology. It strengthens the claim that concordant and discordant genes have different association architectures, not the claim that either architecture is active in both diseases.

KnockTF and miRTarBase add further hypotheses rather than a completed regulatory mechanism. Perturbational target-set direction cannot be read as transcription-factor activity, and curated interaction presence cannot establish regulation in an unmeasured tissue context.
"""
    if discussion_old not in text:
        raise ValueError("Could not locate cellular/regulatory Discussion section")
    text = text.replace(discussion_old, discussion_new)
    text = text.replace(
        "five single-cell adapters, a gene-cell-function matrix, complete MR provenance, focused perturbational/curated regulatory matrices, exact figure source data, and a claim-evidence registry.",
        "five single-cell adapters, CCSS and sample-aware UCell, a gene-cell-function matrix, sample-resolved CellChat, a bounded NicheNet prior overlay, complete MR provenance, focused perturbational/curated regulatory matrices, exact figure source data, and a claim-evidence registry.",
    )
    text = text.replace(
        "Single-cell datasets had unequal controls, annotations, and replication. KnockTF target-set enrichment does not measure transcription-factor activity, and miRNA abundance was unavailable.",
        "Single-cell datasets had unequal controls, annotations, and replication. CellChat is expression- and prior-based, lacks spatial contact information, and was eligible in only two atlases; sample-consensus thresholds and cell-type-balanced subsampling do not remove database, annotation, or sampling bias. NicheNet ligand activity was not estimated, and the prior overlay cannot establish regulation or mediation. STRING is database-derived and not a tissue-specific interaction assay. KnockTF target-set enrichment does not measure transcription-factor activity, and miRNA abundance was unavailable.",
    )
    text = text.replace(
        "The ten genes are an interpretable evidence summary, and the focused pathway/regulatory analyses define testable contexts rather than diagnostic markers, therapeutic targets, or an established shared mechanism.",
        "The ten genes are an interpretable evidence summary. Direction-aware STRING, sample-aware single-cell scoring, sample-consensus CellChat, and the NicheNet prior overlay define distinct and testable contexts, but they do not establish diagnostic markers, therapeutic targets, intercellular causation, or a shared disease mechanism.",
    )
    text = text.replace(
        "including `run_submission_v24.ps1`",
        "including `run_submission_v30.ps1`",
    )
    text = text.replace(
        "(Figure 3; Figure S1; Table S4)",
        "(Figure S8A-B; Figure S1; Table S4)",
    )
    text = text.replace(
        "(Figure 3C; Figure S8; Table S5)",
        "(Figure S8C; Table S5)",
    )
    text = text.replace(
        "As a secondary reproducibility assessment, Figure 4 presents",
        "As a secondary reproducibility assessment, Figure 5 presents",
    )
    text = text.replace(
        "(Figure 4D; Tables S6 and S21)",
        "(Figure 5D; Tables S6 and S21)",
    )
    text = text.replace(
        "(Figure 6; Figure S1; Table S7)",
        "(Figure S1; Table S7)",
    )
    text = text.replace(
        "(Figure S2; Table S12). MR was performed",
        "(Figure 6A; Figure S2; Table S12). MR was performed",
    )

    reference_anchor = "35. Huang HY, Lin YCD, Cui S, et al. miRTarBase update 2022: an informative resource for experimentally validated miRNA-target interactions. *Nucleic Acids Res*. 2022;50:D222-D230. doi:10.1093/nar/gkab1079."
    if reference_anchor not in text:
        raise ValueError("Reference anchor not found")
    text = text.replace(reference_anchor, reference_anchor + "\n" + REFERENCES.rstrip())

    legends = args.figure_legends.read_text(encoding="utf-8").strip()
    text = replace_section(text, "## Figure legends", "## Supplementary table index", FIGURE_LEGENDS_HEADER + "\n\n" + legends)
    text = re.sub(r"## Supplementary table index\n.*\Z", SUPPLEMENTARY_INDEX.rstrip() + "\n", text, flags=re.S)
    return re.sub(r"\n{3,}", "\n\n", text).rstrip() + "\n"


def response_matrix() -> str:
    return """# OC-OA manuscript V3.0 strengthening matrix

| Recommendation | V3.0 action | Location | Status |
|---|---|---|---|
| Add direction-aware PPI rather than a ten-gene decorative network | Queried all 286 shared DEGs in STRING v12.0; retained mapped isolates, physical graph as primary, functional sensitivity, and 10,000 fixed-size direction-label permutations. Topology did not select candidates. | Figure 2D; Figure S12; Tables S25a-e | Implemented |
| Upgrade single-cell localization | Added exact-label CCSS at 50/100/200-cell thresholds, sample-aware UCell scoring, and eligible pseudobulk evidence across five atlases. | Figure 4; Tables S24a-c | Implemented |
| Add CellChat without patient pseudoreplication | Fitted CellChat independently to 5 OC samples and 11 OA/control donors, then applied prespecified sample-consensus thresholds. | Figure S13; Tables S26a-d | Implemented |
| Add NicheNet if scientifically identifiable | Added official v2 human prior overlay from consensus ligands to the fixed ten-gene targets. | Figure S14; Table S27 | Implemented as prior consistency |
| Avoid unsupported NicheNet activity claims | Did not run ligand-activity inference because OC lacks a symmetric receiver-cell disease/reference contrast and the fixed target set is too small. | Methods; Results; Table S28 | Preserved boundary |
| Keep communication context supplementary | CellChat and NicheNet remain Figures S13-S14 and are absent from prediction/causality claims. | Main narrative; Discussion | Implemented |
| Do not compare OA and OC CellChat probabilities | Reported each atlas/condition descriptively; no numerical cross-disease communication comparison was made. | Methods; Figure S13 legend | Implemented |
| Preserve candidate identities and source labels | Retained the fixed ten genes; did not introduce TNF, RELA, or IL6; did not relabel Fibroblast or Ovarian.cancer.cell. | Tables S16, S24, S27 | Implemented |
| Keep AUC secondary | Preserved PCA, effect sizes, permutation, sample omission, and ROC-last presentation. | Figure 5; Figure S10 | Preserved |
| Keep MR as a causality boundary | Limited wording to disease-to-disease genetic liability under the selected GWAS, instruments, and assumptions. | Figure 6; Figure S2; Discussion | Preserved |
| Avoid clinical-model template expansion | No DCA, nomogram, additional ML algorithm, drug prediction, or therapeutic-target claim was added. | Full project | Preserved |

## Author-controlled items still required

- Authors, affiliations, and corresponding-author details.
- Target journal and journal-specific formatting.
- Funding, competing interests, and CRediT contributions.
- Public repository URL and archival DOI after deposition.
- Final citation, statistical, language, and scientific approval.
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--candidate-pathways", required=True, type=Path)
    parser.add_argument("--validation-context", required=True, type=Path)
    parser.add_argument("--calibration", required=True, type=Path)
    parser.add_argument("--tf-context", required=True, type=Path)
    parser.add_argument("--mirna-context", required=True, type=Path)
    parser.add_argument("--ppi-topology", required=True, type=Path)
    parser.add_argument("--cellchat-audit", required=True, type=Path)
    parser.add_argument("--cellchat-interactions", required=True, type=Path)
    parser.add_argument("--anchored", required=True, type=Path)
    parser.add_argument("--nichenet", required=True, type=Path)
    parser.add_argument("--figure-legends", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--response-output", required=True, type=Path)
    args = parser.parse_args()
    text = build(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")
    args.response_output.write_text(response_matrix(), encoding="utf-8")
    print(f"Wrote {args.output}")
    print(f"Wrote {args.response_output}")


if __name__ == "__main__":
    main()
