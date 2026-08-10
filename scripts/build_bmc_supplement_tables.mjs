import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";


const [inputDir, outputPath, previewDir, auditPath] = process.argv.slice(2);
if (!inputDir || !outputPath || !previewDir || !auditPath) {
  throw new Error("Usage: node build_bmc_supplement_tables.mjs <inputDir> <outputPath> <previewDir> <auditPath>");
}

const mappings = [
  ["S1_Cohorts", "Table_S1_tissue_data_sources_and_cohorts.csv", "Table S1", "Tissue cohorts, platforms, sample groups and analysis roles"],
  ["S2_Shared_DEGs", "Table_S2_shared_tissue_DEGs.csv", "Table S2", "Shared tissue differentially expressed genes and direction classes"],
  ["S3a_Thresholds", "Table_S3a_DEG_threshold_summary.csv", "Table S3a", "Differential-expression threshold sensitivity summary"],
  ["S3b_Membership", "Table_S3b_DEG_threshold_membership.csv", "Table S3b", "Gene membership across prespecified thresholds"],
  ["S4_Ext_effects", "Table_S4_external_tissue_gene_effects.csv", "Table S4", "External tissue gene-level effects"],
  ["S5_Direction", "Table_S5_external_tissue_direction_summary.csv", "Table S5", "External tissue direction-agreement summary"],
  ["S6_Hallmark", "Table_S6_external_tissue_Hallmark_GSEA.csv", "Table S6", "External tissue Hallmark enrichment"],
  ["S7_Illustrative", "Table_S7_candidate_evidence_summary.csv", "Table S7", "Illustrative-gene evidence summary"],
  ["S8a_scRNA_QC", "Table_S8a_single_cell_QC_and_status.csv", "Table S8a", "Single-cell compatibility, quality control and permitted analysis layers"],
  ["S8b_scRNA_detect", "Table_S8b_illustrative_gene_single_cell_detection.csv", "Table S8b", "Illustrative-gene detection by exact source label"],
  ["S9_Blood_audit", "Table_S9_blood_dataset_audit.csv", "Table S9", "Peripheral-blood cohort audit"],
  ["S10_Blood_FDR", "Table_S10_FDR_supported_blood_component.csv", "Table S10", "FDR-supported blood-persistent association"],
  ["S11_Attrition", "Table_S11_blood_screen_attrition.csv", "Table S11", "Prespecified blood-screen attrition"],
  ["S12_GO", "Table_S12_GO_shared_genes.csv", "Table S12", "Gene Ontology enrichment of shared tissue genes"],
  ["S13_KEGG", "Table_S13_KEGG_shared_genes.csv", "Table S13", "Kyoto Encyclopedia of Genes and Genomes enrichment"],
  ["S14_Hallmark_dir", "Table_S14_discovery_Hallmark_direction.csv", "Table S14", "Discovery-cohort Hallmark direction matrix"],
];

await fs.mkdir(path.dirname(outputPath), { recursive: true });
await fs.mkdir(previewDir, { recursive: true });

const workbook = Workbook.create();
const readme = workbook.worksheets.add("README");
readme.showGridLines = false;
const readmeRows = [
  ["Additional file 1", "Supplementary tables for the OA-OC BMC Medical Genomics submission", ""],
  ["Manuscript title", "Directionally heterogeneous transcriptomic overlap between osteoarthritis and ovarian cancer across tissue, cellular and blood contexts", ""],
  ["Purpose", "Machine-readable source tables supporting the main manuscript and Additional file 2", ""],
  ["Version", "BMC submission package, 9 August 2026", ""],
  ["Data source", "Public NCBI GEO datasets; accessions and cohort roles are listed in Table S1", ""],
  ["Code archive", "AUTHOR INPUT REQUIRED: insert public repository URL and archived release DOI", ""],
  ["Table", "Worksheet", "Description"],
  ...mappings.map(([sheetName, _fileName, tableName, description]) => [tableName, sheetName, description]),
];
readme.getRangeByIndexes(0, 0, readmeRows.length, 3).values = readmeRows;
readme.getRange("A1:C1").merge();
readme.getRange("A1").values = [["Additional file 1: Supplementary tables"]];
readme.getRange("A1:C1").format = {
  font: { name: "Arial", size: 14, bold: true, color: "#000000" },
  borders: { bottom: { style: "medium", color: "#000000" } },
};
readme.getRange("A2:C6").format = { font: { name: "Arial", size: 10, color: "#000000" }, wrapText: true };
readme.getRange("A7:C7").format = {
  font: { name: "Arial", size: 10, bold: true, color: "#000000" },
  borders: { bottom: { style: "thin", color: "#000000" } },
};
readme.getRangeByIndexes(7, 0, mappings.length, 3).format = {
  font: { name: "Arial", size: 10, color: "#000000" },
  borders: { bottom: { style: "hair", color: "#D9D9D9" } },
  wrapText: true,
};
readme.getRange("A1:C24").format.autofitRows();
readme.getRange("A1:A24").format.columnWidth = 22;
readme.getRange("B1:B24").format.columnWidth = 30;
readme.getRange("C1:C24").format.columnWidth = 75;
readme.freezePanes.freezeRows(7);

const sheetAudits = [];
for (const [sheetName, fileName, tableName, description] of mappings) {
  const csvPath = path.join(inputDir, fileName);
  const csvText = await fs.readFile(csvPath, "utf8");
  const imported = await Workbook.fromCSV(csvText, { sheetName });
  const importedSheet = imported.worksheets.getItem(sheetName);
  const importedValues = importedSheet.getUsedRange(true).values;
  const sheet = workbook.worksheets.add(sheetName);
  const importedCols = Math.max(...importedValues.map((row) => row.length));
  sheet.getRangeByIndexes(0, 0, importedValues.length, importedCols).values = importedValues;
  sheet.showGridLines = false;
  sheet.freezePanes.freezeRows(1);
  const used = sheet.getUsedRange(true);
  const values = used.values;
  const rowCount = values.length;
  const colCount = Math.max(...values.map((row) => row.length));
  used.format = {
    font: { name: "Arial", size: 9, color: "#000000" },
    borders: { bottom: { style: "hair", color: "#E6E6E6" } },
    wrapText: false,
  };
  sheet.getRangeByIndexes(0, 0, 1, colCount).format = {
    font: { name: "Arial", size: 9, bold: true, color: "#000000" },
    borders: { bottom: { style: "thin", color: "#000000" } },
    wrapText: true,
  };
  used.format.autofitColumns();
  sheet.getRangeByIndexes(0, 0, Math.min(rowCount, 80), colCount).format.autofitRows();
  for (let col = 0; col < colCount; col += 1) {
    const sample = values.slice(0, Math.min(rowCount, 150)).map((row) => String(row[col] ?? ""));
    const longest = Math.max(8, ...sample.map((value) => value.length));
    const width = Math.min(32, Math.max(10, Math.ceil(longest * 0.92)));
    sheet.getRangeByIndexes(0, col, rowCount, 1).format.columnWidth = width;
  }
  sheetAudits.push({ tableName, sheetName, description, source: fileName, rowCount, colCount });
}

const inspectSheets = await workbook.inspect({ kind: "sheet", include: "id,name", maxChars: 12000 });
const inspectFormulas = await workbook.inspect({ kind: "formula", maxChars: 4000, options: { maxResults: 100 } });

for (const sheetName of ["README", ...mappings.map((entry) => entry[0])]) {
  try {
    const preview = await workbook.render({ sheetName, autoCrop: "all", scale: 0.35, format: "png" });
    const bytes = new Uint8Array(await preview.arrayBuffer());
    await fs.writeFile(path.join(previewDir, `${sheetName}.png`), bytes);
  } catch (error) {
    await fs.writeFile(path.join(previewDir, `${sheetName}.render_error.txt`), String(error), "utf8");
  }
}

const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(outputPath);
await fs.writeFile(
  auditPath,
  JSON.stringify(
    {
      generated: new Date().toISOString(),
      sheetCount: 1 + mappings.length,
      sheets: sheetAudits,
      inspectSheets: inspectSheets.ndjson,
      inspectFormulas: inspectFormulas.ndjson,
      excludedFiles: [
        "Table_S15a_STRING_mapping_audit.csv",
        "Table_S15b_STRING_edges.csv",
        "Table_S15c_STRING_node_topology.csv",
      ],
    },
    null,
    2,
  ),
  "utf8",
);
console.log(`Wrote ${outputPath}`);
