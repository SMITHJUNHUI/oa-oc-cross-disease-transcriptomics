param(
    [string]$RscriptPath = "",
    [string]$PythonPath = ""
)

$ErrorActionPreference = "Stop"
$RscriptPath = if ($RscriptPath) { $RscriptPath } else { (Get-Command Rscript -ErrorAction Stop).Source }
$PythonPath = if ($PythonPath) { $PythonPath } else { (Get-Command python -ErrorAction Stop).Source }
$ProjectRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$OutputRoot = Split-Path -Parent $ProjectRoot
$Submission = Join-Path $ProjectRoot "results\submission_v42"
$ReferenceV40 = Join-Path $ProjectRoot "results\submission_v40\reference_audit"
$ReferenceV41 = Join-Path $ProjectRoot "results\submission_v41\reference_audit"
$Manuscript = Join-Path $ProjectRoot "manuscript\OC_OA_manuscript_revision_v42.md"
$BaseDocx = Join-Path $ProjectRoot "manuscript\OC_OA_manuscript_revision_v42.docx"
$FinalDocx = Join-Path $OutputRoot "OC_OA_manuscript_revision_v42_with_figures.docx"
$FinalPdf = Join-Path $OutputRoot "OC_OA_manuscript_revision_v42_with_figures.pdf"

foreach ($Required in @($RscriptPath, $PythonPath)) {
    if (-not (Test-Path -LiteralPath $Required)) { throw "Required runtime was not found: $Required" }
}

Push-Location $ProjectRoot
try {
    & $RscriptPath "run_submission_v42.R"
    if ($LASTEXITCODE -ne 0) { throw "V4.2 figure build failed with exit code $LASTEXITCODE." }

    & $PythonPath "scripts\build_reference_packet_v40.py" --output-dir $ReferenceV40
    if ($LASTEXITCODE -ne 0) { throw "Reference verification failed with exit code $LASTEXITCODE." }

    & $PythonPath "scripts\build_reference_packet_v41.py" --source-dir $ReferenceV40 --output-dir $ReferenceV41
    if ($LASTEXITCODE -ne 0) { throw "V4.1 reference compression failed with exit code $LASTEXITCODE." }

    & $PythonPath "scripts\build_reference_packet_v42.py" --source-dir $ReferenceV41 --output-dir "$Submission\reference_audit"
    if ($LASTEXITCODE -ne 0) { throw "V4.2 reference compression failed with exit code $LASTEXITCODE." }

    & $PythonPath "scripts\build_manuscript_v42.py" `
        --references "$Submission\reference_audit\references_v42.md" `
        --figure-legends "manuscript\figure_legends_v42.md" `
        --template-output "manuscript\OC_OA_manuscript_revision_v42.template.md" `
        --output $Manuscript `
        --response-output "manuscript\OC_OA_revision_response_matrix_v42.md"
    if ($LASTEXITCODE -ne 0) { throw "V4.2 manuscript build failed with exit code $LASTEXITCODE." }

    & $PythonPath "scripts\create_manuscript_docx.py" `
        --input $Manuscript `
        --output $BaseDocx `
        --revision-label "V4.2" `
        --subtitle "Directionally heterogeneous overlap, cohort-specific replication, cellular localization and blood persistence"
    if ($LASTEXITCODE -ne 0) { throw "V4.2 DOCX build failed with exit code $LASTEXITCODE." }

    & $PythonPath "scripts\build_manuscript_with_figures_v42.py" `
        --source $BaseDocx `
        --figure-dir "$Submission\figures" `
        --legends "manuscript\figure_legends_v42.md" `
        --output $FinalDocx
    if ($LASTEXITCODE -ne 0) { throw "V4.2 figure embedding failed with exit code $LASTEXITCODE." }

    powershell.exe -NoProfile -ExecutionPolicy Bypass -File "scripts\export_docx_pdf.ps1" -Document $FinalDocx -Pdf $FinalPdf
    if ($LASTEXITCODE -ne 0) { throw "V4.2 PDF export failed with exit code $LASTEXITCODE." }

    & $PythonPath "scripts\audit_submission_v42.py" --project-root $ProjectRoot --docx $FinalDocx --pdf $FinalPdf
    if ($LASTEXITCODE -ne 0) { throw "V4.2 submission audit failed with exit code $LASTEXITCODE." }

    Write-Host "V4.2 submission package completed and audited."
    Write-Host "DOCX: $FinalDocx"
    Write-Host "PDF:  $FinalPdf"
}
finally {
    Pop-Location
}
