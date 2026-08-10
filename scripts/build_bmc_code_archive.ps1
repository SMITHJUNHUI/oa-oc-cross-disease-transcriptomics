param(
    [Parameter(Mandatory = $true)]
    [string]$ProjectRoot,

    [Parameter(Mandatory = $true)]
    [string]$PackageRoot
)

$ErrorActionPreference = "Stop"
$project = (Resolve-Path -LiteralPath $ProjectRoot).Path
$package = (Resolve-Path -LiteralPath $PackageRoot).Path
$workspace = (Resolve-Path -LiteralPath (Join-Path $project "..\..")).Path

if (-not $project.StartsWith($workspace, [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "Project root is outside the workspace."
}
if (-not $package.StartsWith($workspace, [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "Package root is outside the workspace."
}

$buildId = [Guid]::NewGuid().ToString("N")
$stageParent = Join-Path $package ("work\code_archive_build_" + $buildId)
$stage = Join-Path $stageParent "Additional_file_3_reproducible_code"
New-Item -ItemType Directory -Path $stage -Force | Out-Null

function Copy-SelectedTree {
    param([string]$Source, [string]$Destination)
    Get-ChildItem -LiteralPath $Source -Recurse -File | Where-Object {
        $_.Extension -ne ".pyc" -and $_.FullName -notmatch "[\\/]__pycache__[\\/]"
    } | ForEach-Object {
        $relative = $_.FullName.Substring($Source.Length).TrimStart("\")
        $target = Join-Path $Destination $relative
        New-Item -ItemType Directory -Path (Split-Path -Parent $target) -Force | Out-Null
        Copy-Item -LiteralPath $_.FullName -Destination $target -Force
    }
}

Copy-Item -LiteralPath (Join-Path $package "code_archive_source\README.md") -Destination (Join-Path $stage "README.md")
Copy-Item -LiteralPath (Join-Path $package "code_archive_source\SOURCE_DATASETS.tsv") -Destination (Join-Path $stage "SOURCE_DATASETS.tsv")

foreach ($directory in @("R", "scripts", "tools", "tests")) {
    Copy-SelectedTree -Source (Join-Path $project $directory) -Destination (Join-Path $stage $directory)
}

foreach ($name in @("config.yml", "config.example.yml", "immune_signatures.yml", "submission_sensitivity.yml")) {
    $destination = Join-Path $stage ("config\" + $name)
    New-Item -ItemType Directory -Path (Split-Path -Parent $destination) -Force | Out-Null
    Copy-Item -LiteralPath (Join-Path $project ("config\" + $name)) -Destination $destination
}

foreach ($name in @("analysis_contract.md", "stage_map.md")) {
    $destination = Join-Path $stage ("docs\" + $name)
    New-Item -ItemType Directory -Path (Split-Path -Parent $destination) -Force | Out-Null
    Copy-Item -LiteralPath (Join-Path $project ("docs\" + $name)) -Destination $destination
}

foreach ($name in @("package_manifest.csv", "pipeline_status.csv", "session_info.txt")) {
    $destination = Join-Path $stage ("manifests\" + $name)
    New-Item -ItemType Directory -Path (Split-Path -Parent $destination) -Force | Out-Null
    Copy-Item -LiteralPath (Join-Path $project ("results\manifests\" + $name)) -Destination $destination
}

foreach ($name in @(
    ".Rprofile", "DESCRIPTION", "OC_OA_reproducible.Rproj", "renv.lock", "setup.R",
    "setup_project.ps1", "setup_project.bat", "run_project.R", "run_project.ps1",
    "run_project.bat", "run_submission_v42.R", "run_submission_v42.ps1"
)) {
    Copy-Item -LiteralPath (Join-Path $project $name) -Destination (Join-Path $stage $name)
}

$temporaryZip = Join-Path $stageParent "Additional_file_3_reproducible_code.zip"
Compress-Archive -LiteralPath $stage -DestinationPath $temporaryZip -CompressionLevel Optimal
$finalZip = Join-Path $package "Additional_file_3_reproducible_code.zip"
Copy-Item -LiteralPath $temporaryZip -Destination $finalZip -Force

$files = Get-ChildItem -LiteralPath $stage -Recurse -File
[PSCustomObject]@{
    Stage = $stage
    Files = $files.Count
    Zip = $finalZip
    Bytes = (Get-Item -LiteralPath $finalZip).Length
    HasFigureRebuilder = Test-Path -LiteralPath (Join-Path $stage "scripts\rebuild_bmc_figures.R")
}
