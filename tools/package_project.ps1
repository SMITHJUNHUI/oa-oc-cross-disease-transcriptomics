param()

$ErrorActionPreference = "Stop"
$projectRoot = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$outputRoot = Split-Path -Parent $projectRoot
$sourceZip = Join-Path $outputRoot "OC_OA_reproducible_project_source.zip"
$resultsZip = Join-Path $outputRoot "OC_OA_verified_results.zip"
$checksumPath = Join-Path $outputRoot "OC_OA_package_SHA256.txt"
$fixedTimestamp = [DateTimeOffset]::new(
    2026, 7, 26, 0, 0, 0,
    [TimeSpan]::Zero
)

Add-Type -AssemblyName System.IO.Compression

function New-DeterministicZip {
    param(
        [Parameter(Mandatory = $true)]
        [string]$ArchivePath,
        [Parameter(Mandatory = $true)]
        [object[]]$Items,
        [Parameter(Mandatory = $true)]
        [string]$Prefix
    )

    $resolvedOutputRoot = [System.IO.Path]::GetFullPath($outputRoot)
    $resolvedArchive = [System.IO.Path]::GetFullPath($ArchivePath)
    if (-not $resolvedArchive.StartsWith(
        $resolvedOutputRoot + [System.IO.Path]::DirectorySeparatorChar,
        [System.StringComparison]::OrdinalIgnoreCase
    )) {
        throw "Archive target must remain inside the output directory."
    }
    if (Test-Path -LiteralPath $resolvedArchive) {
        Remove-Item -LiteralPath $resolvedArchive -Force
    }

    $fileStream = [System.IO.File]::Open(
        $resolvedArchive,
        [System.IO.FileMode]::CreateNew,
        [System.IO.FileAccess]::ReadWrite,
        [System.IO.FileShare]::None
    )
    try {
        $archive = [System.IO.Compression.ZipArchive]::new(
            $fileStream,
            [System.IO.Compression.ZipArchiveMode]::Create,
            $false
        )
        try {
            foreach ($item in $Items) {
                $entryName = (
                    $Prefix.TrimEnd("/") + "/" + $item.RelativePath.Replace("\", "/")
                )
                $entry = $archive.CreateEntry(
                    $entryName,
                    [System.IO.Compression.CompressionLevel]::Optimal
                )
                $entry.LastWriteTime = $fixedTimestamp
                $input = [System.IO.File]::OpenRead($item.FullName)
                $output = $entry.Open()
                try {
                    $input.CopyTo($output)
                }
                finally {
                    $output.Dispose()
                    $input.Dispose()
                }
            }
        }
        finally {
            $archive.Dispose()
        }
    }
    finally {
        $fileStream.Dispose()
    }
}

$sourceFiles = Get-ChildItem -LiteralPath $projectRoot -Recurse -Force -File |
    Where-Object {
        $relative = $_.FullName.Substring($projectRoot.Length + 1)
        $inResults = $relative.StartsWith(
            "results\",
            [System.StringComparison]::OrdinalIgnoreCase
        )
        $keepResultsPlaceholder = $relative -in @(
            "results\README.md",
            "results\.gitkeep"
        )
        $excluded = (
            $relative -eq "config\local.yml" -or
            $relative.StartsWith("renv\library\", [System.StringComparison]::OrdinalIgnoreCase) -or
            $relative.StartsWith("renv\staging\", [System.StringComparison]::OrdinalIgnoreCase) -or
            $relative -eq "renv\.restored" -or
            $relative.StartsWith(".Rproj.user\", [System.StringComparison]::OrdinalIgnoreCase) -or
            $relative.StartsWith(".git\", [System.StringComparison]::OrdinalIgnoreCase)
        )
        (-not $excluded) -and ((-not $inResults) -or $keepResultsPlaceholder)
    } |
    ForEach-Object {
        [pscustomobject]@{
            FullName = $_.FullName
            RelativePath = $_.FullName.Substring($projectRoot.Length + 1)
        }
    } |
    Sort-Object RelativePath

$resultsFiles = Get-ChildItem -LiteralPath (Join-Path $projectRoot "results") -Recurse -File |
    Where-Object {
        $relative = $_.FullName.Substring($projectRoot.Length + 1)
        $relative.StartsWith("results\tables\", [System.StringComparison]::OrdinalIgnoreCase) -or
        $relative.StartsWith("results\figures\", [System.StringComparison]::OrdinalIgnoreCase) -or
        $relative.StartsWith("results\reports\", [System.StringComparison]::OrdinalIgnoreCase) -or
        $relative.StartsWith("results\single_cell\", [System.StringComparison]::OrdinalIgnoreCase) -or
        $relative -in @(
            "results\manifests\pipeline_status.csv",
            "results\manifests\package_manifest.csv",
            "results\manifests\single_cell_qc_parameters.yml",
            "results\manifests\single_cell_runtime_manifest.csv",
            "results\manifests\session_info.txt"
        )
    } |
    ForEach-Object {
        [pscustomobject]@{
            FullName = $_.FullName
            RelativePath = $_.FullName.Substring($projectRoot.Length + 1)
        }
    }

$resultsFiles += [pscustomobject]@{
    FullName = Join-Path $projectRoot "VALIDATION_REPORT.md"
    RelativePath = "VALIDATION_REPORT.md"
}
$resultsFiles = $resultsFiles | Sort-Object RelativePath

New-DeterministicZip -ArchivePath $sourceZip -Items $sourceFiles -Prefix "OC_OA_reproducible_project"
New-DeterministicZip -ArchivePath $resultsZip -Items $resultsFiles -Prefix "OC_OA_verified_results"

$checksumLines = foreach ($path in @($sourceZip, $resultsZip)) {
    $hash = (Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash.ToLowerInvariant()
    "$hash  $([System.IO.Path]::GetFileName($path))"
}
[System.IO.File]::WriteAllLines(
    $checksumPath,
    $checksumLines,
    [System.Text.UTF8Encoding]::new($false)
)

Get-Item -LiteralPath $sourceZip, $resultsZip, $checksumPath |
    Select-Object Name, Length, LastWriteTime
