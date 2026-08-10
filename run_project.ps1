param(
    [ValidateSet("preflight", "core", "full", "tests")]
    [string]$Mode = "full",
    [string]$Config = "",
    [switch]$Force,
    [string]$From = "",
    [string]$To = ""
)

$ErrorActionPreference = "Stop"
$projectRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$runner = Join-Path $projectRoot "run_project.R"
if (-not $Config) {
    $localConfig = Join-Path $projectRoot "config\local.yml"
    if (Test-Path -LiteralPath $localConfig) {
        $Config = "config/local.yml"
    }
    else {
        $Config = "config/config.yml"
    }
}

$renvMarker = Join-Path $projectRoot "renv\.restored"
if ((Test-Path -LiteralPath $renvMarker) -and -not $env:OC_OA_USE_RENV) {
    $env:OC_OA_USE_RENV = "true"
}

$candidates = @()
if ($env:R_SCRIPT) {
    $candidates += $env:R_SCRIPT
}

foreach ($registryPath in @(
    "HKLM:\SOFTWARE\R-core\R",
    "HKLM:\SOFTWARE\WOW6432Node\R-core\R",
    "HKCU:\SOFTWARE\R-core\R"
)) {
    $install = Get-ItemProperty -LiteralPath $registryPath -ErrorAction SilentlyContinue
    if ($install -and $install.InstallPath) {
        $candidates += (Join-Path $install.InstallPath "bin\Rscript.exe")
    }
}

$command = Get-Command "Rscript.exe" -ErrorAction SilentlyContinue
if ($command) {
    $candidates += $command.Source
}

$programFilesR = Join-Path $env:ProgramFiles "R"
if (Test-Path -LiteralPath $programFilesR) {
    $candidates += Get-ChildItem -LiteralPath $programFilesR -Directory -Filter "R-4.5*" |
        Sort-Object Name -Descending |
        ForEach-Object { Join-Path $_.FullName "bin\Rscript.exe" }
}

$rscript = $candidates |
    Where-Object { $_ -and (Test-Path -LiteralPath $_) } |
    Select-Object -First 1

if (-not $rscript) {
    throw "Rscript.exe was not found. Set R_SCRIPT or install R 4.5.x."
}

$arguments = @(
    $runner,
    "--mode=$Mode",
    "--config=$Config"
)

if ($Force) {
    $arguments += "--force"
}
if ($From) {
    $arguments += "--from=$From"
}
if ($To) {
    $arguments += "--to=$To"
}

Push-Location $projectRoot
try {
    foreach ($localeVariable in @("LANG", "LC_ALL", "LC_CTYPE", "LC_COLLATE", "LC_TIME", "LC_MONETARY")) {
        Remove-Item -LiteralPath "Env:$localeVariable" -ErrorAction SilentlyContinue
    }
    & $rscript @arguments
    if ($LASTEXITCODE -ne 0) {
        throw "R exited with code $LASTEXITCODE."
    }
}
finally {
    Pop-Location
}
