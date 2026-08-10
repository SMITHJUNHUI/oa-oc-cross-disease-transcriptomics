param()

$ErrorActionPreference = "Stop"
$projectRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$setupScript = Join-Path $projectRoot "setup.R"

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
    throw "Rscript.exe was not found. Install R 4.5.x or set R_SCRIPT."
}

Push-Location $projectRoot
try {
    foreach ($localeVariable in @("LANG", "LC_ALL", "LC_CTYPE", "LC_COLLATE", "LC_TIME", "LC_MONETARY")) {
        Remove-Item -LiteralPath "Env:$localeVariable" -ErrorAction SilentlyContinue
    }
    & $rscript $setupScript "--restore"
    if ($LASTEXITCODE -ne 0) {
        throw "Dependency restore exited with code $LASTEXITCODE."
    }
}
finally {
    Pop-Location
}
