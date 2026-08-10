param(
    [Parameter(Mandatory = $true)]
    [string]$Document,
    [Parameter(Mandatory = $true)]
    [string]$Pdf
)

$ErrorActionPreference = "Stop"
$word = New-Object -ComObject Word.Application
$word.Visible = $false
$word.DisplayAlerts = 0
try {
    $opened = $word.Documents.Open($Document, $false, $true)
    try {
        $opened.ExportAsFixedFormat($Pdf, 17)
    }
    finally {
        $opened.Close($false)
    }
}
finally {
    $word.Quit()
    [System.Runtime.InteropServices.Marshal]::ReleaseComObject($word) | Out-Null
}
