param(
    [Parameter(Mandatory=$true)][string]$StataExe,
    [Parameter(Mandatory=$true)][string]$RscriptExe,
    [string]$PythonExe = 'python',
    [string]$OutputDirectory = (Join-Path ([System.IO.Path]::GetTempPath()) ('summclust-mixed-' + [guid]::NewGuid().ToString('N')))
)
$ErrorActionPreference = 'Stop'
$packageRoot = Split-Path $PSScriptRoot -Parent
foreach ($exe in @($StataExe,$RscriptExe,$PythonExe)) {
    if (-not (Get-Command $exe -ErrorAction SilentlyContinue)) { throw "Executable not found: $exe" }
}
if (Test-Path -LiteralPath $OutputDirectory) { throw 'Choose a new output directory so stale logs cannot count as a pass.' }
New-Item -ItemType Directory -Path $OutputDirectory | Out-Null
$scratch = (Resolve-Path -LiteralPath $OutputDirectory).Path
foreach ($folder in @('Stata','R','python','tests','validation')) {
    Copy-Item -LiteralPath (Join-Path $packageRoot $folder) -Destination $scratch -Recurse
}
New-Item -ItemType Directory -Path (Join-Path $scratch 'output') | Out-Null
Push-Location $scratch
try {
    & $PythonExe validation/python_certification.py
    if ($LASTEXITCODE -ne 0) { throw 'Python certification failed.' }
    & $RscriptExe --vanilla R/certify.R
    if ($LASTEXITCODE -ne 0) { throw 'R certification failed.' }
    $process = Start-Process -FilePath $StataExe -ArgumentList '/e do tests/test_stata_mixed.do' -WorkingDirectory $scratch -WindowStyle Hidden -PassThru
    if (-not $process.WaitForExit(300000)) { throw "Stata mixed certification timed out; process $($process.Id) remains available." }
    $stataLog = Join-Path $scratch 'output/stata_mixed_tests.log'
    if (-not (Test-Path -LiteralPath $stataLog)) { throw 'Missing Stata mixed certification log.' }
    if ((Get-Content -LiteralPath $stataLog -Raw) -notmatch '(?m)^STATA MIXED CERTIFICATION PASSED\s*$') { throw 'Stata mixed certification failed.' }
    & $PythonExe validation/run_mixed.py
    if ($LASTEXITCODE -ne 0) { throw 'Python fixture export failed.' }
    & $RscriptExe --vanilla R/run_mixed_fixtures.R
    if ($LASTEXITCODE -ne 0) { throw 'R fixture export failed.' }
    & $PythonExe validation/compare_mixed.py
    if ($LASTEXITCODE -ne 0) { throw 'Cross-language mixed comparison failed.' }
    Write-Output "All mixed certification and cross-language checks passed. Logs: $scratch"
}
finally { Pop-Location }
