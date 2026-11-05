param(
    [Parameter(Mandatory=$true)][string]$StataExe,
    [string]$OutputDirectory = (Join-Path ([System.IO.Path]::GetTempPath()) ('summclust-tests-' + [guid]::NewGuid().ToString('N')))
)
$ErrorActionPreference = 'Stop'
$packageRoot = Split-Path $PSScriptRoot -Parent
if (-not (Test-Path -LiteralPath $StataExe -PathType Leaf)) { throw 'Stata executable not found.' }
if (Test-Path -LiteralPath $OutputDirectory) { throw 'Choose a new output directory so stale logs cannot count as a pass.' }
New-Item -ItemType Directory -Path $OutputDirectory | Out-Null
$testRoot = (Resolve-Path -LiteralPath $OutputDirectory).Path
Copy-Item -Path (Join-Path $packageRoot 'Stata/*') -Destination $testRoot
Copy-Item -Path (Join-Path $PSScriptRoot '*.do') -Destination $testRoot
$suites = @(
    'certify_summclust',
    'certify_summclust_twoway',
    'certify_summclust_twoway_api',
    'certify_twowayjack_delegate',
    'certify_integration',
    'certify_output_numerics',
    'certify_output_format',
    'certify_string_fixed_effects',
    'verify_net_install_twoway'
)
foreach ($suite in $suites) {
    $process = Start-Process -FilePath $StataExe -ArgumentList "/e do $suite.do" -WorkingDirectory $testRoot -WindowStyle Hidden -PassThru
    if (-not $process.WaitForExit(300000)) { throw "Suite timed out; Stata process $($process.Id) remains available for inspection: $suite" }
    $logPath = Join-Path $testRoot ($suite + '_results.log')
    if (-not (Test-Path -LiteralPath $logPath)) { throw "Missing certification log: $suite" }
    $logText = Get-Content -LiteralPath $logPath -Raw
    # Stata batch exit codes do not reliably expose do-file failures. Require
    # the actual success output line, not the echoed display command.
    if ($logText -notmatch '(?m)^(SUMMCLUST .*CERTIFICATION PASSED|THIN TWOWAYJACK DELEGATE CERTIFICATION PASSED|TWO-WAY API CERTIFICATION PASSED|LOCAL NET INSTALL TWO-WAY VERIFICATION PASSED)\r?$') {
        throw "Certification failed: $suite. Inspect $logPath"
    }
    Write-Output "PASS $suite"
}
& (Join-Path $PSScriptRoot 'check_output_format.ps1') -LogDirectory $testRoot
Write-Output "All $($suites.Count) suites and text/SMCL formatting checks passed. Logs: $testRoot"
