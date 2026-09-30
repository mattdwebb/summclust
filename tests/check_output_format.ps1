param([Parameter(Mandatory=$true)][string]$LogDirectory)
$ErrorActionPreference='Stop'
function Require([bool]$Condition,[string]$Message) { if (-not $Condition) { throw $Message } }
$scenarios=@('ordinary','options','long80','long120','negative','few','identical','suppressed')
foreach($scenario in $scenarios) {
    $texts=@()
    foreach($suffix in @('','_translated')) {
        $file=Join-Path $LogDirectory "format_${scenario}${suffix}.log"
        $raw=Get-Content -LiteralPath $file -Raw
        $match=[regex]::Match($raw,'(?ms)^OUTPUT BEGIN\r?\n(.*?)^OUTPUT END\r?$')
        Require $match.Success "Missing output markers in $file"
        # Loop-generated commands are not echoed inside the marker block.
        $body=$match.Groups[1].Value
        $texts+=,$body
        Require ($body.Contains('Regression Output')) "Missing inference title: $file"
        Require ($body.Contains('    s.e. |      Coeff   Sd. Err.   t-stat  P value    CI-lower    CI-upper')) "Historical inference headings/alignment changed: $file"
        Require ($body.Contains('Coefficients of Variation, G, and G*')) "Missing summary title: $file"
        Require ($body -match 'dimension\s+\|\s+Ng\s+Leverage\s+Partial L\.\s+beta no g\s+G\s+Gstar') "Missing summary headings: $file"
        Require ($body -notmatch 'Dimension-specific diagnostics|(?m)^\s+CV[13] max\s+\||> stic') "Old presentation/wrapping returned: $file"
        Require ($body -notmatch '(?m)^> ') "Stata hard-wrapped the output in $file"
        $width=80
        if($scenario -eq 'long120') { $width=120 }
        foreach($line in ($body -split '\r?\n')) { Require ($line.TrimEnd().Length -le $width) "Line exceeds $width columns in $file : $line" }
        foreach($method in @('CV1max','CV3max')) {
            Require ($body -match "(?m)^\s+$method \|\s+-?\d+\.\d{6}\s+-?\d+\.\d{6}\s+-?\d+\.\d{4}\s+\d+\.\d{4}\s+-?\d+\.\d{6}\s+-?\d+\.\d{6}\s*$") "Incorrect inference precision in $file"
        }
        $summary=[regex]::Match($body,'(?s)Coefficients of Variation, G, and G\*(.*?)Inference:').Groups[1].Value
        $numericRows=[regex]::Matches($summary,'(?m)^.*\|\s+(?:\d+\.\d{4}|\.)\s+(?:\d+\.\d{4}|\.)\s+(?:\d+\.\d{4}|\.)\s+(?:\d+\.\d{4}|\.)\s+\d+\s+(?:\d+\.\d{2}|\.)\s*$')
        Require ($numericRows.Count -eq 3) "Diagnostic precision/integer counts changed in $file"
        Require ($summary -match '(?m)^\s*intersect\s+\|') "Missing intersection row: $file"
        if($scenario -like 'long*') {
            Require ($summary.Contains('clustering_dimension_alpha_12345') -and $summary.Contains('clustering_dimension_beta_123456')) "Long cluster name truncated: $file"
            Require ($body.Contains('control_variable_long_name_12345')) "Long coefficient name truncated: $file"
        } else {
            Require ($summary -match '(?m)^\s*firm\s+\|') "Missing actual cluster name: $file"
            if($scenario -ne 'identical') { Require ($summary -match '(?m)^\s*year\s+\|') "Missing second cluster name: $file" }
        }
        if($scenario -eq 'negative') { Require ($body.Contains('variance is negative')) "Negative-variance warning lost: $file" }
        if($scenario -eq 'few') { Require ($body.Contains('only two clusters')) "Few-cluster warning lost: $file" }
        if($scenario -eq 'identical') { Require ($body.Contains('dimensions are identical') -and $body.Contains('components are tied')) "Identical/tied note lost: $file" }
        if($scenario -eq 'suppressed') { Require ($body.Contains('more than 52 clusters') -and $body.Contains('r(table_g)')) "Suppression note lost: $file" }
        if($scenario -eq 'options') {
            Require ($body.Contains('addmeans has no additional effect') -and $body.Contains('rank-deficient omit-one')) "Option note/warning lost: $file"
            Require ($body -match 'cluster\s+\|\s+Ng\s+Leverage\s+Partial L\.\s+beta no g') "Cluster table headings lost: $file"
            Require ($body -match 'variable\s+\|\s+Coeff\s+CV1max\s+CV3max') "Regression table headings lost: $file"
        }
    }
    # SMCL is independently rendered by Stata's translator. Ignore only end-of-line spaces.
    $plain=($texts[0] -split '\r?\n' | ForEach-Object { $_.TrimEnd() }) -join "`n"
    $smcl=($texts[1] -split '\r?\n' | ForEach-Object { $_.TrimEnd() }) -join "`n"
    Require ($plain -eq $smcl) "Text and rendered SMCL differ: $scenario"
    Write-Output "PASS formatting $scenario (text and SMCL)"
}
