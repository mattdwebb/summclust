param([Parameter(Mandatory=$true)][string]$LogDirectory)
$ErrorActionPreference='Stop'
function Require([bool]$Condition,[string]$Message) { if (-not $Condition) { throw $Message } }
$scenarios=@('ordinary','options','long80','long120','negative','few','identical','suppressed','mixed','mixed_swapped','twoway','mixed_twoway','long_mixed','mixed_tied')
$sources=Import-Csv -LiteralPath (Join-Path $LogDirectory 'format_sources.csv')
Require ($sources.Count -eq 2*$scenarios.Count) 'Missing independent source-variance records'
foreach($scenario in $scenarios) {
    $isMixed=$scenario -in @('mixed','mixed_swapped','mixed_twoway','long_mixed','mixed_tied')
    $secondMethod=if($isMixed) {'CV31max'} else {'CV3max'}
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
        foreach($method in @('CV1max',$secondMethod)) {
            Require ($body -match "(?m)^\s+$method \|\s+-?\d+\.\d{6}\s+-?\d+\.\d{6}\s+(?:-?\d+\.\d{4}|\.)\s+(?:\d+\.\d{4}|\.)\s+-?\d+\.\d{6}\s+-?\d+\.\d{6}\s*$") "Incorrect inference precision in $file"
        }
        # Source wording must agree with the full-precision component variances.
        # Flatten whitespace only for wrapped prose; table layout is checked above.
        $flat=[regex]::Replace($body,'\s+',' ')
        Require (-not $flat.Contains('Max selectors:')) "Obsolete selector shorthand remains: $file"
        $scenarioSources=@($sources | Where-Object { $_.scenario -eq $scenario })
        Require ($scenarioSources.Count -eq 2) "Missing source records: $scenario"
        foreach($source in $scenarioSources) {
            $method=$source.method.ToUpperInvariant()+'max'
            $estimator=$source.method.ToUpperInvariant()
            if($source.method -eq 'cv31') { $estimator='CV31 (mixed)' }
            $component=if($source.method -eq 'cv1') {'CV1'} else {'CV3'}
            $names=@('two-way','G','H')
            $labels=@("conventional two-way $estimator", "one-way $component clustered by $($source.cluster_g) (G)", "one-way $component clustered by $($source.cluster_h) (H)")
            $values=@($source.two_way_var,$source.g_var,$source.h_var) | ForEach-Object { [double]::Parse($_,[Globalization.CultureInfo]::InvariantCulture) }
            $maximum=($values | Measure-Object -Maximum).Maximum
            $scale=[Math]::Max([Math]::Abs($maximum),1e-300)
            $selected=[Array]::IndexOf($names,$source.selected)
            Require ($selected -ge 0) "Invalid source selector: $scenario $method"
            Require ([Math]::Abs($values[$selected]-$maximum) -le 1e-14*$scale) "Selector does not attain the maximum: $scenario $method"
            $expected="$method uses $($labels[$selected])."
            Require ($flat.Contains($expected)) "Wrong or missing source note: $file : $expected"
            Require ($flat.IndexOf($expected) -lt $flat.IndexOf('Coefficients of Variation')) "Source note is buried below diagnostics: $file"
            $tied=@()
            $mask=0
            for($i=0;$i -lt 3;$i++) {
                if([Math]::Abs($values[$i]-$maximum) -le 1e-10*$scale) {
                    $tied+=$labels[$i]
                    $mask+=[Math]::Pow(2,$i)
                }
            }
            Require ($mask -eq [int]$source.ties) "Tie mask disagrees with component variances: $scenario $method"
            if($tied.Count -gt 1) {
                $expectedTie="$method components are tied within relative tolerance 1e-10: $($tied -join '; ')."
                Require ($flat.Contains($expectedTie)) "Incomplete or incorrect named tie list: $file : $expectedTie"
            } else {
                Require (-not $flat.Contains("$method components are tied")) "Spurious tie note: $file"
            }
        }
        # These fixed seeds ensure the scenarios actually exercise each branch.
        $focalSource=$scenarioSources | Where-Object { $_.method -ne 'cv1' }
        if($scenario -in @('mixed','long_mixed')) { Require ($focalSource.selected -eq 'H') "H-winner fixture changed: $scenario" }
        if($scenario -eq 'mixed_swapped') { Require ($focalSource.selected -eq 'G') 'Swapped G-winner fixture changed' }
        if($scenario -in @('twoway','mixed_twoway')) { Require ($focalSource.selected -eq 'two-way') "Conventional-winner fixture changed: $scenario" }
        if($scenario -eq 'mixed_tied') { Require ([int]$focalSource.ties -eq 7) 'Mixed three-component tie fixture changed' }
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

# Render the actual help file through Stata, rather than merely balancing braces.
$help=Get-Content -LiteralPath (Join-Path $LogDirectory 'help_rendered.log') -Raw
Require ($help -notmatch '\{(?:pstd|p_end|phang|cmd:|title:|synopt|opt:|it:)') 'Unrendered SMCL directives in help'
Require ($help -match '(?m)^\s*\. keep if inrange\(age,20,40\)\s*$') 'First sample-restriction example is not a separate line'
Require ($help -match '(?m)^\s*\. summclust ln_wage msp union race,') 'Command after sample restriction did not render separately'
Require ($help -match '(?m)^\s*\. keep if inrange\(age,25,35\)\s*$') 'Mixed example sample restriction did not render separately'
Require ($help -match '(?m)^\s*\. generate vismin = inrange\(race,2,3\)\s*$') 'Mixed example generation command did not render separately'
Require ($help.Contains('Author') -and $help.Contains('Citation') -and $help.Contains('matt.webb@carleton.ca')) 'Help rendering lost sections after examples'
Write-Output 'PASS help rendering (examples, later sections, and SMCL directives)'
