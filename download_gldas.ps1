<#
  download_gldas.ps1
  Downloads the GLDAS .nc4 files listed in a GES DISC "subset" link file and
  files each one into  <Root>\YYYY\DDD\  (the layout ParBal's make_ldas_filelist.m expects).

  PREREQUISITES (one time):
    1. A free NASA Earthdata account:           https://urs.earthdata.nasa.gov/
    2. Authorize "NASA GESDISC DATA ARCHIVE":   https://disc.gsfc.nasa.gov/earthdata-login
    3. A credentials file at $env:USERPROFILE\_netrc containing ONE line:
           machine urs.earthdata.nasa.gov login YOUR_USERNAME password YOUR_PASSWORD

  USAGE (from PowerShell):
    cd E:\ucsb\code\ParBal
    .\download_gldas.ps1
  or point it at a specific list / output folder:
    .\download_gldas.ps1 -List "$env:USERPROFILE\Downloads\subset_...txt" -Root "E:\ucsb\data\ParBal\Shasta\GLDAS"

  It is RESUMABLE: re-run it any time; files already present are skipped.
  Failed downloads are logged to <Root>\_failed.txt so you can see what to retry.
#>
param(
    [string]$List   = "",
    [string]$Root   = "E:\ucsb\data\ParBal\Shasta\GLDAS",
    [string]$Netrc  = "$env:USERPROFILE\_netrc"
)

$ErrorActionPreference = 'Stop'

# --- locate curl + credentials -------------------------------------------------
$curl = (Get-Command curl.exe -ErrorAction SilentlyContinue).Source
if (-not $curl) { throw "curl.exe not found (it ships with Windows 11 - check your PATH)." }
if (-not (Test-Path $Netrc)) {
    throw "Credentials file not found: $Netrc`n" +
          "Create it with ONE line:`n" +
          "  machine urs.earthdata.nasa.gov login YOUR_USERNAME password YOUR_PASSWORD"
}

# --- find the newest subset link file if none was given ------------------------
if (-not $List) {
    $List = (Get-ChildItem "$env:USERPROFILE\Downloads\subset_GLDAS_*.txt" -ErrorAction SilentlyContinue |
             Sort-Object LastWriteTime -Descending | Select-Object -First 1).FullName
    if (-not $List) { throw "No subset_GLDAS_*.txt found in Downloads - pass -List explicitly." }
}
Write-Host "List   : $List"
Write-Host "Output : $Root"
Write-Host ""

$cookies = Join-Path $Root "_urs_cookies.txt"
$failed  = Join-Path $Root "_failed.txt"
New-Item -ItemType Directory -Force -Path $Root | Out-Null
if (Test-Path $failed) { Remove-Item $failed }

# --- only the actual data URLs (skip the README pdf line) ----------------------
$urls = Get-Content $List | Where-Object { $_ -match 'HTTP_services\.cgi' }
$n = $urls.Count
Write-Host "Found $n files to fetch.`n"

$i = 0; $got = 0; $skip = 0; $fail = 0
foreach ($u in $urls) {
    $i++
    # year / day-of-year from the FILENAME path:  ...2.1%2FYYYY%2FDDD%2F...
    if ($u -match '2\.1%2F(\d{4})%2F(\d{3})%2F') { $yr = $Matches[1]; $doy = $Matches[2] }
    else { Add-Content $failed "NO_YR_DOY  $u"; $fail++; continue }

    # output filename from the LABEL= parameter
    if ($u -match 'LABEL=([^&]+)') { $name = [uri]::UnescapeDataString($Matches[1]) }
    else { $name = "gldas_$i.nc4" }

    $dir = Join-Path $Root "$yr\$doy"
    New-Item -ItemType Directory -Force -Path $dir | Out-Null
    $out = Join-Path $dir $name

    if ((Test-Path $out) -and (Get-Item $out).Length -gt 0) {
        $skip++
        if ($i % 200 -eq 0) { Write-Host "[$i/$n] ... (skipping existing)" }
        continue
    }

    & $curl -sS -L -b $cookies -c $cookies --netrc-file $Netrc -o $out $u
    if ($LASTEXITCODE -eq 0 -and (Test-Path $out) -and (Get-Item $out).Length -gt 0) {
        $got++
        if ($got % 50 -eq 0) { Write-Host "[$i/$n] downloaded $got so far ... ($yr\$doy\$name)" }
    } else {
        $fail++
        Add-Content $failed $u
        if (Test-Path $out) { Remove-Item $out }   # drop empty/partial file
    }
}

Write-Host ""
Write-Host "DONE.  downloaded=$got  skipped(existing)=$skip  failed=$fail  of $n"
if ($fail -gt 0) { Write-Host "Failed URLs logged to: $failed  (re-run this script to retry them)" }
