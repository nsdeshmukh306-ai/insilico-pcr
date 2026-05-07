# publish_to_github.ps1
# Run from E:\insilico_pcr in PowerShell

$ErrorActionPreference = "Stop"
$ProjectDir  = Split-Path -Parent $MyInvocation.MyCommand.Path
$GithubUser  = "nsdeshmukh306-ai"
$RepoName    = "insilico-pcr"
$RepoDesc    = "Computational PCR simulation: SantaLucia 1998 NN thermodynamics, FM-index genome search, and an interactive analysis dashboard."

Write-Host ""
Write-Host "======================================" -ForegroundColor Cyan
Write-Host "  In-Silico PCR  ->  GitHub Publisher" -ForegroundColor Cyan
Write-Host "======================================" -ForegroundColor Cyan
Write-Host ""

# Step 1: Remove old broken .git
$gitDir = Join-Path $ProjectDir ".git"
if (Test-Path $gitDir) {
    Write-Host "[1/5] Removing old .git directory..." -ForegroundColor Yellow
    Remove-Item -Recurse -Force $gitDir
    Write-Host "      Done." -ForegroundColor Green
} else {
    Write-Host "[1/5] No existing .git directory - skipping." -ForegroundColor Gray
}

# Step 2: Init fresh repo and commit
Write-Host ""
Write-Host "[2/5] Initialising git repo..." -ForegroundColor Yellow
Set-Location $ProjectDir
git init -b main
git config user.name $GithubUser
git config user.email "niraj_20254009@students.iisertirupati.ac.in"
git add .
git commit -m "Initial commit: In-Silico PCR v1.2.0"
Write-Host "      Committed." -ForegroundColor Green

# Step 3: GitHub PAT
Write-Host ""
Write-Host "[3/5] GitHub authentication" -ForegroundColor Yellow
Write-Host "      You need a Personal Access Token with 'repo' scope."
Write-Host "      Get one at: https://github.com/settings/tokens/new" -ForegroundColor Cyan
Write-Host "      (tick the 'repo' box, set expiry, click Generate Token)"
Write-Host ""
$SecureToken = Read-Host "      Paste your GitHub PAT here" -AsSecureString
$PlainToken  = [Runtime.InteropServices.Marshal]::PtrToStringAuto(
    [Runtime.InteropServices.Marshal]::SecureStringToBSTR($SecureToken)
)

# Step 4: Create repo via GitHub API
Write-Host ""
Write-Host "[4/5] Creating https://github.com/$GithubUser/$RepoName ..." -ForegroundColor Yellow

$Headers = @{
    Authorization          = "Bearer $PlainToken"
    Accept                 = "application/vnd.github+json"
    "X-GitHub-Api-Version" = "2022-11-28"
}
$Body = @{ name = $RepoName; description = $RepoDesc; private = $false; auto_init = $false } | ConvertTo-Json

try {
    $Response = Invoke-RestMethod -Uri "https://api.github.com/user/repos" -Method POST -Headers $Headers -Body $Body -ContentType "application/json"
    Write-Host "      Repository created: $($Response.html_url)" -ForegroundColor Green
} catch {
    $Status = $_.Exception.Response.StatusCode.value__
    if ($Status -eq 422) {
        Write-Host "      Repository already exists - continuing." -ForegroundColor Yellow
    } else {
        Write-Host "      API error ($Status): $($_.Exception.Message)" -ForegroundColor Red
        Write-Host "      If this keeps failing, create the repo manually at https://github.com/new" -ForegroundColor Yellow
    }
}

# Step 5: Push
Write-Host ""
Write-Host "[5/5] Pushing to GitHub..." -ForegroundColor Yellow
$RemoteUrl = "https://$PlainToken@github.com/$GithubUser/$RepoName.git"
git remote remove origin 2>$null
git remote add origin $RemoteUrl
git push -u origin main

Write-Host ""
Write-Host "======================================" -ForegroundColor Green
Write-Host "  Done! Repo is live at:" -ForegroundColor Green
Write-Host "  https://github.com/$GithubUser/$RepoName" -ForegroundColor Cyan
Write-Host "======================================" -ForegroundColor Green
Write-Host ""
$PlainToken = $null
