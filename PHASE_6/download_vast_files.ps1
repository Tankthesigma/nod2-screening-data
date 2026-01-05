# Download all files from Vast.ai instances
# Run this AFTER simulations are complete

$localDir = "C:\Users\vasud\nod2-screening-data\PHASE_6\trajectories"

Write-Host "=== Downloading from 5090 instance ===" -ForegroundColor Cyan

# 5090 instance files
$files5090 = @(
    "febuxostat_rep2.dcd",
    "febuxostat_rep2_solvated.pdb",
    "febuxostat_rep3.dcd",
    "febuxostat_rep3_solvated.pdb",
    "ursodiol_rep1.dcd",
    "ursodiol_rep1_solvated.pdb",
    "ursodiol_rep2.dcd",
    "ursodiol_rep2_solvated.pdb",
    "ursodiol_rep3.dcd",
    "ursodiol_rep3_solvated.pdb"
)

foreach ($file in $files5090) {
    Write-Host "Downloading $file from 5090..." -ForegroundColor Yellow
    scp -P 38103 "root@ssh4.vast.ai:/workspace/nod2-screening-data/PHASE_6/trajectories/$file" "$localDir\$file"
}

Write-Host ""
Write-Host "=== Downloading from other GPU instance ===" -ForegroundColor Cyan

# Other GPU instance files
$filesOther = @(
    "febuxostat_rep1.dcd",
    "febuxostat_rep1_solvated.pdb",
    "budesonide_rep1.dcd",
    "budesonide_rep1_solvated.pdb",
    "budesonide_rep2.dcd",
    "budesonide_rep2_solvated.pdb",
    "budesonide_rep3.dcd",
    "budesonide_rep3_solvated.pdb",
    "natural_top_rep1.dcd",
    "natural_top_rep1_solvated.pdb"
)

foreach ($file in $filesOther) {
    Write-Host "Downloading $file from other GPU..." -ForegroundColor Yellow
    scp -P 10619 "root@ssh4.vast.ai:/workspace/nod2-screening-data/PHASE_6/trajectories/$file" "$localDir\$file"
}

Write-Host ""
Write-Host "=== Downloading logs ===" -ForegroundColor Cyan
scp -P 38103 "root@ssh4.vast.ai:/workspace/nod2-screening-data/PHASE_6/logs/*" "C:\Users\vasud\nod2-screening-data\PHASE_6\logs\"
scp -P 10619 "root@ssh4.vast.ai:/workspace/nod2-screening-data/PHASE_6/logs/*" "C:\Users\vasud\nod2-screening-data\PHASE_6\logs\"

Write-Host ""
Write-Host "=== DOWNLOAD COMPLETE ===" -ForegroundColor Green
Write-Host "Files saved to: $localDir"
