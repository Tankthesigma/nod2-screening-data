# Auto-download DCD files as they complete
# Run this in background - it checks every 5 minutes

$localDir = "C:\Users\vasud\nod2-screening-data\PHASE_6\trajectories"
$checkInterval = 300  # 5 minutes

# Track what we've downloaded
$downloaded = @{}

while ($true) {
    Write-Host "`n[$(Get-Date -Format 'HH:mm:ss')] Checking for completed files..." -ForegroundColor Cyan

    # Check 5090 instance
    Write-Host "  Checking 5090..." -ForegroundColor Yellow
    $files5090 = ssh -o StrictHostKeyChecking=no -o ConnectTimeout=10 -p 38103 root@ssh4.vast.ai "ls -la /workspace/nod2-screening-data/PHASE_6/trajectories/*.dcd 2>/dev/null" 2>$null

    foreach ($line in $files5090 -split "`n") {
        if ($line -match "(\S+\.dcd)$") {
            $filename = $matches[1]
            $size = ($line -split "\s+")[4]

            # Check if file is complete (>4GB for most, or not growing)
            if ($size -gt 4000000000 -and -not $downloaded[$filename]) {
                Write-Host "  Downloading $filename ($([math]::Round($size/1GB, 2)) GB)..." -ForegroundColor Green
                scp -P 38103 "root@ssh4.vast.ai:/workspace/nod2-screening-data/PHASE_6/trajectories/$filename" "$localDir\$filename"
                scp -P 38103 "root@ssh4.vast.ai:/workspace/nod2-screening-data/PHASE_6/trajectories/$($filename -replace '\.dcd$','_solvated.pdb')" "$localDir\" 2>$null
                $downloaded[$filename] = $true
                Write-Host "  Done: $filename" -ForegroundColor Green
            }
        }
    }

    # Check other GPU instance
    Write-Host "  Checking other GPU..." -ForegroundColor Yellow
    $filesOther = ssh -o StrictHostKeyChecking=no -o ConnectTimeout=10 -p 10619 root@ssh4.vast.ai "ls -la /workspace/nod2-screening-data/PHASE_6/trajectories/*.dcd 2>/dev/null" 2>$null

    foreach ($line in $filesOther -split "`n") {
        if ($line -match "(\S+\.dcd)$") {
            $filename = $matches[1]
            $size = ($line -split "\s+")[4]

            if ($size -gt 4000000000 -and -not $downloaded[$filename]) {
                Write-Host "  Downloading $filename ($([math]::Round($size/1GB, 2)) GB)..." -ForegroundColor Green
                scp -P 10619 "root@ssh4.vast.ai:/workspace/nod2-screening-data/PHASE_6/trajectories/$filename" "$localDir\$filename"
                scp -P 10619 "root@ssh4.vast.ai:/workspace/nod2-screening-data/PHASE_6/trajectories/$($filename -replace '\.dcd$','_solvated.pdb')" "$localDir\" 2>$null
                $downloaded[$filename] = $true
                Write-Host "  Done: $filename" -ForegroundColor Green
            }
        }
    }

    Write-Host "  Downloaded so far: $($downloaded.Count) files"
    Write-Host "  Sleeping $checkInterval seconds..." -ForegroundColor Gray
    Start-Sleep -Seconds $checkInterval
}
