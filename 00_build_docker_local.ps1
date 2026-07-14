# =============================================================================
# 00_build_docker_local.ps1
# Build TR_Bear_Connectivity Docker image (ENMTML pipeline) and save as tar.gz.
# Run from PowerShell:  .\00_build_docker_local.ps1
# =============================================================================

$ErrorActionPreference = 'Stop'

$SCRIPT_DIR  = $PSScriptRoot
$OUT_DIR     = (Join-Path $PSScriptRoot "programs")
$DOCKERFILE  = Join-Path $SCRIPT_DIR "00_Dockerfile_trbear"
$IMAGE_TAG   = "trbear:latest"
$TAR_PATH    = Join-Path $OUT_DIR "trbear.tar"
$TARGZ_PATH  = Join-Path $OUT_DIR "trbear.tar.gz"
$BUILD_LOG   = Join-Path $SCRIPT_DIR "docker_build.log"
$BUILD_ERR   = Join-Path $SCRIPT_DIR "docker_build.err"

if (-not (Test-Path $OUT_DIR)) {
    New-Item -ItemType Directory -Force -Path $OUT_DIR | Out-Null
}

Write-Host "[INFO] Checking Docker availability..." -ForegroundColor Cyan
docker --version
if ($LASTEXITCODE -ne 0) { throw "Docker not available. Start Docker Desktop." }

# Reset logs
Remove-Item -Force -ErrorAction Ignore $BUILD_LOG, $BUILD_ERR

Write-Host "[INFO] Building image $IMAGE_TAG" -ForegroundColor Cyan
Write-Host "[INFO]   stdout -> $BUILD_LOG" -ForegroundColor Cyan
Write-Host "[INFO]   stderr -> $BUILD_ERR" -ForegroundColor Cyan
Write-Host "[INFO] BuildKit progress goes to stderr. To watch live, open a 2nd PowerShell:" -ForegroundColor DarkGray
Write-Host "       Get-Content $BUILD_ERR -Wait -Tail 20" -ForegroundColor DarkGray

$env:DOCKER_BUILDKIT = "1"

$dockerArgs = @(
    "build",
    "--progress=plain",
    "-f", $DOCKERFILE,
    "-t", $IMAGE_TAG,
    $SCRIPT_DIR
)

$proc = Start-Process -FilePath "docker" `
    -ArgumentList $dockerArgs `
    -NoNewWindow -Wait -PassThru `
    -RedirectStandardOutput $BUILD_LOG `
    -RedirectStandardError  $BUILD_ERR

if ($proc.ExitCode -ne 0) {
    Write-Host ""
    Write-Host "[FAIL] Docker build failed (exit $($proc.ExitCode))." -ForegroundColor Red
    Write-Host "[FAIL] Last 150 lines of stderr below:" -ForegroundColor Red
    Write-Host ("-" * 80) -ForegroundColor DarkGray
    if (Test-Path $BUILD_ERR) {
        Get-Content $BUILD_ERR -Tail 150
    }
    Write-Host ("-" * 80) -ForegroundColor DarkGray
    Write-Host "[FAIL] Full logs:" -ForegroundColor Red
    Write-Host "       $BUILD_LOG" -ForegroundColor Red
    Write-Host "       $BUILD_ERR" -ForegroundColor Red
    throw "Docker build failed."
}

Write-Host "[INFO] Build succeeded. Sanity-check ENMTML..." -ForegroundColor Green

# Sanity test inside the container (single-line R command to avoid PS parser issues)
$rCheck = "cat(sprintf('%-12s %s\n', 'ENMTML', packageVersion('ENMTML'))); cat(sprintf('%-12s %s\n', 'terra', packageVersion('terra'))); cat(sprintf('%-12s %s\n', 'dismo', packageVersion('dismo'))); cat(sprintf('%-12s %s\n', 'maxnet', packageVersion('maxnet'))); cat(sprintf('%-12s %s\n', 'rgdal', packageVersion('rgdal')))"
docker run --rm $IMAGE_TAG R --no-save -e $rCheck

docker run --rm $IMAGE_TAG python -c "import numpy, rasterio, geopandas, networkx; print('py deps OK')"
docker run --rm $IMAGE_TAG ls /opt/programs/UNICOR | Select-Object -First 5
docker images $IMAGE_TAG

Write-Host "[INFO] Saving image to tar..." -ForegroundColor Cyan
if (Test-Path $TAR_PATH)   { Remove-Item $TAR_PATH -Force }
if (Test-Path $TARGZ_PATH) { Remove-Item $TARGZ_PATH -Force }
docker save -o $TAR_PATH $IMAGE_TAG
if ($LASTEXITCODE -ne 0) { throw "docker save failed." }

Write-Host "[INFO] Compressing tar to tar.gz..." -ForegroundColor Cyan
Add-Type -AssemblyName System.IO.Compression.FileSystem
$inStream  = [System.IO.File]::OpenRead($TAR_PATH)
$outStream = [System.IO.File]::Create($TARGZ_PATH)
$gzStream  = New-Object System.IO.Compression.GZipStream($outStream, [System.IO.Compression.CompressionMode]::Compress)
$inStream.CopyTo($gzStream)
$gzStream.Close()
$outStream.Close()
$inStream.Close()
Remove-Item $TAR_PATH -Force

$gzSize = (Get-Item $TARGZ_PATH).Length / 1MB
$gzMsg  = "[INFO] Saved: $TARGZ_PATH  ({0:N1} MB)" -f $gzSize
Write-Host $gzMsg -ForegroundColor Green

Write-Host "" -ForegroundColor Yellow
Write-Host "[INFO] Next steps:" -ForegroundColor Yellow
Write-Host "  1. scp '$TARGZ_PATH' <user>@<cluster-login-host>:<project-dir>/" -ForegroundColor Yellow
Write-Host "  2. ssh to cluster, cd to scripts dir, sbatch 00_convert_to_sif.slurm" -ForegroundColor Yellow
