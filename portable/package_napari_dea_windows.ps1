# portable napari environment with napari-spatialdata and napari-spatialdata-dea
# University Medicine Greifswald, 2025
# please cite napari and spatialdata if you find this script useful

Param(
    [string]$AppName = "napari-dea-windows"
)

$ErrorActionPreference = "Stop"

$EnvPrefix = ".\tmp_env"
$AppDir    = ".\$AppName"
$EnvDir    = Join-Path $AppDir "env"
$BinDir    = Join-Path $AppDir "bin"
$ZipName   = "$AppName-win64.zip"

Write-Host "[0/6] Cleaning previous build..." -ForegroundColor Cyan
if (Test-Path $EnvPrefix) { Remove-Item -Recurse -Force $EnvPrefix }
if (Test-Path $AppDir)    { Remove-Item -Recurse -Force $AppDir }
if (Test-Path $ZipName)   { Remove-Item -Force $ZipName }
if (Test-Path "env-win64.zip") { Remove-Item -Force "env-win64.zip" }

Write-Host "[1/6] Creating conda env at $EnvPrefix ..." -ForegroundColor Cyan
# Requires mamba/conda on PATH from e.g. Miniforge/Anaconda Prompt
mamba create -y -p $EnvPrefix -c conda-forge `
    python=3.11 `
    pip `

Write-Host "[2/6] Installing your plugin into that env..." -ForegroundColor Cyan
& "$EnvPrefix\python.exe" -m pip install --no-warn-script-location `
    git+https://github.com/Yasas1994/napari-spatialdata-dea.git@main `
    "napari[all]" `
    spatialdata `
    napari-spatialdata `
    spatialdata-io `
    spatialdata-plot

Write-Host "[3/6] Ensuring conda-pack is available..." -ForegroundColor Cyan
mamba install -y -c conda-forge conda-pack

Write-Host "[4/6] Packing environment to env-win64.zip ..." -ForegroundColor Cyan
conda-pack -p $EnvPrefix -o env-win64.zip

Write-Host "[4.5/6] Creating app layout in $AppDir ..." -ForegroundColor Cyan
New-Item -ItemType Directory -Path $AppDir, $EnvDir, $BinDir -Force | Out-Null

Write-Host "[4.6/6] Unpacking env into $EnvDir ..." -ForegroundColor Cyan
Expand-Archive -Path "env-win64.zip" -DestinationPath $EnvDir -Force

# NOTE: we *could* run conda-unpack here, but we'll also call it in the launcher
#       so it relocates on first run if needed.

Write-Host "[5/6] Creating launcher script bin\napari-dea.bat ..." -ForegroundColor Cyan
$LauncherPath = Join-Path $BinDir "napari-dea.bat"
@'
@echo off
setlocal enabledelayedexpansion

REM Resolve APP_DIR = parent of this script directory
set "SCRIPT_DIR=%~dp0"
REM Remove trailing backslash
if "%SCRIPT_DIR:~-1%"=="\" set "SCRIPT_DIR=%SCRIPT_DIR:~0,-1%"
set "APP_DIR=%SCRIPT_DIR%\.."
set "ENV_DIR=%APP_DIR%\env"

REM Activate the packed env
call "%ENV_DIR%\Scripts\activate.bat"

REM Try to run conda-unpack once (idempotent)
if exist "%ENV_DIR%\Scripts\conda-unpack.exe" (
    "%ENV_DIR%\Scripts\conda-unpack.exe" 2>NUL
)

REM Launch napari with any extra args passed through
python -m napari %*
endlocal
'@ | Out-File -FilePath $LauncherPath -Encoding ASCII -Force

Write-Host "[5.5/6] Creating README.txt ..." -ForegroundColor Cyan
$ReadmePath = Join-Path $AppDir "README.txt"
@"
$AppName (Windows, win-64)

This directory contains:
- A relocatable conda environment in .\env
- A launcher script in .\bin\napari-dea.bat

Requirements:
- 64-bit Windows
- A normal GPU/graphics stack (napari uses Qt/ OpenGL)

Usage:
1. Unzip $ZipName
2. Navigate to bin Folder
3. Double click on napari-dea to launch napari

You do NOT need conda installed on the target machine.
Everything needed is bundled in .\env.
"@ | Out-File -FilePath $ReadmePath -Encoding UTF8 -Force

Write-Host "[6/6] Creating distributable zip $ZipName ..." -ForegroundColor Cyan
Compress-Archive -Path $AppDir -DestinationPath $ZipName -Force

Write-Host "Done!"
Write-Host "Created: $ZipName" -ForegroundColor Green

# powershell -ExecutionPolicy Bypass -File .\package_napari_dea_windows.ps1