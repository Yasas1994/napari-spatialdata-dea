## Packaging napari with spatialdata dependencies for portability

### Requirements
1. conda or mamba environment with conda-pack

```
pip install conda-pack
```
2. run the the appropriate script for your OS

**for linux and osx**


```
chmod + ./package_napari_dea_unix.sh & ./package_napari_dea_unix.sh
```
**for windows**

```
powershell -ExecutionPolicy Bypass -File .\package_napari_dea_windows.ps1
```