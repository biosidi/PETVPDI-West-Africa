## Processing and PETVPDI Calculation Codes

This repository contains the MATLAB scripts used for data preprocessing, spatial and temporal harmonization, normalization, and PETVPDI calculation over West Africa for the period 2001–2022.

### 1. `Compute_GPM.m` — GPM Precipitation Preprocessing

This MATLAB script processes the GPM precipitation dataset. It performs the main preprocessing steps, including extraction of the West African study area, spatial resampling to the target grid of 1 km, and pixel-wise normalization of precipitation.

Before running the script:
- Update the input and output directory paths according to the location of the datasets on your computer.
- Install Generic Mapping Tools (GMT) and add the GMT installation directory to the MATLAB path.

### 2. `Compute_NDVI.m` — NDVI Normalization

This MATLAB script processes and normalizes the MOD13A3 NDVI dataset. The NDVI values are normalized independently for each pixel using the temporal minimum and maximum values over the study period (2001–2022).

### 3. `Compute_PET.m` — PET Preprocessing, Temporal Aggregation, Resampling, and Normalization

This MATLAB script processes the MOD16A2GF potential evapotranspiration (PET) dataset. The original 8-day PET composites are quality controlled and aggregated to monthly PET values. The monthly PET data are then resampled to the common analysis grid and normalized using pixel-wise temporal minimum and maximum values over 2001–2022.

### 4. `PETVPDI_Calculation.m` — PETVPDI Calculation

This MATLAB script calculates the Potential Evapotranspiration–Vegetation–Precipitation Drought Index (PETVPDI) over West Africa from January 2001 to December 2022.

The calculation integrates three normalized variables:
- Normalized potential evapotranspiration (NPET)
- Normalized vegetation condition (NNDVI)
- Normalized precipitation (NP)

PETVPDI is calculated as:

PETVPDI = sqrt[(1 − NPET)^2 + NNDVI^2 + NP^2]

The resulting monthly PETVPDI maps are exported as GeoTIFF files.

## Important Note

Users should modify all input and output directory paths in the MATLAB scripts according to their local data-storage structure before running the codes. The scripts should be executed in the following order:

`Compute_GPM.m` → `Compute_NDVI.m` → `Compute_PET.m` → `PETVPDI_Calculation.m`