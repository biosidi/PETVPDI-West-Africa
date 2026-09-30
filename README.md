Processing code associated with the study:
"Development and Application of the Potential Evapotranspiration
Vegetation Precipitation Drought Index (PETVPDI) in West Africa
using Remote Sensing Datasets"

Study period: 2001–2022

Main PETVPDI inputs:
- Precipitation (GPM)
- Potential evapotranspiration (MOD16A2GF)
- NDVI (MOD13A3)

The repository contains scripts for Data preprocessing,
PETVPDI calculation, validation, and statistical analyses.

PETVPDI = √([(〖〖NPET〗_max-NPET)〗^2+(NNDVI-N〖NDVI〗_min )^2+(NP-〖NP〗_min )^2])
