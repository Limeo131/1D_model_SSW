# Input Data

The 1-D model requires three input NetCDF files. These are not included in the
repository due to their size and licensing. Instructions for obtaining or
generating each file are provided below.

## Required Files

| File | Size | Description |
|------|------|-------------|
| `fawa_interp.nc` | ~1 MB | Interpolated FAWA climatology from ERA5 |
| `urad_below32_fit_above32_HMshear_501_1000.nc` | ~1.9 GB | Radiative-equilibrium wind (1000yr) |
| `init_501_1000.nc` | ~6 MB | Initial conditions (1000yr) |
| `forcing_total_1000_ctl.nc` | ~4 MB | Control stochastic forcing (1000yr) |

## ERA5 Source Data

The model is calibrated against ERA5 reanalysis (1979–2021).

- **Access**: https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5
- **Variables needed**: zonal-mean zonal wind, EP flux, FAWA at 60°N
- **Levels**: 1000–1 hPa (37 levels)
- **Temporal resolution**: 6-hourly
- **Spatial**: 60°N zonal mean

## Generating Input Files

1. **fawa_interp.nc**: Generated from ERA5 EP-flux data. See `scripts/preprocess/`.
2. **urad**: Run `scripts/preprocess/construct_urad.ipynb` using ERA5-derived profiles.
3. **init**: Run `scripts/preprocess/generate_init_conditions.ipynb` from ERA5 early-winter states.
4. **forcing**: Run `scripts/preprocess/generate_stochastic_forcing.ipynb` to create 1000yr IAAFT surrogate.

## Environment Variables

The model reads input file paths from environment variables:

```bash
export FAWA_INTERP_FILE="data/fawa_interp.nc"
export URAD_FILE="data/urad_below32_fit_above32_HMshear_501_1000.nc"
export INIT_FILE="data/init_501_1000.nc"
```

If unset, the model defaults to `data/` relative paths.
