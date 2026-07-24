# 1-D Semi-Empirical Model of Sudden Stratospheric Warmings

A computationally efficient one-dimensional model of the winter polar stratosphere
for studying the frequency of sudden stratospheric warmings (SSWs).

## Associated Paper

> **What Controls the Frequency of Sudden Stratospheric Warmings? Insights from a Semi-Empirical 1-D Model**
>
> Siming Liu and Noboru Nakamura
>
> Department of the Geophysical Sciences, University of Chicago

## Overview

This model represents the coupled evolution of finite-amplitude wave activity (FAWA)
and zonal-mean zonal wind in the Northern Hemisphere winter stratosphere at 60°N.
Key features:

- Vertical domain: 10–110 km (kmax=501, Δz=200 m, Δt=120 s)
- Parameters calibrated against 42 years of ERA5 reanalysis (1979–2021)
- Simulates ~1000 years per hour, enabling robust SSW frequency estimates
- Sensitivity experiments for wave forcing amplitude and radiative-equilibrium wind

## Quick Start

### Requirements

- Fortran compiler (gfortran recommended)
- NetCDF-Fortran library
- Python 3.9+ with NumPy, SciPy, xarray, netCDF4, Matplotlib

### Build

```bash
make
```

Or manually:
```bash
gfortran -O3 -ffree-line-length-none -I$(nc-config --includedir) \
  src/fortran/1d_model.f90 $(nf-config --flibs) -o d.out
```

### Run

```bash
# Single 1000-year control run (using default parameters from Table 1)
./d.out 0.233 0.5 5 7 16 0.5 2.2e-3 4.0 0.3 control forcing.nc .
```

Command-line arguments:
```
beta s_tau damp1_days u_thres flux_extra alpha const1 ppow aexp tag forcing_file output_dir
```

### Environment Variables for Input Data

```bash
export FAWA_INTERP_FILE="data/fawa_interp.nc"
export URAD_FILE="data/urad_below32_fit_above32_HMshear_501_1000.nc"
export INIT_FILE="data/init_501_1000.nc"
```

## Repository Structure

```
├── src/fortran/          # Fortran model source
├── scripts/
│   ├── preprocess/       # Input data preparation (uRAD, forcing, init)
│   ├── run/              # Compile and execute scripts
│   ├── figures/          # Manuscript figure generation
│   └── analysis/        # SSW detection and diagnostics (TODO)
├── calibration/          # Gray-box parameter optimization
├── config/               # Reference parameters and experiment configs
├── data/                 # Input data (not tracked; see data/README.md)
├── docs/                 # Documentation and provenance
├── manuscript/           # Submitted manuscript (read-only reference)
└── output/               # Model output (not tracked)
```

## Experiments in the Paper

1. **ERA5 Calibration** (Sec 2.4–2.5): 42-winter integrations with ERA5 forcing
2. **Stochastic Control** (Sec 3): 1000-year ensemble (10 members)
3. **Forcing Trends** (Sec 3.1): CTL, LINEAR, SATUR growth experiments
4. **uRAD Perturbations** (Sec 3.2): Top, Top+Bottom, Top−Bottom
5. **Parameter Sweep** (Sec 3.3): Forcing amplitude × uRAD perturbation

## Reference Parameters (Table 1)

| Parameter | Value | Description |
|-----------|-------|-------------|
| α | 0.5 | Coupling coefficient |
| C | 2.2×10⁻³ | E-P flux propagation coefficient |
| χ_A | 4.0 | FAWA damping curvature exponent |
| a_exp | 0.30 day⁻¹ | Upper-level damping amplitude |
| β | 0.233 | Nonlinear mixing strength |
| s_τ | 0.5 | Mixing timescale scaling |
| d₁ | 5 days | Mesospheric relaxation timescale |
| u* | 7 m/s | Weak-wind threshold |
| Δu | 16 m/s | Weak-wind enhancement |

## Data Availability

ERA5 reanalysis: https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5

See `data/README.md` for instructions on preparing input files.

## Acknowledgments

This work is supported by NSF Award No. 2501292.

## License

See LICENSE file.
