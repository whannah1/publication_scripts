# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Publication scripts for the **SOHIP** study — analysis of June 2023 stratospheric gravity wave observations from the ISS, compared against SCREAM RRM (Regionally Refined Model) simulations at 256x2 and 256x3 refinement levels. Observation regions: Patagonia, South/Equatorial Indian Ocean, South/Southeast Pacific.

## Running Scripts

All scripts are run directly from the `code/` directory or the repo root:

```bash
# Production figures
python code/FXX-curtain-v2.py        # vertical curtain plots (current version)
python code/FXX-map.py               # 2D map plots
python code/FXX-map-precip.py        # precipitation maps
python code/FXX-hovmoller-v1.py      # time-height Hovmoller diagrams
python code/FXX-stockwell-v1.py      # Stockwell transform gravity wave analysis

# Data preparation (run before plotting if curtain files don't exist)
python code/mask_data.v1.py          # spatial subsetting of model output
python code/calc.curtain.v1.py       # interpolate model data to observation paths
python code/calc.curtain.v2.py       # same but uses pre-masked (reduced) data

# Utility
python convert_IMERG.hourly.py       # convert IMERG HDF5 → NetCDF4
```

Figures are saved to `figs/`. Long runs are typically launched with `nohup`.

## Architecture

### Data Flow

```
SCREAM model output (10-min NetCDF snapshots, unstructured ncol grid)
  → [optional] mask_data.v1.py         # crop to 20°×20° box, save *reduced.nc
  → calc.curtain.v{1,2}.py             # interpolate to ISS limb-scan paths
  → curtain_data/*.curtain.*.nc        # pre-computed curtain files
  → FXX-curtain-v2.py / FXX-hovmoller / FXX-stockwell  → figs/
```

Maps bypass the curtain step and read model output directly.

### Key Paths

| Resource | Path |
|---|---|
| Model output | `/global/cfs/cdirs/m4842/whannah/cases/<case>/run/` |
| Reduced (masked) model output | `.../cases/<case>/data_reduced/` |
| Pre-computed curtain files | `/global/cfs/cdirs/m4842/whannah/curtain_data/` |
| Grid descriptor files | `/global/cfs/cdirs/m4842/whannah/files_grid/` |
| ERA5 validation data | `/global/cfs/projectdirs/m4842/whannah/ERA5/` |
| IMERG precipitation | `/global/cfs/cdirs/m4842/whannah/IMERG/` |

### `sohip_methods.py` — Shared Utility Library

All plotting and calc scripts do `from sohip_methods import *`. Key functions:

- **`add_case()` / `add_var()`** — register case and variable configurations into module-level dicts used by all scripts
- **`calculate_path()`** — define the 1200 km limb-scan path from ISS position to tangent point
- **`find_closest_cells()` / `find_path_ncol_wgt()`** — nearest-neighbor search + inverse-distance weights on the unstructured SCREAM grid
- **`interpolate_to_path()` / `interpolate_to_path_numba()`** — interpolate 3D model fields along observation paths (numba-accelerated inner loop)
- **`reduce_file_list_to_target_time()`** — select model snapshot files nearest to an observation time
- **`calc_gw_ep()`** — compute gravity wave potential energy via Stockwell transform: background removal (polynomial fit, order configurable), S-transform, buoyancy frequency N, returns (g/N)² × power (multiply by 0.5 for J/kg)
- **`get_coi_mask()`** — cone-of-influence mask for wavelet edge effects
- **`get_grid_file()`** — maps case name strings to SCRIP grid files

### Script Configuration Pattern

Each plotting/calc script configures cases and variables at the top, then runs the analysis:

```python
from sohip_methods import *

add_case('256x2-sw-ind-v1', n='grid_name', xtime='2023-06-12 19:00',
         xlat=-49.0, xlon=60.0,   # ISS position
         tlat=-52.0, tlon=62.0,   # tangent point
         slat=-49.0, slon=60.0)   # scan center

add_var('T_mid', vstr='Temperature', htype='i', ...)

path_len_km = 1200   # total path length
path_spc_km = 2      # along-path point spacing
path_ncells = 2      # nearest model cells for inv-distance averaging
```

Multiple cases are looped over; to focus on one region, comment out the others.

## Code Style Notes

Follow the separator convention in `CLAUDE.md` (global):
- Major section dividers: dashed line at column 100, description on next line
- Minor section dividers: dashed line at column 80, description on next line
- Never embed text inside dashes

## Case Name Conventions

`{resolution}-{region}-v{N}` e.g. `256x2-sw-ind-v1`

Regions: `ptgnia` (Patagonia), `eq-ind` (Equatorial Indian), `sc-ind` (South Central Indian), `sc-pac` (South Central Pacific), `se-pac` (Southeast Pacific), `sw-ind` (Southwest Indian).

Script versions: `v0`/`v1` = older/experimental, `v2` = current production. `code_test/` holds archived experiments.
