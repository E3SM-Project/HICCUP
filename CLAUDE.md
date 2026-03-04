# HICCUP - Claude Code Guide

## Project Overview

**HICCUP** (Hindcast Initial Condition Creation Utility/Processor) creates E3SM atmospheric initial condition files from observational/reanalysis data (ERA5, CAMS, NOAA). Key capabilities: horizontal/vertical regridding, surface T/P adjustment for new topography, SST/sea-ice processing, support for EAM and EAMxx target models.

## Directory Structure

```
hiccup/                    # Main Python package
  hiccup.py                # Public API: create_hiccup_data() factory
  hiccup_data_class.py     # Core base class (~3800 lines)
  hiccup_data_class_common.py        # Shared imports (__all__)
  hiccup_data_class_grid_methods.py  # Grid creation/manipulation
  hiccup_data_class_sstice_methods.py
  hiccup_data_class_memory_methods.py
  hiccup_data_class_timer_methods.py
  hiccup_data_subclass_ERA5.py       # ERA5 source handler
  hiccup_data_subclass_CAMS.py
  hiccup_data_subclass_NOAA.py
  hiccup_data_subclass_EAM.py
  hiccup_state_adjustment.py         # Surface P/T adjustment algorithms
  hiccup_utilities.py                # tcolor, run_cmd, version checks
  hiccup_constants.py                # Physical constants
test_scripts/              # Unit and integration tests
  unit_test_all.py         # Main test runner
  unit_test_data_class.py
  unit_test_state_adjustment.py
  unit_test_memory_methods.py
  end2end_test.*.py
template_hiccup_scripts/   # Example user scripts
test_data/                 # Low-res test input files
files_vert/                # Vertical coordinate files
```

## Installation

```bash
conda create --name hiccup_env -c conda-forge xarray dask pandas numpy scipy \
  netcdf4 hdf5 cdsapi tempest-remap "nco>=5.3.1"
source activate hiccup_env
pip install ./
```

**System dependencies:** NCO (>=5.3.1), TempestRemap

## Running Tests

Tests must be run from the `test_scripts/` directory:

```bash
cd test_scripts
python unit_test_all.py           # all tests
python unit_test_data_class.py    # specific module
```

New test modules must be added to both `unit_test_all.py` (import + `suite_list.append`).

## Coding Conventions

- **Indentation:** 2 spaces
- **Naming:** snake_case for functions/methods/variables; classes use lowercase with underscores (e.g., `hiccup_data`, `ERA5`)
- **Constants:** ALL_CAPS (e.g., `Rdair`, `gravit`)
- **Section separators:** `# ----------...----------` (78 dashes)
- **Imports:** wildcard `from hiccup.hiccup_data_class_common import *` used internally; all shared imports go through `hiccup_data_class_common.py`
- **Docstrings:** short, one-line descriptions in triple quotes
- **Error handling:** `ValueError`, `KeyError`, `OSError` with descriptive messages using `tcolor`
- **Tests:** use `unittest`, mock `psutil.Process` for memory tests, call `print_timer` at end of each test method

## Architecture

- **Factory pattern:** `create_hiccup_data(src_data_name=...)` returns appropriate subclass
- **Mixin pattern:** large method groups split into separate files (grid, sstice, memory, timer methods) and composed into main class
- **Subclass selection:** `classmethod is_name_for()` on each subclass
- **Variable mapping:** dictionaries `atm_var_name_dict`, `sfc_var_name_dict` per data source
- **Environment:** `HDF5_USE_FILE_LOCKING="FALSE"` set in `hiccup.py`

## Key Dependencies

| Package | Purpose |
|---------|---------|
| xarray | Primary data manipulation |
| numpy | Numerical operations |
| pandas | Time coordinates |
| scipy | Pole-filling |
| psutil | Memory monitoring |
| NCO | Shell-based NetCDF operations |
| TempestRemap | Unstructured grid regridding |
| cdsapi | ERA5 downloads (requires `~/.cdsapirc`) |
