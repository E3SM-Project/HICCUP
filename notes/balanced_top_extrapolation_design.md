# Design brief: balanced top-of-model extrapolation in HICCUP

Status: design proposal, not yet implemented. Written 2026-09-22.

This brief is meant to be self-contained, so that someone (or some agent) with no
prior context can pick up the implementation.

## 1. Background: why this is needed

HICCUP's vertical remap fills target levels that lie **above the source data's top**
by clamping: every extrapolated level gets a copy of the source's top value, for T
and winds alike. When the target grid extends well above the source, this produces an
initial state that is badly out of thermal-wind balance, and the model adjusts
violently in the first hours.

**Concrete case that motivated this (E3SM v4 development, Sep 2026):**

- The target grid is the new E3SM v4 L128 (`vertical_coordinates_L128v4_c20260820.nc`),
  with its top interface at **0.499 hPa**.
- The IC `screami_ne256np4L128v4_ifs-20200120_20260825.nc` was made by *vertically
  regridding an older L128 IC* whose top interface is **2.255 hPa**. The target's top
  5 levels (midpoints 0.59, 0.82, 1.12, 1.51, 2.04 hPa, ~1.5 scale heights) are above
  the source.
- In the resulting file, T and u on those 5 levels are **exact copies** of the source's
  top level, in every column. For example, at ~64N 60E: T = 228.1 K and u = 74.1 m/s
  on all five levels.
- Mechanism: each column keeps its own T, so horizontal T gradients persist unchanged
  through 1.5 scale heights, while u stops shearing. Thermal wind then demands tens of
  m/s more shear than exists.
- Observed result (EAMxx, ne256, all gravity-wave drag off): the lid-level max wind
  goes from the IC's ~108 m/s to 124 -> 136 -> 158 -> 172 m/s over the first 3 hours,
  then oscillates with a ~12-13 h period, the inertial period at 60-65N. It reached
  >200 m/s on day 1. This is classic geostrophic adjustment from an unbalanced state.

## 2. Current HICCUP code (as of Sep 2026)

- `hiccup/hiccup_vertical_remap.py`: `remap_vertical_py()` is the entry point. The
  per-column kernel `_interp_column()` interpolates in log-pressure; `extrap='constant'`
  (the default) is `np.interp`'s clamp, and `'linear'` extends the end slope. Columns
  are processed with `xr.apply_ufunc(..., vectorize=True, dask='parallelized')`, so the
  kernel is strictly column-local.
- `hiccup_data_class.py`: `remap_vertical()` wraps `remap_vertical_py()`;
  `remap_vertical_multifile()` loops over files. The NCO path (`remap_vertical_nco`)
  is legacy.
- Source layouts are auto-detected in `_compute_input_pressure()`: EAM/EAMxx hybrid
  (hyam/hybm/P0), ECMWF IFS hybrid (hyam in Pa + lnsp), or pure pressure levels.
- Workflow order in the templates (e.g.
  `template_hiccup_scripts/create_EAMxx_IC_from_ERA5-NOAA.py`): **horizontal remap
  first, then vertical**. So the vertical remap runs on the target (usually
  unstructured GLL np4) grid.
- `get_hindcast_data.ERA5.py` downloads ERA5 **pressure levels**, whose top is
  **1 hPa**. So even a fresh ERA5-based v4 IC extrapolates 2 levels (~0.5 scale height).
- Repo conventions (see `CLAUDE.md`): 2-space indent, snake_case, 78-dash section
  separators, errors via `ValueError`/`KeyError`/`OSError` with `tcolor`, `unittest`
  tests in `test_scripts/` that must be registered in `unit_test_all.py`.
- **EAMxx variable names:** EAMxx ICs *recently* switched from a packed
  `horiz_winds(time,ncol,dim2,lev)` to separate `U`/`V`. Temperature is `T_mid` in
  EAMxx and `T` in EAM. Older IC files still carry `horiz_winds`; support it only if
  backward compatibility is explicitly wanted.

## 3. Goal and non-goals

**Goal:** when target levels lie above the source top, fill T, U and V there so that
(a) T is climatologically reasonable for the date and latitude, (b) winds are in
thermal-wind balance with T, and (c) everything is continuous with the source at the
junction.

**Non-goals:** changing interpolation *within* the source range; changing
bottom/surface extrapolation (handled elsewhere by HICCUP's state adjustment); an
accurate mesospheric state (above ~1 hPa the model spins up its own circulation within
days, so the aim there is "no shock", not accuracy).

## 4. Design

### 4.1 Climatology data (committed to the repo)

- A **monthly, zonal-mean temperature climatology** on latitude x pressure, 12 months.
  - Source: ERA5 monthly means on **all 37 pressure levels** (1000-1 hPa), e.g.
    1991-2020, zonally averaged, at ~1 degree latitude.
  - Extended above 1 hPa with **NRLMSIS 2.x** (e.g. `pymsis`, zonal mean at mid-month,
    moderate constant F10.7/Ap) up to ~1e-4 hPa. Convert MSIS altitude to pressure by
    hydrostatic integration anchored to ERA5 at 1 hPa, and blend over an overlap range
    so T is continuous at the join.
- A matching **climatological zonal-wind shear field** `U_clim(phi, p, month)`, computed
  offline from `T_clim` by thermal wind (or gradient wind, see section 7), stored
  relative to a reference level. At runtime only differences
  `U_clim(p) - U_clim(p_a)` are used.
- **Size is negligible:** ~181 lat x ~60 levels x 12 months x 2 fields ~ 1 MB.
- Commit the **generation script** alongside the data (ERA5 years, MSIS version,
  parameters). MSIS/pymsis is a **generation-time dependency only**, never a runtime one.
- Keep all levels, not just the upper stratosphere. This covers future sources with low
  tops (anywhere above the tropopause) at no extra cost.

### 4.2 Thermal wind (pressure coordinates; y = a*phi, x = a*cos(phi)*lambda, R = R_d)

- `du/dln(p) = (R/f) * dT/dy`
- `dv/dln(p) = -(R/f) * dT/dx`
- Offline climatological shear:
  `U_clim(phi,p) - U_clim(phi,p_ref) = (R/(f a)) * integral_{ln p_ref}^{ln p} dT_clim/dphi dln(p')`
- **Equatorial taper:** thermal wind fails as f -> 0. Taper all shear terms smoothly to
  zero within |phi| <~ 10-15 degrees. Make the taper latitude a parameter and document
  the choice.

### 4.3 Runtime algorithm (inside `remap_vertical_py`, after normal interpolation)

Everything happens at the vertical remap stage, on whatever grid the data is on.
**No intermediate files and no separate pass over the data.**

1. Determine the **anchor pressure** `p_a`: the source top by default, or a few source
   levels below the top for model-generated sources, whose top levels are contaminated
   by the source model's sponge layer. For each column, find the target levels with
   `p < p_a`. If there are none anywhere, it's a no-op.
2. Interpolate all variables as today.
3. Get the climatology for the IC date: interpolate between months (section 4.4) and in
   latitude to each column's `lat`.
4. Define the blend weight `w(p)`: 1 at `p_a`, decaying smoothly to 0 over a blend depth
   `H_b` in ln p (default ~0.5). Suggested form: `w = exp(-(ln p_a - ln p)/H_b)` or a
   cosine taper to exactly 0. Also define `W(p) = integral_{ln p}^{ln p_a} w dln(p')`
   (1D, >= 0 above the anchor).
5. For target levels above the anchor:
   - **T:** `T = T_clim(phi,p) + T'_a * w(p)`, where `T'_a = T_src(p_a) - T_clim(phi,p_a)`.
   - **U (option A, baseline):** `U = U_src(p_a) + [U_clim(phi,p) - U_clim(phi,p_a)]`.
   - **V (option A):** `V = V_src(p_a)`, since the zonal-mean climatology has no
     meridional shear.
   - **Option B (exact; add later as a separate, additive term):** the anomaly is
     separable (`T'_a(phi,lambda) * w(p)`), so its thermal-wind shear is too:
     - `U += -(R/(f a)) * dT'_a/dphi * W(p)`
     - `V += +(R/(f a cos(phi))) * dT'_a/dlambda * W(p)`
     - This needs **one** horizontal gradient field (grad T'_a on the anchor level, 2D
       only). On the unstructured target grid, use `uxarray` with the control-volume
       SCRIP grid for GLL np4, and verify it at cube edges and poles. Apply the same
       equatorial taper.
6. All other variables keep the current behavior (clamp).
7. **Always print an extrapolation report**, whether or not the feature is on: the
   number of target levels above the source top (or anchor), the depth in scale heights
   (ln-p units), and the min/max source top pressure. This alone would have caught the
   v4 IC problem.
8. Raise an error if the target top is above the climatology's top.

Why not an abrupt switch (`w` = 0 immediately above the anchor)? Daily T anomalies in
the polar vortex region can be 20-40 K, so an abrupt switch puts a jump of that size
into one layer and produces a large spike or collapse in N^2, which the model's
gravity-wave schemes also use. A smooth blend avoids that.

### 4.4 Dates and calendars

- Take the date from the source's `time` variable. Handle **cftime calendars, including
  `noleap` and year 0001**, which are common in model-generated ICs.
- Interpolate linearly between adjacent monthly climatologies by day of year.
- Provide an override (`clim_date=`) for files with missing or meaningless dates, such
  as perpetual or spun-up states.

### 4.5 Proposed API additions to `remap_vertical_py` / `remap_vertical`

- `extrap_top='constant' | 'climatology'`. The default stays `'constant'` initially
  (opt-in rollout); flip it after validation.
- `anchor_pressure=None` or `anchor_levels_below_top=0` (use one or the other).
- `blend_depth=0.5` (ln-p units)
- `clim_file=None` (defaults to the file shipped in the repo)
- `clim_date=None`
- `equator_taper_lat=15.0`
- `balance_anomaly=False` (option B)
- Variable-name mapping for T/U/V per target model (EAM: `T`, `U`, `V`; EAMxx:
  `T_mid`, `U`, `V`), with a clear error if T/U/V can't be identified when
  `extrap_top='climatology'`.

## 5. Use cases (all run through the same vertical-remap code path)

| use case | steps | what's new |
|---|---|---|
| ERA5 -> EAM/EAMxx (horizontal + vertical) | download ERA5 PL -> horizontal remap -> vertical remap | vertical remap with `extrap_top='climatology'`; date from ERA5 time; anchor = source top |
| Model IC -> new vertical grid only (e.g. L128 -> L128v4) | vertical remap | same option, plus `anchor_levels_below_top` (e.g. 2-3) to skip the source sponge levels; date from file or override |
| Model IC -> new horizontal and vertical grid | horizontal remap -> vertical remap | same as above |
| Target top <= source top (the common model-to-model case) | unchanged | no-op; the report says 0 levels extrapolated |
| **Repair mode** | run the correction on an already-generated IC | treat the old source top (e.g. 2.255 hPa for the v4 IC above) as the anchor and rewrite the levels above it; lets a bad IC be fixed without regenerating it |

A complementary option is to download **ERA5/IFS model-level data** (137 levels, top
0.01 hPa). That avoids top extrapolation entirely for ERA5-based ICs, and
`_compute_input_pressure()` already supports the IFS hybrid/`lnsp` layout. It doesn't
help model-to-model remaps.

## 6. Validation

**Unit tests** (`test_scripts/`, registered in `unit_test_all.py`):

- No-op when the target top <= the source top: bit-for-bit with the current output.
- Continuity at the anchor: T, U, V equal the source values there.
- Far above the anchor (w ~ 0), T equals `T_clim`.
- On a synthetic zonal-mean T with a known analytic latitude gradient, the computed U
  shear matches thermal wind; the equatorial taper behaves as specified.
- Option B on a synthetic separable anomaly: U/V increments match the analytic formula.
- Date handling: `noleap`, year 0001, month interpolation, and the override.
- The report reports the correct level count and depth.
- An error when the target top is above the climatology top.

**End-to-end acceptance test:** regenerate or repair the ne256 L128v4 IC and run EAMxx
for 1 day with hourly `U`, `V`, `T` output on the top ~20 levels.

- **Before** (current IC): the lid-level max wind jumps from ~108 m/s to 170+ within
  3 h and oscillates with a ~13 h period.
- **Pass:** no inertial-period ringing; lid-level winds stay near their initial values
  over the first day. Later drift due to model physics (e.g. missing gravity-wave drag)
  is a separate issue.

## 7. Open decisions and risks

- **Equatorial treatment:** the taper latitude and shape.
- **Thermal vs gradient wind.** For strong polar jets the curvature term matters (at
  100 m/s and 60N, `u^2 tan(phi) / a` ~ 20% of `f u`). Gradient-wind balance is still 1D
  and still offline, so consider it for the stored `U_clim`, especially above 1 hPa.
- **Default anchor for model sources:** how many top levels to distrust.
- **Blend depth** default.
- **Option B:** whether the uxarray gradient on the GLL/SCRIP grid is accurate enough at
  cube edges and poles. Option A is the baseline, and option B is a purely additive term
  on top of it.
- **Rollout:** opt-in first, validated on both an ERA5 IC and a model-to-model remap,
  then change the default.
- **Other variables** (ozone, tracers) keep the clamp. A later extension could use the
  same climatology machinery for ozone, but that's out of scope here.
