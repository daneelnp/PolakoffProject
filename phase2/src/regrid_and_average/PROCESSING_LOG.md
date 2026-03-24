# Processing Log — Regrid and Annual Average
**Date run:** 2026-03-23
**Script:** `src/regrid_and_average/regrid_and_average.jl`

---

## Input Files

| File | Variable | Original Grid | Time Dimension | Units |
|------|----------|---------------|----------------|-------|
| `TS_Vcmax25.nc` | `variable` → `vcmax25` | 0.5° (lon=720, lat=360) | z=1..35 (1982–2016) | µmol CO₂ m⁻² s⁻¹ (assumed) |
| `TS_LMA_decidudous.nc` | `variable` → `lma_deciduous` | 0.5° (lon=720, lat=360) | z=1..35 (1982–2016) | g m⁻² (assumed) |
| `TS_LMA_evergreen.nc` | `variable` → `lma_evergreen` | 0.5° (lon=720, lat=360) | z=1..35 (1982–2016) | g m⁻² (assumed) |
| `clm_refl_tran_1m_weighted.nc` | `par_refl, par_tran, nir_refl, nir_tran` | 2.5°×1.875° (lon=144, lat=96) | ind=1..12 (Jan–Dec) | dimensionless |

**Note:** `TS_LMA_decidudous.nc` has a typo in the original filename ("decidudous"); this is preserved to match the file on disk.

---

## Target Grid (Wang & Braghiere 2025)

- **Longitude:** 1.25° to 358.75°, step 2.5°, N=144 (0–360 convention)
- **Latitude:** -89.0625° to 89.0625°, step 1.875°, N=96
- **Important:** The latitude step is 1.875°, not 2.5°. This is a reduced Gaussian-type grid. A fixed 5×5 block window would be incorrect in latitude.

---

## Processing Steps

### Dong et al. files (Vcmax25, LMA deciduous, LMA evergreen)

**Step 1 — Temporal average (3b)**
- Averaged over all 35 z-slices (years 1982–2016) to produce a single 2D climatological mean
- Missing values (`Union{Missing, Float32}` from NCDatasets) excluded before averaging
- A pixel with no valid data across all 35 years → NaN in the output

**Step 2 — Spatial regrid (3a)**
- Each Wang output cell is assigned the mean of all source (0.5°) cells whose centre falls within the Wang cell bounding box
- Lon bounding box: target_lon ± 1.25° (= half of 2.5° Wang cell width)
- Lat bounding box: target_lat ± 0.9375° (= half of 1.875° Wang cell height)
- This gives approximately 5 source cells in longitude and 3–4 in latitude per output cell (~15–20 source values)
- Cells where the entire block is NaN/missing → NaN in the output; otherwise finite values are averaged
- Source longitude was in -180/180 convention; Wang longitudes >180° were converted to -180/180 for index matching, but saved in 0–360 convention to match the source Wang file

### Wang & Braghiere optics

- Already on the target grid — no spatial regridding needed
- Fill value −1.0 (not NaN) replaced with NaN before averaging
- Annual mean computed as the mean over ind=1..12, ignoring NaN
- A pixel with no valid months → NaN in the output

---

## Output Files

| File | Variable(s) | Grid | Description |
|------|-------------|------|-------------|
| `dong_vcmax25_2p5deg_annual.nc` | `vcmax25` | lon=144, lat=96 | Vcmax25 annual mean 1982–2016 |
| `dong_lma_deciduous_2p5deg_annual.nc` | `lma_deciduous` | lon=144, lat=96 | LMA deciduous annual mean 1982–2016 |
| `dong_lma_evergreen_2p5deg_annual.nc` | `lma_evergreen` | lon=144, lat=96 | LMA evergreen annual mean 1982–2016 |
| `wang_optics_2p5deg_annual.nc` | `par_refl, par_tran, nir_refl, nir_tran` | lon=144, lat=96 | Annual mean leaf optics (Jan–Dec) |

All files use the same lon/lat coordinate arrays as the Wang & Braghiere source file.

---

## Output Statistics

| Output key | Variable | min | max | mean | N valid |
|------------|----------|-----|-----|------|---------|
| `lma_deciduous` | `lma_deciduous` | 0.5894 | 146.0983 | 47.1781 | 3861 |
| `lma_evergreen` | `lma_evergreen` | 41.8204 | 253.5716 | 132.8810 | 3861 |
| `vcmax25` | `vcmax25` | 12.9164 | 143.1638 | 53.5298 | 3861 |
| `wang_nir_refl` | `nir_refl` | 0.3384 | 0.3966 | 0.3709 | 6247 |
| `wang_nir_tran` | `nir_tran` | 0.4134 | 0.4781 | 0.4499 | 6247 |
| `wang_par_refl` | `par_refl` | 0.0493 | 0.3353 | 0.1146 | 6247 |
| `wang_par_tran` | `par_tran` | 0.0249 | 0.3685 | 0.1153 | 6247 |

---

## Assumptions

1. **Units** for Dong files are assumed (not verified from file metadata): Vcmax25 in µmol CO₂ m⁻² s⁻¹, LMA in g m⁻².
2. **Fill value** for Dong files is handled via Julia's `missing` (NCDatasets native); no manual fill detection needed.
3. **Fill value** for Wang file is −1.0 (confirmed from task specification; not stored as NaN in the source file).
4. **Temporal averaging** treats each year/month equally (unweighted mean). No area-weighting applied.
5. **Spatial regridding** uses an unweighted mean over source cells in the bounding box. No area-weighting applied (acceptable for the relatively small source cells).
6. **Longitude convention:** output files use 0–360 (matching the Wang source file) so all four processed files share identical coordinate arrays.
