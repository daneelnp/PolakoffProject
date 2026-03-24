using NCDatasets
using Statistics
using Printf
using Dates

# ─── Paths ────────────────────────────────────────────────────────────────────
const DATA_DIR  = joinpath(@__DIR__, "..", "..", "data")
const OUT_DIR   = joinpath(DATA_DIR, "processed")
const LOG_PATH  = joinpath(@__DIR__, "PROCESSING_LOG.md")
mkpath(OUT_DIR)

# Collect stats for each output variable so we can write the log at the end
log_stats = Dict{String, NamedTuple}()

# ─── STEP 0: Load the target grid from Wang & Braghiere ───────────────────────
# This file defines the authoritative output grid.
# lon: 1.25–358.75, step 2.5°, N=144  (0–360 convention)
# lat: -89.0625–89.0625, step 1.875°, N=96
wang_path = joinpath(DATA_DIR, "clm_refl_tran_1m_weighted.nc")

ds_wang = NCDataset(wang_path, "r")
wang_lon_raw = Float64.(ds_wang["lon"][:])   # 0–360 convention; saved as-is in output
wang_lat     = Float64.(ds_wang["lat"][:])   # south-to-north
close(ds_wang)

nlon = length(wang_lon_raw)   # 144
nlat = length(wang_lat)       # 96

# Convert Wang lons to -180/180 so they align with the Dong source grid (-179.75–179.75)
wang_lon = copy(wang_lon_raw)
wang_lon[wang_lon .> 180.0] .-= 360.0

# Half-widths of each Wang cell (used to define bounding box for block averaging)
dlon_half = (wang_lon_raw[2] - wang_lon_raw[1]) / 2.0   # 1.25°
dlat_half = abs(wang_lat[2]  - wang_lat[1])  / 2.0      # 0.9375°

println("=" ^ 65)
println("Target grid (Wang & Braghiere)")
println("=" ^ 65)
println("  lon (0–360): $(wang_lon_raw[1]) to $(wang_lon_raw[end])  N=$nlon  step=$(wang_lon_raw[2]-wang_lon_raw[1])°")
println("  lat:         $(wang_lat[1]) to $(wang_lat[end])  N=$nlat  step=$(abs(wang_lat[2]-wang_lat[1]))°")
println()

# ─── Helper: time-average a 3-D (lon×lat×time) array ignoring missing/NaN ────
# Works on Union{Missing, Float32} arrays from NCDatasets.
# Returns a 2-D Float64 array; cells with no valid data → NaN.
function time_avg(src::AbstractArray)
    nlon_s, nlat_s, nt = size(src)
    out = fill(NaN, nlon_s, nlat_s)
    for i in 1:nlon_s, j in 1:nlat_s
        s = 0.0; n = 0
        for t in 1:nt
            v = src[i, j, t]
            if !ismissing(v) && isfinite(Float64(v))
                s += Float64(v)
                n += 1
            end
        end
        n > 0 && (out[i, j] = s / n)
    end
    return out
end

# ─── Helper: block-average a 2-D (lon×lat) Float64 source onto the Wang grid ─
# For each Wang cell, finds all source cells whose centre falls within the
# Wang cell bounding box (±dlon_half in lon, ±dlat_half in lat) and averages
# the finite values.
# NOTE: source and target lons must both be in the same convention (-180/180).
function spatial_regrid(src_2d::Matrix{Float64},
                        src_lon::Vector{Float64},
                        src_lat::Vector{Float64},
                        wang_lon::Vector{Float64},
                        wang_lat::Vector{Float64},
                        dlon_half::Float64,
                        dlat_half::Float64)
    out = fill(NaN, length(wang_lon), length(wang_lat))
    for i in 1:length(wang_lon)
        # Boolean mask of source lons within this Wang cell's lon band
        lon_mask = abs.(src_lon .- wang_lon[i]) .< dlon_half
        for j in 1:length(wang_lat)
            lat_mask = abs.(src_lat .- wang_lat[j]) .< dlat_half
            block = src_2d[lon_mask, lat_mask]   # ~5×4 = ~20 values
            s = 0.0; n = 0
            for v in block
                if isfinite(v)
                    s += v
                    n += 1
                end
            end
            n > 0 && (out[i, j] = s / n)
        end
    end
    return out
end

# ─── Helper: save a single-variable 2-D result to NetCDF ─────────────────────
function save_nc(path, varname, units, long_name, data, wang_lon_raw, wang_lat)
    NCDataset(path, "c") do ds
        defDim(ds, "lon", size(data, 1))
        defDim(ds, "lat", size(data, 2))

        vlon = defVar(ds, "lon", Float64, ("lon",))
        vlon[:] = wang_lon_raw
        vlon.attrib["units"]     = "degrees_east"
        vlon.attrib["long_name"] = "longitude"
        vlon.attrib["note"]      = "0–360 convention matching Wang & Braghiere source file"

        vlat = defVar(ds, "lat", Float64, ("lat",))
        vlat[:] = wang_lat
        vlat.attrib["units"]     = "degrees_north"
        vlat.attrib["long_name"] = "latitude"

        vd = defVar(ds, varname, Float64, ("lon", "lat"))
        vd.attrib["units"]       = units
        vd.attrib["long_name"]   = long_name
        vd.attrib["description"] = "Annual climatological mean regridded to Wang & Braghiere 2.5° grid"
        vd[:] = data
    end
end

# ─── Helper: print & record stats ────────────────────────────────────────────
function report(label, varname, data)
    valid = filter(isfinite, vec(data))
    mn  = minimum(valid)
    mx  = maximum(valid)
    mu  = mean(valid)
    nv  = length(valid)
    @printf("  %-16s  min=%8.4f  max=%8.4f  mean=%8.4f  N_valid=%d\n",
            varname, mn, mx, mu, nv)
    log_stats[label] = (varname=varname, min=mn, max=mx, mean=mu, n_valid=nv)
end

# ═══════════════════════════════════════════════════════════════════════════════
# TASK 3a + 3b: Regrid and time-average Dong et al. files
# ═══════════════════════════════════════════════════════════════════════════════
# File list: (source filename, output varname, output nc name, units, log key)
dong_files = [
    ("TS_Vcmax25.nc",        "vcmax25",       "dong_vcmax25_2p5deg_annual.nc",        "umol CO2 m-2 s-1", "vcmax25"),
    ("TS_LMA_decidudous.nc", "lma_deciduous", "dong_lma_deciduous_2p5deg_annual.nc",  "g m-2",            "lma_deciduous"),
    ("TS_LMA_evergreen.nc",  "lma_evergreen", "dong_lma_evergreen_2p5deg_annual.nc",  "g m-2",            "lma_evergreen"),
]
# Note: "TS_LMA_decidudous.nc" is the actual filename on disk (typo in source).

println("=" ^ 65)
println("TASK 3a+3b — Dong files: regrid 0.5°→2.5° + annual mean")
println("=" ^ 65)

for (fname, varname, outname, units, log_key) in dong_files
    println("\nProcessing: $fname")

    # ── Load source data ──────────────────────────────────────────────────────
    src_path = joinpath(DATA_DIR, fname)
    ds = NCDataset(src_path, "r")
    src_lon  = Float64.(ds["longitude"][:])    # -179.75 to  179.75, step 0.5°, N=720
    src_lat  = Float64.(ds["latitude"][:])     #   89.75 to  -89.75, step -0.5° (N→S)
    src_data = ds["variable"][:, :, :]         # lon×lat×z, Union{Missing, Float32}
    close(ds)

    println("  Source grid: lon $(src_lon[1]) to $(src_lon[end])  lat $(src_lat[1]) to $(src_lat[end])")
    println("  Step 1: averaging over $(size(src_data, 3)) years (z=1..35, 1982–2016)...")

    # ── Step 1: time-average 720×360×35 → 720×360 ────────────────────────────
    src_2d = time_avg(src_data)

    # ── Step 2: spatial block-average 720×360 → 144×96 ───────────────────────
    println("  Step 2: spatial regrid 0.5°→2.5° (bounds-based block average)...")
    result = spatial_regrid(src_2d, src_lon, src_lat, wang_lon, wang_lat, dlon_half, dlat_half)

    # ── Verify output grid ────────────────────────────────────────────────────
    # (the output array is indexed by wang_lon/wang_lat, so the grid matches by construction)
    println("  Output grid: lon $(wang_lon_raw[1]) to $(wang_lon_raw[end])  lat $(wang_lat[1]) to $(wang_lat[end])")
    report(log_key, varname, result)

    # ── Save ──────────────────────────────────────────────────────────────────
    out_path = joinpath(OUT_DIR, outname)
    save_nc(out_path, varname, units,
            "$(varname) annual mean 1982-2016 regridded to Wang 2.5° grid",
            result, wang_lon_raw, wang_lat)
    println("  Saved: data/processed/$outname")
end

# ═══════════════════════════════════════════════════════════════════════════════
# TASK 3b: Annual-average Wang & Braghiere optics (already on target grid)
# ═══════════════════════════════════════════════════════════════════════════════
println()
println("=" ^ 65)
println("TASK 3b — Wang & Braghiere: annual mean over 12 months")
println("=" ^ 65)

const WANG_FILL = -1.0   # fill value in Wang file (not NaN — must mask manually)

optical_vars = [
    ("par_refl", "PAR reflectance",    "dimensionless"),
    ("par_tran", "PAR transmittance",  "dimensionless"),
    ("nir_refl", "NIR reflectance",    "dimensionless"),
    ("nir_tran", "NIR transmittance",  "dimensionless"),
]

wang_annual = Dict{String, Matrix{Float64}}()   # stores results for saving

ds_wang = NCDataset(wang_path, "r")
println("\nProcessing: clm_refl_tran_1m_weighted.nc")

for (vname, long_name, units) in optical_vars
    raw = Float64.(ds_wang[vname][:, :, :])   # lon×lat×ind (144×96×12)

    # Replace fill value (-1.0) with NaN before averaging
    raw[raw .≈ WANG_FILL] .= NaN

    # Annual mean: average over 12 months, ignoring NaN
    result = fill(NaN, nlon, nlat)
    for i in 1:nlon, j in 1:nlat
        vals = filter(isfinite, raw[i, j, :])
        isempty(vals) || (result[i, j] = mean(vals))
    end

    wang_annual[vname] = result
    report("wang_$vname", vname, result)
end

close(ds_wang)

# Save all four optical variables into one NetCDF file
wang_out_path = joinpath(OUT_DIR, "wang_optics_2p5deg_annual.nc")
NCDataset(wang_out_path, "c") do ds
    defDim(ds, "lon", nlon)
    defDim(ds, "lat", nlat)

    vlon = defVar(ds, "lon", Float64, ("lon",))
    vlon[:] = wang_lon_raw
    vlon.attrib["units"]     = "degrees_east"
    vlon.attrib["long_name"] = "longitude"
    vlon.attrib["note"]      = "0–360 convention"

    vlat = defVar(ds, "lat", Float64, ("lat",))
    vlat[:] = wang_lat
    vlat.attrib["units"]     = "degrees_north"
    vlat.attrib["long_name"] = "latitude"

    for (vname, long_name, units) in optical_vars
        vd = defVar(ds, vname, Float64, ("lon", "lat"))
        vd.attrib["units"]       = units
        vd.attrib["long_name"]   = long_name
        vd.attrib["description"] = "Annual climatological mean (avg over ind=1..12) from Wang & Braghiere 2025"
        vd[:] = wang_annual[vname]
    end
end
println("  Saved: data/processed/wang_optics_2p5deg_annual.nc")

# ═══════════════════════════════════════════════════════════════════════════════
# Write PROCESSING_LOG.md
# ═══════════════════════════════════════════════════════════════════════════════
run_date = Dates.format(now(), "yyyy-mm-dd")

open(LOG_PATH, "w") do io
    write(io, """
# Processing Log — Regrid and Annual Average
**Date run:** $run_date
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
""")

    for (key, s) in sort(collect(log_stats), by=first)
        @printf(io, "| `%s` | `%s` | %.4f | %.4f | %.4f | %d |\n",
                key, s.varname, s.min, s.max, s.mean, s.n_valid)
    end

    write(io, """

---

## Assumptions

1. **Units** for Dong files are assumed (not verified from file metadata): Vcmax25 in µmol CO₂ m⁻² s⁻¹, LMA in g m⁻².
2. **Fill value** for Dong files is handled via Julia's `missing` (NCDatasets native); no manual fill detection needed.
3. **Fill value** for Wang file is −1.0 (confirmed from task specification; not stored as NaN in the source file).
4. **Temporal averaging** treats each year/month equally (unweighted mean). No area-weighting applied.
5. **Spatial regridding** uses an unweighted mean over source cells in the bounding box. No area-weighting applied (acceptable for the relatively small source cells).
6. **Longitude convention:** output files use 0–360 (matching the Wang source file) so all four processed files share identical coordinate arrays.
""")
end

println()
println("Processing log written to: src/regrid_and_average/PROCESSING_LOG.md")
println("\nAll done.")
