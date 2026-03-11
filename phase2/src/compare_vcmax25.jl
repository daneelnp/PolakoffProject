"""
    compare_vcmax25.jl

Compare Vcmax25 between CliMA Land simulation output and
Dong et al. 2023 optimality-based predictions.

  • vcmax25_1M_average.nc  — CliMA Land (1°, 24 months: 2008-03 → 2010-02)
  • TS_Vcmax25.nc          — Dong et al. (0.5°, z=1..35 = years 1982–2016)

Decisions (confirmed by user):
  1. CliMA units left as-is: mol CO2 m⁻² s⁻¹ (flagged for later confirmation)
  2. Dong: average over all 35 z-slices → single climatological annual mean
  3. CliMA: use only last 12 months (indices 13–24, i.e. 2009-03 → 2010-02)

Usage:
    julia --project=. src/compare_vcmax25.jl
"""

using NCDatasets
using Statistics
using CairoMakie

# ── Paths ────────────────────────────────────────────────────────────────────
const DATA_DIR = joinpath(@__DIR__, "..", "data")
const PLOT_DIR = joinpath(@__DIR__, "..", "output_plots")
mkpath(PLOT_DIR)

const FILE_CLIMA = joinpath(DATA_DIR, "vcmax25_1M_average.nc")
const FILE_DONG  = joinpath(DATA_DIR, "TS_Vcmax25.nc")

# ═════════════════════════════════════════════════════════════════════════════
#  STEP 1 — Load and prepare data
# ═════════════════════════════════════════════════════════════════════════════

println("─── Loading CliMA Land data ───")

# Variable name: "vcmax25"
# Dims: (time=24, lon=360, lat=180) — 1° grid
# Time: months 1–12 = 2008-03..2009-02, months 13–24 = 2009-03..2010-02
# Units: mol CO2 m⁻² s⁻¹ (not converted; flagged for confirmation)

clima_lon = Float64[]
clima_lat = Float64[]
clima_vcmax_monthly = Array{Float64,3}(undef, 0, 0, 0)  # placeholder

NCDataset(FILE_CLIMA, "r") do ds
    global clima_lon = ds["lon"][:]                  # (360,)
    global clima_lat = ds["lat"][:]                  # (180,)
    # Variable name may need update ──▼
    raw = ds["vcmax25"][:, :, :]                     # (time=24, lon=360, lat=180)

    # Use only last 12 months (indices 13–24) to skip spin-up
    # global clima_vcmax_monthly = raw[13:24, :, :]    # (12, 360, 180)

    # Use all months (1–24) for now, including spin-up ---▼
    global clima_vcmax_monthly = raw[:, :, :]        # (24, 360, 180)
end

println("  CliMA grid:  lon $(length(clima_lon))  lat $(length(clima_lat))")
println("  Monthly data shape (last 12 mo): ", size(clima_vcmax_monthly))

# Annual mean over 12 months, treating NaN as missing
function nanmean_dim1(A::AbstractArray{Float64,3})
    nt, nx, ny = size(A)
    out = fill(NaN, nx, ny)
    for j in 1:ny, i in 1:nx
        vals = filter(isfinite, A[:, i, j])
        if !isempty(vals)
            out[i, j] = mean(vals)
        end
    end
    return out
end

clima_annual = nanmean_dim1(clima_vcmax_monthly)   # (360, 180)
println("  CliMA annual mean shape: ", size(clima_annual))
valid_c = filter(isfinite, clima_annual)
println("  CliMA annual mean — valid: $(length(valid_c)), min: $(minimum(valid_c)), max: $(maximum(valid_c)), mean: $(round(mean(valid_c), digits=6))")


println("\n─── Loading Dong et al. 2023 data ───")

# Variable name: "variable"
# Dims: (longitude=720, latitude=360, z=35) — 0.5° grid
# z = 1..35 → years 1982–2016
# Units: assumed µmol CO2 m⁻² s⁻¹ (no units attr in file)

dong_lon = Float64[]
dong_lat = Float64[]
dong_annual = Array{Float64,2}(undef, 0, 0)

NCDataset(FILE_DONG, "r") do ds
    global dong_lon = ds["longitude"][:]             # (720,)
    global dong_lat = ds["latitude"][:]              # (360,)
    # Variable name may need update ──▼
    raw = ds["variable"][:, :, :]                    # (720, 360, 35)

    # Average over all 35 z-slices (years) → climatological mean
    # Replace missing with NaN first
    nlon, nlat, nz = size(raw)
    data = fill(NaN, nlon, nlat, nz)
    for k in 1:nz, j in 1:nlat, i in 1:nlon
        v = raw[i, j, k]
        if !ismissing(v)
            data[i, j, k] = Float64(v)
        end
    end
    # NaN-aware mean over z (dim 3)
    out = fill(NaN, nlon, nlat)
    for j in 1:nlat, i in 1:nlon
        vals = filter(isfinite, data[i, j, :])
        if !isempty(vals)
            out[i, j] = mean(vals)
        end
    end
    global dong_annual = out                         # (720, 360)
end

println("  Dong grid:  lon $(length(dong_lon))  lat $(length(dong_lat))")
println("  Dong annual mean shape: ", size(dong_annual))
valid_d = filter(isfinite, dong_annual)
println("  Dong annual mean — valid: $(length(valid_d)), min: $(minimum(valid_d)), max: $(maximum(valid_d)), mean: $(round(mean(valid_d), digits=4))")

# NOTE: Grids differ in resolution (CliMA = 1°, Dong = 0.5°).
#       No regridding is performed per user instructions.

println("\n─── Data loading complete ───")


# ═════════════════════════════════════════════════════════════════════════════
#  STEP 2 — Plots
# ═════════════════════════════════════════════════════════════════════════════

# ─────────────────────────────────────────────────────────────────────────────
#  Plot 1 — Side-by-side global maps of annual mean Vcmax25
# ─────────────────────────────────────────────────────────────────────────────
println("\n─── Plot 1: Side-by-side global maps ───")

# Determine shared color range from both datasets
# CliMA: mol CO2 m⁻² s⁻¹; Dong: µmol CO2 m⁻² s⁻¹ — different units!
# Plot each on its own natural scale with its own colorbar.
# (A shared scale is not meaningful until units are reconciled.)

fig1 = Figure(size=(1600, 600))

# CliMA panel (left)
ax1a = Axis(fig1[1, 1],
    title  = "CliMA Land — Annual Mean Vcmax25\n(2009-03 → 2010-02)",
    xlabel = "Longitude",
    ylabel = "Latitude")
# dims: clima_annual is (lon=360, lat=180), heatmap wants (x, y)
hm1a = heatmap!(ax1a, clima_lon, clima_lat, clima_annual;
    colormap = :viridis, nan_color = :gray90)
Colorbar(fig1[1, 2], hm1a, label="Vcmax25 [mol CO₂ m⁻² s⁻¹]")

# Dong panel (right)
ax1b = Axis(fig1[1, 3],
    title  = "Dong et al. 2023 — Climatological Mean Vcmax25\n(1982–2016 average)",
    xlabel = "Longitude",
    ylabel = "Latitude")
hm1b = heatmap!(ax1b, dong_lon, dong_lat, dong_annual;
    colormap = :viridis, nan_color = :gray90)
Colorbar(fig1[1, 4], hm1b, label="Vcmax25 [µmol CO₂ m⁻² s⁻¹] (assumed)")

save(joinpath(PLOT_DIR, "plot1_side_by_side_maps.png"), fig1, px_per_unit=3)
println("  Saved: plot1_side_by_side_maps.png")


# ─────────────────────────────────────────────────────────────────────────────
#  Plot 2 — Histogram
# ─────────────────────────────────────────────────────────────────────────────
println("\n─── Plot 2: Histogram ───")

vals_clima = filter(isfinite, vec(clima_annual))
vals_dong  = filter(isfinite, vec(dong_annual))

fig4 = Figure(size=(1200, 500))

# Panel A: CliMA
ax4a = Axis(fig4[1, 1],
    xlabel = "Vcmax25 [mol CO₂ m⁻² s⁻¹]",
    ylabel = "Count",
    title  = "CliMA Land — Vcmax25 distribution")
hist!(ax4a, vals_clima; bins = 50, color = (:blue, 0.6))

# Panel B: Dong
ax4b = Axis(fig4[1, 2],
    xlabel = "Vcmax25 [µmol CO₂ m⁻² s⁻¹] (assumed)",
    ylabel = "Count",
    title  = "Dong et al. — Vcmax25 distribution")
hist!(ax4b, Float64.(vals_dong); bins = 50, color = (:red, 0.6))

save(joinpath(PLOT_DIR, "plot2_histogram.png"), fig4, px_per_unit=3)
println("  Saved: plot2_histogram.png")

println("\n─── All plots complete ───")