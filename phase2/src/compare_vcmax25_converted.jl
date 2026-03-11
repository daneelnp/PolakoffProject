"""
    compare_vcmax25_converted.jl

Same as compare_vcmax25.jl but with CliMA values converted from
mol CO2 m⁻² s⁻¹ → µmol CO2 m⁻² s⁻¹ (×10⁶) so both datasets
are on the same scale.

Outputs:
  plot1_side_by_side_maps_converted.png
  plot2_histogram_converted.png

Usage:
    julia --project=. src/compare_vcmax25_converted.jl
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

clima_lon = Float64[]
clima_lat = Float64[]
clima_vcmax_monthly = Array{Float64,3}(undef, 0, 0, 0)

NCDataset(FILE_CLIMA, "r") do ds
    global clima_lon = ds["lon"][:]
    global clima_lat = ds["lat"][:]
    raw = ds["vcmax25"][:, :, :]           # (24, 360, 180)
    global clima_vcmax_monthly = raw[:, :, :]
end

println("  CliMA grid:  lon $(length(clima_lon))  lat $(length(clima_lat))")
println("  Monthly data shape: ", size(clima_vcmax_monthly))

# Annual mean, NaN-aware
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

# ── Convert mol → µmol ──
clima_annual .*= 1e6
println("  ✓ Converted CliMA: mol → µmol (×10⁶)")

valid_c = filter(isfinite, clima_annual)
println("  CliMA annual mean — valid: $(length(valid_c)), min: $(minimum(valid_c)), max: $(maximum(valid_c)), mean: $(round(mean(valid_c), digits=2)) µmol")


println("\n─── Loading Dong et al. 2023 data ───")

dong_lon = Float64[]
dong_lat = Float64[]
dong_annual = Array{Float64,2}(undef, 0, 0)

NCDataset(FILE_DONG, "r") do ds
    global dong_lon = ds["longitude"][:]
    global dong_lat = ds["latitude"][:]
    raw = ds["variable"][:, :, :]          # (720, 360, 35)

    nlon, nlat, nz = size(raw)
    data = fill(NaN, nlon, nlat, nz)
    for k in 1:nz, j in 1:nlat, i in 1:nlon
        v = raw[i, j, k]
        if !ismissing(v)
            data[i, j, k] = Float64(v)
        end
    end
    out = fill(NaN, nlon, nlat)
    for j in 1:nlat, i in 1:nlon
        vals = filter(isfinite, data[i, j, :])
        if !isempty(vals)
            out[i, j] = mean(vals)
        end
    end
    global dong_annual = out
end

valid_d = filter(isfinite, vec(dong_annual))
println("  Dong annual mean — valid: $(length(valid_d)), min: $(minimum(valid_d)), max: $(maximum(valid_d)), mean: $(round(mean(valid_d), digits=2)) µmol")

println("\n─── Data loading complete ───")


# ═════════════════════════════════════════════════════════════════════════════
#  STEP 2 — Plots (both in µmol, shared scales)
# ═════════════════════════════════════════════════════════════════════════════

# ─────────────────────────────────────────────────────────────────────────────
#  Plot 1 — Side-by-side global maps (shared colorscale)
# ─────────────────────────────────────────────────────────────────────────────
println("\n─── Plot 1: Side-by-side global maps (converted) ───")

# Each panel on its own colorscale (independent colorbars)
fig1 = Figure(size=(1600, 600))

ax1a = Axis(fig1[1, 1],
    title  = "CliMA Land — Annual Mean Vcmax25\n(mol→µmol converted)",
    xlabel = "Longitude",
    ylabel = "Latitude")
hm1a = heatmap!(ax1a, clima_lon, clima_lat, clima_annual;
    colormap = :viridis, nan_color = :gray90)
Colorbar(fig1[1, 2], hm1a, label="Vcmax25 [µmol CO₂ m⁻² s⁻¹]")

ax1b = Axis(fig1[1, 3],
    title  = "Dong et al. 2023 — Climatological Mean Vcmax25\n(1982–2016 average)",
    xlabel = "Longitude",
    ylabel = "Latitude")
hm1b = heatmap!(ax1b, dong_lon, dong_lat, dong_annual;
    colormap = :viridis, nan_color = :gray90)
Colorbar(fig1[1, 4], hm1b, label="Vcmax25 [µmol CO₂ m⁻² s⁻¹]")

save(joinpath(PLOT_DIR, "plot1_side_by_side_maps_converted.png"), fig1, px_per_unit=3)
println("  Saved: plot1_side_by_side_maps_converted.png")


# ─────────────────────────────────────────────────────────────────────────────
#  Plot 2 — Histogram (shared x-axis, both in µmol)
# ─────────────────────────────────────────────────────────────────────────────
println("\n─── Plot 2: Histogram (converted) ───")

vals_clima = filter(x -> isfinite(x) && x > 0, vec(clima_annual))
vals_dong  = filter(isfinite, vec(dong_annual))

fig2 = Figure(size=(1200, 500))

# Panel A: CliMA (now in µmol) — independent x-axis
ax2a = Axis(fig2[1, 1],
    xlabel = "Vcmax25 [µmol CO₂ m⁻² s⁻¹]",
    ylabel = "Count",
    title  = "CliMA Land — Vcmax25 distribution (converted)")
hist!(ax2a, vals_clima; bins = 50, color = (:blue, 0.6))

# Panel B: Dong — independent x-axis
ax2b = Axis(fig2[1, 2],
    xlabel = "Vcmax25 [µmol CO₂ m⁻² s⁻¹]",
    ylabel = "Count",
    title  = "Dong et al. — Vcmax25 distribution")
hist!(ax2b, Float64.(vals_dong); bins = 50, color = (:red, 0.6))

save(joinpath(PLOT_DIR, "plot2_histogram_converted.png"), fig2, px_per_unit=3)
println("  Saved: plot2_histogram_converted.png")

println("\n═══ All converted plots complete ═══")
