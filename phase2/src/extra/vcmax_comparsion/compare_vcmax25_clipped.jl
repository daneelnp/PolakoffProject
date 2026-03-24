"""
    compare_vcmax25_clipped.jl

Same as converted version but with CliMA map colorbar clipped to 0–150 µmol
to match Dong's range and reveal spatial structure.

Outputs:
  plot1_side_by_side_maps_clipped.png
  plot2_histogram_clipped.png

Usage:
    julia --project=. src/compare_vcmax25_clipped.jl
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

clima_annual = nanmean_dim1(clima_vcmax_monthly)

# Convert mol → µmol
clima_annual .*= 1e6
println("  ✓ Converted CliMA: mol → µmol (×10⁶)")

valid_c = filter(x -> isfinite(x) && x > 0, vec(clima_annual))
println("  CliMA annual mean (positive) — n=$(length(valid_c)), min=$(minimum(valid_c)), max=$(maximum(valid_c)), mean=$(round(mean(valid_c), digits=2))")

# ── Percentiles ──
sorted_c = sort(valid_c)
for p in [25, 50, 75, 95, 99]
    idx = max(1, round(Int, p / 100 * length(sorted_c)))
    println("  p$p: $(round(sorted_c[idx], digits=2)) µmol")
end
n_above_150 = count(x -> x > 150, valid_c)
println("  Values > 150 µmol: $n_above_150 / $(length(valid_c)) ($(round(100 * n_above_150 / length(valid_c), digits=2))%)")


println("\n─── Loading Dong et al. 2023 data ───")

dong_lon = Float64[]
dong_lat = Float64[]
dong_annual = Array{Float64,2}(undef, 0, 0)

NCDataset(FILE_DONG, "r") do ds
    global dong_lon = ds["longitude"][:]
    global dong_lat = ds["latitude"][:]
    raw = ds["variable"][:, :, :]

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
println("  Dong annual mean — n=$(length(valid_d)), min=$(minimum(valid_d)), max=$(maximum(valid_d)), mean=$(round(mean(valid_d), digits=2))")

sorted_d = sort(valid_d)
for p in [25, 50, 75, 95, 99]
    idx = max(1, round(Int, p / 100 * length(sorted_d)))
    println("  p$p: $(round(sorted_d[idx], digits=2)) µmol")
end

println("\n─── Data loading complete ───")


# ═════════════════════════════════════════════════════════════════════════════
#  STEP 2 — Plots
# ═════════════════════════════════════════════════════════════════════════════

# ─────────────────────────────────────────────────────────────────────────────
#  Plot 1 — Side-by-side maps, CliMA clipped to 0–150
# ─────────────────────────────────────────────────────────────────────────────
println("\n─── Plot 1: Side-by-side maps (clipped) ───")

fig1 = Figure(size=(1600, 600))

ax1a = Axis(fig1[1, 1],
    title  = "CliMA Land — Annual Mean Vcmax25 (mol→µmol)\n⚠ Colorbar clipped to 0–150; values above 150 are saturated",
    xlabel = "Longitude",
    ylabel = "Latitude")
hm1a = heatmap!(ax1a, clima_lon, clima_lat, clima_annual;
    colormap = :viridis, colorrange = (0, 150), nan_color = :gray90)
Colorbar(fig1[1, 2], hm1a, label="Vcmax25 [µmol CO₂ m⁻² s⁻¹]")

ax1b = Axis(fig1[1, 3],
    title  = "Dong et al. 2023 — Climatological Mean Vcmax25\n(1982–2016 average)",
    xlabel = "Longitude",
    ylabel = "Latitude")
hm1b = heatmap!(ax1b, dong_lon, dong_lat, dong_annual;
    colormap = :viridis, colorrange = (0, 150), nan_color = :gray90)
Colorbar(fig1[1, 4], hm1b, label="Vcmax25 [µmol CO₂ m⁻² s⁻¹]")

save(joinpath(PLOT_DIR, "plot1_side_by_side_maps_clipped.png"), fig1, px_per_unit=3)
println("  Saved: plot1_side_by_side_maps_clipped.png")


# ─────────────────────────────────────────────────────────────────────────────
#  Plot 2 — Histogram (shared x-axis, both in µmol)
# ─────────────────────────────────────────────────────────────────────────────
println("\n─── Plot 2: Histogram (clipped) ───")

vals_clima_hist = filter(x -> isfinite(x) && x > 0, vec(clima_annual))
vals_dong_hist  = filter(isfinite, vec(dong_annual))

fig2 = Figure(size=(1200, 500))

ax2a = Axis(fig2[1, 1],
    xlabel = "Vcmax25 [µmol CO₂ m⁻² s⁻¹]",
    ylabel = "Count",
    title  = "CliMA Land — Vcmax25 distribution (converted)")
hist!(ax2a, vals_clima_hist; bins = 50, color = (:blue, 0.6))

ax2b = Axis(fig2[1, 2],
    xlabel = "Vcmax25 [µmol CO₂ m⁻² s⁻¹]",
    ylabel = "Count",
    title  = "Dong et al. — Vcmax25 distribution")
hist!(ax2b, Float64.(vals_dong_hist); bins = 50, color = (:red, 0.6))

save(joinpath(PLOT_DIR, "plot2_histogram_clipped.png"), fig2, px_per_unit=3)
println("  Saved: plot2_histogram_clipped.png")

println("\n═══ All clipped plots complete ═══")
