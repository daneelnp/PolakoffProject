"""
    compare_lma_overlap.jl

Load the Dong et al. TS_LMA_evergreen.nc and TS_LMA_decidudous.nc files,
plot them on maps to see whether their spatial coverage overlaps or is disjoint.

Outputs:
  plot1_lma_side_by_side.png    — individual maps
  plot2_lma_overlap.png         — overlap / classification map

Usage:
    julia --project=. src/lma_comparison/compare_lma_overlap.jl
"""

using NCDatasets
using Statistics
using CairoMakie

# ── Paths ────────────────────────────────────────────────────────────────────
const DATA_DIR = joinpath(@__DIR__, "..", "..", "data")
const PLOT_DIR = joinpath(@__DIR__, "..", "..", "output_plots")
mkpath(PLOT_DIR)

const FILE_EV  = joinpath(DATA_DIR, "TS_LMA_evergreen.nc")
const FILE_DC  = joinpath(DATA_DIR, "TS_LMA_decidudous.nc")

# ═════════════════════════════════════════════════════════════════════════════
#  STEP 0 — Quick inspection
# ═════════════════════════════════════════════════════════════════════════════
println("═══ Inspecting files ═══")
for (label, path) in [("Evergreen", FILE_EV), ("Deciduous", FILE_DC)]
    NCDataset(path, "r") do ds
        println("\n─── $label: $(basename(path)) ───")
        println("  Dimensions: ", ds.dim)
        for vname in keys(ds)
            v = ds[vname]
            println("  Variable '$vname': dims=$(dimnames(v)), size=$(size(v)), type=$(eltype(v))")
            atts = Dict(k => v.attrib[k] for k in keys(v.attrib))
            if !isempty(atts)
                println("    Attributes: $atts")
            end
        end
    end
end

# ═════════════════════════════════════════════════════════════════════════════
#  STEP 1 — Load data
# ═════════════════════════════════════════════════════════════════════════════
println("\n═══ Loading data ═══")

function load_dong_lma(path)
    lon = Float64[]
    lat = Float64[]
    annual = Array{Float64,2}(undef, 0, 0)
    NCDataset(path, "r") do ds
        lon = ds["longitude"][:]
        lat = ds["latitude"][:]
        raw = ds["variable"][:, :, :]   # (lon, lat, z=35 years)

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
        annual = out
    end
    return lon, lat, annual
end

ev_lon, ev_lat, ev_annual = load_dong_lma(FILE_EV)
dc_lon, dc_lat, dc_annual = load_dong_lma(FILE_DC)

valid_ev = filter(isfinite, vec(ev_annual))
valid_dc = filter(isfinite, vec(dc_annual))
println("  Evergreen — n=$(length(valid_ev)), min=$(round(minimum(valid_ev), digits=4)), max=$(round(maximum(valid_ev), digits=4)), mean=$(round(mean(valid_ev), digits=4))")
println("  Deciduous — n=$(length(valid_dc)), min=$(round(minimum(valid_dc), digits=4)), max=$(round(maximum(valid_dc), digits=4)), mean=$(round(mean(valid_dc), digits=4))")

# ── Check overlap ──
has_ev = .!isnan.(ev_annual)
has_dc = .!isnan.(dc_annual)
both   = has_ev .& has_dc
only_ev = has_ev .& .!has_dc
only_dc = .!has_ev .& has_dc
neither = .!has_ev .& .!has_dc

n_both    = count(both)
n_only_ev = count(only_ev)
n_only_dc = count(only_dc)
n_neither = count(neither)
n_total   = length(ev_annual)

println("\n─── Spatial overlap ───")
println("  Both evergreen & deciduous: $n_both  ($(round(100*n_both/n_total, digits=2))%)")
println("  Evergreen only:             $n_only_ev  ($(round(100*n_only_ev/n_total, digits=2))%)")
println("  Deciduous only:             $n_only_dc  ($(round(100*n_only_dc/n_total, digits=2))%)")
println("  Neither (ocean/missing):    $n_neither  ($(round(100*n_neither/n_total, digits=2))%)")


# ═════════════════════════════════════════════════════════════════════════════
#  STEP 2 — Plots
# ═════════════════════════════════════════════════════════════════════════════

# ─────────────────────────────────────────────────────────────────────────────
#  Plot 1 — Side-by-side maps
# ─────────────────────────────────────────────────────────────────────────────
println("\n─── Plot 1: Side-by-side LMA maps ───")

# Use a shared colorrange for direct comparison
cmin = min(minimum(valid_ev), minimum(valid_dc))
cmax = max(maximum(valid_ev), maximum(valid_dc))

fig1 = Figure(size=(1600, 600))

ax1a = Axis(fig1[1, 1],
    title  = "Evergreen LMA (Dong et al.)\nClimatological mean (1982–2016)",
    xlabel = "Longitude",
    ylabel = "Latitude")
hm1a = heatmap!(ax1a, ev_lon, ev_lat, ev_annual;
    colormap = :YlGn, colorrange = (cmin, cmax), nan_color = :gray90)
Colorbar(fig1[1, 2], hm1a, label="LMA [g m⁻²]")

ax1b = Axis(fig1[1, 3],
    title  = "Deciduous LMA (Dong et al.)\nClimatological mean (1982–2016)",
    xlabel = "Longitude",
    ylabel = "Latitude")
hm1b = heatmap!(ax1b, dc_lon, dc_lat, dc_annual;
    colormap = :YlOrBr, colorrange = (cmin, cmax), nan_color = :gray90)
Colorbar(fig1[1, 4], hm1b, label="LMA [g m⁻²]")

save(joinpath(PLOT_DIR, "plot1_lma_side_by_side.png"), fig1, px_per_unit=3)
println("  Saved: plot1_lma_side_by_side.png")


# ─────────────────────────────────────────────────────────────────────────────
#  Plot 2 — Overlap classification map
#   0 = neither, 1 = evergreen only, 2 = deciduous only, 3 = both
# ─────────────────────────────────────────────────────────────────────────────
println("\n─── Plot 2: Overlap classification map ───")

overlap_map = zeros(Float64, size(ev_annual))
overlap_map[only_ev] .= 1.0
overlap_map[only_dc] .= 2.0
overlap_map[both]    .= 3.0
# Set ocean/neither to NaN so it renders as background
overlap_map[neither] .= NaN

fig2 = Figure(size=(1000, 500))
ax2 = Axis(fig2[1, 1],
    title  = "LMA Spatial Coverage — Evergreen vs Deciduous\n(Dong et al. 2023)",
    xlabel = "Longitude",
    ylabel = "Latitude")

hm2 = heatmap!(ax2, ev_lon, ev_lat, overlap_map;
    colormap = cgrad([:forestgreen, :orange, :purple], 3, categorical=true),
    colorrange = (0.5, 3.5),
    nan_color = :gray90)

# Custom legend instead of colorbar
Legend(fig2[1, 2],
    [MarkerElement(marker=:rect, color=:forestgreen, markersize=20),
     MarkerElement(marker=:rect, color=:orange, markersize=20),
     MarkerElement(marker=:rect, color=:purple, markersize=20)],
    ["Evergreen only ($n_only_ev px)",
     "Deciduous only ($n_only_dc px)",
     "Both ($n_both px)"],
    "Coverage")

save(joinpath(PLOT_DIR, "plot2_lma_overlap.png"), fig2, px_per_unit=3)
println("  Saved: plot2_lma_overlap.png")


println("\n═══ All LMA plots complete ═══")
