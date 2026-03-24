using NCDatasets
using CairoMakie
using Statistics
using Printf

# ─── Paths ────────────────────────────────────────────────────────────────────
data_path  = joinpath(@__DIR__, "..", "data", "clm_refl_tran_1m_weighted.nc")
output_dir = joinpath(@__DIR__, "..", "output_plots", "leaf_optical")
mkpath(output_dir)

# Fill value used in this file (not NaN — must be masked manually)
const FILL = -1.0

# ─── STEP 0: Inspect the file ─────────────────────────────────────────────────
println("=" ^ 60)
println("STEP 0 — File Inspection")
println("=" ^ 60)

ds = NCDataset(data_path, "r")

println("\nDimensions:")
for (name, dim) in ds.dim
    println("  $name = $dim")
end

println("\nVariables:")
for (name, var) in ds
    println("  $name  dims=$(dimnames(var))  size=$(size(var))")
end

optical_vars = ["par_refl", "par_tran", "nir_refl", "nir_tran"]

println("\nSummer (July, ind=7) statistics after masking fill value ($FILL):")
for vname in optical_vars
    raw  = Float64.(ds[vname][:, :, 7])   # lon × lat slice for July
    raw[raw .≈ FILL] .= NaN
    valid = filter(!isnan, vec(raw))
    @printf("  %-10s  min=%.4f  max=%.4f  mean=%.4f  N_valid=%d\n",
            vname, minimum(valid), maximum(valid), mean(valid), length(valid))
end

close(ds)
println("=" ^ 60, "\n")

# ─── Helpers ──────────────────────────────────────────────────────────────────

# Load one month's lon×lat slice and replace fill with NaN.
# month_ind is 1-based (1=Jan … 12=Dec).
function load_slice(ds, vname, month_ind)
    data = Float64.(ds[vname][:, :, month_ind])
    data[data .≈ FILL] .= NaN
    return data
end

# Read lon/lat from file; roll longitude from [0,360] to [-180,180]
# so the map is centred on the prime meridian.
function get_coords(ds)
    lon = Float64.(ds["lon"][:])
    lat = Float64.(ds["lat"][:])
    # Shift any longitudes > 180 by -360
    lon[lon .> 180] .-= 360
    # Re-sort so the array is monotonically increasing west→east
    order = sortperm(lon)
    return lon[order], lat, order
end

# Re-order the data columns to match the sorted longitude
function reorder_lon(data, order)
    return data[order, :]
end

# ─── Fixed color limits (same across all seasons for comparability) ────────────
# Update these if your data range differs.
clims = Dict(
    "par_refl" => (0.0, 0.5),
    "par_tran" => (0.0, 0.5),
    "nir_refl" => (0.0, 0.5),
    "nir_tran" => (0.0, 0.5),
)

var_labels = Dict(
    "par_refl" => "PAR Reflectance",
    "par_tran" => "PAR Transmittance",
    "nir_refl" => "NIR Reflectance",
    "nir_tran" => "NIR Transmittance",
)

# ─── Layout constants ─────────────────────────────────────────────────────────
FIG_W   = 1600   # figure width  (px)
FIG_H   =  480   # figure height (px)
FONTSIZE = 13

# ─── STEP 1: One figure per season, 4 panels per figure ───────────────────────
seasons = [
    ("Spring", "April",   4,  "seasonal_spring_april.png"),
    ("Summer", "July",    7,  "seasonal_summer_july.png"),
    ("Fall",   "October", 10, "seasonal_fall_october.png"),
    ("Winter", "January", 1,  "seasonal_winter_january.png"),
]

ds  = NCDataset(data_path, "r")
lon, lat, lon_order = get_coords(ds)

for (season_name, month_name, month_ind, fname) in seasons

    fig = Figure(size = (FIG_W, FIG_H),
                 figure_padding = (10, 20, 10, 10))

    # ── Suptitle ──────────────────────────────────────────────────────────────
    Label(fig[0, 1:4],
        "Leaf Optical Properties — $season_name ($month_name)";
        fontsize  = FONTSIZE + 5,
        font      = :bold,
        tellwidth = false,
        padding   = (0, 0, 8, 0),
    )

    for (col, vname) in enumerate(optical_vars)
        data = reorder_lon(load_slice(ds, vname, month_ind), lon_order)

        ax = Axis(fig[1, col];
            title      = var_labels[vname],
            titlesize  = FONTSIZE,
            xlabel     = "Longitude (°)",
            ylabel     = col == 1 ? "Latitude (°)" : "",
            xlabelsize = FONTSIZE - 1,
            ylabelsize = FONTSIZE - 1,
            xticklabelsize = FONTSIZE - 2,
            yticklabelsize = FONTSIZE - 2,
            xticks     = -180:60:180,
            yticks     = -90:30:90,
            # 2:1 geographic aspect (360° wide × 180° tall)
            aspect     = 2.0,
        )

        hm = heatmap!(ax, lon, lat, data;
            colormap   = :viridis,
            colorrange = clims[vname],
            nan_color  = :lightgray,  # ocean / missing land = gray
        )

        Colorbar(fig[2, col], hm;
            vertical       = false,
            width          = Relative(0.88),
            ticklabelsize  = FONTSIZE - 2,
            label          = var_labels[vname],
            labelsize      = FONTSIZE - 1,
        )
    end

    rowsize!(fig.layout, 2, Fixed(55))
    colgap!(fig.layout, 12)

    save(joinpath(output_dir, fname), fig; px_per_unit = 2)
    println("Saved: $fname")
end

# ─── STEP 2: Summer minus Winter difference, 4 panels ─────────────────────────
# Symmetric diverging limits; red = summer > winter, blue = winter > summer
diff_clims = (-0.3, 0.3)

fig_diff = Figure(size = (FIG_W, FIG_H),
                  figure_padding = (10, 20, 10, 10))

Label(fig_diff[0, 1:4],
    "Leaf Optical Properties — Summer (July) minus Winter (January)";
    fontsize  = FONTSIZE + 5,
    font      = :bold,
    tellwidth = false,
    padding   = (0, 0, 8, 0),
)

for (col, vname) in enumerate(optical_vars)
    july    = reorder_lon(load_slice(ds, vname, 7), lon_order)   # ind=7  July
    january = reorder_lon(load_slice(ds, vname, 1), lon_order)   # ind=1  January
    diff    = july .- january

    ax = Axis(fig_diff[1, col];
        title      = var_labels[vname],
        titlesize  = FONTSIZE,
        xlabel     = "Longitude (°)",
        ylabel     = col == 1 ? "Latitude (°)" : "",
        xlabelsize = FONTSIZE - 1,
        ylabelsize = FONTSIZE - 1,
        xticklabelsize = FONTSIZE - 2,
        yticklabelsize = FONTSIZE - 2,
        xticks     = -180:60:180,
        yticks     = -90:30:90,
        aspect     = 2.0,
    )

    hm = heatmap!(ax, lon, lat, diff;
        colormap   = :RdBu,           # diverging: red positive, blue negative
        colorrange = diff_clims,
        nan_color  = :lightgray,
    )

    Colorbar(fig_diff[2, col], hm;
        vertical       = false,
        width          = Relative(0.88),
        ticklabelsize  = FONTSIZE - 2,
        label          = "Δ $(var_labels[vname])  (Jul − Jan)",
        labelsize      = FONTSIZE - 1,
    )
end

rowsize!(fig_diff.layout, 2, Fixed(55))
colgap!(fig_diff.layout, 12)

save(joinpath(output_dir, "seasonal_summer_minus_winter_difference.png"),
     fig_diff; px_per_unit = 2)
println("Saved: seasonal_summer_minus_winter_difference.png")

close(ds)
println("\nDone.")
