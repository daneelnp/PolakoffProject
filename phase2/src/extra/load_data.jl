"""
    load_data.jl

Read and inspect the four NetCDF datasets for the leaf-optics project.
Uses NCDatasets.jl to load each file, print metadata, and return key arrays.

Datasets (all from Dong et al. 2023 / Wang & Braghiere et al. 2025):
  • clm_refl_tran_1m_weighted.nc — leaf α & τ (PAR, NIR), monthly climatology
        Grid: 144×96 (2.5°), time dim: `ind` (1–12 = Jan–Dec)
  • TS_LMA_decidudous.nc — deciduous LMA (g/m²), annual 1982–2016
        Grid: 720×360 (0.5°), time dim: `z` (1–35 = years 1982–2016)
  • TS_LMA_evergreen.nc  — evergreen LMA (g/m²), annual 1982–2016
        Grid: 720×360 (0.5°), time dim: `z` (1–35 = years 1982–2016)
  • TS_Vcmax25.nc — Vcmax25 (µmol/m²/s), annual 1982–2016
        Grid: 720×360 (0.5°), time dim: `z` (1–35 = years 1982–2016)

Usage:
    julia --project=. src/load_data.jl
"""

using NCDatasets
using Statistics

# ── Paths ────────────────────────────────────────────────────────────────────
const DATA_DIR = joinpath(@__DIR__, "..", "data")

const FILE_REFL_TRAN = joinpath(DATA_DIR, "clm_refl_tran_1m_weighted.nc")
const FILE_LMA_DECID = joinpath(DATA_DIR, "TS_LMA_decidudous.nc")   # note spelling in filename
const FILE_LMA_EVER  = joinpath(DATA_DIR, "TS_LMA_evergreen.nc")
const FILE_VCMAX25   = joinpath(DATA_DIR, "TS_Vcmax25.nc")

# ── Helpers ──────────────────────────────────────────────────────────────────

"""
    print_file_summary(filepath)

Open a NetCDF file and print all dimensions, variables (with their dims,
shape, dtype, and attributes). Returns nothing; purely diagnostic.
"""
function print_file_summary(filepath::String)
    println("\n", "="^80)
    println("FILE: ", basename(filepath))
    println("="^80)

    NCDataset(filepath, "r") do ds
        # Global attributes
        gatts = keys(ds.attrib)
        if !isempty(gatts)
            println("\n── Global Attributes ──")
            for k in gatts
                println("  $k = ", ds.attrib[k])
            end
        end

        # Dimensions
        println("\n── Dimensions ──")
        for (name, dim) in ds.dim
            println("  $name : $(length(dim))")
        end

        # Variables
        println("\n── Variables ──")
        for varname in keys(ds)
            v = ds[varname]
            println("\n  Variable: $varname")
            println("    dimensions : $(dimnames(v))")
            println("    size       : $(size(v))")
            println("    eltype     : $(eltype(v))")

            # Variable attributes
            for aname in keys(v.attrib)
                println("    @$aname = ", v.attrib[aname])
            end
        end
    end
end

"""
    load_array_summary(filepath, varnames)

For each variable name in `varnames`, load the full array and print basic
statistics: shape, min, max, mean, and count of missing values.
Returns a Dict{String, Array} of the loaded data.
"""
function load_array_summary(filepath::String, varnames::Vector{String})
    data = Dict{String, Any}()

    NCDataset(filepath, "r") do ds
        for vn in varnames
            if !haskey(ds, vn)
                @warn "Variable '$vn' not found in $(basename(filepath)), skipping."
                continue
            end

            raw = ds[vn][:]          # load full array (Union{Missing, T})
            data[vn] = raw

            n_total   = length(raw)
            n_miss    = count(ismissing, raw)
            n_valid   = n_total - n_miss
            valid     = collect(skipmissing(raw))

            println("\n  ▸ $vn  (shape=$(size(raw)), eltype=$(eltype(raw)))")
            println("    valid values : $n_valid / $n_total")
            println("    missing      : $n_miss")
            if n_valid > 0
                println("    min          : $(minimum(valid))")
                println("    max          : $(maximum(valid))")
                println("    mean         : $(mean(valid))")
            end
        end
    end

    return data
end

"""
    load_coords(filepath)

Load dimension coordinate arrays (time, lat, lon) from a NetCDF file.
Also handles non-standard time dimensions:
  • `ind` 1–12  → monthly climatology (Jan–Dec)
  • `z`  1–35   → annual timeseries (1982–2016, from Dong et al. 2023)
Returns a NamedTuple (time, time_label, lat, lon) — any missing coordinate is `nothing`.
"""
function load_coords(filepath::String)
    time_vals  = nothing
    time_label = nothing
    lat_vals   = nothing
    lon_vals   = nothing

    NCDataset(filepath, "r") do ds
        # Standard time coordinates
        for tname in ("time", "Time", "t")
            if haskey(ds, tname)
                time_vals  = ds[tname][:]
                time_label = "time"
                break
            end
        end

        # Non-standard: `ind` = monthly index (1–12)
        if time_vals === nothing && haskey(ds, "ind")
            time_vals  = ds["ind"][:]
            time_label = "month index (ind)"
        end

        # Non-standard: `z` = year index (1–35 → 1982–2016)
        if time_vals === nothing && haskey(ds, "z")
            z = ds["z"][:]
            time_vals  = 1981 .+ z          # map 1→1982, …, 35→2016
            time_label = "year (z → 1982–2016)"
        end

        for latname in ("lat", "latitude", "Lat", "y")
            if haskey(ds, latname)
                lat_vals = ds[latname][:]
                break
            end
        end
        for lonname in ("lon", "longitude", "Lon", "x")
            if haskey(ds, lonname)
                lon_vals = ds[lonname][:]
                break
            end
        end
    end

    if time_vals !== nothing
        println("    time ($time_label) : $(first(time_vals))  →  $(last(time_vals))  ($(length(time_vals)) steps)")
    end
    if lat_vals !== nothing
        println("    lat  range   : $(minimum(lat_vals)) → $(maximum(lat_vals))  ($(length(lat_vals)) pts)")
    end
    if lon_vals !== nothing
        println("    lon  range   : $(minimum(lon_vals)) → $(maximum(lon_vals))  ($(length(lon_vals)) pts)")
    end

    return (; time=time_vals, time_label=time_label, lat=lat_vals, lon=lon_vals)
end

# ── Per-file loaders ─────────────────────────────────────────────────────────

"""
Load the CLM reflectance / transmittance dataset.
Returns (coords, data_dict).
"""
function load_refl_tran()
    println("\n", "#"^80)
    println("# LOADING: clm_refl_tran_1m_weighted.nc")
    println("#"^80)

    print_file_summary(FILE_REFL_TRAN)

    println("\n── Coordinates ──")
    coords = load_coords(FILE_REFL_TRAN)

    # Discover non-coordinate data variables (skip dim coord vars)
    data_vars = String[]
    NCDataset(FILE_REFL_TRAN, "r") do ds
        skip = Set(["time","Time","t","lat","latitude","lon","longitude","Lat","Lon","x","y",
                     "ind","z","crs"])
        for vn in keys(ds)
            vn ∉ skip && push!(data_vars, vn)
        end
    end

    println("\n── Data Variable Statistics ──")
    data = load_array_summary(FILE_REFL_TRAN, data_vars)

    return coords, data
end

"""
Load the LMA deciduous dataset.
Returns (coords, data_dict).
"""
function load_lma_deciduous()
    println("\n", "#"^80)
    println("# LOADING: TS_LMA_decidudous.nc")
    println("#"^80)

    print_file_summary(FILE_LMA_DECID)

    println("\n── Coordinates ──")
    coords = load_coords(FILE_LMA_DECID)

    data_vars = String[]
    NCDataset(FILE_LMA_DECID, "r") do ds
        skip = Set(["time","Time","t","lat","latitude","lon","longitude","Lat","Lon","x","y",
                     "ind","z","crs"])
        for vn in keys(ds)
            vn ∉ skip && push!(data_vars, vn)
        end
    end

    println("\n── Data Variable Statistics ──")
    data = load_array_summary(FILE_LMA_DECID, data_vars)

    return coords, data
end

"""
Load the LMA evergreen dataset.
Returns (coords, data_dict).
"""
function load_lma_evergreen()
    println("\n", "#"^80)
    println("# LOADING: TS_LMA_evergreen.nc")
    println("#"^80)

    print_file_summary(FILE_LMA_EVER)

    println("\n── Coordinates ──")
    coords = load_coords(FILE_LMA_EVER)

    data_vars = String[]
    NCDataset(FILE_LMA_EVER, "r") do ds
        skip = Set(["time","Time","t","lat","latitude","lon","longitude","Lat","Lon","x","y",
                     "ind","z","crs"])
        for vn in keys(ds)
            vn ∉ skip && push!(data_vars, vn)
        end
    end

    println("\n── Data Variable Statistics ──")
    data = load_array_summary(FILE_LMA_EVER, data_vars)

    return coords, data
end

"""
Load the Vcmax25 dataset.
Returns (coords, data_dict).
"""
function load_vcmax25()
    println("\n", "#"^80)
    println("# LOADING: TS_Vcmax25.nc")
    println("#"^80)

    print_file_summary(FILE_VCMAX25)

    println("\n── Coordinates ──")
    coords = load_coords(FILE_VCMAX25)

    data_vars = String[]
    NCDataset(FILE_VCMAX25, "r") do ds
        skip = Set(["time","Time","t","lat","latitude","lon","longitude","Lat","Lon","x","y",
                     "ind","z","crs"])
        for vn in keys(ds)
            vn ∉ skip && push!(data_vars, vn)
        end
    end

    println("\n── Data Variable Statistics ──")
    data = load_array_summary(FILE_VCMAX25, data_vars)

    return coords, data
end

# ── Main ─────────────────────────────────────────────────────────────────────

function main()
    println("Loading all datasets from: $DATA_DIR\n")

    for f in [FILE_REFL_TRAN, FILE_LMA_DECID, FILE_LMA_EVER, FILE_VCMAX25]
        if !isfile(f)
            error("File not found: $f")
        end
    end

    refl_coords,  refl_data  = load_refl_tran()
    decid_coords, decid_data = load_lma_deciduous()
    ever_coords,  ever_data  = load_lma_evergreen()
    vcmax_coords, vcmax_data = load_vcmax25()

    println("\n", "="^80)
    println("All four datasets loaded successfully.")
    println("="^80)

    return (
        refl  = (coords=refl_coords,  data=refl_data),
        decid = (coords=decid_coords, data=decid_data),
        ever  = (coords=ever_coords,  data=ever_data),
        vcmax = (coords=vcmax_coords, data=vcmax_data),
    )
end

# Run when executed as a script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
