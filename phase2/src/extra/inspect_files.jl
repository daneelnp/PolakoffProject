using NCDatasets

data_dir = joinpath(@__DIR__, "..", "data")

files = [
    ("vcmax25_1M_average.nc", "CliMA Land simulation output"),
    ("TS_Vcmax25.nc",         "Dong et al. 2023 optimality-based predictions"),
]

for (fname, label) in files
    fpath = joinpath(data_dir, fname)
    println("=" ^ 80)
    println("FILE: $fname")
    println("LABEL: $label")
    println("=" ^ 80)

    ds = NCDataset(fpath, "r")

    # --- All variable names ---
    println("\n📦 Variables:")
    for vname in keys(ds)
        println("  • $vname")
    end

    # --- Dimensions and sizes ---
    println("\n📐 Dimensions:")
    for (dname, dobj) in ds.dim
        println("  • $dname : $(dobj)")
    end

    # --- Per-variable details ---
    println("\n📋 Variable details:")
    for vname in keys(ds)
        v = ds[vname]
        dims = dimnames(v)
        sz = size(v)
        units_str = haskey(v.attrib, "units") ? v.attrib["units"] : "(no units)"
        long_name = haskey(v.attrib, "long_name") ? v.attrib["long_name"] : "(no long_name)"
        println("  ─ $vname")
        println("      dims:      $dims")
        println("      size:      $sz")
        println("      units:     $units_str")
        println("      long_name: $long_name")
    end

    # --- Time axis info ---
    time_vars = filter(v -> any(d -> lowercase(d) in ["time", "month", "t"], dimnames(ds[v])) || lowercase(v) in ["time", "month", "t"], keys(ds))
    # Also check if "time" or "month" is a variable directly
    for tname in ["time", "Time", "month", "Month", "t"]
        if haskey(ds, tname)
            tvar = ds[tname]
            tdata = tvar[:]
            println("\n🕐 Time axis variable: '$tname'")
            println("      length:  $(length(tdata))")
            println("      first:   $(first(tdata))")
            println("      last:    $(last(tdata))")
            println("      type:    $(eltype(tdata))")
        end
    end

    # --- Vcmax25 stats ---
    # Try common variable names
    vcmax_candidates = ["vcmax25", "Vcmax25", "VCMAX25", "vcmax_25", "Vcmax_25", "vcmax25_annual_mean"]
    found = false
    for cand in keys(ds)
        if lowercase(cand) in ["vcmax25", "vcmax_25", "vcmax25_annual_mean", "vcmax"]
            data = ds[cand][:]
            valid = filter(!ismissing, data)
            valid_f = Float64.(valid)
            println("\n📊 Stats for '$cand':")
            println("      shape:     $(size(data))")
            println("      # missing: $(count(ismissing, data))")
            println("      # valid:   $(length(valid_f))")
            if length(valid_f) > 0
                println("      min:       $(minimum(valid_f))")
                println("      max:       $(maximum(valid_f))")
                println("      mean:      $(round(mean(valid_f), digits=4))")
            end
            found = true
        end
    end

    if !found
        # If not found, just print stats for every non-coordinate variable
        println("\n⚠ Could not auto-detect Vcmax25 variable. Printing stats for all data variables:")
        for vname in keys(ds)
            v = ds[vname]
            if length(dimnames(v)) >= 2  # likely a data variable
                data = v[:]
                valid = filter(!ismissing, data)
                if length(valid) > 0
                    valid_f = Float64.(valid)
                    println("  📊 Stats for '$vname':")
                    println("      shape:     $(size(data))")
                    println("      # missing: $(count(ismissing, data))")
                    println("      # valid:   $(length(valid_f))")
                    println("      min:       $(minimum(valid_f))")
                    println("      max:       $(maximum(valid_f))")
                    println("      mean:      $(round(mean(valid_f), digits=4))")
                end
            end
        end
    end

    close(ds)
    println("\n")
end
