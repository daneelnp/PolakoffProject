#=
forcing_comparison.jl
=====================
Side-by-side validation of CRU-JRA vs ERA5 forcing data + model outputs.

This script does three things:
  1. Evaluates the RAW forcing fields (T, P, q, u, precip, SW, LW) from both
     datasets at every hour over the simulation window and compares them directly.
  2. Runs a standalone CanopyModel with each forcing set and collects diagnostics
     that include both model outputs AND driver fields as the model sees them.
  3. Produces overlay plots and prints summary statistics (mean, max, RMSE of
     the differences).

Usage (on Clima cluster with GPU):
    srun --gpus=1 --mpi=none -t 02:00:00 --pty bash -l
    cd /home/dpolakof/dp_CliMAProject/ClimaLand.jl
    julia --project=. forcing_comparison.jl

Usage (CPU, e.g. local laptop — will just be slower):
    julia --project=. forcing_comparison.jl
=#

# ── 0. Imports ──────────────────────────────────────────────────────────────────

import ClimaParams as CP
import ClimaLand
using ClimaLand: PrescribedGroundConditions, evaluate!
using ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
using ClimaUtilities.TimeManager: ITime, date
using ClimaLand.Domains
using ClimaLand.Canopy
import ClimaLand.Simulations: LandSimulation, solve!
import ClimaLand.Parameters as LP
using Dates
using Statistics
import ClimaDiagnostics
using CairoMakie
using ClimaCore

# ── 1. Common setup ─────────────────────────────────────────────────────────────

const FT        = Float32
const toml_dict  = LP.create_toml_dict(FT)

# Point location: Amazon tropical forest (lon, lat)
const longlat      = FT.((-60.0, -3.0))
const domain       = Domains.Point(; z_sfc = FT(0.0), longlat)
const surface_space = domain.space.surface

const start_date = DateTime(2008)
const stop_date  = start_date + Day(14)   # 14-day comparison window
const dt         = 900.0                  # 15-minute timestep

# Output directory for all plots
const outdir = joinpath(@__DIR__, "forcing_comparison_output_amazon")
mkpath(outdir)

# ── 2. Build both forcing sets ──────────────────────────────────────────────────

println("┌─ Building ERA5 (high-res) forcing …")
atmos_era5, rad_era5 = ClimaLand.prescribed_forcing_era5(
    start_date, stop_date, surface_space, toml_dict, FT;
    use_lowres_forcing = false,
)
println("└─ done.")

println("┌─ Building CRU-JRA forcing …")
atmos_cru, rad_cru = ClimaLand.prescribed_forcing_crujra(
    start_date, stop_date, surface_space, toml_dict, FT,
)
println("└─ done.")

# ── 3. Raw forcing comparison (evaluate! loop) ──────────────────────────────────
#
# We sample every hour for 7 days and compare the raw interpolated values
# from each dataset before they ever enter the model.

println("\n══════ Part 1: Raw Forcing Field Comparison ══════")

# Time vector in seconds from start_date: every hour for 14 days
const n_hours    = 14 * 24
const eval_times = Float64.(0:3600:(n_hours * 3600))   # length n_hours+1
const eval_dates = [start_date + Second(Int(t)) for t in eval_times]

# Helper: evaluate a TimeVaryingInput at every sample time, return Float64 vector
function eval_timeseries(tvi, times, buf)
    out = Vector{Float64}(undef, length(times))
    for (i, t) in enumerate(times)
        evaluate!(buf, tvi, t)
        out[i] = Float64(Array(parent(buf))[1])
    end
    return out
end

buf = ClimaCore.Fields.zeros(FT, surface_space)

# Define which forcing fields to compare
# (label, ERA5 field, CRU-JRA field, units)
forcing_fields = [
    ("Temperature",         atmos_era5.T,            atmos_cru.T,            "K"),
    ("Pressure",            atmos_era5.P,            atmos_cru.P,            "Pa"),
    ("Specific Humidity",   atmos_era5.q,            atmos_cru.q,            "kg/kg"),
    ("Wind Speed",          atmos_era5.u,            atmos_cru.u,            "m/s"),
    ("Liquid Precip",       atmos_era5.liquid_precip, atmos_cru.liquid_precip, "m/s"),
    ("Snow Precip",         atmos_era5.snow_precip,  atmos_cru.snow_precip,  "m/s"),
    ("SW Down",             rad_era5.SW_d,           rad_cru.SW_d,           "W/m²"),
    ("LW Down",             rad_era5.LW_d,           rad_cru.LW_d,           "W/m²"),
]

# Collect all timeseries
era5_series = Dict{String, Vector{Float64}}()
cru_series  = Dict{String, Vector{Float64}}()

for (label, era5_field, cru_field, _) in forcing_fields
    era5_series[label] = eval_timeseries(era5_field, eval_times, buf)
    cru_series[label]  = eval_timeseries(cru_field,  eval_times, buf)
end

# ── 3a. Print summary statistics ────────────────────────────────────────────────

println("\n  Variable                | ERA5 mean  | CRU mean   | Mean Diff  | Max |Diff|  | RMSE")
println("  " * "─"^95)
for (label, _, _, units) in forcing_fields
    e = era5_series[label]
    c = cru_series[label]
    d = e .- c
    mn_e  = round(mean(e); sigdigits=5)
    mn_c  = round(mean(c); sigdigits=5)
    mn_d  = round(mean(d); sigdigits=4)
    mx_d  = round(maximum(abs.(d)); sigdigits=4)
    rmse  = round(sqrt(mean(d .^ 2)); sigdigits=4)
    println("  $(rpad(label, 24))| $(rpad(mn_e, 11))| $(rpad(mn_c, 11))| $(rpad(mn_d, 11))| $(rpad(mx_d, 11))| $rmse  [$units]")
end

# ── 3b. NaN / missing-value check ───────────────────────────────────────────────

println("\n  NaN check:")
any_nan = false
for (label, _, _, _) in forcing_fields
    n_era = count(isnan, era5_series[label])
    n_cru = count(isnan, cru_series[label])
    if n_era > 0 || n_cru > 0
        println("    ⚠  $label: ERA5 has $n_era NaNs, CRU-JRA has $n_cru NaNs")
        global any_nan = true
    end
end
any_nan || println("    ✓  No NaNs detected in any forcing field.")

# ── 3c. Plot raw forcing overlays ────────────────────────────────────────────────

println("\n  Saving raw forcing overlay plots to $outdir …")

for (label, _, _, units) in forcing_fields
    fig = Figure(size = (900, 350))
    ax  = Axis(fig[1, 1];
        xlabel = "Date (UTC)",
        ylabel = "$label [$units]",
        title  = "Raw Forcing: $label",
    )
    lines!(ax, eval_dates, era5_series[label]; label = "ERA5", color = :blue)
    lines!(ax, eval_dates, cru_series[label];  label = "CRU-JRA", color = :red)
    axislegend(ax; position = :lt)

    fname = replace(lowercase(label), " " => "_")
    save(joinpath(outdir, "raw_$(fname).png"), fig)
end

# Difference plots
for (label, _, _, units) in forcing_fields
    fig = Figure(size = (900, 350))
    ax  = Axis(fig[1, 1];
        xlabel = "Date (UTC)",
        ylabel = "ERA5 − CRU [$units]",
        title  = "Forcing Difference: $label",
    )
    diff = era5_series[label] .- cru_series[label]
    lines!(ax, eval_dates, diff; color = :black)
    hlines!(ax, [0.0]; color = :gray, linestyle = :dash)

    fname = replace(lowercase(label), " " => "_")
    save(joinpath(outdir, "diff_$(fname).png"), fig)
end

println("  ✓  Raw forcing plots saved.\n")


# ══════════════════════════════════════════════════════════════════════════════════
# ── 4. Run canopy simulations with both forcings ─────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════════

println("══════ Part 2: Canopy Model Simulation Comparison ══════\n")

# Diagnostic variables: model outputs + driver fields the model actually uses
output_vars_model   = ["gpp", "ct", "lai", "trans", "er", "sif", "lhf", "shf"]
output_vars_drivers = ["swd", "lwd", "rain", "snow", "ws", "airp"]
output_vars = vcat(output_vars_model, output_vars_drivers)

ground = PrescribedGroundConditions{FT}()
LAI    = TimeVaryingInput((t) -> FT(1.0))

"""
    run_canopy_sim(atmos, radiation; label)

Build and run a standalone canopy simulation, returning the DictWriter
containing all diagnostics.
"""
function run_canopy_sim(atmos, radiation; label::String)
    println("  ┌─ Building CanopyModel ($label) …")
    model = Canopy.CanopyModel{FT}(
        domain, (; atmos, radiation, ground), LAI, toml_dict,
    )

    function set_ic!(Y, p, t0, model)
        ψ_leaf_0 = FT(-2e5 / 9800)
        (; retention_model, ν, S_s) = model.hydraulics.parameters
        S_l_ini = Canopy.PlantHydraulics.inverse_water_retention_curve(
            retention_model, ψ_leaf_0, ν, S_s,
        )
        Y.canopy.hydraulics.ϑ_l.:1 .=
            Canopy.PlantHydraulics.augmented_liquid_fraction.(ν, S_l_ini)
        evaluate!(Y.canopy.energy.T, atmos.T, t0)
    end

    writer = ClimaDiagnostics.Writers.DictWriter()
    diagnostics = ClimaLand.Diagnostics.default_diagnostics(
        model, start_date;
        output_vars  = output_vars,
        output_writer = writer,
        reduction_period = :hourly,
    )

    simulation = LandSimulation(
        start_date, stop_date, dt, model;
        set_ic!,
        updateat      = Second(dt),
        user_callbacks = (),
        diagnostics,
    )

    println("  │  Solving …")
    solve!(simulation)
    println("  └─ $label simulation complete.")
    return writer, simulation
end

writer_era5, sim_era5 = run_canopy_sim(atmos_era5, rad_era5; label = "ERA5")
writer_cru,  sim_cru  = run_canopy_sim(atmos_cru,  rad_cru;  label = "CRU-JRA")


# ── 5. Extract & compare simulation diagnostics ─────────────────────────────────

println("\n  Extracting diagnostic timeseries …")

"""
    extract_diag(writer, short_name)

Pull (times, values) vectors from a DictWriter for a given diagnostic short_name.
Looks up the actual output key (e.g. "gpp_1h_average") automatically.
"""
function extract_diag(writer, short_name)
    # DictWriter keys look like "gpp_1h_average" — find the one starting with our short_name
    matching_keys = filter(k -> startswith(k, short_name * "_"), collect(keys(writer.dict)))
    if isempty(matching_keys)
        # Try exact match
        matching_keys = filter(k -> k == short_name, collect(keys(writer.dict)))
    end
    isempty(matching_keys) && error("No diagnostic key found for '$short_name'. Available: $(keys(writer.dict))")
    key = first(matching_keys)
    times_raw, vals = ClimaLand.Diagnostics.diagnostic_as_vectors(writer, key)
    # Convert ITime objects to DateTime
    dates = [date(t) for t in times_raw]
    return dates, vals
end


# ── 5a. Print model-output comparison statistics ────────────────────────────────

println("\n  Diagnostic          | ERA5 mean  | CRU mean   | Mean Diff  | Max |Diff|  | RMSE")
println("  " * "─"^90)
for sn in output_vars
    try
        t_e, v_e = extract_diag(writer_era5, sn)
        t_c, v_c = extract_diag(writer_cru,  sn)
        # Align lengths (should be identical, but be safe)
        n = min(length(v_e), length(v_c))
        d = v_e[1:n] .- v_c[1:n]
        mn_e = round(mean(v_e[1:n]); sigdigits=5)
        mn_c = round(mean(v_c[1:n]); sigdigits=5)
        mn_d = round(mean(d); sigdigits=4)
        mx_d = round(maximum(abs.(d)); sigdigits=4)
        rmse = round(sqrt(mean(d .^ 2)); sigdigits=4)
        println("  $(rpad(sn, 20))| $(rpad(mn_e, 11))| $(rpad(mn_c, 11))| $(rpad(mn_d, 11))| $(rpad(mx_d, 11))| $rmse")
    catch e
        println("  $(rpad(sn, 20))| ⚠ skipped: $(sprint(showerror, e))")
    end
end


# ── 5b. Overlay diagnostic timeseries plots ──────────────────────────────────────

println("\n  Saving diagnostic overlay plots to $outdir …")

for sn in output_vars
    try
        dates_e, v_e = extract_diag(writer_era5, sn)
        dates_c, v_c = extract_diag(writer_cru,  sn)

        # ── Overlay plot ──
        fig = Figure(size = (900, 350))
        ax  = Axis(fig[1, 1];
            xlabel = "Date (UTC)",
            ylabel = sn,
            title  = "Diagnostic: $sn",
        )
        lines!(ax, dates_e, v_e; label = "ERA5",    color = :blue)
        lines!(ax, dates_c, v_c; label = "CRU-JRA", color = :red)
        axislegend(ax; position = :lt)
        save(joinpath(outdir, "diag_$(sn).png"), fig)

        # ── Difference plot ──
        n = min(length(v_e), length(v_c))
        fig2 = Figure(size = (900, 300))
        ax2  = Axis(fig2[1, 1];
            xlabel = "Date (UTC)",
            ylabel = "ERA5 − CRU",
            title  = "Diagnostic Difference: $sn",
        )
        lines!(ax2, dates_e[1:n], v_e[1:n] .- v_c[1:n]; color = :black)
        hlines!(ax2, [0.0]; color = :gray, linestyle = :dash)
        save(joinpath(outdir, "diag_diff_$(sn).png"), fig2)
    catch e
        println("    ⚠  Skipping $sn: $(sprint(showerror, e))")
    end
end

println("  ✓  Diagnostic plots saved.\n")


# ── 6. Final summary ────────────────────────────────────────────────────────────

println("══════ Summary ══════")
println("  Output directory : $outdir")
println("  Raw forcing plots: raw_*.png, diff_*.png")
println("  Diagnostic plots : diag_*.png, diag_diff_*.png")
println("  Location         : lon=$(longlat[1])°, lat=$(longlat[2])°")
println("  Period           : $start_date → $stop_date")
println("")
println("  Review the plots and statistics above.")
println("  Key things to check before a global run:")
println("    1. No NaNs in CRU-JRA fields")
println("    2. Forcing magnitudes are physically reasonable")
println("    3. Model outputs (GPP, ET, etc.) are in the same ballpark")
println("    4. Differences are explainable by dataset resolution / source")
println("══════════════════════════════════════════════════════════")
