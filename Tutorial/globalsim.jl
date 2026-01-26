import ClimaComms
ClimaComms.@import_required_backends
using ClimaUtilities
import Interpolations
import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import ClimaParams as CP
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!
using Dates
using CairoMakie, ClimaAnalysis, GeoMakie, Printf, StatsBase
import ClimaLand.LandSimVis as LandSimVis;

const FT = Float64;
context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
root_path = "land_longrun_$(device_suffix)"
diagnostics_outdir = joinpath(root_path, "global_diagnostics")
outdir =
    ClimaUtilities.OutputPathGenerator.generate_output_path(diagnostics_outdir);
toml_dict = LP.create_toml_dict(FT)

parameter_log_file = joinpath(root_path, "parameters.toml")
CP.log_parameter_information(toml_dict, parameter_log_file)

Δt = 450.0
start_date = DateTime(2008)
stop_date = DateTime(2009);

nelements = (20, 7)
domain = ClimaLand.Domains.global_domain(FT; context, nelements);

forcing = ClimaLand.prescribed_forcing_era5(
    start_date,
    stop_date,
    domain.space.surface,
    toml_dict,
    FT;
    use_lowres_forcing = true,
    max_wind_speed = 25.0,
    time_interpolation_method = LinearInterpolation(PeriodicCalendar()),
    regridder_type = :InterpolationsRegridder,
    context,
);

LAI = ClimaLand.Canopy.prescribed_lai_modis(
    domain.space.surface,
    start_date,
    stop_date,
);

model = ClimaLand.LandModel{FT}(forcing, LAI, toml_dict, domain, Δt);

simulation = ClimaLand.Simulations.LandSimulation(
    start_date,
    stop_date,
    Δt,
    model;
    outdir,
    user_callbacks = (),
);

ClimaLand.Simulations.solve!(simulation)
LandSimVis.make_annual_timeseries(simulation; savedir = root_path)
LandSimVis.make_heatmaps(simulation;date = stop_date, savedir = root_path)
LandSimVis.make_leaderboard_plots(simulation, savedir = root_path)