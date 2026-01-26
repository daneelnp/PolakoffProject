import ClimaParams as CP
using ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
using ClimaLand
using ClimaLand.Domains
using ClimaLand.Canopy
import ClimaLand.Simulations: LandSimulation, solve!
import ClimaLand.Parameters as LP
using Dates
import ClimaDiagnostics
using CairoMakie, ClimaAnalysis, GeoMakie, Printf, StatsBase
import ClimaLand.LandSimVis as LandSimVis

FT = Float32
toml_dict = LP.create_toml_dict(FT);

longlat = FT.((-110.6, 44.6))
domain = Domains.Point(; z_sfc = FT(0.0), longlat);
surface_space = domain.space.surface;

start_date = DateTime(2008);
stop_date = start_date + Second(60 * 60 * 72);
dt = 900.0;

atmos, radiation = ClimaLand.prescribed_forcing_era5(
    start_date,
    stop_date,
    surface_space,
    toml_dict,
    FT;
    use_lowres_forcing = true,
);
ground = PrescribedGroundConditions{FT}();
LAI = TimeVaryingInput((t) -> FT(1.0));

model = Canopy.CanopyModel{FT}(domain, (; atmos, radiation, ground), LAI, toml_dict);

function set_ic!(Y, p, t0, model)
    ψ_leaf_0 = FT(-2e5 / 9800)
    (; retention_model, ν, S_s) = model.hydraulics.parameters
    S_l_ini = Canopy.PlantHydraulics.inverse_water_retention_curve(
        retention_model,
        ψ_leaf_0,
        ν,
        S_s,
    )
    Y.canopy.hydraulics.ϑ_l.:1 .=
        Canopy.PlantHydraulics.augmented_liquid_fraction.(ν, S_l_ini)
    evaluate!(Y.canopy.energy.T, atmos.T, t0)
end

set_ic!

diag_writer = ClimaDiagnostics.Writers.DictWriter();
diagnostics = ClimaLand.Diagnostics.default_diagnostics(
    model,
    start_date;
    output_vars = ["ct", "trans"],
    output_writer = diag_writer,
    reduction_period = :hourly,
);

diag_writer = ClimaDiagnostics.Writers.DictWriter();
diagnostics = ClimaLand.Diagnostics.default_diagnostics(
    model,
    start_date;
    output_vars = ["ct", "trans"],
    output_writer = diag_writer,
    reduction_period = :hourly,
);

simulation = LandSimulation(
    start_date,
    stop_date,
    dt,
    model;
    set_ic!,
    updateat = Second(dt),
    user_callbacks = (),
    diagnostics,
);

solve!(simulation);

parameter_log_file = "default_canopy_parameters.toml"
CP.log_parameter_information(toml_dict, parameter_log_file)

LandSimVis.make_diurnal_timeseries(
    simulation;
    short_names = ["lai"],
    plot_stem_name = "default_canopy",
);