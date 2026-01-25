import ClimaParams as CP
using ClimaLand
using ClimaLand.Domains
using ClimaLand.Soil
import ClimaLand.Simulations: LandSimulation, solve!
import ClimaLand.Parameters as LP
import ClimaLand.LandSimVis as LandSimVis;
root_path = ""
using Dates

FT = Float32
toml_dict = LP.create_toml_dict(FT);

zmax = FT(0)
zmin = FT(-1.0)
longlat = FT.((-118.1, 34.1))
domain = Domains.Column(; zlim = (zmin, zmax), nelements = 10, longlat);
surface_space = domain.space.surface;

start_date = DateTime(2008);
stop_date = start_date + Second(60 * 60 * 72);
dt = 1000.0; # time step in seconds

atmos, radiation = ClimaLand.prescribed_forcing_era5(
    start_date,
    stop_date,
    surface_space,
    toml_dict,
    FT;
    use_lowres_forcing = true,
);

model = Soil.EnergyHydrology{FT}(
    domain,
    (; atmos, radiation),
    toml_dict,
);

function set_ic!(Y, p, t0, model)
    Y.soil.ϑ_l .= FT(0.24);
    Y.soil.θ_i .= FT(0.0);
    T = FT(290.15);
    ρc_s = Soil.volumetric_heat_capacity.(
        Y.soil.ϑ_l,
        Y.soil.θ_i,
        model.parameters.ρc_ds,
        model.parameters.earth_param_set,
    );
    Y.soil.ρe_int .=
        Soil.volumetric_internal_energy.(
            Y.soil.θ_i,
            ρc_s,
            T,
            model.parameters.earth_param_set,
        );
end

simulation = LandSimulation(start_date, stop_date, dt, model; set_ic!, user_callbacks = ());

solve!(simulation);

# LandSimVis.make_annual_timeseries(simulation)