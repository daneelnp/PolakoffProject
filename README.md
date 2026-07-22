# Predicting Leaf Optical Properties from Plant Traits

Random Forest models that predict leaf optical properties (PAR/NIR reflectance
and transmittance) from plant traits, integrated into [ClimaLand](https://github.com/CliMA/ClimaLand.jl).

## Parts

- **Phase 1** — Set up and validate ClimaLand canopy simulations and CRUJRA forcing as a baseline.
- **Phase 2** — Compute global plant-trait predictors and train the Random Forest optical-property models.
- **Phase 3** — Integrate the trained models into ClimaLand's TwoStream radiation scheme and run column and global experiments.

## Paper

`Predicting_Leaf_Optical_Properties_from_Plant_Traits_with_Random_Forest.pdf` writes up the methods and results of this work.
