"""
    plot_results.jl

Generate diagnostic plots for the leaf-optics log-linear fit (Option B).
Includes the fit pipeline from fit_leaf_optics.jl, then produces:

  1. Predictor scatter:  log(Vcmax25) vs log(LMA), colored by α_PAR
  2. Observed vs Fitted:  4-panel (α_PAR, τ_PAR, α_NIR, τ_NIR)
  3. Residual histograms: 4-panel
  4. Residual maps:       4-panel global maps of residuals on the 2.5° grid
  5. Predictor marginals: partial-dependence-style scatter for each target vs each predictor

Outputs saved to phase2/plots/

Usage:
    julia --project=. src/plot_results.jl
"""

using CairoMakie
using Statistics

# Include the fit pipeline
include(joinpath(@__DIR__, "fit_leaf_optics.jl"))

# ── Output directory ─────────────────────────────────────────────────────────
const PLOT_DIR = joinpath(@__DIR__, "..", "plots")
mkpath(PLOT_DIR)

# ═════════════════════════════════════════════════════════════════════════════
#  Plot 1 — Predictor space: log(Vcmax) vs log(LMA)
# ═════════════════════════════════════════════════════════════════════════════

function plot_predictor_space(tbl)
    fig = Figure(size=(700, 550))
    ax = Axis(fig[1,1],
        xlabel = "log(Vcmax25)",
        ylabel = "log(LMA)",
        title  = "Predictor Space — colored by α_PAR")

    sc = scatter!(ax, tbl.log_vcmax, tbl.log_lma;
        color = tbl.α_par,
        colormap = :viridis,
        markersize = 3,
        alpha = 0.6)

    Colorbar(fig[1,2], sc, label="α_PAR")

    # Add correlation annotation
    r = cor(tbl.log_vcmax, tbl.log_lma)
    text!(ax, 0.05, 0.95; text="r = $(round(r; digits=3))",
        space=:relative, fontsize=14, align=(:left, :top))

    save(joinpath(PLOT_DIR, "01_predictor_space.png"), fig, px_per_unit=3)
    println("  Saved: 01_predictor_space.png")
    return fig
end

# ═════════════════════════════════════════════════════════════════════════════
#  Plot 2 — Observed vs Fitted (4-panel)
# ═════════════════════════════════════════════════════════════════════════════

function plot_obs_vs_fitted(tbl, results)
    fig = Figure(size=(1000, 900))

    targets = [
        ("α_PAR", tbl.α_par),
        ("τ_PAR", tbl.τ_par),
        ("α_NIR", tbl.α_nir),
        ("τ_NIR", tbl.τ_nir),
    ]

    for (i, (name, obs)) in enumerate(targets)
        row = (i - 1) ÷ 2 + 1
        col = (i - 1) % 2 + 1
        res = results[name]

        ax = Axis(fig[row, col],
            xlabel = "Fitted $name",
            ylabel = "Observed $name",
            title  = "$name  (R² = $(round(res.R²; digits=3)), RMSE = $(round(res.rmse; digits=4)))",
            aspect = DataAspect())

        scatter!(ax, res.fitted, obs;
            markersize = 2, alpha = 0.4, color = :steelblue)

        # 1:1 line
        lo = min(minimum(obs), minimum(res.fitted))
        hi = max(maximum(obs), maximum(res.fitted))
        lines!(ax, [lo, hi], [lo, hi]; color=:red, linewidth=1.5, linestyle=:dash)
    end

    save(joinpath(PLOT_DIR, "02_obs_vs_fitted.png"), fig, px_per_unit=3)
    println("  Saved: 02_obs_vs_fitted.png")
    return fig
end

# ═════════════════════════════════════════════════════════════════════════════
#  Plot 3 — Residual histograms (4-panel)
# ═════════════════════════════════════════════════════════════════════════════

function plot_residual_histograms(results)
    fig = Figure(size=(1000, 900))

    for (i, name) in enumerate(["α_PAR", "τ_PAR", "α_NIR", "τ_NIR"])
        row = (i - 1) ÷ 2 + 1
        col = (i - 1) % 2 + 1
        res = results[name]
        resid = res.residuals

        ax = Axis(fig[row, col],
            xlabel = "Residual",
            ylabel = "Count",
            title  = "$name residuals (std = $(round(std(resid); digits=4)))")

        hist!(ax, resid; bins=50, color=(:steelblue, 0.7))
        vlines!(ax, [0.0]; color=:red, linewidth=1.5, linestyle=:dash)
    end

    save(joinpath(PLOT_DIR, "03_residual_histograms.png"), fig, px_per_unit=3)
    println("  Saved: 03_residual_histograms.png")
    return fig
end

# ═════════════════════════════════════════════════════════════════════════════
#  Plot 4 — Residual maps (4-panel global)
# ═════════════════════════════════════════════════════════════════════════════

"""
Rebuild a 2D map from flat residuals. We need to track which pixels were
used in the fit, so we re-run the masking logic from fit_leaf_optics.jl.
"""
function residuals_to_map(par_r, par_t, nir_r, nir_t, lma, vcmax, qc_mask, residuals_vec)
    α_par_map = absorptance(par_r, par_t)
    α_nir_map = absorptance(nir_r, nir_t)

    idx = findall(eachindex(par_r)) do k
        qc_mask[k] &&
        !isnan(lma[k]) && !isnan(vcmax[k]) &&
        par_r[k] >= 0 && par_t[k] >= 0 &&
        nir_r[k] >= 0 && nir_t[k] >= 0 &&
        α_par_map[k] >= 0 && α_nir_map[k] >= 0 &&
        lma[k] > 0 && vcmax[k] > 0
    end

    resid_map = fill(NaN, size(par_r))
    for (j, k) in enumerate(idx)
        resid_map[k] = residuals_vec[j]
    end
    return resid_map
end

function plot_residual_maps(lon, lat, par_r, par_t, nir_r, nir_t, lma, vcmax, qc_mask, results)
    fig = Figure(size=(1200, 900))

    for (i, name) in enumerate(["α_PAR", "τ_PAR", "α_NIR", "τ_NIR"])
        row = (i - 1) ÷ 2 + 1
        col = (i - 1) % 2 + 1
        res = results[name]

        resid_map = residuals_to_map(par_r, par_t, nir_r, nir_t, lma, vcmax, qc_mask, res.residuals)

        # Shift lon from [0,360) → [-180,180) for display
        lon_display = [l > 180 ? l - 360 : l for l in lon]
        perm = sortperm(lon_display)
        lon_sorted = lon_display[perm]
        data_sorted = resid_map[perm, :]

        ax = Axis(fig[row, col],
            xlabel = "Longitude",
            ylabel = "Latitude",
            title  = "$name residuals")

        clim = maximum(abs.(filter(!isnan, res.residuals))) * 0.8

        hm = heatmap!(ax, lon_sorted, lat, data_sorted;
            colormap = :RdBu_10,
            colorrange = (-clim, clim),
            nan_color = :gray90)

        Colorbar(fig[row, col][1,2], hm; width=12)
    end

    save(joinpath(PLOT_DIR, "04_residual_maps.png"), fig, px_per_unit=3)
    println("  Saved: 04_residual_maps.png")
    return fig
end

# ═════════════════════════════════════════════════════════════════════════════
#  Plot 5 — Marginal scatter: each target vs each predictor
# ═════════════════════════════════════════════════════════════════════════════

function plot_marginal_scatter(tbl)
    fig = Figure(size=(1200, 1000))

    targets = [
        ("α_PAR", tbl.α_par),
        ("τ_PAR", tbl.τ_par),
        ("α_NIR", tbl.α_nir),
        ("τ_NIR", tbl.τ_nir),
    ]

    predictors = [
        ("log(Vcmax25)", tbl.log_vcmax),
        ("log(LMA)",     tbl.log_lma),
    ]

    for (i, (tname, tvals)) in enumerate(targets)
        for (j, (pname, pvals)) in enumerate(predictors)
            ax = Axis(fig[i, j],
                xlabel = pname,
                ylabel = tname,
                title  = "$tname vs $pname (r=$(round(cor(pvals, tvals); digits=3)))")

            scatter!(ax, pvals, tvals;
                markersize = 2, alpha = 0.3, color = :steelblue)
        end
    end

    save(joinpath(PLOT_DIR, "05_marginal_scatter.png"), fig, px_per_unit=3)
    println("  Saved: 05_marginal_scatter.png")
    return fig
end

# ═════════════════════════════════════════════════════════════════════════════
#  Main
# ═════════════════════════════════════════════════════════════════════════════

function plot_main()
    println("Running fit pipeline...")
    out = main()  # from fit_leaf_optics.jl

    tbl     = out.tbl
    results = out.results
    lon     = out.lon
    lat     = out.lat

    println("\nGenerating plots → $(PLOT_DIR)/\n")

    # Need to reconstruct qc_mask and grids for residual maps
    rt = load_refl_tran_annual()
    qc_mask = quality_check_optics(rt.par_refl, rt.par_tran, rt.nir_refl, rt.nir_tran)

    plot_predictor_space(tbl)
    plot_obs_vs_fitted(tbl, results)
    plot_residual_histograms(results)
    plot_residual_maps(lon, lat, rt.par_refl, rt.par_tran, rt.nir_refl, rt.nir_tran,
                       out.lma_2p5, out.vcmax_2p5, qc_mask, results)
    plot_marginal_scatter(tbl)

    println("\nAll plots saved to $(PLOT_DIR)/")
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_main()
end
