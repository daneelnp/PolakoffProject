"""
    fit_leaf_optics.jl

Option B fit: 35-year mean traits + annual-mean optics on a 2.5° grid.

Model (log-linear, per Claude's recommendation):
    target = β₀ + β₁·log(Vcmax25) + β₂·log(LMA)

for each target ∈ {α_PAR, τ_PAR, α_NIR, τ_NIR}.

Pipeline:
  1. Load clm_refl_tran monthly climatology → annual mean (144×96)
     • Mask fill value (-1.0) → NaN
     • Quality check: refl + tran ≤ 1 per waveband
  2. Load LMA_deciduous, LMA_evergreen, Vcmax25 → 35-year temporal mean (720×360)
  3. Combine deciduous + evergreen LMA (mean where both exist)
  4. Shift Dong longitudes [-180,180) → [0,360)
  5. Upscale traits from 0.5° → 2.5° via nearest-neighbor block averaging
  6. Mask to valid pixels in ALL variables, apply physical filters
  7. Fit log-linear model, report R², RMSE, VIF

Usage:
    julia --project=. src/fit_leaf_optics.jl
"""

using NCDatasets
using Statistics
using LinearAlgebra

# ── Paths ────────────────────────────────────────────────────────────────────
const DATA_DIR = joinpath(@__DIR__, "..", "data")

const FILE_REFL_TRAN = joinpath(DATA_DIR, "clm_refl_tran_1m_weighted.nc")
const FILE_LMA_DECID = joinpath(DATA_DIR, "TS_LMA_decidudous.nc")
const FILE_LMA_EVER  = joinpath(DATA_DIR, "TS_LMA_evergreen.nc")
const FILE_VCMAX25   = joinpath(DATA_DIR, "TS_Vcmax25.nc")

const YEAR_START = 1982
const YEAR_END   = 2016
const N_YEARS    = YEAR_END - YEAR_START + 1   # 35

# ═════════════════════════════════════════════════════════════════════════════
#  Step 1 — Load clm_refl_tran → annual mean
# ═════════════════════════════════════════════════════════════════════════════

"""
Load clm_refl_tran, replace -1.0 fill with NaN, average across 12 months.
Returns NamedTuple with (lon, lat, par_refl, par_tran, nir_refl, nir_tran),
each 2D array (144, 96).
"""
function load_refl_tran_annual()
    lon = lat = nothing
    par_r = par_t = nir_r = nir_t = nothing

    NCDataset(FILE_REFL_TRAN, "r") do ds
        lon   = ds["lon"][:]
        lat   = ds["lat"][:]
        par_r = ds["par_refl"][:,:,:]   # (144, 96, 12)
        par_t = ds["par_tran"][:,:,:]
        nir_r = ds["nir_refl"][:,:,:]
        nir_t = ds["nir_tran"][:,:,:]
    end

    # Replace fill (-1.0) with NaN
    for arr in (par_r, par_t, nir_r, nir_t)
        arr[arr .== -1.0] .= NaN
    end

    # NaN-aware mean over months (dim 3)
    function nanmean3(A)
        out = similar(A, size(A,1), size(A,2))
        for j in axes(A,2), i in axes(A,1)
            vals = filter(!isnan, A[i,j,:])
            out[i,j] = isempty(vals) ? NaN : mean(vals)
        end
        return out
    end

    return (
        lon      = lon,
        lat      = lat,
        par_refl = nanmean3(par_r),
        par_tran = nanmean3(par_t),
        nir_refl = nanmean3(nir_r),
        nir_tran = nanmean3(nir_t),
    )
end

# ═════════════════════════════════════════════════════════════════════════════
#  Step 2 — Load Dong datasets → 35-year temporal mean
# ═════════════════════════════════════════════════════════════════════════════

"""
Load a Dong 0.5° dataset and compute the temporal mean across all 35 years
(z dimension). Returns (lon, lat, data_2d).
"""
function load_dong_mean(filepath::String)
    lon = lat = nothing
    data = nothing

    NCDataset(filepath, "r") do ds
        lon  = ds["longitude"][:]       # (720,)
        lat  = ds["latitude"][:]        # (360,)
        data = ds["variable"][:,:,:]    # (720, 360, 35)
    end

    # Convert to Float64, Missing → NaN
    nlon, nlat, ntime = size(data)
    arr = Array{Float64}(undef, nlon, nlat, ntime)
    for i in eachindex(data)
        arr[i] = ismissing(data[i]) ? NaN : Float64(data[i])
    end

    # Temporal mean (dim 3), NaN-aware
    out = Array{Float64}(undef, nlon, nlat)
    for j in 1:nlat, i in 1:nlon
        vals = filter(!isnan, arr[i,j,:])
        out[i,j] = isempty(vals) ? NaN : mean(vals)
    end

    return (lon=lon, lat=lat, data=out)
end

# ═════════════════════════════════════════════════════════════════════════════
#  Step 3 — Combine deciduous + evergreen LMA
# ═════════════════════════════════════════════════════════════════════════════

"""
Combine deciduous and evergreen LMA maps.
Where both valid → mean; where one valid → use it; else NaN.
"""
function combine_lma(decid::Matrix{Float64}, ever::Matrix{Float64})
    @assert size(decid) == size(ever)
    out = similar(decid)
    for i in eachindex(decid)
        d = !isnan(decid[i])
        e = !isnan(ever[i])
        if d && e
            out[i] = 0.5 * (decid[i] + ever[i])
        elseif d
            out[i] = decid[i]
        elseif e
            out[i] = ever[i]
        else
            out[i] = NaN
        end
    end
    return out
end

# ═════════════════════════════════════════════════════════════════════════════
#  Step 4 — Shift longitudes [-180, 180) → [0, 360)
# ═════════════════════════════════════════════════════════════════════════════

function shift_lon_0_360(lon::Vector{Float64}, data::Matrix{Float64})
    new_lon = mod.(lon, 360.0)
    perm    = sortperm(new_lon)
    return new_lon[perm], data[perm, :]
end

# ═════════════════════════════════════════════════════════════════════════════
#  Step 5 — Upscale 0.5° → 2.5° (nearest-neighbor binning + block average)
# ═════════════════════════════════════════════════════════════════════════════

function upscale_to_2p5(lon_05, lat_05, data_05, lon_2p5, lat_2p5)
    nlon = length(lon_2p5)
    nlat = length(lat_2p5)

    lon_idx = assign_bins(lon_05, lon_2p5)
    lat_idx = assign_bins(lat_05, lat_2p5)

    sums   = zeros(nlon, nlat)
    counts = zeros(Int, nlon, nlat)

    for j in eachindex(lat_05)
        tj = lat_idx[j]
        tj == 0 && continue
        for i in eachindex(lon_05)
            ti = lon_idx[i]
            ti == 0 && continue
            v = data_05[i, j]
            if !isnan(v)
                sums[ti, tj]   += v
                counts[ti, tj] += 1
            end
        end
    end

    out = fill(NaN, nlon, nlat)
    for k in eachindex(out)
        out[k] = counts[k] > 0 ? sums[k] / counts[k] : NaN
    end
    return out
end

function assign_bins(src::Vector{<:Real}, tgt::Vector{<:Real})
    idx = zeros(Int, length(src))
    Δ = length(tgt) > 1 ? 0.5 * median(diff(tgt)) : Inf
    for i in eachindex(src)
        best_d = Inf
        best_j = 0
        for j in eachindex(tgt)
            d = abs(src[i] - tgt[j])
            if d < best_d
                best_d = d
                best_j = j
            end
        end
        idx[i] = best_d <= Δ ? best_j : 0
    end
    return idx
end

# ═════════════════════════════════════════════════════════════════════════════
#  Step 6 — Quality checks + build pixel table
# ═════════════════════════════════════════════════════════════════════════════

absorptance(refl, tran) = 1.0 .- refl .- tran

"""
Quality checks on α/τ data (per Claude's Q4 recommendations).
Returns a boolean mask (true = pass).
"""
function quality_check_optics(par_r, par_t, nir_r, nir_t; tol=0.05)
    nlon, nlat = size(par_r)
    mask = trues(nlon, nlat)
    n_refl_tran_fail = 0
    n_zero_spike     = 0

    for k in eachindex(par_r)
        if isnan(par_r[k]) || isnan(par_t[k]) || isnan(nir_r[k]) || isnan(nir_t[k])
            mask[k] = false
            continue
        end

        par_sum = par_r[k] + par_t[k]
        nir_sum = nir_r[k] + nir_t[k]
        if par_sum > 1.0 + tol || nir_sum > 1.0 + tol
            mask[k] = false
            n_refl_tran_fail += 1
            continue
        end

        if par_r[k] == 0.0 || par_t[k] == 0.0 || nir_r[k] == 0.0 || nir_t[k] == 0.0
            n_zero_spike += 1
        end
    end

    println("    QC: $(n_refl_tran_fail) pixels failed refl+tran ≤ 1 test")
    println("    QC: $(n_zero_spike) pixels have exact-zero (kept)")
    println("    QC: $(count(mask)) pixels pass (of $(length(mask)))")
    return mask
end

"""
Build paired pixel table with log-transformed predictors.
"""
function build_pixel_table(par_r, par_t, nir_r, nir_t, lma, vcmax, qc_mask)
    @assert size(par_r) == size(lma) == size(vcmax) == size(qc_mask)

    α_par = absorptance(par_r, par_t)
    α_nir = absorptance(nir_r, nir_t)

    idx = findall(eachindex(par_r)) do k
        qc_mask[k] &&
        !isnan(lma[k]) && !isnan(vcmax[k]) &&
        par_r[k] >= 0 && par_t[k] >= 0 &&
        nir_r[k] >= 0 && nir_t[k] >= 0 &&
        α_par[k] >= 0 && α_nir[k] >= 0 &&
        lma[k] > 0 && vcmax[k] > 0
    end

    return (
        α_par     = α_par[idx],
        τ_par     = par_t[idx],
        α_nir     = α_nir[idx],
        τ_nir     = nir_t[idx],
        ρ_par     = par_r[idx],
        ρ_nir     = nir_r[idx],
        lma       = lma[idx],
        vcmax     = vcmax[idx],
        log_lma   = log.(lma[idx]),
        log_vcmax = log.(vcmax[idx]),
        npix      = length(idx),
    )
end

# ═════════════════════════════════════════════════════════════════════════════
#  Step 7 — Log-linear fit + diagnostics
# ═════════════════════════════════════════════════════════════════════════════

ols(X, y) = X \ y

function r_squared(y, ŷ)
    ss_res = sum((y .- ŷ).^2)
    ss_tot = sum((y .- mean(y)).^2)
    return 1.0 - ss_res / ss_tot
end

"""
Variance Inflation Factor for each predictor column (excluding intercept).
"""
function compute_vif(X_pred)
    p = size(X_pred, 2)
    vifs = zeros(p)
    for j in 1:p
        y_j = X_pred[:, j]
        X_j = hcat(ones(size(X_pred, 1)), X_pred[:, setdiff(1:p, j)])
        ŷ_j = X_j * ols(X_j, y_j)
        R²_j = r_squared(y_j, ŷ_j)
        vifs[j] = 1.0 / max(1.0 - R²_j, 1e-12)
    end
    return vifs
end

"""
Fit log-linear models and report diagnostics.
"""
function fit_log_linear(tbl)
    n = tbl.npix
    println("\n  Fitting log-linear model on $n pixels")
    println("  Model: target = β₀ + β₁·log(Vcmax25) + β₂·log(LMA)\n")

    X      = hcat(ones(n), tbl.log_vcmax, tbl.log_lma)
    X_pred = hcat(tbl.log_vcmax, tbl.log_lma)

    # VIF
    vifs = compute_vif(X_pred)
    println("  ── Collinearity Check ──")
    println("    VIF(log_Vcmax25) = $(round(vifs[1]; digits=2))")
    println("    VIF(log_LMA)     = $(round(vifs[2]; digits=2))")
    if maximum(vifs) > 5
        println("    ⚠  VIF > 5 — predictors are collinear!")
    else
        println("    ✓  VIF < 5 — no severe collinearity")
    end
    println("    Pearson r(log_Vcmax, log_LMA) = $(round(cor(tbl.log_vcmax, tbl.log_lma); digits=4))")
    println()

    targets = [
        ("α_PAR", tbl.α_par),
        ("τ_PAR", tbl.τ_par),
        ("α_NIR", tbl.α_nir),
        ("τ_NIR", tbl.τ_nir),
    ]

    results = Dict{String, NamedTuple}()

    for (name, y) in targets
        β    = ols(X, y)
        ŷ    = X * β
        R²   = r_squared(y, ŷ)
        rmse = sqrt(mean((y .- ŷ).^2))
        resid = y .- ŷ

        println("  ── $name ──")
        println("    β₀ (intercept)   = $(round(β[1]; digits=6))")
        println("    β₁ (log Vcmax25) = $(round(β[2]; sigdigits=4))")
        println("    β₂ (log LMA)     = $(round(β[3]; sigdigits=4))")
        println("    R²               = $(round(R²; digits=4))")
        println("    RMSE             = $(round(rmse; digits=6))")
        println("    y range          : $(round(minimum(y); digits=4)) – $(round(maximum(y); digits=4))")
        println("    residual range   : $(round(minimum(resid); digits=4)) – $(round(maximum(resid); digits=4))")
        println("    residual std     : $(round(std(resid); digits=6))")
        println()

        results[name] = (β=β, R²=R², rmse=rmse, residuals=resid, fitted=ŷ)
    end

    return results
end

# ═════════════════════════════════════════════════════════════════════════════
#  Main
# ═════════════════════════════════════════════════════════════════════════════

function main()
    println("="^70)
    println("  Leaf Optics Fit — Option B (35-year mean, log-linear)")
    println("="^70)

    # Step 1
    println("\n[1/6] Loading clm_refl_tran → annual mean...")
    rt = load_refl_tran_annual()
    println("  Grid: $(length(rt.lon)) lon × $(length(rt.lat)) lat")

    # Step 2
    println("\n[2/6] Loading Dong trait data (35-year means, 1982–2016)...")
    decid = load_dong_mean(FILE_LMA_DECID)
    ever  = load_dong_mean(FILE_LMA_EVER)
    vcmax = load_dong_mean(FILE_VCMAX25)
    println("  Dong grid: $(length(decid.lon)) lon × $(length(decid.lat)) lat")
    println("  LMA deciduous — valid: $(count(!isnan, decid.data))")
    println("  LMA evergreen — valid: $(count(!isnan, ever.data))")
    println("  Vcmax25       — valid: $(count(!isnan, vcmax.data))")

    # Step 3
    println("\n[3/6] Combining deciduous + evergreen LMA...")
    lma_raw = combine_lma(decid.data, ever.data)
    println("  Combined LMA: $(count(!isnan, lma_raw)) valid (of $(length(lma_raw)))")

    # Step 4
    println("\n[4/6] Shifting Dong longitudes to [0, 360)...")
    lon_shift, lma_shift   = shift_lon_0_360(decid.lon, lma_raw)
    _,         vcmax_shift = shift_lon_0_360(vcmax.lon, vcmax.data)
    println("  Lon range: $(round(minimum(lon_shift);digits=2)) → $(round(maximum(lon_shift);digits=2))")

    # Step 5
    println("\n[5/6] Upscaling traits 0.5° → 2.5°...")
    lma_2p5   = upscale_to_2p5(lon_shift, decid.lat, lma_shift,   rt.lon, rt.lat)
    vcmax_2p5 = upscale_to_2p5(lon_shift, decid.lat, vcmax_shift, rt.lon, rt.lat)
    println("  LMA   on 2.5°: $(count(!isnan, lma_2p5)) valid / $(length(lma_2p5))")
    println("  Vcmax on 2.5°: $(count(!isnan, vcmax_2p5)) valid / $(length(vcmax_2p5))")

    # Step 6
    println("\n[6/6] Quality checks & building pixel table...")
    qc_mask = quality_check_optics(rt.par_refl, rt.par_tran, rt.nir_refl, rt.nir_tran)
    tbl = build_pixel_table(rt.par_refl, rt.par_tran, rt.nir_refl, rt.nir_tran,
                            lma_2p5, vcmax_2p5, qc_mask)
    println("  Final valid pixels: $(tbl.npix)")

    # Data summary
    println("\n  ── Data Summary ──")
    for (name, vals) in [("α_PAR",      tbl.α_par),    ("τ_PAR",      tbl.τ_par),
                          ("α_NIR",      tbl.α_nir),    ("τ_NIR",      tbl.τ_nir),
                          ("LMA",        tbl.lma),      ("Vcmax25",    tbl.vcmax),
                          ("log(LMA)",   tbl.log_lma),  ("log(Vcmax)", tbl.log_vcmax)]
        println("    $name : min=$(round(minimum(vals);digits=4)), max=$(round(maximum(vals);digits=4)), mean=$(round(mean(vals);digits=4)), std=$(round(std(vals);digits=4))")
    end

    # Step 7
    results = fit_log_linear(tbl)

    println("="^70)
    println("  Done.")
    println("="^70)

    return (; tbl, results, lon=rt.lon, lat=rt.lat, lma_2p5, vcmax_2p5)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end