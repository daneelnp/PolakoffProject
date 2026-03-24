# Phase 2 Data Report — Leaf Optics Project
## For: Claude (Research/Theory Assistant)
## Date: March 4, 2026

---

## Project Goal

Learn the function `f` such that:

```
α_leaf(PAR, NIR), τ_leaf(PAR, NIR) = f(Vcmax25, LMA)
```

Where:
- **α** = leaf absorptance, **τ** = leaf transmittance (reflectance ρ = 1 − α − τ)
- **LMA** = Leaf Mass per Area (g/m²)
- **Vcmax25** = maximum carboxylation rate at 25°C (µmol/m²/s)
- Constraints: all outputs ∈ [0, 1]; α + τ ≤ 1 per waveband

This function will replace the current hardcoded leaf optics in CliMA Land.

---

## What We Did

1. **Loaded and inspected all four NetCDF datasets** (`load_data.jl`)
2. **Identified the structure mismatch** between datasets (see below)
3. **Wrote a fitting pipeline** (`fit_leaf_optics.jl`) that picks one year, upscales, aligns, and fits a linear model — not yet executed

---

## Dataset Summary

| Dataset | File | Grid | Lon Convention | Time Dimension | Source |
|---------|------|------|----------------|----------------|--------|
| Leaf α & τ (PAR, NIR) | `clm_refl_tran_1m_weighted.nc` | 144×96 (2.5°) | [0, 360) | `ind` 1–12 = **monthly climatology** (no year info) | Wang & Braghiere et al. 2025 |
| Deciduous LMA | `TS_LMA_decidudous.nc` | 720×360 (0.5°) | [−180, 180) | `z` 1–35 = **annual, 1982–2016** | Dong et al. 2023 |
| Evergreen LMA | `TS_LMA_evergreen.nc` | 720×360 (0.5°) | [−180, 180) | `z` 1–35 = **annual, 1982–2016** | Dong et al. 2023 |
| Vcmax25 | `TS_Vcmax25.nc` | 720×360 (0.5°) | [−180, 180) | `z` 1–35 = **annual, 1982–2016** | Dong et al. 2023 |

### Key Variable Details

| Variable | Shape | Missing | Range | Fill value |
|----------|-------|---------|-------|------------|
| `par_refl` | (144, 96, 12) | 0 nominal, but -1.0 = ocean/fill | valid: 0 – 0.43 | -1.0 |
| `par_tran` | (144, 96, 12) | same | valid: 0 – 0.47 | -1.0 |
| `nir_refl` | (144, 96, 12) | same | valid: 0 – 0.40 | -1.0 |
| `nir_tran` | (144, 96, 12) | same | valid: 0 – 0.48 | -1.0 |
| LMA (deciduous) | (720, 360, 35) | 79% missing | 0.56 – 203 g/m² | NaN |
| LMA (evergreen) | (720, 360, 35) | 79% missing | 39.9 – 291 g/m² | NaN |
| Vcmax25 | (720, 360, 35) | 79% missing | 11.1 – 170 µmol/m²/s | NaN |

---

## The Core Problem: Temporal Mismatch

The **target variables** (α, τ) have **monthly** resolution but **no interannual variability** — they are a single climatological seasonal cycle.

The **predictor variables** (LMA, Vcmax25) have **annual** resolution across **35 years** but **no monthly variability**.

There is no dataset that has both monthly AND interannual resolution. This constrains how we can pair the data.

---

## Options for Aligning the Data

### Option A: Single Annual Snapshot (simplest — currently implemented)

- **Pick one year** (e.g., 2000) from LMA/Vcmax
- **Average** clm_refl_tran across all 12 months → one annual-mean α/τ map
- **Upscale** LMA/Vcmax from 0.5° → 2.5° to match α/τ grid
- **Fit** on ~4,000–6,000 valid land pixels

**Pros:** Simple, clean, no assumptions  
**Cons:** Throws away seasonal info, uses only 1/35th of predictor data, clm_refl_tran seasonal cycle is lost  

### Option B: Multi-Year Annual Mean (climatological average)

- **Average LMA/Vcmax across all 35 years** → one mean trait map each
- **Average α/τ across 12 months** → one annual mean map
- **Upscale** and fit

**Pros:** Uses all years of trait data (more robust mean), still simple  
**Cons:** Still loses seasonal + interannual variability  

### Option C: Monthly Climatology Matching

- **Average LMA/Vcmax across all 35 years** → one mean trait map each
- **Keep** clm_refl_tran as 12 monthly maps
- **Replicate** the single LMA/Vcmax map across 12 months
- **Fit** with 12× more data points (same predictors, different targets per month)

**Pros:** Preserves the seasonal cycle of α/τ, 12× more training points  
**Cons:** LMA/Vcmax are constant across months — the model must learn seasonality purely from spatial variation. This may artificially inflate R² due to repeated predictors. The model can't explain seasonal variation in α/τ because the predictors don't vary seasonally.

### Option D: Year-by-Year with Repeated Climatology

- **For each year** (1982–2016), pair that year's LMA/Vcmax with the same annual-mean α/τ
- **Stack** all 35 years → 35× more points than Option A
- **Fit** on the combined pool

**Pros:** Uses all interannual trait variation, large sample size  
**Cons:** Target (α/τ) is identical across years — the model is learning whether trait *level differences* across years map to the *same* spatial pattern of α/τ. Risk of pseudoreplication (same target, slightly different predictors).

### Option E: Monthly × Yearly (maximum data, most complex)

- **Replicate** each year's LMA/Vcmax across 12 months
- **Pair** with the 12-month climatological α/τ
- **Fit** on 35 × 12 = 420 maps' worth of pixels

**Pros:** Maximum data volume  
**Cons:** Severe pseudoreplication — same 12-month α/τ target repeated 35 times. Statistical degrees of freedom are much lower than the point count suggests. Likely misleading R².

---

## Spatial Resolution Options

| Strategy | Direction | Method | Result |
|----------|-----------|--------|--------|
| **Upscale traits** (recommended) | 0.5° → 2.5° | 5×5 block average | 144×96 grid, matches α/τ natively |
| **Downscale α/τ** | 2.5° → 0.5° | Bilinear interpolation | 720×360 grid, more points but interpolated targets |
| **Intermediate grid** | Both → e.g. 1° | Regrid both | Compromise, requires two regridding steps |

Upscaling the traits is preferred because:
- Averaging is physically defensible (area-mean trait values)
- Avoids creating artificial spatial detail in the target via interpolation
- Fewer pixels = faster fitting, less spatial autocorrelation

---

## Longitude Convention

- clm_refl_tran uses **[0, 360)**: lon = 1.25, 3.75, ..., 358.75
- Dong datasets use **[−180, 180)**: lon = −179.75, −179.25, ..., 179.75

Our pipeline shifts Dong → [0, 360) via `mod(lon, 360)` before upscaling.

---

## Combining Deciduous + Evergreen LMA

Each pixel may have deciduous LMA, evergreen LMA, both, or neither.  
Current approach: **mean where both exist, use whichever is available otherwise**.

Alternative approaches to consider:
- **Weighted average** by PFT fraction (requires a land cover map)
- **Treat separately** — fit two models, one per vegetation type
- **Use only one** (e.g., dominant PFT per pixel)

---

## Recommendation for First Fit

**Option A** (single year, annual mean, upscaled to 2.5°) is ready to run.  
It's the cleanest baseline. After evaluating R² and residual patterns, we can decide whether to:
- Expand to Option B or D for more data
- Try Option C to explore seasonal signal
- Switch to nonlinear models (polynomial, neural net) if linear R² is poor

---

## Open Questions for Claude

1. Given that α/τ is climatological and LMA/Vcmax is annual — what is the most statistically defensible way to pair these for regression? Is there a risk of ecological fallacy in any of the options?

2. Should we combine deciduous + evergreen LMA, or treat them as separate models? The physics of leaf optics differs between the two types.

3. Is a linear model `α = β₀ + β₁·Vcmax + β₂·LMA` a reasonable starting point, or should we go directly to something nonlinear given the known biophysics (e.g., PROSPECT model relates leaf structure parameter N and pigment content to spectra)?

4. The clm_refl_tran data uses -1.0 as fill for ocean pixels (not standard CF missing). After masking these, the valid range of reflectance/transmittance is physically reasonable (0–0.48). Should we apply additional quality filters?

5. For the final CliMA Land implementation, `f` needs to run at every timestep. What functional form is most practical — lookup table, polynomial, or a small neural net?

---
---

## Claude's Response & Agreed Decisions (March 4, 2026)

### Reframing the Temporal Mismatch

The temporal mismatch is less troubling than it appears. We are not trying to learn how optics *change over time* — we are learning whether **spatial variation** in traits explains **spatial variation** in optics across the globe. The Wang 2025 climatology gives the "typical" optical property per location; Dong's multi-year mean gives the "typical" trait value per pixel. This is a coherent **cross-sectional regression**. The mismatch only becomes problematic if we later try to predict *temporal changes* in optics from year-to-year trait changes — which is not the current goal.

### Decision: Temporal Alignment → Option B

**Use Option B**: 35-year mean traits + annual-mean optics, upscaled to 2.5°.

Rationale for rejecting the alternatives:
- **Options D & E** (multi-year stacking): reuse the same target across years → pseudoreplication, inflated confidence intervals. Avoid.
- **Option C** (monthly climatology): predictors are constant across months but targets vary seasonally → model is forced to explain seasonal α/τ variation from spatial trait variation only, which is spurious. Avoid for main fit (worth running as a diagnostic).
- **Option A** (single year): wastes 34 years of trait data for no reason.

**Ecological fallacy risk**: mitigated since we work at the pixel level, not biome aggregates. The real risk is **spatial autocorrelation** — neighboring pixels share similar traits AND optics because they're in the same biome, inflating R² and shrinking standard errors. Recommended diagnostics:
- Moran's I test on residuals
- Spatially-blocked cross-validation (hold out entire regions, not random pixels)

### Decision: LMA Treatment → Combined Baseline First

Start with **combined deciduous + evergreen LMA** (mean where both exist), then check residuals colored by vegetation type. If deciduous and evergreen residuals are systematically offset in opposite directions → split into separate models.

Physical justification for eventually splitting: evergreen leaves have high LMA from structural carbon investment (thick cuticle, dense mesophyll); deciduous leaves have high LMA for different reasons (nutrient storage). The chlorophyll-per-LMA relationship differs, and since PAR absorptance is driven by chlorophyll, lumping them risks fitting a bimodal distribution with a single line.

Note: check whether Wang 2025 `clm_refl_tran_1m_weighted.nc` already blends PFTs in its "weighted" composite — if so, combining LMA types is *more consistent* with how the target was constructed.

### Decision: Model Form → Log-Linear

Start with **log-transformed predictors** in a linear regression:

```
α_PAR = β₀ + β₁·log(Vcmax25) + β₂·log(LMA)
```

Physical justification:
- Chlorophyll absorption follows Beer-Lambert law → transmittance is roughly exponential in pigment content → absorptance saturates at high chlorophyll. Since Vcmax25 proxies for chlorophyll/nitrogen, `α_PAR ~ Vcmax25` is concave. `log(Vcmax25)` linearizes this.
- LMA spans nearly three orders of magnitude (0.56–291 g/m²) → log transform makes the distribution more symmetric and the fit more stable.

This is still linear in coefficients (easy to fit and interpret) but encodes physically reasonable curvature. If residuals remain structured → add `log(Vcmax25)·log(LMA)` interaction term or degree-2 polynomial in log space before going to random forest.

**Collinearity check**: Vcmax25 and LMA are correlated (leaf economics spectrum). Compute VIF after fitting — if VIF > 5, the two predictors are collinear and may not identify independent effects.

### Decision: Quality Filters

After masking -1.0 fills in Wang 2025:
1. Verify `refl + tran ≤ 1` per waveband everywhere. Mask pixels where this fails by > 0.05.
2. Check for spikes at exactly 0.0 in α or τ — could be edge artifacts from the weighting procedure.

### Decision: Functional Form for CliMA Land → Polynomial in Log Space

A **polynomial in log-transformed predictors** (degree ≤ 2) is the target:

```julia
α_PAR = f(log(Vcmax25), log(LMA))   # 5–6 coefficients
```

Runs in nanoseconds per cell, is differentiable (critical for adjoint-based calibration), trivially portable to Julia, no external dependencies. Neural nets introduce serialization overhead, extrapolation issues, and reviewer friction. Lookup tables require interpolation logic and don't extrapolate gracefully.

### Agreed Next Steps

1. ~~Compute 35-year mean of each Dong variable~~ ✅
2. ~~Compute annual mean of Wang 2025 (masking -1.0 fills)~~ ✅
3. ~~Upscale Dong to 2.5°~~ ✅
4. ~~Plot scatter of `log(Vcmax25)` vs. `α_PAR`~~ ✅
5. ~~Fit the log-linear model, report R², RMSE, VIF~~ ✅
6. Check residuals by vegetation type — partially done (see residual maps below)

---
---

## Option B Fit Results (March 4, 2026)

### Pipeline Summary

| Step | Detail |
|------|--------|
| Target (α/τ) | Wang 2025 `clm_refl_tran_1m_weighted.nc`, annual mean of 12-month climatology |
| Predictors | Dong 2023: 35-year mean (1982–2016) of LMA and Vcmax25 |
| LMA combination | Mean of deciduous + evergreen where both exist; single PFT where only one valid |
| Spatial resolution | Traits upscaled 0.5° → 2.5° via nearest-neighbor block averaging |
| Longitude alignment | Dong [-180, 180) shifted to [0, 360) to match Wang grid |
| Quality filter | refl + tran ≤ 1.05 per waveband (0 pixels failed) |
| Valid pixels | **3,776** (out of 13,824 grid cells) |

### Data Summary

| Variable | Min | Max | Mean | Std |
|----------|-----|-----|------|-----|
| α_PAR | 0.296 | 0.926 | 0.739 | 0.162 |
| τ_PAR | 0.025 | 0.369 | 0.132 | 0.088 |
| α_NIR | 0.125 | 0.248 | 0.178 | 0.022 |
| τ_NIR | 0.413 | 0.478 | 0.451 | 0.012 |
| LMA (g/m²) | 21.2 | 197.4 | 90.2 | 27.7 |
| Vcmax25 (µmol/m²/s) | 12.9 | 143.2 | 53.7 | 19.5 |
| log(LMA) | 3.05 | 5.29 | 4.45 | 0.31 |
| log(Vcmax25) | 2.56 | 4.96 | 3.91 | 0.42 |

### Collinearity Check

| Metric | Value |
|--------|-------|
| VIF(log Vcmax25) | 2.44 |
| VIF(log LMA) | 2.44 |
| Pearson r(log Vcmax, log LMA) | 0.769 |

**Verdict:** VIF < 5 — no severe collinearity. The predictors are correlated (leaf economics spectrum) but each carries independent information.

### Log-Linear Regression Results

Model: `target = β₀ + β₁·log(Vcmax25) + β₂·log(LMA)`

| Target | β₀ | β₁ (log Vcmax) | β₂ (log LMA) | R² | RMSE |
|--------|-----|-----------------|---------------|-----|------|
| α_PAR | 0.305 | −0.404 | +0.452 | **0.453** | 0.119 |
| τ_PAR | 0.365 | +0.222 | −0.247 | **0.462** | 0.064 |
| α_NIR | 0.081 | −0.050 | +0.065 | **0.409** | 0.017 |
| τ_NIR | 0.501 | +0.026 | −0.034 | **0.407** | 0.009 |

### Key Observations

1. **R² ≈ 0.41–0.46** — the two-predictor log-linear model explains ~45% of the spatial variance. This is a meaningful signal but substantial residual structure remains.

2. **Coefficient signs are physically consistent:**
   - Higher LMA → higher absorptance, lower transmittance (denser leaves absorb more light)
   - Higher Vcmax25 → lower absorptance, higher transmittance in PAR — this seems counterintuitive (more photosynthetic capacity = less absorption?) and may reflect that high-Vcmax leaves in this dataset tend to be thinner/less pigmented per unit area, or it may be an artifact of the LMA-Vcmax collinearity

3. **PAR has much wider dynamic range than NIR** — α_PAR ranges from 0.30 to 0.93 while α_NIR only ranges from 0.13 to 0.25. This is expected: PAR absorptance is driven by chlorophyll (highly variable), while NIR is dominated by leaf internal structure (less variable).

4. **RMSE in context:**
   - α_PAR RMSE = 0.119 on a range of 0.63 → relative error ~19%
   - τ_PAR RMSE = 0.064 on a range of 0.34 → relative error ~19%
   - α_NIR RMSE = 0.017 on a range of 0.12 → relative error ~14%
   - τ_NIR RMSE = 0.009 on a range of 0.065 → relative error ~14%

5. **Quality checks passed cleanly** — zero pixels had refl + tran > 1.05, zero exact-zero spikes.

---

## Diagnostic Plots

All plots saved to `phase2/plots/`.

### Plot 1: Predictor Space (`01_predictor_space.png`)

**What it shows:** Scatter of log(Vcmax25) vs log(LMA) for all 3,776 valid pixels, colored by α_PAR.

**What to look for:**
- The correlation between predictors (r = 0.77) — visible as a positive trend
- Whether α_PAR varies smoothly across the predictor space or has clusters/discontinuities
- Whether the color gradient is aligned with one predictor more than the other (would indicate which predictor dominates)

### Plot 2: Observed vs Fitted (`02_obs_vs_fitted.png`)

**What it shows:** 4-panel scatter of observed vs model-predicted values for each target. Red dashed line = perfect 1:1 fit.

**What to look for:**
- How tightly points cluster around the 1:1 line
- Whether scatter is homoscedastic (uniform spread) or heteroscedastic (spread varies with fitted value)
- Whether there are systematic biases at extremes (curve away from 1:1 at high/low values → nonlinearity needed)
- Outliers far from the line

### Plot 3: Residual Histograms (`03_residual_histograms.png`)

**What it shows:** 4-panel histograms of residuals (observed − fitted) for each target.

**What to look for:**
- Are residuals centered at zero? (unbiased)
- Are they approximately normal? (important for confidence intervals)
- Heavy tails or skewness (would suggest outliers or model misspecification)
- The α_PAR residuals have a notable right tail (residual up to +1.0) — some pixels have much higher absorptance than predicted

### Plot 4: Residual Maps (`04_residual_maps.png`)

**What it shows:** 4-panel global maps of residuals on the 2.5° grid. Blue = model over-predicts, Red = model under-predicts. Gray = no data.

**What to look for:**
- **Geographic clustering of residuals** — if entire biomes are systematically red or blue, that signals a missing predictor (e.g., vegetation type, soil background, canopy structure)
- Whether tropical, temperate, and boreal regions behave differently
- Whether the residual pattern differs between PAR and NIR (would indicate different missing physics)
- This is the key plot for deciding whether to split by vegetation type (Claude's Q2)

### Plot 5: Marginal Scatter (`05_marginal_scatter.png`)

**What it shows:** 4×2 panel grid — each target vs each predictor individually, with Pearson r annotated.

**What to look for:**
- The shape of each relationship — linear in log space, or does it curve?
- Whether the scatter is tight or diffuse (tight = predictor is informative)
- Comparing r values: which predictor matters more for which target?
- If curves are visible → justification for adding quadratic terms or the interaction term `log(Vcmax)·log(LMA)`

---

## Interpretation & Next Steps for Discussion

### What the R² = 0.45 means

The model explains nearly half the spatial variance using only two leaf traits. For a first-pass cross-sectional regression with global data, this is encouraging — it confirms that LMA and Vcmax25 carry real information about leaf optical properties. However, ~55% of variance is unexplained, which could come from:

- **Missing predictors**: chlorophyll content, leaf water content, leaf thickness (distinct from LMA), canopy clumping effects contaminating the "leaf-level" α/τ estimates
- **Vegetation type effects**: deciduous vs evergreen have different LMA-optics relationships (Claude's point). The residual maps should reveal if this is the dominant source of error.
- **Scale mismatch**: the optics are from a radiative transfer inversion (remote sensing derived), while traits are from an optimality model — both have their own uncertainties
- **Nonlinearity**: the log-linear model may not capture saturation or threshold effects

### Suggested next steps (in priority order)

1. **Examine the residual maps** — if biome-scale patterns are visible, split by vegetation type (deciduous vs evergreen) and refit separately
2. **Try adding the interaction term**: `target = β₀ + β₁·log(V) + β₂·log(L) + β₃·log(V)·log(L)`
3. **Try degree-2 polynomial**: add `log(V)²` and `log(L)²` terms (6-parameter model)
4. If R² improves substantially with the polynomial → that's the production model for CliMA Land
5. If R² plateaus ~0.5 → the remaining variance is likely from missing predictors, and we should document this as a known limitation
