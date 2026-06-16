# Selective Borrowing of External Controls via Uniform Confidence Bands

**Band-gated borrowing** fuses a randomized controlled trial (RCT) with an external-control (EC)
data source by certifying and borrowing only the external units that are exchangeable with the
trial — and then quantifying exactly what the borrowing did.

The method re-points the doubly-robust uniform confidence band of Lee, Okui & Whang (2017) at a
**source-bias function**, inverts it into an equivalence region where the bias is certifiably
negligible, gates per-unit borrowing to that region, and feeds the admitted controls into a
transport-weighted AIPW estimator.


## The problem

We observe i.i.d. `O = (X, S, A, Y)`: covariates `X`, source indicator `S` (1 = trial, 0 = external
control), treatment `A` (randomized in the trial, fixed to 0 externally), and outcome `Y`. External
controls can sharpen the trial's control arm — but only where they are exchangeable with it. Pooling
indiscriminately injects bias; discarding everything forfeits power. This method borrows
*selectively*, with guarantees.

## Core idea

Target the **conditional source bias**

```
δ(x) = ν(x) − μ₀₁(x) = E[Y | X=x, S=0] − E[Y | X=x, S=1, A=0]
```

Under trial randomization and external-control consistency this equals
`E[Y(0) | X, S=0] − E[Y(0) | X, S=1]`, so `{x : δ(x) = 0}` is exactly the region where external
controls agree with the trial. We borrow where `|δ| ≤ ε`, for a substantive tolerance `ε` on the
outcome scale, fixed before looking at the data.

## Pipeline

1. **DR pseudo-outcome** — a cross-fitted, doubly-robust score `ψ_δ` whose conditional mean is
   `δ(x)` (two AIPW estimators differenced).
2. **Undersmoothed local-linear fit** — smooth `ψ_δ` along a low-dimensional axis `V = v(X)` to
   estimate the bias profile `δ̄(v) = E[δ(X) | V=v]`, undersmoothing so the smoothing bias is
   negligible relative to the band's width.
3. **Uniform band** — a simultaneous confidence band for `δ̄(·)` by Gaussian multiplier bootstrap,
   valid over the whole axis (and over a family of axes that share the bootstrap draws).
4. **Certify (invert the band)** — keep grid points whose entire band lies inside `[−ε, ε]`:
   `T̂_ε = { v : |ĝ(v)| + c·s(v) ≤ ε }`. This is a *positive* certificate, not "failure to reject."
5. **Gate per unit** — a unit is borrowed only if its axis value is in range and both flanking grid
   points are certified.
6. **Fuse** — a transport-weighted AIPW estimator targets the trial-population ATE, using the trial's
   treated and control arms plus the admitted external controls. Two weights carry distinct jobs:
   `ρ` decides who enters, `w` transports admitted controls to the trial population.

## Guarantee stack

- **G1 — Certificate.** With probability ≥ 1−α, every gated grid point has `|δ̄(v)| ≤ ε`,
  simultaneously and across all axes sharing the bootstrap supremum.
- **G2 — Bias budget.** A computable worst-case ceiling `B̂_max` on the leaked bias, plus a
  bias-aware confidence interval that widens to absorb it.
- **G3 — Audit.** An assumption-free, Neyman-orthogonal estimate of the bias *actually injected*,
  checked against the budget.

A single-number **power readout** (the median band half-height) reports how demanding the data are.
If the defensible `ε` falls below it, the gate declines to borrow and the estimator collapses to the
trial-only AIPW — no coverage is bought that wasn't earned.

## Assumptions

Trial randomization (known `π_A`); external-control consistency (`Y = Y(0)` for `S=0`); overlap on
the borrowing region; harmonized covariates across sources; and partial exchangeability (agreement
holds on some unknown set). G2 adds a bias-index assumption along the gated axis; **G3 drops it**.

## Repository contents

| File | What it is |
| --- | --- |
| `rct.md` | Full method writeup: estimand, assumptions, the DR pseudo-outcome, the uniform band, the equivalence region and per-unit gate, the fused estimator, and the G1 / G2 / G3 derivations. |
| `band_gated_borrowing.png` | Schematic of the band → equivalence region → gate → fused-estimator pipeline. |

## Status

This is a research repository. It currently provides the method and its derivations (`rct.md`).

## Reference

Lee, Okui & Whang (2017), *Doubly robust uniform confidence band for the conditional average
treatment effect function* — the band construction this method re-points at the source-bias function.
