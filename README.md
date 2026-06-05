# HTML6 — Frequentist Network Meta-Analysis (offline)

A single-page, fully offline tool for **frequentist network meta-analysis (NMA)** using
the graph-theoretical (Rücker 2012) approach, with P-scores/SUCRA ranking, league tables,
node-splitting and design-by-treatment consistency checks, Q-profile τ² confidence
intervals, a contribution-matrix approximation, a comparison-adjusted funnel plot, CINeMA-style
confidence ratings, and threshold analysis. It runs entirely in the browser with no network
access — Plotly is vendored locally.

## Files

| File | Purpose |
|------|---------|
| `index.html` | The app UI + inline application/DOM/Plotly controller. |
| `engine.js`  | Pure NMA / linear-algebra engine (matrix ops, GLS estimator, τ² estimators, P-scores/SUCRA, consistency, contribution matrix, etc.) plus the `Validation` self-test module. Single source of truth — the HTML loads this; there are no duplicate inline copies. |
| `tests.js`   | Pure Node test harness (`node tests.js`) with hand-derived expected values. |
| `plotly.min.js` | Vendored Plotly 2.24.1 (offline). |

## Run

Open `index.html` in any modern browser (no server needed), or visit the GitHub Pages
deployment. Everything is offline.

## Tests

```
node tests.js
```

Pure Node, no dependencies. Prints per-check PASS/FAIL and a final `N passed, M failed`
line; exits non-zero on any failure. Expected values are hand-derived (triangle network
contrast recovery, inverse-variance pooling, Q statistic, P-score sum = K/2, Moore-Penrose
property `A·A⁺·A = A`, OR log-scale back-transform).

## Method notes / scope

- **Engine model:** contrast-level input. Each row is a contrast `(t1, t2, effect, SE)`; the
  Laplacian `L = Bᵀ W B` (with `W = diag(1/(SE² + τ²))`) is assembled and the treatment
  effects are `θ = L⁺ Bᵀ W y` (Moore-Penrose pseudoinverse via Jacobi SVD).
- For ratio measures (OR/RR) enter **log-transformed** effects; the app pools on the log
  scale and back-transforms with `exp()` for display.
- **τ²:** REML (Fisher scoring, default), DerSimonian-Laplace (Jackson 2014), Paule-Mandel.
  Q-profile CI for τ² (Viechtbauer 2007).
- **Ranking:** P-scores (Rücker & Schwarzer 2015) and simulation-based SUCRA.
- **Consistency:** design-by-treatment Q decomposition and node-splitting.
- **Multi-arm trials (limitation):** this engine treats each contrast row independently with a
  diagonal weight matrix; it does **not** construct the within-study contrast covariance that a
  full arm-based / shared-control multi-arm correction (off-diagonal covariance `τ²/2`) would
  require. For networks with multi-arm trials, validate against R `netmeta`, which builds the
  full multi-arm covariance.

## Fixes applied during revival (2026-06-05)

- **Offline:** vendored Plotly 2.24.1 locally as `plotly.min.js` (was a `cdn.plot.ly`
  `<script src>`); removed the Google Fonts `<link>`. No remaining external `src`/`href`/`@import`.
- **Engine extracted** verbatim from the inline `<script>` into `engine.js` as the single
  source of truth; the HTML now loads it via `<script src="engine.js">` and the duplicated
  inline copy was deleted.
- **Correctness bug fixed (SVD pseudoinverse):** the one-sided Jacobi SVD used
  `Math.sign(zeta)`, and `Math.sign(0) === 0`. When two columns had equal norm (`zeta == 0`,
  e.g. the rank-deficient `[[1,-1],[-1,1]]` Laplacian of a 2-treatment network, or any balanced
  network with repeated singular values) the rotation angle collapsed to `t = 0`, the columns
  were never orthogonalised, and the pseudoinverse came out **2× too large**. This produced
  doubled pairwise contrasts and inflated variances on 2-treatment and symmetric networks
  (`A·A⁺·A` returned `2A` instead of `A`). Fixed by treating `zeta == 0` as `+1`. Verified:
  `pinv([[1,-1],[-1,1]]) = [[0.25,-0.25],[-0.25,0.25]]`, `A·A⁺·A = A`, and the netmeta-validated
  4-treatment example (τ²≈0, Q≈0.038, P-score sum≈2.0) is unchanged.
- Added `tests.js`, `.nojekyll`, `.gitignore`, `README.md`, `E156-PROTOCOL.md`.

Note: the in-app header retains its original "NMA Ultimate" branding; "Ultimate" is legacy UI
copy, not a capability claim — see the scope/limitation notes above for what the tool actually does.

## License

Apache License 2.0 — see `LICENSE`.
