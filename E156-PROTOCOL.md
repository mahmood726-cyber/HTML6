# E156 Protocol — HTML6

- **Project:** HTML6 — offline frequentist Network Meta-Analysis tool
- **Revived:** 2026-06-05
- **Type:** Single-file HTML app + extracted pure JS engine + Node test harness
- **Dashboard:** https://mahmood726-cyber.github.io/HTML6/

## What changed

- Vendored Plotly 2.24.1 locally (`plotly.min.js`); removed Google Fonts link; zero external `src`/`href`/`@import` references — fully offline.
- Extracted the pure NMA / linear-algebra engine (matrix ops, GLS estimator, τ² estimators, P-scores/SUCRA, consistency, contribution matrix, plus the `Validation` module) verbatim into `engine.js` as the single source of truth; deleted the duplicated inline copy and wired `<script src="engine.js">`.
- Fixed a real SVD-pseudoinverse correctness bug: `Math.sign(0) === 0` left equal-norm columns un-rotated, yielding a 2×-too-large Moore-Penrose inverse on rank-deficient / balanced Laplacians (doubled contrasts and variances on 2-treatment and symmetric networks). Verified `A·A⁺·A = A` and that the netmeta-validated example is unchanged.
- Added `tests.js` (50 checks, hand-derived expectations), `.nojekyll`, `.gitignore`, `README.md`.
- Renamed `NMAHTML` → `index.html`.

## Body (E156 draft — CURRENT BODY)

When a clinician faces several competing treatments compared only indirectly across separate trials, which option ranks best and how trustworthy is that ranking? HTML6 ingests contrast-level network data — each row a comparison with its log-effect and standard error — for any connected set of treatments. It fits a frequentist graph-theoretical network meta-analysis, assembling the weighted Laplacian and solving treatment effects through its Moore-Penrose pseudoinverse, with REML, DerSimonian-Laplace, or Paule-Mandel heterogeneity. On a hand-derived triangle and inverse-variance pooling benchmarks it reproduces pooled contrasts, the Q statistic, and P-score sums exactly, matching R's netmeta on the bundled example. During revival a one-sided-Jacobi SVD defect that doubled the pseudoinverse on balanced networks was identified and corrected, restoring the defining `A·A⁺·A = A` identity. The tool runs entirely offline in one browser file, with a fifty-check Node suite guarding the engine. It does not model within-study multi-arm covariance, so multi-arm networks still warrant netmeta confirmation. SUBMITTED: [ ]
