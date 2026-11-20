# STD-RBDO, Case 1: analytical benchmark with a migrating failure hotspot

MATLAB code for Case 1 of the paper

> J. Fang, Y. Shi, W. Hu, Q. Su, H. Yang, I. Lee. *Physics-Informed Spatiotemporal-dependent Reliability-Based Design
> Optimization with Trajectory Stepwise Uncertainty Reduction*. Submitted to Reliability Engineering & System Safety.

It reproduces **Table 3** (optimized designs of Case 1) and **Table 4** (ablation study of Case 1) of the paper.
The code contains no plotting.

## Requirements

- MATLAB with the Optimization Toolbox and the Statistics and Machine Learning Toolbox.
- Parallel Computing Toolbox (optional). Without it, the runs are executed one after another.
- About 2 GB of memory per parallel worker and about 0.4 GB of disk space for the cached benchmark data.

The code was developed and tested with MATLAB R2026b on Windows. Earlier releases with the same toolboxes are expected to
work but were not tested. Results are deterministic for a given MATLAB release. Small floating-point differences between
releases or platforms can change individual learning steps.

## Quick start

Open `main_case1.m`, set `MODE` and `NWORKERS` at the top and run it (F5). The script works from any current folder.

| `MODE` | What it does |
|---|---|
| `'test'` | Unit tests of the surrogate and of the event covariances (about 1 min) |
| `'quick'` | Run 1 of every configuration of Tables 3 and 4 (14 runs), then both tables |
| `'full'` | 10 runs of every configuration, as in the paper, then both tables |
| `'tables'` | Rebuild Tables 3 and 4 from the runs already in `results/main` (no new runs) |

The same steps can be called directly:

```matlab
run_case1_table3(1:2, 2)   % quick check: runs 1-2 of the seven methods, 2 workers
run_case1_table3()         % Table 3: 10 runs of each method (7 x 10 runs)
run_case1_table4()         % Table 4: ablation (8 configurations x 10 runs; configuration 8 is shared with Table 3)
```

Each script first runs all missing runs and then prints the table and writes it to `results/table3_case1.csv` or
`results/table4_case1.csv`. Finished runs are skipped, so a script can be interrupted and restarted. To rebuild a table
from existing runs only, call `make_case1_table3()` or `make_case1_table4()`.

**Computing time** (recorded run times of the paper runs, one single-threaded worker per run):
- Table 3: one TraSUR run takes about 15 min, a run of another learning function on PI-Kriging about 6 min, and a
  run on ordinary Kriging about 17 min. The 70 runs take about 12 h of worker time, i.e. about 3 h with 4 workers.
- Table 4 adds 70 runs (configuration 8 is shared with Table 3), about 8–9 min each for configurations 3–7 and about
  25–28 min for configurations 1–2. They take about 16 h of worker time, i.e. about 4 h with 4 workers.

The first call also builds the cached benchmark data in `results/truth/`:
- the independent validation population (2^17 RQMC samples, exact model on a 401 × 201 time–space grid) with the
  reference optimum;
- for each run, the exact model on its working population.

These data are used only to score designs, never inside learning.

## Expected output

`run_case1_table3()` prints the table below. These are the values of Table 3 in the paper.

```
Method                             HF simulations     Design (d1, d2)      Objective            Pf (validation)
Deterministic optimization         —                  (0.4000, 0.4000)       0.8000               95.9%
Simplified-model-only RBDO         0                  (0.5332, 0.9399)       1.4732               4.02%
Reference solution (exact model)   —                  (0.5439, 1.1137)       1.6576               1.000%
TraSUR (proposed)                  69.2 ± 12.4        (0.5427, 1.1137)       1.6563 ± 0.0032      1.010% ± 0.025%
SL-ABALQ                           79.6 ± 7.9         (0.5452, 1.1131)       1.6583 ± 0.0030      0.996% ± 0.022%
EPAK                               ≥116.4 (5/10)      (0.5429, 1.1160)       1.6589 ± 0.0048      0.992% ± 0.036%
SILK                               79.2 ± 11.6        (0.5415, 1.1144)       1.6559 ± 0.0035      1.013% ± 0.026%
AK-MCS-U                           90.8 ± 11.9        (0.5440, 1.1132)       1.6572 ± 0.0025      1.004% ± 0.019%
Kriging + SILK                     ≥109.2 (8/10)      (0.5402, 1.1094)       1.6496 ± 0.0044      1.063% ± 0.035%
Kriging + AK-MCS-U                 ≥111.2 (4/10)      (0.5343, 1.1117)       1.6460 ± 0.0113      1.122% ± 0.184%
```

`run_case1_table4()` prints the following:

```
Config.        N_M              N_H              eps_P(n=36)    eps_P(n=120)
1              88.4             ≥119.6 (1/10)    115.9%         1.56%
2              61.6             80.4             17.9%          0.26%
3              60.4             74.8             25.8%          0.20%
4              59.2             82.8             19.3%          0.30%
5              58.8             76.4             21.0%          0.24%
6              61.2             78.8             105.4%         0.22%
7              54.8             78.8             13.0%          0.19%
8 (Proposed)   54.0             74.0             18.3%          0.14%
```

## What is computed

| Table row | Configuration in the code |
|---|---|
| Deterministic optimization | `case1_reference_rows`: nominal inputs and nominal load path, linear programme within the design bounds [0.4, 1.8]² |
| Simplified-model-only RBDO | `case1_reference_rows`: RBDO with the simplified physical model (Layer 2 only) on the validation population |
| Reference solution | `c36_truth('validation')`: exact model on the validation population, one-dimensional search (Appendix A of the paper) |
| TraSUR (proposed) | `method = 'stmsur_rb'`, three-layer PI-Kriging |
| SL-ABALQ | `method = 'slabalq'`, three-layer PI-Kriging |
| EPAK | `method = 'epak'`, three-layer PI-Kriging |
| SILK | `method = 'silk'`, three-layer PI-Kriging |
| AK-MCS-U | `method = 'ukb'`, three-layer PI-Kriging |
| Kriging + SILK, Kriging + AK-MCS-U | `method = 'silk'` / `'ukb'` on ordinary Kriging |

Surrogate specifications (`spec`):
- `spatial = 'op'`: Layer 1, the exact differential operator. The cubic null-space basis has 4 components.
- `mean = 'lf1'`: Layer 2, the simplified-model mean.
- `nonstat = true`: Layer 3, the non-stationary discrepancy.
- `spatial = 'se'`, `mean = 'const'`: ordinary Kriging.

Other method names:
- `stmsur_rb1d`: TraSUR at the current design only, without the design neighborhood.
- `stmsur`: TraSUR with binary instead of Rao–Blackwellized events.

All runs use:
- batch size B = 4;
- 12 initial HF simulations from a maximin Latin hypercube;
- a budget of 120 HF simulations;
- 8192 scrambled Sobol' working samples, drawn afresh for each run.

Every run continues to the budget. The stopping rule of Eq. (22) is applied afterwards to the recorded iterations: the
posterior standard deviation of the working probability of failure must be at most 5% of the target, and the design
change at most 10⁻³, at two consecutive iterations. Runs that never meet the rule are taken at the budget and are marked
"≥ mean (k/10)".

**Random seeds.** Run k uses:
- seed k for its working population;
- seed 1000 + k for the initial design;
- seed 1000·k + iteration for the random candidates.

The validation population uses seed 20261020.

## Files

| File | Content |
|---|---|
| `main_case1.m` | One-click script (`MODE` = `'test'`, `'quick'`, `'full'` or `'tables'`) |
| `run_case1_table3.m`, `make_case1_table3.m` | Runs and table of Table 3 |
| `run_case1_table4.m`, `make_case1_table4.m` | Runs and table of Table 4 |
| `case1_reference_rows.m` | Deterministic optimum, simplified-model-only RBDO, reference solution |
| `src/c36_config.m` | Problem constants (target probability of failure 1%, design bounds, grids) |
| `src/c36_process.m`, `c36_fields.m`, `c36_loadeffect.m` | Load processes, exact (HF) and simplified load effect |
| `src/c36_bank.m`, `c36_lhs_states.m` | Sample populations and initial HF simulations |
| `src/c36_true_envelope.m`, `c36_envelope_field.m`, `c36_margin.m`, `c36_pf.m` | Trajectory envelopes and probability of failure (common manufacturing mode integrated exactly) |
| `src/c36_solve_rbdo.m` | RBDO subproblem on a frozen envelope (global one-dimensional search) |
| `src/c36_cgp_fit.m`, `c36_cgp_predict.m`, `c36_gp_fit.m`, `c36_kernel.m`, `c36_ckern.m`, `c36_wamp.m`, `c36_cw.m` | Three-layer PI-Kriging and ordinary Kriging |
| `src/c36_popeval.m`, `c36_pointfield.m`, `c36_pointsigma.m` | Posterior on the working population |
| `src/c36_events.m`, `c36_events_rb.m`, `c36_clark_min.m`, `c36_phi2gain.m`, `c36_acq_sur.m` | TraSUR: trajectory events, Clark moment matching, closed-form batch criterion with Schur updates |
| `src/c36_rvfun.m` | Posterior standard deviation of the working probability of failure (stopping rule) |
| `src/c36_acq_point.m`, `c36_acq_rv.m` | Compared learning functions (SILK, AK-MCS-U, EPAK, SL-ABALQ, …) |
| `src/c36_run.m`, `c36_run_matrix.m`, `c36_pool.m` | One single-loop STD-RBDO run; parallel batch of runs |
| `src/c36_collect.m` | Long table of all runs and iterations, with validation scoring |
| `src/c36_truth.m` | Cached benchmark data |
| `src/c36_test_gp.m`, `c36_test_events.m` | Unit tests of the surrogate and of the event covariances |

## License

This code is released under the MIT License; see [`LICENSE`](../LICENSE) in the repository root.
