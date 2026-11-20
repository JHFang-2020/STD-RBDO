# STD-RBDO, Case 2: dynamic RBDO of a variable-depth beam under moving loads

MATLAB code for Case 2 of the paper

> J. Fang, Y. Shi, W. Hu, Q. Su, H. Yang, I. Lee. *Physics-Informed Spatiotemporal-dependent Reliability-Based Design
> Optimization with Trajectory Stepwise Uncertainty Reduction*. Submitted to Reliability Engineering & System Safety.

It reproduces **Table 5** (optimized designs at stopping in Case 2) of the paper. The code contains no plotting.

## Requirements

- MATLAB with the Optimization Toolbox.
- Parallel Computing Toolbox (optional). Without it, the runs and the HF validations are executed one after another.
- About 1.5 GB of memory per parallel worker and 2 GB for the MATLAB client. Disk space is small (about 10 MB of
  results for the 60 runs).

The code was developed and tested with MATLAB R2026b on Windows. Earlier releases with the same toolboxes are expected to
work but were not tested. Results are deterministic for a given MATLAB release. Small floating-point differences between
releases or platforms can change individual learning steps.

## Quick start

Open `main_case2.m`, set `MODE` and `NWORKERS` at the top and run it (F5). The script works from any current folder.

| `MODE` | What it does |
|---|---|
| `'test'` | Unit tests of the surrogate, the margins and the acquisition (about 1 min) |
| `'quick'` | Run 1 of every configuration of Table 5 (6 runs), the reference rows and the validations, then the table |
| `'full'` | 10 runs of every configuration, as in the paper, then the table |
| `'tables'` | Rebuild Table 5 from the runs already in `results/main` (no new runs) |

The same steps can be called directly:

```matlab
run_case2_table5(1, 2)     % quick check: run 1 of the six configurations, 2 workers
run_case2_table5()         % Table 5: 10 runs of each configuration (6 x 10 runs)
```

The script first runs all missing runs, then computes the three reference rows, validates the stopped designs, and
prints the table. It writes the table to `results/table5_case2.csv` and the stopped design of every run to
`results/table5_case2_runs.csv`. Finished runs and validations are skipped, so the script can be interrupted and
restarted. To rebuild the table from existing runs only, call `make_case2_table5()`. To compute only the first three
rows, call `case2_reference_rows()`. The unit tests take about 30 s: `addpath('src'); b36_test_core(); b36_test_acq();`.

**Computing time** (measured with 2 single-threaded workers on an Intel Core i7-12700H laptop):

| Step | Time |
|---|---|
| Reference solution `b36_truth_ref` (built once, before the first run) | 33 min on 2 workers |
| One learning run to the budget of 80 HF simulations (one worker) | 14 min (TraSUR, run 1), 16 min (SILK, run 1); about 15–25 min on average, longest for "Layer 1 only" |
| Validation of one stopped design (`b36_validate`, about 1000–2100 full-grid HF simulations, one worker) | about 3 min |
| Deterministic optimum `b36_ddo` (including about 8500 full-grid HF simulations for its validation) | 9 min on 2 workers |
| Simplified-model-only RBDO `b36_lf_only` and its validation | about 2 min |

A large part of each run is the scoring of every iterate with the exact model on the working population (fields `qTrue`,
`nFailW`, `costGapW` of the run record). This scoring is never used inside learning. The 60 runs of Table 5 take about
17 h on one worker, i.e. about 4–5 h with 4 workers; the reference rows and the 60 validations add about 1 h with
4 workers.

The first call builds the cached benchmark data:
- `results/truth/truth_ref.mat`: the reference solution (exact HF optimum of the 4096 working histories);
- `results/C2_val_cache.mat`: the validation probability of failure of every validated design;
- `results/C2_ddo.mat` and `results/C2_lf_only_designs.csv`: the deterministic and simplified-model-only designs.

These data are used only to score designs, never inside learning.

## Expected output

`run_case2_table5()` prints the table below. These are the values of Table 5 in the paper.

```
Method                             HF simulations   Design (d1, d2, d3)          Objective            Pf (validation)
Deterministic optimization         —                (0.1600, 0.2796, 0.16000)     0.59962              51.7%
Simplified-model-only RBDO         0                (0.1600, 0.3000, 0.19139)     0.65139              0.677%
Reference solution (exact model)   —                (0.1600, 0.3000, 0.18916)     0.64916              0.934%
TraSUR (proposed)                  38.0 ± 3.4       (0.1601, 0.3000, 0.18919)     0.64933 ± 0.00024    0.928% ± 0.015%
SL-ABALQ                           56.4 ± 8.3       (0.1602, 0.3000, 0.18954)     0.64969 ± 0.00019    0.903% ± 0.017%
EPAK                               57.2 ± 15.6      (0.1603, 0.2999, 0.18923)     0.64946 ± 0.00019    0.926% ± 0.008%
SILK                               38.0 ± 5.4       (0.1601, 0.3000, 0.18933)     0.64935 ± 0.00019    0.925% ± 0.007%
AK-MCS-U                           40.8 ± 3.7       (0.1601, 0.3000, 0.18935)     0.64942 ± 0.00020    0.925% ± 0.010%
Layer 1 only (ablation)            51.2 ± 9.8       (0.1740, 0.3000, 0.19069)     0.66465 ± 0.08269    2.364% ± 0.836%
```

The HF simulations include the 16 initial ones. The objective is C = d1 + d2 + d3 (m); the expected mass is
3768 kg/m × C. Values after ± are standard deviations over the 10 runs. In run 6 of the "Layer 1 only" ablation the
design jumps to the upper bounds (0.30, 0.30, 0.30); this run causes the large standard deviation of C.

## What is computed

| Table row | Configuration in the code |
|---|---|
| Deterministic optimization | `b36_ddo`: nominal inputs (ξ = 0), min C subject to a nominal HF peak stress of at most 160 MPa (full grid), SQP within the design bounds [0.16, 0.30]³ m |
| Simplified-model-only RBDO | `b36_lf_only({'Efield'})`: d1 = 0.16 m and d2 = 0.30 m at their bounds (as at the reference solution); smallest d3 for which the 41st-largest simplified-model peak stress of the 4096 working histories equals 160 MPa; no HF simulation |
| Reference solution | `b36_truth_ref`: exact (HF) optimum of the 4096 working histories on the working grid (Pf ≤ 1% ⇔ 41st-largest HF peak stress < 160 MPa); d1 and d2 at their bounds, checked by perturbing both bounds; d3 by a one-dimensional root search |
| TraSUR (proposed) | `method = 'stmsur'`, `spec.lf = 'Efield'`, `spec.nonstat = true` |
| SL-ABALQ | `method = 'slabalq'`, same surrogate |
| EPAK | `method = 'epak'`, same surrogate |
| SILK | `method = 'silk'`, same surrogate |
| AK-MCS-U | `method = 'ukb'`, same surrogate |
| Layer 1 only (ablation) | `method = 'stmsur'`, `spec.lf = 'qs'`, `spec.nonstat = false` |

The run files are named `<method>_<spec>_B4_s<run>.mat`, for example `stmsur_EfieldNS_B4_s001.mat` (TraSUR, run 1) and
`stmsur_qs_B4_s001.mat` (Layer 1 only, run 1).

Surrogate specifications (`spec`):
- Layer 1, the exact operator: the beam is statically determinate, so the quasi-static bending moment follows from the
  Green's function of the load (`b36_loads`). It needs no FE model and does not depend on the design.
- Layer 2, `lf = 'Efield'`: the simplified physical model adds two assumed sine modes on each history's own section and
  modulus field (mode-acceleration form, `b36_lf`).
- Layer 3, `nonstat = true`: the Gaussian-process discrepancy y = log(W_HF/W_LF) of the 0.1 s window maxima has a
  non-stationary amplitude driven by the dynamic amplification λ = log(W_LF/W_QS) (`b36_gp_fit`).
- `lf = 'qs'`, `nonstat = false`: the exact quasi-static moment alone with a stationary discrepancy (Layer 1 only).

Other method names: `stmsur_1d` (TraSUR at the current design only, without the design neighborhood), `pbalc`, `efflp`
and `random`.

All runs use:
- batch size B = 4;
- 16 initial HF simulations: designs from a Latin hypercube in the design box, histories drawn at random;
- a budget of 80 HF simulations; one HF simulation is one load history at one design and returns the 0.1 s window
  maxima of the stress on the working grid (41 sections, 5 ms, 1–4 s);
- the same 4096 working histories (17 standard normal inputs each) for every method and run.

Every run continues to the budget. The stopping rule of Eq. (22) is applied afterwards to the recorded iterations: the
posterior standard deviation of the working probability of failure must be at most 5% of the target, and the largest
design change at most 10⁻³ m, at two consecutive iterations. Runs that never meet the rule are taken at the budget and
are marked "≥ mean (k/10)". All 60 runs of Table 5 meet the rule.

The probability of failure in the last column is evaluated with the HF model on the full grid (400 points, 1 ms,
0–6 s) on an independent validation population of 16,384 histories (`b36_validate`). Every history whose
simplified-model peak stress lies within a log-margin of 0.15 of min(160 MPa, 164th-largest simplified-model peak) is
simulated with the HF model (at least 204 histories). The others are counted as safe; this is certified when the largest
observed HF/simplified-model error cannot carry them across 160 MPa. The deterministic and simplified-model-only designs,
which are far from the target, use `b36_pf_val`: every history whose simplified-model peak stress lies within a
log-margin of 0.15 of 160 MPa is simulated with the HF model, and the histories below and above this band are counted as
safe and failed, respectively.

**Random seeds.** Run k uses:
- seed 2000 + k for the initial designs and histories;
- seed 1000·k + iteration for the random candidates.

The working population uses seed 1222001 and the validation population seed 1229001 (generator `mt19937ar`).

## Files

| File | Content |
|---|---|
| `main_case2.m` | One-click script (`MODE` = `'test'`, `'quick'`, `'full'` or `'tables'`) |
| `run_case2_table5.m`, `make_case2_table5.m` | Runs and table of Table 5 |
| `case2_reference_rows.m` | Deterministic optimum, simplified-model-only RBDO, reference solution |
| `src/b36_config.m`, `src/c2_config.m` | Problem constants (target probability of failure 1%, allowable stress 160 MPa, design bounds, grids) |
| `src/c2_load.m`, `c2_section.m`, `c2_hermite.m`, `c2_assemble.m`, `c2_field_operators.m`, `b36_hf.m` | Moving loads, random section and modulus fields, HF model (80 Hermite beam elements, Newmark, 1 ms) |
| `src/b36_loads.m`, `b36_lf.m`, `b36_cache.m` | Layer 1 (exact quasi-static moment) and Layer 2 (two-mode simplified model); load cache |
| `src/b36_pop.m`, `b36_screen.m` | Working and validation populations; screening with the simplified model |
| `src/b36_feat.m`, `b36_block.m`, `b36_gp_fit.m`, `b36_gp_k.m`, `b36_gp_kpair.m`, `b36_gp_pred.m` | Layer 3: discrepancy Gaussian process on the critical windows |
| `src/b36_popeval.m`, `b36_clark_min.m` | Posterior of the trajectory margins on the working population (Clark moment matching) |
| `src/b36_solve.m` | RBDO subproblem on the frozen surrogate (SQP) |
| `src/b36_events.m`, `b36_phi2cov.m`, `b36_phi2gain.m`, `b36_cands.m`, `b36_catcands.m`, `b36_subcands.m`, `b36_acq_sur.m` | TraSUR: trajectory events, posterior standard deviation of the probability of failure, candidates, closed-form batch criterion over the design neighborhood |
| `src/b36_acq_point.m` | Compared learning functions (SL-ABALQ, EPAK, SILK, AK-MCS-U, …) |
| `src/b36_run.m`, `b36_run_matrix.m`, `b36_pool.m` | One single-loop STD-RBDO run; parallel batch of runs |
| `src/b36_collect.m` | Long table of all runs and iterations |
| `src/b36_truth_at.m`, `b36_truth_ref.m`, `b36_validate.m`, `b36_pf_val.m` | Exact (HF) scoring, reference solution and validation |
| `src/b36_ddo.m`, `b36_lf_only.m` | Deterministic and simplified-model-only designs |
| `src/b36_test_core.m`, `b36_test_acq.m` | Unit tests of the surrogate, the margins and the acquisition against brute-force posteriors |

## License

This code is released under the MIT License; see [`LICENSE`](../LICENSE) in the repository root.
