# STD-RBDO: contents of the repository

This repository contains the MATLAB code for Cases 1 and 2 of the paper *Physics-Informed Spatiotemporal-dependent
Reliability-Based Design Optimization with Trajectory Stepwise Uncertainty Reduction*, submitted to Reliability
Engineering & System Safety.

- Case 1 reproduces **Table 3** and **Table 4** of the paper. Case 2 reproduces **Table 5**.
- The code contains no plotting, and no result files are shipped. Every number is recomputed by the code.
- The code and result files of Case 3 (wind-turbine tower, COMSOL finite element model) are available on request.

---

## 1 Directory structure

```
STD-RBDO/
├── README.md                    overview of the repository
├── CONTENTS.md                  this file
├── LICENSE                      MIT License
├── .gitignore                   excludes results/ and *.asv
│
├── Case1/                       Case 1: analytical benchmark with a migrating failure hotspot
│   ├── main_case1.m             ★ one-click script (MODE = 'test' / 'quick' / 'full' / 'tables')
│   ├── run_case1_table3.m       runs the 7 methods of Table 3 (10 runs each), then prints the table
│   ├── make_case1_table3.m      builds Table 3 from existing runs
│   ├── case1_reference_rows.m   first three rows of Table 3: deterministic optimization,
│   │                            simplified-model-only RBDO, reference solution
│   ├── run_case1_table4.m       runs the 8 configurations of the Table 4 ablation study (10 runs each),
│   │                            then prints the table
│   ├── make_case1_table4.m      builds Table 4 from existing runs
│   ├── README.md                details of Case 1 (expected output, random seeds)
│   ├── .gitignore
│   └── src/                     36 functions c36_*.m (see Section 6)
│
└── Case2/                       Case 2: dynamic RBDO of a variable-depth beam under moving loads
    ├── main_case2.m             ★ one-click script (MODE = 'test' / 'quick' / 'full' / 'tables')
    ├── run_case2_table5.m       runs the 6 configurations of Table 5 (10 runs each), then prints the table
    ├── make_case2_table5.m      builds Table 5 from existing runs (missing validations are computed)
    ├── case2_reference_rows.m   first three rows of Table 5: deterministic optimization,
    │                            simplified-model-only RBDO, reference solution
    ├── README.md                details of Case 2 (expected output, random seeds, validation procedure)
    ├── .gitignore
    └── src/                     42 functions b36_*.m and c2_*.m (see Section 6)
```

Running a case creates a `results/` folder inside the case folder. It holds the run records, caches and output tables.
`.gitignore` excludes it, so it is not uploaded.

---

## 2 Tables of the paper and the scripts that produce them

| Table | Content | One-click script | Step functions | Output file |
|---|---|---|---|---|
| Table 3 | Optimized designs of Case 1 | `Case1/main_case1.m` | `run_case1_table3`, `make_case1_table3`, `case1_reference_rows` | `Case1/results/table3_case1.csv` |
| Table 4 | Ablation study of Case 1 | `Case1/main_case1.m` | `run_case1_table4`, `make_case1_table4` | `Case1/results/table4_case1.csv` |
| Table 5 | Optimized designs at stopping in Case 2 | `Case2/main_case2.m` | `run_case2_table5`, `make_case2_table5`, `case2_reference_rows` | `Case2/results/table5_case2.csv` (table), `Case2/results/table5_case2_runs.csv` (stopped design and validation probability of failure of every run) |

The two cases are independent. Either one can be downloaded and run on its own.

---

## 3 How to run

**Requirements**

| | Case 1 | Case 2 |
|---|---|---|
| Required | MATLAB, Optimization Toolbox, Statistics and Machine Learning Toolbox | MATLAB, Optimization Toolbox |
| Optional | Parallel Computing Toolbox (without it, the runs are executed one after another) | Same |
| Memory | About 2 GB per parallel worker | About 1.5 GB per parallel worker and 2 GB for the MATLAB client |
| Disk | About 0.4 GB (cached benchmark data) | About 10 MB |

The code was developed and tested with MATLAB R2026b on Windows.

**Steps**

1. Open `Case1/main_case1.m` or `Case2/main_case2.m`.
2. Set the two parameters at the top of the file:
   ```matlab
   MODE     = 'full';   % 'test' | 'quick' | 'full' | 'tables'
   NWORKERS = 4;        % number of parallel workers
   ```
3. Run the script (F5). It works from any current folder and adds its own paths.

| `MODE` | What it does |
|---|---|
| `'test'` | Unit tests only, about 1 min. Checks that the environment works |
| `'quick'` | Run 1 of every configuration, then the tables. A trial run |
| `'full'` | 10 runs of every configuration, as in the paper, then the tables |
| `'tables'` | No new runs. Rebuilds the tables from the runs already in `results/main` |

Finished runs are skipped. A script can be interrupted at any time and continues where it stopped when it is run again.

---

## 4 Code configuration of each table row

### Case 1, Table 3

Settings shared by all runs:
- batch size B = 4;
- 12 initial HF simulations;
- a budget of 120 HF simulations;
- 8192 Sobol' working samples, drawn afresh for each run.

| Table row | Code configuration |
|---|---|
| Deterministic optimization | `case1_reference_rows`: nominal inputs and nominal load path, linear programme within the design bounds [0.4, 1.8]² |
| Simplified-model-only RBDO | `case1_reference_rows`: RBDO with the simplified physical model alone (Layer 2) on the validation population |
| Reference solution | `c36_truth('validation')`: one-dimensional search with the exact model on the validation population (Appendix A of the paper) |
| TraSUR (proposed) | `method = 'stmsur_rb'`, three-layer PI-Kriging |
| SL-ABALQ / EPAK / SILK / AK-MCS-U | `method = 'slabalq'` / `'epak'` / `'silk'` / `'ukb'`, three-layer PI-Kriging |
| Kriging + SILK / Kriging + AK-MCS-U | `method = 'silk'` / `'ukb'`, ordinary Kriging |

### Case 1, Table 4 (ablation study, all with TraSUR)

| Configuration | Code configuration |
|---|---|
| 1 | Ordinary Kriging (SE kernel, constant mean) |
| 2 | Layer 2 only: SE kernel and the calibrated simplified-model trend |
| 3 | Layer 1 only: operator constraint, constant mean |
| 4 | Layers 1 and 2, stationary discrepancy |
| 5 | Three layers, nested scheme, inner tolerance 10% (at most 3 batches per design) |
| 6 | Three layers, nested scheme, inner tolerance 5% |
| 7 | Three layers, single loop, learning at the current design only |
| 8 (Proposed) | Three layers, single loop, with the design neighborhood (the proposed setting; same runs as TraSUR in Table 3) |

### Case 2, Table 5

Settings shared by all runs:
- batch size B = 4;
- 16 initial HF simulations;
- a budget of 80 HF simulations;
- the same 4096 working load histories for every method.

| Table row | Code configuration |
|---|---|
| Deterministic optimization | `b36_ddo`: nominal inputs, nominal HF peak stress at most 160 MPa, solved by SQP |
| Simplified-model-only RBDO | `b36_lf_only({'Efield'})`: simplified physical model alone, no HF simulation |
| Reference solution | `b36_truth_ref`: exact (HF) optimum of the 4096 working histories (Appendix A of the paper) |
| TraSUR (proposed) | `method = 'stmsur'`, `spec.lf = 'Efield'`, `spec.nonstat = true` (with the design neighborhood) |
| SL-ABALQ / EPAK / SILK / AK-MCS-U | `method = 'slabalq'` / `'epak'` / `'silk'` / `'ukb'`, same surrogate |
| Layer 1 only (ablation) | `method = 'stmsur'`, `spec.lf = 'qs'`, `spec.nonstat = false` |

The probability of failure in the last column is computed by `b36_validate`. It evaluates the full-grid HF model on an
independent validation population of 16,384 histories.

---

## 5 Computing time

The times below are estimated from the recorded run times of the paper runs. Each run used one single-threaded
worker on an Intel Core i7-12700H laptop.

| Step | Total worker time | Time with 4 parallel workers |
|---|---|---|
| Case 1, Table 3 (70 runs) | About 12 h | About 3 h |
| Case 1, the 70 additional runs of Table 4 | About 16 h | About 4 h |
| Case 1, first call of `MODE = 'tables'` (builds the exact-model scoring data) | — | About 7 min |
| Case 2, Table 5 (60 runs) | About 17 h | About 4–5 h |
| Case 2, reference solution `b36_truth_ref` | — | About 33 min (2 workers) |
| Case 2, reference rows and validation of the 60 stopped designs | — | About 1 h |
| `MODE = 'test'`, either case | — | About 1 min each |

The exact model is used only to score designs and build the tables, never inside learning. In Case 2, a large part
of each run is spent scoring every iterate with the exact model.

---

## 6 Files in `src/`

### Case 1: `Case1/src/` (36 files)

| Files | Content |
|---|---|
| `c36_config.m` | Problem constants (target probability of failure 1%, design bounds, time–space grids) |
| `c36_process.m`, `c36_fields.m`, `c36_loadeffect.m` | Load processes; exact (HF) and simplified load effects |
| `c36_bank.m`, `c36_lhs_states.m` | Sample populations; initial HF simulations (maximin Latin hypercube) |
| `c36_true_envelope.m`, `c36_envelope_field.m`, `c36_margin.m`, `c36_pf.m` | Trajectory envelopes and probability of failure (common manufacturing mode integrated exactly) |
| `c36_solve_rbdo.m` | RBDO subproblem on a frozen envelope (global one-dimensional search) |
| `c36_cgp_fit.m`, `c36_cgp_predict.m`, `c36_gp_fit.m`, `c36_kernel.m`, `c36_ckern.m`, `c36_wamp.m`, `c36_cw.m` | Three-layer PI-Kriging and ordinary Kriging |
| `c36_popeval.m`, `c36_pointfield.m`, `c36_pointsigma.m` | Posterior on the working population |
| `c36_events.m`, `c36_events_rb.m`, `c36_clark_min.m`, `c36_phi2gain.m`, `c36_acq_sur.m` | TraSUR: trajectory events, Clark moment matching, closed-form batch criterion with Schur updates |
| `c36_rvfun.m` | Posterior standard deviation of the working probability of failure (stopping rule) |
| `c36_acq_point.m`, `c36_acq_rv.m` | Compared learning functions (SILK, AK-MCS-U, EPAK, SL-ABALQ, …) |
| `c36_run.m`, `c36_run_matrix.m`, `c36_pool.m` | One single-loop STD-RBDO run; parallel batch of runs; parallel pool (serial without the Parallel Computing Toolbox) |
| `c36_collect.m` | Long table of all runs and iterations, scored on the validation population |
| `c36_truth.m` | Cached benchmark data (validation population and the working population of each run) |
| `c36_test_gp.m`, `c36_test_events.m` | Unit tests of the surrogate and of the event covariances |

### Case 2: `Case2/src/` (42 files)

| Files | Content |
|---|---|
| `b36_config.m`, `c2_config.m` | Problem constants (target probability of failure 1%, allowable stress 160 MPa, design bounds, grids) |
| `c2_load.m`, `c2_section.m`, `c2_hermite.m`, `c2_assemble.m`, `c2_field_operators.m`, `b36_hf.m` | Moving loads; random section and modulus fields; HF model (80 Hermite beam elements, Newmark integration, 1 ms step) |
| `b36_loads.m`, `b36_lf.m`, `b36_cache.m` | Layer 1 (exact quasi-static moment) and Layer 2 (two-mode simplified model); load cache |
| `b36_pop.m`, `b36_screen.m` | Working and validation populations; screening with the simplified model |
| `b36_feat.m`, `b36_block.m`, `b36_gp_fit.m`, `b36_gp_k.m`, `b36_gp_kpair.m`, `b36_gp_pred.m` | Layer 3: discrepancy Gaussian process on the critical windows |
| `b36_popeval.m`, `b36_clark_min.m` | Posterior of the trajectory margins on the working population (Clark moment matching) |
| `b36_solve.m` | RBDO subproblem on the frozen surrogate (SQP) |
| `b36_events.m`, `b36_phi2cov.m`, `b36_phi2gain.m`, `b36_cands.m`, `b36_catcands.m`, `b36_subcands.m`, `b36_acq_sur.m` | TraSUR: trajectory events, posterior standard deviation of the probability of failure, candidates, closed-form batch criterion over the design neighborhood |
| `b36_acq_point.m` | Compared learning functions (SL-ABALQ, EPAK, SILK, AK-MCS-U, …) |
| `b36_run.m`, `b36_run_matrix.m`, `b36_pool.m` | One single-loop STD-RBDO run; parallel batch of runs; parallel pool (serial without the Parallel Computing Toolbox) |
| `b36_collect.m` | Long table of all runs and iterations |
| `b36_truth_at.m`, `b36_truth_ref.m`, `b36_validate.m`, `b36_pf_val.m` | Exact (HF) scoring, reference solution, validation probability of failure |
| `b36_ddo.m`, `b36_lf_only.m` | Deterministic and simplified-model-only designs |
| `b36_test_core.m`, `b36_test_acq.m` | Unit tests of the surrogate, the margins and the acquisition against brute-force posteriors |

---

## 7 Reproducibility checks

All checks were run with MATLAB R2026b on Windows.

- **Tables rebuilt from the archived paper runs** (`MODE = 'tables'`): every number of Tables 3, 4 and 5 equals the
  paper.
- **Runs recomputed from scratch with this code only**: the results are bit-identical to the archived runs.
  - Case 1: run 1 of each of the 7 methods of Table 3.
  - Case 2: the reference solution, the deterministic optimum, the simplified-model-only RBDO, and run 1 of TraSUR
    and SILK.
- **Unit tests** (`MODE = 'test'`): all pass in both cases.

Small floating-point differences between MATLAB releases or platforms can change individual learning steps, so the
final statistics may differ slightly.
