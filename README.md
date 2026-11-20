# STD-RBDO

MATLAB code for Cases 1 and 2 of the paper

> J. Fang, Y. Shi, W. Hu, Q. Su, H. Yang, I. Lee. *Physics-Informed Spatiotemporal-dependent Reliability-Based Design
> Optimization with Trajectory Stepwise Uncertainty Reduction*. Submitted to Reliability Engineering & System Safety.

The code contains no plotting, and no result files are shipped; every number is recomputed by the code. The code and
result files of Case 3 (wind-turbine tower, COMSOL model) are available on request.

## Contents

| Folder | Case | Reproduces | One-click script |
|---|---|---|---|
| [`Case1/`](Case1/README.md) | Analytical benchmark with a migrating failure hotspot | Table 3 (optimized designs), Table 4 (ablation study) | `Case1/main_case1.m` |
| [`Case2/`](Case2/README.md) | Dynamic RBDO of a variable-depth beam under moving loads | Table 5 (optimized designs at stopping) | `Case2/main_case2.m` |

A file-by-file description, the code configuration of every table row and the computing times are given in
[`CONTENTS.md`](CONTENTS.md).

## Quick start

Open `Case1/main_case1.m` or `Case2/main_case2.m`, set the two parameters at the top and run the script (F5):

```matlab
MODE     = 'full';   % 'test' | 'quick' | 'full' | 'tables'
NWORKERS = 4;        % number of parallel workers
```

| `MODE` | What it does |
|---|---|
| `'test'` | Unit tests only (about 1 min) |
| `'quick'` | Run 1 of every configuration, then the tables |
| `'full'` | 10 runs of every configuration, as in the paper, then the tables |
| `'tables'` | Rebuild the tables from the runs already in `results/main` |

Finished runs are skipped, so a script can be interrupted and restarted. Results are written to `Case1/results/` and
`Case2/results/`.

## Requirements

- MATLAB with the Optimization Toolbox. Case 1 also needs the Statistics and Machine Learning Toolbox.
- Parallel Computing Toolbox (optional). Without it, the runs are executed one after another.

The code was developed and tested with MATLAB R2026b on Windows.

**Computing time** with 4 workers:
- Case 1: about 3 h for Table 3 and about 4 h more for Table 4.
- Case 2: about 5–6 h for Table 5, including the reference solution and the validations.

See the README of each case for details.

## License

This code is released under the MIT License; see [`LICENSE`](LICENSE).
