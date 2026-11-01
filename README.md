# BIOS6011 Midterm Project

## Motivation

My project uses simulation to compare two different methods for filling left-censored zeros in a compositional dataset. 
Compositional data exists on the simplex meaning values exist between 0 and 1, inclusive, and all values for an observation must sum to 1. 
As such, when one variable goes up, another goes down. 
The compositional constraint leads to spurious negative correlations which we want to avoid. 
One method of avoiding this is by transforming the simplex, AKA Aitchison space, to Euclidean space, which works well with most statistical tests and tools, like the t-test and PCA.

One issue in most biological compositional data is that values are sparse, most features, e.g. bacterial taxa, are so rarely observed, the majority of their values are zero. 
This means they either weren't there, or that the technical process was not sensitive enough to record them--this is left-censoring.
To convert the simplex into euclidean space, the most common method is to use the centered log-ratio (CLR) transformation, which essentially centers the log abundances on their geometric mean: log(abundance) - mean(log(abundance)).
One issue here though is that we can't calculate log(0) and the data is sparse so that would be a lot of Nan values. 

Two methods exist to deal with these zeros. 
First, a small value, pseudo-count, simply replaces the zeros, a value of min(abundance\[abundance!=0]\) / 2 is often used, but the choice may have large downstream effects. 
Second, the zeros can be imputed using soft imputation which infers values the zeros can take such that the overall variance of the matrix is not altered.
This method iteratively calculates the singular value decomposition of all the data, and fills the zero cells each time. 
When the reconstruction starts to look essentially the same each iteration, it has converged.
Additionally, the singular values are regularized so that only the most informative singular vectors impute the zero cells.

Performing statistical tests on the transformed matrix is the ultimate goal, whether one method is truer to the original data structure is unknown to me.
Soft imputation uses the structure of the data, but it has the potential to introduce arbitrary values which on a whole may reconstruct well, but for individual features may introduce noise and artificial structure. 
It also uses more parameters which should be well tuned. 
Pseudo-counts are simple, but also introduce artificial variance. 


_I hypothesize that soft imputation will better recapitulate the original multi-variate structure than pseudo-count._

## Approach

To test which method is better, I simulate sparse compositional data and introduce an effect. 
I choose to simulate with zeros built in instead of adding them artificially as this is more accurate of biological data.
I then evaluate FDR, sensitivity, and specificity to determine if the method has introduced some unwanted structure.
The method that has higher sensitivity (power) at an FDR < .05 is superior. 

Both methods use lots of parameters, and the simulation step also uses a lot, so I use a grid search to identify the best ones for each method.
The grid searches run at a fold change of 5 with 1% of features differential, averaged over five seeds (120, 100, 113, 1111, 245): 30 pseudocounts from 1e-8 to 1e-2, and 10 values each of `lam` and `rank_frac` from 0.01 to 0.9 (log-spaced).
The best pseudocount is the one with the highest sensitivity at FDR < .05; the best soft imputation setting is the one with the lowest FDR, then the highest sensitivity.
The two methods are then compared at their best settings over 100 fold changes from 0 to 10, one simulated dataset per fold change.
The pseudocount replaces only the zeros of the relative abundances; soft imputation runs in log space with imputed values capped below each feature's smallest observed value.

# Code

Python 3.12 managed with [uv](https://docs.astral.sh/uv/).
Dependencies are in `pyproject.toml` and pinned in `uv.lock`.

## Setup

```bash
git clone https://github.com/canderson318/BIOS6011_midterm.git && cd BIOS6011_midterm
uv sync --locked      # creates .venv
```

`code/main.py` sets the working directory with `wd = Path.home() / 'Documents/school/BIOS6011/assignments/midterm'`.
Change `wd` to the path of your clone (the repository root) before running, since outputs are written to `out/` relative to it.

## Layout

| File | Contents |
|---|---|
| `code/sim.py` | `simulate()`: Dirichlet-multinomial counts with a group effect (`fc`) on a subset of features (`effect_frac`) |
| `code/impute.py` | `soft_impute()`: iterative soft-thresholded SVD imputation of zeros in log space |
| `code/process.py` | `clr`, `standardize`, per-feature t-tests (`test_features`), and sensitivity/specificity/FDR (`metrics`) |
| `code/pipeline.py` | `run()` for one simulation, and `sweep()` to score a grid of settings over seeds in parallel |
| `code/main.py` | Tuning sweeps for the pseudocount and for `lam` and `rank_frac`, the sweep over fold changes, and the plots and tables |
<!--| `code/other/compare_lrsvd.py` | Not used in the report. Compares `soft_impute` to R's `zCompositions::lrSVD` through rpy2 (needs R and `zCompositions`; slow) |-->
| `report/main.typ` | Typst source of the report; reads the results from `out/` |

## Running

Reproduce the report's results (from the repository root, after setting `wd`):

```bash
uv run python code/main.py
```

This runs the grid sweeps in parallel with `multiprocess` and writes to `out/`:
Compile the report with `typst compile --root . report/main.typ`.

One simulation scored for both methods:

```bash
uv run python -c "
import sys; sys.path.insert(0, 'code')
from pipeline import run
*_, scores = run(fc=3, col='q')
print(scores)
"
```

Change the settings through arguments, for example `run(fc=2, zero_target=.7, lam=.05)` or `sweep('imp', {'lam': lams, 'rank_frac': ranks}, fc=2, seeds=range(20))`.

`main.py` keeps its work under `if __name__ == "__main__":`; scripts that call `sweep()` need the same guard, otherwise the spawned workers re-import and re-run the script.
