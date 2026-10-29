# BIOS6011 Midterm Project

## Motivation

My project uses simulation to compare two different methods for filling left-censored zeros in a compositional dataset. 
Compositional data exists on the simplex meaning values exist between 0 and 1, inclusive, and all values for an observation must sum to 1. 
As such, when one variable goes up, another goes down, e.g. for the following vector if the first value goes down, the remainder increases \[0.8,0.1,0.1\] -> \[0.6, 0.2, 0.2\]. 
This compositional constraint leads to spurios negative correlations which we want to avoid. 
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
Soft imputation is unbiased and uses the structure of the data, but it has the potential to introduce arbitrary values which on a whole may reconstruct well, but for individual features may introduce noise and artificial structure. 
It also uses more parameters which should be well tuned. 
Pseudo-counts are simple, but also introduce artificial variance, although unbiasedly. 


_I hypothesize that soft imputation will better recapitulate the original multi-variate structure than pseudo-count._

## Approach

To test which method is better, I simulate sparse compositional data and introduce an effect. 
I choose to simulate with zeros built in instead of adding them artificially as this is more accurate of biological data.
I then evaluate FDR, sensitivity, and specificity to determine if the method has introduced some unwanted structure.
The method that has higher sensitivity (power) at an FDR < .05 is superior. 

Both methods use lots of parameters, and the simulation step also uses a lot, so i use a grid search to identify the best ones for each method.

# Code

Python 3.12 managed with [uv](https://docs.astral.sh/uv/). 
Dependencies are in `uv.lock`.

## Setup

```bash
git clone <repo-url> && cd midterm
# create .venv 
uv sync --locked      
```

## Layout

| File | Contents |
|---|---|
| `code/sim.py` | `simulate()`: Dirichlet-multinomial counts with a group effect (`fc`) on a subset of features |
| `code/impute.py` | `soft_impute()`: iterative soft-thresholded SVD imputation of zeros |
| `code/process.py` | `clr`, `standardize`, per-feature t-tests (`test_features`), and sensitivity/specificity/FDR (`metrics`) |
| `code/pipeline.py` | `run()` for one simulation, and `sweep()` to score a grid of settings over seeds in parallel |
| `code/main.py` | Tuning sweeps for the pseudocount and for `lam` and `rank_frac`, with plots |
| `code/compare_lrsvd.py` | Compares `soft_impute` to R's `zCompositions::lrSVD` through rpy2 (slow) |

## Running

One simulation, scored for both methods:

```bash
uv run python -c "
import sys; sys.path.insert(0, 'code')
from pipeline import run
*_, scores = run(fc=3, col='q')
print(scores)
"
```

`code/main.py` is written to be run interactively (Jupyter or an editor REPL), where the plots display:

```bash
uv run jupyter lab
```

Sweeps run in parallel with `multiprocess` (spawned workers, one BLAS thread each). To run `main.py` as a plain script, put its sweeps under `if __name__ == "__main__":`, otherwise the workers re-import and re-run it.

Change the settings through arguments, for example `run(fc=2, zero_target=.7, lam=.05)` or `sweep('imp', {'lam': lams, 'rank_frac': ranks}, fc=2, seeds=range(20))`.
