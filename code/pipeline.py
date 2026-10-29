import itertools
import os

import numpy as np
import pandas as pd
from multiprocess import cpu_count, get_context

from sim import simulate
from impute import soft_impute
from process import clr, standardize, test_features, metrics

METRICS = ['spec', 'sens', 'fdr']


def run(fc=3, lam=.1, pcount=1e-5, rank_frac=.2, maxit=100, tol=1e-3, col='p', alpha=.05, verbose=True, **sim_kw):
    """Simulate, impute zeros two ways, test each feature, and score.
    fc and any other simulate() argument (n, p, zero_target, alpha_sigma, ...) go to the simulation.
    lam, rank_frac, maxit, tol configure soft_impute; pcount is the pseudocount for the comparison.
    col is 'p' or 'q' and alpha the cutoff used to call a feature detected.
    Returns (sim, delta, results, scores).
    """
    s = simulate(fc=fc, **sim_kw)
    W = s.P == 0
    delta, imp = soft_impute(s.P, W, rank=int(min(s.P.shape) * rank_frac), maxit=maxit, tol=tol, lam=lam)
    if verbose:
        print(f"Imputation: iter = {len(delta)}; SSE = {delta[-1]:.2e}; median imputed zero = {np.median(imp[W]):.2e}")
    res = {
        'imp': test_features(standardize(clr(imp)), s.groupA, s.is_effect),
        'pcount': test_features(standardize(clr(s.P.clip(min=pcount))), s.groupA, s.is_effect),
    }
    scores = pd.DataFrame({k: metrics(v, col, alpha) for k, v in res.items()}).T
    return s, delta, res, scores


def pmap(fn, jobs, n_jobs=None):
    """Map fn over jobs in worker processes (n_jobs=1 runs serially).
    Uses `multiprocess` (dill) so it works from Jupyter as well as scripts. Workers are spawned, not forked:
    forking after the parent has run threaded BLAS can deadlock. BLAS threads are pinned to 1 in the
    workers so the SVDs don't oversubscribe the cores.
    """
    n_jobs = n_jobs or max(cpu_count() - 1, 1)
    if n_jobs == 1:
        return [fn(j) for j in jobs]
    names = ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"]
    old = {k: os.environ.get(k) for k in names}
    os.environ.update({k: "1" for k in names})  # inherited by workers started below
    try:
        with get_context("spawn").Pool(n_jobs) as pool:
            return pool.map(fn, jobs, chunksize=1)
    finally:
        for k, v in old.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


def _task(args):
    method, fc, fixed, kw = args
    *_, scores = run(**{'fc': fc, 'col': 'q', 'verbose': False, **fixed, **kw})
    return scores.loc[[method], METRICS].assign(**kw)


def sweep(method, grid, fc=3, seeds=range(5), n_jobs=None, **fixed):
    """Score one method over a grid of arguments and seeds in parallel.
    method: 'imp' or 'pcount'.
    grid:  e.g. {'lam': lams, 'rank_frac': ranks}.
    fixed: other run() arguments held constant (calls default to col='q').
    Returns METRICS averaged over seeds, one row per grid point.
    """
    keys = list(grid)
    jobs = [(method, fc, fixed, {**dict(zip(keys, vals)), 'seed': seed})
            for seed in seeds for vals in itertools.product(*grid.values())]
    out = pd.concat(pmap(_task, jobs, n_jobs))
    return out.groupby(keys, as_index=False)[METRICS].mean()
