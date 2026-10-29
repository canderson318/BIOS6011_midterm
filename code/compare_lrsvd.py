"""Compare soft_impute to zCompositions::lrSVD (called from R via rpy2) on the simulated data.
lrSVD is called as in the metagenome pipeline's utils/processing.py (label=0, z_delete=True, z_warning=.9),
without importing it: that module needs scanpy/anndata.
"""
import time
import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

from sim import simulate
from impute import soft_impute
from process import clr, standardize, test_features, metrics

zCompositions = importr('zCompositions')


def lrsvd(P, ncp, z_warning=.9):
    """Impute the zeros of P (samples x features) with zCompositions::lrSVD; returns an array like P.

    lrSVD drops rows/columns with a zero fraction above z_warning. The sparsity filter in simulate()
    keeps columns below 0.9, so nothing should be dropped, and a changed shape raises an error.
    """
    with localconverter(ro.default_converter + numpy2ri.converter + pandas2ri.converter):
        out = zCompositions.lrSVD(ro.conversion.py2rpy(P), label=0, ncp=ncp, z_delete=True, z_warning=z_warning)
        out = np.asarray(ro.conversion.rpy2py(out))
    if out.shape != P.shape:
        raise ValueError(f"lrSVD returned {out.shape} for input {P.shape}; rows/columns were dropped")
    return out


def compare(fc=3, ncps=(2, 5, 10, 20), lam=.1, rank_frac=.2, pcount=1e-5, col='q', alpha=.05, **sim_kw):
    """Impute one simulated dataset with soft_impute, lrSVD at each ncp, and a pseudocount; test and score each.

    fc and any other simulate() argument go to the simulation; lam and rank_frac configure soft_impute;
    col is 'p' or 'q' and alpha the detection cutoff. Returns (sim, imputed, scores).
    """
    s = simulate(fc=fc, **sim_kw)
    W = s.P == 0

    imputed, secs = {}, {}
    def timed(name, f):
        t = time.time()
        imputed[name] = f()
        secs[name] = time.time() - t

    timed('pcount', lambda: np.where(W, pcount, s.P))
    timed('soft_impute', lambda: soft_impute(s.P, W, rank=int(min(s.P.shape) * rank_frac), maxit=100,
                                             tol=1e-3, lam=lam)[1])
    for k in ncps:
        timed(f'lrSVD ncp={k}', lambda k=k: lrsvd(s.P, k))

    res = {name: test_features(standardize(clr(Z)), s.groupA, s.is_effect) for name, Z in imputed.items()}
    scores = pd.DataFrame({name: metrics(r, col, alpha) for name, r in res.items()}).T
    scores['secs'] = pd.Series(secs)
    return s, imputed, scores


def agreement(s, imputed, ref='soft_impute'):
    """Correlation of log imputed values at the zero entries between `ref` and every other method."""
    W = s.P == 0
    lref = np.log(imputed[ref][W])
    return pd.Series({k: np.corrcoef(lref, np.log(Z[W]))[0, 1] for k, Z in imputed.items() if k != ref})


for fc in (1, 3):
    s, imputed, scores = compare(fc=fc)
    print(f"fc = {fc}")
    print(scores.round(3), "\n")
    print("log-imputed correlation with soft_impute at zeros:")
    print(agreement(s, imputed).round(3), "\n")
