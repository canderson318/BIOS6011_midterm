from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq


@dataclass
class Sim:
    X: np.ndarray             # (n, p_retained) counts
    P: np.ndarray          # (n, p_retained) relative abundances
    is_effect: np.ndarray     # (p_retained,) truth labels, aligned with columns of X
    groupA: np.ndarray        # (n,) True for samples with the effect
    N: np.ndarray          # (n,) library sizes
    zero_frac: float       # zero fraction before the sparsity filter
    alpha_concentration: float


def zero_frac(conc, m, N_med):
    """Gamma-Poisson approximation of the mean P(X=0) across features."""
    return np.mean((1 + N_med / conc) ** (-m * conc))


def simulate(n=100, p=5000, alpha_sigma=1.5, zero_target=0.5, alpha_concentration=None, N_median=1e6, N_sigma=np.log(2)/2, N_clip=(1e5, 1e7), fc=3, effect_frac=0.01, max_zero=0.9, seed=21233):
    """Simulate an n x p count matrix with a group effect on a subset of features.
    m ~ Lognormal(0, alpha_sigma), normalized: mean composition.
    alpha = m * alpha_concentration; if alpha_concentration is None it is solved so the expected zero fraction is zero_target.
    theta ~ Dirichlet(alpha) per sample; group A's theta on effect features is multiplied by fc, then renormalized. N ~ Lognormal(log(N_median), N_sigma), clipped to N_clip.
    X ~ Multinomial(N, theta). Features with zero fraction >= max_zero are dropped.
    is_effect is all False when fc == 1 (null simulation).
    """
    rng = np.random.default_rng(seed)

    m = rng.lognormal(mean=0, sigma=alpha_sigma, size=p)
    m /= m.sum()
    if alpha_concentration is None:
        alpha_concentration = brentq(lambda c: zero_frac(c, m, N_median) - zero_target, 1, 1e6)
    alpha = m * alpha_concentration

    theta = rng.dirichlet(alpha, size=n)

    feats = rng.choice(p, int(round(p * effect_frac)), replace=False)
    groupA = np.arange(n) < (n / 2)
    theta[np.ix_(groupA, feats)] *= fc
    theta /= theta.sum(1, keepdims=True)

    N = rng.lognormal(mean=np.log(N_median), sigma=N_sigma, size=n).clip(*N_clip).astype(int)
    X = rng.multinomial(N, theta).astype(float)
    zf = (X == 0).mean()

    is_effect = np.zeros(p, bool)
    if fc != 1:
        is_effect[feats] = True
    keep = (X == 0).mean(axis=0) < max_zero
    X, is_effect = X[:, keep], is_effect[keep]
    P = X / X.sum(axis=1, keepdims=True)

    return Sim(X, P, is_effect, groupA, N, zf, alpha_concentration)
