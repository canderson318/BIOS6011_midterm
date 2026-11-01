import pandas as pd
import numpy as np
from numpy import log
from sklearn.datasets import load_iris
import matplotlib.pyplot as plt
import os
import scipy.stats as stats
import statsmodels.api as sm
from scipy.optimize import brentq

# \\\\
# \\\\
# ———— Functions
# \\\\
# \\\\

def soft_impute(X, W,maxit=50, rank=10, lam = 1e-3,tol = 1e-3, floor = 0.0, early_stop = True):
    """Impute masked entries of X with iterative soft-thresholded SVD.
    Fills entries where W is True with column means, then repeatedly takes the SVD of the filled matrix,
    shrinks singular values by lam, keeps the top `rank`, and reconstructs.
    Reconstructions are clipped at `floor` to keep values non-negative. Observed entries are never changed.

    Parameters:
        X: (n, p) data matrix.
        W: (n, p) boolean mask, True = missing/to impute.
        maxit: max iterations. rank: max singular values kept.
        lam: shrinkage.
        tol: stop when the change in delta is below this if early_stop is True.
        floor: lower clip value.
        early_stop: stop when difference between iterations is < tol
    Returns:
        (delta, X_imputed):  where delta is the squared change in imputed entries per iteration.
    """
    def _reconstruct(X, rank, floor, lam):
        u, s, vt = np.linalg.svd(X, full_matrices=False)
        s_ = np.maximum(s - lam, 0)
        return np.maximum(u[:, :rank] @ np.diag(s_[:rank]) @ vt[:rank, :], floor)
    col_means = np.nanmean(np.where(W, np.nan, X), axis=0) # colmeans ignoring where W
    M = np.where(W, col_means, X)  # initial fill
    Z = _reconstruct(M, rank, floor, lam)
    delta = []
    for _ in range(maxit):
        M = np.where(W, Z, X)  # blend: imputed where 0, original elsewhere
        prev = Z
        Z = _reconstruct(M, rank, floor, lam)
        delta.append(np.sum((prev[W] - Z[W]) ** 2))
        if len(delta) > 1 and early_stop:
            if np.abs(delta[-1] - delta[-2]) < tol:
                return delta, np.where(W, Z, X)
    if len(delta) > 1:
        if np.abs(delta[-1] - delta[-2]) >= tol:
            print(f'No convergence to `tol` at {maxit} iterations')
    return delta, np.where(W, Z, X)

def kde(x,r = None,**kwargs):
    if r is None:
        min,max  = x.min(), x.max()
    else:
        min,max = r
    range = np.linspace(min,max,500)
    return range,stats.gaussian_kde(x, **kwargs)(range)

def clr(X):
    L = np.log(X)
    return L - L.mean(1, keepdims = True)

def standardize(X):
    return (X - X.mean(0, keepdims = True)) / X.std(0, keepdims = True)

# \\\
# \\\
# ––– Evaluate
# \\\
# \\\


"""
X ~ Multinomial(N,theta) # n x p counts matrix; rows sum to N_i
theta ~ Dirichlet(alpha) # alpha is mean composition length p; low value; sums to 1
alpha ~ Lognormal(0,1) # how abundant each feature is on average
N ~ LogNormal(log_mu, sigma) # log_mu is log(mu); clip and round to integers; length n; sample library size
"""

# ——— Hyper-Parameters ———

# seed
rng = np.random.default_rng(21233)
# n samples
n = 100
# num features
p = 5000


# library size pars: mean is in log space (median) and sd is in ± fold changes from mean
N_log_mu, N_sigma = log(1e6), log(2)/2


# —— alpha parameter for dirichlet ———
alpha_sigma = 1.5
# mean abundance per feature
m = rng.lognormal(mean = 0, sigma = alpha_sigma, size = p)
m /= m.sum()

## sigma is on log scale so not in unit deviations, in fold deviations so 1 sigma = e^sigma. 2 fold spread around mean is ln(2)/2
alpha_concentration = 2000
def zero_frac(conc, m, N_med):
    return np.mean((1 + N_med / conc) ** (-m * conc))
target = 0.50
alpha_concentration = brentq(lambda c: zero_frac(c, m, 1e6) - target, 1, 1e6)
# alpha controls how prevalent each feature is (relative to its total); and values control how much each feature's between-sample variability
alpha = m * alpha_concentration # zero fraction falls as alpha goes up
# plt.hist(alpha, bins = 30, ); plt.show()

# ——— Dirichlet samples proportions of counts ———
theta = rng.dirichlet(alpha, size = n) # vectors of feature probabilities across samples

# ——— add effects to theta ———
feats = rng.choice(p, p//100, replace = False)
is_effect = np.zeros(p, bool)
is_effect[feats] = True
groupA = np.arange(n) < (n/2)
groupB = ~groupA
# fc = 3
fc = 1
theta[np.ix_(groupA, feats)] *= fc
theta/=theta.sum(1, keepdims=True) # re-simplex

# —— Sample Ns (library sizes) ———
N = rng.lognormal(mean=N_log_mu, sigma = N_sigma, size = n).clip(1e5, 1e7).astype(int)

# ——— sample counts: c ~ count_i * theta ————
X = rng.multinomial(N, theta).astype(float)
xdescr = pd.Series((X == 0).mean(axis = 0)).describe()
print(f"{((X == 0).mean(axis = 0) > .9).sum()} features with > 90% zeros")

# remove cols that are too sparse
keep = (X==0).mean(axis = 0) < .9
X = X[:, keep]
is_effect = is_effect[keep]  # truth labels, aligned with columns of X after filtering
kept_feats = np.where(is_effect)[0]  # column indices (filtered space) of planted features

print(f"P(X=0) = {xdescr['mean']:.2f} ± {xdescr['std']:.2f}")

# ——— relative abundances ———
P = X / X.sum(axis = 1, keepdims = True)

# ——— Imputed zeros ———
W = P==0
floor = P[P>0].min() / 2
# np.linalg.svd(P, compute_uv=False)[:15]
delta, imp = soft_impute(P, W, rank = min(P.shape) // 5, maxit = 100, floor = floor, tol = 1e-7, lam = .01)
plt.figure()
plt.plot(delta);
print(f"iter = {len(delta)}; SSE = {delta[-1]:.2e}; count filled = {(imp[W] <= floor).mean():.2%}")
clr_imp = clr(imp)
clr_imp = standardize(clr_imp)

# ——— Filled zeros ———
# pcount = imp.min() / 2
pcount = 1e-5
clr_pcount = clr(np.where(P==0, pcount, P))
clr_pcount = standardize(clr_pcount)

plt.figure()
plt.plot(*kde(clr_imp[:,kept_feats[0]]))
plt.plot(*kde(clr_pcount[:,kept_feats[0]]))

y = groupA.astype(int)
i = kept_feats[0]

res_pcount = []
for i in range(clr_pcount.shape[1]):
    a,b = clr_pcount[groupA,i],clr_pcount[groupB,i]
    t, pval, *_ = stats.ttest_ind(a,b)
    diff = a.mean() - b.mean()
    res_pcount.append({'i':i, 't':t,'p':pval, 'a-b':diff})

res_pcount = pd.DataFrame(res_pcount).set_index("i")
res_pcount['q'] = stats.false_discovery_control(res_pcount['p'], method = 'BH')
res_pcount['has_effect'] = is_effect
res_pcount['effect_detected'] = res_pcount['p'] < .05
# res_pcount['effect_detected'] = res_pcount['q'] < .05
confusion_pcount = pd.crosstab(res_pcount['effect_detected'], res_pcount['has_effect']).iloc[::-1, ::-1]

res_imp = []
for i in range(clr_imp.shape[1]):
    a,b = clr_imp[groupA,i],clr_imp[groupB,i]
    t, pval, *_ = stats.ttest_ind(a,b)
    diff = a.mean() - b.mean()
    res_imp.append({'i':i, 't':t,'p':pval, 'a-b':diff})

res_imp = pd.DataFrame(res_imp).set_index("i")
res_imp['q'] = stats.false_discovery_control(res_imp['p'], method = 'BH')
res_imp['has_effect'] = is_effect
res_imp['effect_detected'] = res_imp['p'] < .05
# res_imp['effect_detected'] = res_imp['q'] < .05
confusion_imp = pd.crosstab(res_imp['effect_detected'], res_imp['has_effect']).iloc[::-1, ::-1]


print(confusion_imp)
print(confusion_pcount)

def spec_sens(x):
    x = np.array(x)
    return  x[1,1] / x[:,1].sum(), x[0,0] / x[:,0].sum()

def fdr(x):
    x = np.array(x)
    return x[0,1] / x[0,:].sum()

def ppv_npv(x):
    x = np.array(x)
    return x[0,0] / x[0,:].sum(), x[1,1] / x[1,:].sum()

print(f"Imputation specificity/sensitivity {spec_sens(confusion_imp)[0]:.2f}/{spec_sens(confusion_imp)[1]:.2f}")
print(f"Pcount specificity/sensitivity {spec_sens(confusion_pcount)[0]:.2f}/{spec_sens(confusion_pcount)[1]:.2f}")

print(f"Imputation FDR {fdr(confusion_imp):.2f}")
print(f"Pcount FDR {fdr(confusion_pcount):.2f}")

print(f"Imputation PPV/NPV {ppv_npv(confusion_imp)[0]:.2f}/{ppv_npv(confusion_imp)[1]:.2f}")
print(f"Pcount PPV/NPV {ppv_npv(confusion_pcount)[0]:.2f}/{ppv_npv(confusion_pcount)[1]:.2f}")

plt.hist(res_imp['p'], alpha = .5)
plt.hist(res_pcount['p'], alpha = .5)
