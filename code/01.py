import pandas as pd
import numpy as np
from sklearn.datasets import load_iris
import matplotlib.pyplot as plt
import os
from sklearn.model_selection import RepeatedStratifiedKFold
import scipy.stats as stats

def soft_impute(X, W,maxit=50, rank=10):
      def _reconstruct(X, rank):
          u, s, vt = np.linalg.svd(X, full_matrices=False)
          return u[:, :rank] @ np.diag(s[:rank]) @ vt[:rank, :]
      # condition, if, else
      col_means = np.nanmean(np.where(W, np.nan, X), axis=0) # colmeans ignoring where W
      M = np.where(W, col_means, X)  # initial fill
      reconstructed = _reconstruct(M, rank)
      delta = []
      for _ in range(maxit):
          M = np.where(W, reconstructed, X)  # blend: imputed where 0, original elsewhere
          prev = reconstructed
          reconstructed = _reconstruct(M, rank)
          delta.append(np.sum((prev[W] - reconstructed[W]) ** 2))
      return delta, np.where(W, reconstructed, X)

# \\\
# \\\
# ––– Evaluate
# \\\
# \\\

def kde(x,r = None,**kwargs):
    if r is None:
        min,max  = x.min(), x.max()
    else:
        min,max = r
    range = np.linspace(min,max,500)
    return range,stats.gaussian_kde(x, **kwargs)(range)





S = 1000
rng = np.random.default_rng(123)
test_inds = (tuple(rng.choice(len(X), size=S)), tuple(rng.choice(X.shape[1], size=S)))

X_test = X.copy()
X_test[test_inds] =  np.nan

X.shape
W = np.isnan(X_test)
delta , X_res = soft_impute(X_test, W, rank = 3,maxit = 20)
plt.figure()
plt.plot(delta, 'o-')

plt.figure()
plt.plot(X[W], X_res[W], 'o')
np.sqrt(np.mean((X[W] - X_res[W])**2))


def log_likelihood(y, X, mod):
	p = mod.predict_proba(X)
	return np.sum(np.log(p[np.arange(len(y)), y]))

cv_ll = RepeatedStratifiedKFold(n_splits=4, n_repeats=15, random_state=123)
LL = {"RCLR": [], "CLR": []}

for train_idx, test_idx in cv_ll.split(RCLR_PCs, y):
    y_train, y_test = y[train_idx], y[test_idx]
    for k, X in zip(["RCLR", "CLR"], [RCLR_PCs, CLR_PCs]):
        mod = LogReg().fit(X.iloc[train_idx], y_train)
        LL[k].append(log_likelihood(y_test, X.iloc[test_idx], mod))

LL = pd.DataFrame(LL)

plt.figure(figsize = (4,5))
plt.boxplot(LL.loc[:,["CLR", "RCLR"]], tick_labels  = ["CLR", "RCLR"], showfliers = False)
plt.scatter(np.ones(len(LL)), LL["CLR"], alpha =.3, color = 'grey')
plt.scatter(np.ones(len(LL)) + 1, LL["RCLR"], alpha =.3, color = 'grey')
plt.title("K-fold cross-validation")
plt.ylabel("Log-likelihood")
plt.xlabel("Transformation Method")
plt.savefig(results_dir / f"LL_KF_CV_RCLR_CLR.pdf")



def soft_impute(X, W,maxit=50, rank=10):
    def _reconstruct(X, rank):
        u, s, vt = np.linalg.svd(X, full_matrices=False)
        return u[:, :rank] @ np.diag(s[:rank]) @ vt[:rank, :]
    # condition, if, else
    col_means = np.nanmean(np.where(W, np.nan, X), axis=0) # colmeans ignoring where W
    M = np.where(W, col_means, X)  # initial fill
    reconstructed = _reconstruct(M, rank)
    delta = []
    for _ in range(maxit):
    M = np.where(W, reconstructed, X)  # blend: imputed where 0, original elsewhere
    prev = reconstructed
    reconstructed = _reconstruct(M, rank)
    delta.append(np.sum((prev[W] - reconstructed[W]) ** 2))
    return delta, np.where(W, reconstructed, X)

# \\\
# \\\
# ––– Evaluate
# \\\
# \\\

X_raw, y, *_ = load_iris(as_frame = True).values()

# X = X_raw.copy()
X = np.column_stack([X_raw,X_raw.iloc[:,:2]])


S = 1000
rng = np.random.default_rng(123)
test_inds = (tuple(rng.choice(len(X), size=S)), tuple(rng.choice(X.shape[1], size=S)))

X_test = X.copy()
X_test[test_inds] =  np.nan

X.shape
W = np.isnan(X_test)
delta , X_res = soft_impute(X_test, W, rank = 3,maxit = 20)
plt.figure()
plt.plot(delta, 'o-')

plt.figure()
plt.plot(X[W], X_res[W], 'o')
np.sqrt(np.mean((X[W] - X_res[W])**2))
