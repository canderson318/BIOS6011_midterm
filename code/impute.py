import numpy as np


def soft_impute(X, W, maxit=50, rank=10, lam=.1, tol=1e-3, early_stop=True, cap=True):
    """Impute masked entries of X with iterative soft-thresholded SVD in log space.
    Masked entries start at half the smallest observed value in their column. Each iteration centers the
    log matrix by column, takes the SVD, shrinks singular values by lam * (largest singular value), keeps
    the top `rank`, and reconstructs. Observed entries are never changed; imputed values are exp() of the
    reconstruction, so always positive.
    Parameters:
        X: (n, p) non-negative data matrix.
        W: (n, p) boolean mask, True = missing/to impute.
        maxit: max iterations. rank: max singular values kept.
        lam: shrinkage as a fraction of the largest centered singular value (0 = none, 1 = everything removed).
        tol: stop when the change in delta is below this if early_stop is True (delta is in log units).
        early_stop: stop when difference between iterations is < tol.
        cap: cap imputed log values at the column's smallest observed log value (a zero is below detection).
    Returns:
        (delta, X_imputed):  where delta is the squared change in imputed log values per iteration.
    """
    def _reconstruct(L):
        mu = L.mean(0, keepdims=True)
        u, s, vt = np.linalg.svd(L - mu, full_matrices=False)
        s_ = np.maximum(s - lam * s[0], 0)
        return mu + (u[:, :rank] * s_[:rank]) @ vt[:rank]
    logX = np.log(np.where(W, 1, X))  # masked entries are placeholders
    col_min = np.nanmin(np.where(W, np.nan, logX), axis=0)  # log of the smallest observed value per column
    Z = _reconstruct(np.where(W, col_min - np.log(2), logX))  # initial fill with min - log(2) -> min/2
    if cap:
        Z = np.minimum(Z, col_min)
    delta = []
    for _ in range(maxit):
        prev = Z
        Z = _reconstruct(np.where(W, Z, logX))  # imputed where 0, original elsewhere
        if cap:
            Z = np.minimum(Z, col_min)
        delta.append(np.sum((prev[W] - Z[W]) ** 2))
        if early_stop and len(delta) > 1 and np.abs(delta[-1] - delta[-2]) < tol:
            break
    else:
        if len(delta) > 1:
            print(f'No convergence to `tol` at {maxit} iterations')
    return delta, np.where(W, np.exp(Z), X)
