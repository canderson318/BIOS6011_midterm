import numpy as np
import pandas as pd
import scipy.stats as stats


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

def test_features(Z, groupA, is_effect):
    """
    two-sample t-test of group A vs B on the columns of Z.
    Returns:
        DataFrame indexed by feature with t, p, BH q, mean difference a-b and the truth label.
    """
    t, p = stats.ttest_ind(Z[groupA], Z[~groupA], axis=0)
    return pd.DataFrame({
        't': t,
        'p': p,
        'q': stats.false_discovery_control(p, method='BH'),
        'a-b': Z[groupA].mean(0) - Z[~groupA].mean(0),
        'has_effect': is_effect,
    })

def metrics(res, col='p', alpha=.05):
    """
    Confusion metrics.
    """
    det, eff = (res[col] < alpha).to_numpy(), res['has_effect'].to_numpy()
    tp, fp = (det & eff).sum(), (det & ~eff).sum()
    fn, tn = (~det & eff).sum(), (~det & ~eff).sum()
    div = lambda a, b: a / b if b else np.nan
    return {'tp': tp, 'fp': fp, 'fn': fn, 'tn': tn,
            'sens': div(tp, tp + fn), 'spec': div(tn, tn + fp),
            'fdr': div(fp, tp + fp), 'ppv': div(tp, tp + fp), 'npv': div(tn, tn + fn)}
