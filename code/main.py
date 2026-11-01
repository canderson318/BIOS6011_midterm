from re import L
import numpy as np
from plotnine import *
from pipeline import run, sweep
import os
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
from dataclasses import dataclass
from tqdm import tqdm
from scipy import stats


if __name__ == '__main__':
    wd = Path.home() / 'Documents/school/BIOS6011/assignments/midterm'
    os.chdir(wd)

    # \\\\\\
    # \\\\\\
    # \\\\\\
    fc = 5
    seeds = [120,100,113,1111,245]

    # \\\\
    # ———— Tune p-count
    # \\\\

    pcounts = np.geomspace(1e-8, 1e-2, 30)
    RES_p = sweep('pcount', {'pcount': pcounts}, fc=fc, seeds=seeds )

    p = (
        ggplot(data = RES_p, mapping = aes(x = 'sens', y = 'fdr'))
        + geom_hline(yintercept = 0.05, linetype = '--', color = 'green')
        + geom_point(aes(color = 'pcount'), alpha = .7, size = 5)
        + scale_color_continuous(trans = 'log10')
        + labs(title="FDR versus Sensitivity at different p-counts", x = 'Sensitivity', y = 'FDR')
        + theme_classic()
    )
    p.save('out/pcount_sens_fdr_scatter.pdf', height = 2.5, width = 6)

    optim_p = RES_p.loc[RES_p['fdr'].lt(.05)].sort_values("sens", ascending = False).iloc[0]
    optim_p.to_frame().T.to_csv("out/optim_p.csv", index = False)

    # \\\\
    # ———— Tune imp
    # \\\\

    lams = np.geomspace(.01, .9, 10)  # fraction of the largest singular value
    ranks = np.geomspace(0.01, .9, 10)
    RES_i = sweep('imp', {'lam': lams, 'rank_frac': ranks}, fc=fc, seeds=seeds)

    def panel(var, trans = 'identity'):
        lab = ''
        if var == 'lam':
            lab =  "Lambda"
        elif var == 'rank_frac':
            lab =  "Rank Fraction"
        return (
            ggplot(data = RES_i, mapping = aes(x = 'sens', y = 'fdr'))
            + geom_hline(yintercept = 0.05, linetype = '--', color = 'green')
            + geom_point(aes(color = var), alpha = .7, size = 5)
            + scale_color_continuous(trans = trans)
            + labs(title = lab, x = 'Sensitivity', y = "FDR")
            + theme_classic()
        )
    p = panel('lam', 'log10') / panel('rank_frac')
    p.save('out/imp_sens_fdr_scatter.pdf', height = 5, width = 6)

    optim_i = RES_i.loc[RES_i.fdr.eq(RES_i.fdr.min())].sort_values("sens", ascending = False).iloc[0]
    optim_i.to_frame().T.to_csv("out/optim_i.csv", index = False)

    X = RES_i.select_dtypes("number")
    corr = X.corr(method="spearman")
    PREC = np.linalg.inv(corr)
    D_inv = np.diag(1 / np.sqrt(np.diag(PREC)))
    partial_corr = pd.DataFrame(-D_inv @ PREC @ D_inv, index=X.columns, columns=X.columns)

    n, k = len(X), X.shape[1] - 2
    df = n - 2 - k

    # table over unique pairs only (upper triangle)
    iu = np.triu_indices(X.shape[1], k=1)
    r = partial_corr.values[iu]
    t = r * np.sqrt(df / (1 - r**2))
    res = pd.DataFrame({"var1": X.columns[iu[0]], "var2": X.columns[iu[1]],"r": r, "p": 2 * stats.t.sf(np.abs(t), df),})
    res["q"] = stats.false_discovery_control(res.p, method="bh")   # scipy >= 1.11
    res["signif"] = res.q < .05
    print(res.sort_values("p").round(3))

    # \\\\
    # ———— Plot fdr/sens vs fc for both
    # \\\\
    # \\\\
    sens_fdrs = []
    FCs = np.linspace(0,10,100)

    for fc in tqdm(FCs):
        s, delta, res, scores =  run(fc = fc, pcount = optim_p.pcount, col = 'q', lam = optim_i.lam, rank_frac = optim_i.rank_frac, \
            seed = seeds, verbose = False)
        sens_fdrs.append(scores[['sens', 'fdr']].assign(fc = fc))
    sens_fdrs_d = pd.concat([x.reset_index() for x in sens_fdrs], axis = 0)

    raw = sens_fdrs_d.copy()
    raw['fc_group'] = (raw['fc']*2).round()/2
    raw.columns = ['index', 'Sensitivity', 'FDR', 'fc', 'fc_group']

    plt.figure(figsize = (6, 4))
    colors = {'pcount': 'C0', 'imp': 'C1'}
    styles = {'Sensitivity': ('-', 'o'), 'FDR': ('--', '^')}
    for meth,lab in zip(['pcount', 'imp'], ['Pseudocount', 'Soft-Imputation']):
        dm = raw.loc[raw['index'].eq(meth)]
        for metric, (ls, mk) in styles.items():
            g = dm.groupby('fc_group')[metric].agg(['mean', 'std', 'count']).dropna()   # drops groups with NaN/1 run
            half = stats.t.ppf(.975, g['count'] - 1) * g['std'] / np.sqrt(g['count'])  # 95% CI half-width
            plt.errorbar(g.index, g['mean'], yerr = half, linestyle = ls, marker = mk,
                        markersize = 3, capsize = 2, color = colors[meth], label = f"{lab} {metric}")
    plt.axhline(.05, color = 'green', linestyle = ':', linewidth = 1)   # nominal FDR
    plt.axhline(0, color = 'black', linestyle = '-', linewidth = 1)   # nominal FDR
    plt.ylabel("Sensitivity (solid) / FDR (dashed)")
    plt.xlabel("Fold Change")
    plt.xscale("linear")
    plt.xticks(range(FCs.min().astype(int), FCs.max().astype(int)+1, 1))
    plt.legend()
    plt.tight_layout()
    plt.savefig("out/sens_fdr_fc.pdf")


    print("*** Optimal parameters ***")
    print(f"lambda = {optim_i.lam:.3f}, rank fraction = {optim_i.rank_frac:.3f}")
    print(f"pcount = {optim_p.pcount:.3e}")
    s, delta, res, scores =  run(fc = .0001, pcount = optim_p.pcount, col = 'q', lam = optim_i.lam, rank_frac = optim_i.rank_frac, \
        seed = np.arange(100), effect_frac = .1, verbose = False)
    scores
    scores[['tp', 'fp', 'fn', 'tn', 'sens', 'spec', 'fdr']].round(2).to_csv("out/table.csv")

    d = sens_fdrs_d.reset_index(drop = True)
    d = d[d.fdr < .05]

    d.loc[d.groupby("index").sens.idxmax()].rename({'index':'method'}, axis= 1).to_csv("out/best_sens_under_fdr.csv")

    optim_i
