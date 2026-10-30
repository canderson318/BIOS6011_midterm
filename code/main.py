from re import L

import numpy as np
from plotnine import *
from pipeline import run, sweep
import os
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
from dataclasses import dataclass
from scipy import stats

os.chdir(Path.home() / 'Documents/school/BIOS6011/assignments/midterm')

# \\\\\\
# \\\\\\
# \\\\\\
fc = 5
seeds = [120,100,113,1111,245]

# \\\\
# ———— Tune p-count
# \\\\

pcounts = np.geomspace(1e-8, 1e-2, 30)
RES_p = sweep('pcount', {'pcount': pcounts}, fc=fc, seeds=seeds)

p = (
    ggplot(data = RES_p, mapping = aes(x = 'sens', y = 'fdr'))
    + geom_hline(yintercept = 0.05, linetype = '--', color = 'green')
    + geom_point(aes(color = 'pcount'), alpha = .7, size = 5)
    + scale_color_continuous(trans = 'log10')
    + labs(title="FDR versus Sensitivity at different p-counts", x = 'Sensitivity', y = 'FDR')
    + theme_classic()
)
p.save('out/pcount_sens_fdr_scatter.pdf', height = 5, width = 6)

optim_p = RES_p.loc[RES_p['fdr'].lt(.05)].sort_values("sens", ascending = False).iloc[0]

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
optim_i

# \\\\
# \\\\
# ———— Plot fdr/sens vs fc for both
# \\\\
# \\\\
sens_fdrs = []
for fc in np.linspace(1,10,100):
    s, delta, res, scores =  run(fc = fc, pcount = optim_p.pcount, col = 'q', lam = optim_i.lam, rank_frac = optim_i.rank_frac, \
        seed = seeds, effect_frac = .1, verbose = False)
    sens_fdrs.append(scores[['sens', 'fdr']].assign(fc = fc))

sens_fdrs_d = pd.concat([x.reset_index() for x in sens_fdrs], axis = 0)

sens_fdrs_d[sens_fdrs_d.fc == 1]
raw = sens_fdrs_d.copy()
raw['fc_group'] = (raw['fc']*2).round()/2

plt.figure(figsize = (6, 4))
colors = {'pcount': 'C0', 'imp': 'C1'}
styles = {'sens': ('-', 'o'), 'fdr': ('--', '^')}
for meth in ['pcount', 'imp']:
    dm = raw.loc[raw['index'].eq(meth)]
    for metric, (ls, mk) in styles.items():
        g = dm.groupby('fc_group')[metric].agg(['mean', 'std', 'count']).dropna()   # drops groups with NaN/1 run
        half = stats.t.ppf(.975, g['count'] - 1) * g['std'] / np.sqrt(g['count'])  # 95% CI half-width
        plt.errorbar(g.index, g['mean'], yerr = half, linestyle = ls, marker = mk,
                    markersize = 3, capsize = 2, color = colors[meth], label = f"{meth} {metric}")
plt.axhline(.05, color = 'green', linestyle = ':', linewidth = 1)   # nominal FDR
plt.ylabel("Sensitivity (solid) / FDR (dashed)")
plt.xlabel("Fold Change")
plt.legend()
plt.tight_layout()
plt.savefig("out/sens_fdr_fc.pdf")
