import numpy as np
from plotnine import *
from pipeline import run, sweep
import os
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
from dataclasses import dataclass

os.chdir(Path.home() / 'Documents/school/BIOS6011/assignments/midterm')

# \\\\\\
# \\\\\\
# \\\\\\
fc = 5
seeds = [120,100,113]

# \\\\
# \\\\
# ———— Tune p-count
# \\\\
# \\\\

pcounts = np.geomspace(1e-8, 1e-2, 30)
RES_p = sweep('pcount', {'pcount': pcounts}, fc=fc, seeds=seeds)

(
    ggplot(data = RES_p.melt(id_vars='pcount'), mapping = aes(x = 'pcount', y = 'value'))
    + geom_point()
    + geom_line()
    + facet_wrap("variable", ncol = 1, scales = 'free_y')
    + theme_bw()
    + scale_x_log10()
)


p = (
    ggplot(data = RES_p, mapping = aes(x = 'sens', y = 'fdr'))
    + geom_hline(yintercept = 0.05, linetype = '--', color = 'green')
    + geom_point(aes(color = 'pcount'), alpha = .7, size = 5)
    + scale_color_continuous(trans = 'log10')
    + labs(title="FDR versus Sensitivity at different p-counts")
    + theme_classic()
)
p.save('out/pcount_sens_fdr_scatter.pdf', height = 5, width = 6)

optim_p = RES_p.loc[RES_p['fdr'].lt(.05)].sort_values("sens", ascending = False).iloc[0]

# \\\\
# \\\\
# ———— Tune imp
# \\\\
# \\\\

lams = np.geomspace(.01, .9, 10)  # fraction of the largest singular value
ranks = np.geomspace(0.01, .9, 10)
RES_i = sweep('imp', {'lam': lams, 'rank_frac': ranks}, fc=fc, seeds=seeds)
(
    ggplot(data = RES_i.melt(id_vars=['lam', 'rank_frac']), mapping = aes(x = 'lam', y = 'value', color = 'factor(rank_frac)', group = 'rank_frac'))
    + geom_point()
    + geom_line()
    + scale_x_log10()
    + facet_wrap("variable", ncol = 1, scales = 'free_y')
    + theme_bw()
)

p = (
    ggplot(data = RES_i.melt(id_vars= ['fdr', 'sens', 'spec']), mapping = aes(x = 'sens', y = 'fdr'))
    + geom_hline(yintercept = 0.05, linetype = '--', color = 'green')
    + scale_color_continuous(trans = 'log10')
    + geom_point(aes(color = 'value'), alpha = .7, size = 5)
    + facet_wrap("variable", ncol = 1)
    + labs(title="FDR versus Sensitivity at different lambda and ranks")
    + theme_classic()
)
p.save('out/imp_sens_fdr_scatter.pdf', height = 5, width = 6)
p

optim_i = RES_i.loc[RES_i.fdr < .05].sort_values("sens", ascending = False).iloc[0]


# \\\\
# \\\\
# ———— Plot fdr/sens vs fc for both
# \\\\
# \\\\
sens_fdrs = []
for fc in np.linspace(1.1,5,100):
    s, delta, res, scores =  run(fc = fc, pcount = optim_p.pcount, col = 'q', lam = optim_i.lam, rank_frac = optim_i.rank_frac, \
        seed = seeds, effect_frac = .1)
    sens_fdrs.append(scores[['sens', 'fdr']].assign(fc = fc))

sens_fdrs_d = pd.concat([x.reset_index() for x in sens_fdrs], axis = 0)

plt.figure(figsize = (6, 4))
colors = {'pcount': 'C0', 'imp': 'C1'}
styles = {'sens': ('-', 'o'), 'fdr': ('--', '^')}   # metric -> (linestyle, marker)
for meth in ['pcount', 'imp']:
    d = sens_fdrs_d.loc[sens_fdrs_d['index'].eq(meth)]
    for metric, (ls, mk) in styles.items():
        plt.plot(d['fc'], d[metric], linestyle = ls, marker = mk,markersize = 3, color = colors[meth],
                label = f"{meth} {metric}")
plt.axhline(.05, color = 'green', linestyle = ':', linewidth = 1)   # nominal FDR
plt.ylabel("Sensitivity (solid) / FDR (dashed)")
plt.xlabel("Fold Change")
plt.ylim(0, 1)
plt.legend()
plt.tight_layout()
plt.savefig("out/sens_fdr_fc.pdf")
