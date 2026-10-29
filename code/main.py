import numpy as np
from plotnine import *
from pipeline import run, sweep
import os
from pathlib import Path

os.chdir(Path.home() / 'Documents/school/BIOS6011/assignments/midterm')

# \\\\\\
# \\\\\\
# \\\\\\
fc = 3
seeds = [100]

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

optim_pcount = RES_p.loc[RES_p.loc[RES_p['fdr'].lt(.05), 'sens'].idxmax(), 'pcount']
optim_pcount

# \\\\
# \\\\
# ———— Tune imp
# \\\\
# \\\\

lams = np.geomspace(.01, .9, 20)  # fraction of the largest singular value
ranks = np.geomspace(0.01, .9, 20)
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


optim_lam, optim_rank_frac, *_ = RES_i.loc[RES_i.fdr < .05].values.ravel()
