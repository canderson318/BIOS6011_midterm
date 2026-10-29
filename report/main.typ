#import "utils.typ": *

#show: conf

//////////////////////
//////////////////////
//////// BODY ////////
//////////////////////
//////////////////////

// Strict maximum of 2 pages (aim for 1.5-2), single spaced, including at least 1 table and 1 figure; references, supplemental material, and the appendix do not count. Use a formal scientific writing style (the README wording is informal in places).

#counter(page).update(0)
#set page(numbering:none)
#title[Comparison Between Pseudocount and Soft-Imputation Left-censored Zero Handling]
#align(center)[By: Christian Anderson\ #today.display("[day] [month repr:long] [year]")]
#underline[*Abstract.*]
Zeros in sparse compositional data are problematic because values lie on a simplex and form spurious corelations which most analyses that assume euclidean geometry do not account for. Thus, non-linear transformations are necessary before data is analyzed. The most common among these is the centered log-ratio, but this only works on non-zero data. I compare two methods to 'impute' zero entries: (1) I fill zeros with a small value, a pseudocount, and (2) I impute zeros using SVD based 'soft' imputation. These approaches make the CLR transform possible and data amenable to statistical testing. To test these methods I generate datasets with a set group effect using a Dirichlet-multinomial model and evaluate pseudocount vs. soft imputation using sensitivity at FDR < 0.05.
#underline[*Keywords.*]
compositional data, zero imputation, centered log-ratio, Dirichlet

#align(center,
	link("https://github.com/canderson318/BIOS6011_midterm",
	image("github.svg", height: .8em))
)



#pagebreak()
#set page(numbering:"1")
#set page(columns:2)

= Introduction
Compositional data lie on the simplex: each value is between 0 and 1, and the values for an observation sum to 1.
As a result, when one variable goes up, another goes down.
For example, if the first value of $[0.8, 0.1, 0.1]$ decreases, the remainder must increase, giving $[0.6, 0.2, 0.2]$.
This compositional constraint leads to spurious negative correlations, which we want to avoid.
One way to avoid them is to transform the simplex, also called Aitchison space, to Euclidean space, which works well with most statistical tests and tools, such as the t-test and PCA.

Most biological compositional data are also sparse.
Most features, for example bacterial taxa, are so rarely observed that the majority of their values are zero.
This means the feature was either not there, or the technical process was not sensitive enough to record it.
The latter is left-censoring.
To convert the simplex into Euclidean space, the most common method is the centered log-ratio (CLR) transformation, which centers the log abundances on their geometric mean:
$ "clr"(x_j) = log(x_j) - "mean"(log(x)) $
However, $log(0)$ is undefined, and because the data are sparse, the result would contain a large number of missing values.
#TODO[Add citations: Aitchison on compositional data analysis, and a reference for zero handling in microbiome data.]

= Motivation
Two methods exist to deal with these zeros.
First, a small value, a pseudocount, simply replaces the zeros.
A value of $min(x_(x != 0)) / 2$ is often used, but the choice may have large downstream effects.
Second, the zeros can be imputed using soft imputation, which infers values the zeros can take such that the overall variance of the matrix is not altered.
This method iteratively calculates the singular value decomposition (SVD) of all the data and fills the zero cells each time.
When the reconstruction is essentially the same from one iteration to the next, it has converged.
Additionally, the singular values are regularized so that only the most informative singular vectors impute the zero cells.

Performing statistical tests on the transformed matrix is the ultimate goal.
Whether one method is truer to the original data structure is unknown.
Soft imputation is unbiased and uses the structure of the data, but it has the potential to introduce arbitrary values that may reconstruct well on the whole while introducing noise and artificial structure for individual features.
It also uses more parameters, which must be well tuned.
Pseudocounts are simple, but also introduce artificial variance, although unbiasedly.
#TODO[Why this comparison matters in practice (e.g. which method is the default in microbiome pipelines, and what a wrong choice costs downstream).]

== Explored Question
_I hypothesize that soft imputation will better recapitulate the original multivariate structure than a pseudocount._

Specifically, which method gives higher sensitivity (power) at a false discovery rate (FDR) below 0.05 when testing for a group difference in the CLR-transformed data?
#TODO[State the question in one sentence, with the exact outcome measures and the hypothesized direction of the difference.]

= Simulation Setup
To test which method is better, I simulate sparse compositional data and introduce an effect.
I simulate with the zeros built in instead of adding them artificially, as this is more representative of biological data.
I then evaluate FDR, sensitivity, and specificity to determine whether a method has introduced unwanted structure.
The method that has higher sensitivity (power) at an FDR below 0.05 is superior.

Both methods use many parameters, and the simulation step also uses many, so I use a grid search to identify the best ones for each method.

$X_i tilde "Multinomial"(N_i, theta_i)$, proportions $theta_i tilde "Dirichlet"(alpha)$, with $alpha = m dot "conc"$ and $m$ lognormal; library sizes $N_i$ lognormal. State $n$, $p$, the sparsity target and the prevalence filter.

$
	X_i tilde "Multinomial"(N_i, theta_i)
	\ theta_i tilde "Dirichlet"(alpha)
	\ N_i tilde "lognormal"(mu, sigma)
	\ P = X_(i,j) / (sum_(j=1)^p  X_(i,j))
 \ "where" X in RR^(n times p), i in {1, dots, n}, "and" j in {1, dots, p}

$


#TODO[Effect. How many features carry a group effect, what the fold change $"fc"$ is, how the effect is applied to $theta$, and the null scenario ($"fc" = 1$).]
#TODO[Methods. Pseudocount value(s) and the soft-imputation settings tuned ($lambda$ and rank), and how zeros are handled before the CLR transform.]
#TODO[Testing and evaluation. Per-feature t-test, Benjamini-Hochberg adjustment, sensitivity, specificity, and FDR, the number of replicate simulations, and how tuning and evaluation seeds are separated.]
#TODO[Table 1: simulation scenarios and parameter values.]

= Results
#TODO[Summarize the grid-search results: which pseudocount and which $(lambda, "rank")$ gave the highest sensitivity with FDR < 0.05, and how the two methods compare at their best settings. Add the null-scenario calibration and the comparison across fold changes. Refer to every table and figure in the text.]
#TODO[Table 2 (journal-ready): sensitivity, specificity, and FDR for each method at its best setting, with an informative title and units. Use meaningful labels, not code names (e.g. "Sensitivity", not `sens`).]
#TODO[Figures: replace code names on the axes, colorbars, and captions (`sens`, `fdr`, `pcount`, `lam`, `rank_frac`, `value`) with readable labels (Sensitivity, FDR, pseudocount, $lambda$, rank fraction). Captions should state what is plotted and the simulation settings (e.g. $"fc"$, seeds). Extra figures and tables go in Supplemental Materials as Figure S1, Table S1, etc.]

#align(center,figure(image("../out/imp_sens_fdr_scatter.pdf", width: 60%),caption:"Soft Imputation grid search for lamdas and ranks."))
#align(center,figure(image("../out/pcount_sens_fdr_scatter.pdf", width: 60%),caption:"Pcount grid search."))
#align(center,figure(image("../out/sens_fdr_fc.pdf", width: 60%),caption:"Pcount and soft imputation FDR versus sensitivity (power)."))

= Discussion
#TODO(title: "TODO (10 pts)")[Interpretation: interpret the results in the context of the original question. Was the hypothesis supported, and why do the methods differ or not? End with a conclusion.]
#TODO(title: "TODO (5 pts)")[Limitations of the simulation or of the question examined (single data-generating model, how the zeros arise, the range of fold changes, the effect size and sparsity values explored, and tuning on the same simulation model), and next steps.]

= References
#TODO[References in one consistent style (not counted in the page limit). It is fine to have none.]

#set page(columns:1)
= Supplemental Materials
#TODO[Optional extra figures and tables (Figure S1, Table S1, ...), each referenced in the text.]

= Appendix
== Code
#TODO(title: "TODO (25 pts)")[Reproducibility: include the code for every simulation, analysis, figure, and table in the report, with clear comments, and exclude unused code. State the random seed(s) explicitly. The code will be run by the grader, so check that it runs without error from a clean environment and gives exactly the numbers in the report. The assignment names R or SAS but allows other code.]
```python

```
