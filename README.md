# Welcome to this project:

The repository contains text, analysis code and figure code for the manuscript **Beyond Overall Effects: A Bayesian Approach to Finding Constraints In Meta-Analysis**

1. To see the fill pipeline behind the manuscript, clone this repository.  Do not just download a single file as files are dependent on one another.

2. After you have cloned the repository, go hunt for the manuscript you wish under papers.  For example, "papers/psyMethRev" is the invited revision from Psychological Methods.

3. Go to the file p.Rmd to start.  It is an R Markdown file.  You should be able to execute it.

4. Our `R` scripts depend on libraries from others.  The first block should call all the needed libraries. The following packages are used for the analyses in the manuscript: papaja, spatialfil, tmvtnorm, MCMCpack, msm, gridBase, ggplot2, grid, RColorBrewer, reshape2, BayesFactor, MCMCpack, diagram. All of these packages have to be installed in order to compile the manuscript using `install.packages("packagename")`. To install papaja, please use the following code: `devtools::install.github("crsh/papaja")`


