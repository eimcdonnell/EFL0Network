# EFL0Network

The `EFL0Network` package fits an Elastic Fused $l_0$ dynamic Gaussian graphical model.  
The proposed method is an extension of the Augmented and Penalized Minimization $l_0$ (APM-$L_0$) method$^1$. At each time point of interest, our model uses a neighborhood selection approach$^2$ to estimate a Gassian graphical model with $p$ nodes. Rather than fitting $p$ ordinary lasso regression models as in a traditional neighborhood selection approach, we add kernel weights to borrow information from neighboring time points, and we apply $l_0$, $l_1$, and $l_2$ penalties to promote sparsity and temporal smoothness.  
Source code for this package has been adapted from CRAN R package `APML0`$^3$.

1. Li, X., Xie, S., Zeng, D., and Wang, Y. (2018). Efficient L0-norm feature selection based on augmented and
  penalized minimization. Statistics in Medicine, 37, 473--486.

2. Meinshausen, N. and B&uuml;hlmann, P. (2006). High-dimensional graphs and variable selection with the lasso. Annals of Statistics, 34, 1436--1462.

3. Li X., Xie S., Zeng D., and Wang Y. (2020). APML0: Augmented and Penalized Minimization Method L0 (R package version 0.10). https://cran.r-project.org/web/packages/APML0/index.html.