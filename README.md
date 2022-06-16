# Incremental Latent Confounder search with VBEM (ILC-V) and Hill-Climbing Latent Confounder search with VBEM (HCLC-V) algorithms

The algorithms for the preprint "Discovery and density estimation of latent confounders in Bayesian networks with evidence lower bound".
instruction


## R Library Requirement
R >= 3.3.0, 
[stringr](https://cran.r-project.org/web/packages/stringr/),
[rJava](https://cran.r-project.org/web/packages/rJava/index.html) 

## Installation
- Install Java development Kit version 13.0.1 
```R
Sys.setenv(JAVA_HOME="C:\\Program Files\\Java\\jdk-13.0.1")
```
- Install rjava version 0.9:
```R
require(devtools)
install.packages("stringr")
install.packages("rJava",type='source',"https://cran.r-project.org/src/contrib/Archive/rJava/rJava_0.9-12.tar.gz") 
```
- Install r-causal from github:

```R
install_github("bd2kccd/r-causal", INSTALL_opts=c("--no-multiarch"),force = TRUE)
```

# Instructions

1) Inputs
- algo list = 'GFCI', 'FCI', 'FCI_tetrad'
- case list = 'PROPERTY','asia' , 'Sports','Alarm'

2) Run the files HCLC-V.R and ILC-V.R for both algorithms
3) Outputs

- '*bestdag_*.csv' is the best-found DAG
- '*bestBIC_*.csv' is the BIC score of the best-found DAG
- '*bestLL_*.csv' is the log-likelihood score of the best-found DAG
- '*bestELBO_*.csv' is the p-ELBO score of the best-found DAG
- '*dim_*.csv' is the free parameters of the best-found DAG
- '*list_confounder_*.csv' is the numbers of latent confounders being searched  in the step 4 of both algorithms
- '*listELBO_*.csv' is the p-ELBO score being searched in the step 4 of both algorithms

## Citation

[https://arxiv.org/abs/2206.05490](https://arxiv.org/abs/2206.05490)




