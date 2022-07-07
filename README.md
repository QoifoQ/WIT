# WIT Estimator

This is the github repo for the WIT estimator proposed by LIN et al. (2022) that tests overidentifying restrictions (or instrumental variable validity) with high-dimensional data, which is robust to heteroskedastic errors. The method is called the Q test since it is based on estimation and inference for quadratic functionals of high-dimensional vectors. The paper is available at .

## Required packages 

You can install the required package with 

``` r
install.pacakges(c("ManyIV,ncvreg,ivmodel"))
```

## Install in R

You can install package WITEstimator by the following command in R.

```R
devtools::install_github("https://github.com/QoifoQ/WIT/WITEstimator")
```

## Example

This is a example 1 show in section 2.2 of our article. 

```{r example}
rm(list = ls())
library(MASS)

source("Lasso.R")
source("Projection.R")
source("QTest.R")

n = 300
px = 150
pz = 100 

p = px + pz
Sigma = diag(p)

phi <- c(seq(0.1,0.5,0.1),rep(0,145))
psi <- c(seq(0.3,0.7,0.1),rep(0,145))

s1 = 10 # num of relevant IV
gamma <- matrix(c(rep(0.5,10),rep(0,pz-10)),ncol = 1)

beta = 1 

set.seed(2022)

## covariates and instruments 
W <- MASS::mvrnorm(n=n, rep(1, p), Sigma)
if (px == 0){
  X <- NULL
  Z <- W
}else{
  X <- W[,1:px]
  Z <- W[,(px+1):p]
}

## error terms 
err <- MASS::mvrnorm(n=n, rep(0, 2), matrix(1.5 * c(1, .5, .5,1),2))
e <- err[,1]
eps2 <- err[,2]
```


Test results with valid instruments 
```{r}
D <-  X%*%psi  + Z %*% gamma + eps2
Y <-   D*beta + X%*%phi + e

QTest(Y,D,Z,X) 
# $invalid
# [1] FALSE
# 
# $sig.level
# [1] 0.05
# 
# $pval
# [1] 0.6518306
```


Test results with invalid instruments 
```{r}
pi <-  c(seq(0.1,1,0.1),rep(0,90))
Y <- D*beta + X%*%phi + Z %*% pi + e

QTest(Y,D,Z,X) 
# $invalid
# [1] TRUE
# 
# $sig.level
# [1] 0.05
# 
# $pval
# [1] 9.769963e-15
```

## Reference 
