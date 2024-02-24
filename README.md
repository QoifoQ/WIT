# WIT Estimator

## Overview

This repository hosts the implementation of the Weak and some Invalid instruments robust Treatment (WIT) estimator, as proposed by Lin et al. (2022). The estimator is designed to address the challenges of estimating treatment effects in the presence of many weak and potentially some invalid instrumental variables (IVs). 

The paper detailing the theory and methodology behind the WIT estimator is publicly available and can be accessed at [https://arxiv.org/abs/2207.03035](https://arxiv.org/abs/2207.03035).

## Required packages 

You can install the required package with 

``` r
install.pacakges(c("ManyIV,ncvreg,ivmodel"))
```

## Repository Contents
This repository contains the R scripts used for the numerical demonstrations presented in our paper. Each script is associated with a specific section of the paper, allowing for replication of the results as described below:

- `Numerical Demonstrations/`
  - `Example1.R`: Generates the results for Example 1 presented in Section 2.3 of the paper.
  - `Simulation.R`: Replicates the Case 1 simulation results from Section 4. To reproduce other simulation scenarios discussed in the paper, you can modify the coefficients as needed within the script.
  - `Application.R`: Provides a toy example that mimics the real data analysis detailed in Section 5. This script also computes the WIT estimate for the toy dataset.

## Example

The following is example 1 considered in Section 2.3 of our article. 

```{r example}
rm(list = ls())
library(MASS)

library(WITEstimator)
n = 500
## TRUE Value 
gamma = c(rep(0.04,3),rep(0.5,2),0.2,rep(0.1,4))
p = length(gamma)
alpha = matrix(c(rep(0,5),1,rep(0.7,4)),p,1)

#####
# IVs: Z1-Z5 are valid; Z6-Z10 are invalid.
#####  
p=length(gamma)
Sigma = 0.8*diag(1,p,p)
    
for(i in 1:p)
{
    for (j in 1:p)
    {
      if(i !=j)
      {
        Sigma[i,j] = 0.3^(abs(i-j))*0.8
      }
    }
}
Z = mvrnorm(n , rep(0,p), Sigma)
    
valid_index = which(alpha ==0)
invalid_index = which(alpha !=0)
QR_invalid = qr(Z[,invalid_index]);
Z_valid_adj = qr.resid(QR_invalid,Z[,valid_index]);
txx = tcrossprod(t(gamma[valid_index])%*%t(Z_valid_adj))/(n)
    
Sigma_e = 1;
Sigma_U = 1
mu = txx = tcrossprod(t(gamma[valid_index])%*%t(Z_valid_adj))/Sigma_U
    
Sigma_2 = matrix(c(Sigma_e,0.6*sqrt(Sigma_U*Sigma_e),0.6*sqrt(Sigma_U*Sigma_e),Sigma_U),2,2)
error = mvrnorm(n,rep(0,2),Sigma_2)
```


Conduct estimation on simulated data using WIT estimator with MCD tuning.
```{r}
## Simulated data
D = Z%*%matrix(gamma,p,1)+error[,2] # The remaining is the intercept
Y =  1*D+Z%*%alpha+error[,1]

## Estimate 
WIT_practice(D,Y,Z,seq(0.1,0.25,0.02),ini_lam = 0.05,num_trail = 4)

## Results
#$WIT
#[1] 0.9976033

#$WIT_robust_std
#[1] 0.06231824

#$WIT_CI
#[1] 0.8754595 1.1197470

#$alpha_MCP
# [1] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 1.0172806 0.6735095 0.7529272 0.6534325 0.6735992

#$WIT_valid
#[1] 1 2 3 4 5

#$alpha_WIT
# [1] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 1.0219745 0.6751280 0.7539500 0.6545832 0.6741253

#$Sargen_p_value
#[1] 0.3611227

#$`Modified-CD`
#[1] 0.3822153
```


## Reference 
Lin, Y., Windmeijer, F., Song, X., & Fan, Q. (2022). On the instrumental variable estimation with many weak and invalid instruments. arXiv preprint arXiv:2207.03035.
