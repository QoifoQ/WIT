# WIT Estimator

This is the github repo for the WIT estimator proposed by LIN et al. (2022) that estimate treatment effect using weak and some invalid IVs. The proposed WIT estimator is robust to many individually weak and some invalid IVs (high or low dimension). This paper is now avaible at: 

## Required packages 

You can install the required package with 

``` r
install.pacakges(c("ManyIV,ncvreg,ivmodel"))
```


## Example

This is a example 1 show in section 2.2 of our article. 

```{r example}
rm(list = ls())
library(MASS)

library(WITEstimator)

    n = 500
    gamma = c(rep(0.04,3),rep(0.5,2),0.2,rep(0.1,4))
    p = length(PI)
    alpha = matrix(c(rep(0,5),1,rep(0.7,4)),p,1)
    n =500
  
    p=length(PI)
    Sigma = 0.8*diag(1,p,p)
    
    for(i in 1:p){
      for (j in 1:p){
        if(i !=j){
          Sigma[i,j] = 0.3^(abs(i-j))*0.8
        }
      }
    }
    Z = mvrnorm(n , rep(0,p), Sigma)
    
    valid_index = which(alpha ==0)
    invalid_index = which(alpha !=0)
    QR_invalid = qr(Z[,invalid_index]);
    Z_valid_adj = qr.resid(QR_invalid,Z[,valid_index]);
    txx = tcrossprod(t(PI[valid_index])%*%t(Z_valid_adj))/(n)
    
    Sigma_e = 1;
    Sigma_U = 1
    mu = txx = tcrossprod(t(PI[valid_index])%*%t(Z_valid_adj))/Sigma_U
    
    Sigma_2 = matrix(c(Sigma_e,0.6*sqrt(Sigma_U*Sigma_e),0.6*sqrt(Sigma_U*Sigma_e),Sigma_U),2,2)
    error = mvrnorm(n,rep(0,2),Sigma_2)
    
    med = rep(0,p)
    med[(alpha ==0)] = rep(1,sum((alpha ==0)))
    true_valid = med #true valid index
    

```


Estimate simulated data with WIT estimator 
```{r}
    
    
D = Z%*%matrix(PI,p,1)+error[,2] # The remaining is the intercept
Y =  1*D+Z%*%alpha+error[,1]

WIT_practice(D,Y,Z,seq(0.1,0.25,0.02),ini_lam = 0.05,num_trail = 4)
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
