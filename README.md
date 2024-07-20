# WIT Estimator

## Overview

This repository hosts the implementation of the Weak and some Invalid instruments robust Treatment (WIT) estimator, as proposed by Lin et al. (2024). The estimator is designed to address the challenges of estimating treatment effects in the presence of many weak and potentially some invalid instrumental variables (IVs). 

The paper detailing the theory and methodology behind the WIT estimator is publicly available and can be accessed at [https://academic.oup.com/jrsssb/advance-article-abstract/doi/10.1093/jrsssb/qkae025/7629972](https://academic.oup.com/jrsssb/advance-article-abstract/doi/10.1093/jrsssb/qkae025/7629972).

## Required packages 

You can install the required package with 

``` r
install.pacakges(c("ManyIV,ncvreg,ivmodel"))
```

## Installation 

``` r
install.packages("devtools")
library(devtools)
source("https://raw.githubusercontent.com/QoifoQ/WIT/master/WIT_Estimator_V1.R")
```

## Repository Contents
This repository contains the R scripts used for the numerical demonstrations presented in our paper. Each script is associated with a specific section of the paper, allowing for replication of the results as described below:

- `Numerical Demonstrations/`
  - `Simulation.R`: Replicates the Case 1 simulation results from Section 4. To reproduce other simulation scenarios discussed in the paper, you can modify the coefficients as needed within the script.
  - `Application.R`: Provides a toy example that mimics the real data analysis detailed in Section 5. This script also computes the WIT estimate for the toy dataset.

## Example of Usage

The following is example 1 considered in Section 2.3 of our article. 

```{r example}
# Clear the environment
rm(list = ls())

# Load necessary libraries
library(MASS)
source("https://raw.githubusercontent.com/QoifoQ/WIT/master/WIT_Estimator_V1.R")

# Set the number of observations
n = 500

# TRUE parameter values
gamma = c(rep(0.04, 3), rep(0.5, 2), 0.2, rep(0.1, 4))
p = length(gamma)
alpha = matrix(c(rep(0, 5), 1, rep(0.7, 4)), p, 1)

# Define the covariance matrix for the IVs
Sigma = 0.8 * diag(1, p, p)
for (i in 1:p) {
  for (j in 1:p) {
    if (i != j) {
      Sigma[i, j] = 0.3^(abs(i - j)) * 0.8
    }
  }
}

# Generate instruments (Z) with multivariate normal distribution
Z = mvrnorm(n, rep(0, p), Sigma)

# Define additional parameters for the error term
Sigma_e = 1
Sigma_U = 1


# Define the covariance matrix for the error term
Sigma_2 = matrix(c(Sigma_e, 0.6 * sqrt(Sigma_U * Sigma_e), 0.6 * sqrt(Sigma_U * Sigma_e), Sigma_U), 2, 2)

# Generate the error term with multivariate normal distribution
error = mvrnorm(n, rep(0, 2), Sigma_2)
```


Conduct estimation on simulated data using WIT estimator with MCD tuning.
```{r}
# Generate simulated data
D = Z%*%matrix(gamma,p,1)+error[,2] 
Y =  1*D+Z%*%alpha+error[,1]

# Estimate 
WIT_Results = WIT_practice(D,Y,Z,seq(0.1,0.25,0.02),ini_lam = 0.05,num_trail = 4)

# Show estimated results
WIT_Results$Final

$WIT
[1] 0.9900749

$WIT_robust_std
[1] 0.06454843

$WIT_CI
[1] 0.86356 1.11659

$alpha_MCP
 [1] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.9877521 0.6450006 0.7149057 0.7898840 0.6759620

$WIT_valid
[1] 1 2 3 4 5

$alpha_WIT
 [1] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.9889044 0.6450930 0.7156186 0.7902349 0.6762512

$Sargen_p_value
[1] 0.3421309

$`Modified-CD`
[1] 0.3650936
```


## Reference 
Lin, Y., Windmeijer, F., Song, X., & Fan, Q. (2024). On the instrumental variable estimation with many weak and invalid instruments. Journal of the Royal Statistical Society Series B: Statistical Methodology, qkae025.
