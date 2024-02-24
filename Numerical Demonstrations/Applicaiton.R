# Clear the environment
rm(list = ls())

# Load Estimator
source("https://raw.githubusercontent.com/QoifoQ/WIT/master/WIT_Estimator_V1.R")

# Load necessary libraries
library(MASS)


# Set the number of observations
n = 105239

# TRUE parameter values
gamma = c(rep(0.001, 96))
p = length(gamma)
alpha = matrix(c( 0,  0 , 0,  2.578878e-06  ,0 , 0 , 0 , 0,  0,  0,  0,  0,
                  9.121386e-05, -1.591784e-04,  0 , 0,  0 , 0 , 0 , 0,  0, -3.296838e-04 , 0 , 0,
                   0 , 0,  0,  0,  0,  0,  0,  3.415639e-04,  0, -5.777665e-04,  0, -2.024845e-04,
                   0,  0,  0,  0,  0, -2.757444e-04,  0,  0 , 0,  0 , 0,  0,
                   0,  0,  0,  0,  0,  0,  0,  5.713716e-05,  0,  0,  0,  0,
                   0,  0,  0,  7.911372e-05,  0,  0,  0,  0,  0,  0,  0,  0,
                   0,  0, 0,  0,  0 , 0,  0,  0,  8.042636e-04,  0,  7.506719e-05,  0,
                   0,  0,  0,  0,  0,  0,  7.971449e-05,  0,  0,  0,  3.665708e-04,  1.427829e-04), p, 1)

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
Sigma_e = 0.05
Sigma_U = 0.05


# Define the covariance matrix for the error term
Sigma_2 = matrix(c(Sigma_e, 0.6 * sqrt(Sigma_U * Sigma_e), 0.6 * sqrt(Sigma_U * Sigma_e), Sigma_U), 2, 2)

# Generate the error term with multivariate normal distribution
error = mvrnorm(n, rep(0, 2), Sigma_2)


# Generate simulated data
D = Z%*%matrix(gamma,p,1)+error[,2] 
Y =  0.123*D+Z%*%alpha+error[,1]

# Estimate 
WIT_Results = WIT_practice(D,Y,Z,exp(seq(-9,-4,0.5)),ini_lam = 0.05,num_trail = 2)

# Show estimated results
WIT_Results$Final
