rm(list=ls())
med_time = 0
set.seed(1234)
options(warn=-1)


#### Simulation 
# loading package and files
library(MASS)
library(sisVIVE)
library(AER)
library(ManyIV)
library("stringr")
source("https://raw.githubusercontent.com/hyunseungkang/invalidIV/master/TSHT.R")
source("https://raw.githubusercontent.com/hyunseungkang/invalidIV/master/JMCode.R")
source("https://raw.githubusercontent.com/xlbristol/CIIV/main/R/CIIV_Functions.R")
source("https://raw.githubusercontent.com/QoifoQ/WIT/master/WIT_Estimator_V1.R")
#####

N_repli = 500
##
Report_Matrix <- list()
Report_Matrix <- matrix(rep(NA,32),nrow=8,ncol = 4)
rownames(Report_Matrix) <- c("TSLS","LIML","oracle_TSLS", "oracle-LIML","TSHT","CIIV", "sisVIVE", "WIT")
colnames(Report_Matrix) <- c("MAD", "CP","FPR","FNR")

## generate data

Const = 1
########
# Case A1(1)
#  
#PI = c(rep(0.5,5),rep(0.6,5))
#p = length(PI)
#alpha = matrix(c(rep(0,5),rep(0.4,3),rep(0.8,2)),p,1)
#n = 200

# Case A1(2)  
# PI = c(rep(0.04,3),rep(0.5,2),0.2,rep(0.1,4))
#  p = length(PI)
# alpha = matrix(c(rep(0,5),1,rep(0.7,4)),p,1)
#  n =500

# Case 1
#PI = c(rep(0.4,21))
#p = length(PI)
#alpha = matrix(c(rep(0,9),rep(0.4,6),rep(0.2,6)),p,1)
#n =1000

# Case 2
PI = c(rep(0.15,21))
p = length(PI)
alpha = matrix(c(rep(0,9),rep(0.4,6),rep(0.2,6)),p,1)
n =1000

########
## Generate Temperate Result

WIT_beta = rep(N_repli);
TSLS_beta = rep(N_repli);
oracle_LIML_beta = rep(N_repli);
sisVIVE_beta = rep(N_repli);
CIIV_beta = rep(N_repli);
TSHT_beta = rep(N_repli);
oracle_TSLS_beta = rep(N_repli);
LIML_beta = rep(N_repli);

WIT_CI_U = rep(N_repli)
WIT_CI_L = rep(N_repli)
WIT_FPR = rep(N_repli)
WIT_FNR = rep(N_repli)


TSHT_CI_U = rep(N_repli)
TSHT_CI_L = rep(N_repli)
TSHT_FPR = rep(N_repli)
TSHT_FNR = rep(N_repli)

CIIV_CI_U = rep(N_repli)
CIIV_CI_L = rep(N_repli)
CIIV_FPR = rep(N_repli)
CIIV_FNR = rep(N_repli)

sisVIVE_FPR = rep(N_repli)
sisVIVE_FNR = rep(N_repli)

TSLS_CI_U = rep(N_repli)
TSLS_CI_L = rep(N_repli) 
oracle_LIML_CI_U = rep(N_repli)
oracle_LIML_CI_L = rep(N_repli)

oracle_TSLS_CI_U = rep(N_repli)
oracle_TSLS_CI_L = rep(N_repli)

LIML_CI_U = rep(N_repli)
LIML_CI_L = rep(N_repli)

WIT_alpha = matrix(rep(0,N_repli*p),nrow = N_repli)
sisVIVE_alpha = matrix(rep(0,N_repli*p),nrow = N_repli)

WIT_valid = matrix(rep(0,N_repli*p),nrow = N_repli)
sisVIVE_valid = matrix(rep(0,N_repli*p),nrow = N_repli)
TSHT_valid = matrix(rep(0,N_repli*p),nrow = N_repli)
CIIV_valid = matrix(rep(0,N_repli*p),nrow = N_repli)

## Start doting replication 
for (replication in 1:N_repli){
  
  #Generate Data

  p=length(PI)
  mu = 30 # Concentration Para.
  Sigma = 0.8*diag(1,p,p)
  
  for(i in 1:p){
    for (j in 1:p){
      if(i !=j){
        Sigma[i,j] = 0.3^(abs(i-j))*0.8
      }
    }
  }
  
  set.seed(med_time+replication)
  Z = mvrnorm(n , rep(0,p), Sigma)
  
  valid_index = which(alpha ==0)
  invalid_index = which(alpha !=0)
  QR_invalid = qr(Z[,invalid_index]);
  Z_valid_adj = qr.resid(QR_invalid,Z[,valid_index]);
  txx = tcrossprod(t(PI[valid_index])%*%t(Z_valid_adj))/(n)
  
  Sigma_U = txx/Const
  Sigma_e = 1;
  Sigma_U = 1
  mu = txx = tcrossprod(t(PI[valid_index])%*%t(Z_valid_adj))/Sigma_U
  
  set.seed(med_time+replication+1)
  Sigma_2 = matrix(c(Sigma_e,0.6*sqrt(Sigma_U*Sigma_e),0.6*sqrt(Sigma_U*Sigma_e),Sigma_U),2,2)
  error = mvrnorm(n,rep(0,2),Sigma_2)
  
  med = rep(0,p)
  med[(alpha ==0)] = rep(1,sum((alpha ==0)))
  true_valid = med
  
  D = Z%*%matrix(PI,p,1)+error[,2] 
  Y =  1*D+Z%*%alpha+error[,1]
  

  
  ## Estimating 
  WIT_est_rough = WIT_practice(D,Y,Z,seq(0.01,0.2,0.02),ini_lam = 0.05,num_trail = 3)
  sis = cv.sisVIVE(Y,D,Z,K = 5,intercept = FALSE,normalize = FALSE)
  ciiv = CIIV(Y,D,Z,robust=FALSE, firststage=T)
  TS = TSHT(Y,D,Z)
  re_TSLS = IVreg(Y~-1+D|Z)
  re_or_TSLS = IVreg(Y~-1+D+Z[,invalid_index]|Z[,valid_index])
  re_LIML = IVreg(Y~-1+D+Z[,invalid_index]|Z[,valid_index],  inference = "re")
  
  ## Recording
  WIT_est = WIT_est_rough[1]$Final
  WIT_beta[replication] = WIT_est$WIT
  WIT_alpha[replication,] = WIT_est$alpha_MCP
  WIT_valid[replication, WIT_est$WIT_valid] = rep(1,length(WIT_est$WIT_valid))
  
  WIT_CI_U[replication] = WIT_est$WIT+1.96*WIT_est$WIT_robust_std
  WIT_CI_L[replication] = WIT_est$WIT-1.96*WIT_est$WIT_robust_std
  WIT_judge = FPR_FNR(true_valid, as.numeric(WIT_valid[replication,]))
  WIT_FPR[replication] = WIT_judge$FPR
  WIT_FNR[replication] = WIT_judge$FNR
  
  TSHT_beta[replication] = TS$betaHat
  TSHT_valid[replication, TS$VHat] = rep(1,length(TS$VHat))
  TSHT_judge = FPR_FNR(true_valid, as.numeric(TSHT_valid[replication,]))
  TSHT_FPR[replication] = TSHT_judge$FPR
  TSHT_FNR[replication] = TSHT_judge$FNR
  TSHT_CI_U[replication] = TS$ci[2]
  TSHT_CI_L[replication] = TS$ci[1]
  
  
  sisVIVE_beta[replication] = sis$beta
  sisVIVE_alpha[replication,] = sis$alpha
  sisVIVE_valid[replication,] = as.numeric(abs(sis$alpha)<0.05)
  
  sisVIVE_judge = FPR_FNR(true_valid, as.numeric(sisVIVE_valid[replication,]))
  sisVIVE_FPR[replication] = sisVIVE_judge$FPR
  sisVIVE_FNR[replication] = sisVIVE_judge$FNR
  
  
  CIIV_beta[replication] = as.numeric(ciiv$Coefficients_CIM)[1]
  ciiv_valid_index = as.numeric(str_sub(ciiv$Valid_Instruments,2))
  CIIV_valid[replication,ciiv_valid_index] = rep(1,length(ciiv_valid_index))
  CIIV_judge = FPR_FNR(true_valid, as.numeric(CIIV_valid[replication,]))
  CIIV_FPR[replication] = CIIV_judge$FPR
  CIIV_FNR[replication] = CIIV_judge$FNR
  CIIV_CI_U[replication] = as.numeric(ciiv$ci_CIM[2])
  CIIV_CI_L[replication] = as.numeric(ciiv$ci_CIM[1])
  
  
  #TSLS_summary = summary(re_TSLS)
  
  TSLS_beta[replication]  = re_TSLS$estimate[2,1]
  TSLS_CI_U[replication] = re_TSLS$estimate[2,1]+1.96*re_TSLS$estimate[2,2]
  TSLS_CI_L[replication]  = re_TSLS$estimate[2,1]-1.96*re_TSLS$estimate[2,2]
  
  LIML_beta[replication]  = re_TSLS$estimate[3,1]
  LIML_CI_U[replication] = re_TSLS$estimate[3,1]+1.96*re_TSLS$estimate[3,2]
  LIML_CI_L[replication]  = re_TSLS$estimate[3,1]-1.96*re_TSLS$estimate[3,2]
  
  oracle_LIML_beta[replication]  = re_LIML$estimate[3,1]
  oracle_LIML_CI_U[replication]  = re_LIML$estimate[3,1]+1.96*re_LIML$estimate[3,2]
  oracle_LIML_CI_L[replication]  = re_LIML$estimate[3,1]-1.96*re_LIML$estimate[3,2]
  
  oracle_TSLS_beta[replication]  = re_or_TSLS$estimate[2,1]
  oracle_TSLS_CI_U[replication]  = re_or_TSLS$estimate[2,1]+1.96*re_or_TSLS$estimate[2,2]
  oracle_TSLS_CI_L[replication]  = re_or_TSLS$estimate[2,1]-1.96*re_or_TSLS$estimate[2,2]
  
}
## Calculate for each replicate


WIT_CP = (WIT_CI_U>=1)*(WIT_CI_L <=1)
TSHT_CP = (TSHT_CI_U>=1)*(TSHT_CI_L <=1)
CIIV_CP = (CIIV_CI_U>=1)*(CIIV_CI_L <=1)
TSLS_CP = (TSLS_CI_U>=1)*(TSLS_CI_L <=1)
oracle_LIML_CP = (oracle_LIML_CI_U>=1)*(oracle_LIML_CI_L <=1)
oracle_TSLS_CP = (oracle_TSLS_CI_U>=1)*(oracle_TSLS_CI_L <=1)
LIML_CP = (LIML_CI_U>=1)*(LIML_CI_L <=1)

Report_Matrix[,"CP"] = c(mean(TSLS_CP),mean(LIML_CP),mean(oracle_TSLS_CP),mean(oracle_LIML_CP),mean(TSHT_CP),mean(CIIV_CP),NA,mean(WIT_CP))
Report_Matrix[,"MAD"] = c(mad(TSLS_beta,center = 1,constant = 1),mad(LIML_beta,center = 1,constant = 1),mad(oracle_TSLS_beta,center = 1,constant = 1),mad(oracle_LIML_beta,center = 1,constant = 1),mad(TSHT_beta,center = 1,constant = 1),mad(CIIV_beta,center = 1,constant = 1),mad(sisVIVE_beta,center = 1,constant = 1),mad(WIT_beta,center = 1,constant = 1)) 
Report_Matrix[,"FPR"] = c(NA,NA,NA,NA,mean(TSHT_FPR),mean(CIIV_FPR),mean(sisVIVE_FPR),mean(WIT_FPR))
Report_Matrix[,"FNR"] = c(NA,NA,NA,NA,mean(TSHT_FNR),NA,mean(CIIV_FNR),mean(WIT_FNR))

Report_Matrix


