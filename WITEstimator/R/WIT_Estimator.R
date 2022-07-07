## WIT-Estimator

library(ncvreg)
library(ManyIV)
library(ivmodel)

#########################################################
##################### Description #######################
#########################################################

# In the following description

# Y is a nx1 vector standards for responds
# Z is a nxp vector means potential p-IVs that contains
# valid and invalid IVs
# D is a nx1 vector means endogenous varibles of interest

#########################################################

## Helper function
## MCP with two parameter lambda and kappa, kappa = 1/alpha, which alpha is usually = 2;
mcp = function(x,lambda,kappa)
{
  (lambda*abs(x)-kappa*x^2/2)*(abs(x)<=lambda/kappa) + (lambda^2/(2*kappa))*(abs(x)>lambda/kappa)
}

## Derivative of MCP penalty
dmcp = function(x,lambda,kappa)
{
  (lambda-kappa*abs(x))*(abs(x)<=(lambda/kappa))
}

## Operatror of SoftThresholding
SoftThresholding = function(X,Lambda)
{
  sign(X)*pmax(0,abs(X)-Lambda)
}


dL = function(nXX,nXY,beta0)
{
  (-nXY+nXX%*%beta0)
}

## Transofmred Loss
L = function(nXX,nXY,nYY,x)
{
  nYY-2*crossprod(x,nXY)+crossprod(x,(nXX%*%x))
}

FPR_FNR = function(true_valid_index,estimate_valid)
{
  p = length(estimate_valid)
  true_invalid_index = rep(1,p)-true_valid_index
  estimate_invalid = rep(1,p)-estimate_valid

  TP = sum(true_valid_index*estimate_valid)
  FP = sum(true_invalid_index*estimate_valid)
  FN = sum(true_valid_index*estimate_invalid)
  TN = sum(true_invalid_index*estimate_invalid)

  result = list("FPR" = (FP/(FP+TN)),"FNR" = (FN/(FN+TP)))
  return(result)
}

ALL_oracle = function(alpha,PI)
{
  uni = unique(alpha/PI)
  number_oracle = length(uni)
  All_oracle_solutions = matrix(rep(0,length(alpha)*number_oracle),number_oracle,length(alpha))
  for (i in 1:number_oracle)
  {
    All_oracle_solutions[i,] = t(alpha+PI*-uni[i])
  }
  return(All_oracle_solutions)
}

######################
## I-LAMM algorithm ##
######################

## LAMM updating (inner iteration)
LAMM = function(x_old,CAP_Lambda,phi,kappa,nYY,nXX,nXY)
{
  x_new = SoftThresholding(x_old-dL(nXX,nXY,x_old)/phi,CAP_Lambda/phi)
  return(x_new)
}

## I-LAMM updating (outter iteration)
ILAMM = function(CAP_Lambda,x0,phi,kappa,nYY,nXX,nXY,error)
{
  x_old = x0
  TMAX = 500;
  for (J in 1:TMAX)
  {
    x_new = LAMM(x_old,CAP_Lambda,phi,kappa,nYY,nXX,nXY)
    judge = (SoftThresholding(dL(nXX,nXY,x_new),abs(CAP_Lambda)))*(x_new ==0)+(dL(nXX,nXY,x_new)+CAP_Lambda*sign(x_new))*(x_new!=0)

    if(norm(judge,"I")<= error)
    {
      print(paste("sub-iteration = ", J))

      return(x_new)
      break
    }
    x_old = x_new
  }
  return(x_new)
}

## FUll-ILAMM algorithm that integrate the previous two terms
FINAL_ILAMM = function(Y,X,beta_old,lambda,kappa,error_c,error_t)
{

  p = ncol(X)
  n = nrow(X)
  nYY = t(Y)%*%Y/n
  nXY = t(X)%*%Y/n
  nXX = t(X)%*%X/n
  iterationMAX = 2000
  phi = max(eigen(nXX)$value)## only compute once to achieve the majorization.
  ## Outer iteration
  for(i in 1: iterationMAX)
  {
    print(paste("iteration = ", i))
    CAP_Lambda = dmcp(abs(beta_old),lambda,kappa)
    if (i==1)
    {
      error = error_c
    }
    else
    {
      error = error_t
    }
    ## complete
    beta_new = ILAMM(CAP_Lambda,beta_old,phi,kappa,nYY,nXX,nXY,error)

    if(norm(beta_new-beta_old,"I") < error )
    {
      print("success I-LAMM algorithm")
      Loss = L(nXX,nXY,nYY,beta_new)
      print(paste("Loss = ", Loss))
      return(beta_new)
      break
    }
    beta_old = beta_new
  }
  return(beta_new)
}

######################
## Main Body of WIT ##
######################

##
#' Title WIT estimator with a predetermined initial points and a fixed tuning parameter lambda
#'
#' @param D is a nx1 vector means endogenous variables of interest
#' @param Y is a nx1 vector standards for responds
#' @param Z is a nxp vector means potential p-IVs that contains valid and invalid IVs
#' @param beta_old is the initial points we choose to use
#' @param lambda is the tuning lambda parameter
#' @param kappa is 1/rho in paper, take 0.5 for default.
#' @param error_c tolerance for contraction region, take 1e-3 for default.
#' @param error_t tolerance for tightening region, take 1e-5 for default.
#'
#' @return WIT estimate
#' @export
#'
#'
WIT = function(D,Y,Z,beta_old,lambda,kappa=0.5,error_c=1e-3,error_t=1e-5)
{
  if( (!is.vector(Y) && !is.matrix(Y)) | !is.numeric(Y)) stop("Y is not a numeric vector.")
  if( (!is.vector(D) && !is.matrix(D)) | !is.numeric(D)) stop("D is not a numeric vector.")
  if( !is.matrix(Z) | !is.numeric(Z)) stop("Z is not a numerical matrix.")
  if(nrow(Z) != length(Y)) stop("The dimension of Y differs from the row dimension of Z")
  if(nrow(Z) != length(D)) stop("The dimension of D differs from the row dimension of Z")
  if(length(D) != length(Y)) stop("The dimension of Y differs from the dimension of D")

  n = nrow(Z); L = ncol(Z);
  p = L

  ## Compute transformed tilde_Z and tilde_Y
  QR = qr(Z); Yhat = qr.fitted(QR,Y); Dhat = qr.fitted(QR,D);
  Z.transformed = Z - Dhat %*% t(Dhat) %*% Z / sum( Dhat^2)
  Y.transformed = Yhat - Dhat * (sum(Dhat * Yhat) / sum( Dhat^2))

  ### SELECTION0-STAGE
  ## Using I-LAMM algorithm to obtain alpha_MCP
  alpha_hat = FINAL_ILAMM(Y.transformed,Z.transformed,beta_old,lambda,kappa,error_c,error_t)
  if (sum(alpha_hat !=0)==p)
  {
    return(NA)
  }
  alphaSuppSize = apply(alpha_hat,1,function(x){sum(x != 0)}); indexProperSupp = which(alphaSuppSize < p);

  invalid_index = which(alphaSuppSize ==1)
  valid_index = which(alphaSuppSize ==0)
  ##

  Z_valid = Z[,valid_index]
  Z_invalid = Z[,invalid_index]

  ### ESTIMATION-STAGE using LIML
  many_est = IVreg(Y~-1+D+Z[,invalid_index]|Z[,valid_index],  inference = "re") ## Depends on ManyIV package
  xxx = ivmodel(Y,D,Z[,valid_index],Z[,invalid_index],intercept=FALSE) # for LIML alpha
  WIT =  many_est$estimate[3,1] #as.numeric(xxx$Fuller$point.est)      #many_est$estimate[3,1] ## liml
  WIT_robost_std = many_est$estimate[3,2]
  test = IVoverid(many_est)

  ### Get LIML result of alpha estimate
  WIT_alpha = rep(0,p);
  WIT_alpha[invalid_index] = xxx$LIML$point.est.other

  WIT_CI = c(WIT-1.96*WIT_robost_std,WIT+1.96*WIT_robost_std)
  ## Summary for Result
  Result = list("WIT" = WIT,"WIT_robust_std" = WIT_robost_std,"WIT_CI" = WIT_CI,"alpha_MCP" = as.numeric(alpha_hat),"WIT_valid" = valid_index,"alpha_WIT" = WIT_alpha,"Sargen_p_value" = test[1,2],"Modified-CD"= test[2,2])
  return(Result)
}

######################################
## Main Body of WIT with MCD tuning ##
######################################

## Tuning Strategy with MCD Test
MCD_tuning = function(D,Y,Z,Lambda_Seq,kappa,error_c,error_t,initial) ##
{
  ## Obtain Transformed Z and Y
  n = nrow(Z); L = ncol(Z);
  p = L
  QR = qr(Z); Yhat = qr.fitted(QR,Y); Dhat = qr.fitted(QR,D);
  Z.transformed = Z - Dhat %*% t(Dhat) %*% Z / sum( Dhat^2)
  Y.transformed = Yhat - Dhat * (sum(Dhat * Yhat) / sum( Dhat^2))

  nYY = t(Y.transformed)%*%Y.transformed/n
  nXY = t(Z.transformed)%*%Y.transformed/n
  nXX = t(Z.transformed)%*%Z.transformed/n


  ## Record MCD value
  MCD_value = rep(0,length(Lambda_Seq))
  Model_size_value = rep(0,length(Lambda_Seq))
  for (ii in (1:length(Lambda_Seq)))
  {
    lambda = Lambda_Seq[ii]

    WIT_Current = WIT(D,Y,Z,initial,lambda,kappa,error_c,error_t)
    if (is.na(WIT_Current))
    {
      next
    }
    Model_size_value[ii] = sum(WIT_Current$alpha_WIT !=0)
    MCD_value[ii] = WIT_Current$`Modified-CD`
  }
  ## Selected the result
  filtered_index = (MCD_value > 0.5/log(n))
  if (sum(is.na(MCD_value)) == length(Lambda_Seq)){
    # not censored
    filtered_index = 1:length(Lambda_Seq)
  }

  filtered_index[is.na(filtered_index)] = rep(T,sum(is.na(filtered_index)))

  if(sum(!filtered_index) == length(Lambda_Seq)) ## all cannot pass MCD
  {
    return(NA)
  }

  Optimal_Model_size_value = which(Model_size_value == min(Model_size_value[filtered_index]))[1]

  WIT_final_estimation = WIT(D,Y,Z,initial,Lambda_Seq[Optimal_Model_size_value],kappa,error_c,error_t)
  return(WIT_final_estimation)
}



#' Title WIT Estimator with automatically MCD tuning strategy.
#'
#' @param D is a nx1 vector means endogenous variables of interest
#' @param Y is a nx1 vector standards for responds
#' @param Z is a nxp vector means potential p-IVs that contains valid and invalid IVs
#' @param Lambda_Seq is a tuning parameter sequence of lambda. Generally take sqrt(p/n)*(0.1,0.2,...,2)
#' @param ini_lam for grouping initial lambda, generally vary from 0.02 to 0.2.
#' @param kappa is 1/rho in paper, take 0.5 for default.
#' @param error_c tolerance for contraction region, take 1e-3 for default.
#' @param error_t tolerance for tightening region, take 1e-5 for default.
#' @param num_trail is upper bound of how many initials we may take
#'
#' @return
#' @export

WIT_practice = function(D,Y,Z,Lambda_Seq,ini_lam = 0.05,kappa=2,error_c=1e-3,error_t=1e-4,num_trail=2)
{
  ## Ratio Estimate
  n= dim(Z)[1]
  p = dim(Z)[2]
  Gamma_est = qr.solve(Z,Y)
  gamma_est = qr.solve(Z,D)

  inidividual_esti = Gamma_est/gamma_est
  ord_ini_est = order (inidividual_esti)
  sorted_ini_est = sort(inidividual_esti)
  XXX = matrix(rep(1,p*p),p,p);
  XXX[lower.tri(XXX)] = 0;

  ## Rough clustering
  fuzed_model = ncvfit(XXX,sorted_ini_est,sorted_ini_est,penalty = "MCP",lambda = ini_lam,penalty.factor = c(0,rep(1,length(inidividual_esti)-1)))
  fuzed_ini_est = XXX%*%fuzed_model$beta
  uni_fuzed_ini_est = unique(fuzed_ini_est)
  number_group = length(uni_fuzed_ini_est)
  group_count = matrix(c(rep(0,2*number_group)),number_group,2)
  group_count[,1] = uni_fuzed_ini_est
  for (i in 1:length(uni_fuzed_ini_est))
  {
    group_count[i,2] = sum(fuzed_ini_est == uni_fuzed_ini_est[i])

  }

  ## Determine how many initials we use
  num_trail = min(num_trail,number_group)
  order_in = order(group_count[,2],decreasing = T)
  group_count = group_count[order_in,]

  ## Compute results
  initials = matrix(rep(0,p*(num_trail+1)),p,num_trail+1)
  for (j in 1:num_trail)
  {
    initials[,j] = Gamma_est-group_count[j,1]*gamma_est
    initials[(ord_ini_est[fuzed_ini_est==group_count[j,1]]),j] = rep(0,group_count[j,2])
  }
  l_0_alpha = p+1

  ## Tuning with MCD
  for (j in 1:(num_trail+1))
  {
    med_model = MCD_tuning(D,Y,Z,Lambda_Seq,kappa,error_c,error_t,initials[,j])
    if (is.na(med_model))
    {
      next
    }
    if (sum(med_model$alpha_WIT!=0)<l_0_alpha)
    {
      l_0_alpha = sum(med_model$alpha_WIT!=0)
      final_model = med_model
    }
  }
  return(final_model)
}




