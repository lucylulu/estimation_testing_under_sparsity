# ETUS Assignment 2
# Xia Zou and Kieran Marray (Depts. of Econometrics and Spatial Economics, VU Amsterdam)
# Simulating some LASSO 

######################################
# Code for a simulation study of ETUS course
# author: @Xia Zou @Kieran Marray
# Date: February 19th 2022 
#######################################  

rm(list = ls())
########Loading packages  
#install.package('RPtests')
library(RPtests) 
library(glmnet)

#####Setting working path 
setwd("/Users/xiazou/Desktop/Tinbergen_Study/Year2/Block 3/advanced time series")

##########Function of Data generating #################### 
# inputs
# N, sigma_0, ??, p, ????? as in the assignment
# nos - number of simulations
# optim - optimiser we use - default value is finite-difference BFGS method
# ?????init - initial values for our betas 
######################################################## 
#####Parameters of data generating 
N=100
p = 200 
sigma_0 = 1 
rho = 0.8  

DGP <- function(N,sigma_0,rho,p){ 
  factor_F = rnorm(N)  
  X =matrix(nrow = N,ncol = p) 
  Y = matrix(nrow = N) 
  eta_all = rnorm(N*p) 
  eta_all = matrix(eta_all,nrow = N,ncol = p) 
  X = sqrt(rho)*factor_F + sqrt(1-rho)*eta_all 
  Y = X[,1]+X[,2] + rnorm(N,sigma_0)
  return(list(X=X,Y=Y))
} 

######################Function for Lasso regression and confidence interval of beta_1 and beta_3 #### 

###Step 1: For the lasso regression, we use the glmnet function with lasso regularization  
###Step 2: Constract de-sparsified lasso estimator 
###Step3: Constract CI for beta_1 and beta_3 using formula from Theorem 5.2 (Van de Geer(2016))

data_x = DGP(N,sigma_0,rho,p)$X 
data_y = DGP(N,sigma_0,rho,p)$Y  
index = c(1,3) 
sig = 0.05
CI_lasso<- function(data_x,data_y,index,sig){   
  ##Index is a tuple about the predictors that you want to have confidence interval
  ##sig: significant level for the confidence interval 
  ####Lasso estimator 
  ####Choose a lambda of sqrt(log(p)/n) for the result of asymptotic
  p = ncol(data_x) 
  n= nrow(data_x)
  #### lamda = sqrt(log(p)/n)
  lasso_model = glmnet(data_x,data_y,alpha = 1,intercept = FALSE,lambda = sqrt(log(p)/n))   
  coef_all = coef(lasso_model) 
  res_lasso = data_y - data_x%*%coef_all[-1] 
  hat_sigma0 = sqrt(mean(res_lasso^2))
  # Construct confidence interval for variables in the index  
  CI_index = matrix(nrow = length(index),ncol = 4) 
  colnames(CI_index) = c('beta_estimator','De-sparsified Estimator','lower_bound','Higher_bound')

  for (j in (1:length(index))){ 
    i = index[j]
    beta = coef_all[i+1] 
    ####The ith column of the surrogate inverse 
    ####Use square root lasso, which is  scale free  
    x_y = data_x[,i] 
    x_x = data_x[,-i]  
    Theta = matrix(nrow = p)
    Theta[i] = 1
    sqrt_lasso_model = sqrt_lasso(x_x,x_y,lam0 = sqrt(log(p)/n),output_all = TRUE,intercept = FALSE)  
    Theta[-i] = -sqrt_lasso_model$beta 
    ####Scale free tau 
    hat_tau2 = sqrt_lasso_model$sigma_hat^2  
    tilde_tau2 = hat_tau2 + sqrt(log(p)/n)*sqrt(hat_tau2)*norm(as.matrix(sqrt_lasso_model$beta),type = 'o')
    Theta = Theta/tilde_tau2 
    
    ###De-sparsified lasso 
    hat_b = beta + (t(Theta)%*%t(data_x)%*%res_lasso)/n  
    critical_value = qnorm(1-sig/2)
    lower_b = hat_b - critical_value* (hat_sigma0)*(hat_tau2)/((tilde_tau2^2)*sqrt(n)) 
    higher_b = hat_b +critical_value* (hat_sigma0)*(hat_tau2)/((tilde_tau2^2)*sqrt(n))  
    
    CI_index[j,] = c(beta,hat_b,lower_b,higher_b)
  } 
  
   
  return(CI_index)
}



#####################Monte Carlo simulation for Lasso regression and confidence interval of beta_1 and beta_3 
#####Number of simulations  
M=200 

rejection_fre = matrix(nrow = 1,ncol = length(index))

rejection_table = matrix(data = 1,nrow = M,ncol = length(index))

#####Simulation ###########
##Step 1: generate data 
##Step 2: calculate the confidence interval using CI_lasso function 
##Step 3: check if the CI includes 0, if yes, accept null hypothesis, if no, reject
##        1: reject , 0: not reject
###end 
###Rejection frequency: mean of the rejection table


for (m in (1:M)){  
  data_dgp = DGP(N,sigma_0,rho,p)
  data_x = data_dgp$X 
  data_y =data_dgp$Y    
  lassoCI = CI_lasso(data_x = data_x,data_y = data_y,index,sig)
  ###Testing for beta 1  
  if ((lassoCI[1,'lower_bound'] <= 0) & (lassoCI[1,'Higher_bound'] >= 0) ){
    rejection_table[m,1] = 0
  } 
  if ((lassoCI[2,'lower_bound'] <= 0) & (lassoCI[2,'Higher_bound'] >= 0) ){
    rejection_table[m,2] = 0
  }
} 

rejection_fre = colMeans(rejection_table)



#####################Monte Carlo simulation for Lasso regression and confidence interval of beta_1 and beta_3 
####By using package 





