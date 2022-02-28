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
library(ggplot2)
#####Setting working path 
setwd("/Users/xiazou/Desktop/Tinbergen_Study/Year2/Block 3/Estimation under sparsity/Assignment2/estimation_testing_under_sparsity")
set.seed(123)
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
s_0 = 2
DGP <- function(N,sigma_0,rho,p,s_0){ 
  factor_F = rnorm(N)  
  X =matrix(nrow = N,ncol = p) 
  Y = matrix(nrow = N) 
  eta_all = rnorm(N*p) 
  eta_all = matrix(eta_all,nrow = N,ncol = p) 
  X = sqrt(rho)*factor_F + sqrt(1-rho)*eta_all  
  Y =  rnorm(N,sigma_0)
  for( i in (1:s_0)){ 
    if(i < 3){
      Y = Y+ X[,i]
    }
    else{
      Y = Y + X[,(i+1)]
    }
  }
  return(list(X=X,Y=Y))
}  


# DGP with latent spatially autoregressive process
# First stage is that individuals form links based on the similarity of their
# Xs. Then we generate the Ys based on the variables plus these 
# spatially autocorrelated errors.
DGP_net <- function(N,sigma_0,rho,p, beta, decay_rate){ 
  factor_F = rnorm(N)  
  X =matrix(nrow = N,ncol = p) 
  Y = matrix(nrow = N)
  
  # generate Xs
  eta_all = rnorm(N*p) 
  eta_all = matrix(eta_all,nrow = N,ncol = p) 
  X = sqrt(rho)*factor_F + eta_all 
  
  # generate endogenous network with logistic-distributed random utility shocks
  
  G =matrix(nrow = N,ncol = N)
  for (i in 1:N){
    p_i = exp(beta * (X[,3] - X[i,3]))/(1+exp(beta * (X[,3] - X[i,3])))
    p_i[p_i<=0.5] = 0
    p_i[p_i>=0.5] = 1
    G[,i] = p_i
  }
  
  G[lower.tri(G)] <- 0 # remove lower triangular part so we get symmetry
  G = t(G) + G
  
  # generate Ys by SAR process (see LeSage and Pace 2001) based on above network
  
  I = diag(N)
  linv = (I- decay_rate* G)^-1
  linv[is.infinite(linv)] <- 0
  
  Y = linv %*%(X[,1]+X[,2] + rnorm(N,sigma_0))
  
 X_with_y = cbind(X, t(G)%*%Y)
  return(list(X=X,X2 = X_with_y, Y=Y))
}

######################Function for Lasso regression and confidence interval of beta_1 and beta_3 #### 

###Step 1: For the lasso regression, we use the glmnet function with lasso regularization  
###Step 2: Constract de-sparsified lasso estimator 
###Step3: Constract CI for beta_1 and beta_3 using formula from Theorem 5.2 (Van de Geer(2016))

data_x = DGP(N,sigma_0,rho,p,s_0)$X 
data_y = DGP(N,sigma_0,rho,p,s_0)$Y  
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
  lasso_sq_model = sqrt_lasso(data_x,data_y,lam0 = sqrt(log(p)/n),output_all = TRUE,intercept = FALSE)  
  
  ####Sigma_hat by sqrt-root lasso  
  hat_sigma0 = lasso_sq_model$sigma_hat
  
  
  ####Sigma_hat by sample standard error of residuals from lasso 
  #hat_sigma0 = sqrt(mean(res_lasso^2)) 
  
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
  data_dgp = DGP(N,sigma_0,rho,p,s_0)
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



#####################Extension by changing values of rho, sigma_0, 
#####################p and |S_0|  


#############################Change rho #########################
N=100
p = 200 
sigma_0 = 1 
#rho = 0.8  
s_0 = 2

rho_seq = seq(0.01,0.99,by= 0.05)

rejection_fre_rho_1 = matrix(nrow = length(rho_seq))
rejection_fre_rho_3 = matrix(nrow = length(rho_seq))

for (i in (1:length(rho_seq))){
  rho = rho_seq[i]  
  rejection_table = matrix(data = 1,nrow = M,ncol = length(index))
  
  #####Simulation ###########
  ##Step 1: generate data 
  ##Step 2: calculate the confidence interval using CI_lasso function 
  ##Step 3: check if the CI includes 0, if yes, accept null hypothesis, if no, reject
  ##        1: reject , 0: not reject
  ###end 
  ###Rejection frequency: mean of the rejection table
  
  
  for (m in (1:M)){  
    data_dgp = DGP(N,sigma_0,rho,p,s_0)
    data_x = data_dgp$X 
    data_xy = data_dgp$X2
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
  rejection_fre_rho_1[i,] = colMeans(rejection_table)[1] 
  rejection_fre_rho_3[i,] = colMeans(rejection_table)[2] 

}

#########Plotting ####################### 
rejection_fre_rho_1 = as.data.frame(rejection_fre_rho_1)
rejection_fre_rho_1[,'type'] = 'beta_1'  
rejection_fre_rho_1[,'rho'] = rho_seq  

rejection_fre_rho_3 = as.data.frame(rejection_fre_rho_3)
rejection_fre_rho_3[,'type'] = 'beta_3'  
rejection_fre_rho_3[,'rho'] = rho_seq  

rejection_fre_rho = rbind(rejection_fre_rho_1,rejection_fre_rho_3)
colnames(rejection_fre_rho)[1] = 'rejection_frequency'
png(filename = 'rejection_fre_rho.png')
ggplot(data = rejection_fre_rho,mapping = aes(x = rho, y= rejection_frequency,color = type))+geom_line()+ggtitle('Rejection frequency under different rho')
dev.off()




#############################Change sigma_0 #########################
N=100
p = 200 
#sigma_0 = 1 
rho = 0.8  
s_0 = 2

sigma0_seq = seq(0.01,2,by= 0.1)

rejection_fre_sigma0_1 = matrix(nrow = length(sigma0_seq))
rejection_fre_sigma0_3 = matrix(nrow = length(sigma0_seq))

for (i in (1:length(sigma0_seq))){
  sigma_0 = sigma0_seq[i]  
  rejection_table = matrix(data = 1,nrow = M,ncol = length(index))
  
  #####Simulation ###########
  ##Step 1: generate data 
  ##Step 2: calculate the confidence interval using CI_lasso function 
  ##Step 3: check if the CI includes 0, if yes, accept null hypothesis, if no, reject
  ##        1: reject , 0: not reject
  ###end 
  ###Rejection frequency: mean of the rejection table
  
  
  for (m in (1:M)){  
    data_dgp = DGP(N,sigma_0,rho,p,s_0)
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
  rejection_fre_sigma0_1[i,] = colMeans(rejection_table)[1] 
  rejection_fre_sigma0_3[i,] = colMeans(rejection_table)[2] 
  
}

#########Plotting ####################### 
rejection_fre_sigma0_1 = as.data.frame(rejection_fre_sigma0_1)
rejection_fre_sigma0_1[,'type'] = 'beta_1'  
rejection_fre_sigma0_1[,'sigma0'] = sigma0_seq  

rejection_fre_sigma0_3 = as.data.frame(rejection_fre_sigma0_3)
rejection_fre_sigma0_3[,'type'] = 'beta_3'  
rejection_fre_sigma0_3[,'sigma0'] = sigma0_seq  

rejection_fre_sigma0 = rbind(rejection_fre_sigma0_1,rejection_fre_sigma0_3)
colnames(rejection_fre_sigma0)[1] = 'rejection_frequency'
png(filename = 'rejection_fre_sigma0.png')
ggplot(data = rejection_fre_sigma0,mapping = aes(x =sigma0, y= rejection_frequency,color = type))+geom_line()+ggtitle('Rejection frequency under different sigma0')
dev.off()




#############################Change p #########################
N=100
#p = 200 
sigma_0 = 1 
rho = 0.8  
s_0 = 2

p_seq = seq(100,200,by= 2)

rejection_fre_p_1 = matrix(nrow = length(p_seq))
rejection_fre_p_3 = matrix(nrow = length(p_seq))

for (i in (1:length(p_seq))){
  p = p_seq[i]  
  rejection_table = matrix(data = 1,nrow = M,ncol = length(index))
  
  #####Simulation ###########
  ##Step 1: generate data 
  ##Step 2: calculate the confidence interval using CI_lasso function 
  ##Step 3: check if the CI includes 0, if yes, accept null hypothesis, if no, reject
  ##        1: reject , 0: not reject
  ###end 
  ###Rejection frequency: mean of the rejection table
  
  
  for (m in (1:M)){  
    data_dgp = DGP(N,sigma_0,rho,p,s_0)
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
  rejection_fre_p_1[i,] = colMeans(rejection_table)[1] 
  rejection_fre_p_3[i,] = colMeans(rejection_table)[2] 
  
}

#########Plotting ####################### 
rejection_fre_p_1 = as.data.frame(rejection_fre_p_1)
rejection_fre_p_1[,'type'] = 'beta_1'  
rejection_fre_p_1[,'p'] = p_seq  

rejection_fre_p_3 = as.data.frame(rejection_fre_p_3)
rejection_fre_p_3[,'type'] = 'beta_3'  
rejection_fre_p_3[,'p'] = p_seq  

rejection_fre_p = rbind(rejection_fre_p_1,rejection_fre_p_3)
colnames(rejection_fre_p)[1] = 'rejection_frequency'
png(filename = 'rejection_fre_p.png')
ggplot(data = rejection_fre_p,mapping = aes(x =p, y= rejection_frequency,color = type))+geom_line()+ggtitle('Rejection frequency under different p')
dev.off()




#############################Change |S_0| #########################
N=100
p = 200 
sigma_0 = 1 
rho = 0.8  
#s_0 = 2
s0_seq = seq(2,102,by= 2)

rejection_fre_s0_1 = matrix(nrow = length(s0_seq))
rejection_fre_s0_3 = matrix(nrow = length(s0_seq))

for (i in (1:length(s0_seq))){
  s_0 = s0_seq[i]  
  rejection_table = matrix(data = 1,nrow = M,ncol = length(index))
  
  #####Simulation ###########
  ##Step 1: generate data 
  ##Step 2: calculate the confidence interval using CI_lasso function 
  ##Step 3: check if the CI includes 0, if yes, accept null hypothesis, if no, reject
  ##        1: reject , 0: not reject
  ###end 
  ###Rejection frequency: mean of the rejection table
  
  
  for (m in (1:M)){  
    data_dgp = DGP(N,sigma_0,rho,p,s_0)
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
  rejection_fre_s0_1[i,] = colMeans(rejection_table)[1] 
  rejection_fre_s0_3[i,] = colMeans(rejection_table)[2] 
  
}

#########Plotting ####################### 
rejection_fre_s0_1 = as.data.frame(rejection_fre_s0_1)
rejection_fre_s0_1[,'type'] = 'beta_1'  
rejection_fre_s0_1[,'s0'] = s0_seq  

rejection_fre_s0_3 = as.data.frame(rejection_fre_s0_3)
rejection_fre_s0_3[,'type'] = 'beta_3'  
rejection_fre_s0_3[,'s0'] = s0_seq  

rejection_fre_s0 = rbind(rejection_fre_s0_1,rejection_fre_s0_3)
colnames(rejection_fre_s0)[1] = 'rejection_frequency'
png(filename = 'rejection_fre_s0.png')
ggplot(data = rejection_fre_s0,mapping = aes(x =s0, y= rejection_frequency,color = type))+geom_line()+ggtitle('Rejection frequency under different |S_0|(s_0)')
dev.off()



######################Function for Lasso regression and confidence interval of beta_1 and beta_3 #### 

###Step 1: For the lasso regression, we use the glmnet function with lasso regularization  
###Step 2: Constract de-sparsified lasso estimator 
###Step3: Constract CI for beta_1 and beta_3 using formula from Theorem 5.2 (Van de Geer(2016))

data_x = DGP(N,sigma_0,rho,p,s_0)$X 
data_y = DGP(N,sigma_0,rho,p,s_0)$Y  
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
  lasso_sq_model = sqrt_lasso(data_x,data_y,lam0 = sqrt(log(p)/n),output_all = TRUE,intercept = FALSE)  
  
  ####Sigma_hat by sqrt-root lasso  
  hat_sigma0 = lasso_sq_model$sigma_hat
  
  
  ####Sigma_hat by sample standard error of residuals from lasso 
  #hat_sigma0 = sqrt(mean(res_lasso^2)) 
  
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



#####################Extension by considering different DGP
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
  data_dgp = DGP_net(N,sigma_0,rho,p,beta = 0.5,decay_rate = 0.5)
  data_x = data_dgp$X 
  data_y =data_dgp$Y[,1] 
  
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



#####################Extension by considering different DGP2
#####Number of simulations  
M=200 

index = c(1,3,201)

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
  data_dgp = DGP_net(N,sigma_0,rho,p,beta = 0.5,decay_rate = 0.5)
  data_x = data_dgp$X2 
  data_y =data_dgp$Y[,1] 
  
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


