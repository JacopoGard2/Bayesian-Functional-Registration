
#--- Creazione dati simulati ---# 

library(splines)
library(MASS)
library(fda)
library(ddalpha)
library(emdbook)
library(RcppArmadillo)
library(splines2)
library(mvnfast)
library(coda)
library(LearnBayes)
library(lattice)



creation_dataset = function( n,N_OBS,p_tilde,q_tilde,true_sigma_eps = 0.01, true_lambda = 5, 
                             true_c0 = 3,true_sigma_c = 1,true_a0 = 1.2,true_sigma_a = 0.5,true_sigma_phi = 2)
{
  p_tilde = 10
  p_true = p_tilde + 4
  q_tilde = 4
  q_true = q_tilde + 4 
  n_knots_m = p_tilde
  n_knots_h = q_tilde 
  
  #-- True sigma epsilon^2 --#
  eps = matrix(NA,n,N_OBS)
  for(i in 1:n){
    eps[i,] =  mvrnorm(n = N_OBS, mu = 0 , Sigma = true_sigma_eps )
  }
  
  #-- True beta --#
  Omega=omega_P(p_true)
  true_beta = mvrnorm(n = 1, mu = rep(0,p_true) , Sigma = solve(Omega/true_lambda) )
  
  #-- True a_i --# 
  cond = TRUE
  while(cond){
    true_a = rnorm(n=n,true_a0,true_sigma_a)
    cond = ! sum(true_a > 0 ) == n
  }
  
  #-- True c_i --# 
  true_c = rnorm(n=n,true_c0,true_sigma_c)
  
  # True phi
  true_phi = matrix(NA,n,q_true)
  for(i in 1:n){
    true_csi_i = list()
    b_g = sample(1:10, size = 1, replace = TRUE)
    A_gi = sample(7:15, size = 1, replace = TRUE)
    true_csi = rep(0,q_true)
    for(j in 2:q_true){
      true_csi[j] = rgamma(n=1,A_gi,b_g)
    }
    true_phi[i,] = cumsum(true_csi) / sum(true_csi)
    
  }
  
  #-- True h_i --#
  time = seq(0,1,length=N_OBS)
  knots_h = seq(0, 1, length.out = q_tilde+2)[-c(1,q_tilde+2)]
  true_Bh <- bs(time, knots = knots_h, intercept = T)
  h_i = true_Bh %*% t(true_phi)
  
  #-- True Bm --#
  time_w = true_Bh %*% t(true_phi)   # time x n 
  knots_m = seq(0, 1, length.out = n_knots_m+2)[-c(1,n_knots_m+2)]
  true_Bm_beta = lapply(1:n, function(i){bs(h_i[,i], knots = knots_m, intercept = T)})
  true_Bm_beta_time_vero = bs(time, knots = knots_m, intercept = T)
  
  y_t = lapply(1:n,function(i){true_a[i]*true_Bm_beta[[i]]%*%true_beta + true_c[i] + eps[i,]})
  y_noise_free = lapply(1:n,function(i){true_a[i]*true_Bm_beta[[i]]%*%true_beta + true_c[i] })
  y_all = lapply(1:n,function(i){true_a[i]*true_Bm_beta_time_vero%*%true_beta+true_c[i]+ eps[i,]})
  y_all_noise_free = lapply(1:n,function(i){true_a[i]*true_Bm_beta_time_vero%*%true_beta+true_c[i]})
  y_true=matrix(NA,nrow=N_OBS,ncol=n)
  y_true_noise_free = matrix(NA,nrow=N_OBS,ncol=n)
  y_true_all=matrix(NA,nrow=N_OBS,ncol=n)
  y_true_all_noise_free=matrix(NA,nrow=N_OBS,ncol=n)
  
  for(i in 1:n){
    y_true[1:N_OBS,i]=y_t[[i]]
    y_true_noise_free[1:N_OBS,i]=y_noise_free[[i]]
    y_true_all[1:N_OBS,i]=y_all[[i]]
    y_true_all_noise_free[1:N_OBS,i]=y_all_noise_free[[i]]
  }
  
  
  y_to_use = y_true
  
  return(list( y = list(y_to_use = y_to_use,y_mis = y_true, y_mis_noise_free = y_true_noise_free,
                        y_al = y_true_all, y_al_noise_free = y_true_all_noise_free),
               params = list(beta = true_beta,  phi = true_phi, a = true_a, c = true_c, sigma_c = true_sigma_c,
                             sigma_a = true_sigma_a, c0 = true_c0, a0 = true_a0,
                             lambda = true_lambda, sigma_phi = true_sigma_phi, sigma_eps = true_sigma_eps)
  ))
}

omega_P <- function(dim){
  K <- Matrix::bandSparse(dim, k=-c(1), diag=list(rep(-1,dim)), symmetric=TRUE)
  diag(K) <- c(rep(2,dim-1),1)
  K <- matrix(as.numeric(K),nrow=dim, byrow=TRUE)
  return (K)
}

n = 30
N_OBS = 300
p_tilde = 10
q_tilde = 4
true_c0 = 3
true_sigma_c = 1
true_sigma_phi = 2 
true_a0 = 1.2
true_sigma_a = 0.5
true_lambda = 1.3
true_sigma_eps = 0.002

seed = -4
set.seed(seed)
telesca_data = creation_dataset(n,N_OBS,p_tilde,q_tilde,true_sigma_eps = true_sigma_eps , true_lambda = true_lambda, 
                                true_c0 = true_c0 ,true_sigma_c = true_sigma_c,true_a0 = true_a0,true_sigma_a = true_sigma_a,
                                true_sigma_phi = true_sigma_phi)

time = seq(0,1,length=N_OBS)
knots_m = seq(0, 1, length.out = p_tilde+2)[-c(1,p_tilde+2)]
true_Bm_beta_time_vero = bs(time, knots = knots_m, intercept = T)


col_pat = rainbow(n)
yplot=c(-2,13)


x11()
par(mfrow = c(1,2))
plot(seq(0,1, len=N_OBS),telesca_data$y$y_mis[1:N_OBS,1], type='l', main='y misaligned',
     col=col_pat[1], xlab='time', ylab='y',ylim=yplot)
for(i in 2:n){
  points(seq(0,1, len=N_OBS),telesca_data$y$y_mis[1:N_OBS,i], type='l', col=col_pat[i])
}
points(seq(0,1, len=N_OBS),true_Bm_beta_time_vero%*%telesca_data$params$beta,col='black',lwd=5,type='l')
plot(seq(0,1, len=N_OBS),telesca_data$y$y_al[1:N_OBS,1], type='l', main='y aligned',
     col=col_pat[1], xlab='time', ylab='y',ylim=yplot)
for(i in 2:n){
  points(seq(0,1, len=N_OBS),telesca_data$y$y_al[1:N_OBS,i], type='l', col=col_pat[i])
}
points(seq(0,1, len=N_OBS),true_Bm_beta_time_vero%*%telesca_data$params$beta,col='black',lwd=5,type='l')


