

omega_P <- function(dim){
  K <- Matrix::bandSparse(dim, k=-c(1), diag=list(rep(-1,dim)), symmetric=TRUE)
  diag(K) <- c(rep(2,dim-1),1)
  K <- matrix(as.numeric(K),nrow=dim, byrow=TRUE)
  return (K)
}

creation_dataset = function( N_groups, N_per_group, N_OBS,
                             p_tilde,k_tilde,q_tilde,
                             true_sigma_eps = 0.01, true_lambda = 5, true_sigma_gamma = 5,
                             RANDOM_NOISE_GROUP_MEAN,  RANDOM_NOISE_GROUP_VAR, q_cub=TRUE,q_lin=FALSE,q_quad = FALSE)
{
  p_true = p_tilde + 3
  n_knots_m = p_tilde
  k_true = k_tilde + 4 
  n_knots_gamma = k_tilde 
  if(q_lin){degree_q = 1}
  if(q_quad){degree_q = 2}
  if(q_cub){degree_q = 3}
  q_true = q_tilde + degree_q + 1 
  n_knots_h = q_tilde 
  
  #-- True sigma epsilon^2 --#
  true_eps_group = list()
  
  for(g in 1:N_groups){
    eps = matrix(NA,N_per_group[g],N_OBS)
    for(i in 1:N_per_group[g]){
      eps[i,] =  mvrnorm(n = N_OBS, mu = 0 , Sigma = true_sigma_eps )
    }
    true_eps_group[[g]] = eps
  }
  
  # True beta #
  Omega=omega_P(p_true)
  true_beta_group  = mvrnorm(n = N_groups, mu = rep(0,p_true) , Sigma = solve(Omega/true_lambda) )
  true_beta_group
  
  for (g in 2:N_groups){
    for(i in 1:p_true){
      true_beta_group[g,i] = true_beta_group[1,i] + rnorm(1, RANDOM_NOISE_GROUP_MEAN, RANDOM_NOISE_GROUP_VAR)
    }
  }
  
  # True gamma # 
  true_gamma_group = list()
  for(g in 1:N_groups){
    true_gamma = mvrnorm(n = N_per_group[g], mu = rep(0,k_true), Sigma = diag(k_true)*true_sigma_gamma )
    true_gamma_group[[g]] = true_gamma
  }
  
  # True phi
  
  true_phi_group = list()
  true_csi_group = list()
  for(g in 1:N_groups){
    true_csi_i = matrix(NA,q_true,10)
    true_phi_i = matrix(NA,q_true,10)
    for(i in 1:N_per_group[g]){
      b_g = sample(1:15, size = 1, replace = TRUE)
      A_gi = sample(5:20, size = 1, replace = TRUE)
      true_csi = rep(0,q_true)
      for(j in 2:q_true){
        true_csi[j] = rgamma(n=1,A_gi,b_g)
      }
      true_csi_i[,i] = true_csi
      true_phi_i[,i] = cumsum(true_csi) / sum(true_csi)
      
    }
    true_csi_group[[g]] = t(true_csi_i)
    true_phi_group[[g]] = t(true_phi_i)
  }
  
  #-- True h_i --#
  time = seq(0,1,length=N_OBS)
  knots_h = seq(0, 1, length.out = q_tilde+2)[-c(1,q_tilde+2)]
  true_Bh <- bs(time, knots = knots_h, intercept = T, degree = degree_q )
  true_h_i_group = list()
  for(g in 1:N_groups){
    h_i = true_Bh %*% t(true_phi_group[[g]])
    true_h_i_group[[g]] = h_i
  }
  
  #-- True Bm --#
  knots_m = seq(0, 1, length.out = n_knots_m+2)[-c(1,n_knots_m+2)]
  knots_gamma = seq(0, 1, length.out = n_knots_gamma+2)[-c(1,n_knots_gamma+2)]
  true_Bm_beta_group = list()
  for(g in 1:N_groups){
    true_Bm_beta = lapply(1:N_per_group[g], function(i){bs(true_h_i_group[[g]][,i], knots = knots_m, intercept = F)})
    true_Bm_beta_group[[g]] = true_Bm_beta 
  }
  true_Bm_gamma_group = list()
  for(g in 1:N_groups){
    true_Bm_gamma = lapply(1:N_per_group[g], function(i){bs(true_h_i_group[[g]][,i], knots = knots_gamma, intercept = T)})
    true_Bm_gamma_group[[g]] = true_Bm_gamma 
  }
  
  true_Bm_beta_time_vero = bs(time, knots = knots_m, intercept = F)
  true_Bm_gamma_time_vero = bs(time, knots = knots_gamma, intercept = T)
  
  #-- True y_i --#
  Y_group_list = list()
  Y_group_noise_free_list = list()
  Y_all_group_list = list()
  Y_all_group_noise_free_list = list()
  for(g in 1:N_groups){
    
    y = lapply(1:N_per_group[g],function(i){true_Bm_beta_group[[g]][[i]]%*%true_beta_group[g,] + 
        true_Bm_gamma_group[[g]][[i]]%*%true_gamma_group[[g]][i,] + true_eps_group[[g]][i,]})
    
    y_noise_free = lapply(1:N_per_group[g],function(i){true_Bm_beta_group[[g]][[i]]%*%true_beta_group[g,] + 
        true_Bm_gamma_group[[g]][[i]]%*%true_gamma_group[[g]][i,] })
    
    y_all = lapply(1:N_per_group[g],function(i){true_Bm_beta_time_vero%*%true_beta_group[g,] + 
        true_Bm_gamma_time_vero%*%true_gamma_group[[g]][i,] + true_eps_group[[g]][i,]})
    
    y_all_noise_free = lapply(1:N_per_group[g],function(i){true_Bm_beta_time_vero%*%true_beta_group[g,] + 
        true_Bm_gamma_time_vero%*%true_gamma_group[[g]][i,]})
    
    Y_group_list[[g]] = y
    Y_group_noise_free_list[[g]] = y_noise_free
    Y_all_group_list[[g]] = y_all
    Y_all_group_noise_free_list[[g]] = y_all_noise_free
    
  }
  
  Y_group = list()
  Y_group_noise_free = list()
  Y_all_group = list()
  Y_all_group_noise_free = list()
  for(g in 1:N_groups){
    y_true=matrix(NA,nrow=N_OBS,ncol=N_per_group[g])
    y_true_noise_free=matrix(NA,nrow=N_OBS,ncol=N_per_group[g])
    y_true_all=matrix(NA,nrow=N_OBS,ncol=N_per_group[g])
    y_true_all_noise_free=matrix(NA,nrow=N_OBS,ncol=N_per_group[g])
    for(i in 1:N_per_group[g]){
      y_true[1:N_OBS,i]=Y_group_list[[g]][[i]]
      y_true_noise_free[1:N_OBS,i]=Y_group_noise_free_list[[g]][[i]]
      y_true_all[1:N_OBS,i]=Y_all_group_list[[g]][[i]]
      y_true_all_noise_free[1:N_OBS,i]=Y_all_group_noise_free_list[[g]][[i]]
    }
    Y_group[[g]] = y_true
    Y_group_noise_free[[g]] = y_true_noise_free
    Y_all_group[[g]] = y_true_all
    Y_all_group_noise_free[[g]] = y_true_all_noise_free
  }
  
  y_to_use = do.call(cbind,Y_group)
  
  return(list( y = list(y_to_use = y_to_use,y_mis = Y_group, y_mis_noise_free = Y_group_noise_free,
                        y_al = Y_all_group, y_al_noise_free = Y_all_group_noise_free),
               params = list(beta = true_beta_group,           phi = true_phi_group, csi = true_csi_group, 
                             gamma = true_gamma_group, lambda = true_lambda, sigma_gamma = true_sigma_gamma, sigma_eps = true_sigma_eps)
  ))
}



N_groups = 2
N_per_group = rep(10,N_groups)
N_OBS = 300
p_tilde = 10
k_tilde = 3
q_tilde = 6
true_sigma_eps = 0.01
true_lambda = 5
true_sigma_gamma = 5
RANDOM_NOISE_GROUP_MEAN = 1
RANDOM_NOISE_GROUP_VAR = 1 

q_cub = TRUE
q_quad = FALSE
q_lin = FALSE


N_data = 3

dati_simulati = list()
seed = -2
set.seed(seed)
for(n in 1:N_data){
  dati_simulati[[n]] = creation_dataset(N_groups, N_per_group, N_OBS,
                                        p_tilde,k_tilde,q_tilde,
                                        true_sigma_eps , true_lambda , true_sigma_gamma ,
                                        RANDOM_NOISE_GROUP_MEAN,  RANDOM_NOISE_GROUP_VAR, q_cub, q_quad, q_lin)
}

dati_sim = dati_simulati[[1]]

Y_all_group = dati_sim$y$y_al_noise_free
Y_group = dati_sim$y$y_mis_noise_free
BETA = dati_sim$params$beta

time = seq(0,1,length=N_OBS)
knots_m = seq(0, 1, length.out = p_tilde+2)[-c(1,p_tilde+2)]
true_Bm_beta_time_vero = bs(time, knots = knots_m, intercept = F)

col_pat = rainbow(sum(N_per_group))
col_beta = c("green","red","blue")

greens <- colorRampPalette(c("lightgreen", "darkgreen"))(10)
reds <- colorRampPalette(c("pink", "darkred"))(10)
blues <- colorRampPalette(c("lightblue", "darkblue"))(10)

colors = list()
colors[[1]] = greens
colors[[2]] = reds 
colors[[3]] = blues

yplot=c(-10,20)

# PLOT THE CURVES # 
x11()
par(mfrow = c(1,2))
plot(seq(0,1, len=N_OBS),Y_all_group[[1]][1:N_OBS,1], type='l', main='y aligned',
     col=colors[[1]][1], xlab='time', ylab='y',ylim=yplot)
for(g in 1:N_groups){
  for(i in 1:N_per_group[g]){
    points(seq(0,1, len=N_OBS),Y_all_group[[g]][1:N_OBS,i], type='l', col=colors[[g]][i])
  }
  points(seq(0,1, len=N_OBS),true_Bm_beta_time_vero%*%prova$params$beta[g,],col=col_beta[g],lwd=5,type='l')
}
plot(seq(0,1, len=N_OBS),Y_group[[1]][1:N_OBS,1], type='l', main='y misaligned',
     col=colors[[1]][1], xlab='time', ylab='y',ylim=yplot)
for(g in 1:N_groups){
  for(i in 1:N_per_group[g]){
    points(seq(0,1, len=N_OBS),Y_group[[g]][1:N_OBS,i], type='l', col=colors[[g]][i])
  }
  points(seq(0,1, len=N_OBS),true_Bm_beta_time_vero%*%prova$params$beta[g,],col=col_beta[g],lwd=5,type='l')
}
