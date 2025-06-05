
# PLOTS: 

load("Our_model_knee_data_result.RData")

n1 <- 32 # Group 1 (green)
n2 <- 32 # Group 2 (red)
n3 <- 36 # Group 3 (blue)

# Generate color vectors with shades
greens <- colorRampPalette(c("lightgreen", "darkgreen"))(n1)
reds <- colorRampPalette(c("pink", "darkred"))(n2)
blues <- colorRampPalette(c("lightblue", "darkblue"))(n3)

colors = list()
colors[[1]] = greens
colors[[2]] = reds 
colors[[3]] = blues
colors2 = c('green3','coral1','blue3')

n_obs_max  = dim(y)[1] 
n_obs=result$n_obs
n_patients = dim(y)[2] # y has patients in the columns
y_list = lapply(seq_len(ncol(y)), function(i) na.omit(y[, i]))
n_obs=sapply(y_list,length)

y_all = list()    
y_val  = list()    
n_obs  = list()

n0 = 1 
nF = n_per_group[1]

for(g in 1:n_groups){
  
  y_temp = y[,n0:nF]
  y_list_temp = lapply(seq_len(ncol(y_temp)), function(i) na.omit(y_temp[, i]))
  n_obs_temp=sapply(y_list_temp,length)
  y_temp = unlist(y_list_temp)
  y_temp = as.numeric(y_temp)
  
  y_all[[g]] = y_list_temp 
  y_val[[g]]  = y_temp 
  n_obs[[g]]  = n_obs_temp
  
  n0 = nF + 1 
  nF = nF + n_per_group[g+1]
}

time_grid = seq(0,1, length.out=1000)
knots_m = seq(0, 1, length.out = n_knots_m+2)[-c(1,n_knots_m+2)]
BM = bs(seq(0,1,len = 100),knots = knots_m, intercept = F)

ylim_plot=c(-5,100)

#-- Aligned curves control --# 
x11()
par(mfrow=c(1,2))
plot(seq(0,1, len=n_obs[[1]][1]),y_control[1:n_obs[[1]][1],1], type='l', main='y_control',
     col=greens[1], xlab='time', ylab='y',ylim=ylim_plot)
for(i in 2:n_per_group[1]){
  points(seq(0,1, len=n_obs[[1]][i]),y_control[1:n_obs[[1]][i],i], type='l', col=greens[i])
}
points(seq(0,1, len=100),BM%*%result$post$beta_group[[1]],col='green3',lwd=5,type='l')
plot(seq(0,1, len=n_obs[[1]][1]),result$y_star[[1]][[1]], type='l', main='',
     col=greens[1], xlab='', ylab='',ylim=ylim_plot)
for(i in 2:n_per_group[1]){
  points(seq(0,1, len=n_obs[[1]][i]),result$y_star[[1]][[i]], type='l', col=greens[i])
}
points(seq(0,1, len=100),BM%*%result$post$beta_group[[1]],col='green3',lwd=5,type='l')

#-- Aligned curves surgery --# 
x11()
par(mfrow=c(1,2))
plot(seq(0,1, len=n_obs[[2]][1]),y_surgery[1:n_obs[[2]][1],1], type='l', main='',
     col=reds[1], xlab='time', ylab='y',ylim=ylim_plot)
for(i in 2:n_per_group[2]){
  points(seq(0,1, len=n_obs[[2]][i]),y_surgery[1:n_obs[[2]][i],i], type='l', col=reds[i])
}
points(seq(0,1, len=100),BM%*%result$post$beta_group[[2]],col='coral1',lwd=5,type='l')
plot(seq(0,1, len=n_obs[[2]][1]),result$y_star[[2]][[1]], type='l', main='',
     col=reds[1], xlab='', ylab='',ylim=ylim_plot)
for(i in 2:n_per_group[2]){
  points(seq(0,1, len=n_obs[[2]][i]),result$y_star[[2]][[i]], type='l', col=reds[i])
}
points(seq(0,1, len=100),BM%*%result$post$beta_group[[2]],col='coral1',lwd=5,type='l')

#-- Aligned curves physiotherapy --# 
x11()
par(mfrow=c(1,2))
plot(seq(0,1, len=n_obs[[3]][1]),y_physio[1:n_obs[[3]][1],1], type='l', main='y_physiotherapy ',
     col=blues[1], xlab='time', ylab='y',ylim=ylim_plot)
for(i in 2:n_per_group[3]){
  points(seq(0,1, len=n_obs[[3]][i]),y_physio[1:n_obs[[3]][i],i], type='l', col=blues[i])
}
points(seq(0,1, len=100),BM%*%result$post$beta_group[[3]],col='blue3',lwd=5,type='l')
plot(seq(0,1, len=n_obs[[3]][1]),result$y_star[[3]][[1]], type='l', main='',
     col=blues[1], xlab='', ylab='',ylim=ylim_plot)
for(i in 2:n_per_group[3]){
  points(seq(0,1, len=n_obs[[3]][i]),result$y_star[[3]][[i]], type='l', col=blues[i])
}
points(seq(0,1, len=100),BM%*%result$post$beta_group[[3]],col='blue3',lwd=5,type='l')

# CREDIBILY BANDS

color_group = c("green3","red1","blue3")
color_bands = c("green4","red4","blue4")
color_region = c("lightgreen",'lightcoral','lightblue')

time_grid = seq(0,1, length.out=1000)
knots_m = seq(0, 1, length.out = n_knots_m+2)[-c(1,n_knots_m+2)]
BM_band = bs(seq(0,1,len = 1000),knots = knots_m, intercept = F)

alfa = 0.05 

mean_curve = list()
lower_bound = list()
upper_bound = list()

for(g in 1:n_groups){
  beta_samples = result$full$beta_group_s[[g]][(nburn+1):(niter+nburn),]
  posterior_curves <- beta_samples %*% t(BM_band)
  mean_curve[[g]] <- apply(posterior_curves, 2, mean)
  lower_bound[[g]] <- apply(posterior_curves, 2, quantile, probs = alfa/2)
  upper_bound[[g]] <- apply(posterior_curves, 2, quantile, probs = 1-alfa/2)
}

# Control 
x11()
plot(seq(0,1,len=length(mean_curve[[1]])), mean_curve[[1]], type = 'l', col = color_group[1], ylim = ylim_plot, 
     main = '', xlab = '', ylab = '')
polygon(c(time_grid, rev(time_grid)), c(upper_bound[[1]], rev(lower_bound[[1]])), col = color_region[1], border = NA)
lines(seq(0,1,len=length(mean_curve[[1]])), mean_curve[[1]], col = color_group[1], lwd = 2)
lines(time_grid, lower_bound[[1]], col = color_bands[1], lty = 2)
lines(time_grid, upper_bound[[1]], col = color_bands[1], lty = 2)

# Surgery 
x11()
plot(seq(0,1,len=length(mean_curve[[2]])), mean_curve[[2]], type = 'l', col = color_group[2], ylim = ylim_plot, 
     main = '', xlab = '', ylab = '')
polygon(c(time_grid, rev(time_grid)), c(upper_bound[[2]], rev(lower_bound[[2]])), col = color_region[2], border = NA)
lines(seq(0,1,len=length(mean_curve[[2]])), mean_curve[[2]], col = color_group[2], lwd = 2)
lines(time_grid, lower_bound[[2]], col = color_bands[2], lty = 2)
lines(time_grid, upper_bound[[2]], col = color_bands[2], lty = 2)

# Physiotherapy 
x11()
plot(seq(0,1,len=length(mean_curve[[3]])), mean_curve[[3]], type = 'l', col = color_group[3], ylim = ylim_plot, 
     main = '', xlab = '', ylab = '')
polygon(c(time_grid, rev(time_grid)), c(upper_bound[[3]], rev(lower_bound[[3]])), col = color_region[3], border = NA)
lines(seq(0,1,len=length(mean_curve[[3]])), mean_curve[[3]], col = color_group[3], lwd = 2)
lines(time_grid, lower_bound[[3]], col = color_bands[3], lty = 2)
lines(time_grid, upper_bound[[3]], col = color_bands[3], lty = 2)


# WARPING FUNCTIONS 

h_p = result$post$h

x11()
par(mfrow=c(1,n_groups))
for(g in 1:n_groups){
  plot(seq(0,1, len=n_obs[[g]][1]),h_p[[g]][[1]], type='l', main='',col=colors[[g]][1], xlab='', ylab='',ylim=c(0,1))
  for(i in 2:n_per_group[g]){
    points(seq(0,1, len=n_obs[[g]][i]),h_p[[g]][[i]], type='l', col=colors[[g]][i])
  }
  points(seq(0,1, len=n_obs[[g]][i]),seq(0,1, len=n_obs[[g]][i]),col='black',type='l',lwd=4)
}

graphics.off()

# Bands together with aligned curves 

# Control 
x11()
plot(seq(0,1, len=n_obs[[1]][1]),result$y_star[[1]][[1]], type='l', main='',
     col=greens[1], xlab='', ylab='',ylim=ylim_plot)
for(i in 2:n_per_group[1]){
  points(seq(0,1, len=n_obs[[1]][i]),result$y_star[[1]][[i]], type='l', col=greens[i])
}
points(seq(0,1, len=100),BM%*%result$post$beta_group[[1]],col='green3',lwd=5,type='l')
polygon(c(time_grid, rev(time_grid)), c(upper_bound[[1]], rev(lower_bound[[1]])), col = color_region[1], border = NA)
lines(seq(0,1,len=length(mean_curve[[1]])), mean_curve[[1]], col = color_group[1], lwd = 2)
lines(time_grid, lower_bound[[1]], col = color_bands[1], lty = 2)
lines(time_grid, upper_bound[[1]], col = color_bands[1], lty = 2)

# Surgery 
x11()
plot(seq(0,1, len=n_obs[[2]][1]),result$y_star[[2]][[1]], type='l', main='',
     col=reds[1], xlab='', ylab='',ylim=ylim_plot)
for(i in 2:n_per_group[2]){
  points(seq(0,1, len=n_obs[[2]][i]),result$y_star[[2]][[i]], type='l', col=reds[i])
}
points(seq(0,1, len=100),BM%*%result$post$beta_group[[2]],col='coral1',lwd=5,type='l')
polygon(c(time_grid, rev(time_grid)), c(upper_bound[[2]], rev(lower_bound[[2]])), col = color_region[2], border = NA)
lines(seq(0,1,len=length(mean_curve[[2]])), mean_curve[[2]], col = color_group[2], lwd = 2)
lines(time_grid, lower_bound[[2]], col = color_bands[2], lty = 2)
lines(time_grid, upper_bound[[2]], col = color_bands[2], lty = 2)

# Physiotherapy 
x11()
plot(seq(0,1, len=n_obs[[3]][1]),result$y_star[[3]][[1]], type='l', main='',
     col=blues[1], xlab='', ylab='',ylim=ylim_plot)
for(i in 2:n_per_group[3]){
  points(seq(0,1, len=n_obs[[3]][i]),result$y_star[[3]][[i]], type='l', col=blues[i])
}
points(seq(0,1, len=100),BM%*%result$post$beta_group[[3]],col='blue3',lwd=5,type='l')
polygon(c(time_grid, rev(time_grid)), c(upper_bound[[3]], rev(lower_bound[[3]])), col = color_region[3], border = NA)
lines(seq(0,1,len=length(mean_curve[[3]])), mean_curve[[3]], col = color_group[3], lwd = 2)
lines(time_grid, lower_bound[[3]], col = color_bands[3], lty = 2)
lines(time_grid, upper_bound[[3]], col = color_bands[3], lty = 2)



