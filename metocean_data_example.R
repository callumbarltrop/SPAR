#loading all required functions and packages

source("master_functions.R")

# Tuning parameters -------------------------------------------------------

thresh_prob = 0.7 #Non-exceedance probability. 

k = 35 #Number of knots for threshold and scale spline models

k_shape = 12 #Number of knots for shape spline model

pred_Q = seq(-2,2,length.out=1001) #Angles at which to evaluate predicted quantities

pred_Q_loc = seq(-2,2,length.out=201) #Angles at which to evaluate predicted quantities for local fit

bw = 50 #Bandwdith parameter for the circular density estimation

num_neigh = 500 #Number of neighbours for local windown estimation

norm_choice = "L1" #Selecting the L1 norm. Can also choose L2 if desired/required

ret_period = 10 #return period for evaluating return level set

# Data -----------------------------------------

which.dataset = "A" #Select either A, B or C

data = readRDS(file=paste0("datafiles/dataset",which.dataset,".rds")) #loading in data from the datafiles folder

data = data[,2:1] #selecting only columns of interest - Tz and Hs

std_data = apply(data,2,function(x){(x-mean(x))/sd(x)}) #normalising the data
#this standardising ensures the mean will be at (0,0), as required for the SPAR model

#computing standard deviations from each variable. Required for transforming back to the original scale
sds_data = apply(data,2,sd)

#computing means of each variable. Required for transforming back to the original scale
mus_data = apply(data,2,mean)

#Frequency of observation. We observe hourly observations 365 days per year
obs_year = 365*24 

#Number of points to simulate from model
nsim = dim(data)[1]

# Fitting SPAR model ------------------------------------------------------------

#We fit the SPAR model smoothly with L1 coordinates
SPAR_smooth_fit = SPAR_smooth(sample_data = std_data,norm_choice = norm_choice,thresh_prob = thresh_prob,k=k,k_shape = k_shape,pred_Q = pred_Q)

#We fit the SPAR model locally (via local windows) with L1 coordinates
SPAR_local_fit = SPAR_local(sample_data = std_data,norm_choice = norm_choice,thresh_prob = thresh_prob,pred_Q = pred_Q_loc,num_neigh = num_neigh)

#We estimate the angular density with L1 coordinates
SPAR_angular_density = SPAR_angular_density(sample_data = std_data,norm_choice = norm_choice,pred_Q = pred_Q,bw = bw)

#Define density levels for which to evaluate equidensity contours 
density_levels = 10^(-c(3,6)) 

#Adjusting density levels to account for Jacobian of standardisation
density_levels = density_levels*prod(sds_data)

if(norm_choice == "L2"){
  density_levels = (pi/2)*density_levels
}

#Estimate equidensity contours for the desired levels
SPAR_equidensity_curves = SPAR_equidensity_contours(density_levels = density_levels,norm_choice=norm_choice,SPAR_GPD=SPAR_smooth_fit,SPAR_ang=SPAR_angular_density)

#Estimate return level set for desired return period
SPAR_RL_set = SPAR_ret_level_sets(ret_period = ret_period,obs_year = obs_year,norm_choice = norm_choice,SPAR_GPD = SPAR_smooth_fit)

#Simulate new data from fitted SPAR model
SPAR_simulated_data = SPAR_simulation(sample_data=std_data,nsim=nsim,norm_choice = norm_choice,thresh_prob = thresh_prob,k=k,k_shape = k_shape,pred_Q = pred_Q,bw=bw)

# Validating model fits ---------------------------------------------------

#Obtaining polar coordinates of dataset, along with points on the norm unit circle
if(norm_choice == "L1"){
  
  L1_rad = function(vec){
    return(abs(vec[1])+abs(vec[2]))
  }
  
  L1_ang = function(vec){
    return(sgn(vec[2])*(1-vec[1]))
  }
  
  #observed radial data
  R = apply(std_data,1,L1_rad)
  
  #observed angular data
  Q = apply( apply(std_data,2,function(x,r){return(x/r)},r=R),1,L1_ang)
  
  #dataframe of angular-radial data
  polar_data = data.frame(R=R,Q=Q)
  
  #Defining points on the unit sphere for the L1 norm
  u_vec = ifelse(pred_Q>=0,(1-pred_Q),(pred_Q+1))
  v_vec = ifelse(pred_Q>=0, 1-abs(u_vec),-1+abs(u_vec))
  
} else {
  
  L2_rad = function(vec){
    return(sqrt(vec[1]^2+vec[2]^2))
  }
  
  L2_ang = function(vec){
    return(sgn(vec[2])*(2/pi)*acos(vec[1]))
  }
  
  #observed radial data
  R = apply(std_data,1,L2_rad)
  
  #observed angular data
  Q = apply( apply(std_data,2,function(x,r){return(x/r)},r=R),1,L2_ang)
  
  #dataframe of angular-radial data
  polar_data = data.frame(R=R,Q=Q)
  
  #Defining points on the unit sphere for the L2 norm
  u_vec = ifelse(pred_Q>=0,cos(pi*pred_Q/2),cos(-pi*pred_Q/2))
  v_vec = ifelse(pred_Q>=0, sqrt(1-u_vec^2),-sqrt(1-u_vec^2))
  
}

#We transform everything from the normalised data back to the original coordinates
#This function is a wrapper for performing this transformation
normalisation_inverse_function = function(y){
  mu = y[1]
  sd = y[2] 
  norm_data = y[3:length(y)]
  return(norm_data*sd + mu)
}

#Transforming density contours back to original coordinates
for(i in 1:length(SPAR_equidensity_curves)){
  SPAR_equidensity_curves[[i]] = apply(rbind(mus_data,sds_data,SPAR_equidensity_curves[[i]]),2,normalisation_inverse_function)
}

SPAR_RL_set = apply(rbind(mus_data,sds_data,SPAR_RL_set),2,normalisation_inverse_function)

SPAR_simulated_data$data_sample = apply(rbind(mus_data,sds_data,SPAR_simulated_data$data_sample),2,normalisation_inverse_function)

pdf(file="plots/metocean_ang_dens_diag.pdf",width=6,height=6)
#Setting plotting parameters
par(mfrow=c(1,1),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

#Extracting empirical density function data
hist_data = hist(polar_data$Q, plot=F) 
#Computing the empirical histogram for angular density
hist(polar_data$Q, freq = FALSE,xlab="Q",ylab=expression(f[Q](q)), main = "Angular density",sub="L1 coordinates",col=NULL,cex.lab=1.2, cex.axis=1.2,cex.main=1.5,ylim=range(SPAR_angular_density,0,hist_data$density))

#Comparing estimated angular density function to empirical
lines(pred_Q,SPAR_angular_density,lwd=4,col="blue")

dev.off()

pdf(file="plots/metocean_local_smooth.pdf",width=12,height=4)
#Setting plotting parameters
par(mfrow=c(1,3),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

#Comparing local and smooth estimates of threshold function
plot(pred_Q,SPAR_smooth_fit$pred_thresh,xlab="Q",ylab=expression(u[gamma]),main="Threshold",sub="L1 coordinates",typ="l",col=2,cex.lab=1.2, cex.axis=1.2,cex.main=1.5,lwd=4,ylim=range(SPAR_smooth_fit$pred_thresh,SPAR_local_fit$pred_thresh))
lines(pred_Q_loc,SPAR_local_fit$pred_thresh,lwd=4,col="red")

#Comparing local and smooth estimates of scale function
plot(pred_Q,SPAR_smooth_fit$pred_para$scale,xlab="Q",ylab=expression(tau),main="Scale",sub="L1 coordinates",typ="l",col=3,cex.lab=1.2, cex.axis=1.2,cex.main=1.5,lwd=4,ylim=range(SPAR_smooth_fit$pred_para$scale,SPAR_local_fit$pred_para$scale))
lines(pred_Q_loc,SPAR_local_fit$pred_para$scale,lwd=4,col="green")

#Comparing local and smooth estimates of shape function
plot(pred_Q,SPAR_smooth_fit$pred_para$shape,xlab="Q",ylab=expression(xi),main="Shape",sub="L1 coordinates",typ="l",col=4,cex.lab=1.2, cex.axis=1.2,cex.main=1.5,lwd=4,ylim=range(SPAR_smooth_fit$pred_para$shape,SPAR_local_fit$pred_para$shape))
lines(pred_Q_loc,SPAR_local_fit$pred_para$shape,lwd=4,col="blue")

dev.off()

#Plotting estimated equidensity contours
png(file="plots/metocean_equidensity_contours.png",width=600,height=600)
par(mfrow=c(1,1),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

plot(data,pch=16,col="grey",main="Equidensity contours",sub="L1 coordinates",xlab="Tz",ylab="Hs",cex.lab=1.2, cex.axis=1.2,cex.main=1.5,ylim=range(SPAR_equidensity_curves),xlim=range(SPAR_equidensity_curves))
lines(SPAR_equidensity_curves[[1]],lwd=3,col="orange")
lines(SPAR_equidensity_curves[[2]],lwd=3,col="cyan")
legend(range(SPAR_equidensity_curves)[1],range(SPAR_equidensity_curves)[2],legend=c(expression(paste("10"^"-3")),expression(paste("10"^"-6"))),lwd=3,col=c("Orange","Cyan"),cex=1.2,bg="white")

dev.off()

#Plotting estimated return level set
png(file="plots/metocean_ret_level_set.png",width=600,height=600)
par(mfrow=c(1,1),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

plot(data,pch=16,col="grey",main="Return level set",sub="L1 coordinates",xlab="Tz",ylab="Hs",cex.lab=1.2, cex.axis=1.2,cex.main=1.5,ylim=range(SPAR_RL_set),xlim=range(SPAR_RL_set))
lines(SPAR_RL_set,lwd=3,col="purple")

dev.off()

#Plotting SPAR model simulations
png(file="plots/metocean_simulated_data.png",width=600,height=600)
par(mfrow=c(1,1),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

plot(data,pch=16,col="grey",main="SPAR simulations",sub="L1 coordinates",xlab="Tz",ylab="Hs",cex.lab=1.2, cex.axis=1.2,cex.main=1.5,ylim=range(data,SPAR_simulated_data$data_sample),xlim=range(data,SPAR_simulated_data$data_sample))
points(SPAR_simulated_data$data_sample,pch=16,col=adjustcolor(3,alpha.f = 0.2))
legend(range(data,SPAR_simulated_data$data_sample)[1],range(data,SPAR_simulated_data$data_sample)[2],legend=c("Observerd","Simulated"),pch=16,col=c("grey",adjustcolor(3,alpha.f = 0.2)),cex=1.2,bg="white")

dev.off()

#Extracting exceedance values on exponential scale
emp_quants = sort(SPAR_smooth_fit$exp_quants)

#Computing length of vector
m = length(emp_quants)

#Computing theoretical exponential quantiles
model_quants = qexp((1:m)/(m+1))

#Bootstrapping the observed exceedances
nboots = 500
emp_quants_boots = matrix(nrow=nboots,ncol=m)
for(j in 1:nboots){
  exceedance_sample = sort(emp_quants[sample(1:m,m,replace=T)])
  emp_quants_boots[j,] = exceedance_sample  
}

#Computing confidence intervals
upper_quants = apply(emp_quants_boots,2,quantile,probs=0.975)
lower_quants = apply(emp_quants_boots,2,quantile,probs=0.025)

#Plotting overall diagnostic
png(file="plots/metocean_global_diagnostic.png",width=500,height=500)
par(mfrow=c(1,1),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

plot(model_quants,emp_quants,xlim=range(model_quants,emp_quants,upper_quants,lower_quants),ylim=range(model_quants,emp_quants,upper_quants,lower_quants),pch=16,col="grey",xlab="Fitted",ylab="Empirical",main="QQ plot on exponential scale",cex.lab=1.2, cex.axis=1.2,cex.main=1.5)
polygon(c(rev(model_quants), model_quants),c(rev(upper_quants),lower_quants),col = 'grey80', border = NA)
points(model_quants,emp_quants,pch=16,col="black")
abline(a=0,b=1,lwd=3,col=2)

dev.off()

#Comparing local GPD fits
ref_Q = seq(-1.5,2,by=0.5)

#Computing local windows for which to evaluate GPD fits
windows_datasets = sapply(ref_Q,SPAR_empirical_windows,polar_data=polar_data,num_neigh=num_neigh,simplify = F)

pdf(file="plots/metocean_local_diagnostic.pdf",width=16,height=8)

par(mfrow=c(2,4),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

#This code computes QQ plots for each angle
#Taking the data in each window, we compare theoretical quantiles obtained using the SPAR model with empirical quantiles
#This assumes the data within each window are stationary, which appears reasonable given the small window size
for(i in 1:length(windows_datasets)){
  
  Q_index = which(pred_Q == ref_Q[i])
  
  window_R_data = windows_datasets[[i]]
  
  window_thresh = SPAR_smooth_fit$pred_thresh[Q_index]
  
  window_R_exc = window_R_data - window_thresh; window_R_exc = window_R_exc[window_R_exc>0]
  
  window_scale = SPAR_smooth_fit$pred_para$scale[Q_index]
  
  window_shape = SPAR_smooth_fit$pred_para$shape[Q_index]
  
  emp_quants = sort(window_R_exc)
  
  m = length(emp_quants)
  
  model_quants = (window_scale/window_shape)*( (1- ((1:m)/(m+1)) )^(-window_shape) - 1 )
  
  nboots = 500
  emp_quants_boots = matrix(nrow=nboots,ncol=m)
  for(j in 1:nboots){
    exceedance_sample = sort(emp_quants[sample(1:m,m,replace=T)])
    emp_quants_boots[j,] = exceedance_sample  
  }
  
  upper_quants = apply(emp_quants_boots,2,quantile,probs=0.975)
  lower_quants = apply(emp_quants_boots,2,quantile,probs=0.025)
  
  plot(model_quants,emp_quants,xlim=range(model_quants,emp_quants,upper_quants,lower_quants),ylim=range(model_quants,emp_quants,upper_quants,lower_quants),pch=16,col="grey",xlab="Fitted",ylab="Empirical",main=paste0("QQ plot for Q = ",ref_Q[i]),cex.lab=1.2, cex.axis=1.2,cex.main=1.5)
  polygon(c(rev(model_quants), model_quants),c(rev(upper_quants),lower_quants),col = 'grey80', border = NA)
  points(model_quants,emp_quants,pch=16,col="black",xlab="Model",ylab="Empirical",main="QQ plot")
  abline(a=0,b=1,lwd=3,col=2)
  
}

dev.off()