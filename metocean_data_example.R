#loading all required functions and packages

source("master_functions.R")

# Tuning parameters -------------------------------------------------------

thresh_prob = 0.67 #Non-exceedance probability. This was selected through a basic sensitivity analysis

k = 25 #Number of knots for threshold and scale spline models

k_shape = 8 #Number of knots for shape spline model

pred_Q = seq(-2,2,length.out=201) #Angles at which to evaluate predicted quantities

bw = 50 #Bandwdith parameter for the circular density estimation

num_neigh = 500 #Number of neighbours for local windown estimation

norm_choice = "L1" #Selecting the L1 norm. Can also choose L2 if desired/required

# Data -----------------------------------------

data = readRDS(file="dataset1.rds") #load data stored in the "dataset1.rds" file

data = data[,2:1] #selecting only columns of interest

std_data = apply(data,2,function(x){(x-mean(x))/sd(x)}) #normalising the data
#this standardising ensures the mean will be at (0,0), as required for the SPAR model

#computing standard deviations from each variable. Required for evaluating equidensity contours
sds_data = apply(data,2,sd)

#computing means of each variable. Required for transforming back to original scale
mus_data = apply(data,2,mean)

# Fitting SPAR model ------------------------------------------------------------

#We fit the SPAR model smoothly with L1 coordinates
SI_smooth_fit = fit_SPAR_model(sample_data = std_data,norm_choice = "L1",thresh_prob = thresh_prob,k=k,k_shape = k_shape,pred_Q = pred_Q)

#We fit the SPAR model locally (via local windows) with L1 coordinates
SI_local_fit = fit_SPAR_empirical(sample_data = std_data,norm_choice = "L1",thresh_prob = thresh_prob,pred_Q = pred_Q)

#We estimate the angular density with L1 coordinates
SI_angular_density = SPAR_angular_density(sample_data = std_data,norm_choice = "L1",pred_Q = pred_Q,bw = bw)

#Define density levels for which to evaluate equidensity contours 
density_levels = 10^(-(3:6)) 

#Adjusting density levels to account for Jacobian of standardisation
density_levels = density_levels*prod(sds_data)

#Estimate equidensity contours for the desired levels
SI_equidensity_curves = sapply(density_levels,SPAR_equidensity_contours,norm_choice="L1",SPAR_GPD=SI_smooth_fit,SPAR_ang=SI_angular_density,simplify = F)

# Validating model fits ---------------------------------------------------

#Obtaining polar coordinates of dataset
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
  
}

#Setting plotting parameters
par(mfrow=c(1,1),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

#Computing the empirical histogram for angular density
hist(polar_data$Q, freq = FALSE,xlab="Q", main = "Angular density",sub="L1 coordinates",col=NULL,cex.lab=1.2, cex.axis=1.2,cex.main=1.5)

#Comparing estimated angular density function to empirical
lines(pred_Q,SI_angular_density,lwd=4,col=2)

#Comparing local and smooth estimates of threshold function
plot(pred_Q,SI_smooth_fit$pred_thresh,xlab="Q",ylab="R",main="Threshold",sub="L1 coordinates",typ="l",col=2,cex.lab=1.2, cex.axis=1.2,cex.main=1.5,lwd=4,ylim=range(SI_smooth_fit$pred_thresh,SI_local_fit$pred_thresh))
lines(pred_Q,SI_local_fit$pred_thresh,lwd=4,col="red")

#Comparing local and smooth estimates of scale function
plot(pred_Q,SI_smooth_fit$pred_para$scale,xlab="Q",ylab="R",main="Scale",sub="L1 coordinates",typ="l",col=3,cex.lab=1.2, cex.axis=1.2,cex.main=1.5,lwd=4,ylim=range(SI_smooth_fit$pred_para$scale,SI_local_fit$pred_para$scale))
lines(pred_Q,SI_local_fit$pred_para$scale,lwd=4,col="green")

#Comparing local and smooth estimates of shape function
plot(pred_Q,SI_smooth_fit$pred_para$shape,xlab="Q",ylab="R",main="Shape",sub="L1 coordinates",typ="l",col=4,cex.lab=1.2, cex.axis=1.2,cex.main=1.5,lwd=4,ylim=range(SI_smooth_fit$pred_para$shape,SI_local_fit$pred_para$shape))
lines(pred_Q,SI_local_fit$pred_para$shape,lwd=4,col="blue")

#To avoid the plots needing a very long time to render, we only plot a subset of 5000 points
set.seed(1)
rand_5000 = sample(1:dim(data)[1],5000,replace=F)

#We compute GPD upper bound for the conditional radial component
upper_bound = SI_smooth_fit$pred_thresh - SI_smooth_fit$pred_para$scale/SI_smooth_fit$pred_para$shape

#We check which coordinate system we are working in. This allows us to define vectors on the corresponding unit sphere
if(norm_choice == "L1"){
  
  #Defining points on the unit sphere for the L1 norm
  u_vec = ifelse(pred_Q>=0,(1-pred_Q),(pred_Q+1))
  v_vec = ifelse(pred_Q>=0, 1-abs(u_vec),-1+abs(u_vec))
  
} else {
  
  #Defining points on the unit sphere for the L2 norm
  u_vec = ifelse(pred_Q>=0,cos(pi*pred_Q/2),cos(-pi*pred_Q/2))
  v_vec = ifelse(pred_Q>=0, sqrt(1-u_vec^2),-sqrt(1-u_vec^2))
  
}

#Computing upper bound contour on Caresian coordinates 
upper_contour = cbind(upper_bound*u_vec,upper_bound*v_vec); upper_contour = rbind(upper_contour,upper_contour[1,])

#Computing threshold contour on Cartesian coordinates
thresh_contour = cbind(SI_smooth_fit$pred_thresh*u_vec,SI_smooth_fit$pred_thresh*v_vec); thresh_contour = rbind(thresh_contour,thresh_contour[1,])

#We transform everything from the normalised data back to the original coordinates
#This function is a wrapper for performing this transformation
normalisation_inverse_function = function(y){
  mu = y[1]
  sd = y[2] 
  norm_data = y[3:length(y)]
  return(norm_data*sd + mu)
}

#Transforming upper contour to original coordinates
upper_contour = apply(rbind(mus_data,sds_data,upper_contour),2,normalisation_inverse_function)

#Transforming threshold contour to original coordinates
thresh_contour = apply(rbind(mus_data,sds_data,thresh_contour),2,normalisation_inverse_function)

#Transforming density contours back to original coordinates
for(i in 1:length(SI_equidensity_curves)){
  SI_equidensity_curves[[i]] = apply(rbind(mus_data,sds_data,SI_equidensity_curves[[i]]),2,normalisation_inverse_function)
}

#Plotting sample of 5000 datapoints with upper bound
plot(data[rand_5000,],pch=16,col="grey",main="Upper bound",sub="Random 5000 points",xlab="Tz",ylab="Hs",cex.lab=1.2, cex.axis=1.2,cex.main=1.5,ylim=range(upper_contour[,2]),xlim=range(upper_contour[,1]))
lines(upper_contour,lwd=3,col=4)

#Plotting sample of 5000 datapoints with equidensity contours and threshold function
plot(data[rand_5000,],pch=16,col="grey",main="Density contours and threshold function",sub="Random 5000 points",xlab="Tz",ylab="Hs",cex.lab=1.2, cex.axis=1.2,cex.main=1.5,ylim=range(SI_equidensity_curves),xlim=range(SI_equidensity_curves))
for(i in 1:length(SI_equidensity_curves)){
  lines(SI_equidensity_curves[[i]],lwd=3,col=i+2)
}
lines(thresh_contour,lwd=3,col=2)
legend(range(SI_equidensity_curves)[1],range(SI_equidensity_curves)[2],legend=c("Threshold",paste0("10^(-",3:(length(SI_equidensity_curves)+2),")")),lwd=3,col=2:(length(SI_equidensity_curves)+2),cex=1.2,bg="white")

#Comparing local GPD fits
ref_Q = seq(-2,1.5,by=0.5)

#Computing local windows for which to evaluate GPD fits
windows_datasets = sapply(ref_Q,SPAR_empirical_windows,polar_data=polar_data,num_neigh=num_neigh,simplify = F)

par(mfrow=c(2,4),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

#This code computes QQ plots for each angle
#Taking the data in each window, we compare theoretical quantiles obtained using the SPAR model with empirical quantiles
#This assumes the data within each window are stationary, which appears reasonable given the small window size
for(i in 1:length(windows_datasets)){
  
  Q_index = which(pred_Q == ref_Q[i])
  
  window_R_data = windows_datasets[[i]]
  
  window_thresh = SI_smooth_fit$pred_thresh[Q_index]
  
  window_R_exc = window_R_data - window_thresh; window_R_exc = window_R_exc[window_R_exc>0]
  
  window_scale = SI_smooth_fit$pred_para$scale[Q_index]
  
  window_shape = SI_smooth_fit$pred_para$shape[Q_index]
  
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
