#loading all required functions

source("master_functions.R")

# Example data - bivariate normal ------------------------------------------------------------------

#We simulate 10,000 points on standard Laplace margins from a bivariate Gaussian copula with correlation coefficient 0.5

n = 10000

d = 2

rho = 0.5

normc = normalCopula(param = rho, dim = d)

set.seed(1)
example_data = apply(rCopula(n, copula = normc),2,Laplace_inverse)

#polar transformation
polar_data = rect2polar(t(example_data))

#altering the angular component to be a vector rather than a matrix
polar_data$phi = as.vector(polar_data$phi)

names(polar_data) = c("R","Phi")

polar_data = as.data.frame(polar_data)

# Tuning parameters -------------------------------------

#The below tuning parameters are recommended values for the SPAR inference. They can be adjusted to alter model fits

#Non-exceedances probability
thresh_prob = 0.8

#Basis dimension/number of spline knots. Higher = more flexibility
k = 25

#Grid of angular points
pred_phi = seq(0,2*pi,length.out=1001)

#Bandwidth for kernel density estimation. Higher = more flexibility
bw = 50

#Number of neighbours for local estimation approach
num_neigh = 300

# Fitting SPAR ------------------------------------------------------------

#The below code fits the SPAR model using both angular systems. The threshold and parameter functions are estimated both locally and smoothly
#The code also estimates isodensity contours for a range of density levels

smooth_fit = SPAR_smooth_polar(sample_data = example_data,thresh_prob = thresh_prob,k=k,pred_phi = pred_phi)

local_fit = SPAR_local_polar(sample_data = example_data,thresh_prob = thresh_prob,pred_phi = pred_phi,num_neigh = num_neigh)

angular_density = SPAR_angular_density_polar(sample_data = example_data,pred_phi = pred_phi,bw = bw)

# Obtaining statistics and simulations ----------------------------------------------

#density levels for which to evaluate isodensity contours
density_levels = 10^(-(3:8))  

#obtaining isodensity contours
isodensity_contours = SPAR_equidensity_contours_polar(density_levels = density_levels,SPAR_GPD=smooth_fit,SPAR_ang=angular_density)

#return period for evaluating return level set
ret_period = 10 

#number of data points per year (treating the simulations as daily data)
obs_year = 365

#Estimate return level set for desired return period
RL_set = SPAR_ret_level_sets_polar(ret_period = ret_period,obs_year = obs_year,SPAR_GPD=smooth_fit)

#Number of observations to simulate from the SPAR model
nsim = n

#obtaining model simulations
simulated_data = SPAR_simulation_polar(nsim=nsim,SPAR_GPD=smooth_fit,SPAR_ang=angular_density)

# Comparing estimates and validation -----------------------------------------------------

pdf(file="plots/gaussian_local_smooth_polar.pdf",width=12,height=4)
#Setting plotting parameters
par(mfrow=c(1,3),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

#Comparing local and smooth estimates of threshold function
plot(pred_phi,smooth_fit$pred_thresh,xlab=expression(phi),ylab=expression(u[gamma]),main="Threshold",typ="l",col=2,cex.lab=1.2, cex.axis=1.2,cex.main=1.55,cex.sub=1.3,lwd=4,ylim=range(smooth_fit$pred_thresh,local_fit$pred_thresh))
lines(pred_phi,local_fit$pred_thresh,lwd=4,col="red")

#Comparing local and smooth estimates of scale function
plot(pred_phi,smooth_fit$pred_para$scale,xlab=expression(phi),ylab=expression(tau),main="Scale",typ="l",col=3,cex.lab=1.2, cex.axis=1.2,cex.main=1.55,cex.sub=1.3,lwd=4,ylim=range(smooth_fit$pred_para$scale,local_fit$pred_para$scale))
lines(pred_phi,local_fit$pred_para$scale,lwd=4,col="green")

#Comparing local and smooth estimates of shape function
plot(pred_phi,smooth_fit$pred_para$shape,xlab=expression(phi),ylab=expression(xi),main="Shape",typ="l",col=4,cex.lab=1.2, cex.axis=1.2,cex.main=1.55,cex.sub=1.3,lwd=4,ylim=range(smooth_fit$pred_para$shape,local_fit$pred_para$shape))
lines(pred_phi,local_fit$pred_para$shape,lwd=4,col="blue")

dev.off()

#Load in data files containing true angular density
true_angdens = readRDS(file="datafiles/gaussian_L2.rds")

true_angdens = rbind(true_angdens[true_angdens$theta >0,],true_angdens[true_angdens$theta < 0,])

true_angdens$theta[true_angdens$theta<0] = 2*pi + true_angdens$theta[true_angdens$theta<0]

pdf(file="plots/gaussian_ang_dens_polar.pdf",width=6,height=6)

#Setting plotting parameters
par(mfrow=c(1,1),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

#Plots comparing angular density estimates to truth for both sets of estimates
plot(true_angdens$theta,true_angdens$ft,type="l",lwd=3,col=2,xlab=expression(phi),ylab=bquote(f[Phi](phi)),main="Angular density",ylim=range(angular_density,true_angdens$ft),cex.lab=1.2, cex.axis=1.2,cex.main=1.55,cex.sub=1.3)
lines(pred_phi,angular_density,lwd=3,col=3)
legend(0,0.1,legend=c("True","Estimated"),lwd=3,col=c(2,3),cex=1.2,bg="white")

dev.off()

pdf(file="plots/gaussian_ang_dens_diag_polar.pdf",width=6,height=6)
#Setting plotting parameters
par(mfrow=c(1,1),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

#Extracting empirical density function data
hist_data = hist(polar_data$Phi, plot=F) 
#Computing the empirical histogram for angular density
hist(polar_data$Phi, freq = FALSE,xlab=expression(Phi),ylab=expression(f[Phi](phi)), main = "Angular density",col=NULL,cex.lab=1.2, cex.axis=1.2,cex.main=1.5,ylim=range(angular_density,0,hist_data$density))

#Comparing estimated angular density function to empirical
lines(pred_phi,angular_density,lwd=4,col="blue")

dev.off()

#Computing true isodensity contours for Gaussian copula. Isodensity contours will be the same for both coordinate systems
gaussian_root_wrapper = function(phi,val,normc){
  gaussian_root = function(r,phi,val,normc){
    u <- cos(phi)
    v <- sin(phi)
    vec = c(r*u,r*v)
    return(dCopula(u=Laplace_cdf(vec),copula=normc,log=F)*(exp(-abs(vec[1])-abs(vec[2])))/4 -val)
  }
  return(uniroot(gaussian_root,interval = c(0,150),phi=phi,val=val,normc=normc)$root)
}

#Define points in L1 unit circle
u_vec = cos(pred_phi)
v_vec = sin(pred_phi)

#Compute isodensity contours for 3 density levels, store in a list
gaussian_density_contours = list()
for(i in 1:length(density_levels)){
  
  gaussian_radii = sapply(pred_phi, gaussian_root_wrapper,val=density_levels[i],normc=normc)
  
  gaussian_density_contours[[i]] = cbind(gaussian_radii*u_vec,gaussian_radii*v_vec); gaussian_density_contours[[i]] = rbind(gaussian_density_contours[[i]],gaussian_density_contours[[i]][1,])
  
}

#Colours for plotting
colfunc_true = colorRampPalette(c( "red","blue"))
colfunc_est = colorRampPalette(c( "red","blue"))

pdf(file="plots/gaussian_isodensity_contours_polar.pdf",width=6,height=6)

#Setting plotting parameters
par(mfrow=c(1,1),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

#Plots comparing estimated and true isodensity contours
plot(gaussian_density_contours[[1]],xlab="X",ylab="Y",main="Isodensity contours for Gaussian copula",type="l",col=1,lwd=3,xlim=range(gaussian_density_contours),ylim=range(gaussian_density_contours),cex.lab=1.2, cex.axis=1.2,cex.main=1.55,cex.sub=1.3)
for(i in 1:length(gaussian_density_contours)){
  lines(gaussian_density_contours[[i]],lwd=3,col=colfunc_true(length(gaussian_density_contours))[i])
  lines(isodensity_contours[[i]],lwd=3,col=colfunc_est(length(gaussian_density_contours))[i],lty=2)
}
legend(-15,15,legend=paste0("10^(",-(1:length(gaussian_density_contours))-2,")"),lwd=3,col=colfunc_true(length(gaussian_density_contours)),cex=1.2,bg="white")

dev.off()

pdf(file="plots/gaussian_rl_set_polar.pdf",width=6,height=6)

#Setting plotting parameters
par(mfrow=c(1,1),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

#Plotting return level sets 
plot(example_data,xlab="X",ylab="Y",main="Return level set",col="grey",pch=16,lwd=3,xlim=c(-10,10),ylim=c(-10,10),cex.lab=1.2, cex.axis=1.2,cex.main=1.55,cex.sub=1.3)
lines(RL_set,lwd=3,col="purple")

dev.off()

pdf(file="plots/gaussian_simulated_data_polar.pdf",width=6,height=6)

#Setting plotting parameters
par(mfrow=c(1,1),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

#Plotting simulated data over original sample with L1 breakdown 
plot(example_data,xlab="X",ylab="Y",main="SPAR model simulations",col="grey",pch=16,lwd=3,xlim=range(example_data,simulated_data$data_sample),ylim=range(example_data,simulated_data$data_sample),cex.lab=1.2, cex.axis=1.2,cex.main=1.55,cex.sub=1.3)
points(simulated_data$data_sample,pch=16,col=adjustcolor(3,alpha.f = 0.2))
legend(range(example_data,simulated_data$data_sample)[1],range(example_data,simulated_data$data_sample)[2],legend=c("Observerd","Simulated"),pch=16,col=c("grey",adjustcolor(3,alpha.f = 0.2)),cex=1.2,bg="white")

dev.off()
