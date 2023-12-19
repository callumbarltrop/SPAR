#loading all required functions

source("master_functions.R")

# Example data - bivariate normal ------------------------------------------------------------------

#We simulatte 10,000 points on standard Laplace margins from a bivariate gaussian copula with correlation coefficient 0.5

n = 10000

d = 2

rho = 0.5

normc = normalCopula(param = rho, dim = d)

set.seed(1)
example_data = apply(rCopula(n, copula = normc),2,Laplace_inverse)


# Tuning parameters -------------------------------------

#The below tuning parameters are recommended values for the SPAR inference. They can be adjusted to alter model fits

#Non-exceedances probability
thresh_prob = 0.8

#Basis dimension/number of spline knots. Higher = more flexibility
k = 25

#Grid of angular points
pred_Q = seq(-2,2,length.out=1001)

#Bandwidth for kernel density estimation. Higher = more flexibility
bw = 50

#Number of neighbours for local estimation approach
num_neigh = 300

#Number of obervsations to simulate from the SPAR model
nsim = 5000

# Fitting SPAR ------------------------------------------------------------

#The below code fits the SPAR model using both angular systems. The threshold and parameter functions are estimated both locally and smoothly
#The code also estimates equidensity contours for a range of density levels

L1_smooth_fit = SPAR_smooth(sample_data = example_data,norm_choice = "L1",thresh_prob = thresh_prob,k=k,pred_Q = pred_Q)

L2_smooth_fit = SPAR_smooth(sample_data = example_data,norm_choice = "L2",thresh_prob = thresh_prob,k=k,pred_Q = pred_Q)

L1_local_fit = SPAR_local(sample_data = example_data,norm_choice = "L1",thresh_prob = thresh_prob,pred_Q = pred_Q,num_neigh = num_neigh)

L2_local_fit = SPAR_local(sample_data = example_data,norm_choice = "L2",thresh_prob = thresh_prob,pred_Q = pred_Q,num_neigh = num_neigh)

L1_angular_density = SPAR_angular_density(sample_data = example_data,norm_choice = "L1",pred_Q = pred_Q,bw = bw)

L2_angular_density = SPAR_angular_density(sample_data = example_data,norm_choice = "L2",pred_Q = pred_Q,bw = bw)

density_levels = 10^(-(3:6)) #density levels for which to evaluate equidensity contours 

L1_equidensity_density_curves = SPAR_equidensity_contours(density_levels = density_levels,norm_choice="L1",SPAR_GPD=L1_smooth_fit,SPAR_ang=L1_angular_density)

L2_equidensity_density_curves = SPAR_equidensity_contours(density_levels = density_levels,norm_choice="L2",SPAR_GPD=L2_smooth_fit,SPAR_ang=L2_angular_density)

L1_simulation = SPAR_simulation(sample_data=example_data,nsim=nsim,norm_choice = "L1",thresh_prob = thresh_prob,k=k,pred_Q = pred_Q,bw=bw)

L2_simulation = SPAR_simulation(sample_data=example_data,nsim=nsim,norm_choice = "L2",thresh_prob = thresh_prob,k=k,pred_Q = pred_Q,bw=bw)

# Comparing estimates and validation -----------------------------------------------------

#Setting plotting parameters
par(mfrow=c(1,3),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

#Comparing local and smooth estimates of threshold function
plot(pred_Q,L1_smooth_fit$pred_thresh,xlab="Q",ylab=expression(u[gamma]),main="Threshold",sub="L1 coordinates",typ="l",col=2,cex.lab=1.2, cex.axis=1.2,cex.main=1.5,lwd=4,ylim=range(L1_smooth_fit$pred_thresh,L1_local_fit$pred_thresh))
lines(pred_Q,L1_local_fit$pred_thresh,lwd=4,col="red")

#Comparing local and smooth estimates of scale function
plot(pred_Q,L1_smooth_fit$pred_para$scale,xlab="Q",ylab=expression(tau),main="Scale",sub="L1 coordinates",typ="l",col=3,cex.lab=1.2, cex.axis=1.2,cex.main=1.5,lwd=4,ylim=range(L1_smooth_fit$pred_para$scale,L1_local_fit$pred_para$scale))
lines(pred_Q,L1_local_fit$pred_para$scale,lwd=4,col="green")

#Comparing local and smooth estimates of shape function
plot(pred_Q,L1_smooth_fit$pred_para$shape,xlab="Q",ylab=expression(xi),main="Shape",sub="L1 coordinates",typ="l",col=4,cex.lab=1.2, cex.axis=1.2,cex.main=1.5,lwd=4,ylim=range(L1_smooth_fit$pred_para$shape,L1_local_fit$pred_para$shape))
lines(pred_Q,L1_local_fit$pred_para$shape,lwd=4,col="blue")

#Comparing local and smooth estimates of threshold function
plot(pred_Q,L2_smooth_fit$pred_thresh,xlab="Q",ylab=expression(u[gamma]),main="Threshold",sub="L2 coordinates",typ="l",col=2,cex.lab=1.2, cex.axis=1.2,cex.main=1.5,lwd=4,ylim=range(L2_smooth_fit$pred_thresh,L2_local_fit$pred_thresh))
lines(pred_Q,L2_local_fit$pred_thresh,lwd=4,col="red")

#Comparing local and smooth estimates of scale function
plot(pred_Q,L2_smooth_fit$pred_para$scale,xlab="Q",ylab=expression(tau),main="Scale",sub="L2 coordinates",typ="l",col=3,cex.lab=1.2, cex.axis=1.2,cex.main=1.5,lwd=4,ylim=range(L2_smooth_fit$pred_para$scale,L2_local_fit$pred_para$scale))
lines(pred_Q,L2_local_fit$pred_para$scale,lwd=4,col="green")

#Comparing local and smooth estimates of shape function
plot(pred_Q,L2_smooth_fit$pred_para$shape,xlab="Q",ylab=expression(xi),main="Shape",sub="L2 coordinates",typ="l",col=4,cex.lab=1.2, cex.axis=1.2,cex.main=1.5,lwd=4,ylim=range(L2_smooth_fit$pred_para$shape,L2_local_fit$pred_para$shape))
lines(pred_Q,L2_local_fit$pred_para$shape,lwd=4,col="blue")

#Computing true equidensity contours for gaussian copula. Equidensity contours will be the same for both coordinate systems
gaussian_root_wrapper = function(q,val,normc){
  gaussian_root = function(r,q,val,normc){
    u <- ifelse(q>=0,(1-q),(q+1))
    v <- ifelse(q>=0, 1-abs(u),-1+abs(u))
    vec = c(r*u,r*v)
    return(dCopula(u=Laplace_cdf(vec),copula=normc,log=F)*(exp(-abs(vec[1])-abs(vec[2])))/4 -val)
  }
  return(uniroot(gaussian_root,interval = c(0,150),q=q,val=val,normc=normc)$root)
}

#Define points in L1 unit circle
u_vec = ifelse(pred_Q>=0,(1-pred_Q),(pred_Q+1))
v_vec = ifelse(pred_Q>=0, 1-abs(u_vec),-1+abs(u_vec))

#Compute equidensity contours for 3 density levels, store in a list
gaussian_density_contours = list()
for(i in 1:length(density_levels)){
  
  gaussian_radii = sapply(pred_Q, gaussian_root_wrapper,val=density_levels[i],normc=normc)
  
  gaussian_density_contours[[i]] = cbind(gaussian_radii*u_vec,gaussian_radii*v_vec); gaussian_density_contours[[i]] = rbind(gaussian_density_contours[[i]],gaussian_density_contours[[i]][1,])
  
}

#Colours for plotting
colfunc_true = colorRampPalette(c( "red","blue"))
colfunc_est = colorRampPalette(c( "orange","cyan"))

#Setting plotting parameters
par(mfrow=c(1,2),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

#Plots comparing estimated and true equidensity contours
plot(gaussian_density_contours[[1]],xlab="X",ylab="Y",main="Equidensity contours for gaussian copula",sub="L1 coordinates",type="l",col=1,lwd=3,xlim=range(gaussian_density_contours),ylim=range(gaussian_density_contours),cex.lab=1.2, cex.axis=1.2,cex.main=1.5)
for(i in 1:length(gaussian_density_contours)){
  lines(gaussian_density_contours[[i]],lwd=3,col=colfunc_true(length(gaussian_density_contours))[i])
  lines(L1_equidensity_density_curves[[i]],lwd=3,col=colfunc_est(length(gaussian_density_contours))[i])
}
legend(-10,10,legend=paste0("10^(",-(1:length(gaussian_density_contours))-2,")"),lwd=3,col=colfunc_true(length(gaussian_density_contours)),cex=1.2,bg="white")

plot(gaussian_density_contours[[1]],xlab="X",ylab="Y",main="Equidensity contours for gaussian copula",sub="L2 coordinates",type="l",col=1,lwd=3,xlim=range(gaussian_density_contours),ylim=range(gaussian_density_contours),cex.lab=1.2, cex.axis=1.2,cex.main=1.5)
for(i in 1:length(gaussian_density_contours)){
  lines(gaussian_density_contours[[i]],lwd=3,col=colfunc_true(length(gaussian_density_contours))[i])
  lines(L2_equidensity_density_curves[[i]],lwd=3,col=colfunc_est(length(gaussian_density_contours))[i])
}
legend(-10,10,legend=paste0("10^(",-(1:length(gaussian_density_contours))-2,")"),lwd=3,col=colfunc_true(length(gaussian_density_contours)),cex=1.2,bg="white")

#Load in data files containing true angular densities for both coordinate systems
true_angdens_L1 = readRDS(file="datafiles/gaussian_L1.rds")
true_angdens_L2 = readRDS(file="datafiles/gaussian_L2.rds")

#Transforming angles from polar to (-2,2] interval
true_angdens_L2$q = 2*true_angdens_L2$theta/pi
true_angdens_L2$fq = (pi/2)*true_angdens_L2$ft

#Setting plotting parameters
par(mfrow=c(1,2),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

#Plots comparing angular density estimates to truth for both sets of estimates
plot(true_angdens_L1$q,true_angdens_L1$fq,type="l",lwd=3,col=2,xlab="q",ylab=expression(f[Q](q)),main="Angular density",sub="L1 coordinates",ylim=range(L1_angular_density,true_angdens_L1$fq),cex.lab=1.2, cex.axis=1.2,cex.main=1.5)
lines(pred_Q,L1_angular_density,lwd=3,col=3)
legend(-2,1,legend=c("True","Estimated"),lwd=3,col=c(2,3),cex=1.2,bg="white")

#Plots comparing angular density estimates to truth for both sets of estimates
plot(true_angdens_L2$q,true_angdens_L2$fq,type="l",lwd=3,col=2,xlab="q",ylab=expression(f[Q](q)),main="Angular density",sub="L2 coordinates",ylim=range(L2_angular_density,true_angdens_L2$fq),cex.lab=1.2, cex.axis=1.2,cex.main=1.5)
lines(pred_Q,L2_angular_density,lwd=3,col=3)
legend(-2,1,legend=c("True","Estimated"),lwd=3,col=c(2,3),cex=1.2,bg="white")

#Setting plotting parameters
par(mfrow=c(1,2),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

#Plotting simulated data over original sample with L1 breakdown 
plot(example_data,xlab="X",ylab="Y",main="SPAR model simulations",sub="L1 coordinates",col="grey",pch=16,lwd=3,xlim=range(example_data,L1_simulation$data_sample),ylim=range(example_data,L1_simulation$data_sample),cex.lab=1.2, cex.axis=1.2,cex.main=1.5)
points(L1_simulation$data_sample,pch=16,col=adjustcolor(3,alpha.f = 0.2))
legend(range(example_data,L1_simulation$data_sample)[1],range(example_data,L1_simulation$data_sample)[2],legend=c("Observerd","Simulated"),pch=16,col=c("grey",adjustcolor(3,alpha.f = 0.2)),cex=1.2,bg="white")

#Plotting simulated data over original sample with L2 breakdown 
plot(example_data,xlab="X",ylab="Y",main="SPAR model simulations",sub="L2 coordinates",col="grey",pch=16,lwd=3,xlim=range(example_data,L2_simulation$data_sample),ylim=range(example_data,L2_simulation$data_sample),cex.lab=1.2, cex.axis=1.2,cex.main=1.5)
points(L2_simulation$data_sample,pch=16,col=adjustcolor(3,alpha.f = 0.2))
legend(range(example_data,L2_simulation$data_sample)[1],range(example_data,L2_simulation$data_sample)[2],legend=c("Observerd","Simulated"),pch=16,col=c("grey",adjustcolor(3,alpha.f = 0.2)),cex=1.2,bg="white")
