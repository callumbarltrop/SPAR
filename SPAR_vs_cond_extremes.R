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

# Fitting SPAR geometric method ------------------------------------------------------------

#The below code fits the SPAR model using both angular systems. The threshold and parameter functions are estimated smoothly with GAMs

smooth_fit = SPAR_smooth_polar(sample_data = example_data,thresh_prob = thresh_prob,k=k,pred_phi = pred_phi)

angular_density = SPAR_angular_density_polar(sample_data = example_data,pred_phi = pred_phi,bw = bw)

#Number of observations to simulate from the SPAR model
nsim = n

#obtaining model simulations
simulated_data = SPAR_simulation_polar(nsim=nsim,SPAR_GPD=smooth_fit,SPAR_ang=angular_density)


# Heffernan Tawn conditional extremes model ----------------------------------------------------------

source("cond_extremes_functions.R")

heff_tawn_model_fit = heff_tawn_function(data_laplace = example_data,cond_var = 1,q_HT = 0.9)

nsim = n

heff_tawn_sim_data = heff_tawn_simulation(fit=heff_tawn_model_fit,nsim = n)

# Comparing ---------------------------------------------------------------

pdf(file="plots/SPAR_vs_cond_extremes_prefit.pdf",width=12,height=6)

#Setting plotting parameters
par(mfrow=c(1,2),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

plot(example_data,pch=16,cex=0.5,col="grey",xlab="X",ylab="Y",main="Conditional extremes simulations",lwd=3,xlim=range(example_data,simulated_data$data_sample),ylim=range(example_data,simulated_data$data_sample),cex.lab=1.2, cex.axis=1.2,cex.main=1.55,cex.sub=1.3)

#Plotting simulated data over original sample with L1 breakdown 
plot(example_data,col="grey",pch=16,cex=.5,xlab="X",ylab="Y",main="SPAR model simulations",lwd=3,xlim=range(example_data,simulated_data$data_sample),ylim=range(example_data,simulated_data$data_sample),cex.lab=1.2, cex.axis=1.2,cex.main=1.55,cex.sub=1.3)

dev.off()


pdf(file="plots/SPAR_vs_cond_extremes_simulation.pdf",width=12,height=6)

#Setting plotting parameters
par(mfrow=c(1,2),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

plot(heff_tawn_sim_data,xlab="X",ylab="Y",main="Conditional extremes simulations",col=adjustcolor(2,alpha.f = 0.2),pch=16,lwd=3,xlim=range(example_data,simulated_data$data_sample),ylim=range(example_data,simulated_data$data_sample),cex.lab=1.2, cex.axis=1.2,cex.main=1.55,cex.sub=1.3)
points(example_data[example_data[,heff_tawn_model_fit$cond_var] > heff_tawn_model_fit$u_cond_var,],pch=16,cex=0.5,col="grey")
abline(v=heff_tawn_model_fit$u_cond_var,lwd=4,col=4)

#Plotting simulated data over original sample with L1 breakdown 
plot(simulated_data$data_sample,pch=16,col=adjustcolor(2,alpha.f = 0.2),xlab="X",ylab="Y",main="SPAR model simulations",lwd=3,xlim=range(example_data,simulated_data$data_sample),ylim=range(example_data,simulated_data$data_sample),cex.lab=1.2, cex.axis=1.2,cex.main=1.55,cex.sub=1.3)
points(example_data[smooth_fit$polar_data$R > smooth_fit$polar_data$thresh_func,],col="grey",pch=16,cex=.5)
lines(smooth_fit$pred_thresh*cos(smooth_fit$pred_phi),smooth_fit$pred_thresh*sin(smooth_fit$pred_phi),lwd=4,col=4)
legend(range(example_data,simulated_data$data_sample)[1],range(example_data,simulated_data$data_sample)[2],legend=c("Observerd","Simulated"),pch=16,col=c("grey",adjustcolor(2,alpha.f = 0.2)),cex=1.2,bg="white")

dev.off()


pdf(file="plots/SPAR_vs_cond_extremes_comparison.pdf",width=12,height=6)

#Setting plotting parameters
par(mfrow=c(1,2),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

plot(heff_tawn_sim_data,xlab="X",ylab="Y",main="Conditional extremes simulations",col=adjustcolor(2,alpha.f = 0.2),pch=16,lwd=3,xlim=range(example_data,simulated_data$data_sample),ylim=range(example_data,simulated_data$data_sample),cex.lab=1.2, cex.axis=1.2,cex.main=1.55,cex.sub=1.3)
points(example_data[example_data[,heff_tawn_model_fit$cond_var] > heff_tawn_model_fit$u_cond_var,],pch=16,cex=0.5,col="grey")
abline(v=heff_tawn_model_fit$u_cond_var,lwd=4,col=4)

#Plotting simulated data over original sample with L1 breakdown 
plot(simulated_data$data_sample,pch=16,col=adjustcolor(2,alpha.f = 0.2),xlab="X",ylab="Y",main="SPAR model simulations",lwd=3,xlim=range(example_data,simulated_data$data_sample),ylim=range(example_data,simulated_data$data_sample),cex.lab=1.2, cex.axis=1.2,cex.main=1.55,cex.sub=1.3)
points(example_data[smooth_fit$polar_data$R > smooth_fit$polar_data$thresh_func,],col="grey",pch=16,cex=.5)
lines(smooth_fit$pred_thresh*cos(smooth_fit$pred_phi),smooth_fit$pred_thresh*sin(smooth_fit$pred_phi),lwd=4,col=4)
legend(range(example_data,simulated_data$data_sample)[1],range(example_data,simulated_data$data_sample)[2],legend=c("Observerd","Simulated"),pch=16,col=c("grey",adjustcolor(2,alpha.f = 0.2)),cex=1.2,bg="white")

u2 = 4

points(simulated_data$data_sample[ simulated_data$data_sample[,1] > u2 ,],pch=16,col=adjustcolor(3,alpha.f = 0.2),lwd=3)
points((example_data[smooth_fit$polar_data$R > smooth_fit$polar_data$thresh_func,])[(example_data[smooth_fit$polar_data$R > smooth_fit$polar_data$thresh_func,])[,1]>u2,],col="grey2",pch=16,cex=.5)
abline(v=u2,lwd=4,col="blue3")

dev.off()
