rm(list=ls())
#Checking for required packages. This function will install any required packages if they are not already installed
packages = c("evgam","copula","evd","parallel","circular")
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

#Standard Laplace quantile function. Transforms data from uniform margins to standard Laplace
Laplace_inverse = function(u){ 
  x = c()
  x[u<=0.5] = log(2*u[u<=0.5])
  x[u>0.5] = -log(2*(1-u[u>0.5]))
  return(x)
}

#Standard Laplace cumulative distribution function. Transforms data from standard Laplace to uniform margins
Laplace_cdf = function(x){ 
  u = c()
  u[x<0] = exp(x[x<0])/2
  u[x>=0] = 1-exp(-x[x>=0])/2
  return(u)
}

#Updated sign function, ensures sgn(0) = 1 
sgn = function(y){
  if(y>=0){
    return(1)
  } else {
    return(-1)
  }
}

#Radial function associated with the L1 norm 
L1_rad = function(vec){
  return(abs(vec[1])+abs(vec[2]))
}

#Angular function associated with the L1 norm
L1_ang = function(vec){
  return(sgn(vec[2])*(1-vec[1]))
}

#Radial function associated with the L2 norm 
L2_rad = function(vec){
  return(sqrt(vec[1]^2+vec[2]^2))
}

#Angular function associated with the L2 norm
L2_ang = function(vec){
  return(sgn(vec[2])*(2/pi)*acos(vec[1]))
}

#Wrapper for fitting the SPAR modelling framework 
SPAR_smooth = function(sample_data,norm_choice="L1",thresh_prob,k,k_shape=NULL,pred_Q=seq(-2,2,length.out=201)){ #wrapper for fitting the SPAR model 
  #sample_data - this denotes the bivariate data centered at (0,0) for which to fit SPAR. Must be an n x 2 matrix, where n denotes number of observations
  #norm_choice - this denotes the choice of norm. Must equal either "L1" or "L2". Defaults to "L1"
  #thresh_prob - non-exceedance probability for which to estimate threshold. Must be a probability in (0,1)
  #k - spline basis dimension for threshold and scale functions. Must be a natural number greater that 3
  #k_shape - spline basis dimension for shape functions. Must be equal to a natural number greater that 3, or NULL. Defaults to NULL, corresponding to a constant shape
  #pred_Q - angular values at which to evaluate estimated functions. This must be a grid of angles in [-2,2]. Defaults to seq(-2,2,length.out=201)
  
  if(!is.matrix(sample_data) | ncol(sample_data) != 2){
    stop("sample_data must be an n x 2 matrix, where n denotes number of observations")
  }
  if(!(norm_choice %in% c("L1","L2"))){
    stop("norm_choice must equal either 'L1' or 'L2'")
  }
  if(thresh_prob >= 1 | thresh_prob <= 0 | length(thresh_prob) != 1){
    stop("thresh_prob must be a single probability in (0,1)")
  }
  if(round(k) != k | k<=3){
    stop("k must be a natural number greater that 3")
  }
  if(!is.null(k_shape)){
    if(round(k_shape) != k_shape | k_shape<=3){
      stop("If not set to NULL (i.e., constant shape), k_shape must be a natural number greater that 3")
    }
  }
  if(min(pred_Q) < -2 | max(pred_Q)>2 | length(pred_Q)<=1){
    stop("pred_Q should be a grid of angles in the interval [-2,2]")
  }
  
  
  #We first transform data to radial-angular coordinates
  #All radial/angular/sgn functions are repeated inside the SPAR_smooth function. This allows for easier parallelisation of the code
  
  sgn = function(y){
    if(y>=0){
      return(1)
    } else {
      return(-1)
    }
  }
  
  if(norm_choice == "L1"){
    
    L1_rad = function(vec){
      return(abs(vec[1])+abs(vec[2]))
    }
    
    L1_ang = function(vec){
      return(sgn(vec[2])*(1-vec[1]))
    }
    
    #observed radial data
    R = apply(sample_data,1,L1_rad)
    
    #observed angular data
    Q = apply( apply(sample_data,2,function(x,r){return(x/r)},r=R),1,L1_ang)
    
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
    R = apply(sample_data,1,L2_rad)
    
    #observed angular data
    Q = apply( apply(sample_data,2,function(x,r){return(x/r)},r=R),1,L2_ang)
    
    #dataframe of angular-radial data
    polar_data = data.frame(R=R,Q=Q)
    
  }
  
  #transform radial component to the log scale. This guarantees positivity in the resulting threshold function
  polar_data$logR = log(polar_data$R)
  
  #Formulation for threshold function. s(Q, bs='cc',k=k) is a cyclic cubic spline over Q with basis dimension k
  fmla_ald = paste0("logR ~ s(Q, bs='cc',k=",k,")")
  
  #Fit asymmetric Laplace distribution to estimate threshold function using evgam
  m_ald = evgam(as.formula(fmla_ald), data=polar_data, family="ald", ald.args=list(tau=thresh_prob)) 
  
  #Estimated threshold functions. The exp function is used to transform back to R from the logR scale
  gpd_thresh_function = exp(predict(m_ald, newdata=polar_data)$location) 
  
  #all exceedances of threshold function
  gpd_thresh_exceedances = polar_data$R - gpd_thresh_function
  
  #checking which exceedances are positive 
  positive_exceedances_indicator = which(gpd_thresh_exceedances>0)
  
  #extracting positive exceedances
  gpd_thresh_exceedances = gpd_thresh_exceedances[positive_exceedances_indicator] 
  
  #dataframe of threshold exceedance data. We only save exceedances - we don't need the rest of the radial data
  polar_exceedance_data = data.frame(Q=polar_data$Q[positive_exceedances_indicator],R_exc=gpd_thresh_exceedances)
  
  #We now fit the non-stationary generalised Pareto distribution (GPD) using evgam
  
  #Checking if the shape is assumed to be constant
  if(is.null(k_shape)){
    
    #formulation for scale function
    spl=paste0("R_exc ~ s(Q, bs='cc',k=",k,")")
    
    #formulation for both scale and shape. ~1 specifies a constant shape
    fmla_gpd = list(as.formula(spl), ~1) 
    
    #fit non-stationary GPD model using evgam 
    m_gpd = evgam(fmla_gpd, data=polar_exceedance_data, family="gpd") 
    
    #obtain GPD parameter function estimates for values in pred_Q interval
    pred_para = predict(m_gpd, newdata=data.frame(Q=pred_Q),type="response")
    
  } else {
    
    #define second angular component that is identical as first. This is so that we can define separate splines for the shape and scale 
    polar_exceedance_data$Q2 = polar_exceedance_data$Q
    
    #formulation for both scale and shape. Note that k and k_shape are unlikely to be the same
    fmla_gpd = list(as.formula(paste0("R_exc ~ s(Q, bs='cc',k=",k,")")), as.formula(paste0("R_exc ~ s(Q2, bs='cc',k=",k_shape,")"))) 
    
    #fit non-stationary GPD model using evgam 
    m_gpd = evgam(fmla_gpd, data=polar_exceedance_data, family="gpd") 
    
    #obtain GPD parameter function estimates for values in pred_Q interval
    pred_para = predict(m_gpd,type = "response",newdata=data.frame(Q=pred_Q,Q2=pred_Q))
    
  }
  
  #obtain threshold function estimate for values in pred_Q interval
  pred_thresh = exp(predict(m_ald, newdata=data.frame(Q=pred_Q))$location) 
  
  exp_quants = qexp(1 - ( 1 + ((predict(m_gpd,type = "response",newdata=polar_exceedance_data))$shape/(predict(m_gpd,type = "response",newdata=polar_exceedance_data))$scale)*(polar_exceedance_data$R_exc)  )^(-1/(predict(m_gpd,type = "response",newdata=polar_exceedance_data))$shape))
  
  #return estimated threshold and parameter functions, along with the angular grid and non-exceedance probability. These are given in a list
  return(list(pred_thresh=pred_thresh,pred_para=pred_para,pred_Q=pred_Q,thresh_prob=thresh_prob,exp_quants=exp_quants))
  
}

#Wrapper for fitting the SPAR modelling framework 
SPAR_angular_density = function(sample_data,norm_choice="L1",pred_Q=seq(-2,2,length.out=201),bw=50){ #wrapper for estimating the angular density in the SPAR framework
  #sample_data - this denotes the bivariate data centered at (0,0) for which to fit SPAR. Must be an n x 2 matrix, where n denotes number of observations
  #norm_choice - this denotes the choice of norm. Must equal either "L1" or "L2". Defaults to "L1"
  #pred_Q - angular values at which to evaluate estimated functions. This must be a grid of angles in [-2,2]. Defaults to seq(-2,2,length.out=201)
  #bw - this is the bandwidth parameter for the kernel density estimation of the circular angular density. Should be a positive number. Defaults to 50
  
  if(!is.matrix(sample_data) | ncol(sample_data) != 2){
    stop("sample_data must be an n x 2 matrix, where n denotes number of observations")
  }
  if(!(norm_choice %in% c("L1","L2"))){
    stop("norm_choice must equal either 'L1' or 'L2'")
  }
  if(min(pred_Q)< -2 | max(pred_Q)>2 | length(pred_Q)<=1){
    stop("pred_Q should be a grid of angles in the interval [-2,2]")
  }
  if(bw < 0 | length(bw) != 1){
    stop("bw should be a positive real number")
  }
  
  #We first transform data to radial-angular coordinates
  #All radial/angular/sgn functions are repeated inside this function. This allows for easier parallelisation of the code
  
  sgn = function(y){
    if(y>=0){
      return(1)
    } else {
      return(-1)
    }
  }
  
  if(norm_choice == "L1"){
    
    L1_rad = function(vec){
      return(abs(vec[1])+abs(vec[2]))
    }
    
    L1_ang = function(vec){
      return(sgn(vec[2])*(1-vec[1]))
    }
    
    #observed radial data
    R = apply(sample_data,1,L1_rad)
    
    #observed angular data
    Q = apply( apply(sample_data,2,function(x,r){return(x/r)},r=R),1,L1_ang)
    
  } else {
    
    L2_rad = function(vec){
      return(sqrt(vec[1]^2+vec[2]^2))
    }
    
    L2_ang = function(vec){
      return(sgn(vec[2])*(2/pi)*acos(vec[1]))
    }
    
    #observed radial data
    R = apply(sample_data,1,L2_rad)
    
    #observed angular data
    Q = apply( apply(sample_data,2,function(x,r){return(x/r)},r=R),1,L2_ang)
    
  }
  
  #we transform the angular component to be within the interval [0,2pi). This is because the density.circular function is defined only for standard polar angles   
  scale_Q = (Q+2)*pi/2
  
  #estimating the density function using the `circular' R package with a von Mises kernel
  dens_est = density.circular(as.circular(scale_Q,type="angles",units="radians",template="none",modulo="2pi",zero=0,rotation="counter"), bw=bw,kernel="vonmises")
  
  #Making sure we equal 2pi at the endpoint 
  dens_est$x[length(dens_est$x)] = 2*pi
  
  #obtaining the estimated density function on the original angular scale. (pi/2) represents the corresponding Jacobian 
  f_Q = approxfun(x = (2*dens_est$x)/pi - 2,y = (pi/2)*dens_est$y) 
  
  #obtain density estimates for values in pred_Q interval
  return(f_Q(pred_Q))
}

#function for computing equidensity contours using the SPAR approximation 
#density of R|Q=q is given by a GPD, while the density of Q is estimated non-parametrically 
#this formulation is treating the angle Q as fixed, i.e., Q=q
joint_dens_root_func = function(r,u,scale,shape,thresh_prob,f_q,dens_val){
  #r - this is the radial value on the equidensity contour that we wish to compute for Q=q
  #u - estimated threshold at Q=q
  #scale - estimated scale parameter at Q=q
  #shape - estimated shape parameter at Q=q
  #thresh_prob - non-exceedance probability for which to estimate threshold. Must be a probability in (0,1)
  #f_q - estimated density value at Q=q
  #dens_val - the density value along the equidensity contours. Note this is constrained by the estimated threshold function - see paper for further details
  
  if(u <= 0 | length(u) != 1){
    stop("u must be a positive real number")
  }
  if(scale <= 0 | length(scale) != 1){
    stop("scale must be a positive real number")
  }
  if(thresh_prob >= 1 | thresh_prob <= 0 | length(thresh_prob) != 1){
    stop("thresh_prob must be a single probability in (0,1)")
  }
  if(f_q <= 0 | length(f_q) != 1){
    stop("f_q must be a positive real number")
  }
  
  #computing gpd density for R|Q=q
  gpd_density = ((scale^(-1)) * pmax((1 + shape * (r - u)/scale),0)^((-1/shape) - 1))
  
  #estimated joint density at (r,q)
  estimated_joint_dens = ((1-thresh_prob)*f_q/r)*gpd_density
  
  #difference between desired density value and estimated density. We use this difference, along with the root function, to find the corresponding radial value 
  return(dens_val - estimated_joint_dens)
}

#Function for estimating equidensity contours from fitted SPAR model
SPAR_equidensity_contours = function(density_levels,norm_choice="L1",SPAR_GPD,SPAR_ang){
  #density_level - level for which equidensity contour is required. Note that this value is bounded above by the fact we do not model data below the threshold function
  #norm_choice - this denotes the choice of norm. Must equal either "L1" or "L2". Defaults to "L1"
  #SPAR_GPD - output from either SPAR_smooth or SPAR_local. Must be a list
  #SPAR_ang - output from SPAR_angular_density. Must be a vector
  
  dens_lb = max( (1-SPAR_GPD$thresh_prob)*SPAR_ang/SPAR_GPD$pred_thresh )
  
  if(sum(density_levels < 0) > 0 | !is.vector(density_levels)){
    stop("density_levels should be a vector of positive real numbers (or a single positive real number)")
  }
  if(min(density_levels)>dens_lb){
    stop("At least one of the density_levels is too high, such that the SPAR model is not valid. Please input values greater than")
  }
  if(!(norm_choice %in% c("L1","L2"))){
    stop("norm_choice must equal either 'L1' or 'L2'")
  }
  if(!is.list(SPAR_GPD)){
    stop("SPAR_GPD must be a list obtained as output from either SPAR_smooth or SPAR_local")
  }
  if(!is.vector(SPAR_ang)){
    stop("SPAR_ang must be a vector obtained as output from SPAR_angular_density")
  }
  
  #We check which coordinate system we are working in
  if(norm_choice == "L1"){
    
    #Defining points on the unit sphere for the L1 norm
    u_vec = ifelse(SPAR_GPD$pred_Q>=0,(1-SPAR_GPD$pred_Q),(SPAR_GPD$pred_Q+1))
    v_vec = ifelse(SPAR_GPD$pred_Q>=0, 1-abs(u_vec),-1+abs(u_vec))
    
  } else {
    
    #Defining points on the unit sphere for the L2 norm
    u_vec = ifelse(SPAR_GPD$pred_Q>=0,cos(pi*SPAR_GPD$pred_Q/2),cos(-pi*SPAR_GPD$pred_Q/2))
    v_vec = ifelse(SPAR_GPD$pred_Q>=0, sqrt(1-u_vec^2),-sqrt(1-u_vec^2))
    
  }
  
  #Computing radii values at density value for each fixed angle
  dens_radii = matrix(NA,ncol = length(density_levels),nrow=length(pred_Q))
  for(i in 1:length(density_levels)){
    for(j in 1:length(pred_Q)){
      dens_radii[j,i] = uniroot(joint_dens_root_func,interval = c(SPAR_GPD$pred_thresh[j],50),u=SPAR_GPD$pred_thresh[j],scale=SPAR_GPD$pred_para$scale[j],shape=SPAR_GPD$pred_para$shape[j],thresh_prob=SPAR_GPD$thresh_prob,f_q=SPAR_ang[j],dens_val=density_levels[i])$root
    }
  }
  
  #Computing the corresponding contours in cartesian coordinates
  estimated_density_contours = list()
  for(i in 1:length(density_levels)){
    estimated_density_contours[[i]] = cbind(dens_radii[,i]*u_vec,dens_radii[,i]*v_vec); estimated_density_contours[[i]] = rbind(estimated_density_contours[[i]],estimated_density_contours[[i]][1,])
  }
  
  if(length(estimated_density_contours)==1){
    estimated_density_contours = estimated_density_contours[[1]]
  }
  
  #returning estimated density contour(s). This is in the form of a list of matrices, or a singular matrix if just one density level is entered.
  return(estimated_density_contours)
}

#Function for estimating return level sets from fitted SPAR model
SPAR_ret_level_sets = function(ret_period,obs_year,norm_choice="L1",SPAR_GPD){
  #ret_period - the return period for evaluating the return level set
  #obs_year - the number of observations per year within the data
  #norm_choice - this denotes the choice of norm. Must equal either "L1" or "L2". Defaults to "L1"
  #SPAR_GPD - output from either SPAR_smooth or SPAR_local. Must be a list
  
  if( ret_period < 0 | length(ret_period)!=1){
    stop("ret_period must be a positive real number")
  }
  if( obs_year < 0 | length(obs_year)!=1){
    stop("obs_year must be a positive real number")
  }
  if(!(norm_choice %in% c("L1","L2"))){
    stop("norm_choice must equal either 'L1' or 'L2'")
  }
  if(!is.list(SPAR_GPD)){
    stop("SPAR_GPD must be a list obtained as output from either SPAR_smooth or SPAR_local")
  }
  
  #We check which coordinate system we are working in
  if(norm_choice == "L1"){
    
    #Defining points on the unit sphere for the L1 norm
    u_vec = ifelse(SPAR_GPD$pred_Q>=0,(1-SPAR_GPD$pred_Q),(SPAR_GPD$pred_Q+1))
    v_vec = ifelse(SPAR_GPD$pred_Q>=0, 1-abs(u_vec),-1+abs(u_vec))
    
  } else {
    
    #Defining points on the unit sphere for the L2 norm
    u_vec = ifelse(SPAR_GPD$pred_Q>=0,cos(pi*SPAR_GPD$pred_Q/2),cos(-pi*SPAR_GPD$pred_Q/2))
    v_vec = ifelse(SPAR_GPD$pred_Q>=0, sqrt(1-u_vec^2),-sqrt(1-u_vec^2))
    
  }
  
  #Corresponding probability for return level set 
  prob = 1/(ret_period*obs_year)
  
  #Checking if the probability is valid
  if(prob > (1 - SPAR_GPD$thresh_prob)){
    stop("Return period too small such that the SPAR model is no longer valid. Please try a larger return period value")
  }
  
  #Computing radii values of return level set for fixed angles
  ret_level_set = SPAR_GPD$pred_thresh + (SPAR_GPD$pred_para$scale/SPAR_GPD$pred_para$shape)*( (prob/(1 - SPAR_GPD$thresh_prob))^(-SPAR_GPD$pred_para$shape) - 1 )
  
  #Computing return level set contour on original scale
  ret_level_contour = cbind(ret_level_set*u_vec,ret_level_set*v_vec); ret_level_contour = rbind(ret_level_contour,ret_level_contour[1,])
  
  #returning estimated return level set
  return(ret_level_contour)
}

#Wrapper for fitting the empirical version of the SPAR model 
SPAR_local = function(sample_data,norm_choice="L1",thresh_prob,pred_Q=seq(-2,2,length.out=201),num_neigh){ #wrapper for fitting the SPAR model 
  #sample_data - this denotes the bivariate data centered at (0,0) for which to fit SPAR. Must be an n x 2 matrix, where n denotes number of observations
  #norm_choice - this denotes the choice of norm. Must equal either "L1" or "L2". Defaults to "L1"
  #thresh_prob - non-exceedance probability for which to estimate threshold. Must be a probability in (0,1)
  #pred_Q - angular values at which to evaluate estimated functions. This must be a grid of angles in [-2,2]. Defaults to seq(-2,2,length.out=201)
  #num_neigh - the number of neighbours in each local windows. This must be a positive integer greater than 50 for numerical stability. Defaults to 500 
  
  if(!is.matrix(sample_data) | ncol(sample_data) != 2){
    stop("sample_data must be an n x 2 matrix, where n denotes number of observations")
  }
  if(!(norm_choice %in% c("L1","L2"))){
    stop("norm_choice must equal either 'L1' or 'L2'")
  }
  if(thresh_prob >= 1 | thresh_prob <= 0 | length(thresh_prob) != 1){
    stop("thresh_prob must be a single probability in (0,1)")
  }
  if(min(pred_Q) < -2 | max(pred_Q)>2 | length(pred_Q)<=1){
    stop("pred_Q should be a grid of angles in the interval [-2,2]")
  }
  if(num_neigh <= 50 | length(num_neigh) != 1){
    stop("num_neigh must be a positive real number greater than 50")
  }
  
  #We first transform data to radial-angular coordinates
  #All radial/angular/sgn functions are repeated inside the SPAR_smooth function. This allows for easier parallelisation of the code
  
  sgn = function(y){
    if(y>=0){
      return(1)
    } else {
      return(-1)
    }
  }
  
  if(norm_choice == "L1"){
    
    L1_rad = function(vec){
      return(abs(vec[1])+abs(vec[2]))
    }
    
    L1_ang = function(vec){
      return(sgn(vec[2])*(1-vec[1]))
    }
    
    #observed radial data
    R = apply(sample_data,1,L1_rad)
    
    #observed angular data
    Q = apply( apply(sample_data,2,function(x,r){return(x/r)},r=R),1,L1_ang)
    
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
    R = apply(sample_data,1,L2_rad)
    
    #observed angular data
    Q = apply( apply(sample_data,2,function(x,r){return(x/r)},r=R),1,L2_ang)
    
    #dataframe of angular-radial data
    polar_data = data.frame(R=R,Q=Q)
    
  }
  
  #wrapper for obtaining empirical radial windows for each of the reference angles
  SPAR_empirical_windows = function(ref_Q,polar_data,num_neigh){
    #ref_q - reference angle 
    #polar_data - angular radial dataframe 
    #num_neigh - number of neighbours for each local angular windows
    
    #finding differences between all observed angles and the reference angle. This metric is adjusted to account for periodicity in angles
    diff = apply(cbind(4 - abs(polar_data$Q - ref_Q),abs(polar_data$Q - ref_Q)),1,min)
    
    #computing indices of num_neigh closest observed angles to ref_Q
    eps.Q = sort(diff)[num_neigh]
    nei.Q = which(diff <= eps.Q)
    
    #returning radial component at angles within local window
    return(polar_data$R[nei.Q])
  }
  
  #function for obtaining local GPD parameter estimates via the ismev package
  #this function should NOT be used outside of the SPAR_local function
  SPAR_local_GPD_fit = function(radial_data,thresh_prob){
    #radial_data - vector containing radial values within the local window
    #thresh_prob - non-exceedance probability for which to estimate threshold.
    
    #empirical threshold estimate
    thresh = quantile(radial_data,thresh_prob)
    
    #estimated GPD parameters for the local window, obtained using the ismev package
    gpd_par = ismev::gpd.fit(radial_data, threshold=thresh, show = F)$mle
    
    #returns vector containaing estimated threshold and GPD parameters for the local window
    return(c(thresh,gpd_par))
  }
  
  #Computing local angular windows for radial component
  windows_datasets = sapply(pred_Q,SPAR_empirical_windows,polar_data=polar_data,num_neigh=num_neigh,simplify = F)
  
  #Obtaining local GPD estimates in each local window
  empirical_estimates = lapply(windows_datasets,SPAR_local_GPD_fit,thresh_prob=thresh_prob);empirical_estimates = do.call(rbind,empirical_estimates);
  
  #storing output from local estimation
  emp_thresh = empirical_estimates[,1]
  emp_para = data.frame(scale = empirical_estimates[,2], shape = empirical_estimates[,3])
  
  #return estimated threshold and parameter functions, along with the angular grid and non-exceedance probability. These are given in a list
  return(list(pred_thresh=emp_thresh,pred_para=emp_para,pred_Q=pred_Q,thresh_prob=thresh_prob))
}

#wrapper for obtaining empirical radial windows for each of the reference angles
SPAR_empirical_windows = function(ref_Q,polar_data,num_neigh){
  #ref_q - reference angle 
  #polar_data - angular radial dataframe 
  #num_neigh - number of neighbours for each local angular windows
  
  #finding differences between all observed angles and the reference angle. This metric is adjusted to account for periodicity in angles
  diff = apply(cbind(4 - abs(polar_data$Q - ref_Q),abs(polar_data$Q - ref_Q)),1,min)
  
  #computing indices of num_neigh closest observed angles to ref_Q
  eps.Q = sort(diff)[num_neigh]
  nei.Q = which(diff <= eps.Q)
  
  #returning radial component at angles within local window
  return(polar_data$R[nei.Q])
}

#Function for simulating from fitted SPAR model. This is only possible with the smooth model fits
SPAR_simulation = function(sample_data,nsim,norm_choice="L1",thresh_prob,k,k_shape=NULL,pred_Q=seq(-2,2,length.out=201),bw=50){ #wrapper for simulating from smooth model fits
  #sample_data - this denotes the bivariate data centered at (0,0) for which to fit SPAR. Must be an n x 2 matrix, where n denotes number of observations
  #nsim - number of simulations to obtain from fitted SPAR model
  #norm_choice - this denotes the choice of norm. Must equal either "L1" or "L2". Defaults to "L1"
  #thresh_prob - non-exceedance probability for which to estimate threshold. Must be a probability in (0,1)
  #k - spline basis dimension for threshold and scale functions. Must be a natural number greater that 3
  #k_shape - spline basis dimension for shape functions. Must be equal to a natural number greater that 3, or NULL. Defaults to NULL, corresponding to a constant shape
  #pred_Q - angular values at which to evaluate estimated functions. This must be a grid of angles in [-2,2]. Defaults to seq(-2,2,length.out=201)
  #bw - this is the bandwidth parameter for the kernel density estimation of the circular angular density. Should be a positive number. Defaults to 50
  
  if(!is.matrix(sample_data) | ncol(sample_data) != 2){
    stop("sample_data must be an n x 2 matrix, where n denotes number of observations")
  }
  if(nsim < 1 | length(nsim) != 1){
    stop("nsim must be a positive real integer")
  }
  if(!(norm_choice %in% c("L1","L2"))){
    stop("norm_choice must equal either 'L1' or 'L2'")
  }
  if(thresh_prob >= 1 | thresh_prob <= 0 | length(thresh_prob) != 1){
    stop("thresh_prob must be a single probability in (0,1)")
  }
  if(round(k) != k | k<=3){
    stop("k must be a natural number greater that 3")
  }
  if(!is.null(k_shape)){
    if(round(k_shape) != k_shape | k_shape<=3){
      stop("If not set to NULL (i.e., constant shape), k_shape must be a natural number greater that 3")
    }
  }
  if(min(pred_Q) < -2 | max(pred_Q)>2 | length(pred_Q)<=1){
    stop("pred_Q should be a grid of angles in the interval [-2,2]")
  }
  if(bw < 0 | length(bw) != 1){
    stop("bw must be a positive real number")
  }
  
  #We first transform data to radial-angular coordinates
  #All radial/angular/sgn functions are repeated inside this function. This allows for easier parallelisation of the code
  
  sgn = function(y){
    if(y>=0){
      return(1)
    } else {
      return(-1)
    }
  }
  
  if(norm_choice == "L1"){
    
    L1_rad = function(vec){
      return(abs(vec[1])+abs(vec[2]))
    }
    
    L1_ang = function(vec){
      return(sgn(vec[2])*(1-vec[1]))
    }
    
    #observed radial data
    R = apply(sample_data,1,L1_rad)
    
    #observed angular data
    Q = apply( apply(sample_data,2,function(x,r){return(x/r)},r=R),1,L1_ang)
    
  } else {
    
    L2_rad = function(vec){
      return(sqrt(vec[1]^2+vec[2]^2))
    }
    
    L2_ang = function(vec){
      return(sgn(vec[2])*(2/pi)*acos(vec[1]))
    }
    
    #observed radial data
    R = apply(sample_data,1,L2_rad)
    
    #observed angular data
    Q = apply( apply(sample_data,2,function(x,r){return(x/r)},r=R),1,L2_ang)
    
  }
  
  #dataframe of angular-radial data
  polar_data = data.frame(R=R,Q=Q)
  
  #transform radial component to the log scale. This guarantees positivity in the resulting threshold function
  polar_data$logR = log(polar_data$R)
  
  #Formulation for threshold function. s(Q, bs='cc',k=k) is a cyclic cubic spline over Q with basis dimension k
  fmla_ald = paste0("logR ~ s(Q, bs='cc',k=",k,")")
  
  #Fit asymmetric Laplace distribution to estimate threshold function using evgam
  m_ald = evgam(as.formula(fmla_ald), data=polar_data, family="ald", ald.args=list(tau=thresh_prob)) 
  
  #Estimated threshold functions. The exp function is used to transform back to R from the logR scale
  gpd_thresh_function = exp(predict(m_ald, newdata=polar_data)$location) 
  
  #all exceedances of threshold function
  gpd_thresh_exceedances = polar_data$R - gpd_thresh_function
  
  #checking which exceedances are positive 
  positive_exceedances_indicator = which(gpd_thresh_exceedances>0)
  
  #extracting positive exceedances
  gpd_thresh_exceedances = gpd_thresh_exceedances[positive_exceedances_indicator] 
  
  #dataframe of threshold exceedance data. We only save exceedances - we don't need the rest of the radial data
  polar_exceedance_data = data.frame(Q=polar_data$Q[positive_exceedances_indicator],R_exc=gpd_thresh_exceedances)
  
  #We now fit the non-stationary generalised Pareto distribution (GPD) using evgam
  
  #Checking if the shape is assumed to be constant
  if(is.null(k_shape)){
    
    #formulation for scale function
    spl=paste0("R_exc ~ s(Q, bs='cc',k=",k,")")
    
    #formulation for both scale and shape. ~1 specifies a constant shape
    fmla_gpd = list(as.formula(spl), ~1) 
    
    #fit non-stationary GPD model using evgam 
    m_gpd = evgam(fmla_gpd, data=polar_exceedance_data, family="gpd") 
    
  } else {
    
    #define second angular component that is identical as first. This is so that we can define separate splines for the shape and scale 
    polar_exceedance_data$Q2 = polar_exceedance_data$Q
    
    #formulation for both scale and shape. Note that k and k_shape are unlikely to be the same
    fmla_gpd = list(as.formula(paste0("R_exc ~ s(Q, bs='cc',k=",k,")")), as.formula(paste0("R_exc ~ s(Q2, bs='cc',k=",k_shape,")"))) 
    
    #fit non-stationary GPD model using evgam 
    m_gpd = evgam(fmla_gpd, data=polar_exceedance_data, family="gpd") 
    
  }
  
  #we transform the angular component to be within the interval [0,2pi). This is because the density.circular function is defined only for standard polar angles   
  scale_Q = (Q+2)*pi/2
  
  #estimating the density function using the `circular' R package with a von Mises kernel
  dens_est = density.circular(as.circular(scale_Q,type="angles",units="radians",template="none",modulo="2pi",zero=0,rotation="counter"), bw=bw,kernel="vonmises")
  
  #Making sure we equal 2pi at the endpoint 
  dens_est$x[length(dens_est$x)] = 2*pi
  
  #obtaining the estimated density function on the original angular scale. (pi/2) represents the corresponding Jacobian 
  f_Q = approxfun(x = (2*dens_est$x)/pi - 2,y = (pi/2)*dens_est$y)
  
  kd_integral = function(u,x,f_Q){
    if(x>0.99999){
      return(1 - u) #numerical error for really large probabilities. This ensures we integrate to 1
    } else {
      return(integrate(f_Q,-2,x)$value - u)
    }
  }
  
  kd_root = function(u,f_Q){
    return(uniroot(f=kd_integral,interval = c(-2,2),u=u,f_Q = f_Q)$root)
  }
  
  unif_sample = runif(nsim)
  
  Q_sample = sapply(unif_sample,kd_root,f_Q=f_Q)
  
  #obtain threshold function estimates for sample angles  
  thresh_sample = exp(predict(m_ald, newdata=data.frame(Q=Q_sample))$location) 
  
  if(is.null(k_shape)){
   
    #obtain GPD parameter function estimates for sample angles
    para_sample = predict(m_gpd, newdata=data.frame(Q=Q_sample),type="response")
    
  } else {
    
    #obtain GPD parameter function estimates for sample angles
    para_sample = predict(m_gpd,type = "response",newdata=data.frame(Q=Q_sample,Q2=Q_sample))
    
  }
  
  unif_sample = runif(nsim)
  
  R_sample = thresh_sample + (para_sample$scale/para_sample$shape)*(unif_sample^(-para_sample$shape)-1)
  
  #We check which coordinate system we are working in
  if(norm_choice == "L1"){
    
    #Defining points on the unit sphere for the L1 norm
    u_vec = ifelse(Q_sample>=0,(1-Q_sample),(Q_sample+1))
    v_vec = ifelse(Q_sample>=0, 1-abs(u_vec),-1+abs(u_vec))
    
  } else {
    
    #Defining points on the unit sphere for the L2 norm
    u_vec = ifelse(Q_sample>=0,cos(pi*Q_sample/2),cos(-pi*Q_sample/2))
    v_vec = ifelse(Q_sample>=0, sqrt(1-u_vec^2),-sqrt(1-u_vec^2))
    
  }
  
  data_sample = cbind(R_sample*u_vec,R_sample*v_vec)
  
  #returning simulated datasets
  return(list(Q_sample = Q_sample,R_sample = R_sample,data_sample = data_sample))
}
