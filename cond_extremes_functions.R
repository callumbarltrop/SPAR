#To impose bounds on the alpha and beta parameters, we use link functions 
#This speeds up optimisation of Heffernan Tawn model fitting. 

inv_link_func = function(y,a,b){
  return( log( (y-a)/(b-y) ) )
}

link_func = function(x,a,b){ 
  return( (a + b*exp(x))/(1+exp(x)) )
}

inv_link_func_2 = function(y){ 
  return( log( 1 - y ) )
}

link_func_2 = function(x){
  return( 1 - exp(x) )
}

heff_tawn_function = function(data_laplace,cond_var,q_HT){ #fits Heffernan Tawn model to data_laplace, conditioning on cond_var above quantile q_HT
  #data_matrix is a matrix of all data, with columns corresponding to different locations/variables
  #con_var is the conditioning variable; must be between 1 and ncol(data_matrix)
  if(!(cond_var %in%1:ncol(data_laplace))){
    stop("cond_var must be between 1 and ncol(data_laplace)")
  }
  #q_HT is the quantile of the conditioning variable we use to fit the Heffernan Tawn model
  #q_pred is the threshold we predict above in our model
  
  u = Laplace_inverse(q_HT) #obtain quantile value from Laplace distribution
  
  data_laplace = data_laplace[data_laplace[,cond_var] > u, ] #Extract rows where conditioning variable exceeds some threshold
  
  non_cond_var = (1:ncol(data_laplace))[-cond_var] #non conditioning sites
  para_estimates = apply(as.matrix(data_laplace[,non_cond_var]),2,HT_fit_wrapper,cond_data = data_laplace[,cond_var]) #for each pair of sites, we get parameter estimates
  
  return(list(cond_var=cond_var, cond_data_laplace = data_laplace, u_cond_var = u, HT_parameter_estimates = para_estimates))
}

HT_fit_wrapper = function(other_data,cond_data){ #function for estimating parameters across pairs of sites. Speeds up fitting to have this in a wrapper
  optim_current = optim(fn=HeffTawnNegLL,cond_site=cond_data,other_site=other_data,par=c(inv_link_func(1/2,-1,1),inv_link_func_2(1/2) ,log(1/2),1/2),method="BFGS",control = list(maxit=100000))
  optim_new = optim(fn=HeffTawnNegLL,cond_site=cond_data,other_site=other_data,par=optim_current$par,method="Nelder-Mead",control = list(maxit=100000))
  while(abs(optim_new$val - optim_current$val)>=10^(-10)){ #this is to ensure we converge to the global minimum (optimum parameter values). Keep putting parameter estimates back into optimiser till it is stable 
    optim_current = optim_new
    optim_new = optim(fn=HeffTawnNegLL,cond_site=cond_data,other_site=other_data,optim_current$par,method = "Nelder-Mead",control = list(maxit=100000))
  }
  optim_outputs = optim_current
  return(c(link_func(optim_outputs$par[1],-1,1),link_func_2(optim_outputs$par[2]))) #link functions make everything quicker
}

HeffTawnNegLL=function(cond_site,other_site,par){ #negative log likelihood for estimation of HT model parameters. We minimise this function. See https://rss.onlinelibrary.wiley.com/doi/full/10.1111/j.1467-9868.2004.02050.x
  #cond_site is conditioning variable, other_site is other variable
  alpha=link_func(par[1],-1,1) #links save us having to use `if` statements for checking validity 
  beta=link_func_2(par[2])
  sig=exp(par[3])
  mu=par[4]
  negloglik = -sum(dnorm(other_site,alpha*cond_site+mu*((cond_site)^beta),sig*((cond_site)^beta),log=T))
  if(is.finite(negloglik)){
    return(negloglik)
  } else {
    return(1e10) #if we diverge, we set the negative log likelihood to a big value and keep going.
  }
}

heff_tawn_simulation=function(fit,nsim){
  
  non_cond_var = (1:ncol(fit$cond_data_laplace))[-fit$cond_var] #non conditioning sites
  non_cond_data = as.matrix(fit$cond_data_laplace[,non_cond_var]) #non conditioning variables 
  
  cond_data = as.matrix(fit$cond_data_laplace[,fit$cond_var]) #conditioning variable, above some high threshold 
  
  z = apply(rbind(fit$HT_parameter_estimates,matrix(rep(cond_data,each=2),ncol=2,byrow=T),non_cond_data),2,function(x){ #this looks complicated, but it is just a clever way of writing the function to work in the general multivariate case
    #we are computing z = (y - alpha*x)/(x^beta), where x is conditioning variable and y is the observed
    return( (x[(length(cond_data) + 3):(2*length(cond_data) + 2) ] - x[1]*x[3:(length(cond_data) + 2) ] )/(x[3:(length(cond_data) + 2) ]^(x[2])))
  })
  
  z_sim = as.matrix(z[sample(1:nrow(z),nsim,replace=T),]) #randomly draw rows from the residual matrix 
  
  x_sim = as.matrix(rexp(nsim) + fit$u_cond_var) #simulate theoretically above high threshold for conditioning variable 
  
  y_sim = sapply(1:ncol(z_sim),function(i){ #calculate sample from non-conditioning variables 
    return( z_sim[,i]*(x_sim^fit$HT_parameter_estimates[2,i]) + x_sim*fit$HT_parameter_estimates[1,i]  )
  })
  
  sim_data = matrix(NA,ncol=ncol(fit$cond_data_laplace),nrow=nsim)
  
  sim_data[,fit$cond_var] = x_sim
  
  sim_data[,non_cond_var] = y_sim
  
  return(sim_data)
  
}
