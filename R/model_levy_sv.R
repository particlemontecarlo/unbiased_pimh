#'@rdname get_levy_sv
#'@title Levy-driven stochastic volatility model as in Barndorff-Nielsen and Shephard (2001)
#'@description This function returns a list with objects such as
#'* rinit, rinit_rand to sample from the initial distribution
#'* rtransition, rtransition_rand to sample from the transition
#'* dtransition to evaluate the transition density
#'* dmeasurement to evaluate the measurement density
#'* dimension, which represents the dimension of the latent process
#'@return A list
#'@export
# SV model
# theta = (mu, beta, xi, w2, lambda)
# transformed so all are in R...
get_levy_sv <- function(){
  ## Levy-SV
  
  rinit = function(nparticles, theta, rand, ...){
    # xi = theta[3]
    # w2 = theta[4]
    # lambda = theta[5]
    # Note: to avoid numerical issues, we artificially set the variance vt to machine epsilon
    # instead of 0 whenever needed
    
    Xs = CoupledPIMH:::rinitial_SVLevy_cpp(nparticles, 1, theta[3], theta[4], theta[5])
    if (all(Xs[,1]<.Machine$double.eps)) {Xs[,1] = rep(.Machine$double.eps, nparticles)}
    if (all(Xs[,2]<.Machine$double.eps)) {Xs[,2] = rep(.Machine$double.eps, nparticles)}
    return (t(Xs))
  }
  
  
  
  rinit_rand <- function(nparticles, theta){
    return(NULL)
  }
  
  rtransition = function(xparticles, theta, t, rand, ...){
    # xi = theta[3]
    # w2 = theta[4]
    # lambda = theta[5]
    # Note: to avoid numerical issues, we artificially set the variance vt to machine epsilon
    # instead of 0 whenever needed
    xparticles <- t(xparticles)
    nparticles <- dim((xparticles))[1]
    
    
    new_Xs = CoupledPIMH:::rtransition_SVLevy_cpp(xparticles, t-1, t, theta[3], theta[4], theta[5])
    if (all(new_Xs[,1]<.Machine$double.eps)) {
      new_Xs[,1] = rep(.Machine$double.eps, nparticles)
      stop('this should not be happening')
    }
    if (all(new_Xs[,2]<.Machine$double.eps)) {
      new_Xs[,2] = rep(.Machine$double.eps, nparticles)
      stop('this should not be happening')
    }
    return (t(new_Xs))
  }
  
  
  rtransition_rand <- function(nparticles, theta){
    return(NULL)
  }
  
  dmeasurement <- function(xparticles, theta, observation, ...) {
    
    return (CoupledPIMH:::dobs_SVLevy_cpp(matrix(observation, ncol = 1), t(xparticles), theta[1], theta[2], TRUE)[1,])
  }
  
  precompute <- function(...){
    return(list())
  }
  
  robs = function(Xt,theta){

    mu = theta[1]
    beta = theta[2]
    return (rnorm(1, mu + beta*Xt[1], sd=sqrt(Xt[1])))
  }
  
  #
  sv <- list(rinit = rinit, rinit_rand = rinit_rand, rtransition = rtransition,
             rtransition_rand = rtransition_rand,robs=robs,
             dmeasurement = dmeasurement, precompute = precompute, dimension = 2)
  return(sv)
}
# get_model_SVLevy_singlefactor <- function(){
#   model = list()
#   # Dimension of parameter, observations, and possibly latent states (int)
#   model$dimtheta = 5
#   model$dimY = 1
#   model$dimX = 2
#   model$dimension <- 2
# 
#   #----------------------------------------------------------------------------------------------------
#   #----------------------------------------------------------------------------------------------------
#   # Note: if no likelihood nor predictive is provided, the method will be SMC2, which requires
#   # specifying the transition kernel and the observation density
#   #----------------------------------------------------------------------------------------------------
#   #----------------------------------------------------------------------------------------------------
#   # Sampler from the initial distribution of the latent states
#   # inputs: theta (single vector), Nx (int)
#   # outputs: matrix (dimX by Nx) of latent states
#   model$generate_randomness <- function(nparticles, nobservations){
#     return(NULL)
#   }
# 
#   model$rinitial = function(Nx, theta, rand, ...){
#     # xi = theta[3]
#     # w2 = theta[4]
#     # lambda = theta[5]
#     # Note: to avoid numerical issues, we artificially set the variance vt to machine epsilon
#     # instead of 0 whenever needed
# 
#     Xs = rinitial_SVLevy_cpp(Nx, 1, theta[3], theta[4], theta[5])
#     if (all(Xs[,1]<.Machine$double.eps)) {Xs[,1] = rep(.Machine$double.eps, Nx)}
#     if (all(Xs[,2]<.Machine$double.eps)) {Xs[,2] = rep(.Machine$double.eps, Nx)}
#     return (Xs)
#   }
#   # Sampler from the transition distribution of the latent states
#   # inputs: current states Xs at time (t-1) (dimX by Nx matrix), time t (int), theta (single vector)
#   # outputs: updated states (dimX by Nx)
#   model$rtransition = function(Xs, theta, t, rand, nparticles, ...){
#     # xi = theta[3]
#     # w2 = theta[4]
#     # lambda = theta[5]
#     # Note: to avoid numerical issues, we artificially set the variance vt to machine epsilon
#     # instead of 0 whenever needed
# 
#     new_Xs = CoupledPIMH:::rtransition_SVLevy_cpp(Xs, t-1, t, theta[3], theta[4], theta[5])
#     if (all(new_Xs[,1]<.Machine$double.eps)) {new_Xs[,1] = rep(.Machine$double.eps, nparticles)}
#     if (all(new_Xs[,2]<.Machine$double.eps)) {new_Xs[,2] = rep(.Machine$double.eps, nparticles)}
#     return (new_Xs)
#   }
#   # observation density
#   # inputs: single observation Yt (dimY by 1), states Xts (dimX by Nx), time t, theta (single vector), log (TRUE by default)
#   # outputs: observation (log)-densities ("vectorized" with respect to the states Xt)
#   # model$dobs = function(Xts, theta, Yt){
#   #   # mu = theta[1]
#   #   # beta = theta[2]
#   #   return (debiasedpmcmc:::dobs_SVLevy_cpp(Yt, Xts, theta[1], theta[2], TRUE)[1,])
#   # }
#   model$dmeasurement <- function(xparticles, theta, observation, ...) {
# 
#     return (CoupledPIMH:::dobs_SVLevy_cpp(matrix(observation, ncol = 1), xparticles, theta[1], theta[2], TRUE)[1,])
#   }
#   #
#   model$robs = function(Xt,theta){
# 
#     mu = theta[1]
#     beta = theta[2]
#     return (rnorm(1, mu + beta*Xt[1], sd=sqrt(Xt[1])))
#   }
#   return(model)
# }
# 
# 



