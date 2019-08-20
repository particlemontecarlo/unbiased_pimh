#'@rdname CPF
#'@title Conditional Particle Filter
#'@description Runs a conditional particle filter, with or without ancestor sampling
#'@param nparticles number of particles
#'@param model a list representing a model, for instance as given by \code{\link{get_ar}}.
#'@param theta a parameter to give to the model functions
#'@param observations a matrix of observations of size datalength x dimension(observation)
#'@param ref_trajectory a reference trajectory, of size dimension(process) x datalength; if missing, runs a standard
#'particle filter.
#'@param with_as whether ancestor sampling should be used (TRUE/FALSE)
#'@return A new trajectory and optionally a trajectory averaged by the weights
#'@export
# perform the conditional CPF refresh for the coupling algorithm
CPF <- function(nparticles, model, theta, observations, ref_trajectory = NULL, with_as = FALSE,all_particles=FALSE){
  datalength <- nrow(observations)
  # create tree representation of the trajectories
  Tree <- new(TreeClass, nparticles, 10*nparticles*model$dimension, model$dimension)
  # initialization
  model_precomputed <- model$precompute(theta)
  init_randomness <- model$rinit_rand(nparticles, theta)
  xparticles <- model$rinit(nparticles, theta, init_randomness, model_precomputed)
  if (!is.null(ref_trajectory)){
    xparticles[,nparticles] <- ref_trajectory[,1]
  }
  Tree$init(xparticles)
  logz <- rep(NA,datalength)
  #  
  normweights <- rep(1/nparticles, nparticles)
  logw <- log(normweights)
  # if ancestor sampling, needs to keep the last generation of particles at hand
  if (with_as){
    x_last <- xparticles
  }
  # step t > 1
  for (time in 1:datalength){
    ancestors <- multinomial_resampling_n(normweights, nparticles)
    # if no observation or first time, no resampling
    if (time == 1 || (time > 1 && is.na(observations[time-1,1]))){
      ancestors <- 1:nparticles
    }
    xparticles <- xparticles[,ancestors]
    if (is.null(dim(xparticles))) xparticles <- matrix(xparticles, nrow = dimension)
    transition_randomness <- model$rtransition_rand(nparticles, theta)
    xparticles <- model$rtransition(xparticles, theta, time, transition_randomness, model_precomputed)
    if (!is.null(ref_trajectory)){
      xparticles[,nparticles] <- ref_trajectory[,time+1]
      if (with_as){
        # Ancestor sampling
        logm <- model$dtransition(ref_trajectory[,time+1], x_last, theta, time, model_precomputed)
        logm <- logw + logm
        w_as <- exp(logm - max(logm))
        w_as <- w_as / sum(w_as)
        ancestors[nparticles] <- systematic_resampling_n(w_as, 1, runif(1))
        x_last <- xparticles
      } else {
        ancestors[nparticles] <- nparticles
      }
    }
    #    
    logw <- model$dmeasurement(xparticles, theta, observations[time,], model_precomputed)
    maxlw <- max(logw)
    w <- exp(logw - maxlw)
    normweights <- w / sum(w)
    #
    Tree$update(xparticles, ancestors - 1)   
    logz[time] <- log(mean(exp(logw - maxlw)))+maxlw
  }
  
  new_trajectory <- Tree$get_path(systematic_resampling_n(normweights, 1, runif(1))-1)
  if(!all_particles){
    return(list(trajectory=new_trajectory,logz=sum(logz)))
  }else{
    path_arr <- sapply((1:nparticles)-1,Tree$get_path,simplify='array')
    exp_path <- pf_get_weighted_average(path_arr,normweights)
    stopifnot(sum(path_arr[1,1,]*normweights)==exp_path[1,1])
    return(list(trajectory=new_trajectory,exp_path=exp_path,logz=sum(logz)))
  }
}
# pf used in pmmh inference
#'@rdname pf_get_weighted_average
#'@title weight average of path
#'@description utility to get weighted average of particle paths
#'@export
pf_get_weighted_average <- function(arr_test,w){
  # takes an array (dimension, observations, particles ) and a vector
  # of weights. Returns array of sum_i ( (dimension , observations) *weight_i  )
  exp_path <- matrix(NA,dim(arr_test)[1],dim(arr_test)[2])
  for(d in 1:dim(arr_test)[1]){
    exp_path[d,] <- colSums(t(arr_test[d,,]) * w)
  }
  return(exp_path)
}



# pf used in pmmh inference
#'@rdname pf
#'@title bootstrap particle filter inference
#'@description runs a particle filter with specified model, model parameters and number of particles
#'
#'@param y observation set
#'@param theta model parameters
#'@param model initial and transition densities for HMM
#'@param nparticles number of particles
#'@export
pf <- function(y, theta, model, nparticles,ess_threshold=1,systematic=FALSE){
  
  nobservations <- nrow(y)
  dimension <- model$dimension
  Tree <- new(TreeClass, nparticles, 10*nparticles*dimension, dimension)
  
  # initial step
  randomness <- model$generate_randomness(nparticles, nobservations)
  xparticles <- model$rinit(nparticles, theta, randomness)
  Tree$init(t(xparticles))
  normweights <- rep(1/nparticles, nparticles)
  #
  loglik <- rep(0, nobservations)
  log_W <- rep(log(1/nparticles),nparticles)
  logw <- rep(0,nparticles)
  pf_means <- matrix(0, nrow = nobservations, ncol = dimension)
  
  for (time in 1:nobservations){
    xparticles <- model$rtransition(xparticles, theta, time, randomness,nparticles)
    logw <- model$dmeasurement(xparticles, theta, y[time,])
    
    log_weights <- logw + log_W
    maxlw <- max(log_weights)
    w <- exp(log_weights - maxlw)
    loglik[time] <- maxlw + log(sum(w))
    normweights <- w / sum(w)
    # filtering means
    pf_means[time,] <- apply(X = xparticles, MARGIN = 2, FUN = function(x) sum(normweights*x))
    #
    if(((1/sum(normweights**2))/nparticles)<ess_threshold){
      if(!systematic){
        ancestors <- sample(x = 1:nparticles, nparticles, replace = TRUE, prob = normweights)
      }else{
        ancestors <- systematic_resampling_n(normweights, nparticles, runif(1))
      }
      
      log_W <- rep(log(1/nparticles),nparticles)
      logw <- rep(0,nparticles)
    }else{
      ancestors = 1:nparticles
      log_W <- log(normweights)
    }
    xparticles <- xparticles[ancestors,,drop=F]
    Tree$update(t(xparticles), ancestors - 1)
    #
  }
  index_path <- sample(1:nparticles, 1, prob = normweights)
  path <- Tree$get_path(index_path - 1)
  return(list(pf_means = pf_means, loglik = sum(loglik), path = t(path),loglik_t =cumsum(loglik)))
}


# smc sampler used in  inference
#'@rdname smc_sampler
#'@title SMC sampler
#'@description runs an SMC sampler with specified model, model parameters and number of particles
#'
#'@param theta model parameters
#'@param model initial and transition densities for HMM
#'@param nparticles number of particles
#'@export
smc_sampler <- function(nparticles, model,theta,systematic=FALSE,h_func=NULL){
  
  ess_threshold <- 1
  ntransitions <- model$ntransitions
  
  # set up arrays
  lognormalising <- rep(0, ntransitions+1)
  
  # # # initial step
  # sample
  xparticles_new <- model$rinit(nparticles, theta)
  # initial incremental weights
  logw <- model$inc_init(nparticles,xparticles_new, theta)
  
  log_W <- rep(log(1/nparticles),nparticles)
  log_weights <- logw + log_W
  maxlw <- max(log_weights)
  w <- exp(log_weights - maxlw)
  lognormalising[1] <- maxlw + log(sum(w))
  normweights <- w / sum(w)
  ancestors = 1:nparticles
  log_W <- log(normweights)
  
  
  # iterate through the kernels
  for (time in 1:ntransitions){
    # if((time%%10)==0){
    #   print(time/ntransitions)
    # }
    xparticles <- xparticles_new[,ancestors,drop=F]
    
    # sample
    xparticles_new <- model$rtransition(xparticles, theta, time ,nparticles)
    # get new incremental weights
    logw <- model$inc_w(xparticles,xparticles_new, theta, time)
    
    # get new weights
    log_weights <- logw + log_W
    maxlw <- max(log_weights)
    w <- exp(log_weights - maxlw)
    lognormalising[time+1] <- maxlw + log(sum(w))
    normweights <- w / sum(w)
    
    
    # resample
    if(((1/sum(normweights**2))/nparticles)<ess_threshold){
      if(!systematic){
        ancestors <- sample(x = 1:nparticles, nparticles, replace = TRUE, prob = normweights)
      }else{
        ancestors <- systematic_resampling_n(normweights, nparticles, runif(1))
      }
      log_W <- rep(log(1/nparticles),nparticles)
      logw <- rep(0,nparticles)
    }else{
      ancestors = 1:nparticles
      log_W <- log(normweights)
    }
    
  }
  
  # if h_func is not specified randomly sample a trajectory - otherwise reuse all the particles with the particular expectation
  if(is.null(h_func)){
    ret_trajectory <- t(xparticles_new[,sample(x = 1:nparticles, 1, prob = normweights),drop=F])
  }else{
    h_vals <- apply(xparticles_new,2,FUN=h_func)
    ret_trajectory <- matrix(sum(h_vals * normweights))
    stopifnot(all(dim(ret_trajectory)==c(1,1)))
  }
  return(list(trajectory = ret_trajectory, logz = sum(lognormalising), loglik_t =cumsum(lognormalising)))
}





