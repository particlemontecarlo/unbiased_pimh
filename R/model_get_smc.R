#'@rdname get_smc
#'@title SMC sampler
#'@description This function returns a list with objects such as
#'* rinit, rinit_rand to sample from the initial distribution
#'* rtransition, rtransition_rand to sample from the transition
#'* dtransition to evaluate the transition density
#'* dmeasurement to evaluate the measurement density
#'* dimension, which represents the dimension of the latent process
#'@return A list
#'@export
get_smc <- function(){
  ## SMC sampler model
  ntransitions <- 2e2
  
  
  gen_data <- function(nobs=NULL){
    if(is.null(nobs)){
      nobs <- 1e2  
    }
    stopifnot(nobs%%dimension==0)
    y_vec <- rep(NA,nobs)
    for(i in 1:dimension){
      n_strat <- nobs/dimension
      
      y_vec[((i-1)*n_strat+1):(i*n_strat)] <- rnorm(n_strat,theta$xstar[i],theta$ss)
    }
    return(y_vec)
  }
  
  logsumexp_pairwse <- function(x1,x2){
    max1 <- pmax(x1,x2)
    min1 <- pmin(x1,x2)
    res <- max1 + log(1+exp(min1-max1))
    return(res)
  }
  
  
  loglik <- function(xparticles,theta){
    
    f_outer <- function(x,y) dnorm(x,y,s=theta$ss,log=T)
    
    # this creates of size n x N x d
    big_array <- outer(yobs,t(xparticles),FUN=f_outer)
    
    # will then want to do log sum exp on each dimension
    # do the pairwise log sum exp
    log_l <- big_array[,,1]
    if(dimension>1){
      for(i in 2:dimension){
        log_li <- big_array[,,i]
        log_l <-logsumexp_pairwse(log_l,log_li)
      }  
      
    }
    
    log_ll <- colSums(log_l)
    return(log_ll)
  }
  
  test_loglik <- function(){
    dimension <- 2
    xparticles <- matrix(1:6,dimension)/1
    yobs <- 1:4/1
    
    ll <- log(colSums(dnorm(xparticles,yobs[1],s=theta$ss)))
    for(i in 2:length(yobs)){
      ll <- ll + log(colSums(dnorm(xparticles,yobs[i],s=theta$ss)))
    }
    loglik(xparticles,theta)
  }
  
  target_t <- function(xparticles,theta,lambda){
    
    ll <- loglik(xparticles,theta)
    
    prior_mask <- ((xparticles<=10) & (xparticles>=-10))
    prior_vals <- colSums(!prior_mask)
    prior_vals[prior_vals>0] <- -Inf
    
    ll_target <- lambda*ll + prior_vals
    
    return(ll_target)
  }
  
  # model:  pi0 = uniform distribution on d-dimensional hypercube
  #         pi = sum (1/4) N(mu_i,s)
  
  # sample particles
  rinit <- function(nparticles, theta, rand, ...){
    return((matrix(runif(nparticles*dimension),nrow=dimension)*20)-10)
  }
  rtransition <- function(xparticles, theta, time ,nparticles, ...){
    xparticles_prop <- xparticles + matrix(rnorm(prod(dim(xparticles))),nrow=dim(xparticles)[1])
    lambda_t <- (time/ntransitions)^2
    
    target_curr <- target_t(xparticles,theta,lambda_t)
    target_prop <- target_t(xparticles_prop,theta,lambda_t)
    log_u <- log(runif(nparticles))
    accept_mask <- (target_prop-target_curr)>log_u
    
    xparticles_new <- xparticles
    xparticles_new[,accept_mask] <- xparticles_prop[,accept_mask]
    return(xparticles_new)
  }
  
  # incremental weight functions
  inc_init <- function(nparticles,xparticles, theta, ...){
    return(rep(0,nparticles))
  }
  inc_w <- function(xparticles,xparticles_new, theta, time, ...){
    lambda_t <- (time/ntransitions)^2
    lambda_tm1 <- ((time-1)/ntransitions)^2
    
    log_w <- target_t(xparticles,theta,lambda_t)-target_t(xparticles,theta,lambda_tm1)
    
    return(log_w)
  }
  
  
  smc <- list(ntransitions=ntransitions,
              rinit=rinit,
              rtransition=rtransition,
              inc_init=inc_init,
              inc_w=inc_w,
              gen_data=gen_data,
              dimension=dimension)
  return(smc)
}
