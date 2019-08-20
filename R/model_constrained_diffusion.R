#'@rdname get_constrained_diffusion
#'@title Constrained sine diffusion
#'@description This function returns a list with objects such as
#'* rinit, rinit_rand to sample from the initial distribution
#'* rtransition, rtransition_rand to sample from the transition
#'* dtransition to evaluate the transition density
#'* dmeasurement to evaluate the measurement density
#'* dimension, which represents the dimension of the latent process
#'@return A list
#'@export
### model & gen data
get_constrained_diffusion <- function(time_T,T_final)
{
  dimension <- 1
  h <- T_final/time_T
  
  rinit <- function(nparticles, theta, rand, precomputed, ...) {
    xparticles <- matrix(0, nrow = dimension,ncol=nparticles)
    xparticles[dimension,] <- theta[2]*(2*(pnorm(rand)-1))
    return(xparticles)
  }
  rinit_rand <- function(nparticles, theta) {
    return(rnorm(nparticles))
  }
  rtransition <- function(xparticles, theta, time, rand, precomputed, 
                          ...) {
    xparticles[1,] <-  xparticles[dimension,] + theta[3]*(sin(xparticles[dimension,]-theta[1]))*h + sqrt(h)*rand[1:nparticles]
    if(dimension>1){
      for(d in 2:dimension){
        xparticles[d,] <-  xparticles[d-1,] + theta[3]*(sin(xparticles[d-1,]-theta[1]))*h + sqrt(h)*rand[(1+nparticles*(d-1)):(nparticles*d)]
      }
    }
    return(xparticles)
  }
  rtransition_rand <- function(nparticles, theta) {
    return(rnorm(nparticles * dimension))
  }
  dtransition <- function(next_x, xparticles, theta, time, 
                          precomputed, ...) {
    
    if(dimension>1){
      like_est_delta_0 <- 0
      for(d in 2:dimension){
        next_mu <- xparticles[d-1,]+ theta[3]*(sin(xparticles[d-1,]-theta[1]))*h
        like_est_delta_0 <- like_est_delta_0 + dnorm(xparticles[d,],next_mu,sd=sqrt(h),log=T)
      }
      
      next_mu <- xparticles[dimension,]+ theta[3]*(sin(xparticles[dimension,]-theta[1]))*h #+ sqrt(h)*rand[1:nparticles]
      like_est_delta <- dnorm(next_x[1],next_mu,sd=sqrt(h),log=T)
      
      # for(d in 2:dimension){
      #   next_mu <- next_x[d-1]+ theta[3]*(sin(next_x[d-1]-theta[1]))*h
      #   like_est_delta <- like_est_delta + dnorm(next_x[d],next_mu,sd=sqrt(h),log=T)
      # }
      
      like_ests <- like_est_delta_0 + like_est_delta
    }else{
      next_mu <- xparticles[1,]+ theta[3]*(sin(xparticles[1,]-theta[1]))*h
      like_ests <- dnorm(next_x[1],next_mu,sd=sqrt(h),log=T)
    }
    
    return(like_ests)
  }
  dmeasurement <- function(xparticles, theta, observation, 
                           precomputed, ...) {
    
    like_est <- log(1*(abs(xparticles[dimension,])<theta[2]))#dnorm(observation,xparticles[dimension,],theta[2],log=T)
    return(like_est)
  }
  precompute <- function(theta) {
    return(list())
  }
  constrained_diffusion_model <- list(rinit = rinit, rinit_rand = rinit_rand, 
                               rtransition = rtransition, rtransition_rand = rtransition_rand, 
                               dtransition = dtransition, dmeasurement = dmeasurement, 
                               precompute = precompute, dimension = dimension)
  return(constrained_diffusion_model)
}

