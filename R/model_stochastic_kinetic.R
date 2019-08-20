#'@rdname get_stochastic_kinetic
#'@title Stochastic kinetic model as in Golightly 2005
#'@description This function returns a list with objects such as
#'* rinit, rinit_rand to sample from the initial distribution
#'* rtransition, rtransition_rand to sample from the transition
#'* dtransition to evaluate the transition density
#'* dmeasurement to evaluate the measurement density
#'* dimension, which represents the dimension of the latent process
#'@return A list
#'@export

### model & gen data
get_stochastic_kinetic <- function() 
{
  
  dimension <- 4
  dimension_y <- 2
  k <- 5
  s1 <- c(0,0,1,0,0,0,-1,0)
  s2 <- c(0,0,0,1,-2,2,0,-1)
  s3 <- c(-1,1,0,0,1,-1,0,0)
  s4 <- c(-1,1,0,0,0,0,0,0)
  S <- rbind(s1,s2,s3,s4)
  c_vec <- c(0.1,0.7,0.35,0.2,0.1,0.9,0.3,0.1)
  
  
  get_h_arr <- function(xparticles){
    # this function returns the hazard functions for the eight possible events
    # for each particle
    h_arr <- matrix(NA,8,dim(xparticles)[2])
    h_arr[1,] <- c_vec[1]*xparticles[4,]*xparticles[3,]
    h_arr[2,] <- c_vec[2]*(k-xparticles[4,])
    h_arr[3,] <- c_vec[3]*xparticles[4,]
    h_arr[4,] <- c_vec[4]*xparticles[1,]
    h_arr[5,] <- c_vec[5]*xparticles[2,]*(xparticles[2,]-1)/2
    h_arr[6,] <- c_vec[6]*xparticles[3,]
    h_arr[7,] <- c_vec[7]*xparticles[1,]
    h_arr[8,] <- c_vec[8]*xparticles[2,]
    return(h_arr)
  }
  
  
  rinit <- function(nparticles, theta, rand, precomputed, ...) {
    X_0 <- c(8,8,8,5)
    init_mat <- replicate(nparticles,X_0)
    
    return(init_mat)
  }
  
  
  rinit_rand <- function(nparticles, theta) {
    return(NA)
  }
  
  
  # this function is used to sample the type of event for each particle
  # as well as update the particles
  update_alive_particles <- function(h_arr_alive,h_0_alive,x_particles_alive,n_alive){
    # for the alive particles we need to sample the event type with a certain
    # probability given by the hazard functions, then update according to
    # the stoichiometry matrix
    
    # I think sampling the event type is probably best by explicitly constructing the
    # cumulative distribution function
    u_rnd <- runif(n_alive)*h_0_alive
    cdf_h <- colCumsums(h_arr_alive)
    event_type <- rowSums(t(cdf_h)<u_rnd)+1
    update_mat <- S[,event_type]
    return(x_particles_alive + update_mat)
    
  }
  
  rtransition <- function(xparticles, theta, time, rand, precomputed, 
                          ...) {
    stopifnot(is.null(rand))
    # need to recursively sample inter-arrival times


    dx <- dim(xparticles)[1]
    nparticles <- dim(xparticles)[2]
    xparticles_updated <- matrix(NA,4,nparticles)
    
    xparticles_alive <- xparticles
    n_alive <- nparticles
    alive_particles <- rep(T,nparticles)
    global_alive_mask <- rep(T,nparticles)
    alive_mask <- rep(T,nparticles)
    t_particles_alive <- rep(0,nparticles)
    
    
    # we need to recursively
    # 1. sample new arrival times
    # 2. work out which particles are still alive
    # 3. for those that are alive update the state at the new arrival times
    counter <- 0
    while(n_alive>0){
      
      counter <- counter + 1
      # 1. get the hazard function and sample interarrival times
      h_arr_alive <- get_h_arr(xparticles_alive)
      h_0_alive <- colSums(h_arr_alive)
      delta_t <- rexp(n_alive,h_0_alive)
      
      # 2. work out which particles are alive at the new times
      # we need the particles that have just been killed to record
      # their value
      t_particles_alive <- t_particles_alive +  delta_t
      alive_particles <- t_particles_alive<0.1
      
      global_alive_mask_updated <- global_alive_mask
      global_alive_mask_updated[global_alive_mask] <- t_particles_alive<0.1
      xparticles_updated[,(!global_alive_mask_updated) & global_alive_mask] <- xparticles_alive[,!alive_particles]
      global_alive_mask <- global_alive_mask_updated
      n_alive <- sum(global_alive_mask)
      
      h_arr_alive <- h_arr_alive[,alive_particles,drop=F]
      
      h_0_alive <- h_0_alive[alive_particles]
      xparticles_alive <- xparticles_alive[,alive_particles,drop=F]
      t_particles_alive <- t_particles_alive[alive_particles]
      
      if(n_alive>0){
        # 3. for those that are alive sample the event type and adjust the state of the particles
        xparticles_alive <- update_alive_particles(h_arr_alive,
                                                   h_0_alive,
                                                   xparticles_alive,
                                                   n_alive)  
      }
    }
    
    
    return(xparticles_updated)
  }
  
  rtransition_rand <- function(nparticles, theta) {
    return(NULL)
  }
  
  dtransition <- function(next_x, xparticles, theta, time, 
                          precomputed, ...) {
    stop('no dtransition available for stochasti_kinetic model')
    return(dmvnorm_transpose_cholesky(precomputed$A %*% xparticles, 
                                      next_x, precomputed$di))
  }
  
  dmeasurement <- function(xparticles, theta, observation, 
                           precomputed, ...) {
    A_obs <- theta$A_obs
    s2_obs <- theta$s2_obs
    mean_obs <- A_obs %*% xparticles
    l_obs <- colSums(dnorm(mean_obs,mean=observation,sd=sqrt(s2_obs),log=T))
    
    return(l_obs)
  }
  precompute <- function(theta) {
    return(list())
  }
  stochastic_kinetic_model <- list(rinit = rinit, rinit_rand = rinit_rand, 
                                   rtransition = rtransition, rtransition_rand = rtransition_rand, 
                                   dtransition = dtransition, dmeasurement = dmeasurement, 
                                   precompute = precompute, dimension = dimension,dimension_y=dimension_y)
  return(stochastic_kinetic_model)
}


