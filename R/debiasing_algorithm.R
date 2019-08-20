# do things properly, recycling the code for unbiased mcmc
#'@rdname coupled_ccpf_chains
#'@title Coupled conditional particle filter chains
#'@description Sample two conditional particle filter chains, each following \code{single_kernel} marginally,
#' and \code{coupled_kernel} jointly, until min(max(tau, K), max_iterations), where tau
#' is the first time at which the two chains meet (i.e. take the same value exactly).
#' Or more precisely, they meet with a delay of one, i.e. X_t = Y_{t-1}. The chains
#' are initialized from the distribution provided in \code{rinit}.
#'
#'  See \code{\link{get_hmc_kernel}}
#' for an example of function returning the appropriate kernels.
#'
#'@param single_kernel function taking a state (in a vector) and an iteration, and returning
#' a list with a key named \code{chain_state} and containing the next state.
#'@param coupled_kernel function taking two states (in two vectors) and an iteration,
#'and returning a list with keys \code{chain_state1} and \code{chain_state2}.
#'@param rinit function taking no arguments are returning an initial state for a Markov chain.
#'@param K number of iterations desired (will be proportional to the computing cost if meeting occurs before \code{K},
#' default to 1).
#'@param max_iterations number of iterations at which the function stops if it is still running  (default to Inf).
#'@param preallocate  expected number of iterations, used to pre-allocate memory (default to 10).
#'@export
coupled_ccpf_chains <- function(single_kernel, coupled_kernel, pimh_init,N,K = 1, max_iterations = Inf, preallocate = 10,verbose=FALSE){
  
  initial_conditions <- pimh_init(N)
  chain_state1 <- initial_conditions$chain_state1
  chain_state2 <- initial_conditions$chain_state2
  log_pdf_state1 <- initial_conditions$log_pdf_state1
  log_pdf_state2 <- initial_conditions$log_pdf_state2
  
  
  p <- length(chain_state1[1,])
  samples1 <- matrix(nrow = K+preallocate+1, ncol = p)
  nrowsamples1 <- K+preallocate+1
  samples2 <- matrix(nrow = K+preallocate, ncol = p)
  log_pdf1 <- matrix(nrow = K+preallocate+1, ncol = 1)
  log_pdf2 <- matrix(nrow = K+preallocate, ncol = 1)
  
  samples1[1,] <- chain_state1[1,]
  samples2[1,] <- chain_state2[1,]
  log_pdf1[1,] <- log_pdf_state1
  log_pdf2[1,] <- log_pdf_state2
  
  
  current_nsamples1 <- 1
  iter <- 1
  mh_step <- single_kernel(chain_state1,log_pdf_state1, iter,N)
  chain_state1 <- mh_step$chain_state
  log_pdf_state1 <- mh_step$log_pdf_state
  
  
  current_nsamples1 <- current_nsamples1 + 1
  samples1[current_nsamples1,] <- chain_state1[1,]
  log_pdf1[current_nsamples1,] <- log_pdf_state1
  
  meet <- FALSE
  finished <- FALSE
  meetingtime <- Inf
  while (!finished && iter < max_iterations){
    if(verbose){
      print(iter)
    }
    iter <- iter + 1
    if (meet){
      mh_step <- single_kernel(chain_state1,log_pdf_state1, iter,N)
      chain_state1 = mh_step$chain_state
      log_pdf_state1 = mh_step$log_pdf_state
      chain_state2 = mh_step$chain_state
      log_pdf_state2 = mh_step$log_pdf_state
    } else {
      mh_step <- coupled_kernel(chain_state1, chain_state2,log_pdf_state1,log_pdf_state2,iter,N)
      chain_state1 <- mh_step$chain_state1
      chain_state2 <- mh_step$chain_state2
      log_pdf_state1 <- mh_step$log_pdf_state1
      log_pdf_state2 <- mh_step$log_pdf_state2
      
      
      
      if (all(chain_state1 == chain_state2) && !meet){
        # recording meeting time tau
        meet <- TRUE
        meetingtime <- iter
      }
    }
    if ((current_nsamples1+1) > nrowsamples1){
      # print('increase nrow')
      new_rows <- nrowsamples1-1
      nrowsamples1 <- nrowsamples1 + new_rows
      samples1 <- rbind(samples1, matrix(NA, nrow = new_rows, ncol = p))
      samples2 <- rbind(samples2, matrix(NA, nrow = new_rows, ncol = p))
      log_pdf1 <- rbind(log_pdf1, matrix(NA, nrow = new_rows, ncol = ncol(log_pdf1)))
      log_pdf2 <- rbind(log_pdf2, matrix(NA, nrow = new_rows, ncol = ncol(log_pdf2)))
    }
    samples1[current_nsamples1+1,] <- chain_state1[1,]
    samples2[current_nsamples1,] <- chain_state2[1,]
    log_pdf1[current_nsamples1+1] <- log_pdf_state1
    log_pdf2[current_nsamples1] <- log_pdf_state2
    
    current_nsamples1 <- current_nsamples1 + 1
    # stop after max(K, tau) steps
    if (iter >= max(meetingtime, K)){
      finished <- TRUE
    }
  }
  samples1 <- samples1[1:current_nsamples1,,drop=F]
  samples2 <- samples2[1:(current_nsamples1-1),,drop=F]
  log_pdf1 <- log_pdf1[1:current_nsamples1,,drop=F]
  log_pdf2 <- log_pdf2[1:(current_nsamples1-1),,drop=F]
  return(list(samples1 = samples1, samples2 = samples2, log_pdf1=log_pdf1, log_pdf2=log_pdf2,
              meetingtime = meetingtime, iteration = iter, finished = finished))
}






# do things properly, recycling the code for unbiased mcmc
#'@rdname coupled_pimh_chains
#'@title Coupled PIMH chains
#'@description Sample two PIMH chains, each following \code{single_kernel} marginally,
#' and \code{coupled_kernel} jointly, until min(max(tau, K), max_iterations), where tau
#' is the first time at which the two chains meet (i.e. take the same value exactly).
#' Or more precisely, they meet with a delay of one, i.e. X_t = Y_{t-1}. The chains
#' are initialized from the distribution provided in \code{rinit}.
#'
#'  See \code{\link{get_hmc_kernel}}
#' for an example of function returning the appropriate kernels.
#'
#'@param single_kernel function taking a state (in a vector) and an iteration, and returning
#' a list with a key named \code{chain_state} and containing the next state.
#'@param coupled_kernel function taking two states (in two vectors) and an iteration,
#'and returning a list with keys \code{chain_state1} and \code{chain_state2}.
#'@param rinit function taking no arguments are returning an initial state for a Markov chain.
#'@param K number of iterations desired (will be proportional to the computing cost if meeting occurs before \code{K},
#' default to 1).
#'@param max_iterations number of iterations at which the function stops if it is still running  (default to Inf).
#'@param preallocate  expected number of iterations, used to pre-allocate memory (default to 10).
#'@export
coupled_pimh_chains <- function(single_kernel, coupled_kernel, pimh_init,N,K = 1, max_iterations = Inf, preallocate = 10,verbose=FALSE){
  
  
  initial_conditions <- pimh_init(N)
  chain_state1 <- initial_conditions$chain_state1
  log_pdf_state1 <- initial_conditions$log_pdf_state1
  
  
  p <- length(chain_state1[1,])
  nrowsamples1 <- K+preallocate+1
  
  samples1 <- matrix(nrow = K+preallocate+1, ncol = p)
  samples2 <- matrix(nrow = K+preallocate, ncol = p)
  log_pdf1 <- matrix(nrow = K+preallocate+1, ncol = 1)
  log_pdf2 <- matrix(nrow = K+preallocate, ncol = 1)
  
  samples1[1,] <- chain_state1[1,]
  log_pdf1[1,] <- log_pdf_state1
  
  
  chain_state2 <- NULL
  log_pdf_state2 <- -Inf
  
  # the function works out whether the expected path is returned and allocated new memory accordingly
  all_particles <- F
  if(!is.null(initial_conditions$exp_path1)){
    all_particles <- T
    exp_path1 <- initial_conditions$exp_path1
    exp_samples1 <- matrix(nrow = K+preallocate+1, ncol = p)
    exp_samples2 <- matrix(nrow = K+preallocate, ncol = p)
    exp_samples1[1,] <- exp_path1[1,]
    exp_path2 <- NULL
  }
  
  
  iter <- 0
  
  meet <- FALSE
  finished <- FALSE
  meetingtime <- Inf
  while (!finished && iter < max_iterations){
    if(verbose){
      print(iter)
    }
    iter <- iter + 1
    if (meet){
      if(!all_particles){
        mh_step <- single_kernel(chain_state1,log_pdf_state1, iter,N)  
      }else{
        mh_step <- single_kernel(chain_state1,log_pdf_state1,exp_path1, iter,N)  
      }
      
      chain_state1 = mh_step$chain_state
      log_pdf_state1 = mh_step$log_pdf_state
      chain_state2 = mh_step$chain_state
      log_pdf_state2 = mh_step$log_pdf_state
      if(all_particles){
        exp_path1 = mh_step$exp_path
        exp_path2 = mh_step$exp_path
      }
    } else {
      if(!all_particles){
        mh_step <- coupled_kernel(chain_state1, chain_state2,log_pdf_state1,log_pdf_state2,iter,N)  
      }else{
        mh_step <- coupled_kernel(chain_state1, chain_state2,log_pdf_state1,log_pdf_state2,exp_path1,exp_path2,iter,N)
      }
      
      chain_state1 <- mh_step$chain_state1
      chain_state2 <- mh_step$chain_state2
      log_pdf_state1 <- mh_step$log_pdf_state1
      log_pdf_state2 <- mh_step$log_pdf_state2
      if(all_particles){
        exp_path1 = mh_step$exp_path1
        exp_path2 = mh_step$exp_path2
      }
      
      if ((log_pdf_state1==log_pdf_state2) && all(chain_state1 == chain_state2) && !meet){
        # recording meeting time tau
        meet <- TRUE
        meetingtime <- iter
      }
    }
    if ((iter+1) > nrowsamples1){
      # print('increase nrow')
      new_rows <- nrowsamples1-1
      nrowsamples1 <- nrowsamples1 + new_rows
      samples1 <- rbind(samples1, matrix(NA, nrow = new_rows, ncol = p))
      samples2 <- rbind(samples2, matrix(NA, nrow = new_rows, ncol = p))
      log_pdf1 <- rbind(log_pdf1, matrix(NA, nrow = new_rows, ncol = ncol(log_pdf1)))
      log_pdf2 <- rbind(log_pdf2, matrix(NA, nrow = new_rows, ncol = ncol(log_pdf2)))
      if(all_particles){
        exp_samples1 <- rbind(exp_samples1, matrix(NA, nrow = new_rows, ncol = p))
        exp_samples2 <- rbind(exp_samples2, matrix(NA, nrow = new_rows, ncol = p))
      }
    }
    samples1[iter+1,] <- chain_state1[1,]
    samples2[iter,] <- chain_state2[1,]
    log_pdf1[iter+1] <- log_pdf_state1
    log_pdf2[iter] <- log_pdf_state2
    if(all_particles){
      exp_samples1[iter+1,] <- exp_path1[1,]
      exp_samples2[iter,] <- exp_path2[1,]
    }
    
    # stop after max(K, tau) steps
    if (iter >= max(meetingtime, K)){
      finished <- TRUE
    }
  }
  samples1 <- samples1[1:(iter+1),,drop=F]
  samples2 <- samples2[1:(iter),,drop=F]
  log_pdf1 <- log_pdf1[1:(iter+1),,drop=F]
  log_pdf2 <- log_pdf2[1:(iter),,drop=F]
  if(!all_particles){
    return(list(samples1 = samples1, samples2 = samples2, log_pdf1=log_pdf1, log_pdf2=log_pdf2,
                meetingtime = meetingtime, iteration = iter, finished = finished))
  }else{
    exp_samples1 <- exp_samples1[1:(iter+1),,drop=F]
    exp_samples2 <- exp_samples2[1:(iter),,drop=F]
    
    return(list(samples1 = samples1, samples2 = samples2, log_pdf1=log_pdf1, log_pdf2=log_pdf2,
                exp_samples1=exp_samples1,
                exp_samples2=exp_samples2,
                meetingtime = meetingtime, iteration = iter, finished = finished))
  }
}








# from h_bar --------------------------------------------------------------

#'@rdname H_bar
#'@title Compute unbiased estimators from coupled chains
#'@description Compute the proposed unbiased estimators, for each of the element
#'in the list 'c_chains'. The integral of interest is that of the function h,
#'which can be multivariate. The estimator uses the variance reduction technique
#'whereby the estimator is the MCMC average between times k and K, with probability
#'going to one as k increases.
#'@export
H_bar <- function(c_chains, h = function(x) x, k = 0, K = 1){
  maxiter <- c_chains$iteration
  if (k > maxiter){
    print("error: k has to be less than the horizon of the coupled chains")
    return(NULL)
  }
  if (K > maxiter){
    print("error: K has to be less than the horizon of the coupled chains")
    return(NULL)
  }
  # test the dimension of h(X)
  p <- length(h(c_chains$samples1[1,]))
  h_of_chain <- apply(X = c_chains$samples1[(k+1):(K+1),,drop=F], MARGIN = 1, FUN = h)
  if (is.null(dim(h_of_chain))){
    h_of_chain <- matrix(h_of_chain, ncol = 1)
  } else {
    h_of_chain <- t(h_of_chain)
  }
  H_bar <- apply(X = h_of_chain, MARGIN = 2, sum)
  if (c_chains$meetingtime <= k + 1){
    # nothing else to add
  } else {
    deltas <- matrix(0, nrow = maxiter - k + 1, ncol = p)
    deltas_term <- rep(0, p)
    for (t in k:min(maxiter-1, c_chains$meetingtime-1)){ # t is as in the report, where the chains start at t=0
      coefficient <- min(t - k + 1, K - k + 1)
      delta_tp1 <- h(c_chains$samples1[t + 1 + 1,]) - h(c_chains$samples2[t+1,]) # the +1's are because R starts indexing at 1
      deltas_term <- deltas_term + coefficient * delta_tp1
    }
    H_bar <- H_bar + deltas_term
  }
  return(H_bar / (K - k + 1))
}




#'@rdname BC
#'@title Compute the bias correction for the unbiased estimators
#'@description Compute the bias correction, for each of the element
#'in the list 'c_chains'. The integral of interest is that of the function h,
#'which can be multivariate. The estimator uses the variance reduction technique
#'whereby the estimator is the MCMC average between times k and K, with probability
#'going to one as k increases. In the case that K=Inf, the returned value is the 
#'almost sure limit of the bias correction
#'@export
BC <- function(c_chains, h = function(x) x, k = 0, K = 1){
  maxiter <- c_chains$iteration
  if (k > maxiter){
    print("error: k has to be less than the horizon of the coupled chains")
    return(NULL)
  }
  if ((K > maxiter)& is.finite(K)){
    print("error: K has to be less than the horizon of the coupled chains")
    return(NULL)
  }
  # test the dimension of h(X)
  p <- length(h(c_chains$samples1[1,]))
  
  
  if (c_chains$meetingtime <= k + 1){
    deltas_term <- rep(0, p)
  } else {
    
    deltas_term <- rep(0, p)
    for (t in k:min(maxiter-1, c_chains$meetingtime-1)){ # t is as in the report, where the chains start at t=0
      coefficient <- min(t - k + 1, K - k + 1)
      delta_tp1 <- h(c_chains$samples1[t + 1 + 1,]) - h(c_chains$samples2[t+1,]) # the +1's are because R starts indexing at 1
      deltas_term <- deltas_term + coefficient * delta_tp1
    }
    
  }
  
  if(!is.finite(K)){
    return(deltas_term)
  }else{
    return(deltas_term / (K - k + 1))  
  }
  
}


