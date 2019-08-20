#'@rdname get_csmc_refresh_kernels
#'@title CSMC refresh kernels
#'@description Gets the kernels for CSMC refresh
#'@param with_as whether ancestor sampling should be used (TRUE/FALSE)
#'@return Two kernels - the single kernel and the coupled transition kernel 
#'@export
get_csmc_refresh_kernels <- function(with_as){
  
  single_kernel <- function(chain_state,log_pdf_state, iter,nparticles){
    # update first guy
    cpf_res <- CPF(nparticles, model, theta, observations,ref_trajectory=chain_state,with_as=with_as)
    chain_state <- cpf_res$trajectory
    log_pdf_state <- cpf_res$logz
    
    # propose new path
    cpf_res_prop <- CPF(nparticles, model, theta, observations,with_as=with_as)
    xref_prop <- cpf_res_prop$trajectory
    logz_prop <- cpf_res_prop$logz
    
    # accept for either
    u <- runif(1)
    if(log(u)<(logz_prop-log_pdf_state)){
      chain_state <- xref_prop
      log_pdf_state <- logz_prop
    }
    return(list(chain_state=chain_state,log_pdf_state=log_pdf_state))
  }
  
  
  coupled_kernel <- function(chain_state1, chain_state2,log_pdf_state1,log_pdf_state2,iter,nparticles){
    
    
    not_coupled <- !all(chain_state1==chain_state2)
    
    # update CSMC for both
    cpf_res <- CPF(nparticles, model, theta, observations,ref_trajectory=chain_state1,with_as=with_as)
    chain_state1 <- cpf_res$trajectory
    log_pdf_state1 <- cpf_res$logz
    
    if(not_coupled){
      cpf_res <- CPF(nparticles, model, theta, observations,ref_trajectory=chain_state2,with_as=with_as)
      chain_state2 <- cpf_res$trajectory
      log_pdf_state2 <- cpf_res$logz
    }else{
      chain_state2 <- chain_state1
      log_pdf_state2 <- log_pdf_state1
    }
    
    # propose for both guys
    cpf_res_prop <- CPF(nparticles, model, theta, observations,with_as=with_as)
    xref_prop <- cpf_res_prop$trajectory
    logz_prop <- cpf_res_prop$logz
    
    # accept for either
    u <- runif(1)
    if(log(u)<(logz_prop-log_pdf_state1)){
      accept_1 <- T
      chain_state1 <- xref_prop
      log_pdf_state1 <- logz_prop
    }else{
      accept_1 <- F
    }
    
    if(log(u)<(logz_prop-log_pdf_state2)){
      accept_2 <- T
      chain_state2 <- xref_prop
      log_pdf_state2 <- logz_prop
    }else{
      accept_2 <- F
    }
    
    
    return(list(chain_state1=chain_state1,chain_state2=chain_state2,log_pdf_state1=log_pdf_state1,log_pdf_state2=log_pdf_state2))
  }
  return(list(single_kernel=single_kernel,coupled_kernel=coupled_kernel))
}


#'@rdname get_pimh_kernels
#'@title PIMH kernels
#'@description Gets the kernels for PIMH transitions
#'@param with_as whether ancestor sampling should be used (TRUE/FALSE)
#'@param all_particles whether to return the expected path using all the particles
#'@return Two kernels - the single kernel and the coupled transition kernel 
#'@export
get_pimh_kernels <- function(with_as,all_particles=F){
  
  single_kernel <- function(chain_state,log_pdf_state, iter,nparticles){
    # update first guy
    pf_prop <- CPF(nparticles, model, theta, observations,with_as=F)
    xref_prop <- pf_prop$trajectory
    logz_prop <- pf_prop$logz
    
    # accept for either
    u <- runif(1)
    if(log(u)<(logz_prop-log_pdf_state)){
      chain_state <- xref_prop
      log_pdf_state <- logz_prop
    }
    return(list(chain_state=chain_state,log_pdf_state=log_pdf_state))
  }
  
  
  coupled_kernel <- function(chain_state1, chain_state2,log_pdf_state1,log_pdf_state2,iter,nparticles){
    
    # propose for both guys
    pf_prop <- CPF(nparticles, model, theta, observations,with_as=F)
    xref_prop <- pf_prop$trajectory
    logz_prop <- pf_prop$logz
    
    # accept for either
    u <- runif(1)
    if(log(u)<(logz_prop-log_pdf_state1)){
      accept_1 <- T
      chain_state1 <- xref_prop
      log_pdf_state1 <- logz_prop
    }else{
      accept_1 <- F
    }
    
    if(log(u)<(logz_prop-log_pdf_state2)){
      accept_2 <- T
      chain_state2 <- xref_prop
      log_pdf_state2 <- logz_prop
    }else{
      accept_2 <- F
    }
    
    
    return(list(chain_state1=chain_state1,chain_state2=chain_state2,log_pdf_state1=log_pdf_state1,log_pdf_state2=log_pdf_state2))
  }
  
  
  pimh_init <- function(nparticles){
    cpf_res <- CPF(nparticles, model, theta, observations,all_particles=all_particles)
    
    chain_state1 <- cpf_res$trajectory
    logz1 <- cpf_res$logz
    
    cpf_res <- CPF(nparticles, model, theta, observations)
    chain_state2 <- cpf_res$trajectory
    logz2 <- cpf_res$logz
    
    return(list(chain_state1=chain_state1,chain_state2=chain_state2,log_pdf_state1=logz1,log_pdf_state2=logz2))
  }
  
  single_kernel_all_particles <- function(chain_state,log_pdf_state, exp_path,iter,nparticles){
    # update first guy
    pf_prop <- CPF(nparticles, model, theta, observations,with_as=with_as,all_particles=T)
    xref_prop <- pf_prop$trajectory
    logz_prop <- pf_prop$logz
    exp_path_prop <- pf_prop$exp_path
    
    # accept for either
    u <- runif(1)
    if(log(u)<(logz_prop-log_pdf_state)){
      chain_state <- xref_prop
      log_pdf_state <- logz_prop
      exp_path <- exp_path_prop
    }
    return(list(chain_state=chain_state,log_pdf_state=log_pdf_state,exp_path=exp_path))
  }
  
  
  coupled_kernel_all_particles <- function(chain_state1, chain_state2,log_pdf_state1,log_pdf_state2,exp_path1,exp_path2,iter,nparticles){
    
    # propose for both guys
    pf_prop <- CPF(nparticles, model, theta, observations,with_as=with_as,all_particles=T)
    xref_prop <- pf_prop$trajectory
    logz_prop <- pf_prop$logz
    exp_path_prop <- pf_prop$exp_path
    
    # accept for either
    u <- runif(1)
    if(log(u)<(logz_prop-log_pdf_state1)){
      accept_1 <- T
      chain_state1 <- xref_prop
      log_pdf_state1 <- logz_prop
      exp_path1 <- exp_path_prop
    }else{
      accept_1 <- F
    }
    
    if(log(u)<(logz_prop-log_pdf_state2)){
      accept_2 <- T
      chain_state2 <- xref_prop
      log_pdf_state2 <- logz_prop
      exp_path2 <- exp_path_prop
    }else{
      accept_2 <- F
    }
    
    return(list(chain_state1=chain_state1,chain_state2=chain_state2,
                log_pdf_state1=log_pdf_state1,log_pdf_state2=log_pdf_state2,
                exp_path1=exp_path1,
                exp_path2=exp_path2))
  }
  
  
  pimh_init_all_particles <- function(nparticles){
    
    pf_res <- CPF(nparticles, model, theta, observations,all_particles=T)
    exp_path1 <- pf_res$exp_path
    chain_state1 <- pf_res$trajectory
    logz1 <- pf_res$logz
    

    pf_res <- CPF(nparticles, model, theta, observations,all_particles=T)
    exp_path2 <- pf_res$exp_path
    chain_state2 <- pf_res$trajectory
    logz2 <- pf_res$logz
    
    return(list(chain_state1=chain_state1,chain_state2=chain_state2,
                log_pdf_state1=logz1,log_pdf_state2=logz2,
                exp_path1=exp_path1,
                exp_path2=exp_path2))
  }
  
  if(all_particles){
    single_kernel <- single_kernel_all_particles
    coupled_kernel <- coupled_kernel_all_particles
    pimh_init <- pimh_init_all_particles
  }
  
  return(list(single_kernel=single_kernel,coupled_kernel=coupled_kernel,pimh_init=pimh_init))
}



#'@rdname get_cpf_kernels
#'@title CCPF kernels
#'@description Gets the kernels for CCPF transitions
#'@param with_as whether ancestor sampling should be used (TRUE/FALSE)
#'@return Two kernels - the single kernel and the coupled transition kernel 
#'@export
get_cpf_kernels <- function(with_as){
  
  coupled_resampling <- CR_indexmatching
  
  single_kernel <- function(chain_state,log_pdf_state, iter,nparticles){
    # update first guy
    cpf_res <- CPF(nparticles, model, theta, observations,ref_trajectory=chain_state,with_as=with_as)
    chain_state <- cpf_res$trajectory
    log_pdf_state <- cpf_res$logz
    
    return(list(chain_state=chain_state,log_pdf_state=log_pdf_state))
  }
  
  
  coupled_kernel <- function(chain_state1, chain_state2,log_pdf_state1,log_pdf_state2,iter,nparticles){
    
    # update CSMC for both
    cpf_res <- CPF_coupled(nparticles, model, theta, observations,
                           ref_trajectory1=chain_state1,
                           ref_trajectory2=chain_state2,
                           coupled_resampling=coupled_resampling,
                           with_as=with_as)
    chain_state1 <- cpf_res$new_trajectory1
    chain_state2 <- cpf_res$new_trajectory2
    
    return(list(chain_state1=chain_state1,chain_state2=chain_state2,log_pdf_state1=log_pdf_state1,log_pdf_state2=log_pdf_state2))
  }
  
  return(list(single_kernel=single_kernel,coupled_kernel=coupled_kernel))
}





#'@rdname get_pimh_smc_sampler_kernels
#'@title PIMH kernels for SMC sampler
#'@description Gets the kernels for PIMH transitions
#'@param h_func If specified, this will be the Rao-Blackwellised test function
#'@return Two kernels - the single kernel and the coupled transition kernel 
#'@export
get_pimh_smc_sampler_kernels <- function(h_func=NULL){
  
  single_kernel_smc <- function(chain_state,log_pdf_state, iter,nparticles){
    # update first guy
    pf_prop <- smc_sampler(nparticles, model, theta,h_func=h_func)
    xref_prop <- pf_prop$trajectory
    logz_prop <- pf_prop$logz
    
    # accept for either
    u <- runif(1)
    if(log(u)<(logz_prop-log_pdf_state)){
      chain_state <- xref_prop
      log_pdf_state <- logz_prop
    }
    return(list(chain_state=chain_state,log_pdf_state=log_pdf_state))
  }
  
  
  coupled_kernel_smc <- function(chain_state1, chain_state2,log_pdf_state1,log_pdf_state2,iter,nparticles){
    
    # propose for both guys
    pf_prop <- smc_sampler(nparticles, model, theta,h_func=h_func)
    xref_prop <- pf_prop$trajectory
    logz_prop <- pf_prop$logz
    
    # accept for either
    u <- runif(1)
    if(log(u)<(logz_prop-log_pdf_state1)){
      accept_1 <- T
      chain_state1 <- xref_prop
      log_pdf_state1 <- logz_prop
    }else{
      accept_1 <- F
    }
    
    if(log(u)<(logz_prop-log_pdf_state2)){
      accept_2 <- T
      chain_state2 <- xref_prop
      log_pdf_state2 <- logz_prop
    }else{
      accept_2 <- F
    }
    
    
    return(list(chain_state1=chain_state1,chain_state2=chain_state2,log_pdf_state1=log_pdf_state1,log_pdf_state2=log_pdf_state2))
  }
  
  
  pimh_init_smc <- function(nparticles){
    cpf_res <- smc_sampler(nparticles, model, theta,h_func=h_func)
    
    chain_state1 <- cpf_res$trajectory
    logz1 <- cpf_res$logz
    
    cpf_res <- smc_sampler(nparticles, model, theta,h_func=h_func)
    chain_state2 <- cpf_res$trajectory
    logz2 <- cpf_res$logz
    
    return(list(chain_state1=chain_state1,chain_state2=chain_state2,log_pdf_state1=logz1,log_pdf_state2=logz2))
  }
  
  
  return(list(single_kernel_smc=single_kernel_smc,coupled_kernel_smc=coupled_kernel_smc,pimh_init_smc=pimh_init_smc))
}

















