# First, some initialization
# remove all objects from R environment
rm(list = ls())
set.seed(17)
if (!require(CoupledPIMH)){
  library(devtools)
  devtools::document()
}

library(doRNG)
setmytheme()
set.seed(17)


### get experiment settings
source("inst/lgssm/lgssm_settings.R")
settings <- lgssm_settings()
dimension <- settings$dimension
theta <- settings$theta
datalength <- settings$datalength
with_as <- settings$with_as
N_arr <- settings$N_arr
n_var <- settings$n_var
R_bound_check <- settings$R_bound_check
results_file <- settings$results_file
data_file <- settings$data_file
bound_check_save_file <- settings$bound_check_save_file
m_bc_bound_check <- settings$m_bc_bound_check


###
cores_requested <- min(80,detectCores()-1)
registerDoMC(cores = cores_requested)

module_tree <<- Module("module_tree", PACKAGE = "CoupledPIMH")
TreeClass <<- module_tree$Tree




### model & data
model <-  get_ar(dimension)
load(data_file)



# ok want 
# 1. test new coupling
# 2. grid of var logz
# 3. compare meeting times and predictions
# 4. predicted optimisation
# 5. optimisation for large m

# ok happy with tests now get grid of var z

pimh_kernels <- get_pimh_kernels(with_as)
single_kernel_pimh <- pimh_kernels$single_kernel
coupled_kernel_pimh <- pimh_kernels$coupled_kernel

pmmh_init <- function(nparticles){
  cpf_res <- CPF(nparticles, model, theta, observations)
  
  chain_state1 <- cpf_res$trajectory
  logz1 <- cpf_res$logz
  
  cpf_res <- CPF(nparticles, model, theta, observations)
  chain_state2 <- cpf_res$trajectory
  logz2 <- cpf_res$logz
  
  return(list(chain_state1=chain_state1,chain_state2=chain_state2,log_pdf_state1=logz1,log_pdf_state2=logz2))
}



# get meeting times and estimators
# we will want to compare the predicted with the actual meeting times
# and compare the variance of the estimators with the optimal provided by the
# bound
estimators_lst <- list()
run_times <- rep(NA,length(N_arr))
for(i in 1:length(N_arr)){
  nparticles <- N_arr[i]
  print(sprintf('getting meeting times N=%i',nparticles))

  t1 <- Sys.time()
  pimh_iterations <- foreach(r = 1:R_bound_check) %dorng% {
    coupled_smoothing_res <- coupled_pimh_chains(single_kernel_pimh, coupled_kernel_pimh, pmmh_init,nparticles,K=m_bc_bound_check,verbose=F)
  }
  
  print('time taken:')
  t2 <- Sys.time()-t1
  print(t2)
  run_times[i] <- t2
  estimators_lst[[i]] <- pimh_iterations
  
  
  save.image(save_file)
  
}


# save results
save.image(bound_check_save_file)



