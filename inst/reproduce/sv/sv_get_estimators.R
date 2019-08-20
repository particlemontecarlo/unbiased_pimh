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
source("inst/reproduce/sv/sv_settings.R")
settings <- sv_settings()
with_as <- settings$with_as
N_arr <- settings$N_arr
n_var <- settings$n_var
R <- settings$R
R_hkm <- settings$R_hkm
K_max <- settings$K_max
varsave_file <- settings$varsave_file
meeting_times_varlogz_save_file <- settings$meeting_times_varlogz_save_file
sv_data_file <- settings$sv_data_file
estimators_res_file <- settings$estimators_res_file
parameter_mle_file <- settings$parameter_mle_file
hkm_save_file <- settings$hkm_save_file
hkm_small_save_file <- settings$hkm_small_save_file


###
cores_requested <- min(30,detectCores()-1)
registerDoMC(cores = cores_requested)

module_tree <<- Module("module_tree", PACKAGE = "CoupledPIMH")
TreeClass <<- module_tree$Tree



# load parameters
load(parameter_mle_file)
theta <- theta_mle

# load data
load(file=sv_data_file)
observations_plt <- observationsDF

# # plot data
# g <- ggplot(data = observations_plt, aes(x = index, y = y)) 
# g <- g + geom_line() + ylab("observations")
# print(g)
observations <- as.matrix(observations_plt['y'])

# get model
model <- get_levy_sv()


# estimate var log z
t1 <- Sys.time()
varlogz_arr <- rep(NA,length(N_arr))
for(i in 1:length(N_arr)){
  nparticles <- N_arr[i]
  print(sprintf('getting var logz N=%i',nparticles))
  log_z_runs <- foreach(r = 1:n_var, .combine = rbind) %dorng% {
    pf_res <- CPF(nparticles, model, theta, observations)
    logz1 <- pf_res$logz
    return(logz1)
  }
  varlogz_arr[i] <- var(log_z_runs)[1]
}

print('N arr :')
print(N_arr)
print('var log z')
print(varlogz_arr)
print('Time taken:')
var_est_timing <- Sys.time()-t1
print(var_est_timing)


save(varlogz_arr,file=varsave_file)








### get estimators
pimh_kernels <- get_pimh_kernels(with_as)
single_kernel_pimh <- pimh_kernels$single_kernel
coupled_kernel_pimh <- pimh_kernels$coupled_kernel
pimh_init <- pimh_kernels$pimh_init


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
  pimh_iterations <- foreach(r = 1:R_hkm) %dorng% {
    coupled_smoothing_res <- coupled_pimh_chains(single_kernel_pimh, coupled_kernel_pimh, pimh_init,nparticles,K=K_max,verbose=F)
  }
  
  print('time taken:')
  t2 <- Sys.time()-t1
  print(t2)
  run_times[i] <- t2
  estimators_lst[[i]] <- pimh_iterations
  
  
  save.image(hkm_save_file)
  
}


# extract meeting times
mt_arr <- matrix(NA,length(N_arr),R_hkm)
for(i in 1:length(N_arr)){
  mt_arr[i,] <- sapply(estimators_lst[[i]],function(x) x$meetingtime)
}

# save results
save.image(hkm_save_file)











### Process the results for the plots
# the results will reduce the size
k <- 20


N_select <- c(1:5)
N_arr <- N_arr[N_select]
varlogz_arr <- varlogz_arr[N_select]

# hkm of state
est_var <- rep(NA,length(N_arr))
est_hsquared_mean <- rep(NA,length(N_arr))
est_hmean <- rep(NA,length(N_arr))
for(i in 1:length(N_arr)){
  nparticles <- N_arr[i]
  print(sprintf('getting estimators N=%i',nparticles))
  
  estimators <- sapply(estimators_lst[[i]],function(x) H_bar(x,k=k,K=K_max,h=function(y) sum(y)))
  est_var[i] <- var(estimators)
  est_hmean[i] <- mean(estimators)
  estimators_hsquared <- sapply(estimators_lst[[i]],function(x) H_bar(x,k=k,K=K_max,h=function(y) sum(y)^2))
  est_hsquared_mean[i] <- mean(estimators_hsquared)
}
var_pi_h <- mean(est_hsquared_mean-est_hmean^2)

# hkm of likelihood
est_var_ll <- rep(NA,length(N_arr))
for(i in 1:length(N_arr)){
  nparticles <- N_arr[i]
  print(sprintf('getting estimators N=%i',nparticles))
  
  estimator_batch <- estimators_lst[[i]]
  for(r in 1:length(estimator_batch)){
    estimator_batch[[r]]$samples1 <- estimator_batch[[r]]$log_pdf1
    estimator_batch[[r]]$samples2 <- estimator_batch[[r]]$log_pdf2
  }
  
  estimators <- sapply(estimator_batch,function(x) H_bar(x,k=k,K=K_max,h=function(y) sum(y)))
  est_var_ll[i] <- var(estimators)
  
}


# produce the condensed file of meeting times
print('saving condensed file of meeting times...')
save(N_arr,mt_arr,varlogz_arr,est_var,est_var_ll,var_pi_h,smoothed_ests,file=hkm_small_save_file)

















