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
R <- settings$R
save_file <- settings$save_file
results_file <- settings$results_file
data_file <- settings$data_file



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
  pimh_iterations <- foreach(r = 1:R) %dorng% {
    coupled_smoothing_res <- coupled_pimh_chains(single_kernel_pimh, coupled_kernel_pimh, pmmh_init,nparticles,K=1,verbose=F)
  }
  
  print('time taken:')
  t2 <- Sys.time()-t1
  print(t2)
  run_times[i] <- t2
  estimators_lst[[i]] <- pimh_iterations
  
  
  save.image(save_file)
  
}


# extract meeting times
mt_arr <- matrix(NA,length(N_arr),R)
for(i in 1:length(N_arr)){
  mt_arr[i,] <- sapply(estimators_lst[[i]],function(x) x$meetingtime)
}



save.image(save_file)




### plot optimal variance vs bound
# we need to weight the variance by the number of particles as well as the time taken to obtain the estimators

# 
# 
# # Finally compare these with an increasing k and m
# # for each value of m we will iterate over the grid of N
# # we will want
# # 1) the acceptance rate
# # 2) the variance of the estimators 
# # 3) the mean and mean squared meeting time
# 
# m_arr <- c(1,seq(25,200,25))
# phi_k <- 0.1 # this is used to get k=phi_k*m
# n_m_arr <- length(m_arr)
# 
# ar_mat1 <- matrix(NA,n_m_arr,length(N_arr))
# var_mat <- matrix(NA,n_m_arr,length(N_arr))
# mean_tau_mat <- matrix(NA,n_m_arr,length(N_arr))
# mean_tau2_mat <- matrix(NA,n_m_arr,length(N_arr))
# 
# for(i_m in 1:length(m_arr)){
#   
#   K <- m_arr[i_m]
#   k <- floor(K*phi_k)
#   
#   print(sprintf('getting estimators on grid of m, m=%i',K))
#   
#   # ok, iterate over N
#   for(i in 1:length(N_arr)){
#     nparticles <- N_arr[i]
#     print(sprintf('getting estimators N=%i',nparticles))
#     
#     # get the unbiased estimators
#     t1 <- Sys.time()
#     pimh_iterations <- foreach(r = 1:R) %dorng% {
#       coupled_smoothing_res <- coupled_pimh_chains(single_kernel_pimh, coupled_kernel_pimh, pmmh_init,nparticles,K=K,verbose=F)
#     }
#     
#     # get the acceptance rates, variance and moments of the meeting times
#     # start with the meeting times
#     mt <- sapply(pimh_iterations,function(x) x$meetingtime)
#     mean_tau_mat[i_m,i] <- mean(mt)
#     mean_tau2_mat[i_m,i] <- mean(mt^2)
#     
#     # next get the variance of the estimators
#     estimators <- sapply(pimh_iterations,function(x) H_bar(x,k=k,K=K,h=h1))
#     var_mat[i_m,i] <- var(estimators)
#     
#     # next get the mean acceptance rate of the samples
#     ar_mat1[i_m,i] <- mean(sapply(pimh_iterations,function(x) acc_rate(x$samples1)))
#   }
#   
# }
# 
# # save a snapshot
# save.image(save_file)
# 
# 
# # ok! running on the servers!! can process results while it's running
# 
# # having done all this we should plot the results
# load(results_file)
# 
# ### plot histogram of meeting times and compare with their predicted values
# 
# 
# # for predicted values we need to do a series of integrations over the range of tau appearing in the histograms
# # define the conditional geometric probability
# p_z <- function(z){return(pnorm(-z/sqrt(s2)) + exp(-z+s2/2)*pnorm((z-s2)/sqrt(s2)))}
# # define the integral
# I_z <- function(z){P_z <- p_z(z)*(1-p_z(z))^(k-1) * dnorm(z,0,sd=sqrt(s2))}
# # estimate the integral on a grid of values 
# k_max <- max(mt_arr)
# k_arr = 1:k_max
# 
# p_k_mat <- matrix(NA,length(N_arr),k_max)
# # ok for each estimate of the variance of the log-likelihood we need to estimate the respective tau probabilities
# for(N_i in 1:length(N_arr)){
#   nparticles <- N_arr[N_i]
#   s2 <- varlogz_arr[N_i]
#   
#   lower = -10*sqrt(s2)
#   upper = 10*sqrt(s2)
# 
#   for(i in 1:length(k_arr)){
#     k <- k_arr[i]
#     p_k_mat[N_i,i] <- integrate(I_z, lower, upper)$value
#     
#   }
# }
# 
# 
# # plot the histogram of empirical meeting times
# mt_max <- max(mt_arr)
# count_arr <- matrix(NA,length(N_arr),mt_max)
# for(n_i in 1:length(N_arr)){
#   for(i in 1:mt_max){
#     count_arr[n_i,i] <- sum(mt_arr[n_i,]==i)
#   }
# }
# 
# prob_arr <- count_arr/R
# sd_err <- sqrt(prob_arr*(1-prob_arr)/R)
# prob_arr_u <- prob_arr + 2*sd_err
# prob_arr_l <- prob_arr - 2*sd_err
# 
# plot_df <- data.frame(t(prob_arr))
# names(plot_df) <- N_arr
# plot_df$tau_val <- 1:mt_max
# plot_df_m <- reshape2::melt(plot_df,id='tau_val')
# names(plot_df_m)[2] <- 'N'
# plot_df_m$N <- as.factor(plot_df_m$N)
# 
# plot_df_u <- data.frame(t(prob_arr_u))
# names(plot_df_u) <- N_arr
# plot_df_u$tau_val <- 1:mt_max
# plot_df_m_u <- reshape2::melt(plot_df_u,id='tau_val')
# names(plot_df_m_u)[2] <- 'N'
# plot_df_m_u$N <- as.factor(plot_df_m_u$N)
# names(plot_df_m_u)[3] <- 'value_u'
# 
# plot_df_l <- data.frame(t(prob_arr_l))
# names(plot_df_l) <- N_arr
# plot_df_l$tau_val <- 1:mt_max
# plot_df_l <- reshape2::melt(plot_df_l,id='tau_val')
# names(plot_df_l)[2] <- 'N'
# plot_df_l$N <- as.factor(plot_df_l$N)
# names(plot_df_l)[3] <- 'value_l'
# 
# plot_df <- cbind(plot_df_m,plot_df_m_u$value_u,plot_df_l$value_l)
# names(plot_df)[4:5] <- c('u','l')
# 
# plot_df_exact <- data.frame(t(p_k_mat))
# names(plot_df_exact) <- N_arr
# plot_df_exact$tau_val <- 1:mt_max
# 
# plot_df_exact_m <- reshape2::melt(plot_df_exact,id='tau_val')
# names(plot_df_exact_m)[2] <- 'N'
# plot_df_exact_m$N <- as.factor(plot_df_exact_m$N)
# 
# ggplot() +
#   geom_bar(data = plot_df_exact_m,aes(x=tau_val,y=value,fill=N),stat='identity',position='dodge')+
#   geom_errorbar(data=plot_df, aes(x=tau_val,y=value,fill=N,ymin=l,ymax=u),position='dodge') + 
#   xlim(0,8)
# 
# 
# ### plot comparison of variance for k=m=1
# var_ests <-  diag(var(t(h2_arr)))
# #rowMeans(mt_arr^2+mt_arr)#
# var_proxy_l <- 1*rowMeans((mt_arr^2+mt_arr))# 
# 
# var_proxy_bound <- rep(NA,length(N_arr))
# for(i in 1:length(N_arr)){
#   mt_i <- mt_arr[i,]
#   var_proxy_bound[i] <- mean((2+2*(mt_i[mt_i>1]-1))^2)
# }
# 
# var_proxy <-rowMeans((2*mt_arr-1)^2)#-rowMeans(mt_arr))#(2*rowMeans(mt_arr)-1)+4*(diag(var(t(mt_arr))))#+2*sqrt(diag(var(t(mt_arr)))))
# comp_time <- rowMeans(mt_arr)
# 
# ylim_arr <- c(0.1*min(c(var_ests,var_proxy)),1.5*max(c(var_ests,var_proxy)))
# plot(varlogz_arr,var_ests,ylim=ylim_arr)
# matplot(varlogz_arr,10*var_proxy_bound,add=T,pch=2)
# matplot(varlogz_arr,0.5*var_proxy,add=T,pch=2)
# matplot(varlogz_arr,var_proxy_l,add=T,pch=3)
# 
# 
# ineff <- var_ests*N_arr
# ineff_proxy <-1*var_proxy*N_arr
# #ineff_proxy_l <- 600*var_proxy_l*comp_time/(varlogz_arr)
# 
# ylim_arr <- c(min(c(ineff,ineff_proxy)),max(c(ineff,ineff_proxy)))
# plot(N_arr,ineff,ylim=ylim_arr)
# #matplot(N_arr,ineff_proxy,add=T)
# matplot(N_arr,ineff_proxy,add=T)
# #matplot(varlogz_arr,ineff_proxy_l,add=T)
# 
# ### plot optimal values of N as m increases
# 
# 
# s2 <- 0.1
# r_Z <- function(){
#   return(rnorm(1,mean=-s2/2,sd=sqrt(s2)))
# }
# 
# 
# 
# M <- 10000
# counters <- rep(NA,M)
# for(m in 1:M){
#   
#   Z_u <- max(r_Z(),r_Z())
#   not_met <- T
#   counter <- 0
#   
#   
#   while(not_met){
#     log_W <- r_Z()
#     
#     if(log(runif(1))<(log_W-Z_u)){
#       not_met <- F
#     }
#     counter <- counter+1
#   }
#   counters[m] <- counter
# }
# hist(counters,freq=F,breaks=1000)
# h <- hist(counters, breaks = 1000, plot=FALSE)
# h$counts=h$counts/sum(h$counts)
# plot(h)
# 
# n_Z <- 10000000
# n_s2 <- 20
# 
# s2_arr <- seq(0.05,0.15,length.out=n_s2)
# var_arr <- rep(NA,n_s2)
# for(i in 1:n_s2){
#   s2 <- s2_arr[i]
#   Z <- rnorm(n_Z)*sqrt(s2)
#   pp <- p_z(Z)
#   tau <- rgeom(n_Z,pp)+1
#   var_arr[i] <- mean(tau^2*(tau+1)^2)
#   
# }
# 
# plot(s2_arr,(var_arr)/s2_arr)
# 
# 
# n_Z <- 100000
# n_s2 <- 20
# 
# s2_arr <- seq(0.2,1,length.out=n_s2)
# var_arr <- rep(NA,n_s2)
# for(i in 1:n_s2){
#   s2 <- s2_arr[i]
#   Z <- rnorm(n_Z)*sqrt(s2)
#   pp <- p_z(Z)
#   tau <- rgeom(n_Z,pp)+1
#   #tau_select <- tau[tau-1>=1]
#   var_arr[i] <- mean((tau+1)^2*(tau+2)^2)*mean(tau>1)
#   
# }
# 
# plot(s2_arr,sqrt(var_arr)/s2_arr)
# 
# 
# 
# 
# 
# 
# n_Z <- 100000
# n_s2 <- 20
# 
# s2_arr <- seq(0.01,0.5,length.out=n_s2)
# var_arr <- rep(NA,n_s2)
# for(i in 1:n_s2){
#   s2 <- s2_arr[i]
#   Z <- rnorm(n_Z)*sqrt(s2)
#   pp <- p_z(Z)
#   tau <- rgeom(n_Z,pp)+1
#   #tau_select <- tau[tau-1>=1]
#   var_arr[i] <- mean((tau+1)^2*(tau+2)^2)
#   
# }
# 
# 
# 
# 
# p_k_mat <- matrix(NA,length(N_arr),k_max)
# # ok for each estimate of the variance of the log-likelihood we need to estimate the respective tau probabilities
# for(N_i in 1:length(N_arr)){
#   nparticles <- N_arr[N_i]
#   s2 <- varlogz_arr[N_i]
#   
#   lower = -10*sqrt(s2)
#   upper = 10*sqrt(s2)
#   
#   for(i in 1:length(k_arr)){
#     k <- k_arr[i]
#     p_k_mat[N_i,i] <- integrate(I_z, lower, upper)$value
#     
#   }
# }
# 
# 
# plot(s2_arr,sqrt(var_arr)/s2_arr)
# 
# 
# 



