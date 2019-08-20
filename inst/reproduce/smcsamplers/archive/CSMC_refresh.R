#### This scripts briefly shows how this package works
# First, some initialization
# remove all objects from R environment
rm(list = ls())
# load package
library(CoupledPIMH)
library(doRNG)

cores_requested <- min(80,detectCores()-1)
registerDoMC(cores = cores_requested)
setmytheme()
set.seed(17)
module_tree <<- Module("module_tree", PACKAGE = "CoupledPIMH")
TreeClass <<- module_tree$Tree

###
# We are going to work with the nonlinear growth model
# first we import it
gss <- get_gss()
model <- gss
save_plots <- T

# let us generate some data from the model
dimension <- 1
theta <- c()
datalength <- 50
x <- matrix(0, nrow = datalength + 1, ncol = 1)
observations <- matrix(0, nrow = datalength, ncol = 1)
x[1,] <- gss$rinit(1, theta, gss$rinit_rand(1, theta))
# equivalent to
#x[1,] <- rnorm(1, 0, sd = sqrt(2))
for (i in 1:datalength){
  x[i+1,] <- gss$rtransition(x[i,], theta, i, gss$rtransition_rand(1, theta))
  # equivalent to 
  # x[i+1,] <- 0.5 * x[i,] + 25 * x[i,] / (1 + x[i,]^2) + 8 * cos(1.2*(i - 1)) + sqrt(10) * rnorm(1)
  observations[i,] <- (x[i+1,]^2) / 20 + rnorm(1)
}
# plot the latent process
plot(0:datalength, x, type = "l")
# plot the observations
plot(1:datalength, observations, type = "l")
#

# now let us estimate the smoothing mean with the proposed method
# choose a number of particles
nparticles <- 2^10
# choose a coupled resampling scheme 
coupled_resampling <- CR_indexmatching
# ancestor sampling or not
with_as <- TRUE
# truncation variable parameter between 0 and 1 (0 means that we don't use any truncation variable)
lambda <- 0
# put all the parameters in a list
algorithmic_parameters <- list(nparticles = nparticles, coupled_resampling = coupled_resampling, lambda = lambda, with_as = with_as)
# how many estimators do we want?
R <- 100


# run the particle filter to obtain a sample trajector
xref <- CPF(nparticles, gss, theta, observations)$trajectory
plot(xref[1,],type='l')

M <- 200
xrefs <- matrix(0,M,datalength+1)
varz <- rep(NA,M)
for(i in 1:M){
  print(i)
  xref <- CPF(nparticles, gss, theta, observations,with_as=T)
  xrefs[i,] <- xref$trajectory[1,]
  varz[i] <- xref$logz
}

matplot(t(xrefs),type='l',col=rgb(1,0,0,0.01 ))
matplot(x[,1],add=T,pch=1)



single_kernel <- function(chain_state,log_pdf_state, iter,nparticles,with_as=T){
  # update first guy
  cpf_res <- CPF(nparticles, model, theta, observations,ref_trajectory=chain_state,with_as=with_as)
  chain_state <- cpf_res$trajectory
  log_pdf_state <- cpf_res$logz
  
  # propose for both guys
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




coupled_kernel <- function(chain_state1, chain_state2,log_pdf_state1,log_pdf_state2,iter,nparticles,with_as=T){
  
  
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



single_kernel_cpf <- function(chain_state,log_pdf_state, iter,nparticles,with_as=T){
  # update first guy
  cpf_res <- CPF(nparticles, model, theta, observations,ref_trajectory=chain_state,with_as=with_as)
  chain_state <- cpf_res$trajectory
  log_pdf_state <- cpf_res$logz
  
  return(list(chain_state=chain_state,log_pdf_state=log_pdf_state))
}


coupled_kernel_cpf <- function(chain_state1, chain_state2,log_pdf_state1,log_pdf_state2,iter,nparticles,with_as=T){
  
  # update CSMC for both
  cpf_res <- CPF_coupled(nparticles, model, theta, observations,
                         ref_trajectory1=chain_state1,
                         ref_trajectory2=chain_state2,
                         coupled_resampling=algorithmic_parameters$coupled_resampling,
                         with_as=with_as)
  chain_state1 <- cpf_res$new_trajectory1
  chain_state2 <- cpf_res$new_trajectory2
  
  return(list(chain_state1=chain_state1,chain_state2=chain_state2,log_pdf_state1=log_pdf_state1,log_pdf_state2=log_pdf_state2))
}


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
n_var <- 1000
log_z_runs <- foreach(r = 1:n_var, .combine = rbind) %dorng% {
  pf_res <- CPF(nparticles, model, theta, observations)
  logz1 <- pf_res$logz
  return(logz1)
}
var(log_z_runs)

R <- 100
csmc_iterations <- foreach(r = 1:R, .combine = rbind) %dorng% {
  coupled_smoothing_res <- coupled_pmmh_chains(single_kernel, coupled_kernel, pmmh_init,nparticles,K=1)
  coupled_smoothing_res$meetingtime
}

hist(csmc_iterations,breaks=50)



# simplest estimator 
k <- 1
K <- 1
R <- 2000

h_bar_estimators_batch <- foreach(r = 1:R) %dorng% {
  ptm <- proc.time()
  coupled_smoothing_res <- coupled_pmmh_chains(single_kernel, coupled_kernel, pmmh_init,nparticles,K=K)
  h_bar_res <- H_bar(coupled_smoothing_res,k=k,K=K)
  iteration_time <- proc.time() - ptm
  return(list(h_bar=h_bar_res,meetingtime=coupled_smoothing_res$meetingtime,iteration_time=iteration_time))
}

# coupled cpf values
h_bar_estimators_coupled_cpf_batch <- foreach(r = 1:R) %dorng% {
  ptm <- proc.time()
  coupled_smoothing_res <- coupled_pmmh_chains(single_kernel_cpf, coupled_kernel_cpf, pmmh_init,nparticles,K=K)
  h_bar_res <- H_bar(coupled_smoothing_res,k=k,K=K)
  iteration_time <- proc.time() - ptm
  return(list(h_bar=h_bar_res,meetingtime=coupled_smoothing_res$meetingtime,iteration_time=iteration_time))
}

# get summary statistics for estimators
h_bar_estimators <- t(sapply(h_bar_estimators_batch,function(x) x$h_bar))
meeting_times <- sapply(h_bar_estimators_batch,function(x) x$meetingtime)

h_bar_mean <- apply(h_bar_estimators,2,mean)
h_bar_std <- apply(h_bar_estimators,2,var)^0.5

h_bar_estimators_coupled_cpf <- t(sapply(h_bar_estimators_coupled_cpf_batch,function(x) x$h_bar))
meeting_times_coupled_cpf <- sapply(h_bar_estimators_coupled_cpf_batch,function(x) x$meetingtime)

h_bar_mean_coupled_cpf <- apply(h_bar_estimators_coupled_cpf,2,mean)
h_bar_std_coupled_cpf <- apply(h_bar_estimators_coupled_cpf,2,var)^0.5

# plot summaries for estimators
matplot(t(h_bar_estimators),col=rgb(1,0,0,0.2),type='l',main='csmc')
matplot(x,col='blue',add=T,type='l')

matplot(t(h_bar_estimators_coupled_cpf),col=rgb(1,0,0,0.2),type='l',main='cpf')
matplot(x,col='blue',add=T,type='l')

# get mean times
csmc_time <- mean(sapply(h_bar_estimators_batch,function(x) x$iteration_time[3]))
coupled_cpf_time <- mean(sapply(h_bar_estimators_coupled_cpf_batch,function(x) x$iteration_time[3]))

# plot time-weigted variances
results_matrix <- rbind(h_bar_std_coupled_cpf*coupled_cpf_time,h_bar_std*csmc_time)
results_matrix_unweighted <- rbind(h_bar_std_coupled_cpf,h_bar_std)

name <- 'inst/exps/gss/low_K_est_var.pdf'
pdf(name)
barplot(results_matrix,beside=T,names.arg=1:(datalength+1),legend=c('Coupling PF','CSMC refresh'),
        args.legend = list(x='topleft'),
        xlab='t',
        ylab='Time weighted standard dev.')
dev.off()
print(sprintf('written to %s',name))

name <- 'inst/exps/gss/low_K_est_var_unweighted.pdf'
pdf(name)
barplot(results_matrix_unweighted,beside=T,names.arg=1:(datalength+1),legend=c('Coupling PF','CSMC refresh'),
        args.legend = list(x='topleft'),
        xlab='t',
        ylab='Estimator standard dev.')
dev.off()
print(sprintf('written to %s',name))



hist_breaks <- seq(from=0,to=max(c(meeting_times,meeting_times_coupled_cpf)),length.out=200)
hist(meeting_times-1+0.2,breaks=hist_breaks,col=rgb(1,0,1,0.5),main='',xlab = 'Meeting times')
hist(meeting_times_coupled_cpf-1,add=T,breaks=hist_breaks,col=rgb(1,0,0,0.5))
legend('topright',legend=c('Meeting times CSMC','Meeting times CPF'),fill=c(rgb(1,0,1,0.5),rgb(1,0,0,0.5)))



# estimator with k,K>1
k_csmc <- floor(as.numeric(quantile(meeting_times,0.95)))
K_csmc <- 10*k_csmc
h_bar_estimators_batch_large_K <- foreach(r = 1:R) %dorng% {
  ptm <- proc.time()
  coupled_smoothing_res <- coupled_pmmh_chains(single_kernel, coupled_kernel, pmmh_init,nparticles,K=K_csmc)
  h_bar_res <- H_bar(coupled_smoothing_res,k=k_csmc,K=K_csmc)
  iteration_time <- proc.time() - ptm
  return(list(h_bar=h_bar_res,meetingtime=coupled_smoothing_res$meetingtime,iteration_time=iteration_time))
}

# coupled cpf values
k_cpf <- floor(as.numeric(quantile(h_bar_estimators_coupled_cpf,0.95)))
K_cpf <- 10*k_cpf
ptm_coupling <- proc.time()
h_bar_estimators_coupled_cpf_batch_large_K <- foreach(r = 1:R) %dorng% {
  ptm <- proc.time()
  coupled_smoothing_res <- coupled_pmmh_chains(single_kernel_cpf, coupled_kernel_cpf, pmmh_init,nparticles,K=K_cpf)
  h_bar_res <- H_bar(coupled_smoothing_res,k=k_cpf,K=K_cpf)
  iteration_time <- proc.time() - ptm
  return(list(h_bar=h_bar_res,meetingtime=coupled_smoothing_res$meetingtime,iteration_time=iteration_time))
}

# get estimators for large K
h_bar_estimators_large_K <- t(sapply(h_bar_estimators_batch_large_K,function(x) x$h_bar))
meeting_times_large_K <- sapply(h_bar_estimators_batch_large_K,function(x) x$meetingtime)

h_bar_mean_large_K <- apply(h_bar_estimators_large_K,2,mean)
h_bar_std_large_K <- apply(h_bar_estimators_large_K,2,var)^0.5

h_bar_estimators_coupled_cpf_large_K <- t(sapply(h_bar_estimators_coupled_cpf_batch_large_K,function(x) x$h_bar))
meeting_times_coupled_cpf_large_K <- sapply(h_bar_estimators_coupled_cpf_batch_large_K,function(x) x$meetingtime)

h_bar_mean_coupled_cpf_large_K <- apply(h_bar_estimators_coupled_cpf_large_K,2,mean)
h_bar_std_coupled_cpf_large_K <- apply(h_bar_estimators_coupled_cpf_large_K,2,var)^0.5

# get mean times
csmc_time_large_K <- mean(sapply(h_bar_estimators_batch_large_K,function(x) x$iteration_time[3]))
coupled_cpf_time_large_K <- mean(sapply(h_bar_estimators_coupled_cpf_batch_large_K,function(x) x$iteration_time[3]))


# plot estimators
matplot(t(h_bar_estimators_large_K[sample(1:R,100),]),col=rgb(1,0,0,0.2),type='l',main='csmc')
matplot(x,col='blue',add=T,type='l')

matplot(t(h_bar_estimators_coupled_cpf_large_K[sample(1:R,100),]),col=rgb(1,0,0,0.2),type='l',main='cpf')
matplot(x,col='blue',add=T,type='l')

# show time weighted performances
name <- 'inst/exps/gss/large_K_est_var.pdf'
pdf(name)
results_matrix <- rbind(h_bar_std_coupled_cpf_large_K*coupled_cpf_time_large_K,h_bar_std_large_K*csmc_time_large_K)
barplot(results_matrix,beside=T,legend=c('Coupling PF','CSMC refresh'),
        args.legend = list(x='topleft'),
        names.arg=1:(datalength+1),
        xlab='t',
        ylab='Time weighted standard dev')
dev.off()
print(sprintf('written to %s',name))


name <- 'inst/exps/gss/large_K_est_var_unweighted.pdf'
pdf(name)
results_matrix <- rbind(h_bar_std_coupled_cpf_large_K,h_bar_std_large_K)
barplot(results_matrix,beside=T,legend=c('Coupling PF','CSMC refresh'),
        args.legend = list(x='topleft'),
        names.arg=1:(datalength+1),
        xlab='t',
        ylab='Estimator standard dev')
dev.off()
print(sprintf('written to %s',name))



name <- 'inst/exps/gss/meeting_times.pdf'
pdf(name)
hist_breaks <- seq(from=0,to=max(c(meeting_times_large_K,meeting_times_coupled_cpf_large_K)),length.out=200)
hist(meeting_times_large_K-1+0.2,breaks=hist_breaks,col=rgb(1,0,1,0.5),main='',xlab = 'Meeting times')
hist(meeting_times_coupled_cpf_large_K-1,add=T,breaks=hist_breaks,col=rgb(1,0,0,0.5))
legend('topright',legend=c('Meeting times CSMC','Meeting times CPF'),fill=c(rgb(1,0,1,0.5),rgb(1,0,0,0.5)))
dev.off()
print(sprintf('written to %s',name))



# # coupling cpf pierre's implementatione
# ptm_coupling_pi <- proc.time()
# estimates <- foreach(r = 1:R, .combine = rbind) %dorng% {
#   res <- rheeglynn_estimator(observations, gss, theta, algorithmic_parameters)
#   data.frame(irep = r, time = 0:datalength, estimate = t(res$estimate), iteration = res$iteration)
# }
# ptm_coupling_pi <- proc.time() - ptm_coupling
# #
# # how many iterations did the estimators take?
# summary(estimates$iteration)
# # now plot the confidence intervals on the smoothing mean
# estimates.df <- estimates %>% select(irep, time, starts_with("estimate")) %>% group_by(time) %>%
#   summarise(m = mean(estimate), s = sd(estimate))
# 
# barplot(rbind(h_bar_std_coupled_cpf,h_bar_std,estimates.df$s),beside=T,legend=c('Coupling PF','CSMC refresh' ,'Coupling CPF (Pierre\'s implementation'),args.legend = list(x='topleft'))
# 
# 
# diff_methods <- estimates.df$s-(1.5*h_bar_std/sqrt(R))
# plot(diff_methods)
# 
# hist(diff_methods,breaks=20)
# 


