#### This scripts briefly shows how this package works
# First, some initialization
# remove all objects from R environment
rm(list = ls())
set.seed(17)

# installs the rcpp - useful if R CMD INSTALL requires more priveleges
if (!require(CoupledPIMH)){
  library(devtools)
  devtools::document()
}

library(doRNG)
library(coda)
library(latex2exp)

cores_requested <- min(40,detectCores()-1)
registerDoMC(cores = cores_requested)
setmytheme()
set.seed(17)
module_tree <<- Module("module_tree", PACKAGE = "CoupledPIMH")
TreeClass <<- module_tree$Tree




### ok want to see how using (1/sigma)*logz compares to using variance of estimators
# it's ok to look at one test function
# want to try 3 values of k and m, then vary sigma2


### get experiment settings
source("inst/reproduce/lgssm/lgssm_settings.R")
settings <- lgssm_settings()
dimension <- settings$dimension
theta <- settings$theta
datalength <- settings$datalength
with_as <- settings$with_as
N_arr_large_m <- settings$N_arr_large_m
n_var <- settings$n_var
data_file <- settings$data_file
large_m <- 1e6

# load data
load(data_file)



# let us generate some data from the model
model <-  get_ar(dimension)
precomputed <- model$precompute(theta)

# ok happy with tests now get grid of var z
all_particles <- T
pimh_kernels <- get_pimh_kernels(with_as,all_particles)
pimh_init <- pimh_kernels$pimh_init
single_kernel_pimh <- pimh_kernels$single_kernel
coupled_kernel_pimh <- pimh_kernels$coupled_kernel


var_smoothing_pimh_lst <- list()
var_smoothing_ll_lst <- list()
mt_mat_lst <- list()
large_m_result <- foreach(i = 1:length(N_arr_large_m)) %dorng% {
  
  nparticles <- N_arr_large_m[i]

  print(sprintf('running experiment for nparticles=%.3f',nparticles))
  
  t1 <- Sys.time()
  coupled_smoothing_res <- coupled_pimh_chains(single_kernel_pimh, coupled_kernel_pimh, pimh_init,nparticles,K=large_m,verbose=F)
  
  print('time taken:')
  ttaken <- Sys.time()-t1
  print(ttaken)
  return(list(smoothing_res=coupled_smoothing_res,ttaken=ttaken))
}

timings <- sapply(large_m_result,function(x) x$ttaken)

# get asymptotic variance
burn_indx <- 1e3
asymp_var1 <- rep(NA,length(N_arr_large_m))
asymp_var2 <- rep(NA,length(N_arr_large_m))
asymp_var3 <- rep(NA,length(N_arr_large_m))
asymp_var4 <- rep(NA,length(N_arr_large_m))

asymp_var1_rb <- rep(NA,length(N_arr_large_m))
asymp_var2_rb <- rep(NA,length(N_arr_large_m))
asymp_var3_rb <- rep(NA,length(N_arr_large_m))



for(i in 1:length(N_arr_large_m)){
  print(sprintf('getting asymp var %i of %i',i,length(N_arr_large_m)))
  if(is.matrix(large_m_result[[i]]$smoothing_res$samples1)){
    # state
    chain <- (large_m_result[[i]]$smoothing_res$samples1[burn_indx:(large_m+1),datalength+1,drop=F])
    asymp_var1[i] <- spectrum0.ar(chain)$spec
    chain <- (large_m_result[[i]]$smoothing_res$samples1[burn_indx:(large_m+1),1,drop=F])
    asymp_var2[i] <- spectrum0.ar(chain)$spec
    chain <- rowSums(large_m_result[[i]]$smoothing_res$samples1[burn_indx:(large_m+1),])
    asymp_var3[i] <- spectrum0.ar(chain)$spec
    chain <- rowSums(large_m_result[[i]]$smoothing_res$samples1[burn_indx:(large_m+1),])+rowSums(large_m_result[[i]]$smoothing_res$samples1[burn_indx:(large_m+1),]^2)
    asymp_var4[i] <- spectrum0.ar(chain)$spec  
    
    
    # rb
    chain <- (large_m_result[[i]]$smoothing_res$exp_samples1[burn_indx:(large_m+1),datalength+1,drop=F])
    asymp_var1_rb[i] <- spectrum0.ar(chain)$spec
    chain <- (large_m_result[[i]]$smoothing_res$exp_samples1[burn_indx:(large_m+1),1,drop=F])
    asymp_var2_rb[i] <- spectrum0.ar(chain)$spec
    chain <- rowSums(large_m_result[[i]]$smoothing_res$exp_samples1[burn_indx:(large_m+1),])
    asymp_var3_rb[i] <- spectrum0.ar(chain)$spec

  }
}


# get IACT
IACT1 <- rep(NA,length(N_arr_large_m))
IACT2 <- rep(NA,length(N_arr_large_m))
IACT3 <- rep(NA,length(N_arr_large_m))
IACT4 <- rep(NA,length(N_arr_large_m))


for(i in 1:length(N_arr_large_m)){
  print(sprintf('getting asymp var %i of %i',i,length(N_arr_large_m)))
  if(is.matrix(large_m_result[[i]]$smoothing_res$samples1)){
    # state
    chain <- (large_m_result[[i]]$smoothing_res$samples1[burn_indx:(large_m+1),datalength+1,drop=F])
    IACT1[i] <- (large_m-burn_indx)/effectiveSize(chain)
    chain <- (large_m_result[[i]]$smoothing_res$samples1[burn_indx:(large_m+1),1,drop=F])
    IACT2[i] <- (large_m-burn_indx)/effectiveSize(chain)
    chain <- rowSums(large_m_result[[i]]$smoothing_res$samples1[burn_indx:(large_m+1),])
    IACT3[i] <- (large_m-burn_indx)/effectiveSize(chain)
    chain <- rowSums(large_m_result[[i]]$smoothing_res$samples1[burn_indx:(large_m+1),])+rowSums(large_m_result[[i]]$smoothing_res$samples1[burn_indx:(large_m+1),]^2)
    IACT4[i] <- (large_m-burn_indx)/effectiveSize(chain)
  }
}





# get variance of log likelihood
N_arr_unique <- unique(N_arr_large_m)
var_logz_res <- rep(NA,length(N_arr_unique))
for(i in 1:length(N_arr_unique)){
  
  nparticles <- N_arr_large_m[i]
  print(sprintf('getting varlogz for nparticles=%.3f',nparticles))
  
  # 1. get the variance of the log-likelihood for each value
  print('getting var logz...')
  t1 <- Sys.time()
  
  log_z_runs <- foreach(r = 1:n_var, .combine = rbind) %dorng% {
    pf_res <- CPF(nparticles, model, theta, observations)
    logz1 <- pf_res$logz
    return(logz1)
  }
  s2_est <- as.numeric(var(log_z_runs))
  var_logz_res[i] <-   s2_est
  
  print('time taken for varlogz:')
  print(Sys.time()-t1)
}





# compare actual and predicted meeting times
# for predicted values we need to do a series of integrations over the range of tau appearing in the histograms
# define the conditional geometric probability
p_z <- function(z){return(1-(pnorm(-z/sqrt(s2)) + exp(-z+s2/2)*pnorm((z-s2)/sqrt(s2))))}
# define the integral
I_z <- function(z){P_z <- (1+p_z(z))/(1-p_z(z)) * dnorm(z,s2,sd=sqrt(s2))}
# estimate the integral on a grid of values


asymp_pitt <- rep(NA,length(N_arr_unique))
# ok for each estimate of the variance of the log-likelihood we need to estimate the respective tau probabilities
for(N_i in 1:length(N_arr_unique)){
  s2 <- var_logz_res[N_i]
  
  lower = -10*sqrt(s2)
  upper = 10*sqrt(s2)
  
  asymp_pitt[N_i] <- integrate(I_z, lower, upper)$value
}


save(settings,asymp_var1,asymp_var2,asymp_var3,asymp_var4,
     asymp_var1_rb,asymp_var2_rb,asymp_var3_rb,
     asymp_pitt,var_logz_res,
     IACT1,IACT2,IACT3,IACT4,
     N_arr_large_m,n_var,large_m,
     N_arr_unique,
     file='lgssmlargem.RData')



load('lgssmlargem.RData')


# go from array of number of particles to unique


asymp_pitt_n <- asymp_pitt/asymp_pitt[length(N_arr_large_m)]
asymp_emp <- rbind(asymp_var1_rb,asymp_var2_rb,asymp_var3_rb)[,1:10]
asymp_emp_n <- asymp_emp/asymp_emp[,10]#length(N_arr_large_m)]
matplot(asymp_pitt,t(asymp_emp_n),pch=1:4)
matplot(var_logz_res,t(asymp_emp_n) /var_logz_res,pch=1:4)

### inefficiency plots
get_plot_df <- function(asymp_var,indx){
  asymp_max <- rowMaxs(matrix(asymp_var,length(N_arr_unique)),na.rm=T)/var_logz_res
  asymp_min <- rowMins(matrix(asymp_var,length(N_arr_unique)),na.rm=T)/var_logz_res
  asymp_mean <- rowMeans(matrix(asymp_var,length(N_arr_unique)),na.rm=T)/var_logz_res
  plt_df <- data.frame(t(rbind(asymp_max,asymp_min,asymp_mean)))
  plt_df <- plt_df/(asymp_mean[4])
  plt_df$s2 <- var_logz_res
  plt_df$h = rep(indx,dim(plt_df)[1])
  return(plt_df)
}
plt_df1 <- get_plot_df(asymp_var1,1)
plt_df2 <- get_plot_df(asymp_var2,2)
plt_df3 <- get_plot_df(asymp_var3,3)
plt_df4 <- get_plot_df(asymp_var3,4)
plt_df <- rbind(plt_df1,plt_df2,plt_df3,plt_df4)
plt_df$h <- as.factor(plt_df$h)

pd <- position_dodge(0.2)
library(latex2exp)
g <- ggplot(data=plt_df) + 
  geom_errorbar(aes(x=s2,ymin=asymp_min,ymax=asymp_max,group=h,color=h),position=pd)+
  geom_point(aes(x=s2,y=asymp_mean,group=h,shape=h),position=pd,size=2) +
  xlab(TeX("$\\sigma^2$")) +
  ylab(TeX("IF(h) / $\\sigma^2}$"))

g
name <- 'lgssm_ineffs2_largem.pdf'
ggsave(name,plot=g,width=7,height=5)



### inefficiency plots
get_plot_df <- function(asymp_var,indx){
  asymp_max <- rowMaxs(matrix(asymp_var,length(N_arr_unique)),na.rm=T)
  asymp_min <- rowMins(matrix(asymp_var,length(N_arr_unique)),na.rm=T)
  asymp_mean <- rowMeans(matrix(asymp_var,length(N_arr_unique)),na.rm=T)
  plt_df <- data.frame(t(rbind(asymp_max,asymp_min,asymp_mean)))
  plt_df <- plt_df/(asymp_mean[length(N_arr_unique)])
  plt_df$s2 <- asymp_pitt
  plt_df$h = rep(indx,dim(plt_df)[1])
  return(plt_df)
}
plt_df1 <- get_plot_df(asymp_var1,1)
plt_df2 <- get_plot_df(asymp_var2,2)
plt_df3 <- get_plot_df(asymp_var3,3)
plt_df4 <- get_plot_df(asymp_var3,4)
plt_df <- rbind(plt_df1,plt_df2,plt_df3,plt_df4)
plt_df$h <- as.factor(plt_df$h)

pd <- position_dodge(2)
library(latex2exp)
g <- ggplot(data=plt_df) + 
  geom_errorbar(aes(x=s2,ymin=asymp_min,ymax=asymp_max,group=h,color=h),position=pd)+
  geom_point(aes(x=s2,y=asymp_mean,group=h,shape=h),position=pd,size=2) +
  xlab(TeX("$IF(\\sigma)$ - for $0.2\\leq\\sigma^2\\leq 2$")) +
  ylab(TeX("IF(h)")) +
  xlim(c(0,20))+
  ylim(c(1,8.5))

g
name <- 'lgssm_asympvar_largem.pdf'
ggsave(name,plot=g,width=7,height=5)











### inefficiency plots with IACT
get_plot_df <- function(asymp_var,indx){
  asymp_max <- rowMaxs(matrix(asymp_var,length(N_arr_unique)),na.rm=T)/var_logz_res
  asymp_min <- rowMins(matrix(asymp_var,length(N_arr_unique)),na.rm=T)/var_logz_res
  asymp_mean <- rowMeans(matrix(asymp_var,length(N_arr_unique)),na.rm=T)/var_logz_res
  plt_df <- data.frame(t(rbind(asymp_max,asymp_min,asymp_mean)))
  plt_df <- plt_df
  plt_df$s <- var_logz_res^0.5
  plt_df$h = rep(indx,dim(plt_df)[1])
  return(plt_df)
}
plt_df1 <- get_plot_df(IACT1,1)
plt_df2 <- get_plot_df(IACT2,2)
plt_df3 <- get_plot_df(IACT3,3)
plt_df4 <- get_plot_df(IACT4,4)
plt_df <- rbind(plt_df1,plt_df2,plt_df3,plt_df4)
plt_df$h <- as.factor(plt_df$h)

pd <- position_dodge(0.07)
g <- ggplot(data=plt_df) + 
  geom_errorbar(aes(x=s,ymin=asymp_min,ymax=asymp_max,group=h,color=h),position=pd)+
  geom_point(aes(x=s,y=asymp_mean,group=h,shape=h),position=pd,size=2) +
  xlab(TeX("$\\sigma$")) +
  ylab(TeX("IF(h) / $\\sigma^2}$"))

g
name <- 'lgssm_IACT_s2_largem.pdf'
ggsave(name,plot=g,width=7,height=5)



### inefficiency plots
get_plot_df <- function(asymp_var,indx){
  asymp_max <- rowMaxs(matrix(asymp_var,length(N_arr_unique)),na.rm=T)
  asymp_min <- rowMins(matrix(asymp_var,length(N_arr_unique)),na.rm=T)
  asymp_mean <- rowMeans(matrix(asymp_var,length(N_arr_unique)),na.rm=T)
  plt_df <- data.frame(t(rbind(asymp_max,asymp_min,asymp_mean)))
  #plt_df <- plt_df/(asymp_mean[length(N_arr_unique)])
  plt_df$s2 <- asymp_pitt
  plt_df$h = rep(indx,dim(plt_df)[1])
  return(plt_df)
}
plt_df1 <- get_plot_df(IACT1,1)
plt_df2 <- get_plot_df(IACT2,2)
plt_df3 <- get_plot_df(IACT3,3)
plt_df4 <- get_plot_df(IACT4,4)
plt_df <- rbind(plt_df1,plt_df2,plt_df3,plt_df4)
plt_df$h <- as.factor(plt_df$h)
pd <- position_dodge(2)
g <- ggplot(data=plt_df) + 
  geom_errorbar(aes(x=s2,ymin=asymp_min,ymax=asymp_max,group=h,color=h),position=pd)+
  geom_point(aes(x=s2,y=asymp_mean,group=h,shape=h),position=pd,size=2) +
  xlab(TeX("$IF(\\sigma)$ - for $0.4\\leq\\sigma\\leq 1.7$")) +
  ylab(TeX("IF(h)"))
  #xlim(c(0,20))+
  #ylim(c(1,8.5))

g
name <- 'lgssm_IACT_largem.pdf'
ggsave(name,plot=g,width=7,height=5)






