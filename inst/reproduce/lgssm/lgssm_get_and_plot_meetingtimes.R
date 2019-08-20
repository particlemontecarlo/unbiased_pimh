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


### this script takes around 30 mins on a 90 core server



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
cores_requested <- min(30,detectCores()-1)
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


### plot meeting times at different quantiles
# ok, we want to look at distribution of meeting times as a function of sigma2  to estimate how to pick m to look at bias correction
print('plotting results...')
pmmh_mts <- mt_arr
mt_arr <- mt_arr
pm_mt_quantiles1 <- apply(pmmh_mts,1,function(x) quantile(x,0.99))
pm_mt_quantiles2 <- apply(pmmh_mts,1,function(x) quantile(x,0.999))
pm_mt_quantiles3 <- apply(pmmh_mts,1,function(x) quantile(x,0.9999))
pm_mt_max <- apply(pmmh_mts,1,function(x) max(x))

name <- paste('lgssm_mt.pdf',sep="")
ggplot_df = data.frame(N=N_arr,q0=pm_mt_quantiles1,q1=pm_mt_quantiles2,q2=pm_mt_quantiles3,q4=pm_mt_max)

ggplot_melt <- reshape2::melt((ggplot_df),id='N')
ggplot_melt$quantile <- as.factor(ggplot_melt$variable)
levels(ggplot_melt$quantile) <- c(0.99,0.999,0.9999,'max')

breaks <- 10^(-10:10)
minor_breaks <- rep(1:9, 21)*(10^rep(-10:10, each=9))

g <- ggplot(data=ggplot_melt) +
  geom_point(aes(x=N,y=value,shape=quantile),cex=2) +
  geom_line(aes(x=N,y=value,linetype=quantile)) +
  scale_y_log10(breaks =breaks ,minor_breaks=minor_breaks) +
  ylab('meeting time')
g
#  scale_y_log10()
ggsave(name,plot=g,width=7,height=5)
# get hkm estimators




### plot meeting times at different quantiles
# ok, we want to look at distribution of meeting times as a function of sigma2  to estimate how to pick m to look at bias correction
N_plt <- N_arr
varlogz_arr <- varlogz_arr
v_plt <- varlogz_arr
name <- 'varlogz_lgssm.pdf'
ggplot_df = data.frame(N=N_plt,q0=v_plt)

ggplot_melt <- reshape2::melt((ggplot_df),id='N')
ggplot_melt$quantile <- as.factor(ggplot_melt$variable)
levels(ggplot_melt$quantile) <- c('var log-likelihood')

# breaks <- 10^(-10:10)
# minor_breaks <- rep(1:9, 21)*(10^rep(-10:10, each=9))
minor_breaks <- c(seq(0.1,0.5,0.1),seq(1,5,1))
g <- ggplot(data=ggplot_melt) +
  geom_point(aes(x=N,y=value,shape=quantile),cex=2) +
  geom_line(aes(x=N,y=value,linetype=quantile)) +
  # scale_y_continuous(breaks =breaks ,minor_breaks=minor_breaks) +
  # scale_y_log10(breaks =minor_breaks)+# ,minor_breaks=minor_breaks) +
  ylab('var log-likelihood') +
  theme(legend.title = element_blank()) +
  guides(shape=FALSE,linetype=FALSE)
g


ggsave(name,plot=g,width=7,height=5)








# plot comparison with predicted values


# compare actual and predicted meeting times
# for predicted values we need to do a series of integrations over the range of tau appearing in the histograms
# define the conditional geometric probability
p_z <- function(z){return(pnorm(-z/sqrt(s2)) + exp(-z+s2/2)*pnorm((z-s2)/sqrt(s2)))}
# define the integral
I_z <- function(z){P_z <- p_z(z)*(1-p_z(z))^(k-1) * dnorm(z,0,sd=sqrt(s2))}
# estimate the integral on a grid of values
k_max <- max(mt_arr)
k_plt_arr = 1:k_max

p_k_mat <- matrix(NA,length(N_arr),k_max)
# ok for each estimate of the variance of the log-likelihood we need to estimate the respective tau probabilities
for(N_i in 1:length(N_arr)){
  nparticles <- N_arr[N_i]
  s2 <- varlogz_arr[N_i]
  
  lower = -10*sqrt(s2)
  upper = 10*sqrt(s2)
  
  for(i in 1:length(k_plt_arr)){
    k <- k_plt_arr[i]
    p_k_mat[N_i,i] <- integrate(I_z, lower, upper)$value
    
  }
}
# plot the histogram of empirical meeting times
mt_max <- max(mt_arr)
count_arr <- matrix(NA,length(N_arr),mt_max)
for(n_i in 1:length(N_arr)){
  for(i in 1:mt_max){
    count_arr[n_i,i] <- sum(mt_arr[n_i,]==i)
  }
}
prob_arr <- count_arr/dim(mt_arr)[2]
sd_err <- sqrt(prob_arr*(1-prob_arr)/R)
prob_arr_u <- prob_arr + 2*sd_err
prob_arr_l <- prob_arr - 2*sd_err

plot_df <- data.frame(t(prob_arr))
plot_df_u <- data.frame(t(prob_arr_u))
plot_df_l <- data.frame(t(prob_arr_l))

names(plot_df) <- N_arr
names(plot_df_u) <- N_arr
names(plot_df_l) <- N_arr

plot_df$tau_val <- 1:mt_max
plot_df_u$tau_val <- 1:mt_max
plot_df_l$tau_val <- 1:mt_max

plot_df_m <- reshape2::melt(plot_df,id='tau_val')
plot_df_m_u <- reshape2::melt(plot_df_u,id='tau_val')
plot_df_m_l <- reshape2::melt(plot_df_l,id='tau_val')


names(plot_df_m)[2] <- 'N'
plot_df <- cbind(plot_df_m,plot_df_m_u$value,plot_df_m_l$value)
names(plot_df)[4:5] <- c('u','l')

plot_df_exact <- data.frame(t(p_k_mat))
names(plot_df_exact) <- N_arr
plot_df_exact$tau_val <- 1:mt_max

plot_df_exact_m <- reshape2::melt(plot_df_exact,id='tau_val')
names(plot_df_exact_m)[2] <- 'N'
plot_df_exact_m$N <- as.factor(plot_df_exact_m$N)

g<-ggplot() +
  geom_bar(data = plot_df_exact_m,aes(x=tau_val,y=value,fill=N),stat='identity',position='dodge')+
  geom_errorbar(data=plot_df, aes(x=tau_val,y=value,fill=N,ymin=l,ymax=u),position='dodge') +
  xlab(TeX('$n$'))+
  scale_x_continuous(breaks=1:100,limits=c(0,10)) +
  ylab(TeX('$P\\[\\tau = n\\]$'))

g

name <- 'lgssm_estimated_mt.pdf'
ggsave(name,plot=g,width=7,height=5)



### compare quantiles

# compare actual and predicted meeting times
# for predicted values we need to do a series of integrations over the range of tau appearing in the histograms
# define the conditional geometric probability
p_z <- function(z){return(pnorm(-z/sqrt(s2)) + exp(-z+s2/2)*pnorm((z-s2)/sqrt(s2)))}
# define the integral
I_z <- function(z){P_z <- (1-p_z(z))^(k-1) * dnorm(z,0,sd=sqrt(s2))}
# estimate the integral on a grid of values
k_max <- max(mt_arr)
k_plt_arr = 1:k_max

p_k_quantile_mat <- matrix(NA,length(N_arr),k_max)
pm_quantile_mat <- matrix(NA,length(N_arr),k_max)
pm_quantile_sd <- matrix(NA,length(N_arr),k_max)
# ok for each estimate of the variance of the log-likelihood we need to estimate the respective tau probabilities
for(N_i in 1:length(N_arr)){
  nparticles <- N_arr[N_i]
  s2 <- varlogz_arr[N_i]
  
  lower = -10*sqrt(s2)
  upper = 10*sqrt(s2)
  
  for(i in 1:length(k_plt_arr)){
    k <- k_plt_arr[i]
    p_k_quantile_mat[N_i,i] <- integrate(I_z, lower, upper)$value
    pm_quantile_mat[N_i,i] <- sum(pmmh_mts[N_i,]>=k)/length(pmmh_mts[N_i,])
    pm_quantile_sd[N_i,i] <- sqrt(pm_quantile_mat[N_i,i]*(1-pm_quantile_mat[N_i,i]))/sqrt(length(pmmh_mts[N_i,]))
  }
}


delta_offset <- 0.1
n_data <- 15
x1 <- c(1:n_data)-delta_offset
x2 <- c(1:n_data)+delta_offset
pm_quantile_u <- t(pm_quantile_mat+2*pm_quantile_sd)[1:n_data,c(2,4,6)]
pm_quantile_l <- t(pm_quantile_mat-2*pm_quantile_sd)[1:n_data,c(2,4,6)]


ggplot_df1 <- data.frame(t(p_k_quantile_mat[c(2,4,6),1:n_data]),x1)
names(ggplot_df1)[1:3] <- c((N_arr[2]),(N_arr[4]),(N_arr[6]))
ggplot_df1_m <- melt(ggplot_df1,id.vars='x1')

erbar_df1 <- data.frame(t(pm_quantile_mat)[1:n_data,2],pm_quantile_l[1:n_data,1],pm_quantile_u[1:n_data,1],x2)
erbar_df2 <- data.frame(t(pm_quantile_mat)[1:n_data,4],pm_quantile_l[1:n_data,2],pm_quantile_u[1:n_data,2],x2)
erbar_df3 <- data.frame(t(pm_quantile_mat)[1:n_data,6],pm_quantile_l[1:n_data,3],pm_quantile_u[1:n_data,3],x2)
names(erbar_df1)[1:3] <- c('q','q_l','q_u')
names(erbar_df2)[1:3] <- c('q','q_l','q_u')
names(erbar_df3)[1:3] <- c('q','q_l','q_u')
ggplot_df2 <- data.frame(t(pm_quantile_mat[c(2,4,6),1:n_data]),x2)
ggplot_df2_m <- melt(ggplot_df2,id.vars='x2')



breaks <- 10^(-10:10)
minor_breaks <- rep(1:9, 21)*(10^rep(-10:10, each=9))

g <- ggplot() + 
  geom_point(data=ggplot_df1_m,aes(x=x1,y=value,shape=variable),size=2) +
  geom_point(data=ggplot_df2_m,aes(x=x2,y=value),size=0.4) +
  geom_errorbar(data=erbar_df1,aes(x=x2,ymin=q_l,ymax=q_u),width=0.2) +
  geom_errorbar(data=erbar_df2,aes(x=x2,ymin=q_l,ymax=q_u),width=0.2) +
  geom_errorbar(data=erbar_df3,aes(x=x2,ymin=q_l,ymax=q_u),width=0.2) +
  scale_x_continuous(breaks = 0:100,minor_breaks = NULL)+
  scale_y_log10(breaks =breaks ,minor_breaks=minor_breaks) + 
  scale_shape_manual(values=c(17, 16, 3))+
  ylab(TeX('$P\\[\\tau\\geq n\\]$')) +
  xlab('n') +
  guides(shape=guide_legend(title="N"))
g
name <- ('tailproblgssm.pdf')
ggsave(name,plot=g,width=7,height=5)




