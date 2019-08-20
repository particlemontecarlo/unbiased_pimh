# First, some initialization
# remove all objects from R environment
rm(list = ls())
set.seed(17)
if (!require(CoupledPIMH)){
  library(devtools)
  devtools::document()
}

library(doRNG)
library(latex2exp)
setmytheme()
set.seed(17)


### get experiment settings
source("inst/smc/smc_settings.R")
settings <- smc_settings()
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
meeting_times_varlogz_save_file <- settings$meeting_times_varlogz_save_file


###
cores_requested <- min(80,detectCores()-1)
registerDoMC(cores = cores_requested)

module_tree <<- Module("module_tree", PACKAGE = "CoupledPIMH")
TreeClass <<- module_tree$Tree




### model & data
gss <- get_gss()
model <- gss
load(data_file)


# load the results - this can be large!!!!
print('loading results...')
load(save_file)

# produce the condensed file of meeting times
print('saving condensed file of meeting times...')
save(N_arr,mt_arr,varlogz_arr,file=meeting_times_varlogz_save_file)

### plot meeting times at different quantiles
# ok, we want to look at distribution of meeting times as a function of sigma2  to estimate how to pick m to look at bias correction
pmmh_mts <- mt_arr
pm_mt_quantiles1 <- round(apply(pmmh_mts,1,function(x) quantile(x,0.99)))
pm_mt_quantiles2 <- round(apply(pmmh_mts,1,function(x) quantile(x,0.999)))
pm_mt_quantiles3 <- round(apply(pmmh_mts,1,function(x) quantile(x,0.9999)))
pm_mt_max <- apply(pmmh_mts,1,function(x) max(x))

name <- paste('gss_mt.pdf',sep="")
ggplot_df = data.frame(N=N_arr,q0=pm_mt_quantiles1,q1=pm_mt_quantiles2,q2=pm_mt_quantiles3)

ggplot_melt <- reshape2::melt((ggplot_df),id='N')
ggplot_melt$quantile <- as.factor(ggplot_melt$variable)
levels(ggplot_melt$quantile) <- c(0.99,0.999,0.9999)

breaks <- 10^(-10:10)
minor_breaks <- rep(1:9, 21)*(10^rep(-10:10, each=9))

g <- ggplot(data=ggplot_melt) +
  geom_point(aes(x=N,y=value,shape=quantile),cex=2) +
  geom_line(aes(x=N,y=value,linetype=quantile)) +
  # scale_y_log10(breaks =breaks ,minor_breaks=minor_breaks) +
  ylab('meeting time')
g
#  scale_y_log10()
ggsave(name,plot=g,width=7,height=5)
# get hkm estimators




### we will want to plot also the variance of the log-likelihood

### plot meeting times at different quantiles
# ok, we want to look at distribution of meeting times as a function of sigma2  to estimate how to pick m to look at bias correction
N_plt <- N_arr
v_plt <- varlogz_arr

name <- paste('varlogz_gss.pdf',sep="")
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
  theme(legend.title = element_blank())
g


ggsave(name,plot=g,width=7,height=5)





# compare actual and predicted meeting times
# for predicted values we need to do a series of integrations over the range of tau appearing in the histograms
# define the conditional geometric probability
p_z <- function(z){return(pnorm(-z/sqrt(s2)) + exp(-z+s2/2)*pnorm((z-s2)/sqrt(s2)))}
# define the integral
I_z <- function(z){P_z <- p_z(z)*(1-p_z(z))^(k-1) * dnorm(z,0,sd=sqrt(s2))}
# estimate the integral on a grid of values
k_max <- max(mt_arr)
k_arr = 1:k_max

p_k_mat <- matrix(NA,length(N_arr),k_max)
# ok for each estimate of the variance of the log-likelihood we need to estimate the respective tau probabilities
for(N_i in 1:length(N_arr)){
  nparticles <- N_arr[N_i]
  s2 <- varlogz_arr[N_i]
  
  lower = -5*sqrt(s2)
  upper = 5*sqrt(s2)
  
  for(i in 1:length(k_arr)){
    k <- k_arr[i]
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
prob_arr <- count_arr/R
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
  geom_bar(data = plot_df_exact_m,aes(x=tau_val,y=value,fill=N),stat='identity',position='dodge') + #,stat='identity',position='dodge')+
  # geom_bar(data=plot_df, aes(x=tau_val,y=value,fill=N),stat='identity',position='dodge')+
  geom_errorbar(data=plot_df, aes(x=tau_val,y=value,fill=N,ymin=l,ymax=u),position='dodge')+#,position='dodge') +
  xlab(TeX('$n$'))+
  scale_x_continuous(breaks=1:100,limits=c(0,9)) +
  scale_y_log10(limits=c(2e-4,1))+
  ylab(TeX('$P\\[\\tau = n\\]$'))
g

name <- 'gss_estimated_mt.pdf'
ggsave(name,plot=g,width=7,height=5)




g<-ggplot() +
  geom_bar(data = plot_df_exact_m,aes(x=tau_val,y=value,fill=N),stat='identity',position='dodge') + #,stat='identity',position='dodge')+
  # geom_bar(data=plot_df, aes(x=tau_val,y=value,fill=N),stat='identity',position='dodge')+
  geom_errorbar(data=plot_df, aes(x=tau_val,y=value,fill=N,ymin=l,ymax=u),position='dodge')+#,position='dodge') +
  xlab(TeX('$n$'))+
  scale_x_continuous(breaks=1:100,limits=c(0,9)) +
  # scale_y_log10(limits=c(2e-4,1))+
  ylab(TeX('$P\\[\\tau = n\\]$'))
g

name <- 'gss_estimated_mt_nonlog.pdf'
ggsave(name,plot=g,width=7,height=5)










# 
# # ok we want to compare the expectation of the squared bias correction with the proposed bound
# BC <- function (c_chains, h = function(x) x, k = 0) 
# {
#   maxiter <- c_chains$iteration
#   if (k > maxiter) {
#     print("error: k has to be less than the horizon of the coupled chains")
#     return(NULL)
#   }
# 
#   p <- length(h(c_chains$samples1[1, ]))
#   h_of_chain <- apply(X = c_chains$samples1[,, drop = F], MARGIN = 1, FUN = h)
#   
#   if (is.null(dim(h_of_chain))) {
#     h_of_chain <- matrix(h_of_chain, ncol = 1)
#   }
#   else {
#     h_of_chain <- t(h_of_chain)
#   }
#   
#   H_bar <- apply(X = h_of_chain, MARGIN = 2, sum)
#   if (c_chains$meetingtime <= k + 1) {
#   }
#   else {
#     deltas <- matrix(0, nrow = maxiter - k + 1, ncol = p)
#     deltas_term <- rep(0, p)
#     for (t in 0:(c_chains$meetingtime - 1)) {
#       coefficient <- min(t - k + 1)
#       delta_tp1 <- h(c_chains$samples1[t + 1 + 1, ]) - 
#         h(c_chains$samples2[t + 1, ])
#       deltas_term <- deltas_term + coefficient * delta_tp1
#     }
#     H_bar <- H_bar + deltas_term
#   }
#   return(H_bar)#/(K - k + 1))
# }
# 
# 
# 
# 
# h1 <- function(x){
#   x[1]+x[1]^2+x[datalength+1]+x[datalength+1]^2
# }
# 
# h2 <- function(x){
#   sum(x)+sum(x^2)
# }
# 
# # want a bounded function
# h3 <- function(x){
#   1*(x[1]<0)
# }
# 
# h1_arr <- matrix(NA,length(N_arr),R)
# for(i in 1:length(N_arr)){
#   h1_arr[i,] <- sapply(estimators_lst[[i]],function(x) BC(x,h=h1))
# }
# 
# 
# h2_arr <- matrix(NA,length(N_arr),R)
# for(i in 1:length(N_arr)){
#   h2_arr[i,] <- sapply(estimators_lst[[i]],function(x) BC(x,h=h2))
# }
# 
# h3_arr <- matrix(NA,length(N_arr),R)
# for(i in 1:length(N_arr)){
#   h3_arr[i,] <- sapply(estimators_lst[[i]],function(x) BC(x,h=h3))
# }
# 
# 
# 
# 
# 
# ### plot comparison of variance for k=m=1
# # var_ests <-  diag(var(t(h2_arr)))
# var_ests <- rowMeans(h3_arr^2)
# 
# 
# var_proxy_bound <- rep(NA,length(N_arr))
# for(i in 1:length(N_arr)){
#   mt_i <- mt_arr[i,]
#   var_proxy_bound[i] <- mean((1+mt_i)^2*(2+mt_i)^2)*mean(mt_i>1)
# }
# 
# #var_proxy <-rowMeans((2*mt_arr-1)^2)#-rowMeans(mt_arr))#(2*rowMeans(mt_arr)-1)+4*(diag(var(t(mt_arr))))#+2*sqrt(diag(var(t(mt_arr)))))
# comp_time <- rowMeans(mt_arr)
# 
# 
# # we consider the bound and three test functions
# ylim_arr <- c(0.1*min(c(var_ests,var_proxy_bound)),2*max(c(var_ests,var_proxy_bound)))
# plot(N_arr,var_ests,ylim=ylim_arr,log='y')
# matplot(N_arr,var_proxy_bound,add=T,pch=2)
# 
# plot(N_arr,sqrt(var_ests)*N_arr,log='y')
# matplot(N_arr,sqrt(var_proxy_bound)*N_arr,add=T,pch=2)
# 
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
# s2_arr <- seq(0.2,2,length.out=n_s2)
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



