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


### this script plots some nice things about the meeting times and var logz
### it also produces a small file the array of meeting times (mt_arr),
### the variance of logz (varlogz_arr)  and the grid of N (N_arr)
### this is needed to do the next experiment which gets the estimators on a 
### grid of m


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
results_file <- settings$results_file
data_file <- settings$data_file
meeting_times_varlogz_save_file <- settings$meeting_times_varlogz_save_file

###
cores_requested <- min(80,detectCores()-1)
registerDoMC(cores = cores_requested)

module_tree <<- Module("module_tree", PACKAGE = "CoupledPIMH")
TreeClass <<- module_tree$Tree


K_arr <- 2^(7:11)
k_arr <- round(K_arr/5)


### model & data
model <-  get_ar(dimension)
load(data_file)

# load the results - this can be large!!!!
print('loading results...')
save_file <- 'lgssm_compare_var.RData'
load(save_file)

# # produce the condensed file of meeting times
# print('saving condensed file of meeting times...')
# save(N_arr,mt_arr,varlogz_arr,file=meeting_times_varlogz_save_file)
# #load(meeting_times_varlogz_save_file)

### plot meeting times at different quantiles
# ok, we want to look at distribution of meeting times as a function of sigma2  to estimate how to pick m to look at bias correction
print('plotting results...')
pmmh_mts <- mt_mat_lst[[1]]
mt_arr <- mt_mat_lst[[1]]
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
varlogz_arr <- var_logz_res
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


### plot comparison of variances
library(latex2exp)
s2_names <-c('m(1)','m(2)','m(3)','m(4)','m(5)')
s2_vec1 <- c(rep(s2_names[1],length(N_arr_mat[1,])),
             rep(s2_names[2],length(N_arr_mat[2,])),
             rep(s2_names[3],length(N_arr_mat[3,])),
             rep(s2_names[4],length(N_arr_mat[4,])),
             rep(s2_names[5],length(N_arr_mat[5,])))
s2_vec <- c(s2_vec1,s2_vec1)
N_df$algorithm <- alg_vec
N_df$s2 <- s2_vec

N_df$ineff <- ineff_df$value
N_df$Var2 <- as.factor(N_df$Var2)
N_df$s2 <- as.factor(N_df$s2)
g <- ggplot(data=N_df) + 
  geom_point(aes(x=value,y=ineff,group=Var2,shape=s2),size=2)+ 
  geom_line(aes(x=value,y=ineff,group=Var2,col=algorithm))+ 
  theme(legend.title=element_blank()) +
  scale_y_continuous(trans='log10') +
  xlab('N') + 
  ylab('var(H)')
g

name <- 'comparevariancez.pdf'
ggsave(name,plot=g,width=7,height=5)



var_smoothing_pimh_mat <- do.call(rbind,var_smoothing_pimh_lst)
var_smoothing_ll_mat <- do.call(rbind,var_smoothing_ll_lst)
ineff_smoothing_pimh_mat <- var_smoothing_pimh_mat*N_arr_mat
ineff_smoothing_ll_mat <- var_smoothing_ll_mat*N_arr_mat


### get optimal values of N
Nplt  <- rbind(N_arr_mat,N_arr_mat)
ineffplt <- rbind(var_smoothing_pimh_mat,var_smoothing_ll_mat*180)*Nplt
matplot(t(Nplt),t(ineffplt),log='y',type='l')

opt_N_pimh <- rep(NA,length(K_arr))
opt_N_ll <- rep(NA,length(K_arr))

for(i in 1:length(K_arr)){
  opt_N_pimh[i] <- N_arr[ineff_smoothing_pimh_mat[i,]==min(ineff_smoothing_pimh_mat[i,])]
  opt_N_ll[i] <- N_arr[ineff_smoothing_ll_mat[i,]==min(ineff_smoothing_ll_mat[i,])]
}


### ok we want first a plot of the inefficiency for each value of m, then some uncertainty for the optimal
# m times IF!!


### plot comparison of variances
library(latex2exp)



Nplt  <- rbind(N_arr_mat)
ineffplt <- var_smoothing_pimh_mat*Nplt
matplot(t(Nplt),t(ineffplt),log='y',type='l')


ineffplt_sd <- sqrt(2/(R_arr-1))*var_smoothing_pimh_mat*Nplt
ineffplt_u <- ineffplt + 2*ineffplt_sd
ineffplt_l <- ineffplt - 2*ineffplt_sd

N_df <- melt(t(Nplt))
ineff_df <- melt(t(ineffplt))
ineff_df_u <- melt(t(ineffplt_u))
ineff_df_l <- melt(t(ineffplt_l))

s2_names <-(K_arr)
s2_vec1 <- c(rep(s2_names[1],length(N_arr_mat[1,])),
             rep(s2_names[2],length(N_arr_mat[2,])),
             rep(s2_names[3],length(N_arr_mat[3,])),
             rep(s2_names[4],length(N_arr_mat[4,])),
             rep(s2_names[5],length(N_arr_mat[5,])))
s2_vec <- c(s2_vec1)
ineff_df$algorithm <- rep('state',length(c(N_arr_mat)))
ineff_df$s2 <- s2_vec
ineff_df$N <- N_df$value
ineff_df$m <- as.factor(ineff_df$s2)
ineff_df$u <- ineff_df_u$value
ineff_df$l <- ineff_df_l$value

g <- ggplot(data=ineff_df) + 
  geom_point(aes(x=N,y=value,group=Var2),size=2)+ 
  geom_line(aes(x=N,y=value,group=Var2,color=m))+ 
  geom_errorbar(aes(x=N,ymin=l,ymax=u,group=Var2),width=5)+ 
  # theme(legend.title=element_blank()) +
  scale_y_continuous(trans='log10') +
  xlab('N') + 
  ylab('var(H) x N')
g

name <- 'lgssm_varhkm.pdf'
ggsave(name,plot=g,width=7,height=5)










