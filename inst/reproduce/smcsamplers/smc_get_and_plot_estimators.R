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
source("inst/reproduce/smcsamplers/smc_settings.R")
settings <- smc_settings()
N <- settings$N
K <- settings$K
k_arr <- settings$k_arr
n_var_arr <- settings$n_var_arr
R <- settings$R
save_file <- settings$save_file
meeting_times_varlogz_save_file <- settings$meeting_times_varlogz_save_file
theta <- settings$theta
dimension  <- settings$dimension


###
cores_requested <- min(40,detectCores()-1)
registerDoMC(cores = cores_requested)

module_tree <<- Module("module_tree", PACKAGE = "CoupledPIMH")
TreeClass <<- module_tree$Tree

model <- get_smc()
yobs <- model$gen_data()



# estimate var log z
N_arr <- N
t1 <- Sys.time()
varlogz_arr <- rep(NA,length(N_arr))
for(i in 1:length(N_arr)){
  nparticles <- N_arr[i]
  n_var <- n_var_arr[i]
  
  print(sprintf('getting var logz N=%i',nparticles))
  
  # seems to have issues with the size of the returned value, lets break things down by 10
  log_z_runs <- foreach(r = 1:n_var) %dorng% {
    pf_res <- smc_sampler(nparticles, model, theta,h_func=h_func)
    logz1 <- pf_res$logz
    return(logz1)
  } 
  varlogz_arr[i] <- var(unlist(log_z_runs))
}

print('N arr :')
print(N_arr)
print('var log z')
print(varlogz_arr)
print('Time taken:')
var_est_timing <- Sys.time()-t1
print(var_est_timing)





# get meeting times and estimators
pimh_kernels <- get_pimh_smc_sampler_kernels(h_func)
single_kernel_pimh <- pimh_kernels$single_kernel_smc
coupled_kernel_pimh <- pimh_kernels$coupled_kernel_smc
pmmh_init <-  pimh_kernels$pimh_init_smc

nparticles <- N

print(sprintf('getting estimators N=%i, K=%i',nparticles,K))

t1 <- Sys.time()
pimh_iterations <- foreach(r = 1:R) %dorng% {
  coupled_smoothing_res <- coupled_pimh_chains(single_kernel_pimh, coupled_kernel_pimh, pmmh_init,nparticles,K=K,verbose=F)
}

print('time taken:')
t2 <- Sys.time()-t1
print(t2)



save.image(save_file)


### process results into data frames for plotting
R_plt <- R
n_bootstrap <- 1000
K_arr <- 2^(2:6)
estimator_lst <- list()
var_lst <- list()

k_arr_mat <- matrix(NA,length(K_arr),K_arr[length(K_arr)])
var_arr_mat <- matrix(0,length(K_arr),K_arr[length(K_arr)])
var_arr_mat_u <- matrix(0,length(K_arr),K_arr[length(K_arr)])
var_arr_mat_l <- matrix(0,length(K_arr),K_arr[length(K_arr)])
for(r in 1:length(K_arr)){
  print(sprintf('%i of %i',r,length(K_arr)))
  K_val <- K_arr[r]
  k_arr <- 0:(K_val-1)#seq(0,K_val,by=1)
  h_arr2 <- matrix(NA,length(k_arr),R_plt)
  for(i in 1:length(k_arr)){
    k <- k_arr[i]
    h_arr2[i,] <- sapply(pimh_iterations[1:R_plt],function(x) H_bar(x,k=k,K=K_val))
  }
  
  k_arr_mat[r,1:length(k_arr)] <- k_arr
  var_h2 <- diag(var(t(h_arr2)))
  var_arr_mat[r,1:length(var_h2)] <- var_h2
  
  # get the uncertainty estimates through estimating V[Hkm] with a number of bootstrap samples
  var_h2_tmp_mat <- matrix(0,n_bootstrap,length(var_h2))
  for(i_bs in 1:n_bootstrap){
    h_arr_tmp <- h_arr2[,sample(1:R_plt,replace=T)]
    var_h2_tmp_mat[i_bs,] <- diag(var(t(h_arr_tmp)))
  }
  
  var_arr_mat_u[r,1:length(var_h2)] <- apply(var_h2_tmp_mat,FUN= function(x) quantile(x,0.99),2)
  var_arr_mat_l[r,1:length(var_h2)] <- apply(var_h2_tmp_mat,FUN= function(x) quantile(x,0.01),2)
}

k_plt = k_arr_mat
var_plt = var_arr_mat
var_plt[var_plt==0] = NA
var_arr_mat_u[var_arr_mat_u==0] = NA
var_arr_mat_l[var_arr_mat_l==0] = NA

# expected cost
mt_arr <- sapply(pimh_iterations,function(x) x$meetingtime)
exp_cost <- rep(NA,length(K_arr))
for(r in 1:length(K_arr)){
  K_val <- K_arr[r]
  exp_cost[r] <- mean(2*mt_arr-1 + pmax(1,K_val-mt_arr+1))
}

pdf('kplt.pdf');matplot(t(k_plt),t(var_plt),pch=1:length(K_arr),log='y');dev.off()
pdf('kpltineff.pdf');matplot(t(k_plt),t(var_plt*exp_cost),pch=1:length(K_arr),log='y');dev.off()



# ok two plots:
# 1. optimal k
# 2. var max K against k=K=0
# 3. plot meeting times

###plot results

k_vals <- data.frame(t(k_plt))
names(k_vals) <- K_arr
k_plt_df <- melt(k_vals)

var_vals <- data.frame(t(var_plt))
names(var_vals) <- K_arr
var_plt_df <- melt(var_vals)

var_vals_u <- data.frame(t(var_arr_mat_u))
names(var_vals_u) <- K_arr
var_plt_df_u <- melt(var_vals_u)

var_vals_l <- data.frame(t(var_arr_mat_l))
names(var_vals_l) <- K_arr
var_plt_df_l <- melt(var_vals_l)




var_col <- var_plt_df[,2]
sd <- var_col*sqrt(2/(R_plt-1))
plt_df <- cbind(k_plt_df,var_col,var_plt_df_u[,2],var_plt_df_l[,2])
names(plt_df) <- c('m','k','var','var_u','var_l')


breaks <- 10^(-10:10)
minor_breaks <- rep(1:9, 21)*(10^rep(-10:10, each=9))

g <- ggplot(data=plt_df) +
  #geom_ribbon(aes(x=k,ymin=var_l,ymax=var_u,group=m), fill = "grey70") +
  geom_ribbon(aes(x=k,ymin=var_l,ymax=var_u,group=m,fill=m),alpha=0.3) +
  # geom_ribbon(aes(x=k,ymin=var_l,ymax=var_u,group=m),fill = "grey70") +
  # geom_point(aes(x=k,y=var,shape=m),cex=2) +
  geom_line(aes(x=k,y=var,col=m),size=1.5) +
  scale_y_log10(breaks =breaks ,minor_breaks=minor_breaks,limits=c(8,500)) +
  # scale_x_log10()+
  ylab(TeX('V\\[$H_{k:m}$\\]'))
  #ylim()
g
name <- ('smcvar.pdf')
ggsave(name,plot=g,width=7,height=5)



### plot second plot - (expected cost x variance) for k=m=0 and k,m>0

# get inefficiency for H0:0
var00 <- var(sapply(pimh_iterations[1:R_plt],function(x) H_bar(x,k=0,K=1)))
exp_cost00 <- mean(2*mt_arr-1)
ineff_00 <- var00*exp_cost00

k_vals <- data.frame(t(k_plt))
names(k_vals) <- K_arr
k_plt_df <- melt(k_vals)

var_vals <- data.frame(t(var_plt*exp_cost))
names(var_vals) <- K_arr
var_plt_df <- melt(var_vals)

var_col <- var_plt_df[,2]
sd <- var_col*sqrt(2/(R_plt-1))
plt_df <- cbind(k_plt_df,var_col,var_col+4*sd,var_col-4*sd)
names(plt_df) <- c('m','k','var','var_u','var_l')
plt_df <- plt_df[plt_df$m==64,]


breaks <- 10^(-10:10)
minor_breaks <- rep(1:9, 21)*(10^rep(-10:10, each=9))

g <- ggplot(data=plt_df) +
  geom_errorbar(aes(x=k,ymin=var_l,ymax=var_u),width=1) +
  geom_line(aes(x=k,y=var,col=m)) +
  geom_hline(yintercept=ineff_00, linetype="dashed", color = "red") + 
  ylab(TeX('Var\\[$H_{k:m}$\\]'))+
  xlim(0,20)+
  ylim(750,2500)
g
name <- ('smcineff.pdf')
ggsave(name,plot=g,width=7,height=5)



### plot meeting times 

hist_res <- hist(mt_arr,breaks=0:100,plot=F)
mt_df <- data.frame(breaks=hist_res$breaks[2:101],density_plot=as.numeric(hist_res$density))
g <- ggplot(data=mt_df) + 
  geom_bar(aes(x=breaks,y=density_plot),stat="identity",width=0.5 ) + 
  # scale_y_log10(breaks =breaks ,minor_breaks=minor_breaks)
  scale_x_continuous(breaks=c(0:100),limits=c(0,10))+
  xlab('n') + 
  ylab(TeX('$P\\[\\tau =n\\]$'))
g
name <- ('smcmt.pdf')
ggsave(name,plot=g,width=7,height=5)




