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
datalength <- dim(observations_plt)[1]
# # plot data
# g <- ggplot(data = observations_plt, aes(x = index, y = y)) 
# g <- g + geom_line() + ylab("observations")
# print(g)
observations <- as.matrix(observations_plt['y'])

# get model
model <- get_levy_sv()



# indicate whether to use all the particles with a Rao Blackwellised estimator
all_particles <- F
nparticles <- N_arr[1]
R_plot <- 1000

# get the kernels
pimh_kernels <- get_pimh_kernels(with_as,all_particles)
pimh_init <- pimh_kernels$pimh_init
single_kernel_pimh <- pimh_kernels$single_kernel
coupled_kernel_pimh <- pimh_kernels$coupled_kernel


# run the coupling algorithm
t1 <- Sys.time()
pimh_iterations <- foreach(r = 1:R_plot) %dorng% {
  coupled_smoothing_res <- coupled_pimh_chains(single_kernel_pimh, coupled_kernel_pimh, pimh_init,nparticles,K=1,verbose=F)
  
}

# get the unbiased estimators from the chains
h_est <- do.call(rbind,lapply(pimh_iterations,function(x) H_bar(x)))
state_mat <-h_est
mean_state_mat <- (colMeans(state_mat))
var_state <- (diag(var(state_mat)))


# plot results
mean_state <- mean_state_mat[2:(datalength+1)]
sd_err <- sqrt(var_state[2:(datalength+1)]/R_plot)

noisy_obs <- observations[,1]
mean_state_u <- mean_state+3*sd_err[2:(datalength+1)]
mean_state_l <- mean_state-3*sd_err[2:(datalength+1)]

plot_df <- data.frame(mean_state,noisy_obs,mean_state_u,sd_err)
plot_df$indx <- c(1:datalength)

g <- ggplot(data=plot_df)+
  geom_ribbon(aes(x=indx,ymin=mean_state_l,ymax=mean_state_u), fill = "grey70")+ #,linetype='dotted') +
  geom_line(aes(x=indx,y=noisy_obs,group=1),lwd=0.1) +
#  geom_line(aes(x=indx,y=true_state,group=2,color='true')) +
  geom_line(aes(x=indx,y=mean_state,group=3,color='smoothing')) + 
  labs(colour="",x="t",y=TeX("E\\[W_t|Y_{1:T}\\]")) +
  guides(color=FALSE)
g
name <- 'svunbiased_est_withobs.pdf'
ggsave(name,plot=g,width=7,height=5)



g <- ggplot(data=plot_df)+
  geom_ribbon(aes(x=indx,ymin=mean_state_l,ymax=mean_state_u), fill = "grey70")+ #,linetype='dotted') +
  #geom_line(aes(x=indx,y=noisy_obs,group=1),lwd=0.1) +
  #  geom_line(aes(x=indx,y=true_state,group=2,color='true')) +
  geom_line(aes(x=indx,y=mean_state,group=3,color='smoothing')) + 
  labs(colour="",x="t",y=TeX("E\\[W_t|Y_{1:T}\\]")) +
  guides(color=FALSE)
g
name <- 'svunbiased_est.pdf'
ggsave(name,plot=g,width=7,height=5)


