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
source("inst/reproduce/stochastic_kinetic/stochastic_kinetic_settings.R")
settings <- stochastic_kinetic_settings()
theta <- settings$theta
datalength <- settings$datalength
with_as <- settings$with_as
n_var <- settings$n_var
results_file <- settings$results_file
data_file <- settings$data_file


nparticles <- 1000
R_plot <- 1e2

# set up parallel environment
cores_requested <- min(40,detectCores()-1)#min(detectCores()-1)#
print(cores_requested)
registerDoMC(cores = cores_requested)

# set up memory for the particle filter
module_tree <<- Module("module_tree", PACKAGE = "CoupledPIMH")
TreeClass <<- module_tree$Tree


### model & data
model <-  get_stochastic_kinetic()
dimension <- model$dimension
load(data_file)


# indicate whether to use all the particles with a Rao Blackwellised estimator
all_particles <- F

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
true_state <- x[2:(datalength+1),1]
mean_state_u <- mean_state+3*sd_err[2:(datalength+1)]
mean_state_l <- mean_state-3*sd_err[2:(datalength+1)]

plot_df <- data.frame(mean_state,noisy_obs,true_state,mean_state_u,sd_err)
plot_df$indx <- c(1:datalength)

g <- ggplot(data=plot_df)+
  geom_ribbon(aes(x=indx,ymin=mean_state_l,ymax=mean_state_u), fill = "grey70")+ #,linetype='dotted') +
  geom_point(aes(x=indx,y=noisy_obs,group=1),shape=3) +
  geom_line(aes(x=indx,y=true_state,group=2,color='true')) +
  geom_line(aes(x=indx,y=mean_state,group=3,color='smoothing')) + 
  labs(colour="",x="index",y="observation")
g
name <- 'skunbiased_est.pdf'
ggsave(name,plot=g,width=7,height=5)



