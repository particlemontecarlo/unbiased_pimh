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


N_arr <- c(200,3000)
R_plot <- 1e2
###
cores_requested <- min(40,detectCores()-1)#min(detectCores()-1)#
print(cores_requested)

registerDoMC(cores = cores_requested)

module_tree <<- Module("module_tree", PACKAGE = "CoupledPIMH")
TreeClass <<- module_tree$Tree


### model & data
model <-  get_stochastic_kinetic()
dimension <- model$dimension
load(data_file)



# indicate whether to use all the particles with a Rao Blackwellised estimator
# If true then we return both with and without Rao Blackwellisation
all_particles <- T

# get the kernels
pimh_kernels <- get_pimh_kernels(with_as,all_particles)
pimh_init <- pimh_kernels$pimh_init
single_kernel_pimh <- pimh_kernels$single_kernel
coupled_kernel_pimh <- pimh_kernels$coupled_kernel


var_state <- matrix(NA,length(N_arr),datalength+1)
var_state_rb <- matrix(NA,length(N_arr),datalength+1)
mean_state_mat <- matrix(NA,length(N_arr),datalength+1)
mean_state_rb_mat <- matrix(NA,length(N_arr),datalength+1)
for(i in 1:length(N_arr)){
  print(i)
  nparticles <- N_arr[i]
  
  # run the coupling kernels
  t1 <- Sys.time()
  pimh_iterations <- foreach(r = 1:R_plot) %dorng% {
    coupled_smoothing_res <- coupled_pimh_chains(single_kernel_pimh, coupled_kernel_pimh, pimh_init,nparticles,K=1,verbose=F)
    
  }
  
  # estimators without RB
  h_est <- do.call(rbind,lapply(pimh_iterations,function(x) H_bar(x)))
  
  # estimators with RB
  for(r in 1:R_plot){
    pimh_iterations[[r]]$samples1 <- pimh_iterations[[r]]$exp_samples1
    pimh_iterations[[r]]$samples2 <- pimh_iterations[[r]]$exp_samples2
    
  }
  h_est_rb <- do.call(rbind,lapply(pimh_iterations,function(x) H_bar(x)))
  
  
  state_mat <-h_est
  mean_state_mat[i,] <- (colMeans(state_mat))
  var_state[i,] <- (diag(var(state_mat)))
  state_mat <-h_est_rb
  mean_state_rb_mat[i,] <- (colMeans(state_mat))
  var_state_rb[i,] <- (diag(var(state_mat)))
}


var_state_plt <- var_state[c(1,length(N_arr)),]
var_df <- data.frame(2:(datalength+1),t(var_state_plt[,2:(datalength+1)]))
names(var_df) <- c('indx',toString(N_arr[1]),toString(N_arr[2]))
var_df <- melt(var_df,id.vars='indx')
var_df$variable <- as.factor(var_df$variable)

var_state_rb_plt <- var_state_rb[c(1,length(N_arr)),]
var_df_rb <- data.frame(2:(datalength+1),t(var_state_rb_plt[,2:(datalength+1)]))
names(var_df_rb) <- c('indx',paste(toString(N_arr[1]),'(RB)'),paste(toString(N_arr[2]),'(RB)'))
var_df_rb <- melt(var_df_rb,id.vars='indx')
var_df_rb$variable <- as.factor(var_df_rb$variable)

var_df$method <- 'standard'
var_df_rb$method <- 'rb'


var_df_plt <- rbind(var_df,var_df_rb)
g <- ggplot(data=var_df_plt) + 
  geom_line(aes(x=indx,y=value,color=variable)) +
  scale_y_log10()+
  ylab(TeX("Variance $E\\[X_{t}|Y_{1:T}\\]$"))+
  xlab("index")+
  theme(legend.title=element_blank())
g
name <- 'sk_rb_comparison.pdf'
ggsave(name,plot=g,width=7,height=5)

