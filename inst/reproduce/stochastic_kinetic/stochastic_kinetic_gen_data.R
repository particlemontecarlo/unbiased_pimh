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





model <-  get_stochastic_kinetic()
dimension <- model$dimension
dimension_y <- model$dimension_y
# generate data


nparticles <- 1
x <- matrix(0, nrow = datalength + 1, ncol = dimension)
observations <- matrix(0, nrow = datalength, ncol = dimension_y)
x[1,] <- model$rinit(1, thdeta, model$rinit_rand(1, theta),precomputed)
for (i in 1:datalength){
  x[i+1,] <- model$rtransition(t(x[i,,drop=F]), theta, i, model$rtransition_rand(1, theta),precomputed)
  observations[i,] <- theta$A_obs %*% t(x[i+1,,drop=F]) + sqrt(theta$s2_obs)*rnorm(dimension_y)
}

# plot the data

plot(1:datalength, observations[,2], pch=4)
matplot(0:datalength, x[,2]+ 2*x[,3], type='l',add=T,lwd=10)
#

# run the particle filter to obtain a sample trajector
nparticles <- 2000

n_var <- 100
logz_arr <- rep(NA,n_var)
for(i in 1:n_var){
  res <- CPF(nparticles, model, theta, observations)
  res$logz
  print(res$logz)
  logz_arr[i] <- res$logz
  matplot(0:datalength,res$trajectory[2,] + 2*res$trajectory[3,],type='l',col=rgb(0,0,1,0.2),add=T)
  
}

var(logz_arr)

# save results
save('x','observations',file=data_file)



