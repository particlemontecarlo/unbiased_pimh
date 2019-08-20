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
cores_requested <- min(80,detectCores()-1)
registerDoMC(cores = cores_requested)

module_tree <<- Module("module_tree", PACKAGE = "CoupledPIMH")
TreeClass <<- module_tree$Tree




### model & gen data
model <-  get_ar(dimension)

# generate data
precomputed <- model$precompute(theta)


x <- matrix(0, nrow = datalength + 1, ncol = dimension)
observations <- matrix(0, nrow = datalength, ncol = dimension)
x[1,] <- model$rinit(1, thdeta, model$rinit_rand(1, theta),precomputed)
for (i in 1:datalength){
  x[i+1,] <- model$rtransition(x[i,], theta, i, model$rtransition_rand(1, theta),precomputed)
  observations[i,] <- x[i+1,]+ theta[2]*rnorm(dimension)
}

# plot the data
plot(0:datalength, x[,1], type = "l")
plot(1:datalength, observations[,1], type = "l")
#

# run the particle filter to obtain a sample trajector
nparticles <- 100
xref <- CPF(nparticles, model, theta, observations)$trajectory
plot(xref[1,],type='l')

# save results
save('x','observations',file=data_file)



