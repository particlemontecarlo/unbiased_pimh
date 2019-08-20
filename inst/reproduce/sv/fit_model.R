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


cores_requested <- min(80,detectCores()-1)
registerDoMC(cores = cores_requested)



sv_data_file <- 'inst/reproduce/sv/sp500.RData'
parameter_est_file <- 'inst/reproduce/sv/data/parameter_est.RData'
load(file=sv_data_file)
observations_plt <- observationsDF

g <- ggplot(data = observations_plt, aes(x = index, y = y)) 
g <- g + geom_line() + ylab("observations")
print(g)

observations <- as.matrix(observations_plt['y'])
# set the smoothing theta to the true params
theta <- c(0.25, -0.25, 0.77, 0.15, 0.1)
delta_theta <- c(0.2,0.2,0.2,0.13,0.09)
theta_u <- theta+delta_theta
theta_l <- theta-delta_theta




sv_model <- get_levy_sv()




# particle filtering
module_tree <<- Module("module_tree", PACKAGE = "CoupledPIMH")
TreeClass <<- module_tree$Tree
nparticles <- 50000

cpf_res <- CPF(nparticles, sv_model, theta, observations)
ll_est <- cpf_res$logz

verbose <- F
n_plt <- 40


D_theta <- length(theta)
ll_theta_mat <- matrix(NA,D_theta,n_plt)
theta_grid_mat <- matrix(NA,D_theta,n_plt)

for(i in 1:D_theta){
  print(i)
  ll_theta_plt <- rep(NA,n_plt)
  theta_grid_plt <- seq(theta_l[i],theta_u [i],length.out=n_plt)
  
  
  ll_ests <- foreach(s = 1:n_plt) %dorng% {

    theta_i <- theta_grid_plt[s]
    theta_val <- theta
    theta_val[i] <- theta_i
    cpf_res <- CPF(nparticles, sv_model, theta_val, observations)
    return(cpf_res$logz)
  }

  ll_theta_mat[i,] <- unlist(ll_ests)
  theta_grid_mat[i,] <- theta_grid_plt
}

# save results
save.image(parameter_est_file)


# load and plot
load(parameter_est_file)

plot(theta_grid_mat[1,],ll_theta_mat[1,])
plot(theta_grid_mat[2,],ll_theta_mat[2,])
plot(theta_grid_mat[3,],ll_theta_mat[3,])
plot(theta_grid_mat[4,],ll_theta_mat[4,])
plot(theta_grid_mat[5,],ll_theta_mat[5,])

