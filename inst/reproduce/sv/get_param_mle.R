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
library(latex2exp)
setmytheme()



cores_requested <- min(80,detectCores()-1)
registerDoMC(cores = cores_requested)




parameter_est_file <- 'inst/reproduce/sv/data/parameter_est.RData'
parameter_plot_file <- 'inst/reproduce/sv/figs/parameter_ll.pdf'
obs_plot_file <- 'inst/reproduce/sv/figs/obs.pdf'
parameter_mle_file <- 'inst/reproduce/sv/data/parameter_mle.RData'


sv_data_file <- 'inst/reproduce/sv/sp500.RData'
load(file=sv_data_file)
observations_plt <- observationsDF

g <- ggplot(data = observations_plt, aes(x = index, y = y)) +
  geom_line() + 
  ylab("S&P log-returns\n (re-scaled)")
print(g)
ggsave(file=obs_plot_file,plot=g,width=7,height = 5)

# load and plot
load(parameter_est_file)

plot(theta_grid_mat[1,],ll_theta_mat[1,])
plot(theta_grid_mat[2,],ll_theta_mat[2,])
plot(theta_grid_mat[3,],ll_theta_mat[3,])
plot(theta_grid_mat[4,],ll_theta_mat[4,])
plot(theta_grid_mat[5,],ll_theta_mat[5,])


theta_names <-c('$\\mu$','$\\beta$','$\\xi$','$\\omega^2$','$\\lambda$')
theta_strs <- c("theta1", "theta2", "theta3", "theta4", "theta5")

theta_df <- data.frame(theta_grid_mat[1,],ll_theta_mat[1,])
names(theta_df) <- c('theta','ll')
theta_df$theta_ind <- theta_strs[1]


for(i in 2:D_theta){
  
  theta_df_i <- data.frame(theta_grid_mat[i,],ll_theta_mat[i,])
  names(theta_df_i) <- c('theta','ll')
  theta_df_i$theta_ind <- theta_strs[i]

  theta_df <- rbind(theta_df,theta_df_i)
}

theta_df$theta_ind <- as.factor(theta_df$theta_ind)

# levels(theta_df$theta_ind) <- theta_names



g <- ggplot(data=theta_df) + 
  geom_line(aes(x=theta,y=ll,color=theta_ind)) +
  guides(color=guide_legend(title=NULL)) + 
  scale_color_discrete(labels=lapply(sprintf('%s', theta_names), TeX)) + 
  xlab(TeX('$\\theta$ values')) + 
  ylab(TeX('log-likelihood'))
g               

theta_mle <- rep(NA,D_theta)

for(i in 1:D_theta){
  theta_mle[i] <- theta_grid_mat[i,ll_theta_mat[i,]==max(ll_theta_mat[i,])]
}

theta_mle
save(file=parameter_mle_file,theta_mle)


ggsave(file=parameter_plot_file,plot=g,width=7,height = 5)

# xi = theta[3]
# w2 = theta[4]
# lambda = theta[5]
# 
# R <- 50
# n_L <- 5
# 
# # seems a good range for 20 obs
# theta3_grid <- seq(1e-4,1e-3,length.out=n_L)
# theta4_grid <- seq(1e-7,1e-6,length.out=n_L)
# theta5_grid <- seq(0.01,0.5,length.out=n_L)
# 
# # 
# # # # seems a good range for 40 obs
# # theta3_grid <- seq(1e-4,1e-3,length.out=n_L)
# # theta4_grid <- seq(1e-7,1e-6,length.out=n_L)
# # theta5_grid <- seq(1e-3,1e-2,length.out=n_L)
# 
# 
# # # 80 obs
# # theta3_grid <- seq(1e-4,1e-3,length.out=n_L)
# # theta4_grid <- seq(1e-7,1e-6,length.out=n_L)
# # theta5_grid <- seq(1e-3,1e-2,length.out=n_L)
# 
# # # 160 obs
# # theta3_grid <- seq(1e-4,1e-3,length.out=n_L)
# # theta4_grid <- seq(1e-7,1e-6,length.out=n_L)
# # theta5_grid <- seq(1e-3,1e-2,length.out=n_L)
# # 
# 
# 
# # 
# # 
# # # 640 obs
# # theta3_grid <- seq(1e-4,1e-3,length.out=n_L)
# # theta4_grid <- seq(1e-7,1e-6,length.out=n_L)
# # theta5_grid <- seq(1e-3,1e-2,length.out=n_L)
# # 
# 
# 
# theta_mat <- matrix(NA,n_L^3,5)
# 
# indx  <- 0
# for(i in 1:length(theta3_grid)){
#   for(j in 1:length(theta4_grid)){
#     for(k in 1:length(theta5_grid)){
#       indx <- indx +1
#       theta <- c(0.,0.,theta3_grid[i],theta4_grid[j],theta5_grid[k])
#       theta_mat[indx,] <- theta
#     }
#   }
# }
# 
# 
# cpf_res <- CPF(nparticles, sv_model, theta, observations)
# ll_est <- cpf_res$logz
# 
# n_theta <- dim(theta_mat)[1]
# 
# ll_arr_mean <- rep(NA,n_theta)
# ll_arr_var <- rep(NA,n_theta)
# 
# for(i in 1:n_theta){
#     print(i)
#     theta <- theta_mat[i,]
#     
#     cpf_res <- CPF(nparticles, sv_model, theta, observations)
#     
#     ll_ests <- foreach(r = 1:R) %dorng% {
#         cpf_res <- CPF(nparticles, sv_model, theta, observations)
#         ll_est <- cpf_res$logz
#     }
#     ll_ests <- unlist(ll_ests)
#     ll_arr_mean[i] <- mean(ll_ests)
#     ll_arr_var[i] <- var(ll_ests)
#     
# }
# 
# 
# mat2arr <- function(mat){
#   arr_res <- array(NA,c(length(theta3_grid),length(theta4_grid),length(theta5_grid)))
#   indx <- 0
#   for(i in 1:length(theta3_grid)){
#     for(j in 1:length(theta4_grid)){
#       for(k in 1:length(theta5_grid)){
#         indx <- indx +1
#         arr_res[i,j,k] <- mat[indx]
#       }
#     }
#   }
#   return(arr_res)
# }
# 
# mean_ll <- mat2arr(ll_arr_mean )
# var_ll <- mat2arr(ll_arr_var)
# 
# ll_mask <- ll_arr_var<10
# valid_ll <- ll_arr_mean[ll_mask]
# 
# ll_min <- min(valid_ll)
# ll_max <- max(valid_ll)
# 
# ll_max_mask <-ll_arr_mean==ll_max
# 
# theta_mat[ll_max_mask,]
# 
# theta_opt <- theta_mat[ll_max_mask,]
# 
# 
# 
# n_plt <- 50
# ll_theta3_plt <- rep(NA,n_plt)
# theta3_grid_plt <- seq(theta3_grid[1],theta3_grid[n_L],length.out=n_plt)
# for(r in 1:n_L){
#   
#   theta3 <- theta3_grid_plt[r]
#   theta <- theta_opt
#   theta[3] <- theta3
#   cpf_res <- CPF(nparticles, sv_model, theta, observations)
#   ll_theta3_plt[r]  <- cpf_res$logz
#   
# }
# plot(ll_theta3_plt)
# 
# 
# n_plt <- 50
# ll_theta4_plt <- rep(NA,n_plt)
# theta4_grid_plt <- seq(theta4_grid[1],theta4_grid[n_L],length.out=n_plt)
# for(r in 1:n_plt){
#   
#   theta4 <- theta4_grid_plt[r]
#   theta <- theta_opt
#   theta[4] <- theta4
#   cpf_res <- CPF(nparticles, sv_model, theta, observations)
#   ll_theta4_plt[r]  <- cpf_res$logz
#   
# }
# plot(ll_theta4_plt)
# 
# 
# n_plt <- 20
# ll_theta5_plt <- rep(NA,n_plt)
# theta5_grid_plt <- seq(theta5_grid[1],theta5_grid[n_L],length.out=n_plt)
# for(r in 1:n_plt){
#   
#   theta5 <- theta5_grid_plt[r]
#   theta <- theta_opt
#   theta[5] <- theta5
#   cpf_res <- CPF(nparticles, sv_model, theta, observations)
#   ll_theta5_plt[r]  <- cpf_res$logz
#   
# }
# plot(ll_theta5_plt)
# 
# 
# 
# 
# require(rgl)
# #  open renderer
# open3d()
# #  plot surface
# rgl.surface( 1:5, 1:5 ,  mean_ll[,,2])
# 
# 
# require(plot3D)
# 
# persp3D(z =  mean_ll[,,2], theta = 120)

