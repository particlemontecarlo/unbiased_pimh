### lgssm experiment settings
smc_settings <- function(){
  ### settings
  N <- 100#10^(1:4)
  K <- 64
  k_arr <- seq(0,K,5)#floor(0.1*K_arr)
  n_var_arr <- 10^3#pmax(10^5/N_arr,20)
  R <- 10000
  save_file <- 'smc_estimators.RData'
  meeting_times_varlogz_save_file <- 'inst/reproduce/smcsamplers/data/smcsamplers_mt_varz.RData'
  theta  <- list(xstar=c(-3,0,3,6),ss=1)
  dimension <- 4
  
  return(list(theta=theta,
              dimension=dimension,
              N=N,
              n_var_arr=n_var_arr,
              k_arr=k_arr,
              K=K,
              R=R,
              save_file=save_file,
              meeting_times_varlogz_save_file=meeting_times_varlogz_save_file))  
}







