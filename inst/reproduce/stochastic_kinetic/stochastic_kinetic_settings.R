### stochastic_kinetic experiment settings
stochastic_kinetic_settings <- function(){
  ### settings
  s2_obs <- 1
  A_obs <- rbind(c(1,0,0,0),c(0,1,2,0))
  theta <- list(s2_obs=s2_obs,A_obs=A_obs)
  datalength <- 1e2
  with_as <- FALSE
  N_arr <- c(1000,2000,3000)
  N_arr_large_m <- c(seq(500,3000,500))
  N_arr_hkm <- c(seq(500,3000,500))
  n_var <- 10000#50000
  R <- 10000#50000
  R_bound_check <-10000
  R_N_opt <-500
  N_arr_large_m <- rep(2^c(9:13))#rep(round(1/seq(1/10,1/150,length.out=10)),5)
  
  save_file <- 'inst/reproduce/stochastic_kinetic/data/stochastic_kinetic_varlogz.RData'
  hkm_save_file <- 'stochastic_kinetic_hkm.RData'
  hkm_small_save_file <- 'stochastic_kinetic_hkm_small.RData'
  bound_check_save_file <- 'bound_check_stochastic_kinetic.RData'
  meeting_times_varlogz_save_file <- 'inst/reproduce/stochastic_kinetic/data/stochastic_kinetic_mt_varz.RData'
  data_file <- sprintf('inst/reproduce/stochastic_kinetic/data/data%iobs.RData',datalength)
  results_file <- '/homes/middleto/Documents/DataDrives/CrossbillData/CSMC_res/varlogz_toy.RData'
  bound_check_results_file <- 'inst/reproduce/stochastic_kinetic/data/stochastic_kinetic_bound_check_res.RData'
  N_opt_results_file <- 'inst/reproduce/stochastic_kinetic/data/stochastic_kinetic_N_opt_res.RData'
  var_k_save_file <-  'stochastic_kinetic_var_k.RData'
  m_bc_bound_check <- 100
  m_N_opt <- 100
  n_m_N_opt <- 10
  m_grid_N_opt <- round(seq(500,5000,length.out = n_m_N_opt))
  k_grid_N_opt <- round(m_grid_N_opt/10)
  R_hkm <- 5e2
  R_var_k <- 2e4#2e4
  K_max <- 1024
  R_hkm_large_m <- 2e5
  K_max_large_m <- 1e5
  K_var_k <- 128
  nparticles_var_k <- 2000
  
  return(list(theta=theta,
              datalength=datalength,
              with_as=with_as,
              N_arr=N_arr,
              N_arr_hkm=N_arr_hkm,
              hkm_small_save_file=hkm_small_save_file,
              N_arr_large_m=N_arr_large_m,
              n_var=n_var,
              R=R,
              R_hkm=R_hkm,
              R_hkm_large_m=R_hkm_large_m,
              R_var_k=R_var_k,
              hkm_save_file=hkm_save_file,
              K_max=K_max,
              K_max_large_m=K_max_large_m,
              save_file=save_file,
              bound_check_save_file=bound_check_save_file,
              bound_check_results_file=bound_check_results_file,
              meeting_times_varlogz_save_file=meeting_times_varlogz_save_file,
              var_k_save_file=var_k_save_file,
              results_file=results_file,
              data_file=data_file,
              m_bc_bound_check=m_bc_bound_check,
              m_N_opt=m_N_opt,
              R_bound_check=R_bound_check,
              R_N_opt=R_N_opt,
              m_grid_N_opt=m_grid_N_opt,
              k_grid_N_opt=k_grid_N_opt,
              K_var_k = K_var_k,
              N_opt_results_file=N_opt_results_file,
              nparticles_var_k=nparticles_var_k))  
}
