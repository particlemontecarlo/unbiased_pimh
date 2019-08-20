### lgssm experiment settings
lgssm_settings <- function(){
  ### settings
  dimension <- 1
  theta <- c(0.5,sqrt(10))
  datalength <- 100
  with_as <- FALSE
  N_arr <- seq(20,220,40)/2
  N_arr_large_m <- rep(round(1/seq(1/10,1/150,length.out=10)),5)#rep(2^c(9:13),5)#
  n_var <- 1e4
  R <- 100000#50000
  R_bound_check <-10000
  R_N_opt <-500
  save_file <- 'varlogz_toy.RData'
  bound_check_save_file <- 'bound_check_lgssm.RData'
  meeting_times_varlogz_save_file <- 'inst/lgssm/data/lgssm_mt_varz.RData'
  data_file <- sprintf('inst/lgssm/data/data%iobs.RData',datalength)
  results_file <- '/homes/middleto/Documents/DataDrives/CrossbillData/CSMC_res/varlogz_toy.RData'
  bound_check_results_file <- 'inst/lgssm/data/lgssm_bound_check_res.RData'
  N_opt_results_file <- 'inst/lgssm/data/lgssm_N_opt_res.RData'
  m_bc_bound_check <- 100
  m_N_opt <- 100
  n_m_N_opt <- 10
  m_grid_N_opt <- round(seq(500,5000,length.out = n_m_N_opt))
  k_grid_N_opt <- round(m_grid_N_opt/10)
  
  
  N_arr_k0 <- 2^c(4:11)
  R_k0 <- 1e5
  k0_save_file <- 'inst/lgssm/data/lgssm_k0_comapare_vinf.RData'
  n_var_k0 <- sapply(round(1e6/N_arr_k0),function(x) min(x,1000))
  
  return(list(dimension=dimension,
              theta=theta,
              datalength=datalength,
              with_as=with_as,
              N_arr=N_arr,
              N_arr_large_m=N_arr_large_m,
              n_var=n_var,
              R=R,
              save_file=save_file,
              bound_check_save_file=bound_check_save_file,
              bound_check_results_file=bound_check_results_file,
              meeting_times_varlogz_save_file=meeting_times_varlogz_save_file,
              results_file=results_file,
              data_file=data_file,
              m_bc_bound_check=m_bc_bound_check,
              m_N_opt=m_N_opt,
              R_bound_check=R_bound_check,
              R_N_opt=R_N_opt,
              m_grid_N_opt=m_grid_N_opt,
              k_grid_N_opt=k_grid_N_opt,
              N_opt_results_file=N_opt_results_file,
              
              N_arr_k0=N_arr_k0,
              R_k0=R_k0,
              k0_save_file=k0_save_file,
              n_var_k0=n_var_k0))  
}
