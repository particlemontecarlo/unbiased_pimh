### sv experiment settings
sv_settings <- function(){
  ### settings
  with_as <- FALSE
  N_arr <- seq(100,500,100)
  n_var <- 1e4#50000
  R <- 1e4#50000
  varsave_file <- 'inst/reproduce/sv/data/varlogz_sv.RData'
  hkm_save_file <- 'inst/reproduce/sv/data/sv_hkm.RData'
  hkm_small_save_file <- 'inst/reproduce/sv/data/sv_hkm_reduced.RData'
  meeting_times_varlogz_save_file <- 'sv_mt_varz.RData'
  sv_data_file <- 'inst/reproduce/sv/sp500.RData'
  estimators_res_file <- 'inst/reproduce/sv/data/estimator_res.RData'
  parameter_mle_file <- 'inst/reproduce/sv/data/parameter_mle.RData'
  R_hkm <- 1e4
  K_max <- 512
  
  
  N_arr_k0 <- 2^c(7:15)
  R_k0 <- sapply(round(1e7/N_arr_k0),function(x) min(x,1000))
  k0_save_file <- 'inst/sv/data/sv_k0.RData'
  n_var_k0 <- round(200*1e4/N_arr_k0)
  
  return(list(with_as=with_as,
              N_arr=N_arr,
              N_arr_k0=N_arr_k0,
              R_k0=R_k0,
              k0_save_file=k0_save_file,
              n_var=n_var,
              n_var_k0=n_var_k0,
              R=R,
              R_hkm=R_hkm,
              K_max=K_max,
              hkm_save_file=hkm_save_file,
              hkm_small_save_file=hkm_small_save_file,
              varsave_file=varsave_file,
              meeting_times_varlogz_save_file=meeting_times_varlogz_save_file,
              sv_data_file=sv_data_file,
              estimators_res_file=estimators_res_file,
              parameter_mle_file=parameter_mle_file))  
}
