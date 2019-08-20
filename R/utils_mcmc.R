# MCMC helpers ------------------------------------------------------------
#'@rdname acc_rate
#'@title acc_rate
#'@description mcmc acceptance rate
#'@export
acc_rate <- function(chain){
  if(is.vector(chain)){
    return (sum(diff(chain)!=0)/length(chain))
  }else{
    return (sum(apply(diff(chain)!=0,1,any))/nrow(chain))
  }
}