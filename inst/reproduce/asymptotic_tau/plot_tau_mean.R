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


### this script takes around 30 mins on a 90 core server
s_arr <- seq(0.01,5,length.out=100)

# compare actual and predicted meeting times
# for predicted values we need to do a series of integrations over the range of tau appearing in the histograms
# define the conditional geometric probability
p_z <- function(z){return(pnorm(-z/sqrt(s2)) + exp(-z+s2/2)*pnorm((z-s2)/sqrt(s2)))}

I_mean_tau <- function(z){P_z <- (1/p_z(z)) * dnorm(z,0,sd=sqrt(s2))}
I_var_tau <- function(z){P_z <- ((1-p_z(z))/p_z(z)^2) * dnorm(z,0,sd=sqrt(s2))}

tail_arr <- c(1)
p_mat <- matrix(NA,length(tail_arr),length(s_arr))
tau_mean_arr <- rep(NA,length(s_arr))
tau_var_arr <- rep(NA,length(s_arr))
# ok for each estimate of the variance of the log-likelihood we need to estimate the respective tau probabilities
for(i in 1:length(s_arr)){
  s2 <- s_arr[i]^2
  
  lower = -50*sqrt(s2)
  upper = 50*sqrt(s2)

  for(j in 1:length(tail_arr)){
    tail_val <- tail_arr[j]
    I_z <- function(z){P_z <- (1-(1-p_z(z))^tail_val) * dnorm(z,0,sd=sqrt(s2))}
    p_mat[j,i] <- integrate(I_z, lower, upper)$value
  }
  tau_mean_arr[i] <- integrate(I_mean_tau, lower, upper)$value
  tau_var_arr[i] <- integrate(I_var_tau, lower, upper)$value
}

matplot(s_arr,t(p_mat),pch=1:length(tail_arr),type='l')
plot(s_arr,tau_var_arr,type='l',col='orange',log='y')
plot(s_arr,tau_mean_arr,log='',type='l',col='orange')




p_mat_df <- data.frame(s=s_arr,p=t(p_mat))
g1<-ggplot(data = p_mat_df) +
   geom_line(aes(x=s,y=p))+
   #geom_errorbar(data=plot_df, aes(x=tau_val,y=value,fill=N,ymin=l,ymax=u),position='dodge') +
   xlab(TeX('$\\sigma$'))+
   scale_x_continuous(breaks=seq(0,10,1)) +
   #ylab()
  ylab('')+
  ggtitle(TeX('$P\\[\\tau = 1\\]$'))+
  theme(plot.title = element_text(size = 20),
        #plot.margin = margin(0, 0, 0, 0),
        axis.title.y = element_blank())#, face = "bold"))

g1

p_mat_df <- data.frame(s=s_arr,p=tau_mean_arr)
g2<-ggplot(data = p_mat_df) +
  geom_line(aes(x=s,y=p))+
  #geom_errorbar(data=plot_df, aes(x=tau_val,y=value,fill=N,ymin=l,ymax=u),position='dodge') +
  scale_x_continuous(breaks=seq(0,10,1)) +
  xlab(TeX('$\\sigma$'))+
  ylab('')+
  ggtitle(TeX('$E\\[\\tau\\]$'))+
  theme(plot.title = element_text(size = 20),
        axis.title.y = element_blank())
g2



library(egg)      
g <- ggarrange(g1,g2,ncol=2)
g
# g <- cowplot::plot_grid(g1, g2, labels = c("",""))
# g
name <- ('ptau_summary.pdf')
ggsave(name,plot=g,width=7,height=4)

