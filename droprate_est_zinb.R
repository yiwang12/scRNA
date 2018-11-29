####### estimate dropour rate using ZINB-WaVE
library(zinbwave)
source("zinbft.R")

inv.logit <- function(a) { # inverse logit function
  ans <- exp(a)/(1+exp(a))
  ans[is.nan(ans)] <- 1  # If a == Inf then ans should be 1
  return(ans)
}

esti_droprate_zinb<-function(raw){
  filter <- rowSums(raw>1)>=1
  raw <- raw[filter,]
  zinb <- zinbFit(raw) 
  X = zinb@X
  V = zinb@V
  W = zinb@W
  a_pi = zinb@alpha_pi
  b_pi = zinb@beta_pi
  g_pi = zinb@gamma_pi
  logit_pi = X %*% b_pi + t(V %*% g_pi) + W %*% a_pi
  pi_pre=t(inv.logit(logit_pi)) # estimated dropout rate of each count Y(i,j)
  droprate_est=array(dim=length(pi_pre[,1])) # mean droprate of each gene(only use non-zero counts)
  for(i in 1:length(droprate_est)){
    d=pi_pre[i,]
    droprate_est[i] = mean(d[which (raw[i, ]!=0)] )
  }
  return(droprate_est)
}




