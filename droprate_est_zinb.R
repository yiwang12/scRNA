####### ZINB-WaVE
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


load("/Users/yiwang/Dropbox/yiwang/project_Kasper/2_16/sim_data/allen_sim_counts6.RData")
allen_zinb_droprate_est = esti_droprate_zinb(allen_sim_counts$drop_counts)
save(allen_zinb_droprate_est,file="/Users/yiwang/Dropbox/yiwang/project_Kasper/4_5/output/allen_zinb_droprate_est.RData")

load("/Users/yiwang/Dropbox/yiwang/project_Kasper/2_16/sim_data/maits_sim_counts6.RData")
maits_zinb_droprate_est = esti_droprate_zinb(maits_sim_counts$drop_counts)
save(maits_zinb_droprate_est,file="/Users/yiwang/Dropbox/yiwang/project_Kasper/4_5/output/maits_zinb_droprate_est.RData")

load("/Users/yiwang/Dropbox/yiwang/project_Kasper/2_16/sim_data/islam_sim_counts6.RData")
islam_zinb_droprate_est = esti_droprate_zinb(islam_sim_counts$drop_counts)
save(islam_zinb_droprate_est,file="/Users/yiwang/Dropbox/yiwang/project_Kasper/4_5/output/islam_zinb_droprate_est.RData")

load("/Users/yiwang/Dropbox/yiwang/project_Kasper/2_16/sim_data/Trapnell_sim_counts6.RData")
trapnell_zinb_droprate_est = esti_droprate_zinb(Trapnell_sim_counts$drop_counts)
save(trapnell_zinb_droprate_est,file="/Users/yiwang/Dropbox/yiwang/project_Kasper/4_5/output/trapnell_zinb_droprate_est.RData")

load("/Users/yiwang/Dropbox/yiwang/project_Kasper/2_16/sim_data/Klein_sim_counts6.RData")
klein_zinb_droprate_est = esti_droprate_zinb(Klein_sim_counts$drop_counts)
save(klein_zinb_droprate_est,file="/Users/yiwang/Dropbox/yiwang/project_Kasper/4_5/output/klein_zinb_droprate_est.RData")

load("/Users/yiwang/Dropbox/yiwang/project_Kasper/2_16/sim_data/Tung_sim_counts6.RData")
tung_zinb_droprate_est = esti_droprate_zinb(Tung_sim_counts$drop_counts)
save(tung_zinb_droprate_est,file="/Users/yiwang/Dropbox/yiwang/project_Kasper/4_5/output/tung_zinb_droprate_est.RData")

load("/Users/yiwang/Dropbox/yiwang/project_Kasper/2_16/sim_data/Zeisel_sim_counts6.RData")
zeisel_zinb_droprate_est = esti_droprate_zinb(Zeisel_sim_counts$drop_counts)
save(zeisel_zinb_droprate_est,file="/Users/yiwang/Dropbox/yiwang/project_Kasper/4_5/output/zeisel_zinb_droprate_est.RData")

load("/Users/yiwang/Dropbox/yiwang/project_Kasper/2_16/sim_data/Eagel_sim_counts6.RData")
eagel_zinb_droprate_est = esti_droprate_zinb(Eagel_sim_counts$drop_counts)
save(eagel_zinb_droprate_est,file="/Users/yiwang/Dropbox/yiwang/project_Kasper/4_5/output/eagel_zinb_droprate_est.RData")

load("/Users/yiwang/Dropbox/yiwang/project_Kasper/2_16/sim_data/patel_sim_counts6.RData")
patel_zinb_droprate_est = esti_droprate_zinb(patel_sim_counts$drop_counts)#error
save(patel_zinb_droprate_est,file="/Users/yiwang/Dropbox/yiwang/project_Kasper/4_5/output/patel_zinb_droprate_est.RData")

load("/Users/yiwang/Dropbox/yiwang/project_Kasper/2_16/sim_data/espresso_sim_counts6.RData")
espresso_zinb_droprate_est = esti_droprate_zinb(espresso_sim_counts$drop_counts)
save(espresso_zinb_droprate_est,file="/Users/yiwang/Dropbox/yiwang/project_Kasper/4_5/output/espresso_zinb_droprate_est.RData")




drop_num = get_est_dropout_num(allen_sim_counts)
allen_sim_counts2 = allen_sim_counts
allen_sim_counts2$dropoutrates3 = droprate_est

drop_num_est2 = get_est_dropout_num(allen_sim_counts2)


out=list()
sim_counts = allen_sim_counts
#drop_num_est = get_est_dropout_num(sim_counts)
sim_counts_n_drop = norm_counts_3(sim_counts$drop_counts )
sim_counts_n_beforedrop = norm_counts_3(sim_counts$real_counts )
mean_disp_knowndrop = mv_knowndrop_norm4(sim_counts,sim_counts_n_beforedrop,sim_counts_n_drop)
mean_disp_unknowndrop = mv_unknowndrop_est_3(sim_counts,sim_counts_n_beforedrop,sim_counts_n_drop,drop_num_est2$drop_num_est)
#mean_disp_unknowndrop_truerate = mv_unknowndrop_est_3(sim_counts,sim_counts_n_beforedrop,sim_counts_n_drop,drop_num$drop_num_est)



