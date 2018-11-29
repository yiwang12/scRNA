

################################################################################################################################################################################################
source("functions_mv_cal_v3.R")
library(matrixStats)

plot_compare_drop_num <- function(drop_num_est_in){
  png("compare_drop_num.png")
  plot(drop_num_est_in$drop_num, drop_num_est_in$drop_num_est,main="dropout number of each gene",xlab="true",ylab="estimated")
  abline(0,1)
  dev.off()
}

plot_compare_m_v_beforedrop_est<-function(mean_disp_in,titlename,filename){
  png(paste(filename,"_mean.png",sep="",collapse=""))
  title = paste(titlename," mean comparison",sep="",collapse="")
  plot(mean_disp_in$mean_beforedrop,mean_disp_in$mean_est,xlab="mean(log(cpm)) before dropout",ylab="mean(log(cpm)) estimated" ,main=title)
  abline(0,1)
  dev.off()
  png(paste(filename,"_disp.png",sep="",collapse=""))
  title = paste(titlename,"dispersion comparison",sep="",collapse="")
  plot(mean_disp_in$disp_beforedrop,mean_disp_in$disp_est,xlab="disp(log(cpm)) before dropout",ylab="disp(log(cpm)) estimated" ,main=title)
  abline(0,1)
  dev.off()
}


plot_mv<-function(mean_disp_in,titlename,filename){
  png(paste(filename,".png",sep="",collapse=""))
  title = paste(titlename," mean-dispersion",sep="",collapse="")
  plot(mean_disp_in$mean_beforedrop,mean_disp_in$disp_beforedrop,xlab="mean(log(cpm))",ylab="dispersion",main=title,col="blue",xlim=range(mean_disp_in$mean_beforedrop,mean_disp_in$mean_est ),ylim=range(mean_disp_in$disp_est,mean_disp_in$disp_beforedrop))
  points(mean_disp_in$mean_est,mean_disp_in$disp_est ,col="red")
  legend("topleft",1.4,legend=c("before dropout", "estimated"),col=c("blue","red"), lty=c(2,2), cex=0.8)
  dev.off()
}


plot_lowess_mv <- function(mean_disp_in,titlename,filename,span){
  png(paste(filename,".png",sep="",collapse=""))
  title = paste(titlename," mean-dispersion",sep="",collapse="")
  plot(lowess(mean_disp_in$mean_beforedrop,mean_disp_in$disp_beforedrop,f=span),main=title,xlab="mean(log(cpm))",ylab="dispersion",col="blue",xlim=range(mean_disp_in$mean_beforedrop,mean_disp_in$mean_est ),ylim=range(mean_disp_in$disp_est,mean_disp_in$disp_beforedrop))
  lines(lowess(mean_disp_in$mean_est,mean_disp_in$disp_est,f=span),col="red")
  legend("topleft",1.4,legend=c("before dropout", "estimated"),col=c("blue","red"), lty=c(2,2), cex=0.8)
  dev.off()
}

###### plot all figures
mv_cal <- function(sim_counts){
  out=list()
  drop_num_est = get_est_dropout_num(sim_counts)
  sim_counts_n_drop = norm_counts_3(sim_counts$drop_counts )
  sim_counts_n_beforedrop = norm_counts_3(sim_counts$real_counts )
  mean_disp_knowndrop = mv_knowndrop_norm4(sim_counts,sim_counts_n_beforedrop,sim_counts_n_drop)
  mean_disp_unknowndrop = mv_unknowndrop_est_3(sim_counts,sim_counts_n_beforedrop,sim_counts_n_drop,drop_num_est$drop_num_est)
  out$drop_num_est = drop_num_est
  out$sim_counts_n_drop=sim_counts_n_drop
  out$sim_counts_n_beforedrop = sim_counts_n_beforedrop
  out$mean_disp_knowndrop=mean_disp_knowndrop
  out$mean_disp_unknowndrop=mean_disp_unknowndrop
  return(out)
}



plot_all_figs<-function(sim_counts,mean_disp_knowndrop,mean_disp_unknowndrop,drop_num_est,span_in){
  plot_mv(mean_disp_in=mean_disp_knowndrop,titlename="(dropout known) ", filename="mv_knowndrop") #known dropout 
  plot_mv(mean_disp_in=mean_disp_unknowndrop,titlename="(dropout unknown) ",filename= "mv_unknowndrop") #unknown dropout 
  plot_lowess_mv(mean_disp_in=mean_disp_knowndrop,titlename="(dropout known) ",filename="mv_lowess_knowndrop",span=span_in) #known dropout 
  plot_lowess_mv(mean_disp_in=mean_disp_unknowndrop,titlename="(dropout unknown) ",filename="mv_lowess_unknowndrop",span=span_in) #unknown dropout 
  plot_compare_m_v_beforedrop_est(mean_disp_in=mean_disp_knowndrop,titlename="(dropout known) ",filename="compare_knowndrop") #known dropout 
  plot_compare_m_v_beforedrop_est(mean_disp_in=mean_disp_unknowndrop,titlename="(dropout unknown) ",filename="compare_unknowndrop") #unknown dropout 
  plot_compare_drop_num(drop_num_est_in=drop_num_est)
}  

cal_plot<-function(sim_counts,span_in){
  out=mv_cal(sim_counts)
  mean_disp_knowndrop=out$mean_disp_knowndrop
  mean_disp_unknowndrop=out$mean_disp_unknowndrop
  drop_num_est = out$drop_num_est
  plot_all_figs(sim_counts,mean_disp_knowndrop,mean_disp_unknowndrop,drop_num_est,span_in)
}


### apply in 10 datasets
setwd("/Users/yiwang/Dropbox/yiwang/project_Kasper/3_29/output/allen")
load("/Users/yiwang/Dropbox/yiwang/project_Kasper/2_16/sim_data/allen_sim_counts4.RData")
cal_plot(allen_sim_counts,0.2)

dir.create("/Users/yiwang/Dropbox/yiwang/project_Kasper/3_29/output/maits")
setwd("/Users/yiwang/Dropbox/yiwang/project_Kasper/3_29/output/maits")
load("/Users/yiwang/Dropbox/yiwang/project_Kasper/2_16/sim_data/maits_sim_counts4.RData")
cal_plot(maits_sim_counts,0.2)

dir.create("/Users/yiwang/Dropbox/yiwang/project_Kasper/3_29/output/islam")
setwd("/Users/yiwang/Dropbox/yiwang/project_Kasper/3_29/output/islam")
load("/Users/yiwang/Dropbox/yiwang/project_Kasper/2_16/sim_data/islam_sim_counts4.RData")
cal_plot(islam_sim_counts,0.2)

dir.create("/Users/yiwang/Dropbox/yiwang/project_Kasper/3_29/output/Trapnell")
setwd("/Users/yiwang/Dropbox/yiwang/project_Kasper/3_29/output/Trapnell")
load("/Users/yiwang/Dropbox/yiwang/project_Kasper/2_16/sim_data/Trapnell_sim_counts4.RData")
cal_plot(Trapnell_sim_counts,0.2)

dir.create("/Users/yiwang/Dropbox/yiwang/project_Kasper/3_29/output/Klein")
setwd("/Users/yiwang/Dropbox/yiwang/project_Kasper/3_29/output/Klein")
load("/Users/yiwang/Dropbox/yiwang/project_Kasper/2_16/sim_data/Klein_sim_counts4.RData")
cal_plot(Klein_sim_counts,0.2)

dir.create("/Users/yiwang/Dropbox/yiwang/project_Kasper/3_29/output/Tung")
setwd("/Users/yiwang/Dropbox/yiwang/project_Kasper/3_29/output/Tung")
load("/Users/yiwang/Dropbox/yiwang/project_Kasper/2_16/sim_data/Tung_sim_counts4.RData")
cal_plot(Tung_sim_counts,0.2)

dir.create("/Users/yiwang/Dropbox/yiwang/project_Kasper/3_29/output/Zeisel")
setwd("/Users/yiwang/Dropbox/yiwang/project_Kasper/3_29/output/Zeisel")
load("/Users/yiwang/Dropbox/yiwang/project_Kasper/2_16/sim_data/Zeisel_sim_counts4.RData")
cal_plot(Zeisel_sim_counts,0.2)

dir.create("/Users/yiwang/Dropbox/yiwang/project_Kasper/3_29/output/Eagel")
setwd("/Users/yiwang/Dropbox/yiwang/project_Kasper/3_29/output/Eagel")
load("/Users/yiwang/Dropbox/yiwang/project_Kasper/2_16/sim_data/Eagel_sim_counts4.RData")
cal_plot(Eagel_sim_counts,0.2)

dir.create("/Users/yiwang/Dropbox/yiwang/project_Kasper/3_29/output/patel")
setwd("/Users/yiwang/Dropbox/yiwang/project_Kasper/3_29/output/patel")
load("/Users/yiwang/Dropbox/yiwang/project_Kasper/2_16/sim_data/patel_sim_counts4.RData")
cal_plot(patel_sim_counts,0.2)

dir.create("/Users/yiwang/Dropbox/yiwang/project_Kasper/3_29/output/espresso")
setwd("/Users/yiwang/Dropbox/yiwang/project_Kasper/3_29/output/espresso")
load("/Users/yiwang/Dropbox/yiwang/project_Kasper/2_16/sim_data/espresso_sim_counts4.RData")
cal_plot(espresso_sim_counts,0.2)










