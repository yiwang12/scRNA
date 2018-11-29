# plot mean-vairance relationship
source("functions_mv_cal_v3.R")
library(matrixStats)

# dropout number estimation vs. true dropout number
plot_compare_drop_num <- function(drop_num_est_in){
  png("compare_drop_num.png")
  plot(drop_num_est_in$drop_num, drop_num_est_in$drop_num_est,main="dropout number of each gene",xlab="true",ylab="estimated")
  abline(0,1)
  dev.off()
}

#compare mean-variance relationship before dropout and after dropout
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

#compare mean-variance relationship before dropout and estimated
plot_mv<-function(mean_disp_in,titlename,filename){
  png(paste(filename,".png",sep="",collapse=""))
  title = paste(titlename," mean-dispersion",sep="",collapse="")
  plot(mean_disp_in$mean_beforedrop,mean_disp_in$disp_beforedrop,xlab="mean(log(cpm))",ylab="dispersion",main=title,col="blue",xlim=range(mean_disp_in$mean_beforedrop,mean_disp_in$mean_est ),ylim=range(mean_disp_in$disp_est,mean_disp_in$disp_beforedrop))
  points(mean_disp_in$mean_est,mean_disp_in$disp_est ,col="red")
  legend("topleft",1.4,legend=c("before dropout", "estimated"),col=c("blue","red"), lty=c(2,2), cex=0.8)
  dev.off()
}

#compare mean-variance relationship lowess curve before dropout and after dropout
plot_lowess_mv <- function(mean_disp_in,titlename,filename,span){
  png(paste(filename,".png",sep="",collapse=""))
  title = paste(titlename," mean-dispersion",sep="",collapse="")
  plot(lowess(mean_disp_in$mean_beforedrop,mean_disp_in$disp_beforedrop,f=span),main=title,xlab="mean(log(cpm))",ylab="dispersion",col="blue",xlim=range(mean_disp_in$mean_beforedrop,mean_disp_in$mean_est ),ylim=range(mean_disp_in$disp_est,mean_disp_in$disp_beforedrop))
  lines(lowess(mean_disp_in$mean_est,mean_disp_in$disp_est,f=span),col="red")
  legend("topleft",1.4,legend=c("before dropout", "estimated"),col=c("blue","red"), lty=c(2,2), cex=0.8)
  dev.off()
}









