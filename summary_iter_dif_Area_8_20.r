root = "/Users/yiwang/Dropbox/project_Kasper/8_17/out/"
setwd("/Users/yiwang/Dropbox/project_Kasper/8_17/plots/")

name="klein"

libSizes_time_test = c(0.1,0.2,0.5,  1,2,5,  10,20)

samples_test = c(70,100,500,1000)#1000000

range=c(1:4)
#samples_test = c(70,100,500,1000, 5000,10000)

#cal area with iteration
root = "/Users/yiwang/Dropbox/project_Kasper/8_17/out/"
setwd("/Users/yiwang/Dropbox/project_Kasper/8_17/plots/")
dif_area =array(dim=c(length(samples_test),length(libSizes_time_test),10)) #edgeR - limmavoom
dif_area_l =array(dim=c(length(samples_test),length(libSizes_time_test),10)) #edgeR - limmavoom
dif_area_e =array(dim=c(length(samples_test),length(libSizes_time_test),10)) #edgeR - limmavoom

# with iter
for(i in range){
  for(j in 1:length(libSizes_time_test)){
        nSample2=samples_test[i]
      libSizes_time =libSizes_time_test[j]
      name3=paste0(name,"_",nSample2,"_",libSizes_time)
      #x_limmavoom=c()
      #y_limmavoom=c()
      #x_edgeR=c()
      #y_edgeR=c()
        for( iter in 1:10){
        name2=paste0(name,"_",nSample2,"_",libSizes_time,"_iter",iter)
        path_file= paste0(root,name2,"_cobraplot.RData")
        load(path_file)
        cobraplot2=cobraplot
      cobraplot2@fdrtprcurve $TPR = (1-cobraplot@fdrtprcurve $FDR) #preision - TDR in y lab
      cobraplot2@fdrtprcurve $FDR = cobraplot@fdrtprcurve $TPR #recall - TPR in x lab
      x_limmavoom = cobraplot2@fdrtprcurve$FDR[which(cobraplot2@fdrtprcurve$method=="ZINB-WaVE_limmavoom2_common")]
      y_limmavoom = cobraplot2@fdrtprcurve$TPR[which(cobraplot2@fdrtprcurve$method=="ZINB-WaVE_limmavoom2_common")]
      x_edgeR = cobraplot2@fdrtprcurve$FDR[which(cobraplot2@fdrtprcurve$method=="ZINB-WaVE_edgeR_common")]
      y_edgeR = cobraplot2@fdrtprcurve$TPR[which(cobraplot2@fdrtprcurve$method=="ZINB-WaVE_edgeR_common")]
      #x_limmavoom=as.matrix(rbind((x_limmavoom),as.matrix(cobraplot2@fdrtprcurve$FDR[which(cobraplot2@fdrtprcurve$method=="limmavoom")] )))
      #y_limmavoom=as.matrix(rbind((y_limmavoom),as.matrix(cobraplot2@fdrtprcurve$TPR[which(cobraplot2@fdrtprcurve$method=="limmavoom")] )))
      #x_edgeR=as.matrix(rbind((x_edgeR),as.matrix(cobraplot2@fdrtprcurve$FDR[which(cobraplot2@fdrtprcurve$method=="edgeR")] )))
      #y_edgeR=as.matrix(rbind((y_edgeR),as.matrix(cobraplot2@fdrtprcurve$TPR[which(cobraplot2@fdrtprcurve$method=="edgeR")] )))
      x_edgeR = x_edgeR[2:length(x_edgeR)]
      y_edgeR = y_edgeR[2:length(y_edgeR)]
      x_limmavoom=x_limmavoom[2:length(x_limmavoom)]
      y_limmavoom = y_limmavoom[2:length(y_limmavoom)]
      x_start_loc = min(which(x_edgeR %in% x_limmavoom))
      x_edgeR_2 = x_edgeR [x_start_loc: length(x_edgeR)]
      y_edgeR_2 = y_edgeR[x_start_loc: length(x_edgeR)]
      x_start_value = x_edgeR[x_start_loc]
      x_limmavoom_2 = x_limmavoom [min(which(x_limmavoom== x_start_value)) : length(x_limmavoom)]
      y_limmavoom_2 = y_limmavoom [min(which(x_limmavoom== x_start_value)) : length(x_limmavoom)]
      dif_area[i,j,iter] =  area_cal_pr2(x_edgeR_2,y_edgeR_2) - area_cal_pr2(x_limmavoom_2,y_limmavoom_2) 
      dif_area_l[i,j,iter] =  area_cal_pr2(x_limmavoom_2,y_limmavoom_2) 
      dif_area_e[i,j,iter] =  area_cal_pr2(x_edgeR_2,y_edgeR_2)
      if(x_start_value<0.8){ 
        dif_area_l[i,j,iter] = dif_area_l[i,j,iter] +x_start_value
        dif_area_e[i,j,iter] = dif_area_e[i,j,iter] +x_start_value

      }
      if(!x_start_value<0.8){ 
        dif_area_l[i,j,iter] = 0.8
        dif_area_e[i,j,iter] = 0.8

      }
  }
  }}

dif_area_iter_2 =array(dim=c(10,length(libSizes_time_test),length(samples_test))) #edgeR - limmavoom


for(i in 1:length(libSizes_time_test)){
  dif_area_iter_2[,i,]=t(dif_area[,i,])

}




dif_area_after_re = dif_area
dif_area_after_re_2  = dif_area_iter_2




dif_area_after_re_2_e =array(dim=c(10,length(libSizes_time_test),length(samples_test))) #edgeR - limmavoom


for(i in 1:length(libSizes_time_test)){
  dif_area_after_re_2_e[,i,]=t(dif_area_e[,i,])

}

dif_area_after_re_2_l =array(dim=c(10,length(libSizes_time_test),length(samples_test))) #edgeR - limmavoom


for(i in 1:length(libSizes_time_test)){
  dif_area_after_re_2_l[,i,]=t(dif_area_l[,i,])

}

mean(dif_area_after_re_2_l[,5,4])
mean(dif_area_after_re_2_e[,5,4])




setwd("/Users/yiwang/Dropbox/project_Kasper/8_17/plots/")


samples_test = c(70,100,500,1000)#1000000

for(nsample in range){
  mat = dif_area_after_re_2[,,nsample]
  lib.size=c("0.1","0.2","0.5","1","2","5","10","20")
  lib.size_data = rep(lib.size,10)
  value_data=array(dim=c(80))
  iter_data=array(dim=c(80))

  count=0
  for(iter in 1:10){
     for( j in 1:length(lib.size)){
      count = count+1
      value_data[count] = mat[iter,j]
     iter_data[count]=iter
   }
  }

  data = data.frame(lib.size_data,value_data,iter_data)
  png(paste0("summary",samples_test[nsample],"_samples.png"))

  print(
    ggerrorplot(data, x = "lib.size_data", y = "value_data", 
         desc_stat = "mean_sd", color = "black",#ylim=c(-0.2,0.3),
            add = "jitter", add.params = list(color = "darkgray"),order=c("0.1","0.2","0.5","1","2","5","10","20"),
            title = "", xlab="* lib.size",ylab="area(edgeR)- area(voom)")
  + font("xlab", size = 30)+
 font("ylab", size = 25)+
 font("xy.text", size = 20)
  ) 

  dev.off()

}




setwd("/Users/yiwang/Dropbox/project_Kasper/8_17/plots/")


samples_test = c(70,100,500,1000)#1000000

for(nsample in range){
  mat_e = dif_area_after_re_2_e[,,nsample]
  lib.size=c("0.1","0.2","0.5","1","2","5","10","20")
  lib.size_data = rep(lib.size,10)
  value_data=array(dim=c(80))
  iter_data=array(dim=c(80))

  count=0
  for(iter in 1:10){
     for( j in 1:length(lib.size)){
      count = count+1
      value_data[count] = mat_e[iter,j]
     iter_data[count]=iter
   }
  }
  method=rep("edgeR",length(iter_data))
  data_e = data.frame(lib.size_data,value_data,iter_data,method)


   mat_l = dif_area_after_re_2_l[,,nsample]
  lib.size=c("0.1","0.2","0.5","1","2","5","10","20")
  lib.size_data = rep(lib.size,10)
  value_data=array(dim=c(80))
  iter_data=array(dim=c(80))

  count=0
  for(iter in 1:10){
     for( j in 1:length(lib.size)){
      count = count+1
      value_data[count] = mat_l[iter,j]
     iter_data[count]=iter
   }
  }
  method=rep("voom",length(iter_data))

  data_l = data.frame(lib.size_data,value_data,iter_data,method)

  data=rbind(data_l,data_e)
  png(paste0("compare",samples_test[nsample],"_samples.png"))

  print(
    ggerrorplot(data, x = "lib.size_data", y = "value_data", color="method",palette=c("blue", "red"),ylim=c(0,0.8),
         desc_stat = "mean_sd", 
            add = "jitter", add.params = list(size = 0.4),order=c("0.1","0.2","0.5","1","2","5","10","20"),
            title = "", xlab="* lib.size",ylab="area")
  + font("xlab", size = 30)+
 font("ylab", size = 25)+
 font("xy.text", size = 20)
  ) 

  dev.off()

}



setwd("/Users/yiwang/Dropbox/project_Kasper/8_17/plots/")



libSizes_time_test = c(0.1,0.2,0.5,  1,2,5,  10,20)


for(nlib in 1:length(libSizes_time_test)){
  mat = dif_area_after_re_2[,nlib,1:4]
  samples_test = c("70","100","500","1000")#10000,1000, 5000,10000)
  samples_test_data = rep(samples_test,10)
  value_data=array(dim=c(40))
  iter_data=array(dim=c(40))

  count=0
  for(iter in 1:10){
     for( j in 1:length(samples_test)){
      count = count+1
      value_data[count] = mat[iter,j]
      iter_data[count]=iter
   }
  }

  data = data.frame(samples_test_data,value_data,iter_data)
  png(paste0("summary",libSizes_time_test[nlib],"_libs.png"))


  print(
    ggerrorplot(data, x = "samples_test_data", y = "value_data", 
         desc_stat = "mean_sd", color = "black",ylim=c(-0.2,0.3),
            add = "jitter", add.params = list(color = "darkgray"),order=c("70","100","500","1000"),
            title = "", xlab="sample size",ylab="area(edgeR)- area(voom)")
  + font("xlab", size = 30)+
 font("ylab", size = 25)+
 font("xy.text", size = 20)
  ) 

  dev.off()

}




#dedgeR




libSizes_time_test = c(0.1,0.2,0.5,  1,2,5,  10,20)


for(nlib in 1:length(libSizes_time_test)){
  mat = dif_area_after_re_2_e[,nlib,1:4]
  samples_test = c("70","100","500","1000")#10000,1000, 5000,10000)
  samples_test_data = rep(samples_test,10)
  value_data=array(dim=c(40))
  iter_data=array(dim=c(40))
  count=0
  for(iter in 1:10){
     for( j in 1:length(samples_test)){
      count = count+1
      value_data[count] = mat[iter,j]
      iter_data[count]=iter
   }
  }
  method=rep("edgeR",length(value_data))
  data_e = data.frame(samples_test_data,value_data,iter_data,method)


  mat = dif_area_after_re_2_l[,nlib,1:4]
  samples_test = c("70","100","500","1000")#10000,1000, 5000,10000)
  samples_test_data = rep(samples_test,10)
  value_data=array(dim=c(40))
  iter_data=array(dim=c(40))
  count=0
  for(iter in 1:10){
     for( j in 1:length(samples_test)){
      count = count+1
      value_data[count] = mat[iter,j]
      iter_data[count]=iter
   }
  }
  method=rep("voom",length(value_data))
  data_l = data.frame(samples_test_data,value_data,iter_data,method)

  data=rbind(data_l,data_e)

  png(paste0("compare",libSizes_time_test[nlib],"_libs.png"))
  print(
    ggerrorplot(data, x = "samples_test_data", y = "value_data", palette=c("blue", "red"),ylim=c(0,0.8),
         desc_stat = "mean_sd", color = "method",#ylim=c(-0.05,0.2),
            add = "jitter",add.params = list(size = 0.4),order=c("70","100","500","1000"),
            title = "", xlab="sample size",ylab="area")
  + font("xlab", size = 30)+
 font("ylab", size = 25)+
 font("xy.text", size = 20)
  ) 

  dev.off()

}




















#replace 2nd group:
      samples_test=5000
      libSizes_time_test =0.5
      iter=8,9,10

setwd("/Users/yiwang/Dropbox/project_Kasper/8_17/plots/")
for(i in 1:length(samples_test)){
  for(j in 1:length(libSizes_time_test)){
        nSample2=samples_test[i]
      libSizes_time =libSizes_time_test[j]
      libSizes_time2=libSizes_time
      if(libSizes_time2==0.1){libSizes_time2="0_1"}
      if(libSizes_time2==0.2){libSizes_time2="0_2"}
      if(libSizes_time2==0.5){libSizes_time2="0_5"}
      for( iter in c(1:10)){
        name2=paste0(name,"_",nSample2,"_",libSizes_time,"_iter",iter,"_v2")
        path_file= paste0(root,name2,"_cobraplot.RData")
        load(path_file)
        name3=paste0(name,"_",nSample2,"_",libSizes_time2,"_iter",iter)
        plots_3way(cobraplot,name3)
    }
  }
}



#iter 9 still error
#So replace original iter=8,9,10 with new iter=3,1,2





#replace 3rd group:
      samples_test=10000
      libSizes_time_test =1
      iter=3,6

setwd("/Users/yiwang/Dropbox/project_Kasper/8_17/plots/")
for(i in 1:length(samples_test)){
  for(j in 1:length(libSizes_time_test)){
        nSample2=samples_test[i]
      libSizes_time =libSizes_time_test[j]
      libSizes_time2=libSizes_time
      if(libSizes_time2==0.1){libSizes_time2="0_1"}
      if(libSizes_time2==0.2){libSizes_time2="0_2"}
      if(libSizes_time2==0.5){libSizes_time2="0_5"}
      for( iter in c(1:10)){
        name2=paste0(name,"_",nSample2,"_",libSizes_time,"_iter",iter,"_v2")
        path_file= paste0(root,name2,"_cobraplot.RData")
        load(path_file)
        name3=paste0(name,"_",nSample2,"_",libSizes_time2,"_iter",iter)
        plots_3way(cobraplot,name3)
    }
  }
}



#So replace original iter=3,6 with new iter=1,2





#replace 4rd group:
      samples_test=10000
      libSizes_time_test =0.5
      iter=6,7,8,9,10

setwd("/Users/yiwang/Dropbox/project_Kasper/8_17/plots/")
for(i in 1:length(samples_test)){
  for(j in 1:length(libSizes_time_test)){
        nSample2=samples_test[i]
      libSizes_time =libSizes_time_test[j]
      libSizes_time2=libSizes_time
      if(libSizes_time2==0.1){libSizes_time2="0_1"}
      if(libSizes_time2==0.2){libSizes_time2="0_2"}
      if(libSizes_time2==0.5){libSizes_time2="0_5"}
      for( iter in c(1:10)){
        name2=paste0(name,"_",nSample2,"_",libSizes_time,"_iter",iter,"_v2")
        path_file= paste0(root,name2,"_cobraplot.RData")
        load(path_file)
        name3=paste0(name,"_",nSample2,"_",libSizes_time2,"_iter",iter)
        plots_3way(cobraplot,name3)
    }
  }
}

#iter 4 still error

#So replace original iter=6,7,8,9,10 with new iter=1,6,2,3,5





#replace 5th group:
      samples_test=10000
      libSizes_time_test =0.2
      iter=8

setwd("/Users/yiwang/Dropbox/project_Kasper/8_17/plots/")
for(i in 1:length(samples_test)){
  for(j in 1:length(libSizes_time_test)){
        nSample2=samples_test[i]
      libSizes_time =libSizes_time_test[j]
      libSizes_time2=libSizes_time
      if(libSizes_time2==0.1){libSizes_time2="0_1"}
      if(libSizes_time2==0.2){libSizes_time2="0_2"}
      if(libSizes_time2==0.5){libSizes_time2="0_5"}
      for( iter in c(1:10)){
        name2=paste0(name,"_",nSample2,"_",libSizes_time,"_iter",iter,"_v2")
        path_file= paste0(root,name2,"_cobraplot.RData")
        load(path_file)
        name3=paste0(name,"_",nSample2,"_",libSizes_time2,"_iter",iter)
        plots_3way(cobraplot,name3)
    }
  }
}

#iter 4 still error

#So replace original iter=8 with new iter=1















