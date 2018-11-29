## performance of voom and edgeR evalated by area under percision recall curve (0-0.8)
root = ""

name="klein"

libSizes_time_test = c(0.1,0.2,0.5,  1,2,5,  10,20)

samples_test = c(70,100,500,1000)#1000000

range=c(1:4)
samples_test = c(70,100,500,1000, 5000,10000)

# calculate area with 10 iterations
root = "/Users/yiwang/Dropbox/project_Kasper/8_17/out/"
setwd("/Users/yiwang/Dropbox/project_Kasper/8_17/plots/")
area_voom =array(dim=c(length(samples_test),length(libSizes_time_test),10)) 
area_edgeR =array(dim=c(length(samples_test),length(libSizes_time_test),10)) 

for( iter in 1:10{ 
 for(i in range){
  for(j in 1:length(libSizes_time_test)){ # vary library size
    nSample2=samples_test[i]  #very sample size
    libSizes_time =libSizes_time_test[j]
    for( iter in 1:10){
      name2=paste0(name,"_",nSample2,"_",libSizes_time,"_iter",iter)
      path_file= paste0(root,name2,"_cobraplot.RData")
      load(path_file)
      cobraplot2=cobraplot
      cobraplot2@fdrtprcurve $TPR = (1-cobraplot@fdrtprcurve $FDR) #preision - TDR in y lab
      cobraplot2@fdrtprcurve $FDR = cobraplot@fdrtprcurve $TPR #recall - TPR in x lab
      range_voom =[which(cobraplot2@fdrtprcurve$method=="ZINB-WaVE_limmavoom_common")]
      x_limmavoom = cobraplot2@fdrtprcurve$FDR[range_voom]
      y_limmavoom = cobraplot2@fdrtprcurve$TPR[range_voom]
      range_edgeR =[which(cobraplot2@fdrtprcurve$method=="ZINB-WaVE_edgeR_common")]
      x_edgeR = cobraplot2@fdrtprcurve$FDR[range_edgeR]
      y_edgeR = cobraplot2@fdrtprcurve$TPR[range_edgeR]
      ## remove the 1st point
      x_edgeR = x_edgeR[2:length(x_edgeR)]
      y_edgeR = y_edgeR[2:length(y_edgeR)]
      x_limmavoom=x_limmavoom[2:length(x_limmavoom)]
      y_limmavoom = y_limmavoom[2:length(y_limmavoom)]
      ## find the first common x-axis value of two curves
      x_start_loc = min(which(x_edgeR %in% x_limmavoom))
      x_edgeR_2 = x_edgeR [x_start_loc: length(x_edgeR)]
      y_edgeR_2 = y_edgeR[x_start_loc: length(x_edgeR)]
      x_start_value = x_edgeR[x_start_loc]
      x_limmavoom_2 = x_limmavoom [min(which(x_limmavoom== x_start_value)) : length(x_limmavoom)]
      y_limmavoom_2 = y_limmavoom [min(which(x_limmavoom== x_start_value)) : length(x_limmavoom)]
      # area calculation
      area_voom[i,j,iter] =  area_cal_pr2(x_limmavoom_2,y_limmavoom_2) 
      area_edgeR[i,j,iter] =  area_cal_pr2(x_edgeR_2,y_edgeR_2)
      if(x_start_value<0.8){ 
        area_voom[i,j,iter] = dif_area_l[i,j,iter] +x_start_value
        area_edgeR[i,j,iter] = dif_area_e[i,j,iter] +x_start_value

      }
      if(!x_start_value<0.8){ 
        area_voom[i,j,iter] = 0.8
        area_edgeR[i,j,iter] = 0.8

      }
    }
  }
}
