norm_counts_3<-function(count_inp){
  library_size = colSums(count_inp)
  y <- t(log( t( (count_inp )/(library_size + 1) * 1e+06  ) + 0.5 ))
  return(y )
}


get_est_dropout_num<-function(sim_counts){
  out=list()
  ngene = length(sim_counts$dropoutrates3)
  dropoutrates3 = sim_counts$dropoutrates3
  dropoutrates3[is.na(dropoutrates3)]=0
  drop_num  = array(dim=ngene)
  for(i in 1: ngene){
    drop_num [i] = length(which(sim_counts$Dropout[i,]=="TRUE" & sim_counts$real_counts[i,]!=0) ) # true dropout numbers
  }
  dropoutrates3 = sim_counts$dropoutrates3
  dropoutrates3[is.na(dropoutrates3)]=0
  drop_num_est  = array(dim=ngene)
  for(i in 1: ngene){
    if(dropoutrates3[i]>0 & dropoutrates3[i]<1)
      {drop_num_est [i] = length(which(sim_counts$drop_counts[i,]!=0 )) / ((1/dropoutrates3[i]) - 1)}   # estimated dropout numbers
    else if (dropoutrates3[i]==0){ drop_num_est [i] =0 }
    else { drop_num_est [i] =length(which(sim_counts$drop_counts[i,]!=0 )) }
  }
  drop_num_est=round(drop_num_est,0)
  out$drop_num_est = drop_num_est
  out$drop_num = drop_num
  return(out)
}



mv_knowndrop_norm4<-function(sim_counts,sim_counts_n_real,sim_counts_n_drop){
  Dropout  = sim_counts$Dropout
  droprate_in =  sim_counts$dropoutrates3
  droprate_in[is.na(droprate_in)]=0
  ngenes = dim(sim_counts_n_drop)[1]
  ncells=dim(sim_counts_n_drop)[2]
  mean=array(dim=ngenes)
  disp=mean
  out=list()
  for( i in 1:ngenes){
    count_ = sim_counts_n_drop[i,]
    count_ = count_[which(count_!= log(0.5))] # remove all zeros in after-drop counts
    real_zero_num = length(which(sim_counts_n_real[i,]==log(0.5))) # real_zero_num
    real_zero_num = round(real_zero_num * (1-droprate_in[i]), 0)
    count_ = c(count_, rep(0,real_zero_num) )
    mean[i] = mean(count_)
    disp[i]=(var(count_))^(1/4)
  }
  mean[is.nan(mean)]=0
  disp[is.nan(disp)]=0
  disp[is.na(disp)]=0
  mean[is.na(mean)]=0
  mean[which(mean == "Inf")]=0
  disp[which(disp == "Inf")]=0
  out$mean_est =mean
  out$disp_est=disp
  out$mean_beforedrop = rowMeans(sim_counts_n_real)
  out$disp_beforedrop = (rowVars(sim_counts_n_real))^(1/4)
  return(out)
}



mv_unknowndrop_est_3 <-function(sim_counts,sim_counts_n_real,sim_counts_n_drop,dropnum_est){ # remove (1-droprate) zeros from non-drop zeros
  droprate_in = sim_counts$dropoutrates3
  droprate_in[is.na(droprate_in)]=0
  ngenes = dim(sim_counts_n_drop)[1]
  ncells=dim(sim_counts_n_drop)[2]
  mean=array(dim=ngenes)
  disp=mean
  out=list()
  sim_counts_n_drop[is.na(sim_counts_n_drop)]=log(0.5)
  for( i in 1:ngenes){
    count_ = sim_counts_n_drop[i,]
    real_zero_num = length( which(sim_counts$drop_counts[i,] == 0))-dropnum_est[i] #estimated real_zero_num
    real_zero_num = round(real_zero_num * (1-droprate_in[i]), 0)
    if(real_zero_num>0){
      count_ = c(count_[which(count_!= log(0.5))], rep(log(0.5), real_zero_num ))}
    else{count_ = count_[which(count_!= log(0.5))]}
    mean[i] = mean(count_)
    disp[i]=(var(count_))^(1/4)
  }
  mean[is.nan(mean)]=0
  disp[is.nan(disp)]=0
  disp[is.na(disp)]=0
  mean[is.na(mean)]=0
  mean[which(mean == "Inf")]=0
  disp[which(disp == "Inf")]=0
  out$mean_est =mean
  out$disp_est=disp
  out$mean_beforedrop = rowMeans(sim_counts_n_real)
  out$disp_beforedrop = (rowVars(sim_counts_n_real))^(1/4)
  return(out)
}

mv_unknowndrop_est_2 <-function(sim_counts,sim_counts_n_real,sim_counts_n_drop,dropnum_est){ #keep all non-drop zeros
  droprate_in = sim_counts$dropoutrates3 
  droprate_in[is.na(droprate_in)]=0
  ngenes = dim(sim_counts_n_drop)[1]
  ncells=dim(sim_counts_n_drop)[2]
  mean=array(dim=ngenes)
  disp=mean
  out=list()
  sim_counts_n_drop[is.na(sim_counts_n_drop)]=log(0.5)
  for( i in 1:ngenes){
    count_ = sim_counts_n_drop[i,]
    real_zero_num = length( which(sim_counts$drop_counts[i,] == 0))-dropnum_est[i] #estimated real_zero_num
    if(real_zero_num>0){
      count_ = c(count_[which(count_!= log(0.5))], rep(log(0.5), real_zero_num ))}
    else{count_ = count_[which(count_!= log(0.5))]}
    mean[i] = mean(count_)
    disp[i]=(var(count_))^(1/4)
  }
  mean[is.nan(mean)]=0
  disp[is.nan(disp)]=0
  disp[is.na(disp)]=0
  mean[is.na(mean)]=0
  mean[which(mean == "Inf")]=0
  disp[which(disp == "Inf")]=0
  out$mean_est =mean
  out$disp_est=disp
  out$mean_beforedrop = rowMeans(sim_counts_n_real)
  out$disp_beforedrop = (rowVars(sim_counts_n_real))^(1/4)
  return(out)
}



