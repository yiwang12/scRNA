#source("/Users/yiwang/Dropbox/project_Kasper/8_6/util_methods.R")
source("/users/ywang/data/8_7/util_methods_zinb.R")

NCORES <- 20
registerDoParallel(NCORES)
register(DoparParam())







name="klein" #error
#source("/Users/yiwang/Dropbox/YiWang_project/code/3_23/function_load_datasets_splatter.R") #loadDataset function from splatter
source("/users/ywang/data/function_load_datasets_splatter.R") #loadDataset function from splatter



library(readr) #read_tsv
#root <- "/Users/yiwang/Dropbox/project_Kasper/DE_and_sc_tools/splatter-paper-master/data"
root <- "/users/ywang/data/"


datasets <- read_tsv(file.path(root, "datasets.txt"),
                     col_types = cols(.default = col_character(),
                                      NumCells = col_integer()
                     )
)
load_spla<-function(num){
  set.seed(1)
  dataset <- datasets[num, ]
  name <- unlist(dataset["Dataset"])
  message("***LOADING ", name, "***")
  counts <- loadDataset(dataset, root)
  na.rows <- which(rowSums(is.na(counts)) > 0)
  if (length(na.rows) > 0) {
    counts <- counts[-na.rows, ]
  }
  counts <- counts[, sample(1:ncol(counts), 200)]
  counts <- counts[rowSums(counts) > 0, ]
}
counts_klein = load_spla(2)

counts_realdata=(counts_klein[1:400,1:70])




params0 = getDatasetMoMPositive(counts = counts_realdata)

name="klein"
samples_test = c(70,100,500,1000)

libSizes_time_test = c(0.1)

test_lib_sample_effect_drop<-function(counts_realdata_in,params0_in, nSamples_in, libSizes_time_in, name){
libSizes1=sample(colSums(counts_realdata_in),nSamples_in,replace=TRUE) #library sizes
libSizes2=libSizes_time_in * libSizes1
grp=as.factor(rep(0:1, each = nSamples_in/2)) #two-group comparison
nTags=400 #nr of genes
set.seed(11)
DEind = sample(1:nTags,floor(nTags*.25),replace=FALSE) #25% DE
fcSim=(2 + rexp(length(DEind), rate = 1/2)) #fold changes
simDataTrapnell2 <- NBsimSingleCell(foldDiff = fcSim, ind = DEind,
                                   dataset = counts_realdata_in, nTags = nTags,
                                   group = grp,
                                   verbose = TRUE, params = params0_in,
                                   lib.size = libSizes2, normalizeLambda=TRUE)

cobraplot=prepare_plot(simDataTrapnell2,grp,nTags)
save(cobraplot,file=paste0(name,"_cobraplot.RData"))
#plots_3way(cobraplot,name)
}


for(i in 1:length(samples_test)){
  for(j in 1:length(libSizes_time_test)){
    for( iter in 1){
      set.seed(iter)
      nSample2 = samples_test[i]
      libSizes_time = libSizes_time_test[j]
      name2=paste0(name,"_",nSample2,"_",libSizes_time,"_iter",iter)
      test_lib_sample_effect_drop(counts_realdata, params0,nSample2, libSizes_time , name2)
     }
  }
}


