source("simulationHelpFunctions_v7_diffInZero.R")
source("function_NBsimSingleCell.R") 

load("params_Klein.RData") #params0: lambda, dispersion and lib.size of Klein dataset

#set additional parameters for simulation
nSamples=70 #number of samples
nTags=400 #number of genes/features
grp=as.factor(rep(0:1, each = nSamples/2)) #two-group comparison
DEind = sample(1:nTags,floor(nTags*.25),replace=FALSE) #25% DE genes
fcSim=(2 + rexp(length(DEind), rate = 1/2)) #fold changes of DE genes

#simulate scRNA data
simData <- NBsimSingleCell(foldDiff = fcSim, ind = DEind,
                                  nTags = nTags,
                                   group = grp,
                                   verbose = TRUE, params = params0,
                                    normalizeLambda=TRUE)

#simData$counts0: counts without dropout
#simData$counts: counts with dropout