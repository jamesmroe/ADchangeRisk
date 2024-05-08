#========================================================================================#
# Author: James M Roe, Ph.D.
# Contributors: Didac Vidal-Piñeiro, Ph.D.
# Center for Lifespan Changes in Brain and Cognition, University of Oslo
#
# Purpose: Cross validation grid search for ML
#          Script requires individual-level data as input and is not executable
#========================================================================================#

#--------inputs
args = commandArgs(TRUE)
eta=as.numeric(args[1])
max_depth=as.numeric(args[2])
gamma=as.numeric(args[3])
min_child_weight=as.numeric(args[4])
nrounds=as.numeric(args[5])
i=as.numeric(args[6])
idx=as.numeric(args[7])
odir=as.character(args[8])
parcellation=as.character(args[9])
INPUT=as.character(args[10]) #rSLOPE/FULLEFF
nfold=10
print(args)
#--------inputs


HyperParamSearch = function(eta, max_depth, gamma, min_child_weight, nrounds, i, idx, odir, parcellation) {

  #---load packages
  tmp.packages = c("dplyr", "stringr", "magrittr","xgboost")
  # tmpnew.packages = tmp.packages[!(tmp.packages %in% installed.packages()[,"Package"])]
  # if(length(tmpnew.packages)) {
  #   install.packages(tmpnew.packages)
  # }
  sapply(tmp.packages, require, character.only = T)
  rm(list=ls(pattern="tmp*"))
  
  
  #---load data
  a="/ess/p23/cluster/workspace/projects/11_MemP/James/ADchangeRisk"
  FULLDATA=file.path(a,"data/DF.Rda")
  print("loading DF.Rda")
  load(FULLDATA)
  SLOPEDATA=file.path(a,"ADNI/results/results.prepSlopes.Rda")
  print("loading results.prepSlopes.Rda")
  load(SLOPEDATA)

  
  #---subset cortical ROIs based on atlas selection
  print(paste("keeping", parcellation, "cortical ROIs"))
  if (parcellation == "aparc") {
    idx1=grepl("schaf200",ROISFINAL)
    idx2=grepl("destr",ROISFINAL)
    idxrm=idx1+idx2
  } else if (parcellation == "dest") {
    idx1=grepl("schaf200",ROISFINAL)
    idx2=grepl("aparc",ROISFINAL)
    idxrm=idx1+idx2 
  } else if (parcellation == "schaf200") {
    idx1=grepl("destr",ROISFINAL)
    idx2=grepl("aparc",ROISFINAL)
    idxrm=idx1+idx2
  } else {
    print(paste("Error. Invalid parcellation provided (schaf200 / aparc / dest). Quitting.")); quit()
  }
  
  
  if (INPUT == "rSLOPE") {
    INPUTDAT = rSLOPE
  } else if (INPUT == "rINTONLY") {
    INPUTDAT = rINTONLY
  } else if (INPUT == "FULLEFF") {
    INPUTDAT = FULLEFF
  }
  
  print(paste("input data is", INPUT))
  print("nfeatures in full dataset")
  print(dim(INPUTDAT))
  
  INPUTDATX = INPUTDAT[,!idxrm]
  ROIS = ROISFINAL[!idxrm]
  print("nfeatures in subsetted dataset")
  print(dim(INPUTDATX))

  
  #---cross validate
  train = INPUTDATX %>% as.matrix()
  #####
  # train = scale(train)
  #####
  labels = df.U %>% mutate(label=ifelse(Status == "AD", 1, 0)) %>% select(label) %>% as.matrix()
  params = list(booster = "gbtree",
              objective = "binary:logistic",
              eta = eta,
              max_depth=max_depth,
              gamma = gamma,
              min_child_weight = min_child_weight)
              
  
  #---define output
  aucTrain = aucTest = aucTestN = c()
  aucTestMin = aucTestNMin = c()
  
  
  #scale variables??
  xgbcv <- xgb.cv( params = params,
                   data = train,
                   label = labels,
                   nrounds = nrounds,
                   nfold = nfold,
                   showsd = T,
                   stratified = T,
                   print_every_n = 10,
                   early_stop_round = 10,
                   maximize = F,
                   metrics=list("logloss","auc","error"),
                   prediction = F)
  
  
  aucTest = c(aucTest, max(xgbcv$evaluation_log$test_auc_mean))
  aucTestN = c(aucTestN, which.max(xgbcv$evaluation_log$test_auc_mean))
  aucTrain = c(aucTrain, max(xgbcv$evaluation_log$train_auc_mean))
  
  
  print("loading searchGrid")
  ofile=file.path(odir, "searchGrid")
  load(paste0(ofile,".Rda"))
  
    
  print("including output in searchGrid")
  searchGrid$aucTest[idx]=aucTest
  searchGrid$aucTestN[idx]=aucTestN
  searchGrid$aucTrain[idx]=aucTrain
  searchGrid$i[idx]=i
  
  print("saving file")
  
  print(searchGrid[idx,])
  # save("searchGrid", 
  #      file = paste0(ofile,".Rda"))
  
  
}


HyperParamSearch(eta, max_depth, gamma, min_child_weight, nrounds, i, idx, odir, parcellation)
