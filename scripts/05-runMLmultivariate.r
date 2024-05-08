#================================================================================================================================#
#
# Author: James M Roe, PhD.
# Center for Lifespan Changes in Brain and Cognition, University of Oslo

# Purpose: run ML analysis in AD-control data (ADNI)
#          apply AD-control weights to healthy adult lifespan brain change estimates conditional on age
#          and test PRS-AD-change associations in healthy adults (LCBC as test data)
#          multivariate analysis using features with accelerated change to test PRS-AD-brain change assocations in healthy adults
#
#          Script requires individual-level data as input and is not executable
#================================================================================================================================#


rm(list=ls())
options(bitmapType = "cairo")
#setwd("/ess/p23/cluster/workspace/projects/11_MemP/James/ADchangeRisk")


#---load packages----
tmp.packages = c("dplyr", "stringr", "magrittr","xgboost","caret","MetBrewer","tidymodels",
                 "here","sgof","cowplot","vip", "asbio",
                 "data.table")
tmpnew.packages = tmp.packages[!(tmp.packages %in% installed.packages()[,"Package"])]
if(length(tmpnew.packages)) {
  #install.packages(tmpnew.packages)
}
sapply(tmp.packages, require, character.only = T)
rm(list=ls())


here()
params1 = list(nrounds=500, eta=0.2, max_depth=5, gamma=1, min_child_weight=1)
params3 = list(nrounds=600, eta=0.01, max_depth=7, gamma=0.5, min_child_weight=2)
parcellation = "aparc"


#---load ADNI Slopes----
a = "/ess/p23/cluster/workspace/projects/11_MemP/James/ADchangeRisk"
FULLDATA = file.path(a, "data/DF.Rda")
print("loading DF.Rda")
load(FULLDATA)
SLOPEDATA = file.path(a, "ADNI/results/results.prepSlopes.Rda")
print("loading results.prepSlopes.Rda")
load(SLOPEDATA)


#---analysis steps-----
runML = 1           #step1
applyMLweights = 1  #step2
combineMLres = 1    #step3
plotMLweights = 1   #step4
runPCA1 = 1         #step5
mergePCA1 = 1       #step6
runPCAroll = 1      #step7
mergePCAroll = 1    #step8
mergeALL = 1        #step9
plotALL = 1         #plots



# ---- STEP1 runML -----

if (runML == 1)  {
  
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
    print(paste("Error. Invalid parcellation provided (schaf200 / aparc / destr). Quitting.")); quit()
  }
  
  print("nfeatures in full dataset")
  print(dim(rSLOPE))
  rSLOPEX = rSLOPE[,!idxrm]
  ABSCHANGEX = FULLEFF[,!idxrm]
  ROIS = ROISFINAL[!idxrm]
  print("nfeatures in subsetted dataset")
  print(dim(rSLOPEX))
  print(dim(ABSCHANGEX))
  # writeLines(ROIS,here("xgbModels/features.txt"))
  
  
  ## MODEL TRAINING-----
  # rSLOPE (age-relative change)
  train = rSLOPEX %>% as.matrix()
  labels = df.U %>% mutate(label=ifelse(Status == "AD", 1, 0)) %>% select(label) %>% as.matrix()
  xgb_train1 = xgb.DMatrix(data = train, label = labels)
  
  params = list(booster = "gbtree",
                objective = "binary:logistic",
                eta = params1$eta,
                max_depth = params1$max_depth,
                gamma = params1$gamma,
                min_child_weight = params1$min_child_weight)
              
  
  bstSlope <- xgboost( params = params,
                   data = xgb_train1,
                   nrounds = params1$nrounds,
                   early_stopping_rounds = 20,
                   metrics=list("logloss","auc","error"),
                   prediction = F)
  
  
  importance_matrix1 = xgb.importance(ROIS, model = bstSlope)
  
  
  
  # ABSCHANGE (absolute change)
  train = ABSCHANGEX %>% as.matrix()
  xgb_train3 = xgb.DMatrix(data = train, label = labels)
  
  params = list(booster = "gbtree",
                objective = "binary:logistic",
                eta = params3$eta,
                max_depth = params3$max_depth,
                gamma = params3$gamma,
                min_child_weight = params3$min_child_weight)
  
  
  bstABSCHANGE <- xgboost( params = params,
                         data = xgb_train3,
                         nrounds = params3$nrounds,
                         early_stopping_rounds = 20,
                         metrics=list("logloss","auc","error"),
                         prediction = F)
  
  
  ## POST TRAINING----
  adnirSLOPEX = rSLOPEX
  adniABSCHANGEX = ABSCHANGEX
  
  
  ## importance matrices----
  # rSlope
  xgb.plot.importance(importance_matrix1[1:54,])
  importance_matrix1 = xgb.importance(ROIS, model = bstSlope)
  
  # absChange
  importance_matrix3 = xgb.importance(ROIS, model = bstABSCHANGE)
  xgb.plot.importance(importance_matrix3[1:10,])
  
}



# ---- STEP2 applyMLweights -----

if (applyMLweights == 1) {
  
  # load LCBC test data----
  load(file.path(a, "LCBC/results/results.prepSlopes30-publish.Rda"))
  if (parcellation == "aparc") {
    idx1=grepl("schaf200",ROISFINAL.LCBC)
    idx2=grepl("destr",ROISFINAL.LCBC)
    idxrm=idx1+idx2
  }
  testDat1 = as.data.frame(rSLOPE.LCBC[,!idxrm])
  testDat3 = as.data.frame(FULLEFF.LCBC[,!idxrm])
  names(testDat1) = names(testDat3) = ROISFINAL.LCBC[!idxrm]
  names(testDat1) == ROIS
  names(testDat1) %in% ROIS
  
  
  # VIP features must be ordered same between train and test data
  testDat1 %<>% select(ROIS)
  testDat3 %<>% select(ROIS)
  
  
  # predict using LCBC data
  xgb.pred1.met1.prob = predict(bstSlope, as.matrix(testDat1), type="prob")
  xgb.pred3.met1.prob = predict(bstABSCHANGE, as.matrix(testDat3), type="prob")
  
   
  # log odds ratio of model weights prediction
  U.LCBC %<>% mutate(xgb.pred1.met1.prob, xgb.pred3.met1.prob,
    pred1class=ifelse(xgb.pred1.met1.prob > .5, "AD", "NC"),
    pred3class=ifelse(xgb.pred3.met1.prob > .5, "AD", "NC"),
    p1 = log( xgb.pred1.met1.prob / (1-xgb.pred1.met1.prob) ),
    p3 = log( xgb.pred3.met1.prob / (1-xgb.pred3.met1.prob) ))
  table(U.LCBC$pred1class)
  table(U.LCBC$pred3class)

  
  ggplot(U.LCBC, aes(x=meanAge, y=xgb.pred1.met1.prob))  + geom_point(aes(col=factor(pred1class))) + geom_smooth(method="lm")
  ptmp1 = ggplot(U.LCBC, aes(x=meanAge, y=p1))  + geom_point(aes(col=factor(pred1class))) + geom_smooth(method="lm", alpha=0.2) + geom_hline(yintercept = 0, linetype=2, size=2) + theme_classic()
  ptmp3 = ggplot(U.LCBC, aes(x=meanAge, y=p3))  + geom_point(aes(col=factor(pred3class))) + geom_smooth(alpha=0.2) + geom_hline(yintercept = 0, linetype=2, size=2) + theme_classic()
  
  

  #---load PGS data ----
  loadPGS = function(siglevel) {
    pgs_path = file.path(a, "LCBC/results/new_PGS_for_LCBC")
    PGS = data.frame(
                   #APOE inclusive
                   fread(paste0(pgs_path, "/AD_Jansen/wAPOE/AD_Jansen.",siglevel,".profile")) %>% select(FID,SCORE) %>% mutate("PGSS_AD_Jansen" = SCORE),
                   fread(paste0(pgs_path, "/AD_Wightman/wAPOE/AD_Wightman.",siglevel,".profile")) %>% select(FID,SCORE) %>% mutate("PGSS_AD_Wightman" = SCORE),
                   fread(paste0(pgs_path, "/AD/wAPOE/AD.",siglevel,".profile")) %>% select(FID,SCORE) %>% mutate("PGSS_AD" = SCORE),
                   fread(paste0(pgs_path, "/AD_kunkle/wAPOE/AD_Kunkle.",siglevel,".profile")) %>% select(FID,SCORE) %>% mutate("PGSS_AD_Kunkle" = SCORE),
                   #APOE exclusive
                   fread(paste0(pgs_path, "/AD_Jansen/woAPOE/AD_Jansen_woAPOE.",siglevel,".profile")) %>% select(FID,SCORE) %>% mutate("PGSS_AD_Jansen_woAPOE" = SCORE),
                   fread(paste0(pgs_path, "/AD_Wightman/woAPOE/AD_Wightman_woAPOE.",siglevel,".profile")) %>% select(FID,SCORE) %>% mutate("PGSS_AD_Wightman_woAPOE" = SCORE),
                   fread(paste0(pgs_path, "/AD/woAPOE/AD.",siglevel,".profile")) %>% select(FID,SCORE) %>% mutate("PGSS_AD_woAPOE" = SCORE),
                   fread(paste0(pgs_path, "/AD_kunkle/woAPOE/AD_Kunkle_woAPOE.",siglevel,".profile")) %>% select(FID,SCORE) %>% mutate("PGSS_AD_Kunkle_woAPOE" = SCORE)
    )
    PGS$FID.1=NULL
    PGS$FID.2=NULL
    PGS$FID.3=NULL
    PGS$FID.4=NULL
    names(PGS)[1]="Genetic_ID"
    PGS
  }
  PGS = loadPGS("S1")
  PGSS = merge(U.LCBC, PGS, by="Genetic_ID")
  
  
  #---link genetic pcs----
  pcs = fread(file.path(a, "LIFEBRAIN_LCBC_PCA.eigenvec"), header = F) %>% 
    setnames(., colnames(.), c("Genetic_ID", "IID", paste0("PC", 1:20))) %>% select(1,3:12)
  
  
  # N229 full LCBC genetic sample (only QC'd European ancestry in pcs file)
  PGSS %<>% select(-starts_with("GAF"))
  names(pcs)[2:11] = paste0("GAF_PC", 1:10)
  PGSS %<>% merge(., pcs)
  
  
  #---full genetic data DF----
  PGSSFULL = PGSS
  
  
  #--- ML weights applied to LCBC change estimates----
  agecuts = c(30, 35, 40, 45, 50, 55, 60, 65, 70)
  for (i in 1:length(agecuts)) {
    if (i == 1) { 
        OUTLIST_ML=list()
      }
    
    #--- raise lower age ---
    agecut2 = agecuts[i]
    PGSS = PGSSFULL %>% filter(meanAge >= agecut2)
    #--- raise lower age ---
    
    
    runLM <- function(PGSS, Yvar, Xvar) {
      
      covariates = c(
        "scale(meanAge)",
        "subject_sex",
        "mri_info_site_name",
        "GAF_PC1",
        "GAF_PC2",
        "GAF_PC3",
        "GAF_PC4",
        "GAF_PC5",
        "GAF_PC6",
        "GAF_PC7",
        "GAF_PC8",
        "GAF_PC9",
        "GAF_PC10",
        "nTimepoints",
        "intervalFirstlast")
      lm_formula <- as.formula(paste(Yvar, "~", Xvar, "+", paste(covariates, collapse = " + ")))
      
      model <- lm(lm_formula, data = PGSS)
      modelsum = summary(model)
      
      return(
        list(
          model = model,
          summary = modelsum
        )
      )
    }
    
    
    getPartialR2 = function(model, score, ap=1) {
      
      if (ap == 1) {
        if (score == "Jansen") { model_drop = update(model, as.formula(paste0(". ~ . -scale(PGSS_AD_Jansen)"))) }
        if (score == "AD") { model_drop = update(model, as.formula(paste0(". ~ . -scale(PGSS_AD)"))) }
        if (score == "Kunkle") { model_drop = update(model, as.formula(paste0(". ~ . -scale(PGSS_AD_Kunkle)"))) }
        if (score == "Wightman") { model_drop = update(model, as.formula(paste0(". ~ . -scale(PGSS_AD_Wightman)"))) }
      } else if (ap == 0) {
        if (score == "Jansen") { model_drop = update(model, as.formula(paste0(". ~ . -scale(PGSS_AD_Jansen_woAPOE)"))) }
        if (score == "AD") { model_drop = update(model, as.formula(paste0(". ~ . -scale(PGSS_AD_woAPOE)"))) }
        if (score == "Kunkle") { model_drop = update(model, as.formula(paste0(". ~ . -scale(PGSS_AD_Kunkle_woAPOE)"))) }
        if (score == "Wightman") { model_drop = update(model, as.formula(paste0(". ~ . -scale(PGSS_AD_Wightman_woAPOE)"))) }
      }
      # ( SSE(reduced) - SSE(full) )
      # ___________________________
      #
      #       SSE(reduced)
      
      partialR=asbio::partial.R2(model_drop, model)
      
      
      #calc CI
      cor_ciiR <- cor.test(residuals(model), residuals(model_drop))$conf.int
      partialR_CIlwr = sort(1 - cor_ciiR^2)[1]
      partialR_CIupr = sort(1 - cor_ciiR^2)[2]
      
      return(list(
        data.frame(
          prsq=partialR, 
          prsq_lwr=partialR_CIlwr, 
          prsq_upr=partialR_CIupr))
      )
    }
    
    
    # PRS-AD relchange weights models -----
    model1list_ml = runLM(PGSS, "scale(p1)", "scale(PGSS_AD_Jansen)")
    model1_ml = model1list_ml$model
    m1_ml = model1list_ml$summary
    
    model2list_ml = runLM(PGSS, "scale(p1)", "scale(PGSS_AD)")
    model2_ml = model2list_ml$model
    m2_ml = model2list_ml$summary
    
    model3list_ml = runLM(PGSS, "scale(p1)", "scale(PGSS_AD_Kunkle)")
    model3_ml = model3list_ml$model
    m3_ml = model3list_ml$summary
    
    model4list_ml = runLM(PGSS, "scale(p1)", "scale(PGSS_AD_Wightman)")
    model4_ml = model4list_ml$model
    m4_ml = model4list_ml$summary
    
    
    # PRS-ADnoAPOE relchange weights models -----
    model1noAPlist_ml = runLM(PGSS, "scale(p1)", "scale(PGSS_AD_Jansen_woAPOE)")
    model1noAP_ml = model1noAPlist_ml$model
    m1noAP_ml = model1noAPlist_ml$summary
    
    model2noAPlist_ml = runLM(PGSS, "scale(p1)", "scale(PGSS_AD_woAPOE)")
    model2noAP_ml = model2noAPlist_ml$model
    m2noAP_ml = model2noAPlist_ml$summary
    
    model3noAPlist_ml = runLM(PGSS, "scale(p1)", "scale(PGSS_AD_Kunkle_woAPOE)")
    model3noAP_ml = model3noAPlist_ml$model
    m3noAP_ml = model3noAPlist_ml$summary
    
    model4noAPlist_ml = runLM(PGSS, "scale(p1)", "scale(PGSS_AD_Wightman_woAPOE)")
    model4noAP_ml = model4noAPlist_ml$model
    m4noAP_ml = model4noAPlist_ml$summary
    
    
    # PRS-AD absolute change weights models ----
    model5list_ml = runLM(PGSS, "scale(p3)", "scale(PGSS_AD_Jansen)")
    model5_ml = model5list_ml$model
    m5_ml = model5list_ml$summary
    
    model6list_ml = runLM(PGSS, "scale(p3)", "scale(PGSS_AD)")
    model6_ml = model6list_ml$model
    m6_ml = model6list_ml$summary
    
    model7list_ml = runLM(PGSS, "scale(p3)", "scale(PGSS_AD_Kunkle)")
    model7_ml = model7list_ml$model
    m7_ml = model7list_ml$summary
    
    model8list_ml = runLM(PGSS, "scale(p3)", "scale(PGSS_AD_Wightman)")
    model8_ml = model8list_ml$model
    m8_ml = model8list_ml$summary
    
    
    # PRS-ADnoAPOE absolute change weights models ----
    model5noAPlist_ml = runLM(PGSS, "scale(p3)", "scale(PGSS_AD_Jansen_woAPOE)")
    model5noAP_ml = model5noAPlist_ml$model
    m5noAP_ml = model5noAPlist_ml$summary
    
    model6noAPlist_ml = runLM(PGSS, "scale(p3)", "scale(PGSS_AD_woAPOE)")
    model6noAP_ml = model6noAPlist_ml$model
    m6noAP_ml = model6noAPlist_ml$summary
    
    model7noAPlist_ml = runLM(PGSS, "scale(p3)", "scale(PGSS_AD_Kunkle_woAPOE)")
    model7noAP_ml = model7noAPlist_ml$model
    m7noAP_ml = model7noAPlist_ml$summary
    
    model8noAPlist_ml = runLM(PGSS, "scale(p3)", "scale(PGSS_AD_Wightman_woAPOE)")
    model8noAP_ml = model8noAPlist_ml$model
    m8noAP_ml = model8noAPlist_ml$summary
    
    
    # --- combine results
    OUTLIST_ML[[i]] = rbind(
      
          # PRS-AD relchange results
          tidy(m1_ml, conf.int = T) %>% mutate(N=nrow(PGSS),
                                            score="Jansen", type="rSlope", 
                                            ap=1, agecut2=agecut2,
                                            getPartialR2(model1_ml, "Jansen")[[1]]) %>% filter(grepl("PGS", term)),
          tidy(m2_ml, conf.int = T) %>% mutate(N=nrow(PGSS),score="AD", type="rSlope", 
                                            ap=1, agecut2=agecut2,
                                            getPartialR2(model2_ml, "AD")[[1]]) %>% filter(grepl("PGS", term)),
          tidy(m3_ml, conf.int = T) %>% mutate(N=nrow(PGSS),score="Kunkle", type="rSlope", 
                                            ap=1, agecut2=agecut2,
                                            getPartialR2(model3_ml, "Kunkle")[[1]]) %>% filter(grepl("PGS", term)),
          tidy(m4_ml, conf.int = T) %>% mutate(N=nrow(PGSS),score="Wightman", type="rSlope", 
                                            ap=1, agecut2=agecut2,
                                            getPartialR2(model4_ml, "Wightman")[[1]]) %>% filter(grepl("PGS", term)),
          
          
          # PRS-ADnoAPOE relchange results
          tidy(m1noAP_ml, conf.int = T) %>% mutate(N=nrow(PGSS),score="Jansen", type="rSlope", 
                                                ap=0, agecut2=agecut2,
                                                getPartialR2(model1noAP_ml, "Jansen", ap=0)[[1]]) %>% filter(grepl("PGS", term)),
          tidy(m2noAP_ml, conf.int = T) %>% mutate(N=nrow(PGSS),score="AD", type="rSlope", 
                                                ap=0, agecut2=agecut2,
                                                getPartialR2(model2noAP_ml, "AD", ap=0)[[1]]) %>% filter(grepl("PGS", term)),
          tidy(m3noAP_ml, conf.int = T) %>% mutate(N=nrow(PGSS),score="Kunkle", type="rSlope", 
                                                ap=0, agecut2=agecut2,
                                                getPartialR2(model3noAP_ml, "Kunkle", ap=0)[[1]]) %>% filter(grepl("PGS", term)),
          tidy(m4noAP_ml, conf.int = T) %>% mutate(N=nrow(PGSS),score="Wightman", type="rSlope", 
                                                ap=0, agecut2=agecut2,
                                                getPartialR2(model4noAP_ml, "Wightman", ap=0)[[1]]) %>% filter(grepl("PGS", term)),
          
          
          # PRS-ADnoAPOE abschange results
          tidy(m5noAP_ml, conf.int = T) %>% mutate(N=nrow(PGSS),score="Jansen", type="absChange", 
                                            ap=0, agecut2=agecut2,
                                                     getPartialR2(model5noAP_ml, "Jansen", ap=0)[[1]]) %>% filter(grepl("PGS", term)),
          tidy(m6noAP_ml, conf.int = T) %>% mutate(N=nrow(PGSS),score="AD", type="absChange", 
                                            ap=0, agecut2=agecut2,
                                                     getPartialR2(model6noAP_ml, "AD", ap=0)[[1]]) %>% filter(grepl("PGS", term)),
          tidy(m7noAP_ml, conf.int = T) %>% mutate(N=nrow(PGSS),score="Kunkle", type="absChange", 
                                            ap=0, agecut2=agecut2,
                                                     getPartialR2(model7noAP_ml, "Kunkle", ap=0)[[1]]) %>% filter(grepl("PGS", term)),
          tidy(m8noAP_ml, conf.int = T) %>% mutate(N=nrow(PGSS),score="Wightman", type="absChange", 
                                            ap=0, agecut2=agecut2,
                                                     getPartialR2(model8noAP_ml, "Wightman", ap=0)[[1]]) %>% filter(grepl("PGS", term)),
          
          
          # PRS-AD abschange results
          tidy(m5_ml, conf.int = T) %>% mutate(N=nrow(PGSS),score="Jansen", type="absChange", 
                                            ap=1, agecut2=agecut2,
                                               getPartialR2(model5_ml, "Jansen")[[1]]) %>% filter(grepl("PGS", term)),
          tidy(m6_ml, conf.int = T) %>% mutate(N=nrow(PGSS),score="AD", type="absChange", 
                                            ap=1, agecut2=agecut2,
                                               getPartialR2(model6_ml, "AD")[[1]]) %>% filter(grepl("PGS", term)),
          tidy(m7_ml, conf.int = T) %>% mutate(N=nrow(PGSS),score="Kunkle", type="absChange", 
                                            ap=1, agecut2=agecut2,
                                               getPartialR2(model7_ml, "Kunkle")[[1]]) %>% filter(grepl("PGS", term)),
          tidy(m8_ml, conf.int = T) %>% mutate(N=nrow(PGSS),score="Wightman", type="absChange", 
                                            ap=1, agecut2=agecut2,
                                               getPartialR2(model8_ml, "Wightman")[[1]]) %>% filter(grepl("PGS", term)))
    
  }
}



# ---- STEP3 combineMLres + FDR ----

if (combineMLres == 1) {
  OUT_ML = bind_rows(OUTLIST_ML) %>% mutate(dotalpha=ifelse(p.value < .05, 1, 0))

  OUT_ML %<>% mutate(agerange=paste0(agecut2,"-89"),
                  agerange=factor(agerange, levels=c("30-89", "35-89", "40-89", "45-89", "50-89", "55-89", "60-89", "65-89", "70-89")),
                  PGS=paste0("AD_", score))
  
  OUT_ML$PGS[OUT_ML$PGS=="AD_AD"]="AD_Lambert"
  
    
  ## FDR MLweights ----
  pvalcor0 =
  c(
    OUT_ML$p.value[OUT_ML$ap == 1 & OUT_ML$type == "rSlope"],
    OUT_ML$p.value[OUT_ML$ap == 1 & OUT_ML$type == "absChange"]
    )
  FDRpML=sgof::BH(pvalcor0, alpha=0.05)
  (FDRthreshpML=max(FDRpML$data[FDRpML$Adjusted.pvalues<.05]))
}
  


# ---- STEP4 plotMLweights ----

if (plotMLweights == 1) {
  for (iii in 1:2) {
    
    if (iii==1) {
      pvalcor0=c()
      pML=pMLrsq=pMLBoth=list()
      
      #apoe inclusive results (for FDR correction)
      OUT1_ML = OUT_ML %>% filter(ap==1, type=="rSlope")
      typeSelect = "rSlope"; selection="apoe"; #just to track plot titles
      
      #apoe exclusive results
      OUT2_ML = OUT_ML %>% filter(ap==0, type=="rSlope")
      
      
    } else if (iii==2) {
     
      #apoe inclusive results (for FDR correction)
      OUT1_ML = OUT_ML %>% filter(ap==1, type=="absChange")
      typeSelect = "absChange"; selection="apoe"; #just to track plot titles
      
      #apoe exclusive results
      OUT2_ML = OUT_ML %>% filter(ap==0, type=="absChange")
      
    }
    
    mytheme = theme(
      plot.background = element_rect(fill = "white"),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      title = element_text(size=17),
      text = element_text(color = "black", size = 18, family="Nimbus Sans Narrow"),
      plot.title = element_text(hjust = 0.5),
      axis.ticks = element_blank(),
      axis.title.y = element_text(color = "black", size = 22, vjust =-1, margin = margin(0,20,0,0)),
      axis.title.x = element_text(color = "black", size = 22, vjust = -2, margin = margin(0,20,20,0)),
      axis.text = element_text(color = "black", size = 18),
      legend.key.size = unit(1,"cm"))
    
    
    if (iii != 2) {
    
      (pML[[iii]] = ggplot(OUT1_ML, aes(x=agerange, y=-log10(p.value), group=factor(PGS), shape=factor(PGS))) +
          scale_x_discrete() +
          scale_alpha_discrete(range = c(0,1)) +
          geom_hline(yintercept = -log10(0.05), linetype=2) +
          geom_hline(yintercept = -log10(FDRthreshpML), linetype=3) +
          geom_point(position = position_dodge(width=0.5),colour="grey") +
          geom_point(aes(alpha=factor(dotalpha)), col="black", position=position_dodge(width=0.5), shape=21, size=4, stroke=0.75) +
          theme_classic() + 
          labs(x = "Age-range") +
          ggtitle(paste0("Model weights ",typeSelect, selection)) +
          mytheme +
          theme(axis.text = element_text(size=8),
                axis.title.x = element_text(size=7, vjust=-3, margin=margin(0,30,0,0)),
                axis.title.y = element_text(size=7, margin = margin(0,10,0,0)),
                legend.text = element_text(size=8),
                legend.title = element_text(size=8),
                plot.title = element_text(size=10),
                plot.subtitle = element_text(size=10, hjust=0.5)
          ) +
         #change shape
         scale_shape_manual(values = c("AD_Jansen"=16, "AD_Kunkle"=15, "AD_Lambert"=17, "AD_Wightman"=4))
       )
      
    } else {
      
      (pML[[iii]] = ggplot(OUT1_ML, aes(x=agerange, y=-log10(p.value), group=factor(PGS), shape=factor(PGS))) +
         scale_x_discrete() +
         scale_alpha_discrete(range = c(0,1)) +
         geom_hline(yintercept = -log10(0.05), linetype=2) +
         geom_hline(yintercept = -log10(FDRthreshpML), linetype=3) +
         geom_point(position = position_dodge(width=0.5),colour="black") +
         geom_point(aes(alpha=factor(dotalpha)), col="black", position=position_dodge(width=0.5), shape=21, size=4, stroke=0.75) +
         theme_classic() + 
         labs(x = "Age-range") +
         # scale_colour_manual(values=met.brewer("Derain", 4)) +
         ggtitle(paste0("Model weights ",typeSelect, selection)) +
         mytheme +
         theme(axis.text = element_text(size=8),
               axis.title.x = element_text(size=7, vjust=-3, margin=margin(0,30,0,0)),
               axis.title.y = element_text(size=7, margin = margin(0,10,0,0)),
               legend.text = element_text(size=8),
               legend.title = element_text(size=8),
               plot.title = element_text(size=10),
               plot.subtitle = element_text(size=10, hjust=0.5)
         ) +
        #change shape
      scale_shape_manual(values = c("AD_Jansen"=16, "AD_Kunkle"=15, "AD_Lambert"=17, "AD_Wightman"=4))
      )
    
    }
    
    
    pML[[iii]] = pML[[iii]] +
      theme(strip.text.x = element_blank()) + 
      scale_x_discrete(labels=c('30+', '35+', '40+', '45+','50+','55+','60+','65+','70+'))
    
    
    # filename=filename
    # print(filename)
    # ggsave(
    #   plot = pML[[iii]],
    # file = filename,
    #   width = 12,
    #   height = 10,
    #   units = "cm",
    #   dpi = 600
    # )
  
  
    # add rsq to plot
    maxscalersq=.1844 #keep rsq on same scale as other analyses (i.e. max left hippo)
    
  
    #which were FDR sig w APOE
    OUT1_ML$FDRsig = factor(ifelse(OUT1_ML$p.value <= FDRthreshpML, 1, 0))
  
    #which FDR-corrected PRS-AD associations were significant wo/APOE
    which(OUT2_ML[OUT1_ML$FDRsig == 1,]$p.value < .05)
    OUT2_ML[OUT1_ML$FDRsig == 1 & OUT2_ML$p.value< .05,]$p.value
    OUT2_ML$dotalphapoe = factor(ifelse(OUT1_ML$FDRsig == 1 & OUT2_ML$p.value< .05, 1, 0))
  
    if (iii == 1) {
      
      (pMLrsq[[iii]] = OUT1_ML %>%
  
          ggplot(., aes(x=agerange, y=prsq*-1,group=factor(PGS)), col="grey", fill="grey") +#color=factor(PGS),fill=factor(PGS))) +
          scale_x_discrete() +
          scale_alpha_discrete(range = c(0, 1)) +
          ylim(c(maxscalersq*-1,0)) +
          geom_col(aes(alpha=as.factor(FDRsig), y=prsq*-1,x=agerange), color=NA, fill="grey", position=position_dodge(width=0.5), width=0.4) +
          geom_errorbar(aes(alpha=as.factor(FDRsig), ymin=prsq_lwr*-1,ymax=prsq_upr*-1),color="grey", fill="#945c6d",position=position_dodge(width=0.5),size=0.5, width=0) +
          geom_point(data=OUT2_ML, aes(alpha=dotalphapoe, y=prsq*-1, x=agerange),col="black",position=position_dodge(width=0.5),shape=3, size=3, stroke=1) +
          geom_hline(yintercept = 0, size=0.5) +
          labs(x="Age-range", y=paste(parse(text="R^2"),"PGS")) +
          theme_classic() +
          mytheme +
          theme(
            panel.background = element_rect(fill = "white"),
            axis.text = element_text(size = 8),
            axis.title.y = element_text(size = 8, margin = margin(0,10,0,0)),
            axis.title.x = element_blank(),
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 8),
            plot.title = element_text(size =10),
            plot.subtitle = element_text(size = 10, hjust=0.5),
            axis.ticks.x = element_blank(),
  
          ) +
          labs(color='PGS',fill='PGS', alpha="p<.05", shape=FALSE)
      )
      } else if (iii == 2) {
        
        (pMLrsq[[iii]] = OUT1_ML %>%
  
           ggplot(., aes(x=agerange, y=prsq*-1,group=factor(PGS)), col="grey", fill="grey") +#color=factor(PGS),fill=factor(PGS))) +
           scale_x_discrete() +
           scale_alpha_discrete(range = c(0, 1)) +
           ylim(c(maxscalersq*-1,0)) +
           geom_col(aes(alpha=as.factor(FDRsig), y=prsq*-1,x=agerange), color=NA, fill="dark grey", position=position_dodge(width=0.5), width=0.4) +
           geom_errorbar(aes(alpha=as.factor(FDRsig), ymin=prsq_lwr*-1,ymax=prsq_upr*-1),color="dark grey", fill="#945c6d",position=position_dodge(width=0.5),size=0.5, width=0) +
           geom_point(data=OUT2_ML, aes(alpha=dotalphapoe, y=prsq*-1, x=agerange),col="black",position=position_dodge(width=0.5),shape=3, size=3, stroke=1) +
           geom_hline(yintercept = 0, size=0.5) +
           labs(x="Age-range", y=paste(parse(text="R^2"),"PGS")) +
           theme_classic() +
           mytheme +
           theme(
             panel.background = element_rect(fill = "white"),
             axis.text = element_text(size = 8),
             axis.title.y = element_text(size = 8, margin = margin(0,10,0,0)),
             axis.title.x = element_blank(),
             legend.text = element_text(size = 8),
             legend.title = element_text(size = 8),
             plot.title = element_text(size =10),
             plot.subtitle = element_text(size = 10, hjust=0.5),
             axis.ticks.x = element_blank(),
  
           ) +
           labs(color='PGS',fill='PGS', alpha="p<.05", shape=FALSE)
        )
      }
  
  
      pMLrsq[[iii]] = pMLrsq[[iii]] +
        theme(strip.text.x = element_blank()) +
        scale_x_discrete(labels=c('30+', '35+', '40+', '45+','50+','55+','60+','65+','70+'))
  
      (pMLBoth[[iii]] = cowplot::plot_grid(pML[[iii]],
                                           pMLrsq[[iii]], nrow=2,rel_heights = c(2, 1),align= T))
  
    
    #   #newrsq
    #   # ggsave(
    #   #   plot = pMLBoth[[iii]],
    #   #   file = filename,
    #   #   width = 12,
    #   #   height = 10,
    #   #   units = "cm",
    #   #   dpi = 600
    #   # )
    # 
      
  }
}
pMLBoth[[1]]
pMLBoth[[2]]



# ---- STEP5 runPCA1: PCA1 FEATURE SELECTION ----

if (runPCA1 == 1) {
  agecuts=c(30,35,40,45,50,55,60,65,70)
  for (agecut in agecuts) {
    
    if (agecut == agecuts[1]) {
      j=3 #step
      w=50 #window size
      windowstart=1
      windowend=w
      
      ii=0
      OUTLIST_PC1=list()
      expl1=expl3=c()
      expl1noHC=expl3noHC=c()
    }
    
    ii=ii+1
    print(paste(ii, "agecut =", agecut))
    
    testDat1PCA = testDat1
    testDat3PCA = testDat3
    
    #agecut 30 data TMP
    U=U.LCBC
    testDat1PCA = testDat1PCA[U$meanAge > agecut,]
    testDat3PCA = testDat3PCA[U$meanAge > agecut,]
    U=U[U$meanAge > agecut,]
    
    
    #incuding HC and Amygdala
    pca1=prcomp(testDat1PCA[,importance_matrix1[windowstart:windowend,]$Feature],center=T,scale.=T)
    (spca1=summary(pca1))
    
    pca3=prcomp(testDat3PCA[,importance_matrix1[windowstart:windowend,]$Feature],center=T,scale.=T)
    (spca3=summary(pca3))
    
    expl1[ii] = (spca1$importance[2,1])
    U$spca1 = (spca1$x[,1])
    expl3[ii] = (spca3$importance[2,1])
    U$spca3 = (spca3$x[,1])
    
    
    #skip hippocampus and amygdala
    idxHippoVol=grepl(glob2rx("*vol*Hippocampus"),importance_matrix1[1:155,]$Feature,
                      ignore.case = T)
    idxAmygdVol=grepl(glob2rx("*vol*Amygdala*"),importance_matrix1[1:155,]$Feature,
                      ignore.case = T)
    idxfeatout=idxHippoVol+idxAmygdVol
    
    tmpout=importance_matrix1[1:155,]$Feature[!idxfeatout][windowstart:windowend]
    
    #works as expected
    # for (i in 1:10){
    #   print(importance_matrix1[1:155,]$Feature[!idxfeatout][windows[i,1]:windows[i,ncol(windows)]])
    # }
    
    pca1noHC=prcomp(testDat1PCA[,tmpout],center=T,scale.=T)
    (spca1noHC=summary(pca1noHC))
    
    
    pca3noHC=prcomp(testDat3PCA[,tmpout],center=T,scale.=T)
    (spca3noHC=summary(pca3noHC))
    
    expl1noHC[ii] = (spca1noHC$importance[2,1])
    expl3noHC[ii] = (spca3noHC$importance[2,1])
    
    U$spca1noHC = (spca1noHC$x[,1])
    U$spca3noHC = (spca3noHC$x[,1])
    
    #log exaplined variance
    (varexpl = data.frame(propvar=rbind(expl1, expl3, expl1noHC ,expl3noHC),hc=c(rep(1,2), rep(0,2))))
    
    ggplot(U, aes(x=meanAge, y=spca3noHC)) + geom_point() + geom_smooth() + geom_vline(xintercept = 50, linetype=2)
    
    
    ## absChange plot (no HC no AMY) ----
    if (agecut == 30) {
      (pAbs=ggplot(U, aes(x=meanAge, y=spca3noHC)) + 
       annotate("rect", xmin=50, xmax=max(U$meanAge), ymin=-Inf, ymax=+Inf, alpha=0.2, fill="grey") +
       geom_point(col="#7a334a") + 
       geom_smooth(col="black",method="loess") +
       geom_vline(xintercept = 50, linetype=2) + 
       geom_vline(xintercept = max(U$meanAge), linetype=2) + 
       theme_classic() +
       labs(y="AbsChange loading PC1", x="Age (mean)") +
       mytheme +
       theme(axis.text = element_text(size=8),
             axis.title.x = element_text(size=7, vjust=-3, margin=margin(0,30,0,0)),
             axis.title.y = element_text(size=7, margin = margin(0,10,0,0)),
             legend.text = element_text(size=8),
             legend.title = element_text(size=8),
             plot.title = element_text(size=10),
             plot.subtitle = element_text(size=10, hjust=0.5)
       ))
    }
    
    #reload PGS data and merge with subset
    SS = "S1"
    PGS = loadPGS(SS)
    
    #reload genetic pcs and merge with subset (ensures QCd + European)
    PGSS = merge(U, PGS,by="Genetic_ID")
    PGSS %<>% select(-starts_with("GAF"))
    names(pcs)[2:11] = paste0("GAF_PC", 1:10)
    PGSS %<>% merge(., pcs)
   
    
  
    ## relative change ---
    model1list_pc1 = runLM(PGSS, "scale(spca1noHC)", "scale(PGSS_AD_Jansen)")
    model1_pc1 = model1list_pc1$model
    (m1_pc1 = model1list_pc1$summary)
    
    model2list_pc1 = runLM(PGSS, "scale(spca1noHC)", "scale(PGSS_AD)")
    model2_pc1 = model2list_pc1$model
    (m2_pc1 = model2list_pc1$summary)
    
    model3list_pc1 = runLM(PGSS, "scale(spca1noHC)", "scale(PGSS_AD_Kunkle)")
    model3_pc1 = model3list_pc1$model
    (m3_pc1 = model3list_pc1$summary)
    
    model4list_pc1 = runLM(PGSS, "scale(spca1noHC)", "scale(PGSS_AD_Wightman)")
    model4_pc1 = model4list_pc1$model
    (m4_pc1 = model4list_pc1$summary)
    
  
    model1noAPlist_pc1 = runLM(PGSS, "scale(spca1noHC)", "scale(PGSS_AD_Jansen_woAPOE)")
    model1noAP_pc1 = model1noAPlist_pc1$model
    (m1noAP_pc1 = model1noAPlist_pc1$summary)
    
    model2noAPlist_pc1 = runLM(PGSS, "scale(spca1noHC)", "scale(PGSS_AD_woAPOE)")
    model2noAP_pc1 = model2noAPlist_pc1$model
    (m2noAP_pc1 = model2noAPlist_pc1$summary)
    
    model3noAPlist_pc1 = runLM(PGSS, "scale(spca1noHC)", "scale(PGSS_AD_Kunkle_woAPOE)")
    model3noAP_pc1 = model3noAPlist_pc1$model
    (m3noAP_pc1 = model3noAPlist_pc1$summary)
    
    model4noAPlist_pc1 = runLM(PGSS, "scale(spca1noHC)", "scale(PGSS_AD_Wightman_woAPOE)")
    model4noAP_pc1 = model4noAPlist_pc1$model
    (m4noAP_pc1 = model4noAPlist_pc1$summary)
    
    
    # --- combine PCA1 results
    OUTLIST_PC1[[ii]] = rbind(
      
      #PRS-AD relChange
      tidy(m1_pc1, conf.int = T) %>% mutate(start=windowstart,nfeat=w,windowend=windowend,
                                            score="Jansen", type="rSlope", 
                                            HC="noHC", w, SS=SS, ap=1, agecut=agecut,
                                        getPartialR2(model1_pc1, "Jansen")[[1]]) %>% filter(grepl("PGS", term)),
      tidy(m2_pc1, conf.int = T) %>% mutate(start=windowstart,nfeat=w,windowend=windowend,
                                            score="AD", type="rSlope", 
                                            HC="noHC", w, SS=SS, ap=1, agecut=agecut,
                                        getPartialR2(model2_pc1, "AD")[[1]]) %>% filter(grepl("PGS", term)),
      tidy(m3_pc1, conf.int = T) %>% mutate(start=windowstart,nfeat=w,windowend=windowend,
                                            score="Kunkle", type="rSlope", 
                                            HC="noHC", w, SS=SS, ap=1, agecut=agecut,
                                        getPartialR2(model3_pc1, "Kunkle")[[1]]) %>% filter(grepl("PGS", term)),
      tidy(m4_pc1, conf.int = T) %>% mutate(start=windowstart,nfeat=w,windowend=windowend,
                                            score="Wightman", type="rSlope", 
                                            HC="noHC", w, SS=SS, ap=1, agecut=agecut,
                                        getPartialR2(model4_pc1, "Wightman")[[1]]) %>% filter(grepl("PGS", term)),
      
      
      
      #PRS-ADnoAPOE relChange
      tidy(m1noAP_pc1, conf.int = T) %>% mutate(start=windowstart,nfeat=w,windowend=windowend,
                                                score="Jansen", type="rSlopeNoApoe", HC="noHC", w, SS=SS, ap=0, agecut=agecut,
                                          getPartialR2(model1noAP_pc1, "Jansen", ap=0)[[1]]) %>% filter(grepl("PGS", term)),
      tidy(m2noAP_pc1, conf.int = T) %>% mutate(start=windowstart,nfeat=w,windowend=windowend,
                                                score="AD", type="rSlopeNoApoe", HC="noHC", w, SS=SS, ap=0, agecut=agecut,
                                          getPartialR2(model2noAP_pc1, "AD", ap=0)[[1]]) %>% filter(grepl("PGS", term)),
      tidy(m3noAP_pc1, conf.int = T) %>% mutate(start=windowstart,nfeat=w,windowend=windowend,
                                                score="Kunkle", type="rSlopeNoApoe", HC="noHC", w, SS=SS, ap=0, agecut=agecut,
                                          getPartialR2(model3noAP_pc1, "Kunkle", ap=0)[[1]]) %>% filter(grepl("PGS", term)),
      tidy(m4noAP_pc1, conf.int = T) %>% mutate(start=windowstart,nfeat=w,windowend=windowend,
                                                score="Wightman", type="rSlopeNoApoe", HC="noHC", w, SS=SS, ap=0, agecut=agecut,
                                          getPartialR2(model4noAP_pc1, "Wightman", ap=0)[[1]]) %>% filter(grepl("PGS", term))
      
    )
  }
}



# ---- STEP6 mergePCA1: PCA1-50 FEATURE SELECTION ----

if (mergePCA1 == 1) {
  OUT_PC1 = bind_rows(OUTLIST_PC1) %>% mutate(dotalpha=ifelse(p.value < .05, 1, 0))
  OUT_PC1 %<>% mutate(agerange=paste0(agecut,"-89"),
                  agerange=factor(agerange, levels=c("30-89", "35-89", "40-89", "45-89", "50-89", "55-89", "60-89", "65-89", "70-89")),
                  PGS=paste0("AD_", score))
  OUT_PC1$PGS[OUT_PC1$PGS=="AD_AD"]="AD_Lambert"
  
  
  #PRS-AD results
  OUT1_PC1 = OUT_PC1 %>% filter(type == "rSlope", ap==1)
  
  #PRS-AD noAPOE results
  OUT2_PC1 = OUT_PC1 %>% filter(type == "rSlopeNoApoe", ap==0)
  
  
  #check out PRS-AD results
  ggplot(OUT1_PC1, aes(x=agerange, y=-log10(p.value), group=factor(PGS), shape=factor(PGS))) +
    scale_x_discrete() +
    scale_alpha_discrete(range = c(0,1)) +
    geom_hline(yintercept = -log10(0.05), linetype=2) +
    geom_point(position = position_dodge(width=0.5),colour="#793349") +
    geom_point(aes(alpha=factor(dotalpha)), col="black", position=position_dodge(width=0.5), shape=21, size=4, stroke=0.75) +
    theme_classic() + 
    labs(x = "Age-range") +
    ggtitle(paste0("PCA", w)) +
    mytheme + theme(legend.position = "none")
  
}



# ---- STEP7 runPCAroll: rolling PCA FEATURE SELECTION ----

if (runPCAroll == 1) {
  
  
  j=3 #step
  w=20 #window size
  
  
  #create feature windows for rolling PCA across first 100 features
  guess=100
  vec=c(1:guess)
  windows=matrix(NA, nrow = guess, ncol=w)
  for (c in 1:guess) {
    if (c==1) {
      minwin=1
    } else {
      minwin=(minwin+j)
    }
    maxwin=(minwin+w-1)
    windows[c,]=vec[minwin:maxwin]
  }
  windows %<>% na.omit()
  # print(windows)
  
  
  for (ww in 1:nrow(windows)) {
    
    if (ww == 1) {
      ii=0
      OUTLIST_PCroll=list()
      agecut=50
    }
    
    print(paste("window",ww))
    ii=ii+1
    OUTLIST_PCroll[[ii]] = list()
    
    testDat1PCA = testDat1
    U=U.LCBC
    testDat1PCA = testDat1PCA[U$meanAge > agecut,]
    U=U[U$meanAge > agecut,]
    
    windowstart=windows[ww,1]
    nfeats=windows[ww,w] #nfeats is row of window
    windowend = nfeats
    
  
    #skip hippocampus and amygdala
    idxHippoVol=grepl(glob2rx("*vol*Hippocampus"),importance_matrix1[1:155,]$Feature,
                      ignore.case = T)
    idxAmygdVol=grepl(glob2rx("*vol*Amygdala*"),importance_matrix1[1:155,]$Feature,
                      ignore.case = T)
    idxfeatout=idxHippoVol+idxAmygdVol
    
    tmpout=importance_matrix1[1:155,]$Feature[!idxfeatout][windowstart:windowend]
    
    #works as expected (check against feature bars in plot)
    # for (i in 1:10){
    #   print(importance_matrix1[1:155,]$Feature[!idxfeatout][windows[i,1]:windows[i,ncol(windows)]])
    # }
    
    pca1noHC=prcomp(testDat1PCA[,tmpout],center=T,scale.=T)
    (spca1noHC=summary(pca1noHC))
    
    U$spca1noHC = (spca1noHC$x[,1])
    
    #reload PGS data and merge with subset
    SS = "S1"
    PGS = loadPGS(SS)
    
    #reload genetic pcs and merge with subset (ensures QCd + European)
    PGSS = merge(U, PGS,by="Genetic_ID")
    PGSS %<>% select(-starts_with("GAF"))
    names(pcs)[2:11] = paste0("GAF_PC", 1:10)
    PGSS %<>% merge(., pcs)
    
    
    
    
    # PRS-AD relChange
    model1list_pc = runLM(PGSS, "scale(spca1noHC)", "scale(PGSS_AD_Jansen)")
    model1_pc = model1list_pc$model
    (m1_pc = model1list_pc$summary)
    
    model2list_pc = runLM(PGSS, "scale(spca1noHC)", "scale(PGSS_AD)")
    model2_pc = model2list_pc$model
    (m2_pc = model2list_pc$summary)
    
    model3list_pc = runLM(PGSS, "scale(spca1noHC)", "scale(PGSS_AD_Kunkle)")
    model3_pc = model3list_pc$model
    (m3_pc = model3list_pc$summary)
    
    model4list_pc = runLM(PGSS, "scale(spca1noHC)", "scale(PGSS_AD_Wightman)")
    model4_pc = model4list_pc$model
    (m4_pc = model4list_pc$summary)
    
    # PRS-ADnoAPOE relChange
    model1noAPlist_pc = runLM(PGSS, "scale(spca1noHC)", "scale(PGSS_AD_Jansen_woAPOE)")
    model1noAP_pc = model1noAPlist_pc$model
    (m1noAP_pc = model1noAPlist_pc$summary)
    
    model2noAPlist_pc = runLM(PGSS, "scale(spca1noHC)", "scale(PGSS_AD_woAPOE)")
    model2noAP_pc = model2noAPlist_pc$model
    (m2noAP_pc = model2noAPlist_pc$summary)
    
    model3noAPlist_pc = runLM(PGSS, "scale(spca1noHC)", "scale(PGSS_AD_Kunkle_woAPOE)")
    model3noAP_pc = model3noAPlist_pc$model
    (m3noAP_pc = model3noAPlist_pc$summary)
    
    model4noAPlist_pc = runLM(PGSS, "scale(spca1noHC)", "scale(PGSS_AD_Wightman_woAPOE)")
    model4noAP_pc = model4noAPlist_pc$model
    (m4noAP_pc = model4noAPlist_pc$summary)
    
    
    # --- combine PCAroll results
    OUTLIST_PCroll[[ii]] = rbind(
      
      # PRS-AD relChange
      tidy(m1_pc, conf.int = T) %>% mutate(start=windowstart,nfeat=w,windowend=windowend,score="Jansen", type="rSlope", 
                                            HC="noHC", ww, SS=SS, ap=1, agecut=agecut,
                                            getPartialR2(model1_pc, "Jansen")[[1]]) %>% filter(grepl("PGS", term)),
      tidy(m2_pc, conf.int = T) %>% mutate(start=windowstart,nfeat=w,windowend=windowend,score="AD", type="rSlope", 
                                            HC="noHC", ww, SS=SS, ap=1, agecut=agecut,
                                            getPartialR2(model2_pc, "AD")[[1]]) %>% filter(grepl("PGS", term)),
      tidy(m3_pc, conf.int = T) %>% mutate(start=windowstart,nfeat=w,windowend=windowend,score="Kunkle", type="rSlope", 
                                            HC="noHC", ww, SS=SS, ap=1, agecut=agecut,
                                            getPartialR2(model3_pc, "Kunkle")[[1]]) %>% filter(grepl("PGS", term)),
      tidy(m4_pc, conf.int = T) %>% mutate(start=windowstart,nfeat=w,windowend=windowend,score="Wightman", type="rSlope", 
                                            HC="noHC", ww, SS=SS, ap=1, agecut=agecut,
                                            getPartialR2(model4_pc, "Wightman")[[1]]) %>% filter(grepl("PGS", term)),
      
      
      
      # PRS-ADnoAPOE relChange
      tidy(m1noAP_pc, conf.int = T) %>% mutate(start=windowstart,nfeat=w,windowend=windowend,score="Jansen", type="rSlopeNoApoe", HC="noHC", ww, SS=SS, ap=0, agecut=agecut,
                                                getPartialR2(model1noAP_pc, "Jansen", ap=0)[[1]]) %>% filter(grepl("PGS", term)),
      tidy(m2noAP_pc, conf.int = T) %>% mutate(start=windowstart,nfeat=w,windowend=windowend,score="AD", type="rSlopeNoApoe", HC="noHC", ww, SS=SS, ap=0, agecut=agecut,
                                                getPartialR2(model2noAP_pc, "AD", ap=0)[[1]]) %>% filter(grepl("PGS", term)),
      tidy(m3noAP_pc, conf.int = T) %>% mutate(start=windowstart,nfeat=w,windowend=windowend,score="Kunkle", type="rSlopeNoApoe", HC="noHC", ww, SS=SS, ap=0, agecut=agecut,
                                                getPartialR2(model3noAP_pc, "Kunkle", ap=0)[[1]]) %>% filter(grepl("PGS", term)),
      tidy(m4noAP_pc, conf.int = T) %>% mutate(start=windowstart,nfeat=w,windowend=windowend,score="Wightman", type="rSlopeNoApoe", HC="noHC", ww, SS=SS, ap=0, agecut=agecut,
                                                getPartialR2(model4noAP_pc, "Wightman", ap=0)[[1]]) %>% filter(grepl("PGS", term))

    )
  }
}



# ---- STEP8 mergePCAroll results ----

if (mergePCAroll == 1) {
  
  OUT_PCroll = bind_rows(OUTLIST_PCroll) %>% mutate(dotalpha=ifelse(p.value < .05, 1, 0))
  OUT_PCroll %<>% mutate(agerange=paste0(agecut,"-89"),
                  agerange=factor(agerange, levels=c("30-89", "35-89", "40-89", "45-89", "50-89", "55-89", "60-89", "65-89", "70-89")),
                  PGS=paste0("AD_", score))
  OUT_PCroll$PGS[OUT_PCroll$PGS=="AD_AD"]="AD_Lambert"
  
  
  #PRS-AD results
  OUT1_PCroll = OUT_PCroll %>% filter(type == "rSlope", ap==1)
  
  #PRS-ADnoAPOE results
  OUT2_PCroll = OUT_PCroll %>% filter(type == "rSlopeNoApoe", ap==0)
  
  
  #check out results
  ggplot(OUT1_PCroll, aes(x=1:nrow(OUT1_PCroll), y=-log10(p.value), col=as.factor(ww), shape=as.factor(score))) +
      scale_alpha_discrete(range=c(0,1)) +
      geom_point() +
      geom_point(aes(alpha=as.factor(dotalpha)), col="black", shape=21, size=4,stroke=0.75) +
      geom_hline(yintercept = -log10(0.05), linetype=2) +
      labs(x="Window") +
      theme_classic() +
      mytheme + theme(legend.position = "none")

}



# ---- STEP9 mergeALL results ----

if (mergeALL == 1) {
  
  #FDR correct across all PRS-AD tests
  FDR_wAPOE = BH(c(OUT1_PC1$p.value, OUT1_PCroll$p.value), .05)
  (FDRthreshAPOE=max(FDR_wAPOE$data[FDR_wAPOE$Adjusted.pvalues<.05]))
  
  
  #which were FDR sig w APOE
  OUT1_PC1$FDRsig = factor(ifelse(OUT1_PC1$p.value < FDRthreshAPOE, 1, 0))
  
  
  #14 tests survived correction with PC1relChange 1-50 (no Hippo no Amy)
  sum(as.numeric(as.character(OUT1_PC1$FDRsig)))
  OUT1_PCroll$FDRsig = factor(ifelse(OUT1_PCroll$p.value <= FDRthreshAPOE, 1, 0))
  
  
  #13 tests survived correction with rolling PC1relChange (no Hippo no Amy)
  sum(as.numeric(as.character(OUT1_PCroll$FDRsig)))
  
  
  #one p.val FDR corrected and sig without APOE
  which(OUT2_PC1[OUT1_PC1$FDRsig == 1,]$p.value < .05)
  OUT2_PC1[OUT1_PC1$FDRsig == 1 & OUT2_PC1$p.value< .05,]$p.value
  OUT2_PC1$dotalphapoe = factor(ifelse(OUT1_PC1$FDRsig == 1 & OUT2_PC1$p.value< .05, 1, 0))
  
  
  #two p.vals FDR corrected and sig without APOE
  which(OUT2_PCroll[OUT1_PCroll$FDRsig == 1,]$p.value < .05)
  OUT2_PCroll[OUT1_PCroll$FDRsig == 1 & OUT2_PCroll$p.value< .05,]$p.value
  OUT2_PCroll$dotalphapoe = factor(ifelse(OUT1_PCroll$FDRsig == 1 & OUT2_PCroll$p.value< .05, 1, 0))

}


#save source files
# figure4B_PRSAD = OUT1_PC1
# figure4B_PRSADnoAPOE = OUT2_PC1
# figure4D_PRSAD = OUT1_PCroll
# figure4D_PRSADnoAPOE = OUT2_PCroll
# save('figure4B_PRSAD', 'figure4B_PRSADnoAPOE', 'figure4D_PRSAD', 'figure4D_PRSADnoAPOE',file = file.path(a, "reproduce/results/PRS-ADmodels_figure4_PCmultivariate.Rda"))



# ---- plotALL results ----

if (plotALL == 1) {
  
  
  # PC1 associations ----
  (pPCA = ggplot(OUT1_PC1, aes(x=agerange, y=-log10(p.value), group=factor(PGS), shape=factor(PGS))) +
     scale_x_discrete() +
     scale_alpha_discrete(range = c(0,1)) +
     geom_hline(yintercept = -log10(0.05), linetype=2) +
     geom_hline(yintercept = -log10(FDRthreshAPOE), linetype=3) +
     geom_point(position = position_dodge(width=0.5),colour="#793349") +
     geom_point(aes(alpha=factor(dotalpha)), col="black", position=position_dodge(width=0.5), shape=21, size=4, stroke=0.75) +
     theme_classic() + 
     labs(x = "Age-range") +
     ggtitle(paste0("PCA", w)) +
     mytheme +
     theme(axis.text = element_text(size=8),
           axis.title.x = element_text(size=7, vjust=-3, margin=margin(0,30,0,0)),
           axis.title.y = element_text(size=7, margin = margin(0,10,0,0)),
           legend.title = element_text(size=8),
           plot.title = element_text(size=10),
           plot.subtitle = element_text(size=10, hjust=0.5)
     ) +
     
   #change shape
   scale_shape_manual(values = c("AD_Jansen"=16, "AD_Kunkle"=15, "AD_Lambert"=17, "AD_Wightman"=4))
  )
  
  
  # PC1 estimates ----
  (pEst = ggplot(OUT1_PC1, aes(x=agerange, y=(abs(estimate)*-1), group=factor(PGS), shape=factor(PGS))) +
      scale_x_discrete() +
      scale_alpha_discrete(range = c(0,1)) +
      geom_point(position = position_dodge(width=0.5),colour="#793349",size=3.5) +
      theme_classic() + 
      labs(x = "Age-range") +
      ggtitle(paste0("PCA1")) +
      theme(axis.text = element_text(size=8),
            axis.title.x = element_text(size=7, vjust=-3, margin=margin(0,30,0,0)),
            axis.title.y = element_text(size=7, margin = margin(0,10,0,0)),
            legend.title = element_text(size=8),
            plot.title = element_text(size=10),
            plot.subtitle = element_text(size=10, hjust=0.5)
      ) +
      #change shape
      scale_shape_manual(values = c("AD_Jansen"=16, "AD_Kunkle"=15, "AD_Lambert"=17, "AD_Wightman"=4)) +
      theme(axis.text = element_text(size=14),
              axis.text.x = element_blank(),
              axis.ticks.x = element_line(),
              panel.border = element_rect(color = "black", fill = NA, size = 0.75))
  )
  
  
  # PC1 rsquared ----
  (pPCArsq = OUT1_PC1 %>%
      
      ggplot(., aes(x=agerange, y=prsq*-1,group=factor(PGS)), col="#945c6d", fill="#945c6d") +#color=factor(PGS),fill=factor(PGS))) +
      scale_x_discrete() +
      scale_alpha_discrete(range = c(0, 1)) +
      ylim(c(maxscalersq*-1,0)) +
      geom_col(aes(alpha=as.factor(FDRsig), y=prsq*-1,x=agerange), color=NA, fill="#945c6d", position=position_dodge(width=0.5), width=0.4) +
      geom_errorbar(aes(alpha=as.factor(FDRsig), ymin=prsq_lwr*-1,ymax=prsq_upr*-1),color="#945c6d",position=position_dodge(width=0.5),size=0.5, width=0) +
      geom_point(data=OUT2_PC1, aes(alpha=dotalphapoe, y=prsq*-1, x=agerange),col="black",position=position_dodge(width=0.5),shape=3, size=3, stroke=1) +
      geom_hline(yintercept = 0, size=0.5) +
      labs(x="Age-range", y=paste(parse(text="R^2"),"PGS")) +
      theme_classic() +
      mytheme +
      theme(
        panel.background = element_rect(fill = "white"),
        axis.text = element_text(size = 8),
        axis.title.y = element_text(size = 8, margin = margin(0,10,0,0)),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        plot.title = element_text(size =10),
        plot.subtitle = element_text(size = 10, hjust=0.5),
        axis.ticks.x = element_blank(),
        
      ) +
      labs(color='PGS',fill='PGS', alpha="p<.05", shape=FALSE)
  ) 
  
  pPCANoLegend = pPCA + theme(legend.position = "none",
                              strip.background = element_blank(),
                              strip.text.x = element_blank()) + scale_x_discrete(labels=c('30+', '35+', '40+', '45+','50+','55+','60+','65+','70+'),position = "bottom")
  
  pPCArsqNoLegend = pPCArsq + theme(legend.position = "none",
                                    strip.background = element_blank(),
                                    strip.text.x = element_blank()) + scale_x_discrete(labels=c('30+', '35+', '40+', '45+','50+','55+','60+','65+','70+'))
  
  (pPCABothnoLegend = cowplot::plot_grid(pPCANoLegend,
                                         pPCArsqNoLegend, nrow=2,rel_heights = c(2, 1),align= T))
  
  pPCA = pPCA +
    theme(strip.text.x = element_blank()) + 
    scale_x_discrete(labels=c('30+', '35+', '40+', '45+','50+','55+','60+','65+','70+'))
  
  pPCArsq = pPCArsq +
    theme(strip.text.x = element_blank()) + 
    scale_x_discrete(labels=c('30+', '35+', '40+', '45+','50+','55+','60+','65+','70+'))
  
  (pPCABoth = cowplot::plot_grid(pPCA,
                                 pPCArsq, nrow=2,rel_heights = c(2, 1),align= T))
  
  # pPCA= pPCA+ theme(legend.position = "none")
  
  #orig
  # ggsave(
  #   plot = pPCA,
  #   file = filename,
  #   width = 13,
  #   height = 10,
  #   units = "cm",
  #   dpi = 600
  # )
  
  #concatenated
  # ggsave(
  #   plot = pPCABoth,
  #   file = filename,
  #   width = 13,
  #   height = 10,
  #   units = "cm",
  #   dpi = 600
  # )
  
  #estimates
  # ggsave(
  #   plot = pEst,
  #   file = filename
  #   width = 12,
  #   height = 10,
  #   units = "cm",
  #   dpi = 600
  # )
  
  
  #every increase in age subset associated with .02sd drop
  summary(lm((abs(OUT1_PC1$estimate)*-1) ~ OUT1_PC1$agecut))
  OUT1_PC1$agegroupindex=c( rep(0, 4), rep(1, 4), rep(2, 4), rep(3, 4), rep(4, 4), rep(5, 4), rep(6, 4), rep(7, 4), rep(8, 4) )
  summary(lm((abs(OUT1_PC1$estimate)*-1) ~ OUT1_PC1$agegroupindex))
  
  
  # PCroll associations ----
  (pRoll=
     ggplot(OUT1_PCroll, aes(x=1:nrow(OUT1_PCroll), y=-log10(p.value), col=as.factor(ww), shape=as.factor(score))) +
     scale_alpha_discrete(range=c(0,1)) +
     geom_point() +
     geom_point(aes(alpha=as.factor(dotalpha)), col="black", shape=21, size=4,stroke=0.75) +
     geom_hline(yintercept = -log10(0.05), linetype=2) +
     geom_hline(yintercept = -log10(FDRthreshAPOE), linetype=3) +
     labs(x="Window") +
     theme_classic() +
     mytheme +
     theme(panel.background = element_rect(fill="#f3f3f3"),
           axis.text = element_text(size=8),
           axis.title.x = element_text(size=7, vjust=-3, margin=margin(0,30,0,0)),
           axis.title.y = element_text(size=7, margin = margin(0,10,0,0)),
           legend.text = element_text(size=8),
           legend.title = element_text(size=8),
           plot.title = element_text(size=10),
           plot.subtitle = element_text(size=10, hjust=0.5)
     ) +
     #change shape
     scale_shape_manual(values = c("Jansen"=16, "Kunkle"=15, "AD"=17, "Wightman"=4))
  )
  
  
  # PCroll rsquared ----
  (pRollRsq = OUT1_PCroll %>%
      
      ggplot(., aes(x=1:nrow(OUT1_PCroll), y=-log10(p.value), col=as.factor(ww), fill=as.factor(ww), shape=as.factor(score))) +
      scale_x_discrete() +
      scale_alpha_discrete(range = c(0, 1)) +
      ylim(c(maxscalersq*-1,0)) +
      geom_col(aes(alpha=as.factor(FDRsig), y=prsq*-1,x=1:nrow(OUT1_PCroll)), color=NA, position=position_dodge(width=.4), width=0.75) +
      geom_errorbar(aes(alpha=as.factor(FDRsig), ymin=prsq_lwr*-1,ymax=prsq_upr*-1),position=position_dodge(width=0.5),size=.3, width=0) +
      geom_point(data=OUT2_PCroll, aes(alpha=dotalphapoe, y=prsq*-1, x=1:nrow(OUT1_PCroll)),col="black",position=position_dodge(width=0.5),shape=3, size=3, stroke=1) +
      geom_hline(yintercept = 0, size=0.5) +
      
      labs(x="Age-range", y=paste(parse(text="R^2"),"PGS")) +
      theme_classic() +
      mytheme +
      theme(
        panel.background = element_rect(fill = "#f3f3f3"),
        axis.text = element_text(size = 8),
        axis.title.y = element_text(size = 8, margin = margin(0,10,0,0)),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        plot.title = element_text(size =10),
        plot.subtitle = element_text(size = 10, hjust=0.5),
        axis.ticks.x = element_blank(),
        
      ) +
      labs(color='PGS',fill='PGS', alpha="p<.05", shape=FALSE)
  ) 
  
  
  pRollNoLegend = pRoll + theme(legend.position = "none",
                                strip.background = element_blank(),
                                strip.text.x = element_blank()) + scale_x_discrete(labels=c('30+', '35+', '40+', '45+','50+','55+','60+','65+','70+'),position = "top")
  
  pRollRsqNoLegend = pRollRsq + theme(legend.position = "none",
                                      strip.background = element_blank(),
                                      strip.text.x = element_blank()) + scale_x_discrete(labels=c('30+', '35+', '40+', '45+','50+','55+','60+','65+','70+'))
  
  (pRollBoth = cowplot::plot_grid(pRollNoLegend,
                                  pRollRsqNoLegend, nrow=2,rel_heights = c(2, 1),align= T))
  
  (pRollBothLegend = cowplot::plot_grid(pRoll, 
                                        pRollRsq, nrow=2,rel_heights = c(2, 1),align= T))
  
  #orig
  # ggsave(
  #   plot = pRoll,
  #   file = filename,
  #   width = 13,
  #   height = 10,
  #   units = "cm",
  #   dpi = 600
  # )
  
  #newrsq
  # ggsave(
  #   plot = pRollBothLegend,
  #   file = filename,
  #   width = 13,
  #   height = 10,
  #   units = "cm",
  #   dpi = 600
  # )

}

