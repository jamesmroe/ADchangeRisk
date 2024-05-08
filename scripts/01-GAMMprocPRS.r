#========================================================================================#
# Written by James M Roe, Ph.D.
# Center for Lifespan Changes in Brain and Cognition, Department of Psychology
# University of Oslo, Norway
#========================================================================================#


#========================================================================================#
## Purpose: Run GAMM trajectory modelling example and PRS-AD association tests in paper
##          Script requires individual-level data as input and is not executable
#========================================================================================#


# rm(list=ls())
options(bitmapType='cairo')


#---load packages
tmp.packages = c("tidyverse","magrittr","gamm4","itsadug","numDeriv","gratia","mgcv","viridis","wesanderson","asbio","broom","cowplot","data.table","stringi")
tmpnew.packages = tmp.packages[!(tmp.packages %in% installed.packages()[,"Package"])]
if(length(tmpnew.packages)) {
  install.packages(tmpnew.packages)
}
sapply(tmp.packages, require, character.only = T)
sapply(tmp.packages, packageVersion)
rm(list=ls(pattern="tmp*"))


#---set dir
wd = "/cluster/projects/p274/projects/p040-ad_change/ADchangeRisk/reproduce"
#---set dir


#---make dirstruct
if (! dir.exists(wd)) { dir.create(wd) }
setwd(wd); library("here"); here()
plotdir = "plots"; if (! dir.exists(plotdir)) { dir.create(plotdir)}
resdir = "results"; if (! dir.exists(resdir)) { dir.create(resdir)}


#---load data
savefigs=0
nTime=2
agecut=30
load(file.path(wd,"data/minDF.Rda"))
DF = minDF
saveres=0


#chrome-extension://oemmndcbldboiebfnladdacbdfmadadm/https://adni.bitbucket.io/reference/docs/UCBERKELEYAV1451/UCBERKELEY_AV1451_Methods_Aug2018.pdf
#ADNI ROIs
# #Braak1
# L_entorhinal
# R_entorhinal

# # Braak2
# L_hippocampus
# R_hippocampus

# # Braak3
# L_parahippocampal
# R_parahippocampal 
# L_fusiform
# R_fusiform 
# L_lingual
# R_lingual
# L_amygdala
# R_amygdala

ROIS = c(
  #stage 2 (L R)
  "Left-Hippocampus.fs71",
  "Right-Hippocampus.fs71",
  #stage 1 (L R)
  "lh_entorhinal_volume.aparcnative71",
  "rh_entorhinal_volume.aparcnative71",
  #stage 3 subcortical (L R)
  "Left-Amygdala.fs71",
  "Right-Amygdala.fs71",
  #stage 3 cortical (L)
  "lh_parahippocampal_volume.aparcnative71",
  "lh_fusiform_volume.aparcnative71",
  "lh_lingual_volume.aparcnative71",
  #stage 4 cortical (R)
  "rh_parahippocampal_volume.aparcnative71",
  "rh_fusiform_volume.aparcnative71",
  "rh_lingual_volume.aparcnative71")

ROIS_CONTROL = c("EstimatedTotalIntraCranialVol.fs71")

#minDF
DF = DF %>% filter(visit_age >= agecut)

# DF$brainvar1=DF[[ROIS[1]]]
# DF$brainvar2=DF[[ROIS[2]]]
# DF$brainvar3=DF[[ROIS[3]]]
# DF$brainvar4=DF[[ROIS[4]]]
# DF$brainvar5=DF[[ROIS[5]]]
# DF$brainvar6=DF[[ROIS[6]]]
# DF$brainvar7=( DF[[ROIS[7]]] + DF[[ROIS[8]]] +
#                 DF[[ROIS[9]]] )/3
# DF$brainvar8=( DF[[ROIS[10]]] + DF[[ROIS[11]]] +
#     DF[[ROIS[12]]] )/3



#derived vars from minDF
#longitudinal only
DF = DF %>% group_by(subject_id) %>% mutate(TimeBsl = visit_age - min(visit_age),
                                            meanAge = mean(visit_age),
                                            nTimepoints= length(unique(visit_age)),
                                            PID = substr(subject_id, 1, 2),
                                            intervalFirstlast = max(visit_age) - min(visit_age)) %>% filter(nTimepoints>=nTime) %>% ungroup()
DF %<>% mutate(ICV = (DF$EstimatedTotalIntraCranialVol.fs71 - mean(EstimatedTotalIntraCranialVol.fs71)) / sd(DF$EstimatedTotalIntraCranialVol.fs71))


# LOOP ROIS  ------------------------------------------------
SS="S1"; iii=1
for (iii in 1) { #1:8
  
  print(iii)
  
  #set outputs
  if (iii == 1) { 
    roivec=0 
    AGERELSLOPE=OUTLISTINT=AGERELSLOPENOAPOE=OUTLISTAPOEINT=ABSCHANGESLOPE=ABSCHANGESLOPENOAPOE=list()
  }
  

  #--- select brainvar
  if (iii == 1) {
    ROI = "Left-Hippocampus.fs71"
    roi = paste(strsplit(ROI,"-")[[1]][1],strsplit(ROI,"-")[[1]][2])
    DF$brainvar = DF[[ROI]]
  } else if (iii == 2) {
    ROI = "Right-Hippocampus.fs71"
    roi = paste(strsplit(ROI,"-")[[1]][1],strsplit(ROI,"-")[[1]][2])
    DF$brainvar = DF[[ROI]]
  } 
  
  if (iii >=3 ) {
    braakaparc = 1
  } else { 
    braakaparc = 0 
  }
  
  
  #composite braak 3 rois
  if (braakaparc == 1) {
    print("braak composite ROI")
    
    
    if (iii == 3) { st="braak1ent"; h="lh"}
    if (iii == 4) { st="braak1ent"; h="rh"}
    if (iii == 5) { st="braak3amyg"; h="lh"}
    if (iii == 6) { st="braak3amyg"; h="rh"}
    if (iii == 7) { st="braak3cortx"; h="lh"}
    if (iii == 8) { st="braak3cortx"; h="rh"}
    if (h == "lh") { H = "Left" } else { H ="Right"}
    ver="71"
    
    if (st == "braak1ent") {
      roi = paste("braak stage 1ent (volume)",h, ver) 
      DF$brainvar = DF[[paste0(h,"_entorhinal_volume.aparcnative",ver)]]
    }
    
    if (st == "braak3amyg") {
      roi = paste("braak stage 3 amygdala (volume)",h, ver) 
      DF$brainvar = DF[[paste0(H,"-Amygdala.fs",ver)]]
      
    }
    
    if (st == "braak3cortx") {
      roi = paste("braak stage 3-amy (volume)",h, ver) 
      DF$brainvar = (DF[[paste0(h,"_parahippocampal_volume.aparcnative",ver)]] + DF[[paste0(h,"_fusiform_volume.aparcnative",ver)]] +
                       DF[[paste0(h,"_lingual_volume.aparcnative",ver)]])/3
      
    }
    
  }
  print(paste(roi,":",iii))
  roivec[iii]=roi
  
  

  
  #gamm w/ full random effects
  g = gamm4(brainvar ~ s(visit_age, k = 8) + subject_sex + mri_info_site_name + ICV, data = DF, random = ~ (1+visit_age | subject_id))
  g.sum = summary(g$gam)
  dug = itsadug::get_predictions(g$gam,
                                 cond = list(visit_age = seq(
                                   min(DF$visit_age, na.rm = T),
                                   max(DF$visit_age, na.rm =
                                         T),
                                   length.out = nrow(DF)
                                 ), subject_sex = "Female",
                                 mri_info_site_name = "ousSkyra",
                                 ICV = 0,
                                 se = T))
  
  predictions <- DF %>% 
    mutate(subject_sex = "Female", mri_info_site_name = "ousSkyra", ICV = 0) %>%
    select(visit_age,subject_sex,mri_info_site_name, ICV) %>% 
    predict(g$gam, newdata = ., se.fit = T)
  
  residualsg <- residuals(g$mer)
  DF$partial_residuals = predictions$fit + residualsg
  DF$fit = predictions$fit
  DF$sefit = predictions$se.fit
  DF$cifit = predictions$se.fit*1.96
  #plot(DF$visit_age, DF$partial_residuals)
  #points(DF$visit_age,DF$fit,col="blue")
  #plot(DF$visit_age,DF$brainvar)
  
  
  #correct for sex and scanner - ousAvanto is reflevel - bring to Skyra
  DF$brainvarcor = DF$brainvar
  DF$brainvarcor = ifelse(DF$subject_sex == "Male", DF$brainvarcor - g.sum$p.coeff["subject_sexMale"], DF$brainvarcor)
  DF$brainvarcor = ifelse(DF$mri_info_site_name == "ousAvanto", DF$brainvarcor + g.sum$p.coeff["mri_info_site_nameousSkyra"], DF$brainvarcor)
  DF$brainvarcor = ifelse(DF$mri_info_site_name == "ousPrisma", (DF$brainvarcor - g.sum$p.coeff["mri_info_site_nameousPrisma"]) + g.sum$p.coeff["mri_info_site_nameousSkyra"], DF$brainvarcor)
  
  
  
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
  
  
  
  #colour palette ---
  pal = wesanderson::wes_palettes$Zissou1
  #colour palette ---
  
  
  #set colour by hemi
  if (braakaparc == 0) {
    if (strsplit(roi," ")[[1]][1] == "Left") { 
      pointcol = pal[1]
    } else if (strsplit(roi," ")[[1]][1] == "Right") {
      pointcol = pal[3]
    }
  } else if (braakaparc == 1) {
      ROI=gsub("\\)","",gsub("\\(", ".", gsub(" ", "-", roivec[[iii]])))
    if (stringi::stri_reverse(substr(stringi::stri_reverse(roi),4,5)) == "lh") { 
      pointcol = pal[1]
    } else if (stringi::stri_reverse(substr(stringi::stri_reverse(roi),4,5)) == "rh") { 
      pointcol = pal[3]
    }
  }
  
  
  
  (fig1=ggplot() +
      geom_line(data=DF,aes(x=visit_age,brainvarcor,group=subject_id),color=pointcol,alpha=0.6, size=1) +
      geom_point(data=DF,aes(x=visit_age,brainvarcor,group=subject_id),color=pointcol,stat="identity",alpha=1, size=0.5) +
      geom_ribbon(data=dug,aes(x=visit_age,ymin=fit-CI,ymax=fit+CI),alpha=.7,show.legend=F,fill="dark grey") +
      geom_line(data=DF,aes(x=visit_age,y=fit),col="black") +
      # geom_line(data=dug,aes(x=visit_age,y=fit),col="red",size=2) + #check correct
      ylab("Volume") +
      ggtitle(roi) +
      theme_classic() + mytheme)
  
  
  
  #random slopes
  rr = ranef(g$mer)$subject_id
  length(unique(DF$subject_id))
  U = DF %>% filter(TimeBsl == 0)
  names(rr) = c("rInt", "rSlope")
  rr$subject_id = row.names(rr)
  U = merge(U, rr, by = "subject_id")
  
  
  #check outliers
  outlierThresh = function(U, outlierthresh, ranef) {

    g1 = gam(get(ranef) ~ s(meanAge, k = 8) + subject_sex + mri_info_site_name, data = U)
    Opt.age = U$meanAge
    Q = scale(g1$residuals)
    db.out = data.frame(U$subject_id,Q) #residuals scaled Z
    db.out$QQ = NA
    db.out$QQ [ (abs(Q)>outlierthresh)] = Q[(abs(Q)>outlierthresh)] #places outliers in new col
    outliers = db.out$U.subject_id[which(!is.na(db.out$QQ))] #indexes
    cat("\noutlier thresh =", outlierthresh, "SD from GAMM model: highlighting", length(outliers), "cases for potential removal\n\n")
    N.outliers = length(outliers)

    if (ranef == "rSlope") {
      U$Qcol1="cyan2"
      U$Qcol1[(U$subject_id %in% outliers)]="black"
    } else {
      U$Qcol2="cyan2"
      U$Qcol2[(U$subject_id %in% outliers)]="black"
    }
    return(list(outliers = as.numeric(as.character(outliers)),
                U))
  }
  ranef="rSlope"
  thr = 6
  check = outlierThresh(U, thr, "rSlope")[[2]] #outliers saved for later PGS tests

  ggplot(data=check,aes(x=meanAge,rSlope)) +
    geom_point(color=check$Qcol1,stat="identity",alpha=0.7) +
    geom_smooth(method="lm",col="black") +
    theme_classic() +
    ylab("rSlope") +
    theme(axis.line = element_blank())
   
  
  
  (fig2=ggplot() +
      geom_point(data=U,aes(x=meanAge,rSlope),color=pointcol,stat="identity",alpha=0.7) +
      # geom_smooth(data=U,aes(x=meanAge,rSlope),method = "gam",col="black") +
      geom_hline(yintercept = 0,linetype=2,size=1) +
      labs(y="rSlope",
           x="Age (mean)") +
      ggtitle("Age-relative change") + 
      ylab(expression("Additional Change (mm"^3*"/year)")) +
      theme_classic() + mytheme +
      theme(
        axis.title.y = element_text(color = "black", size = 19, vjust =-1, margin = margin(0,10,0,0)),
        axis.title.x = element_text(color = "black", size = 19, vjust = -2, margin = margin(0,20,20,0)),
        # axis.title.x = element_text(vjust = -4),
        plot.margin = unit(c(1,2,1,1), "lines")
      
      ) )

  
  
  #calculate derivatives
  dd=derivatives(object=g$gam,
                 newdata=U,
                 order = 1L,
                 n = nrow(U),
                 interval = c("confidence")
  )
  
  #add random effect to curve derivative to get vals of slope
  U$absChange=U$rSlope+dd$derivative
  
  U$dd=dd$derivative
  
  
  (fig3=ggplot() +
      geom_point(data=U,aes(x=meanAge,absChange),color=pointcol,stat="identity",alpha=0.7) +
      geom_smooth(data=U,aes(x=meanAge,absChange),method="loess",col="black", size=0.5) +
      geom_hline(yintercept = 0,linetype=2,size=1) +
      labs(
           x="Age (mean)") +
      ggtitle("Absolute change") +
      ylab(expression("Change (mm"^3*"/year)")) +
      theme_classic() + mytheme +
      theme(
        axis.title.y = element_text(color = "black", size = 19, vjust =-1, margin = margin(0,10,0,0)),
        axis.title.x = element_text(color = "black", size = 19, vjust = -2, margin = margin(0,20,20,0)),
        # axis.title.x = element_text(vjust = -4),
        plot.margin = unit(c(1,2,1,1), "lines")
        
      ) )
  

  # describe sample
  descSample = function(U, DF) {
    nU=nrow(U)
    cat("Unique N =",nU)
    nF = U %>% select(subject_id,subject_sex) %>% filter(subject_sex=="Female") %>% nrow()
    nM = U %>% select(subject_sex) %>% filter(subject_sex=="Male") %>% nrow()
    cat("Female/Male =",nF,"/",nM)

    #calc Nscans per person
    sample = DF %>% group_by(subject_id) %>% mutate(
      AgeBsl = min(visit_age),
      nScans = length(unique(visit_age)),
      TimeBsl = visit_age - min(visit_age)) %>% ungroup()

    sample %>% count(nScans)
    maxL=max(DF$TimeBsl)
    mL=mean(DF$TimeBsl)
    rangeL=range(DF$TimeBsl)
    sdL=sd(DF$TimeBsl)

    rangeTP=range(DF$nTimepoints)
    cat("TPs range from",rangeTP)
    mTP=mean(DF$nTimepoints)
    sdTP=sd(DF$nTimepoints)
    cat("mean(sd) N TPs =",mTP,paste0("(",sdTP,")"))


    mA=mean(sample$visit_age)
    sdA=sd(sample$visit_age)
    rangeA=range(sample$visit_age)

    tmpL=DF %>% filter(TimeBsl!=0)
    rangeL=range(tmpL$TimeBsl)

    sampleDesc = data.frame("Sample"= "LCBC",
                     "Unique N" = nU,
                     "N_obs" = nrow(sample),
                     "M_follow" = paste0(round(mL,1),"(",round(sdL,1),")"),
                     "Max_follow" = paste0(round(maxL,1), "(",round(sdL,1),")"),
                     "Range_follow" = paste0(round(rangeL[1],1),"-",round(rangeL[2],1)),
                     "M_TPs" = paste0(round(mTP,1),"(",round(sdTP,1),")"),
                     "Range_TPs" = paste0(rangeTP[1],"-",rangeTP[2]),
                     "M_Age" = paste0(round(mA,1), "(",paste0(round(rangeA[1],1),"-",round(rangeA[2],1)),")"),
                     "SD_Age" = sdA,
                     "Sex_FM" = paste0(nF,"/",nM)
    )

  return(sampleDesc)
  }
  sampleDesc = descSample(U,DF)

  
  PGS = data.frame(fread(paste0("/cluster/projects/p274/projects/p040-ad_change/GWASsumstats/new_PGS_for_LCBC/AD_Jansen/wAPOE/AD_Jansen.",SS,".profile")) %>% select(FID,SCORE) %>% mutate("PGSS_AD_Jansen"=SCORE),
                   fread(paste0("/cluster/projects/p274/projects/p040-ad_change/GWASsumstats/new_PGS_for_LCBC/AD_Wightman/wAPOE/AD_Wightman.",SS,".profile")) %>% select(FID,SCORE) %>% mutate("PGSS_AD_Wightman"=SCORE),
                   fread(paste0("/cluster/projects/p274/projects/p040-ad_change/GWASsumstats/new_PGS_for_LCBC/AD/wAPOE/AD.",SS,".profile")) %>% select(FID,SCORE) %>% mutate("PGSS_AD"=SCORE),
                   fread(paste0("/cluster/projects/p274/projects/p040-ad_change/GWASsumstats/new_PGS_for_LCBC/AD_kunkle/wAPOE/AD_Kunkle.",SS,".profile")) %>% select(FID,SCORE) %>% mutate("PGSS_AD_Kunkle"=SCORE),

                   fread(paste0("/cluster/projects/p274/projects/p040-ad_change/GWASsumstats/new_PGS_for_LCBC/AD_Jansen/woAPOE/AD_Jansen_woAPOE.",SS,".profile")) %>% select(FID,SCORE) %>% mutate("PGSS_AD_Jansen_woAPOE"=SCORE),
                   fread(paste0("/cluster/projects/p274/projects/p040-ad_change/GWASsumstats/new_PGS_for_LCBC/AD_Wightman/woAPOE/AD_Wightman_woAPOE.",SS,".profile")) %>% select(FID,SCORE) %>% mutate("PGSS_AD_Wightman_woAPOE"=SCORE),
                   fread(paste0("/cluster/projects/p274/projects/p040-ad_change/GWASsumstats/new_PGS_for_LCBC/AD/woAPOE/AD.",SS,".profile")) %>% select(FID,SCORE) %>% mutate("PGSS_AD_woAPOE"=SCORE),
                   fread(paste0("/cluster/projects/p274/projects/p040-ad_change/GWASsumstats/new_PGS_for_LCBC/AD_kunkle/woAPOE/AD_Kunkle_woAPOE.",SS,".profile")) %>% select(FID,SCORE) %>% mutate("PGSS_AD_Kunkle_woAPOE"=SCORE)
  )
  missing=PGS$FID[which(!PGS$FID %in% DF$Genetic_ID)]
  PGS$FID.1=NULL
  PGS$FID.2=NULL
  PGS$FID.3=NULL
  PGS$FID.4=NULL
  names(PGS)[1]="Genetic_ID"
  PGSS=merge(U,PGS,by="Genetic_ID")

  
  pcs <- fread("/tsd/p274/data/durable/projects/genetics/GAF/LCBC/LIFEBRAIN_LCBC_PCA.eigenvec", header=F) %>%
    setnames(., colnames(.), c("Genetic_ID", "IID", paste0("PC",1:20)) ) %>% select(1,3:12)


  PGSS %<>% select(-starts_with("GAF"))
  names(pcs)[2:11] = paste0("GAF_PC",1:10)
  PGSS = PGSS %>% merge(.,pcs)
  
  #add genetic data indicator to plot
  tmpU = U
  tmpU$strokealpha = factor(ifelse(U$subject_id %in% PGSS$subject_id, 1, 0),levels=c(0, 1))
  fig2gen = fig2 +
    scale_alpha_discrete(range = c(0, 1)) + geom_point(data = tmpU,aes(alpha=strokealpha,x=meanAge,rSlope),col="black",shape=21, size=2, stroke=0.2) +
    theme(legend.position = "none")
  

  
  (cp1gen=cowplot::plot_grid(fig1,fig3,fig2gen,ncol=3))
  
  
  
  if (savefigs == 1 ) {
    print("saving fig")
    filename = file.path(plotdir,(paste0("pAll-",ROI,".png")))
    ggsave(filename = filename,
    plot = cp1gen,
      width = 30,
      height = 13,
      dpi = 600,
      units = "cm"
    )
  }
  
  
  
  #--- full genetic data DF ----
  PGSSFULL = PGSS
  
  
  # LOOP increase lower age bound ------------------------------------------------
  for (agecut2 in seq(30,70, by=5)) {
    print(agecut2)
    
    if (agecut2 == 30) { jjj=1 }
    if (agecut2 == 35) { jjj=2 }
    if (agecut2 == 40) { jjj=3 }
    if (agecut2 == 45) { jjj=4 }
    if (agecut2 == 50) { jjj=5 }
    if (agecut2 == 55) { jjj=6 }
    if (agecut2 == 60) { jjj=7 }
    if (agecut2 == 65) { jjj=8 }
    if (agecut2 == 70) { jjj=9 }
    
    AGERELSLOPE[[jjj]]=list()
    AGERELSLOPENOAPOE[[jjj]]=list()
    ABSCHANGESLOPE[[jjj]]=list()
    ABSCHANGESLOPENOAPOE[[jjj]]=list()
    
    
    #--- raise lower age ---
    PGSS = PGSSFULL %>% filter(meanAge >= agecut2)
    #--- raise lower age ---

        
    runLM <- function(PGSS, Yvar, Xvar) {
  
      covariates = c(
        "scale(meanAge)" ,
        "subject_sex" ,
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
        "GAF_PC10" ,
        "nTimepoints",
        "intervalFirstlast"
      )
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
    
    #calc partial R sq
    getPartialR2 = function(model, score, apoe=1) {
      
      if (apoe == 1) {
        if (score == "Jansen") { model_drop = update(model, as.formula(paste0(". ~ . -scale(PGSS_AD_Jansen)"))) }
        if (score == "AD") { model_drop = update(model, as.formula(paste0(". ~ . -scale(PGSS_AD)"))) }
        if (score == "Kunkle") { model_drop = update(model, as.formula(paste0(". ~ . -scale(PGSS_AD_Kunkle)"))) }
        if (score == "Wightman") { model_drop = update(model, as.formula(paste0(". ~ . -scale(PGSS_AD_Wightman)"))) }
      } else if (apoe == 0) {
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
    
    
    # PRS-AD relchange models ----
    model1list = runLM(PGSS, "scale(rSlope)", "scale(PGSS_AD_Jansen)")
    model1 = model1list$model
    m1 = model1list$summary
    
    model2list = runLM(PGSS, "scale(rSlope)", "scale(PGSS_AD)")
    model2 = model2list$model
    m2 = model2list$summary
    
    model3list = runLM(PGSS, "scale(rSlope)", "scale(PGSS_AD_Kunkle)")
    model3 = model3list$model
    m3 = model3list$summary
    
    model4list = runLM(PGSS, "scale(rSlope)", "scale(PGSS_AD_Wightman)")
    model4 = model4list$model
    m4 = model4list$summary
      
      
    ## save to list ----
    AGERELSLOPE[[jjj]][[1]]=tidy(m1,
                                  conf.int = T) %>% mutate(roi=roi, score = "Jansen", agecut=agecut, N=nrow(PGSS), Nlife=nrow(DF), 
                                                              agecut2=agecut2,
                                                              getPartialR2(model1, "Jansen")[[1]],
                                                              apoe=1,
                                                              SS=SS)
    AGERELSLOPE[[jjj]][[2]]=tidy(m2,
                                  conf.int = T) %>% mutate(roi=roi, score = "AD", agecut=agecut, N=nrow(PGSS), Nlife=nrow(DF), 
                                                              agecut2=agecut2,
                                                              getPartialR2(model2, "AD")[[1]],
                                                              apoe=1,
                                                              SS=SS)
    AGERELSLOPE[[jjj]][[3]]=tidy(m3,
                                  conf.int = T) %>% mutate(roi=roi, score = "Kunkle", agecut=agecut, N=nrow(PGSS), Nlife=nrow(DF),
                                                              agecut2=agecut2,
                                                              getPartialR2(model3, "Kunkle")[[1]],
                                                              apoe=1,
                                                              SS=SS)
    AGERELSLOPE[[jjj]][[4]]=tidy(m4,
                                  conf.int = T) %>% mutate(roi=roi, score = "Wightman", agecut=agecut, N=nrow(PGSS), Nlife=nrow(DF),
                                                              agecut2=agecut2,
                                                              getPartialR2(model4, "Wightman")[[1]],
                                                              apoe=1,
                                                              SS=SS)
    
    
    # PRS-AD abschange models ----
    model1flist = runLM(PGSS, "scale(absChange)", "scale(PGSS_AD_Jansen)")
    model1f = model1flist$model
    m1f = model1flist$summary
    
    model2flist = runLM(PGSS, "scale(absChange)", "scale(PGSS_AD)")
    model2f = model2flist$model
    m2f = model2flist$summary
    
    model3flist = runLM(PGSS, "scale(absChange)", "scale(PGSS_AD_Kunkle)")
    model3f = model3flist$model
    m3f = model3flist$summary
    
    model4flist = runLM(PGSS, "scale(absChange)", "scale(PGSS_AD_Wightman)")
    model4f = model4flist$model
    m4f = model4flist$summary
    
    
    ## save to list ----
    ABSCHANGESLOPE[[jjj]][[1]]=tidy(m1f,
                                    conf.int = T) %>% mutate(roi=roi, score = "Jansen", agecut=agecut, N=nrow(PGSS), Nlife=nrow(DF),
                                                                    agecut2=agecut2,
                                                                    getPartialR2(model1f, "Jansen")[[1]],
                                                                    apoe=1,
                                                                    SS=SS)
    ABSCHANGESLOPE[[jjj]][[2]]=tidy(m2f,
                                    conf.int = T) %>% mutate(roi=roi, score = "AD", agecut=agecut, N=nrow(PGSS), Nlife=nrow(DF),
                                                                    agecut2=agecut2,
                                                                    getPartialR2(model2f, "AD")[[1]],
                                                                    apoe=1,
                                                                    SS=SS)
    ABSCHANGESLOPE[[jjj]][[3]]=tidy(m3f,
                                    conf.int = T) %>% mutate(roi=roi, score = "Kunkle", agecut=agecut, N=nrow(PGSS), Nlife=nrow(DF),
                                                                    agecut2=agecut2,
                                                                    getPartialR2(model3f, "Kunkle")[[1]],
                                                                    apoe=1,
                                                                    SS=SS)
    ABSCHANGESLOPE[[jjj]][[4]]=tidy(m4f,
                                    conf.int = T) %>% mutate(roi=roi, score = "Wightman", agecut=agecut, N=nrow(PGSS), Nlife=nrow(DF),
                                                                    agecut2=agecut2,
                                                                    getPartialR2(model4f, "Wightman")[[1]],
                                                                    apoe=1,
                                                                    SS=SS)
    
    
    # PRS-ADnoAPOE relchange models ----
    model1noAPlist = runLM(PGSS, "scale(rSlope)", "scale(PGSS_AD_Jansen_woAPOE)")
    model1noAP = model1noAPlist$model
    m1ap = model1noAPlist$summary
    
    model2noAPlist = runLM(PGSS, "scale(rSlope)", "scale(PGSS_AD_woAPOE)")
    model2noAP = model2noAPlist$model
    m2ap = model2noAPlist$summary
    
    model3noAPlist = runLM(PGSS, "scale(rSlope)", "scale(PGSS_AD_Kunkle_woAPOE)")
    model3noAP = model3noAPlist$model
    m3ap = model3noAPlist$summary
    
    model4noAPlist = runLM(PGSS, "scale(rSlope)", "scale(PGSS_AD_Wightman_woAPOE)")
    model4noAP = model4noAPlist$model
    m4ap = model4noAPlist$summary
  
      
    
    ## save to list ----
    AGERELSLOPENOAPOE[[jjj]][[1]]=tidy(m1ap,conf.int = T) %>% mutate(roi=roi, score = "Jansen", agecut=agecut, N=nrow(PGSS), Nlife=nrow(DF),
                                                                    agecut2=agecut2,
                                                                    getPartialR2(model1noAP, "Jansen",ap=0)[[1]],
                                                                    apoe=0,
                                                                    SS=SS)
    AGERELSLOPENOAPOE[[jjj]][[2]]=tidy(m2ap,conf.int = T) %>% mutate(roi=roi, score = "AD", agecut=agecut, N=nrow(PGSS), Nlife=nrow(DF),
                                                                    agecut2=agecut2,
                                                                    getPartialR2(model2noAP, "AD",ap=0)[[1]],
                                                                    apoe=0,
                                                                    SS=SS)
    AGERELSLOPENOAPOE[[jjj]][[3]]=tidy(m3ap,conf.int = T) %>% mutate(roi=roi, score = "Kunkle", agecut=agecut, N=nrow(PGSS), Nlife=nrow(DF),
                                                                    agecut2=agecut2,
                                                                    getPartialR2(model3noAP, "Kunkle",ap=0)[[1]],
                                                                    apoe=0,
                                                                    SS=SS)
    AGERELSLOPENOAPOE[[jjj]][[4]]=tidy(m4ap,conf.int = T) %>% mutate(roi=roi, score = "Wightman", agecut=agecut, N=nrow(PGSS), Nlife=nrow(DF),
                                                                    agecut2=agecut2,
                                                                    getPartialR2(model4noAP, "Wightman",ap=0)[[1]],
                                                                    apoe=0,
                                                                    SS=SS)
    
    
    
    
    # PRS-ADnoAPOE abschange models ----
    model1fnoAPlist = runLM(PGSS, "scale(absChange)", "scale(PGSS_AD_Jansen_woAPOE)")
    model1fnoAP = model1fnoAPlist$model
    m1fap = model1fnoAPlist$summary
    
    model2fnoAPlist = runLM(PGSS, "scale(absChange)", "scale(PGSS_AD_woAPOE)")
    model2fnoAP = model2fnoAPlist$model
    m2fap = model2fnoAPlist$summary
    
    model3fnoAPlist = runLM(PGSS, "scale(absChange)", "scale(PGSS_AD_Kunkle_woAPOE)")
    model3fnoAP = model3fnoAPlist$model
    m3fap = model3fnoAPlist$summary
    
    model4fnoAPlist = runLM(PGSS, "scale(absChange)", "scale(PGSS_AD_Wightman_woAPOE)")
    model4fnoAP = model4fnoAPlist$model
    m4fap = model4fnoAPlist$summary
    
      
    ## save to list ----
    ABSCHANGESLOPENOAPOE[[jjj]][[1]]=tidy(m1fap,conf.int = T) %>% mutate(roi=roi, score = "Jansen", agecut=agecut, N=nrow(PGSS), Nlife=nrow(DF),
                                                                            agecut2=agecut2,
                                                                            getPartialR2(model1fnoAP, "Jansen",ap=0)[[1]],
                                                                            apoe=0,
                                                                            SS=SS)
    ABSCHANGESLOPENOAPOE[[jjj]][[2]]=tidy(m2fap,
                                          conf.int = T) %>% mutate(roi=roi, score = "AD", agecut=agecut, N=nrow(PGSS), Nlife=nrow(DF),
                                                                            agecut2=agecut2,
                                                                            getPartialR2(model2fnoAP, "AD",ap=0)[[1]],
                                                                            apoe=0,
                                                                            SS=SS)
    ABSCHANGESLOPENOAPOE[[jjj]][[3]]=tidy(m3fap,
                                          conf.int = T) %>% mutate(roi=roi, score = "Kunkle", agecut=agecut, N=nrow(PGSS), Nlife=nrow(DF),
                                                                            agecut2=agecut2,
                                                                            getPartialR2(model3fnoAP, "Kunkle",ap=0)[[1]],
                                                                            apoe=0,
                                                                            SS=SS)
    ABSCHANGESLOPENOAPOE[[jjj]][[4]]=tidy(m4fap,
                                          conf.int = T) %>% mutate(roi=roi, score = "Wightman", agecut=agecut, N=nrow(PGSS), Nlife=nrow(DF),
                                                                            agecut2=agecut2,
                                                                            getPartialR2(model4fnoAP, "Wightman",ap=0)[[1]],
                                                                            apoe=0,
                                                                            SS=SS)
    
    
    # plot predictions ----
    predlm1 = predict(
      model1,
      newdata = data.frame(
        PGSS_AD_Jansen = (PGSS$PGSS_AD_Jansen),
        meanAge = mean(PGSS$meanAge),
        subject_sex = "Female",
        mri_info_site_name = "ousSkyra", 
        GAF_PC1 = mean(PGSS$GAF_PC1),
        GAF_PC2 = mean(PGSS$GAF_PC2),
        GAF_PC3 = mean(PGSS$GAF_PC3), 
        GAF_PC4 = mean(PGSS$GAF_PC4), 
        GAF_PC5 = mean(PGSS$GAF_PC5), 
        GAF_PC6 = mean(PGSS$GAF_PC6), 
        GAF_PC7 = mean(PGSS$GAF_PC7),
        GAF_PC8 = mean(PGSS$GAF_PC8), 
        GAF_PC9 = mean(PGSS$GAF_PC9),
        GAF_PC10 = mean(PGSS$GAF_PC10),
        nTimepoints=median(PGSS$nTimepoints),
        intervalFirstlast=mean(PGSS$intervalFirstlast)
      ),
      se.fit = T
    )
    
    PGSS %<>% mutate(predY1=predlm1$fit,
                    predSE1=predlm1$se.fit,
                    predCI1=predlm1$se.fit*1.96,
                    resid1=residuals(model1),
                    presid1=predY1 + resid1,
                    X1=PGSS_AD_Jansen
    )
    
    ggplot(data=PGSS,aes(x=scale(X1), y=(presid1))) +
      geom_point(size=1,alpha=0.8,col="#66ff00") +
      geom_smooth(method="lm") +
      theme_classic() +
      mytheme
    
  
    savePGSplot=0
    if (savePGSplot == 1) {
      if (ROI == "Left-Hippocampus.fs71" & agecut2 == 30) {
      
        (plotPred1 = ggplot(data=PGSS,aes(x=scale(X1), y=scale(presid1),col=meanAge)) +
          geom_point(aes(size=meanAge),alpha=0.5) + #,col=pointcol) +
          scale_color_viridis(option="A") +
          scale_size_continuous(limits=c(30, 90), breaks=seq(30, 90, by=20)) +
          geom_smooth(method="lm",alpha=0.5,col="black",se=T) +
          coord_fixed() +
          labs(x="PGS_AD_Jansen (Z)",
               y= paste("Additional Change (Z)")) +
          theme_classic() +
          mytheme) +
        theme(
          axis.title.y = element_text(color = "black", size = 19, vjust =-1, margin = margin(0,10,0,0)),
          axis.title.x = element_text(color = "black", size = 19, vjust = -2, margin = margin(0,20,20,0)),
          plot.margin = unit(c(1,2,1,1), "lines")
        )
      
      
      }
      if (ROI == "Right-Hippocampus.fs71" & agecut2 == 30) {
        (plotPred2 = ggplot(data=PGSS,aes(x=scale(X1), y=scale(presid1),col=meanAge)) +
           geom_point(aes(size=meanAge),alpha=0.5) +#,col=pointcol) +
           scale_color_viridis(option="A") +
           scale_size_continuous(limits=c(30, 90), breaks=seq(30, 90, by=20)) +
           geom_smooth(method="lm",alpha=0.5,col="black",se=T) +
           coord_fixed() +
           labs(x="PGS_AD_Jansen (Z)",
                y= paste("Additional Change (Z)")) +
           theme_classic() +
           mytheme +
           theme(
             axis.title.y = element_text(color = "black", size = 19, vjust =-1, margin = margin(0,10,0,0)),
             axis.title.x = element_text(color = "black", size = 19, vjust = -2, margin = margin(0,20,20,0)),
             plot.margin = unit(c(1,2,1,1), "lines")
           )) +
          labs(color='Age (mean)',
               size=NULL)
        
        #changes previous loop plot
        plotPred1 = plotPred1 + theme(plot.margin = unit(c(1,2,1,1), "lines")) + labs(color='Age (mean)',
                                                                          size=NULL)
        
        # common legend
        (multi_plot = ggarrange(plotPred1,plotPred2,
                                ncol = 2, nrow = 1,
                                common.legend = 2
        ))
        
    
        # ggsave(
        #   filename = filename,
        #   plot = multi_plot,
        #   width = 30,
        #   height = 12,
        #   dpi = 600,
        #   units = "cm"
        # )
      
      }
    }
  } 
  # LOOPEND increase lower age bound ----
  
  
  # concatenate save list ----
  roiname=gsub("\\)","",gsub("\\(", ".", gsub(" ", "-", roivec[[iii]])))
  OUT=list()
  OUT[[1]]=bind_rows(AGERELSLOPE, .id = "model")
  OUT[[2]]=bind_rows(AGERELSLOPENOAPOE, .id = "model")
  OUT[[3]]=bind_rows(ABSCHANGESLOPE, .id = "model")
  OUT[[4]]=bind_rows(ABSCHANGESLOPENOAPOE, .id = "model")
  
  
  # save results ----
  if (saveres == 1) {
    print("saving")
    save("roivec","OUT","AGERELSLOPE","AGERELSLOPENOAPOE","ABSCHANGESLOPE","ABSCHANGESLOPENOAPOE",file=paste0(resdir,"/PRSADmodels_",SS,"_",roiname,".publication.Rda"))
  }
}
# LOOPEND ROIS  ------------------------------------------------
