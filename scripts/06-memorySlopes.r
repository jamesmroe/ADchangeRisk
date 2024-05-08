#=======================================================================================================#
#
# Author: James M Roe, PhD.
# Center for Lifespan Changes in Brain and Cognition, University of Oslo

# Purpose: estimate CVLT slopes
#          define quadrant risk groups
#          test whether high PRS-AD high change individuals show more adult lifespan memory decline
#          Script requires individual-level data as input and is not executable
#
#=======================================================================================================#


rm(list=ls())
options(bitmapType = "cairo")
options(show.error.locations = TRUE)
b="/cluster/projects/p274/projects/p040-ad_change/ADchangeRisk"
setwd(b)


#---load packages----
tmp.packages = c("tidyverse", "magrittr", "data.table","broom","cowplot","ggpubr", "here",
                 "gamm4","itsadug","numDeriv","gratia","mgcv","viridis","wesanderson","asbio","RColorBrewer","PupillometryR")
tmpnew.packages = tmp.packages[!(tmp.packages %in% installed.packages()[,"Package"])]
if(length(tmpnew.packages)) {
  #install.packages(tmpnew.packages)
}
sapply(tmp.packages, require, character.only = T)
rm(list=ls(pattern="tmp*"))


here()
load(here("reproduce/data/minCogDF.Rda"))
saveplots=0


#---analysis steps-----
runGAMM = 1           #step1
loadGenIntersect = 1  #step2  
defineQuadGroups = 1  #step3
mergeData = 1  #step4



# ---- STEP1 runGAMM -----

if (runGAMM == 1) {
  
  #713 memory obs of 267 subs
  length(unique(cvlt$subject_id)); dim(cvlt)[1]
  PCA1 = prcomp(data.frame(cvlt$CVLT_A_Total,
                           cvlt$CVLT_5min_Free,
                           cvlt$CVLT_30min_Free), center = T,scale. = T)
  summary(PCA1)
  if (sign(cor(cvlt$visit_age,PCA1$x[,1])) == 1) {
    print("inversing PC1")
    cvlt$Memory = PCA1$x[,1]*-1
    plot(cvlt$visit_age,cvlt$Memory)
  } else {
    cvlt$Memory = PCA1$x[,1]
    plot(cvlt$visit_age,cvlt$Memory)
  }
  
  
  scalefac=max(cvlt$Memory)
  cvlt$Memory = (cvlt$Memory/scalefac)
  
  
  #longitudinal cog only
  cvlt = cvlt %>% group_by(subject_id) %>% mutate(TimeBsl = visit_age - min(visit_age),
                                                  meanAgeCog = mean(visit_age),
                                                  nTimepoints = length(unique(visit_age))) %>% filter(nTimepoints>=2)
  substr(cvlt$mri_info_folder,1,2) %>% unique() #subproject overview
  DF %<>% mutate(PID = substr(subject_id, 1, 2))
  
  
  #longitudinal brain only
  DF %<>% group_by(subject_id) %>% mutate(TimeBsl = visit_age - min(visit_age),
                                          meanAge = mean(visit_age),
                                          nTimepoints= length(unique(visit_age)),
                                          intervalFirstlast = max(visit_age) - min(visit_age)) %>% filter(nTimepoints>=2)
  DF$ICV = (DF$EstimatedTotalIntraCranialVol.fs71-mean(DF$EstimatedTotalIntraCranialVol.fs71)) / sd(DF$EstimatedTotalIntraCranialVol.fs71)
  
  
  # random slope on CVLT factor ----
  gcvlt = gamm4(Memory ~ s(visit_age, k = 8) + subject_sex,
                data = cvlt,
                random = ~ (1 + visit_age | subject_id))
  gcvlt.sum = summary(gcvlt$gam)
  rrcvlt = ranef(gcvlt$mer)$subject_id
  UC = cvlt %>% filter(TimeBsl == 0)
  names(rrcvlt) = c("rIntMemory", "rSlopeMemory")
  rrcvlt$subject_id = row.names(rrcvlt)
  UC = merge(UC, rrcvlt)
  


  # predictions for cog GAMMs ----
  dug = itsadug::get_predictions(gcvlt$gam,
                                  cond = list(visit_age = seq(
                                    min(cvlt$visit_age, na.rm = T),
                                    max(cvlt$visit_age, na.rm = T),
                                    length.out = nrow(cvlt)
                                  ), se = T))
  
  predictions <- cvlt %>%
    mutate(subject_sex = "Female") %>%
    select(subject_sex,visit_age) %>%
    predict(gcvlt$gam, newdata = .)
  residualsg <- residuals(gcvlt$mer)
  cvlt$partial_residuals = predictions + residualsg
  
  
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
  
  
  (fig1 = ggplot() +
      geom_line(data=cvlt,aes(x=visit_age,Memory,group=subject_id),color="light grey",alpha=0.4, size=1) +
      geom_point(data=cvlt,aes(x=visit_age,Memory,group=subject_id),color="gold2",stat="identity",alpha=1, size=0.6) +
      geom_ribbon(data=dug,aes(x=visit_age,ymin=fit-CI,ymax=fit+CI),alpha=.4,show.legend=F,fill="dark grey") +
      geom_line(data=dug,aes(x=visit_age,y=fit),col="black") +
      ylab("loading PC1") +
      ggtitle("CVLT memory performance (PC1)") +
      theme_classic() + mytheme +
      labs(x="Age") + theme(
        axis.title.y = element_text(color = "black", size = 19, vjust =-1, margin = margin(0,10,0,0)),
        axis.title.x = element_text(color = "black", size = 19, vjust = -2, margin = margin(0,20,20,0)),
        plot.margin = unit(c(1,2,1,1), "lines")
      ))

}


# ---- STEP2 loadGenIntersect -----

if (loadGenIntersect == 1) {

  #load PCA on PRS-AD scores + brain change sample
  load(here("reproduce/data/PCA_PGS.Rda"))
  dim(PCA_PGS)
  
  
  #plot slopes with genetic indicator
  UC$strokealpha=factor(ifelse(UC$subject_id %in% PCA_PGS$subject_id, 1, 0),levels=c(0, 1))
  sum(UC$strokealpha==1)
  
  
  (fig2 = ggplot() +
      geom_point(data=UC,aes(x=meanAge,rSlopeMemory),color="gold2",stat="identity",alpha=0.7) +
      geom_hline(yintercept = 0,linetype=2,size=1) +
      labs(y="Additional Memory Change (a.u.)",
           x="Age (mean)") +
      ggtitle("Age-relative change") +
      theme_classic() + mytheme +
      theme(
        axis.title.y = element_text(color = "black", size = 19, vjust =-1, margin = margin(0,10,0,0)),
        axis.title.x = element_text(color = "black", size = 19, vjust = -2, margin = margin(0,20,20,0)),
        # axis.title.x = element_text(vjust = -4),
        plot.margin = unit(c(1,2,1,1), "lines")) +
      
      scale_alpha_discrete(range = c(0, 1)) + 
      geom_point(data=UC,aes(alpha=strokealpha,x=meanAge,rSlopeMemory),col="black",shape=21, size=2, stroke=0.2) +
      theme(legend.position = "none")
    
  )
  
  
  #calculate derivatives - cvlt
  cvltdd = derivatives(object=gcvlt$gam,
                       newdata=UC,
                       order = 1L,
                       n = nrow(UC),
                       interval = c("confidence")
  )
  
  
  
  #add random effect to curve derivative to get vals of slope - cvlt
  UC$absChangeMemory = UC$rSlopeMemory + cvltdd$derivative
  
  
  (fig3 = ggplot() +
      geom_point(data=UC,aes(x=meanAge,absChangeMemory),col="gold2",stat="identity",alpha=0.7) +
      geom_smooth(data=UC,aes(x=meanAge,absChangeMemory),col="black",size=0.5) +
      geom_hline(yintercept = 0,linetype=2,size=1) +
      
      ylab("Memory change (a.u.)") +
      ggtitle("Absolute change") +
      
      theme_classic() + mytheme + 
      theme(axis.title.y = element_text(color = "black", size = 19, vjust =-1, margin = margin(0,10,0,0)),
            axis.title.x = element_text(color = "black", size = 19, vjust = -2, margin = margin(0,20,20,0)),
            # axis.title.x = element_text(vjust = -4),
            plot.margin = unit(c(1,2,1,1), "lines")) +
      scale_alpha_discrete(range = c(0, 1)) +
      theme(legend.position = "none")
  )
  
  (cp1=cowplot::plot_grid(fig1,fig3,fig2,ncol=3))

  
  # if (saveplots == 1) {
    # ggsave(
    # filename = filename,
    #   plot = cp1,
    #   width = 30,
    #   height = 13,
    #   dpi = 600,
    #   units = "cm"
    # )
  # }

}


# ---- STEP3 defineQuadGroups -----

if (defineQuadGroups == 1) {

  #link memory change with conjunction of PGS and brain change (spca1 = PC1relChange1-50)
  head(PCA_PGS)
  
  
  #we separated individuals into discrete groups based on the conjunction of brain and genetic risk
  #used partial association between PC1 of age-relative change (spca1)
  #and PC1 calculated across the four PRS-AD scores 
  modpcaPGS = (
    lm((spca1) ~ scale(PCA_PGS) + scale(meanAge) + subject_sex + mri_info_site_name +
         nTimepoints + intervalFirstlast +
         GAF_PC1 + GAF_PC2 + GAF_PC3 + GAF_PC4 + GAF_PC5 + GAF_PC6 + GAF_PC7 +
         GAF_PC8 + GAF_PC9 + GAF_PC10,
       data = PCA_PGS,
       na.action = na.exclude
    ))
  
  
  # plot predictions ----
  predlm1 = predict(
    modpcaPGS,
    newdata = data.frame(
      PCA_PGS = (PCA_PGS$PCA_PGS),
      meanAge = mean(PCA_PGS$meanAge),
      subject_sex = "Female",
      mri_info_site_name = "ousSkyra", 
      GAF_PC1 = mean(PCA_PGS$GAF_PC1,na.rm=T),
      GAF_PC2 = mean(PCA_PGS$GAF_PC2,na.rm=T),
      GAF_PC3 = mean(PCA_PGS$GAF_PC3,na.rm=T), 
      GAF_PC4 = mean(PCA_PGS$GAF_PC4,na.rm=T), 
      GAF_PC5 = mean(PCA_PGS$GAF_PC5,na.rm=T), 
      GAF_PC6 = mean(PCA_PGS$GAF_PC6,na.rm=T), 
      GAF_PC7 = mean(PCA_PGS$GAF_PC7,na.rm=T),
      GAF_PC8 = mean(PCA_PGS$GAF_PC8,na.rm=T), 
      GAF_PC9 = mean(PCA_PGS$GAF_PC9,na.rm=T),
      GAF_PC10 = mean(PCA_PGS$GAF_PC10,na.rm=T),
      nTimepoints=median(PCA_PGS$nTimepoints),
      intervalFirstlast=mean(PCA_PGS$intervalFirstlast)
    ),
    se.fit = T
  )
  
  PCA_PGS %<>% mutate(predY1 = predlm1$fit,
                      predSE1 = predlm1$se.fit,
                      predCI1 = predlm1$se.fit*1.96,
                      resid1 = residuals(modpcaPGS),
                      presid1 = predY1 + resid1,
                      X1 = PCA_PGS)
  
  
  # define quadrant groups ----
  x_mid_scale <- mean(c(max(scale(PCA_PGS$PCA_PGS), na.rm = TRUE), 
                        min(scale(PCA_PGS$PCA_PGS), na.rm = TRUE)))
  
  y_mid_scale <- mean(c(max(scale(PCA_PGS$presid1), na.rm = TRUE), 
                        min(scale(PCA_PGS$presid1), na.rm = TRUE)))
  
  PCA_PGS %<>% 
    mutate(quadrantpresid = case_when(scale(PCA_PGS) > x_mid_scale & scale(presid1) > y_mid_scale   ~ "Q1",
                                      scale(PCA_PGS) <= x_mid_scale & scale(presid1) > y_mid_scale  ~ "Q2",
                                      scale(PCA_PGS) <= x_mid_scale & scale(presid1) <= y_mid_scale ~ "Q3",
                                      TRUE                                         ~ "Q4"))
  
  # quadrant group plot ----
  (pPCA = PCA_PGS %>% 
      ggplot(data=.,aes(x=scale(X1), y=scale(presid1),size=meanAge,col=quadrantpresid)) +
      geom_point(alpha=0.5) +
      scale_size_continuous(limits=c(30, 90), breaks=seq(30, 90, by=20)) +
      geom_vline(xintercept = x_mid_scale) + # plot vertical line
      geom_hline(yintercept = y_mid_scale) + # plot horizontal line
      coord_fixed() +
      labs(x="PCA_PGS_AD (Z)",
           y= paste("PC1 Age-relative change (features 1-50)")) +
      # ggtitle("Geetic AD risk\n(age corrected)") +
      theme_classic() +
      scale_color_brewer(palette="Dark2") +
      # scale_color_manual(values=met.brewer("Derain", 4)) +
      mytheme +
      theme(legend.position = "top") +
      theme(
        axis.title.y = element_text(color = "black", size = 19, vjust =-1, margin = margin(0,10,0,0)),
        axis.title.x = element_text(color = "black", size = 19, vjust = -2, margin = margin(0,20,20,0)),
        plot.margin = unit(c(1,2,1,1), "lines")
      )
  )
}


# ---- STEP4 mergeData -----

if (mergeData == 1) {

  # link brain and memory change ----
  tmp = UC %>% select(-starts_with("GAF"), -"nTimepoints")
  tmp = merge(PCA_PGS,tmp)
  BOTH = merge(PCA_PGS,tmp)
  
  
  # link APOE ----
  APOE = rbind(fread("/tsd/p274/data/durable/projects/genetics/APOE_imputed/LCBC/LCBC_b1_APOE_count.txt"),
               fread("/tsd/p274/data/durable/projects/genetics/APOE_imputed/LCBC/LCBC_b2_APOE_count.txt"))
  APOE %<>% mutate(carrier = ifelse(E4>=1,1,0))
  APOE$IID = NULL
  names(APOE)[1] = "Genetic_ID"
  BOTH = merge(BOTH,APOE)
  
  
  #trace back groups to baseline memory data
  traceMem = cvlt[cvlt$subject_id %in% BOTH$subject_id,]
  traceMem = merge(traceMem, BOTH %>% select(subject_id, quadrantpresid))
  traceMemBase = traceMem %>% filter(TimeBsl==0)
  traceMemBase$baseline_CVLT_A_Total=traceMemBase$CVLT_A_Total
  traceMemBase$baseline_CVLT_5min_Free=traceMemBase$CVLT_5min_Free
  traceMemBase$baseline_CVLT_30min_Free=traceMemBase$CVLT_30min_Free
  traceMemBase$baseline_Memory=traceMemBase$Memory
  BOTH = merge(BOTH, traceMemBase %>% select(subject_id, contains("baseline_")))

}


# group tests on lifespan memory change ----
# variable input for age-range
memSlopesCompare = function(x, filter_by, max_age, min_age, yvar, type, scale = F) {
  
  if (type == "older") {
    #set lower age bound (raise lower age)
    x %<>%
      filter({{filter_by}} >= min_age)
    
  } else if (type == "younger") {
    #set upper and lower age bound ()
    x %<>%
      filter({{filter_by}} >= min_age)
    x %<>%
      filter({{filter_by}} <= max_age)
    
  }
  
  x$yvar=x[[yvar]]
  if (scale == T) {
    x$yvar=scale(x[[yvar]])
  }
  
  #set Q4 to intercept
  x$quadrantpresid = relevel(factor(x$quadrantpresid),"Q4")
  modelquad=lm((yvar) ~ quadrantpresid +
                 scale(meanAge) + subject_sex +
                 # scale(spca1) +
                 # scale(baseline_Memory) +
                 # E4 +
                 carrier +
                 # scale(PCA_PGS) +
                 nTimepoints + intervalFirstlast,
               data = x,
               na.action = na.exclude
  )
  smodquad=summary(modelquad)
  
  
  #track group sizes and ages
  print(dim(x))
  N = dim(x)[1]
  N1 = sum(x$quadrantpresid == "Q1")
  N4 = sum(x$quadrantpresid == "Q4")
  minage = min(x$meanAge)
  maxage = max(x$meanAge)
  
  
  #predict at set covs
  predlm1quad = predict(
    modelquad,
    newdata = data.frame(
      quadrantpresid = x$quadrantpresid,
      meanAge = mean(x$meanAge),
      subject_sex = "Female",
      mri_info_site_name = "ousSkyra",
      PCA_PGS = mean(x$PCA_PGS),
      E4 = 0,
      carrier = 0,
      nTimepoints = median(x$nTimepoints),
      intervalFirstlast = mean(x$intervalFirstlast),
      baseline_Memory = mean(x$baseline_Memory),
      spca1 = mean(x$spca1)
    ),
    se.fit = T
  )
  
  
  #overwrite data with residualized effects for plotting
  x %<>% mutate(predY1quad = predlm1quad$fit,
                predSE1quad = predlm1quad$se.fit,
                predCI1quad = predlm1quad$se.fit*1.96,
                resid1quad = residuals(modelquad),
                presid1quad = predY1quad + resid1quad,
                X1 = x$quadrantpresid)
  
  
  return(list(x, smodquad, modelquad, minage, maxage, N, N1, N4))
}


#dims going into analysis
dim(BOTH); dim(cvlt); dim(UC)


# MAIN RESULTS ----
# on average across the adult lifespan
memSlopesCompare(x = BOTH, filter_by = meanAge, max_age = 89, min_age = 30, yvar = "rSlopeMemory", type = "older", scale = F)[[2]]
memSlopesCompare(x = BOTH, filter_by = meanAge, max_age = 89, min_age = 30, yvar = "absChangeMemory", type = "older", scale = F)[[2]]


#write residualized effects to df and plot main results
res1 = memSlopesCompare(x = BOTH, filter_by = meanAge, max_age = 89, min_age = 30, yvar = "rSlopeMemory", type = "older", scale = F)[[1]]
res2 = memSlopesCompare(x = BOTH, filter_by = meanAge, max_age = 89, min_age = 30, yvar = "absChangeMemory", type = "older", scale = F)[[1]]
res1$newvar = as.numeric(substr(res1$quadrantpresid, 2, 2))
res2$newvar = as.numeric(substr(res2$quadrantpresid, 2, 2))
res1 %<>% mutate(quadrantpresid = factor(quadrantpresid, levels=c("Q1", "Q2", "Q3", "Q4")))
res2 %<>% mutate(quadrantpresid = factor(quadrantpresid, levels=c("Q1", "Q2", "Q3", "Q4")))



recalcPlot = function(res1, res2) {
  
  #recalc distribution plots
  pal = brewer.pal(n=4,"Dark2")
  means1 = res1 %>% group_by(quadrantpresid) %>% summarise(mean=mean(presid1quad))
  means2 = res2 %>% group_by(quadrantpresid) %>% summarise(mean=mean(presid1quad))
  

  set.seed(123)
  (p1 = ggplot(res1,
               aes(x=(factor(quadrantpresid)), y=presid1quad, col=factor(quadrantpresid))) +
      geom_jitter(alpha=0.7, position = position_jitter(seed = 123,width = 0.1)) +
      geom_segment(aes(x=0.75, y=means1$mean[1], xend=1.25, yend=means1$mean[1]), size=.75, color="black") +#color=pal[1]) +
      geom_segment(aes(x=1.75, y=means1$mean[2], xend=2.25, yend=means1$mean[2]), size=.75, color="black") +#color=pal[2]) +
      geom_segment(aes(x=2.75, y=means1$mean[3], xend=3.25, yend=means1$mean[3]), size=.75, color="black") +#color=pal[3]) +
      geom_segment(aes(x=3.75, y=means1$mean[4], xend=4.25, yend=means1$mean[4]), size=.75, color="black") +#color=pal[4]) +
      theme_classic() + mytheme +
      scale_color_brewer(palette="Dark2") +
      theme(legend.position = "top") +
      labs(y="memory change (corrected)",
           x=NULL)
  )
  
  (p2 = ggplot(res2,
               aes(x=factor(quadrantpresid), y=presid1quad, col=factor(quadrantpresid))) +
      geom_jitter(alpha=0.7,position = position_jitter(seed = 123,width = 0.1)) +
      geom_segment(aes(x=0.75, y=means2$mean[1], xend=1.25, yend=means2$mean[1]), size=.75, color="black") +#color=pal[1]) +
      geom_segment(aes(x=1.75, y=means2$mean[2], xend=2.25, yend=means2$mean[2]), size=.75, color="black") +#color=pal[2]) +
      geom_segment(aes(x=2.75, y=means2$mean[3], xend=3.25, yend=means2$mean[3]), size=.75, color="black") +#color=pal[3]) +
      geom_segment(aes(x=3.75, y=means2$mean[4], xend=4.25, yend=means2$mean[4]), size=.75, color="black") +#color=pal[4]) +
      theme_classic() + mytheme +
      scale_color_brewer(palette="Dark2") +
      theme(legend.position = "top") +
      labs(y="memory change (corrected)",
           x=NULL)
  )
  
  
  Q1res1 = res1[res1$quadrantpresid=="Q1",]
  Q4res1 = res1[res1$quadrantpresid=="Q4",]
  Q2res1 = res1[res1$quadrantpresid=="Q2",]
  Q3res1 = res1[res1$quadrantpresid=="Q3",]
  Q1res1$zero=0
  Q4res1$zero=0
  Q2res1$zero=0
  Q3res1$zero=0
  
  Q1res2 = res2[res2$quadrantpresid=="Q1",]
  Q4res2 = res2[res2$quadrantpresid=="Q4",]
  Q2res2 = res2[res2$quadrantpresid=="Q2",]
  Q3res2 = res2[res2$quadrantpresid=="Q3",]
  Q1res2$zero=0
  Q4res2$zero=0
  Q2res2$zero=0
  Q3res2$zero=0
  
  
  #NB! to align plots q3 and q4 must be transparent
  (dens3 <- ggplot() +
      PupillometryR::geom_flat_violin(data=Q1res1,aes(x=zero,y=presid1quad),alpha = 0.3, col = pal[1], fill = pal[1]) +
      PupillometryR::geom_flat_violin(data=Q4res1,aes(x=zero,y=presid1quad),alpha = 0.3, col = pal[4], fill = pal[4]) +
      
      PupillometryR::geom_flat_violin(data=Q2res1,aes(x=zero,y=presid1quad),alpha = 0.0, col = NA, fill = NA) + #col = pal[2]
      PupillometryR::geom_flat_violin(data=Q3res1,aes(x=zero,y=presid1quad),alpha = 0.0, col = NA, fill = NA) + #col = pal[3]
      theme_void())
  
  (dens4 <- ggplot() +
      PupillometryR::geom_flat_violin(data=Q1res2,aes(x=zero,y=presid1quad),alpha = 0.3, col = pal[1], fill = pal[1]) +
      PupillometryR::geom_flat_violin(data=Q4res2,aes(x=zero,y=presid1quad),alpha = 0.3, col = pal[4], fill = pal[4]) +
      
      PupillometryR::geom_flat_violin(data=Q2res2,aes(x=zero,y=presid1quad),alpha = 0.0, col = NA, fill = NA) + #col = pal[2]
      PupillometryR::geom_flat_violin(data=Q3res2,aes(x=zero,y=presid1quad),alpha = 0.0, col = NA, fill = NA) + #col = pal[3]
      theme_void())
  
  
  (pAll1 = p1 + dens3 +
      patchwork::plot_layout(ncol = 2, nrow = 1, widths = c(2, 0.75), heights = c(1, 4))
  )
  (pAll2 = p2 + dens4 +
      patchwork::plot_layout(ncol = 2, nrow = 1, widths = c(2, 0.75), heights = c(1, 4))
  )
  
  test1 = memSlopesCompare(x = BOTH, filter_by = meanAge, max_age = 89, min_age = 30, yvar = "rSlopeMemory", type = "older", scale = F) [[2]]
  test2 = memSlopesCompare(x = BOTH, filter_by = meanAge, max_age = 89, min_age = 30, yvar = "absChangeMemory", type = "older", scale = F) [[2]]
  
  return(list(pAll1 = pAll1, 
              pAll2 = pAll2,
              test1 = test1,
              test2 = test2))
}
pOut = recalcPlot(res1,res2)
pOut$pAll1
pOut$pAll2
pOut$test1
pOut$test2


# change tests in progressively older ----
runAgeSubsetOlder = function(lastage=89, scaleLM = F) {
  c=0
  for (agecut in (seq(30, 60, by=5))) {
    print(agecut)
    c=c+1
    if (c == 1) {
      
      #set outputs
      N = N1 = N4 = c()
      lastages = agecuts = c()
      tmp1slope = tmp2abs = list()
      pvals_mod1_1 = pvals_mod1_2 = pvals_mod1_3 = pvals_mod2_1 = pvals_mod2_2 = pvals_mod2_3 = c()
      bvals_mod1_1 = bvals_mod1_2 = bvals_mod1_3 = bvals_mod2_1 = bvals_mod2_2 = bvals_mod2_3 = c()
      civals_mod1_1 = civals_mod1_2 = civals_mod1_3 = civals_mod2_1 = civals_mod2_2 = civals_mod2_3 = c()
      minages = maxages = c()
      
    }
    
    #call lm function
    memSlopes1 = memSlopesCompare(x = BOTH, filter_by = meanAge, max_age = lastage, min_age = agecut, yvar = "rSlopeMemory", type = "older", scale = scaleLM)
    memSlopes2 = memSlopesCompare(x = BOTH, filter_by = meanAge, max_age = lastage, min_age = agecut, yvar = "absChangeMemory", type = "older", scale = scaleLM)
    
    #store lm results
    tmp1slope[[c]] = tidy(memSlopes1[[2]])
    tmp2abs[[c]] = tidy(memSlopes2[[2]])
    
    #extract lm results
    pvals_mod1_1[c] =  tmp1slope[[c]]$p.value[tmp1slope[[c]]$term == "quadrantpresidQ1"] 
    pvals_mod1_2[c] =  tmp1slope[[c]]$p.value[tmp1slope[[c]]$term == "quadrantpresidQ2"] 
    pvals_mod1_3[c] =  tmp1slope[[c]]$p.value[tmp1slope[[c]]$term == "quadrantpresidQ3"] 
    pvals_mod2_1[c] = tmp2abs[[c]]$p.value[tmp2abs[[c]]$term == "quadrantpresidQ1"]
    pvals_mod2_2[c] = tmp2abs[[c]]$p.value[tmp2abs[[c]]$term == "quadrantpresidQ2"]
    pvals_mod2_3[c] = tmp2abs[[c]]$p.value[tmp2abs[[c]]$term == "quadrantpresidQ3"]
    
    bvals_mod1_1[c] =  tmp1slope[[c]]$estimate[tmp1slope[[c]]$term == "quadrantpresidQ1"] 
    bvals_mod1_2[c] =  tmp1slope[[c]]$estimate[tmp1slope[[c]]$term == "quadrantpresidQ2"] 
    bvals_mod1_3[c] =  tmp1slope[[c]]$estimate[tmp1slope[[c]]$term == "quadrantpresidQ3"] 
    bvals_mod2_1[c] = tmp2abs[[c]]$estimate[tmp2abs[[c]]$term == "quadrantpresidQ1"]
    bvals_mod2_2[c] = tmp2abs[[c]]$estimate[tmp2abs[[c]]$term == "quadrantpresidQ2"]
    bvals_mod2_3[c] = tmp2abs[[c]]$estimate[tmp2abs[[c]]$term == "quadrantpresidQ3"]
    
    civals_mod1_1[c] =  tmp1slope[[c]]$std.error[tmp1slope[[c]]$term == "quadrantpresidQ1"] *1.96
    civals_mod1_2[c] =  tmp1slope[[c]]$std.error[tmp1slope[[c]]$term == "quadrantpresidQ2"] *1.96
    civals_mod1_3[c] =  tmp1slope[[c]]$std.error[tmp1slope[[c]]$term == "quadrantpresidQ3"] *1.96
    civals_mod2_1[c] = tmp2abs[[c]]$std.error[tmp2abs[[c]]$term == "quadrantpresidQ1"] *1.96
    civals_mod2_2[c] = tmp2abs[[c]]$std.error[tmp2abs[[c]]$term == "quadrantpresidQ2"] *1.96
    civals_mod2_3[c] = tmp2abs[[c]]$std.error[tmp2abs[[c]]$term == "quadrantpresidQ3"] *1.96
    
    #track group sizes and ages
    minages[c] = memSlopes1[[4]]
    maxages[c] = memSlopes1[[5]]
    N[c] = memSlopes1[[6]]
    N1[c] = memSlopes1[[7]]
    N4[c] = memSlopes1[[8]]
    lastages[c] = lastage
    agecuts[c] = agecut
  }
  
  df_out = data.frame(pvals_mod1_1, pvals_mod2_1, 
                      bvals_mod1_1, bvals_mod2_1,
                      civals_mod1_1, civals_mod2_1,
                      lastages, agecuts, minages, maxages,
                      agerange = paste0(agecuts, "-", lastages),
                      N, N1, N4)
  
  
  df_out$agerange = factor(df_out$agerange, levels= c("30-89", "35-89", "40-89", "45-89", "50-89", "55-89", "60-89", "65-89", "70-89"))
  
  
  #colour code signifcance
  df_out$colourcode1 = ifelse(pvals_mod1_1 < .05, "pink", "grey")
  df_out$dotalpha1 = factor(ifelse(pvals_mod1_1 < .05, 1, 0), levels = c(0,1))
  df_out$colourcode2 = ifelse(pvals_mod2_1 < .05, "pink", "grey")
  df_out$dotalpha2 = factor(ifelse(pvals_mod2_1 < .05, 1, 0), levels = c(0,1))
  

  #fix effect direction to show decline effect
  #because Q4 is set to intercept, betas are Q1 > Q4 (i.e. positive)
  df_out$bvals_mod1_1 = df_out$bvals_mod1_1*-1
  df_out$civals_mod1_1 = df_out$civals_mod1_1*-1
  df_out$bvals_mod2_1 = df_out$bvals_mod2_1*-1
  df_out$civals_mod2_1 = df_out$civals_mod2_1*-1
  
  
  return(df_out)
}


df_outolder = runAgeSubsetOlder(scaleLM = T)
df_outolder$agerange


# change estimates in different age subsets
runAgeSubsets = function(firstage=30, scaleLM = F) {
  c=0
  for (agecut in rev(seq(65, 85, by=5))) {
    c=c+1
    if (c == 1) {
      
      #set outputs
      N = N1 = N4 = c()
      firstages = agecuts = c()
      tmp1slope = tmp2abs = list()
      pvals_mod1_1 = pvals_mod1_2 = pvals_mod1_3 = pvals_mod2_1 = pvals_mod2_2 = pvals_mod2_3 = c()
      bvals_mod1_1 = bvals_mod1_2 = bvals_mod1_3 = bvals_mod2_1 = bvals_mod2_2 = bvals_mod2_3 = c()
      civals_mod1_1 = civals_mod1_2 = civals_mod1_3 = civals_mod2_1 = civals_mod2_2 = civals_mod2_3 = c()
      minages = maxages = c()
      
    }
    
    #call lm function
    memSlopes1 = memSlopesCompare(x = BOTH, filter_by = meanAge, max_age = agecut, min_age = firstage, yvar = "rSlopeMemory", type = "younger", scale = scaleLM)
    memSlopes2 = memSlopesCompare(x = BOTH, filter_by = meanAge, max_age = agecut, min_age = firstage, yvar = "absChangeMemory", type = "younger", scale = scaleLM)
    
    #store lm results
    tmp1slope[[c]] = tidy(memSlopes1[[2]])
    tmp2abs[[c]] = tidy(memSlopes2[[2]])
    
    #extract lm results
    pvals_mod1_1[c] =  tmp1slope[[c]]$p.value[tmp1slope[[c]]$term == "quadrantpresidQ1"] 
    pvals_mod1_2[c] =  tmp1slope[[c]]$p.value[tmp1slope[[c]]$term == "quadrantpresidQ2"] 
    pvals_mod1_3[c] =  tmp1slope[[c]]$p.value[tmp1slope[[c]]$term == "quadrantpresidQ3"] 
    pvals_mod2_1[c] = tmp2abs[[c]]$p.value[tmp2abs[[c]]$term == "quadrantpresidQ1"]
    pvals_mod2_2[c] = tmp2abs[[c]]$p.value[tmp2abs[[c]]$term == "quadrantpresidQ2"]
    pvals_mod2_3[c] = tmp2abs[[c]]$p.value[tmp2abs[[c]]$term == "quadrantpresidQ3"]
    
    bvals_mod1_1[c] =  tmp1slope[[c]]$estimate[tmp1slope[[c]]$term == "quadrantpresidQ1"] 
    bvals_mod1_2[c] =  tmp1slope[[c]]$estimate[tmp1slope[[c]]$term == "quadrantpresidQ2"] 
    bvals_mod1_3[c] =  tmp1slope[[c]]$estimate[tmp1slope[[c]]$term == "quadrantpresidQ3"] 
    bvals_mod2_1[c] = tmp2abs[[c]]$estimate[tmp2abs[[c]]$term == "quadrantpresidQ1"]
    bvals_mod2_2[c] = tmp2abs[[c]]$estimate[tmp2abs[[c]]$term == "quadrantpresidQ2"]
    bvals_mod2_3[c] = tmp2abs[[c]]$estimate[tmp2abs[[c]]$term == "quadrantpresidQ3"]
    
    civals_mod1_1[c] =  tmp1slope[[c]]$std.error[tmp1slope[[c]]$term == "quadrantpresidQ1"] *1.96
    civals_mod1_2[c] =  tmp1slope[[c]]$std.error[tmp1slope[[c]]$term == "quadrantpresidQ2"] *1.96
    civals_mod1_3[c] =  tmp1slope[[c]]$std.error[tmp1slope[[c]]$term == "quadrantpresidQ3"] *1.96
    civals_mod2_1[c] = tmp2abs[[c]]$std.error[tmp2abs[[c]]$term == "quadrantpresidQ1"] *1.96
    civals_mod2_2[c] = tmp2abs[[c]]$std.error[tmp2abs[[c]]$term == "quadrantpresidQ2"] *1.96
    civals_mod2_3[c] = tmp2abs[[c]]$std.error[tmp2abs[[c]]$term == "quadrantpresidQ3"] *1.96
    
    #track group sizes and ages
    N[c] = memSlopes1[[6]]
    N1[c] = memSlopes1[[7]]
    N4[c] = memSlopes1[[8]]
    minages[c] = memSlopes1[[4]]
    maxages[c] = memSlopes1[[5]]
    firstages[c] = firstage
    agecuts[c] = agecut
  }
  
  df_out = data.frame(pvals_mod1_1, pvals_mod2_1, 
                      bvals_mod1_1, bvals_mod2_1,
                      civals_mod1_1, civals_mod2_1,
                      firstages, agecuts, minages, maxages,
                      agerange = paste0(firstages, "-", agecuts),
                      N, N1, N4)
  
  if (firstage == 30) {
    df_out$agerange = factor(df_out$agerange, levels= c("30-65", "30-70", "30-75", "30-80", "30-85"))
  } else if (firstage == 40) {
    df_out$agerange = factor(df_out$agerange, levels= c("40-65", "40-70", "40-75", "40-80", "40-85"))
  } else if (firstage == 50) {
    df_out$agerange = factor(df_out$agerange, levels= c("50-65", "50-70", "50-75", "50-80", "50-85"))
  }
  
  
  #colour code significance
  df_out$colourcode1 = ifelse(pvals_mod1_1<.05, "pink", "grey")
  df_out$dotalpha1 = factor(ifelse(pvals_mod1_1<.05, 1, 0), levels = c(0,1))
  df_out$colourcode2 = ifelse(pvals_mod2_1<.05, "pink", "grey")
  df_out$dotalpha2 = factor(ifelse(pvals_mod2_1<.05, 1, 0), levels = c(0,1))
  
  
  #fix effect direction to show decline effect
  #because Q4 is set to intercept, betas are Q1 > Q4 (i.e. positive)
  df_out$bvals_mod1_1 = df_out$bvals_mod1_1*-1
  df_out$civals_mod1_1 = df_out$civals_mod1_1*-1
  df_out$bvals_mod2_1 = df_out$bvals_mod2_1*-1
  df_out$civals_mod2_1 = df_out$civals_mod2_1*-1
  
  
  return(df_out)
}


df_younger = rbind(runAgeSubsets(firstage = 30, scaleLM = T), 
                   runAgeSubsets(firstage = 40, scaleLM = T), 
                   runAgeSubsets(firstage = 50, scaleLM = T))


df_younger$agerange = factor(df_younger$agerange, levels=c("30-65", "30-70", "30-75", "30-80", "30-85",
                                                       "40-65", "40-70", "40-75", "40-80", "40-85",
                                                       "50-65", "50-70", "50-75", "50-80", "50-85"))


#remove redundant age ranges
df_younger %<>% filter(agerange!="50-65", agerange!="30-80", agerange!="30-85", agerange!="50-85",
                       agerange != "40-85", agerange != "50-80",
                       agerange != "40-80")



pal = brewer.pal("Dark2", n=4)
(p1 = df_outolder %>%
    ggplot(., aes(x=agerange, y=bvals_mod1_1,group=factor(colourcode1),color=factor(colourcode1))) +
    scale_x_discrete(position = "top") +
    geom_hline(yintercept = 0, linetype=2, size=1) +
    scale_alpha_discrete(range = c(0, 1)) +
    geom_pointrange(aes(ymin=bvals_mod1_1-civals_mod1_1,ymax=bvals_mod1_1+civals_mod1_1), fill="black",position=position_dodge(width=0.3)) +
    geom_point(aes(alpha=dotalpha1),col="black",position=position_dodge(width=0.3),shape=21, size=2.5, stroke=0.75) +
    labs(x="Age-range", y="Beta") +
    theme_classic() +
    
    scale_color_manual(values=rev(c(pal[4], "light grey"))) +
    scale_fill_manual(values=rev(c(pal[4], "light grey"))) +
    
    mytheme +
    theme(
      axis.text = element_text(size = 8),
      axis.text.x = element_text(angle = 0),
      axis.line.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 8, margin = margin(0,10,0,0)),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      plot.title = element_text(size = 10),
      plot.subtitle = element_text(size = 10, hjust=0.5)
    ) +
    labs(color='PGS') +
    labs(color='PGS') + scale_x_discrete(labels=c('30+', '35+', '40+', '45+','50+','55+','60+'),position = "top")
)


(p2 = df_outolder %>% 
    ggplot(., aes(x=agerange, y=bvals_mod2_1,group=factor(colourcode2),color=factor(colourcode2))) +
    scale_x_discrete(position = "top") +
    geom_hline(yintercept = 0, linetype=2, size=1) +
    scale_alpha_discrete(range = c(0, 1)) +
    geom_pointrange(aes(ymin=bvals_mod2_1-civals_mod2_1,ymax=bvals_mod2_1+civals_mod2_1), fill="black",position=position_dodge(width=0.3)) +
    geom_point(aes(alpha=dotalpha2),col="black",position=position_dodge(width=0.3),shape=21, size=2.5, stroke=0.75) +
    labs(x="Age-range", y="Beta") +
    theme_classic() +
    
    scale_color_manual(values=rev(c(pal[4], "light grey"))) +
    scale_fill_manual(values=rev(c(pal[4], "light grey"))) +
    
    mytheme +
    theme(
      axis.text = element_text(size = 8),
      axis.text.x = element_text(angle = 0),
      axis.line.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 8, margin = margin(0,10,0,0)),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      plot.title = element_text(size = 10),
      plot.subtitle = element_text(size = 10, hjust=0.5)
    ) +
    labs(color='PGS') + scale_x_discrete(labels=c('30+', '35+', '40+', '45+','50+','55+','60+'),position = "top")
)


(p3 = df_younger %>% 
    ggplot(., aes(x=agerange, y=bvals_mod1_1,group=factor(colourcode1),color=factor(colourcode1))) +
    scale_x_discrete(position = "top") +
    geom_hline(yintercept = 0, linetype=2, size=1) +
    scale_alpha_discrete(range = c(0, 1)) +
    geom_pointrange(aes(ymin=bvals_mod1_1-civals_mod1_1,ymax=bvals_mod1_1+civals_mod1_1), fill="black",position=position_dodge(width=0.3)) +
    geom_point(aes(alpha=dotalpha1),col="black",position=position_dodge(width=0.3),shape=21, size=2.5, stroke=0.75) +
    labs(x="Age-range", y="Beta") +
    theme_classic() +
    
    scale_color_manual(values=rev(c(pal[4], "light grey"))) +
    scale_fill_manual(values=rev(c(pal[4], "light grey"))) +
    
    mytheme +
    theme(
      axis.text = element_text(size = 8),
      axis.text.x = element_text(angle = 0),
      axis.line.x = element_blank(),
      # axis.title.x = element_text(size = 8, vjust =-3, margin = margin(0,30,0,0)),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 8, margin = margin(0,10,0,0)),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      plot.title = element_text(size = 10),
      plot.subtitle = element_text(size = 10, hjust=0.5)
    ) +
    labs(color='PGS') 
)


(p4 = df_younger %>% 
    ggplot(., aes(x=agerange, y=bvals_mod2_1,group=factor(colourcode2),color=factor(colourcode2))) +
    scale_x_discrete(position = "top") +
    geom_hline(yintercept = 0, linetype=2, size=1) +
    scale_alpha_discrete(range = c(0, 1)) +
    geom_pointrange(aes(ymin=bvals_mod2_1-civals_mod2_1,ymax=bvals_mod2_1+civals_mod2_1), fill="black",position=position_dodge(width=0.3)) +
    geom_point(aes(alpha=dotalpha2),col="black",position=position_dodge(width=0.3),shape=21, size=2.5, stroke=0.75) +
    labs(x="Age-range", y="Beta") +
    theme_classic() +
    scale_color_manual(values=rev(c(pal[4], "light grey"))) +
    scale_fill_manual(values=rev(c(pal[4], "light grey"))) +
    
    mytheme +
    theme(
      axis.text = element_text(size = 8),
      axis.text.x = element_text(angle = 0),
      axis.line.x = element_blank(),
      # axis.title.x = element_text(size = 8, vjust =-3, margin = margin(0,30,0,0)),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 8, margin = margin(0,10,0,0)),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      plot.title = element_text(size = 10),
      plot.subtitle = element_text(size = 10, hjust=0.5)
    ) +
    labs(color='PGS') 
)



# ggsave(
#   filename = "/ess/p274/cluster/projects/p040-ad_change/ADchangeRisk/reproduce/p1-SIFig17.png",
#   plot = p1,
#   width = 9,
#   height = 6,
#   dpi = 600,
#   units = "cm"
# )
# ggsave(
#   filename = "/ess/p274/cluster/projects/p040-ad_change/ADchangeRisk/reproduce/p2-SIFig17.png",
#   plot = p2,
#   width = 9,
#   height = 6,
#   dpi = 600,
#   units = "cm"
# )
# ggsave(
#   filename = "/ess/p274/cluster/projects/p040-ad_change/ADchangeRisk/reproduce/p3-SIFig17_longer.png",
#   plot = p3 + theme(legend.position = "none"),
#   width = 15, #13
#   height = 6,
#   dpi = 600,
#   units = "cm"
# )
# ggsave(
#   filename = "/ess/p274/cluster/projects/p040-ad_change/ADchangeRisk/reproduce/p4-SIFig17_longer.png",
#   plot = p4 + theme(legend.position = "none"),
#   width = 15,
#   height = 6, #13
#   dpi = 600,
#   units = "cm"
# )




# get better text axis for smaller plots
recalcPlot2 = function(res1, res2, agecut, firstage) {
  
  #recalc distribution plots
  pal = brewer.pal(n=4,"Dark2")
  means1 = res1 %>% group_by(quadrantpresid) %>% summarise(mean=mean(presid1quad))
  means2 = res2 %>% group_by(quadrantpresid) %>% summarise(mean=mean(presid1quad))
  
  
  set.seed(123)
  (p1 = ggplot(res1,
               aes(x=(factor(quadrantpresid)), y=presid1quad, col=factor(quadrantpresid))) +
      geom_jitter(alpha=0.7, position = position_jitter(seed = 123,width = 0.1)) +
      geom_segment(aes(x=0.75, y=means1$mean[1], xend=1.25, yend=means1$mean[1]), size=.75, color="black") +#color=pal[1]) +
      geom_segment(aes(x=1.75, y=means1$mean[2], xend=2.25, yend=means1$mean[2]), size=.75, color="black") +#color=pal[2]) +
      geom_segment(aes(x=2.75, y=means1$mean[3], xend=3.25, yend=means1$mean[3]), size=.75, color="black") +#color=pal[3]) +
      geom_segment(aes(x=3.75, y=means1$mean[4], xend=4.25, yend=means1$mean[4]), size=.75, color="black") +#color=pal[4]) +
      theme_classic() + mytheme +
      scale_color_brewer(palette="Dark2") +
      theme(legend.position = "top",
            axis.ticks.y = element_line(),
            axis.text.y = element_text(size=10)) +
      labs(y="memory change (corrected)",
           x=NULL)
    
  )
  
  (p2 = ggplot(res2,
               aes(x=factor(quadrantpresid), y=presid1quad, col=factor(quadrantpresid))) +
      geom_jitter(alpha=0.7,position = position_jitter(seed = 123,width = 0.1)) +
      geom_segment(aes(x=0.75, y=means2$mean[1], xend=1.25, yend=means2$mean[1]), size=.75, color="black") +#color=pal[1]) +
      geom_segment(aes(x=1.75, y=means2$mean[2], xend=2.25, yend=means2$mean[2]), size=.75, color="black") +#color=pal[2]) +
      geom_segment(aes(x=2.75, y=means2$mean[3], xend=3.25, yend=means2$mean[3]), size=.75, color="black") +#color=pal[3]) +
      geom_segment(aes(x=3.75, y=means2$mean[4], xend=4.25, yend=means2$mean[4]), size=.75, color="black") +#color=pal[4]) +
      theme_classic() + mytheme +
      scale_color_brewer(palette="Dark2") +
      theme(legend.position = "top",
            axis.ticks.y = element_line(),
            axis.text.y = element_text(size=10)) +
      labs(y="memory change (corrected)",
           x=NULL)
  )
  
  Q1res1 = res1[res1$quadrantpresid=="Q1",]
  Q4res1 = res1[res1$quadrantpresid=="Q4",]
  Q2res1 = res1[res1$quadrantpresid=="Q2",]
  Q3res1 = res1[res1$quadrantpresid=="Q3",]
  Q1res1$zero=0
  Q4res1$zero=0
  Q2res1$zero=0
  Q3res1$zero=0
  Q1res2 = res2[res2$quadrantpresid=="Q1",]
  Q4res2 = res2[res2$quadrantpresid=="Q4",]
  Q2res2 = res2[res2$quadrantpresid=="Q2",]
  Q3res2 = res2[res2$quadrantpresid=="Q3",]
  Q1res2$zero=0
  Q4res2$zero=0
  Q2res2$zero=0
  Q3res2$zero=0
  
  
  #NB! to align plots q3 and q4 must be transparent
  (dens3 <- ggplot() +
      PupillometryR::geom_flat_violin(data=Q1res1,aes(x=zero,y=presid1quad),alpha = 0.3, col = pal[1], fill = pal[1]) +
      PupillometryR::geom_flat_violin(data=Q4res1,aes(x=zero,y=presid1quad),alpha = 0.3, col = pal[4], fill = pal[4]) +
      PupillometryR::geom_flat_violin(data=Q2res1,aes(x=zero,y=presid1quad),alpha = 0.0, col = NA, fill = NA) + #col = pal[2]
      PupillometryR::geom_flat_violin(data=Q3res1,aes(x=zero,y=presid1quad),alpha = 0.0, col = NA, fill = NA) + #col = pal[3]
      theme_void())
  
  (dens4 <- ggplot() +
      PupillometryR::geom_flat_violin(data=Q1res2,aes(x=zero,y=presid1quad),alpha = 0.3, col = pal[1], fill = pal[1]) +
      PupillometryR::geom_flat_violin(data=Q4res2,aes(x=zero,y=presid1quad),alpha = 0.3, col = pal[4], fill = pal[4]) +
      PupillometryR::geom_flat_violin(data=Q2res2,aes(x=zero,y=presid1quad),alpha = 0.0, col = NA, fill = NA) + #col = pal[2]
      PupillometryR::geom_flat_violin(data=Q3res2,aes(x=zero,y=presid1quad),alpha = 0.0, col = NA, fill = NA) + #col = pal[3]
      theme_void())
  
  
  (pAll1 = p1 + dens3 +
      patchwork::plot_layout(ncol = 2, nrow = 1, widths = c(2, 0.75), heights = c(1, 4))
  )
  (pAll2 = p2 + dens4 +
      patchwork::plot_layout(ncol = 2, nrow = 1, widths = c(2, 0.75), heights = c(1, 4))
  )
  
  #this iis the only difference
  test1 = memSlopesCompare(x = BOTH, filter_by = meanAge, max_age = agecut, min_age = firstage, yvar = "rSlopeMemory", type = "younger", scale = T)[[2]]
  test2 = memSlopesCompare(x = BOTH, filter_by = meanAge, max_age = agecut, min_age = firstage, yvar = "absChangeMemory", type = "younger", scale = T)[[2]]
  # test1=tmpfunc3(res1, meanAge, agecut, "rSlopeMemory", firstage)[[2]]
  # test2=tmpfunc3(res2, meanAge, agecut, "absChangeMemory", firstage)[[2]]
  # test1=tmpfunc2(res1, meanAge, agecut, "rSlopeMemory", 89)[[2]]
  # test2=tmpfunc2(res2, meanAge, agecut, "absChangeMemory", 89)[[2]]
  
  return(list(pAll1 = pAll1, 
              pAll2 = pAll2,
              test1 = test1,
              test2 = test2,
              minage=min(res1$meanAge),
              maxage=max(res1$meanAge)))
}




firstage=30
lastage=89
runplot=1
iii=0
pOut = list()
if (runplot == 1) {
  for (agecut in seq(65, 80, by=5)) {
    iii = iii+1
    res1 = memSlopesCompare(x = BOTH, filter_by = meanAge, max_age = agecut, min_age = firstage, yvar = "rSlopeMemory", type = "younger", scale = T)[[1]]
    res2 = memSlopesCompare(x = BOTH, filter_by = meanAge, max_age = agecut, min_age = firstage, yvar = "absChangeMemory", type = "younger", scale = T)[[1]]
    res1$newvar=as.numeric(substr(res1$quadrantpresid,2,2))
    res2$newvar=as.numeric(substr(res2$quadrantpresid,2,2))
    res1 %<>% mutate(quadrantpresid = factor(quadrantpresid, levels=c("Q1", "Q2", "Q3", "Q4")))
    res2 %<>% mutate(quadrantpresid = factor(quadrantpresid, levels=c("Q1", "Q2", "Q3", "Q4")))
    
    pOut[[iii]]=recalcPlot2(res1,res2, agecut, firstage)
    
    # filename1=file.path(
    #   paste0("/cluster/projects/p274/projects/p040-ad_change/ADchangeRisk/reproduce/memChangeRange/p",firstage,"-", agecut, ".rSlope_smallertext_Z")
    # )
    # filename2=file.path(
    #   paste0("/cluster/projects/p274/projects/p040-ad_change/ADchangeRisk/reproduce/memChangeRange/p",firstage,"-", agecut, ".full_smallertext_Z")
    # )
    
    # ggsave(plot=pOut[[iii]]$pAll1,
    #        file = paste0(filename1,".png"),
    #        width=4.75, height=6, units="cm", dpi=600)
    # ggsave(plot=pOut[[iii]]$pAll2,
    #        file = paste0(filename2,".png"),
    #        width=4.75, height=6, units="cm", dpi=600)
  }
}

#prints plots, lm results, tracked ages 
pOut[[1]]
pOut[[2]]
pOut[[3]]
pOut[[4]]


#age distributions
BOTH %>% 
  filter(meanAge > 30, meanAge <90) %>% 
  ggplot(., aes(x=meanAge, y=absChangeMemory)) + geom_point(aes(col=quadrantpresid))
BOTH %>% 
  filter(meanAge > 30, meanAge <90) %>% 
  ggplot(., aes(x=meanAge)) + geom_density(aes(col=quadrantpresid))
