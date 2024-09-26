#========================================================================================#
# Written by James M Roe, Ph.D.
# Center for Lifespan Changes in Brain and Cognition, Department of Psychology
# University of Oslo, Norway
#========================================================================================#


#========================================================================================#
## Purpose: plot PRS-AD model results
## Instructions:
## Set simulated = 0 to load the provided summary-level source data underlying Fig. 1E-F and Fig. 2 (PRS-AD association tests),
##                   run multiple testing correction, and reproduce plots
## Set simulated = 1 to visualize the simulated PRS-AD association tests from the simulated GAMM analysis
##
## Script is fully executable
#========================================================================================#


rm(list=ls())


# select analysis stream -----------------
simulated = 0
saveplots = 0
load_from = "excel" #/ "Rdat"
# load_from = "Rdat"
# select analysis stream -----------------


# resdir="/cluster/p/p274/cluster/projects/p040-ad_change/ADchangeRisk/reproduce/results"
# setwd(resdir)
library("here")


# load packages ----
loadPackages = function() {
  packages = c("magrittr","tidyverse", "sgof","MetBrewer","openxlsx")
  new.packages = packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) {
    install.packages(new.packages)
  }
  sapply(packages, require, character.only = T)
}
loadPackages()


# ROIs ----
roinames=c("Left-Hippocampus.fs71", #L hippo (stage II)
           "Right-Hippocampus.fs71", #R hippo (stage (II)
           "braak-stage-1ent-.volume-lh-71", #L entorhinal (stage I)
           "braak-stage-1ent-.volume-rh-71", #R entorhinal (stage I)
           "braak-stage-3-amy-.volume-lh-71", #cortical stage III ROI excluding amygdala L
           "braak-stage-3-amy-.volume-rh-71", #cortical stage III ROI excluding amygdala R
           "braak-stage-3-amygdala-.volume-lh-71", #L amygdala (stage II)
           "braak-stage-3-amygdala-.volume-rh-71") #R amygdala (stage III)

nTime=2


if (simulated == 1) {
  #change loop to be 1:1
  roinames="Left-Hippocampus.fs71"
  roiname=roinames
  roi_sheets = length(roinames)
} else {
  roi_sheets = seq(1, length(roinames), by = 2)
}


if (load_from == "excel") {

  # load and merge results from excel file ----
  for (sheet in roi_sheets) {
    if (sheet==1) {
      SLOPESALL1=c()
      SLOPESALL2=c()
      SLOPESALL3=c()
      SLOPESALL4=c()
    }
    
    if (simulated) {
      print("loading simulated results")
      PRSAD <- read.xlsx("simulated.xlsx", sheet = sheet)
      PRSAD_NOAPOE <- read.xlsx("simulated.xlsx", sheet = sheet+1)
      
    } else if (!simulated) {
      print("loading published results")
      PRSAD <- read.xlsx("source_data.xlsx", sheet = sheet)
      PRSAD_NOAPOE <- read.xlsx("source_data.xlsx", sheet = sheet+1)
    }
  
    SLOPESALL1 %<>% rbind(., PRSAD[PRSAD$change == "ageRelChange",]) #AGERELSLOPE (PRS-AD age-relative change)
    SLOPESALL2 %<>% rbind(., PRSAD_NOAPOE[PRSAD_NOAPOE$change == "ageRelChange",]) #AGERELSLOPENOAPOE (PRS-ADnoAPOE age-relative change)
    SLOPESALL3 %<>% rbind(., PRSAD[PRSAD$change == "absChange",]) #ABSCHANGESLOPE (PRS-AD absolute change)
    SLOPESALL4 %<>% rbind(., PRSAD_NOAPOE[PRSAD_NOAPOE$change == "absChange",]) #ABSCHANGESLOPENOAPOE (PRS-ADnoAPOE absolute change)
    
  }
  
} else if (load_from == "Rdat") {
  
  # load and merge results from Rdata files ----
  for (roiname in roinames) {
    if (roiname==roinames[1]) {
      SLOPESALL1=c()
      SLOPESALL2=c()
      SLOPESALL3=c()
      SLOPESALL4=c()
    }
    
    if (simulated) {
      print("loading simulated results")
      load(here("simulate/PRSADmodels_simulated.Rda"))
    
    } else if (!simulated) {
      #main results
      print("loading published results")
      filename = paste0(here(paste0("results/PRSADmodels_S1_",roiname,".publication.Rda")))
      print(paste0("loading ",filename))
      load(filename)  
    }
    
    OUTMERGE=list()
    OUTMERGE[[1]]=OUT[[1]] #AGERELSLOPE (PRS-AD age-relative change)
    OUTMERGE[[2]]=OUT[[2]] #AGERELSLOPENOAPOE (PRS-ADnoAPOE age-relative change)
    OUTMERGE[[3]]=OUT[[3]] #ABSCHANGESLOPE (PRS-AD absolute change)
    OUTMERGE[[4]]=OUT[[4]] #ABSCHANGESLOPENOAPOE (PRS-ADnoAPOE absolute change)
    
  
    SLOPES1 = OUTMERGE[[1]] %>% filter(grepl("PGS",term))
    SLOPES2 = OUTMERGE[[2]] %>% filter(grepl("PGS",term))
    SLOPES3 = OUTMERGE[[3]] %>% filter(grepl("PGS",term))
    SLOPES4 = OUTMERGE[[4]] %>% filter(grepl("PGS",term))
    
    
    SLOPESALL1 %<>% rbind(., SLOPES1)
    SLOPESALL2 %<>% rbind(., SLOPES2)
    SLOPESALL3 %<>% rbind(., SLOPES3)
    SLOPESALL4 %<>% rbind(., SLOPES4)
    
  }
  
  rename_cols = function(dat) {
    dat %<>% rename(N_genetic = N, N_lifespan = Nlife,	lower_agerange_genetic = agecut2, lower_agerange_life = agecut)
  }
  
  SLOPESALL1 = rename_cols(SLOPESALL1)
  SLOPESALL2 = rename_cols(SLOPESALL2)
  SLOPESALL3 = rename_cols(SLOPESALL3)
  SLOPESALL4 = rename_cols(SLOPESALL4)
  
  SLOPESALL1 %<>% mutate(change = "ageRelChange")
  SLOPESALL2 %<>% mutate(change = "ageRelChange")
  SLOPESALL3 %<>% mutate(change = "absChange")
  SLOPESALL4 %<>% mutate(change = "absChange")
}


#DF of all 576 PRS-AD tests
allPRSAD = rbind(SLOPESALL1,
              SLOPESALL3)


#144 FDR-corrected PRS-AD tests
(ALLFDR=BH(allPRSAD$p.value, alpha=0.05))
FDRthresh = max(ALLFDR$data[ALLFDR$Adjusted.pvalues<.05])


if (!simulated & load_from == "excel") {
  # confirm FDR-corrected results
  ifelse(allPRSAD$p.value <=FDRthresh, 1, 0) == allPRSAD$FDRsig_PRSAD
}


#add FDR indicator to DF for plotting
SLOPESALL1$FDRsig=0
SLOPESALL1$FDRsig[SLOPESALL1$p.value<=FDRthresh]=1
SLOPESALL3$FDRsig=0
SLOPESALL3$FDRsig[SLOPESALL3$p.value<=FDRthresh]=1

#all FDR-corrected associations with random slopes are negaative
(sum(SLOPESALL1$FDRsig==1 & SLOPESALL1$estimate<0) )/ sum(SLOPESALL1$FDRsig==1)

#all FDR-corrected associations with abschange are negaative
(sum(SLOPESALL3$FDRsig==1 & SLOPESALL3$estimate<0) )/ sum(SLOPESALL3$FDRsig==1)

#add into noAPOE DF's for plotting
SLOPESALL2$FDRsig = SLOPESALL1$FDRsig
SLOPESALL4$FDRsig = SLOPESALL3$FDRsig


# agerelchange improvement ----
if (!simulated) {
  
  percIncr = function(b1, b2) {
    (t = t.test(b1, b2))
    (t_paired = t.test(b1, b2, paired = T))
    (inc = ((mean(b1) - mean(b2)) / mean(b2)) * 100)

    return(
      list(t = t, t_paired, inc = inc)
      )
  }
  beta_agerel = SLOPESALL1$estimate[SLOPESALL1$change == "ageRelChange"]
  beta_abs = SLOPESALL3$estimate[SLOPESALL3$change == "absChange"]
  
  beta_agerel_bothFDR = SLOPESALL1$estimate[SLOPESALL1$FDRsig == 1 & SLOPESALL3$FDRsig == 1]
  beta_abs_bothFDR = SLOPESALL3$estimate[SLOPESALL1$FDRsig == 1 & SLOPESALL3$FDRsig == 1]
  
  mean(beta_agerel); mean(beta_abs)
  mean(beta_agerel_bothFDR); mean(beta_abs_bothFDR)
  
  print("Across the 576 PRS-AD tests, conditioning brain change estimates on age conferred a 91% average strengthening in the observed PRS-AD regression coefficients, compared to using absolute change estimates and correcting for mean age (t = -26.9, p = 2.0e-80; Figures 1-2)")
  (beta_percIncr = percIncr(beta_agerel, beta_abs))
  print(beta_percIncr)
  beta_percIncr[[2]]$p.value
  
  print("This was also evident in the subset of results where both change estimates were FDR-corrected significant (54%, t = -20.6, p = 1.1e-28)")
  (beta_FDR_percIncr = percIncr(beta_agerel_bothFDR, beta_abs_bothFDR))
  print(beta_FDR_percIncr)
  beta_FDR_percIncr[[2]]$p.value
}


#add better ROI names
SLOPESALL1$roiname = SLOPESALL2$roiname = SLOPESALL3$roiname = SLOPESALL4$roiname = 
  gsub("\\)", "", gsub("\\(", ".", gsub(" ", "-", SLOPESALL1$roi)))
  

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


SLOPESALL1$ci = SLOPESALL1$std.error*1.96


#plot ROI results ----
if (simulated == 1) {
  nrois=1
} else if (simulated == 0) {
  nrois=8
}

pSlopeAllfacet = list()
for (roinum in c(1:nrois)) {
  
  plotProc=1
  if (plotProc == 1) {
    
    
    (roiname=roinames[roinum])
    SLOPES1 = SLOPESALL1 %>% filter(roiname %in%  roinames[roinum]) #%in% roiname
    SLOPES2 = SLOPESALL2 %>% filter(roiname %in%  roinames[roinum])
    SLOPES3 = SLOPESALL3 %>% filter(roiname %in%  roinames[roinum])
    SLOPES4 = SLOPESALL4 %>% filter(roiname %in%  roinames[roinum])
    
    
    #make agecut plot
    addCols = function(dat) {
      dat %<>% mutate(
                      agerange = paste0(lower_agerange_genetic,"-89"),
                      agerange = factor(agerange, levels = rev(c("70-89","65-89","60-89","55-89","50-89","45-89","40-89","35-89","30-89"))),
                      sig = ifelse(p.value <= .05, "sig", ifelse(p.value <= .1 & p.value >.05, "trend", "nonsig")),
                      ci = std.error*1.96,
                      dotalpha = ifelse(sig == "sig",1,0),
                      dotalpha = factor(dotalpha, levels=c(0, 1)),
                      PRS = paste0("AD_",score))
      
      
      #rename lambert scores
      dat$PRS[dat$PRS=="AD_AD"] = "AD_Lambert"
      
      
      #colour coding
      dat %<>% mutate(colourcode = case_when(
        score == "Jansen" & dotalpha==1   ~ "Jansen_sig",
        score == "Jansen" & dotalpha!=1   ~ "Jansen_nosig",
        score == "Kunkle" & dotalpha==1   ~ "Kunkle_sig",
        score == "Kunkle" & dotalpha!=1   ~ "Kunkle_nosig",
        score == "AD"     & dotalpha==1   ~ "Lambert_sig",
        score == "AD"     & dotalpha!=1   ~ "Lambert_nosig",
        score == "Wightman" & dotalpha==1   ~ "Wightman_sig",
        score == "Wightman" & dotalpha!=1   ~ "Wightman_nosig",
      ))
      
      return(dat)
    }
    
    #colour mapping
    colours <- c('Jansen_sig' = '#eec96d',
                 'Jansen_nosig'  = 'light grey',
                 'Kunkle_sig'   = '#96c683',
                 'Kunkle_nosig'   = 'light grey',
                 'Lambert_sig' = '#7f92e2',
                 'Lambert_nosig'   = 'light grey',
                 'Wightman_sig' = '#5b68a7',
                 'Wightman_nosig'   = 'light grey'
    )
    
    
    SLOPES1 = addCols(SLOPES1)
    SLOPES2 = addCols(SLOPES2)
    SLOPES3 = addCols(SLOPES3)
    SLOPES4 = addCols(SLOPES4)
    
    SLOPES1 %<>% mutate(type="rSlope")
    SLOPES2 %<>% mutate(type="rSlopeNoApoe")
    SLOPES3 %<>% mutate(type="absChange")
    SLOPES4 %<>% mutate(type="absChangeNoApoe")
    

    #report reduction in negative beta coeficient for each increasing age subset
    SLOPES_BOTH = rbind(SLOPES1, SLOPES3)
    SLOPES_BOTH %<>% mutate(age_recode = case_when(
      lower_agerange_genetic == 30 ~ 0,
      lower_agerange_genetic == 35 ~ 1,
      lower_agerange_genetic == 40 ~ 2,
      lower_agerange_genetic == 45 ~ 3,
      lower_agerange_genetic == 50 ~ 4,
      lower_agerange_genetic == 55 ~ 5,
      lower_agerange_genetic == 60 ~ 6,
      lower_agerange_genetic == 65 ~ 7,
      lower_agerange_genetic == 70 ~ 8
      
    ))
    aov1 = summary(aov(estimate ~ age_recode*type, data = SLOPES_BOTH))
    checkaov = summary(lm(estimate ~ age_recode*type, data = SLOPES_BOTH))
    checkaov$coefficients
    #plot(SLOPES_BOTH$lower_agerange_genetic,SLOPES_BOTH$estimate); abline(lm(SLOPES_BOTH$estimate ~ SLOPES_BOTH$lower_agerange_genetic))
    #ggplot(SLOPES_BOTH, aes(x=lower_agerange_genetic, y=estimate, col=type)) + geom_point()
    
    
    #plot in same space
    maxscale = max(SLOPESALL1$ci) + .02
    minscale = min(SLOPESALL1$ci) - .02
    maxscalersq = max(SLOPESALL1$prsq) + SLOPESALL1$prsq_upr[SLOPESALL1$estimate==max(SLOPESALL1$estimate)]+.05
    
    
    # where PRS-AD association survived FDR-correction indicate if PRS-ADnoAPOE association was also signifcant (rSlopes)
    SLOPES1$dotalphapoe = as.factor(
      ifelse(
        SLOPES1$FDRsig == 1 &
          SLOPES2$dotalpha == 1 &
          sign(SLOPES1$estimate) == -1 & sign(SLOPES2$estimate) == -1,
        1,
        0
      )
    )
    #add into other df
    SLOPES2$dotalphapoe=SLOPES1$dotalphapoe
    
    
    # where PRS-AD association survived FDR-correction indicate if PRS-ADnoAPOE association was also signifcant (absChange)
    SLOPES3$dotalphapoe = as.factor(
      ifelse(
        SLOPES3$FDRsig == 1 &
          SLOPES4$dotalpha == 1 &
          sign(SLOPES3$estimate) == -1 & sign(SLOPES4$estimate) == -1,
        1,
        0
      )
    )
    #add into other df
    SLOPES4$dotalphapoe = SLOPES3$dotalphapoe
    
    
    #combine PRS-AD associations for facet plot (rSlopes and absChange)
    SLOPESBOTHAP = rbind(SLOPES1 %>% select(names(SLOPES1)[names(SLOPES1) %in% names(SLOPES3)] ),
                       SLOPES3 )
    SLOPESBOTHAP$type = factor(SLOPESBOTHAP$type, levels = c("rSlope","absChange"))
    
    
    #combine PRS-ADnoAPOE associations for facet plot (rSlopes and absChange)
    SLOPESBOTHnoAP = rbind(SLOPES2 %>% select(names(SLOPES2)[names(SLOPES2) %in% names(SLOPES4)] ),
                         SLOPES4 )
    SLOPESBOTHnoAP$type = factor(SLOPESBOTHnoAP$type, levels = c("rSlopeNoApoe","absChangeNoApoe"))
    
    
  
    
    (pSlopeFacet = SLOPESBOTHAP %>%
        
        ggplot(., aes(x=agerange, y=estimate,group=factor(colourcode),color=factor(colourcode))) +
        scale_x_discrete(position = "top") +
        geom_hline(yintercept = 0, linetype=2, size=1) +
        scale_alpha_discrete(range = c(0, 1)) +
        # ylim(c(minscale[1], maxscale[1])) +
        
      geom_pointrange(aes(ymin=estimate-ci,ymax=estimate+ci), fill="black",position=position_dodge(width=0.3)) +
        geom_point(aes(alpha=dotalpha),col="black",position=position_dodge(width=0.3),shape=21, size=2.5, stroke=0.75) +
        labs(x="Age-range", y="Beta") +
        theme_classic() +
        ggtitle(roiname,
                subtitle = "rSlopes / absChange") +
        scale_color_manual(values=colours) +
        scale_fill_manual(values=colours) +
        
        mytheme +
        theme(
          axis.text = element_text(size = 8),
          axis.line.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 8, margin = margin(0,10,0,0)),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          plot.title = element_text(size = 10),
          plot.subtitle = element_text(size = 10, hjust=0.5)
        ) +
        labs(color='PRS') + facet_grid(. ~ type) 
      
    )
    
    
    #alpha mapping will not work if all significant
    if (! 0 %in% unique(SLOPESBOTHAP$dotalpha)) {
      
      (pSlopeFacet = pSlopeFacet + 
         geom_point(aes(),alpha=1,col="black",position=position_dodge(width=0.3),shape=21, size=2.5, stroke=0.75))
      
    }
    
    
    
    #rsq must be from SLOPES1
    #alpha and rsq no APOE from others
    SLOPESBOTHAP$prsqnoAP = SLOPESBOTHnoAP$prsq
    
    (pSlopeRsqFacet = SLOPESBOTHAP %>%
        
        ggplot(., aes(x=agerange, y=prsq*-1,group=factor(PRS),color=factor(PRS),fill=factor(PRS))) +
        scale_x_discrete() +
        scale_alpha_discrete(range = c(0, 1)) +
        ylim(c(maxscalersq*-1,0)) +
       
      geom_col(aes(alpha=as.factor(FDRsig), y=prsq*-1,x=agerange), color=NA, position=position_dodge(width=0.5), width=0.4) +
        geom_errorbar(aes(alpha=as.factor(FDRsig), ymin=prsq_lwr*-1,ymax=prsq_upr*-1),position=position_dodge(width=0.5),size=0.3, width=0) +
        geom_point(aes(alpha=dotalphapoe, y=prsqnoAP*-1, x=agerange),col="black",position=position_dodge(width=0.5),shape=3, size=2, stroke=1) +
        geom_hline(yintercept = 0, size=0.5) +
        
        labs(x="Age-range", y=paste(parse(text="R^2"),"PRS")) +
        theme_classic() +
        scale_color_manual(values=met.brewer("Derain", 4)) + #colours match colourcode mapping in pSlopeFacet
        scale_fill_manual(values=met.brewer("Derain", 4)) + #colours match colourcode mapping in pSlopeFacet
        
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
        labs(color='PRS',fill='PRS', alpha="p<.05", shape=FALSE) + facet_grid(. ~type)
    ) 
    
    
    (pSlopeFacetNoLegend = pSlopeFacet + theme(legend.position = "none",
                                              strip.background = element_blank(),
                                              strip.text.x = element_blank()) + scale_x_discrete(labels=c('30+', '35+', '40+', '45+','50+','55+','60+','65+','70+'),
                                                                                                 position = "top") +
      theme(axis.ticks.y = element_line(),
            axis.text = element_text(size=9.5),
            # plot.title = element_text(colour="transparent"),
            # plot.subtitle = element_text(colour="transparent")
      ))
    
    (pSlopeRsqFacetNoLegend = pSlopeRsqFacet + theme(legend.position = "none",
                                                    strip.background = element_blank(),
                                                    strip.text.x = element_blank()) + scale_x_discrete(labels=c('30+', '35+', '40+', '45+','50+','55+','60+','65+','70+')) +
    theme(axis.ticks.y = element_line(),
          axis.text = element_text(size=9.5),
          # plot.title = element_text(colour="transparent"),
          # plot.subtitle = element_text(colour="transparent")
    ))
    
    (pSlopeAllfacet[[roinum]] = cowplot::plot_grid(pSlopeFacetNoLegend + scale_y_continuous(labels = scales::number_format(accuracy = 0.01)), 
                                         pSlopeRsqFacetNoLegend, nrow=2,rel_heights = c(2, 1)))
    
    
    #SAVE PLOTS ----
    if (simulated == 0) {
      
      if (saveplots == 1) {
        plotdir=file.path(here("plots/BetaAgeRange"))
        if (!dir.exists(plotdir)) { dir.create(plotdir)}
        
        # main facet plot
        ggsave(
          filename = file.path(plotdir, paste0(
            "pBetaAgeRange.facet.rsq.", roiname, ".-nTime",nTime,"_text9.5.png"
          )),
          plot = pSlopeAllfacet[[roinum]],
          width = 26,#13
          bg = "transparent",
          height = 20,#10
          dpi = 600,
          units = "cm",
          device = "png"
        )
      }
    }
    
  }
}

# view final plots ----
if (simulated) {
  pSlopeAllfacet[[1]] + ggtitle("SIMULATED")
} else {
  pSlopeAllfacet[[1]]
  pSlopeAllfacet[[2]]
  pSlopeAllfacet[[3]]
  pSlopeAllfacet[[4]]
  pSlopeAllfacet[[5]]
  pSlopeAllfacet[[6]]
  pSlopeAllfacet[[7]]
  pSlopeAllfacet[[8]]
}
