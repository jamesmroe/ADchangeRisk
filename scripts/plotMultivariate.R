#========================================================================================#
# Written by James M Roe, Ph.D.
# Center for Lifespan Changes in Brain and Cognition, Department of Psychology
# University of Oslo, Norway
#========================================================================================#


#========================================================================================#
## Purpose: plot multivariate PRS-AD model results for main (figure 4) or replication sample (figure 5)
## Script is fully executable
## Instructions: set to analysis_sample to "main" or "replication" to plot the multivariate results for that sample
#========================================================================================#


#---load packages----
tmp.packages = c("dplyr", "stringr", "magrittr","here","sgof")
tmpnew.packages = tmp.packages[!(tmp.packages %in% installed.packages()[,"Package"])]
if(length(tmpnew.packages)) {
  install.packages(tmpnew.packages)
}
sapply(tmp.packages, require, character.only = T)
rm(list=ls())
here()


# set to "main" or "replication"
analysis_sample = "main" #"replication"
# analysis_sample = "replication"


load_from = "excel"
# load_from = "Rdat"


if (analysis_sample == "main") {
  sheets = c("Fig_4b-c PRS-AD",
           "Fig_4b-c noAPOE",
           "Fig_4e PRS-AD",
           "Fig_4e noAPOE")
} else {
  sheets = c("Fig_5c-d replication")
}


if (load_from == "excel") {
  
  if (analysis_sample == "main") {
    OUT1_PC1 = read.xlsx("source_data.xlsx", sheet = sheets[1])
    OUT2_PC1 = read.xlsx("source_data.xlsx", sheet = sheets[2])
    OUT1_PCroll = read.xlsx("source_data.xlsx", sheet = sheets[3])
    OUT2_PCroll = read.xlsx("source_data.xlsx", sheet = sheets[4])
    
    # significance indicators for plotting
    OUT1_PC1$FDRsig = OUT1_PC1$FDRsig_PRSAD
    OUT1_PCroll$FDRsig = OUT1_PCroll$FDRsig_PRSAD
    OUT1_PC1$dotalpha = factor(ifelse(OUT1_PC1$p.value < .05, 1, 0))
    OUT2_PC1$dotalphapoe = factor(ifelse(OUT2_PC1$p.value < .05, 1, 0))
    OUT1_PCroll$dotalpha = factor(ifelse(OUT1_PCroll$p.value < .05, 1, 0))
    OUT2_PCroll$dotalphapoe = factor(ifelse(OUT2_PCroll$p.value < .05, 1, 0))
    OUT1_PC1$agecut = OUT1_PC1$lower_agerange_genetic
    OUT1_PCroll$ww = OUT1_PCroll$window
    OUT2_PCroll$ww = OUT2_PCroll$window
    
  } else if (analysis_sample == "replication") {
    OUT1_PC1 = read.xlsx("source_data.xlsx", sheet = sheets[1])
    OUT1_PC1$dotalpha = factor(ifelse(OUT1_PC1$p.value < .05, 1, 0))
  }
  
} else if (load_from == "Rdat") {
  
  if (analysis_sample == "main") {
    load(here("results/PRS-ADmodels_figure4_PCmultivariate.Rda"))
    OUT1_PC1 = figure4B_PRSAD
    OUT2_PC1 = figure4B_PRSADnoAPOE
    OUT1_PCroll = figure4D_PRSAD
    OUT2_PCroll = figure4D_PRSADnoAPOE
    
    # significance indicators for plotting
    OUT1_PC1$FDRsig_PRSAD = OUT1_PC1$FDRsig
    OUT1_PCroll$FDRsig_PRSAD = OUT1_PCroll$FDRsig
    OUT2_PC1$dotalphapoe = factor(OUT2_PC1$dotalphapoe)
    OUT2_PCroll$dotalphapoe = factor(OUT2_PCroll$dotalphapoe)
    
  } else if (analysis_sample == "replication") {
    load(here("results/PRS-ADmodels_figure5_PCmultivariate_replication.Rda"))
    OUT1_PC1 = figure5C_PRSAD
  }
}

if (analysis_sample == "main") {
  #confirm FDR-corrected results
  FDR_wAPOE = BH(c(OUT1_PC1$p.value, OUT1_PCroll$p.value), .05)
  (FDRthreshAPOE=max(FDR_wAPOE$data[FDR_wAPOE$Adjusted.pvalues<.05]))
  ifelse(OUT1_PC1$p.value <= FDRthreshAPOE, 1, 0) == OUT1_PC1$FDRsig_PRSAD
  ifelse(OUT1_PCroll$p.value <= FDRthreshAPOE, 1, 0) == OUT1_PCroll$FDRsig_PRSAD
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



# ---- plotALL results ----
plotALL = 1
if (plotALL == 1) {
  
  # PC1 associations ----
  if (analysis_sample == "main") {
    
    (pPCA = ggplot(OUT1_PC1, aes(x=agerange, y=-log10(p.value), group=factor(PGS), shape=factor(PGS))) +
       scale_x_discrete() +
       scale_alpha_discrete(range = c(0,1)) +
       geom_hline(yintercept = -log10(0.05), linetype=2) +
       geom_hline(yintercept = -log10(FDRthreshAPOE), linetype=3) +
       geom_point(position = position_dodge(width=0.5),colour="#793349") +
       geom_point(aes(alpha=factor(dotalpha)), col="black", position=position_dodge(width=0.5), shape=21, size=4, stroke=0.75) +
       theme_classic() + 
       labs(x = "Age-range") +
       ggtitle(paste0("PCA 1-50noHippnoAmy")) +
       mytheme +
       theme(axis.text = element_text(size=8),
             axis.title.x = element_text(size=7, vjust=-3, margin=margin(0,30,0,0)),
             axis.title.y = element_text(size=7, margin = margin(0,10,0,0)),
             legend.title = element_text(size=8),
             plot.title = element_text(size=10),
             plot.subtitle = element_text(size=10, hjust=0.5),
             axis.ticks = element_line()
       ) +
       
       #change shape
       scale_shape_manual(values = c("AD_Jansen"=16, "AD_Kunkle"=15, "AD_Lambert"=17, "AD_Wightman"=4))
    )
    
  } else if (analysis_sample == "replication") {
   
     (pPCA = ggplot(OUT1_PC1, aes(x=agerange, y=-log10(p.value), group=factor(PGS), shape=factor(PGS))) +
       scale_x_discrete() +
       scale_alpha_discrete(range = c(0,1)) +
       geom_hline(yintercept = -log10(0.05), linetype=2) +
       # geom_hline(yintercept = -log10(FDRthreshAPOE), linetype=3) +
       geom_point(position = position_dodge(width=0.5),colour="#793349") +
       geom_point(aes(alpha=factor(dotalpha)), col="black", position=position_dodge(width=0.5), shape=21, size=4, stroke=0.75) +
       theme_classic() + 
       labs(x = "Age-range") +
       ggtitle(paste0("PCA 1-50noHippnoAmy")) +
       mytheme +
       theme(axis.text = element_text(size=8),
             axis.title.x = element_text(size=7, vjust=-3, margin=margin(0,30,0,0)),
             axis.title.y = element_text(size=7, margin = margin(0,10,0,0)),
             legend.title = element_text(size=8),
             plot.title = element_text(size=10),
             plot.subtitle = element_text(size=10, hjust=0.5),
             axis.ticks = element_line()
       ) +
       
       #change shape
       scale_shape_manual(values = c("AD_Jansen"=16, "AD_Kunkle"=15, "AD_Lambert"=17, "AD_Wightman"=4))
    ) 
  }
  
  
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
           plot.subtitle = element_text(size=10, hjust=0.5),
           axis.ticks = element_line()
     ) +
     #change shape
     scale_shape_manual(values = c("AD_Jansen"=16, "AD_Kunkle"=15, "AD_Lambert"=17, "AD_Wightman"=4)) +
     theme(axis.text = element_text(size=14),
           axis.text.x = element_blank(),
           axis.ticks.x = element_line(),
           panel.border = element_rect(color = "black", fill = NA, size = 0.75))
  )
  
  
  if (analysis_sample == "main") {
    maxscalersq = 0.1844 #scale of rsquared plot
    
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
         axis.ticks = element_line(),
         
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
             plot.subtitle = element_text(size=10, hjust=0.5),
             axis.ticks = element_line()
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
         axis.ticks = element_line()
         
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
}

