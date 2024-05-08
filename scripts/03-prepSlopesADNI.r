#========================================================================================#
# Author: James M Roe, Ph.D.
# Center for Lifespan Changes in Brain and Cognition, University of Oslo
#
# Purpose: prepare individual-specific slopes in longitudinal AD-control data (ADNI) for ML
#          Script requires individual-level data as input and is not executable
#========================================================================================#


#============INPUTS=================
args = commandArgs(TRUE)
jj=as.integer(args[1])
subset.size = as.integer(args[2])
print(paste("split =",jj))
print(paste("subset =",subset.size))
#==================================#


packages <- c("dplyr","tidyverse", "magrittr","gamm4","itsadug","gratia","here")
sapply(packages, require, character.only = T)


# rm(list=ls())
a="/ess/p23/cluster/workspace/projects/11_MemP/James/ADchangeRisk"
DF = read.csv(file.path(a, "ADNI_export.tsv"), sep="\t", header = T, stringsAsFactors = F)
savefigs = 0
U = DF %>% filter(timeBsl == 0)
dim(DF)
head(DF)


#plot longitudinal trajectories
sub_trajectories = DF %>% drop_na(DX) %>% filter(DX!="Na")
sub_trajectories$group = ifelse(sub_trajectories$includeLong=="D", "AD_long", "NC_long")
which(is.na(sub_trajectories$DX))
unique(sub_trajectories$DX)
sub_trajectories$DXtmp = as.character(sub_trajectories$DX)
sub_trajectories$DXtmp = ifelse(sub_trajectories$DXtmp == "CN","NC",
                   ifelse(sub_trajectories$DXtmp == "Dementia", "AD",
                          ifelse(sub_trajectories$DXtmp == "MCI", "MCI",NA)))
unique(sub_trajectories$DXtmp)
sub_trajectories$DXtmp = factor(sub_trajectories$DXtmp, levels=rev(c("NC", "MCI", "AD")))
sub_trajectoriesADlong = sub_trajectories %>% filter(includeLong == "D")
sub_trajectoriesNClong = sub_trajectories %>% filter(includeLong == "NC")


(traj=ggplot() +
  geom_line(data=sub_trajectories,aes(x=TP,y=DXtmp,group=RID,col=group),position=position_jitter(w=0.15, h=0.15),alpha=0.2) +
  scale_x_continuous(limits = c(0, max(sub_trajectories$TP)+1), breaks = 0:50) +
  labs(x = "Timepoint",  y="Diagnosis") +
  ggtitle("ADNI (Longitudinal)") +
  theme(
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(color = "black", size = 18, family = "Nimbus Sans Narrow"),
    axis.ticks = element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_blank(),
    axis.title.y = element_text(color = "black", size = 22, vjust =-1, margin = margin(0,20,0,0)),
    axis.title.x = element_text(color = "black", size = 22, vjust = 1, margin = margin(0,20,0,0)),
    axis.text = element_text(color = "black", size = 18))
)
(traj = traj + scale_color_manual(values=(c("#194e6a","#d8b707"))))


plotdir=here("plots")
if (! dir.exists(plotdir)) {
  dir.create(plotdir)
}
if (savefigs == 1) {
  ggsave(plot=traj,filename=file.path(plotdir,"trajectoriesADNI_ADC.png"),width=20,height=14,dpi=600,units="cm")
  ggsave(plot=traj,filename=file.path(plotdir,"trajectoriesADNI_ADCbigger.png"),width=15,height=10,dpi=600,units="cm")
}



DF$ICV = (DF$vol.EstimatedTotalIntraCranialVol-mean(DF$vol.EstimatedTotalIntraCranialVol)) / sd(DF$vol.EstimatedTotalIntraCranialVol)
DF$magstren = ifelse(DF$MagStrength == 3.0,1,0)
DF$csex = ifelse(DF$PTGENDER == "Female",0,1)
DF$csex = DF$csex - mean(DF$csex)
DF$cscan=DF$magstren-mean(DF$magstren)
DF$PTGENDER = as.factor(DF$PTGENDER)


BRAIN = readLines(file.path(a, "allFeatures364.txt"))
BRAIN[grepl("w-g", BRAIN)]
BRAIN = gsub("w-g", "w.g", BRAIN)
nrois=length(BRAIN)


ROIs = BRAIN
start = 1
N = 364
end = nrois
loopend = length(ROIs[start:end])
print(ROIs)


pb = txtProgressBar(min=2, max=end, style=3)
Usubs = nrow(U) #N unique subs


for (i in 1:loopend) {
  print(paste(i,"/",nrois))
  if (i == 1) {
    derivMat = absChangeMat = outlierMat1 = outlierMat2 = rIntonlyMat = rIntMat = rSlopeMat = matrix(NA, nrow = Usubs, ncol = loopend)
    larger=0
    predMat = seMat = matrix(NA, nrow = nrow(DF), ncol = loopend)
    RR = list()
  }
  setTxtProgressBar(pb,i)
  
  
  ROI = ROIs[i]
  print(ROI)
  # ROI = "vol.Left.Hippocampus"
  # ROI = "vol.Left.Inf.Lat.Vent"
  DF$brainvar = DF[[ROI]]
  
  
  # --- GAMM model ----
  g = gamm4(brainvar ~ s(Age, k = 5) + PTGENDER + magstren + ICV, data = DF, random = ~ (1+Age | RID))
  g.sum = summary(g$gam)
  dug = itsadug::get_predictions(g$gam,
                                 cond = list(Age = seq(
                                   min(DF$Age, na.rm = T),
                                   max(DF$Age, na.rm = T),
                                   length.out = nrow(DF)),
                                             PTGENDER="Female", magstren=1,se = T))

  
  
  #name and merge random effects
  rr = ranef(g$mer)$RID
  names(rr) = c("rInt", "rSlope")
  rr$RID = row.names(rr)
  U = DF %>% filter(timeBsl == 0)
  U = merge(U, rr, by = "RID")
  
  
  predictions <- DF %>% 
    mutate(PTGENDER = "Female", magstren = 1, ICV=0) %>% 
    select(Age,PTGENDER,magstren,ICV) %>% 
    predict(g$gam, newdata = ., se.fit = T)
  
  
  residualsg <- residuals(g$mer)
  DF$partial_residuals = predictions$fit + residualsg
  DF$fit = predictions$fit
  DF$sefit = predictions$se.fit
  DF$cifit = predictions$se.fit*1.96
  # plot(DF$Age,DF$brainvar)
  # points(DF$Age,DF$fit,col="blue")
   
   
  #correct for sex and field strength
  DF$brainvarcor = DF$brainvar
  DF$brainvarcor = ifelse(DF$PTGENDER == "Male", DF$brainvarcor - g.sum$p.coeff[2], DF$brainvarcor)
  DF$brainvarcor = ifelse(DF$magstren == 1, DF$brainvarcor - g.sum$p.coeff[3], DF$brainvarcor)
  DF$brainvarcor = ( DF$brainvarcor-(DF$ICV * g.sum$p.coeff[4]) )
  
  
  
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
  
  DF$group=ifelse(DF$includeLong=="D", "AD_long", "NC_long")
  
  (pf1=ggplot() +
      geom_line(data=DF,aes(x=Age,brainvarcor,group=RID,col=as.factor(group)),alpha=0.6, linewidth=0.5) +
      geom_point(data=DF,aes(x=Age,brainvarcor,group=RID,col=as.factor(group)),stat="identity",alpha=1, size=0.5) +
      geom_ribbon(data=dug,aes(x=Age,ymin=fit-CI,ymax=fit+CI),alpha=.7,show.legend=F,fill="dark grey") +
      geom_line(data=DF,aes(x=Age,y=fit),col="black") +
      # geom_line(data=dug,aes(x=Age,y=fit),col="red",size=1) +
      ylab("Volume") +
      ggtitle(ROI) +
      theme_classic() + mytheme + scale_color_manual(values=(c("#194e6a","#d8b707"))) +
      theme(
        axis.title.y = element_text(color = "black", size = 19, vjust =-1, margin = margin(0,10,0,0)),
        axis.title.x = element_text(color = "black", size = 19, vjust = -2, margin = margin(0,20,20,0)),
        legend.position = "top"
      ))
  
  
  #detect outliers
  outlierThresh = function(U, outlierthresh, ranef) {
    
    g1 = gam(get(ranef) ~ meanAge + csex + magstren, data = U, na.action = na.exclude)
    Opt.age = U$meanAge
    Q = scale(g1$residuals)
    db.out = data.frame(U$RID,Q) #residuals scaled Z
    db.out$QQ = NA
    db.out$QQ [ (abs(Q)>outlierthresh)] = Q[(abs(Q)>outlierthresh)]
    outliers = db.out$U.RID[which(!is.na(db.out$QQ))]
    cat("\noutlier thresh =", outlierthresh, "SD from lm: highlighting", length(outliers), "cases for potential removal\n\n")
    N.outliers = length(outliers)
    
    if (ranef == "rSlope") {
      cat("ranef is",ranef, "appending Qcol1")
      U$Qcol1="cyan2"
      U$Qcol1[(U$RID %in% outliers)]="black" 
    } else if (ranef == "rIntonly") {
      cat("ranef is",ranef, "appending Qcol2")
      U$Qcol2="cyan2"
      U$Qcol2[(U$RID %in% outliers)]="black" 
    }
    return(list(outliers = as.numeric(as.character(outliers)),
                U,
                Q))
  }
  thr = 5
  # U = outlierThresh(U, thr, "rSlope")[[2]] #returns colored outliers
  outlierMat1[,i] = outlierThresh(U, thr, "rSlope")[[3]] #rSlope
  
  # ggplot(data=U,aes(x=meanAge,rSlope)) +
  #   geom_point(color="blue",stat="identity",alpha=0.7) +
  #   geom_smooth(method="gam",col="black") +
  #   theme_bw()
  
  # ggplot(data=U,aes(x=meanAge,rSlope,color=Status)) +
  #   geom_point(stat="identity",alpha=0.7) +
  #   geom_smooth(method="lm",col="black") +
  #   theme_bw()
  
  
  U$Status = factor(U$Status,levels = c("CN", "AD") )
  t2=t.test(U$rSlope[U$Status=="CN"],U$rSlope[U$Status=="AD"])
  print(t2)
  RR$lm2[[i]] = summary(lm(rSlope ~ Status, U))
  RR$tt2[[i]]=t2
  RR$tt2stat[i]=t2$statistic
  RR$tt2p[i]=t2$p.value  
  RR$lm2stat[i]=coef(RR$lm2[[i]])[2,1]
  RR$lm2p[i]=coef(RR$lm2[[i]])[2,4]
  
  
  print(paste(ROI, "Slope effect = ", t2$statistic))
  
  
  U$group=ifelse(U$includeLong=="D", "AD_long", "NC_long")
  
  
  (pf2 = ggplot(U, aes(x=Status,y=rSlope,col=group)) +
    geom_jitter(width=0.1) +
    geom_boxplot(alpha=0) +
    labs(x = "group", title = ROI) + #, subtitle = "SLOPES") +
    scale_color_manual(values=(c("#194e6a","#d8b707"))) +
    theme_classic() + mytheme + ggtitle("Age-relative change") + labs(y="Additional Change (mm/year)") + theme(
    axis.title.y = element_text(color = "black", size = 19, vjust =-1, margin = margin(0,10,0,0)),
    axis.title.x = element_text(color = "black", size = 19, vjust = -2, margin = margin(0,20,20,0)),
    plot.margin = unit(c(1,2,1,1), "lines"),
    legend.position = "top"
  ))
  
  traj = traj + theme(legend.position = "top",
               axis.title.x = element_text(color = "black",size = 19,vjust = 6),
               axis.text.x = element_text(color = "black",size = 19,vjust = 8))
                                           
  # cp1=cowplot::plot_grid(traj,pf1,pf2,ncol=3)
  
  # ggsave(
  #   filename = file.path(plotdir, paste0(
  #     "pADNIMethod.",ROI,".png"
  #   )),
  # 
  #   plot = cp1,
  #   width = 30,
  #   height = 13,
  #   dpi = 600,
  #   units = "cm"
  # )

  
  #original model derivatives
  Udat = U %>% select(RID, Age, PTGENDER, magstren,ICV)
  dd=derivatives(object=g$gam,
                 newdata=Udat,
                 order = 1L,
                 n = nrow(U),
                 interval = c("confidence")
  )
  
  
  # plot(dd$data,dd$derivative,col="blue")
  #add random effect to curve derivative to get absolute change vals of slope
  U$absChange=U$rSlope+dd$derivative
  plot(U$meanAge, U$absChange)
  
  ggplot(data=U,aes(x=meanAge,absChange,color=Status)) +
    geom_point(stat="identity",alpha=0.7) +
    # geom_smooth(method="lm",col="black") +
    theme_bw()
  
  
  #save outputs
  derivMat[,i] = dd$derivative
  absChangeMat[,i] = U$absChange
  rSlopeMat[,i] = U$rSlope
  predMat[,i] = DF$fit
  seMat[,i] = DF$sefit
  RR$sTab[[i]] = g.sum$s.table
  RR$pTab[[i]] = g.sum$p.table
  
  
  # ggplot() +
  #   geom_point(data=U,aes(x=meanAge,absChange,col=as.factor(Status)),stat="identity",alpha=0.4) +
  #   # geom_smooth(data=U,aes(x=meanAge,absChange,group=as.factor(Status),col=as.factor(Status))) +
  #   geom_smooth(data=U,aes(x=meanAge,absChange,group=as.factor(Status),col=as.factor(Status)),method="lm",alpha=0.2) +
  #   ylab("brainvar change (mm per year)") +
  #   theme_classic()
  
  
}


#save
resdir=file.path(a,"ADNI/results/prepSlopes")
if (! dir.exists(resdir)) {
  system(paste("mkdir -p", resdir))
}


save(
  "ROIs",
  'U',
  "derivMat",
  "absChangeMat",
  "rSlopeMat",
  "predMat",
  "seMat",
  "RR",
  file = file.path(resdir,paste0("results.prepSlopesADNI.Rda")))

quit()
