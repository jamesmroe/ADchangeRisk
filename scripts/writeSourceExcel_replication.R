rm(list=ls())

#make source data file

# options(bitmapType = "cairo")
# resdir="/cluster/p/p274/cluster/projects/p040-ad_change/ADchangeRisk/reproduce/results"
# setwd(resdir)
library("here")


# load packages ----
library("magrittr")
library("tidyverse")
library("sgof")
library("openxlsx")



# ROIs ----
roinames = c("vol.Left.Hippocampus", #L hippo (stage II)
             "vol.Right.Hippocampus", #R hippo (stage (II)
             "vol.Left.Amygdala", #L amygdala (stage II)
             "vol.Right.Amygdala") #R amygdala (stage III)


nTime=2

simulated=0
if (simulated == 1) {
  #change loop to be 1:1
  roinames="Left-Hippocampus.fs71"
  roiname=roinames
}


# load and merge results ----
# load and merge results from Rdata files ----
for (roiname in roinames) {
  if (roiname==roinames[1]) {
    SLOPESALL1=c()
    SLOPESALL2=c()
    SLOPESALL3=c()
    SLOPESALL4=c()
  }
  
  #main results
  print("loading replication results")
  filename = paste0(here(paste0("results/PRS-ADmodels_S1_",roiname,"_replication.publication.Rda")))
  print(paste0("loading ",filename))
  load(filename)
  
  
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


allPRSADout = rbind(SLOPESALL1 %>% mutate(change = "ageRelChange"),
                    SLOPESALL3 %>% mutate(change = "absChange"))
allPRSADnoAPOEout = rbind(SLOPESALL2 %>% mutate(change = "ageRelChange"),
                          SLOPESALL4 %>% mutate(change = "absChange"))


allPRSADnoAPOEout[allPRSADout$p.value > .05,c(2:6, 13:15)] = NA
allPRSADout %<>% mutate(sample ="REPLICATION")
allPRSADnoAPOEout %<>% mutate(sample ="REPLICATION")
allPRSADout %<>% mutate(sig_PRSAD = ifelse(p.value < .05, 1, 0)) %>% dplyr::select(model, sig_PRSAD, everything())
allPRSADnoAPOEout$sig_PRSAD = allPRSADout$sig_PRSAD
allPRSADnoAPOEout %<>% dplyr::select(model, sig_PRSAD, everything())


rename_cols = function(dat) {
  dat %<>% dplyr::select(-agecut) %>% rename(
    lower_agerange_genetic = agecut2,
    n_genetic = N,
    n_lifespan = Nlife
  )
  return(dat)
}



# write to excel source file ----
makeExcel = 1
if (makeExcel) {
  
  wb <- createWorkbook()
  header_style <- createStyle(textDecoration = "bold")
  
  # Add sheets
  sheets = c("Fig_5a PRS-AD rep",
             "Fig_5a noAPOE rep",
             
             "Fig_5b PRS-AD rep",
             "Fig_5b noAPOE rep",
             
             "Fig_5c-d PRS-AD rep"
  )
  for (sheetnum in 1:length(sheets)) {
    addWorksheet(wb, sheets[sheetnum])
  }
  
  df1 = allPRSADout
  df2 = allPRSADnoAPOEout
  df1$model = as.numeric(df1$model)
  df2$model = as.numeric(df2$model)
  df1$SS = NULL; df2$SS = NULL
  df1$filtBet = NULL; df2$filtBet = NULL; df2$`filtBet == filtBet` = NULL
  df1$cohort = NULL; df2$cohort = NULL
  
  header_style <- createStyle(textDecoration = "bold")
  
  writeData(
    wb,
    sheet = sheets[1],
    rename_cols(rbind(
      df1[df1$roi == unique(df1$roi)[1],],
      df1[df1$roi == unique(df1$roi)[2],]))
  )
  
  writeData(
    wb,
    sheet = sheets[2],
    rename_cols(rbind(
      df2[df2$roi == unique(df2$roi)[1],],
      df2[df2$roi == unique(df2$roi)[2],]))
  )
  
  writeData(
    wb,
    sheet = sheets[3],
    rename_cols(rbind(
      df1[df1$roi == unique(df1$roi)[3],],
      df1[df1$roi == unique(df1$roi)[4],]))
  )
  
  writeData(
    wb,
    sheet = sheets[4],
    rename_cols(rbind(
      df2[df2$roi == unique(df2$roi)[3],],
      df2[df2$roi == unique(df2$roi)[4],]))
  )
  
  df1 = rename_cols(df1)
  df2 = rename_cols(df2)
  new_names = names(df1)
  
  # add multivariate results from replication sample
  load(here("results/PRS-ADmodels_figure5_PCmultivariate_replication.Rda"))
  
  #harmonize colnames
  common_names = names(figure5C_PRSAD)[names(figure5C_PRSAD) %in% new_names]
  names(figure5C_PRSAD)[!names(figure5C_PRSAD) %in% new_names]
  new_names[!new_names %in% names(figure5C_PRSAD)]
  
  munge_cols = function(dat, roi) {
    dat$roi = roi
    dat$roiname = dat$roi
    dat %<>% rename(lower_agerange_genetic = agecut, apoe = ap)
    dat %<>% mutate(change = "ageRelChange",
                    model = row_number(),
                    HC = "noHCv_noAmyV"
    )
    dat %<>% select(-contains("conf"), -type)
  }
  
  figure5C_PRSAD = munge_cols(figure5C_PRSAD, "PC1relChange")
  figure5C_PRSAD %<>% select(-absestimateneg)
  new_names[!new_names %in% names(figure5C_PRSAD)]
  names(figure5C_PRSAD)[!names(figure5C_PRSAD) %in% new_names]
  
  
  df7 = figure5C_PRSAD
  
  
  common_names = new_names[new_names %in% names(df7)]
  df7 %<>% dplyr::select(all_of(common_names), everything()) %>% select(-SS, -dotalpha, -ww, -HC) %>% mutate(sample = "REPLICATION")
  
  
  writeData(wb, sheet = sheets[5], df7)
  
  df_list = list(df1, df2, df7)
  maxcols = max(sapply(df_list, ncol))
  unique(unlist(sapply(df_list, names)))
  
  
  for (num in 1:length(sheets)) {
    addStyle(wb, sheet = sheets[num], style = header_style, rows = 1, cols = 1:maxcols, gridExpand = TRUE)
    setColWidths(wb, sheet = sheets[num], cols = 1:maxcols, widths = "auto")
  }
  saveWorkbook(wb, file = "source_data_replication.xlsx", overwrite = TRUE)
  
    
}
