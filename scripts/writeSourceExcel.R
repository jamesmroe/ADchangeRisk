#rm(list=ls())

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
roinames=c("Left-Hippocampus.fs71", #L hippo (stage II)
           "Right-Hippocampus.fs71", #R hippo (stage (II)
           "braak-stage-1ent-.volume-lh-71", #L entorhinal (stage I)
           "braak-stage-1ent-.volume-rh-71", #R entorhinal (stage I)
           "braak-stage-3-amy-.volume-lh-71", #cortical stage III ROI excluding amygdala L
           "braak-stage-3-amy-.volume-rh-71", #cortical stage III ROI excluding amygdala R
           "braak-stage-3-amygdala-.volume-lh-71", #L amygdala (stage II)
           "braak-stage-3-amygdala-.volume-rh-71") #R amygdala (stage III)

nTime=2

simulated=0
if (simulated == 1) {
  #change loop to be 1:1
  roinames="Left-Hippocampus.fs71"
  roiname=roinames
}


# load and merge results ----
for (roiname in roinames) {
  if (roiname==roinames[1]) {
    SLOPESALL1=c()
    SLOPESALL2=c()
    SLOPESALL3=c()
    SLOPESALL4=c()
  }
  
  if (simulated == 1) {
    print("loading simulated results")
    load(here("simulate/PRSADmodels_simulated.Rda"))
    
  } else if (simulated == 0) {
    #main results
    filename=paste0(here(paste0("results/PRSADmodels_S1_",roiname,".publication.Rda")))
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



#DF of all 576 PRS-AD tests
allPRSAD = rbind(SLOPESALL1 %>% mutate(change = "ageRelChange"),
                 SLOPESALL3 %>% mutate(change = "absChange"))


#144 FDR-corrected PRS-AD tests
(ALLFDR=BH(allPRSAD$p.value, alpha=0.05))
FDRthresh = max(ALLFDR$data[ALLFDR$Adjusted.pvalues<.05])


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


#add better ROI names
SLOPESALL1$roiname = SLOPESALL2$roiname = SLOPESALL3$roiname = SLOPESALL4$roiname = 
  gsub("\\)", "", gsub("\\(", ".", gsub(" ", "-", SLOPESALL1$roi)))


allPRSADout = rbind(SLOPESALL1 %>% mutate(change = "ageRelChange"),
                    SLOPESALL3 %>% mutate(change = "absChange"))
allPRSADnoAPOEout = rbind(SLOPESALL2 %>% mutate(change = "ageRelChange"),
                          SLOPESALL4 %>% mutate(change = "absChange"))


allPRSADnoAPOEout[allPRSADnoAPOEout$FDRsig != 1,c(2:6, 13:15)] = NA
allPRSADout %<>% dplyr::select(model, FDRsig, everything()) %>% 
  rename(FDRsig_PRSAD = FDRsig)
allPRSADnoAPOEout %<>% dplyr::select(model, FDRsig, everything()) %>% 
  rename(FDRsig_PRSAD = FDRsig)


if (simulated) {
  allPRSADout %<>% mutate(SIMULATED ="SIMULATED_RESULTS")
  allPRSADnoAPOEout %<>% mutate(SIMULATED ="SIMULATED_RESULTS")
}


rename_cols = function(dat) {
  dat %<>% dplyr::select(-agecut) %>% rename(
    lower_agerange_genetic = agecut2,
    n_genetic = N,
    n_lifespan = Nlife
  )
}


# write to excel source file ----
makeExcel = 1
if (makeExcel) {
  
  wb <- createWorkbook()
  header_style <- createStyle(textDecoration = "bold")
  
  if (!simulated) {
    
    # Add sheets
    sheets = c("Fig_1e-f PRS-AD",
               "Fig_1e-f noAPOE",
               
               "Fig_2a PRS-AD",
               "Fig_2a noAPOE",
               
               "Fig_2b PRS-AD",
               "Fig_2b noAPOE",
               
               "Fig_2c PRS-AD",
               "Fig_2c noAPOE",
               
               "Fig_4b-c PRS-AD",
               "Fig_4b-c noAPOE",
               
               "Fig_4e PRS-AD",
               "Fig_4e noAPOE",
               
               "Fig_5c-d replication"
    )
    for (sheetnum in 1:length(sheets)) {
      addWorksheet(wb, sheets[sheetnum])
    }
    
    df1 = allPRSADout
    df2 = allPRSADnoAPOEout
    df1$model = as.numeric(df1$model)
    df2$model = as.numeric(df2$model)
    df1$SS = NULL; df2$SS = NULL
    
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
    
    writeData(
      wb,
      sheet = sheets[5],
      rename_cols(rbind(
        df1[df1$roi == unique(df1$roi)[7],],
        df1[df1$roi == unique(df1$roi)[8],]))
    )  
    
    writeData(
      wb,
      sheet = sheets[6],
      rename_cols(rbind(
        df2[df2$roi == unique(df2$roi)[7],],
        df2[df2$roi == unique(df2$roi)[8],]))
    )
    
    writeData(
      wb,
      sheet = sheets[7],
      rename_cols(rbind(
        df1[df1$roi == unique(df1$roi)[5],],
        df1[df1$roi == unique(df1$roi)[6],]))
    )
    
    writeData(
      wb,
      sheet = sheets[8],
      rename_cols(rbind(
        df2[df2$roi == unique(df2$roi)[5],],
        df2[df2$roi == unique(df2$roi)[6],]))
    )
    
    # for (sheetnum in 1:length(sheets)) {
    #   addStyle(wb, sheet = sheets[sheetnum], style = header_style, rows = 1, cols = 1:ncol(df1), gridExpand = TRUE)
    #   setColWidths(wb, sheet = sheets[sheetnum], cols = 1:ncol(df1), widths = "auto")
    # }
    # saveWorkbook(wb, file = "source_data.xlsx", overwrite = TRUE)
  
    
    # add multivariate results from main sample
    load(here("results/PRS-ADmodels_figure4_PCmultivariate.Rda"))
    df1 = rename_cols(df1)
    df2 = rename_cols(df2)
    new_names = names(df1)

    #harmonize colnames
    common_names = names(figure4B_PRSAD)[names(figure4B_PRSAD) %in% new_names]
    names(figure4B_PRSAD)[!names(figure4B_PRSAD) %in% new_names]
    new_names[!new_names %in% names(figure4B_PRSAD)]

    munge_cols = function(dat, roi) {
      dat$roi = roi
      dat$roiname = dat$roi
      dat %<>% rename(lower_agerange_genetic = agecut, apoe = ap, FDRsig_PRSAD = FDRsig)
      dat %<>% mutate(change = "ageRelChange",
                      model = row_number(),
                      hippocampus = "noHippocampus_noAmygdala",
                      FDRsig_PRSAD = as.numeric(as.character(FDRsig_PRSAD)))
      dat %<>% select(-contains("conf"), -type, HC)
    }

    figure4B_PRSAD = munge_cols(figure4B_PRSAD, "PC1relChange")
    figure4D_PRSAD = munge_cols(figure4D_PRSAD, "PC1relChange")
    figure4B_PRSADnoAPOE$FDRsig = figure4B_PRSAD$FDRsig_PRSAD
    figure4D_PRSADnoAPOE$FDRsig = figure4D_PRSAD$FDRsig_PRSAD
    figure4B_PRSADnoAPOE$dotalphapoe = as.numeric(as.character(figure4B_PRSADnoAPOE$dotalphapoe))
    figure4D_PRSADnoAPOE$dotalphapoe = as.numeric(as.character(figure4D_PRSADnoAPOE$dotalphapoe))
    figure4B_PRSADnoAPOE = munge_cols(figure4B_PRSADnoAPOE, "PC1relChangenoAPOE")
    figure4D_PRSADnoAPOE = munge_cols(figure4D_PRSADnoAPOE, "PC1relChangenoAPOE")
    

    new_names[!new_names %in% names(figure4B_PRSAD)]
    names(figure4B_PRSAD)[!names(figure4B_PRSAD) %in% new_names]

    new_names[!new_names %in% names(figure4B_PRSADnoAPOE)]
    names(figure4B_PRSADnoAPOE)[!names(figure4B_PRSADnoAPOE) %in% new_names]

    new_names[!new_names %in% names(figure4D_PRSAD)]
    names(figure4D_PRSAD)[!names(figure4D_PRSAD) %in% new_names]

    new_names[!new_names %in% names(figure4D_PRSADnoAPOE)]
    names(figure4D_PRSADnoAPOE)[!names(figure4D_PRSADnoAPOE) %in% new_names]

    link = df1 %>% select(lower_agerange_genetic, n_genetic, n_lifespan) %>% distinct()


    df3 = figure4B_PRSAD
    df4 = figure4B_PRSADnoAPOE
    df5 = figure4D_PRSAD
    df6 = figure4D_PRSADnoAPOE

    df3 = left_join(df3, link)
    df4 = left_join(df4, link)
    df5 = left_join(df5, link)
    df6 = left_join(df6, link)

    common_names = new_names[new_names %in% names(df3)]
    df3 %<>% dplyr::select(all_of(common_names), everything()) %>% select(-SS, -HC, -contains("dotalpha"))
    df4 %<>% dplyr::select(all_of(common_names), everything()) %>% select(-SS, -HC, -contains("dotalpha"))
    df5 %<>% dplyr::select(all_of(common_names), everything()) %>% select(-SS, -HC, -contains("dotalpha")) %>% rename(window = ww)
    df6 %<>% dplyr::select(all_of(common_names), everything()) %>% select(-SS, -HC, -contains("dotalpha")) %>% rename(window = ww)
    
    

    df4[df3$FDRsig_PRSAD != 1,c(2:6, 13:15)] = NA
    df6[df5$FDRsig_PRSAD != 1,c(2:6, 13:15)] = NA

    # newsheets = c("Fig_4b-c PRS-AD",
    #               "Fig_4b-c noAPOE",
    #               "Fig_4e PRS-AD",
    #               "Fig_4e noAPOE"
    # )

    # for (num in 1:length(newsheets)) {
    #   addWorksheet(wb, newsheets[num])
    # }

    writeData(wb, sheet = sheets[9], df3)
    writeData(wb, sheet = sheets[10], df4)
    writeData(wb, sheet = sheets[11], df5)
    writeData(wb, sheet = sheets[12], df6)

    # for (num in 1:length(newsheets)) {
    #   addStyle(wb, sheet = sheets[num], style = header_style, rows = 1, cols = 1:ncol(df6), gridExpand = TRUE)
    #   setColWidths(wb, sheet = sheets[num], cols = 1:ncol(df6), widths = "auto")
    # }
    # saveWorkbook(wb, file = "source_data.xlsx", overwrite = TRUE)
    
    
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
    df7 %<>% dplyr::select(all_of(common_names), everything()) %>% select(-SS, -dotalpha, -ww, -HC)

    # newsheets = c("Fig_5c-d replication"
    # )
    
    # for (num in 1:length(newsheets)) {
    #   addWorksheet(wb, newsheets[num])
    # }
    
    writeData(wb, sheet = sheets[13], df7)
    
    df_list = list(df1, df2, df3, df4, df5, df6, df7)
    maxcols = max(sapply(df_list, ncol))
    unique(unlist(sapply(df_list, names)))
    
    
    for (num in 1:length(sheets)) {
      addStyle(wb, sheet = sheets[num], style = header_style, rows = 1, cols = 1:maxcols, gridExpand = TRUE)
      setColWidths(wb, sheet = sheets[num], cols = 1:maxcols, widths = "auto")
    }
    saveWorkbook(wb, file = "source_data.xlsx", overwrite = TRUE)
    
    
  } else {
    
    # Add sheets
    sheets = c("SIMULATED PRS-AD",
               "SIMULATED noAPOE"
    )
    for (sheetnum in 1:length(sheets)) {
      addWorksheet(wb, sheets[sheetnum])
    }
    
    df1 = allPRSADout
    df2 = allPRSADnoAPOEout
    
    writeData(wb, sheet = sheets[1], rename_cols(df1))
    writeData(wb, sheet = sheets[2], rename_cols(df2))
    saveWorkbook(wb, file = "simulated.xlsx", overwrite = TRUE)
  }
}
