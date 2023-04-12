rm(list = ls())
suppressPackageStartupMessages({
  library(partitions)
  library(ggplot2)
})
readDir <- ""
f <- read.table(paste0(readDir, "NbVarGenes_MOI1_upsBC_AllSurveys_Weight.txt"), header = T)
hist(rep(f$DBLa_upsBC_rep_size, f$n))

# good sample
# f <- read.table(paste0(readDir, "NbVarGenes_3D7_upsBC_Weight.txt"), header = T)
# hist(rep(f$DBLa_upsBC_rep_size, f$n))

a <- min(f$DBLa_upsBC_rep_size)
b <- max(f$DBLa_upsBC_rep_size)
s_single <- max(45, b)

wd <- ""
s_givenMOI_list <- list()

p <- f$n/sum(f$n)
names(p) <- as.character(f$DBLa_upsBC_rep_size)

s_givenMOI_list[[1]] <- p
MOI_max <- 20

for (MOI in 2:MOI_max) {
  temp1 <- s_givenMOI_list[[1]]
  temp2 <- s_givenMOI_list[[MOI-1]]
  
  s_givenMOI_test <- rep(NA, s_single*MOI-a*MOI+1)
  names(s_givenMOI_test) <- as.character(seq(a*MOI, s_single*MOI))
  
  for (s in c(a*MOI):(s_single*MOI)) {
    tempAll <- 0
    for (bb in a:s_single) {
      temp <- temp1[as.character(bb)] * temp2[as.character(s-bb)]
      if(!is.na(temp)) {
        tempAll <- tempAll + temp
      }
    }
    s_givenMOI_test[[as.character(s)]] <- tempAll
  }
  s_givenMOI_list[[MOI]] <- s_givenMOI_test
}

saveDir <- ""
save(s_givenMOI_list, file = paste0(saveDir, "s_givenMOI_list"))






rm(list=ls())
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(cowplot)
  library(stringr)
})

readDirEpi <- ""
epi <- read.csv(paste0(readDirEpi, "Ghana_Survey_Merged_Epi_MOI_S1_S7_070721_UChicago_080822.csv"), header = T, row.names = 1)
epi$SeqID <- str_replace(epi$SeqID, "-", ".")

readDir <- ""
filesSub <- c("survey_1.csv", "survey_2.csv", "survey_3.csv", "survey_4.csv", "survey_5.csv", "survey_6.csv", "survey_7.csv")
# years <- c(2012, 2014, 2015, 2017)

load("s_givenMOI_list")
p_s_givenMOI <- s_givenMOI_list
MOI_max <- 20
names(p_s_givenMOI) <- as.character(1:MOI_max)

MOI_priors <- rep(1/MOI_max, MOI_max)
names(MOI_priors) <- as.character(1:MOI_max)

f <- read.table(paste0("NbVarGenes_MOI1_upsBC_AllSurveys_Weight.txt"), header = T)

a <- min(f$DBLa_upsBC_rep_size)
b <- max(f$DBLa_upsBC_rep_size)

p <- f$n/sum(f$n)
p <- f$weigth_n
names(p) <- as.character(f$DBLa_upsBC_rep_size)

MOI_max <- 20
s_single <- max(45, b)

for (i in 1:length(filesSub)) {
  file <- filesSub[i]
  # year <- years[i]
  
  file2 <- strsplit(file, split = "\\.")
  survey <- file2[[1]][1]
  
  print(survey)
  file_read <- read.csv(paste0(readDir, file), header = T, row.names = 1)
  
  s_temp <- rowSums(file_read == 1)
  # MOI <- file_read$MOI
  names <- names(s_temp)
  names(s_temp) <- NULL
  # names(MOI) <- NULL
  dfAll <- NULL
  for (j in 1:length(s_temp)) {
    s <- s_temp[j]
    MOIs_s <- rep(NA, MOI_max)
    numerator <- 0
    for (b in 1:MOI_max) {
      if (!is.na(p_s_givenMOI[[as.character(b)]][as.character(s)])) {
        numerator <- numerator + p_s_givenMOI[[as.character(b)]][as.character(s)]*MOI_priors[as.character(b)]
      }
    }
    
    if (numerator == 0) {
      p_c_givens_all <- rep(0, length(1:MOI_max))
    } else {
      p_c_givens_all <- rep(NA, length(1:MOI_max))
      # names(p_c_givens_all) <- as.character(1:MOI_max)
      for (c in 1:MOI_max) {
        if (!is.na(p_s_givenMOI[[as.character(c)]][as.character(s)])) {
          temp <- p_s_givenMOI[[as.character(c)]][as.character(s)]*MOI_priors[as.character(c)]/numerator
          names(temp) <- NULL
        } else {
          temp <- 0
        }
        p_c_givens_all[c] <- temp
      }
    }
    
    age <- epi %>% filter(SeqID == names[j])
    age <- age$AgeGroups2
    # plot(x=as.numeric(names(p_c_givens_all)), y=p_c_givens_all)
    # df <- data.frame("p" = max(p_c_givens_all), "MOI" = c(1:MOI_max)[which.max(p_c_givens_all)], "numDBLaTypes" = s, "host_id" = names[j],
    #                  "survey" = survey, "age" = age)
    df <- data.frame("p" = p_c_givens_all, "MOI" = 1:MOI_max, "numDBLaTypes" = s, "host_id" = names[j],
                     "survey" = survey, "age" = age)
    dfAll <- rbind(dfAll, df)
  }
  save(dfAll, file = paste0("surveysMOI/", survey))
}

