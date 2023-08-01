rm(list = ls())
readDir <- "" # the directory where the repertoire size distribution is stored.
f <- read.table(paste0(readDir, "NbVarGenes_MOI1_upsBC_AllSurveys_Weight.txt"), header = T) 
hist(rep(f$DBLa_upsBC_rep_size, f$n)) 

a <- min(f$DBLa_upsBC_rep_size)
b <- max(f$DBLa_upsBC_rep_size)
s_single <- max(45, b)

s_givenMOI_list <- list() # the list which records the probability distribution of s, i.e., the number of DBLa types sequenced and typed, given a certain MOI. The first element of the list corresponds to MOI = 1, and so on forth. 

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
