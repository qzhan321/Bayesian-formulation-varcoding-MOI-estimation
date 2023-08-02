rm(list=ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
})

option_list = list(
  make_option(c("--input"), type="character", default=NULL,
              help="path and name of the input file", metavar="character"),
  make_option(c("--prior"), type="character", default="uniform",
              help="user-specified prior for MOI distribution [default= %default]; two options, either a uniform or a negative binomial (negBinom)", metavar="character"),
  make_option(c("--params"), type="character", default="medium",
              help="when specifying a negative binomial prior, providing a parameter range for it [default= %default]; three options, low or medium or high, corresponding to a negative binomial with a low (~1), medium (~4), and high (~7) mean value", metavar="character"),
  make_option(c("--verbose"), type="logical", default=TRUE,
              help="whether output the prior distribution for MOI estimation or not [default= %default]", metavar="character"),
  make_option(c("--aggregate"), type="character", default="pool",
              help="how to obtain the MOI distribution at the population level from individual MOI estimates; two options, pool or mixtureDist, either pooling the maximum a posteriori MOI estimate for each individual host together, or using mixture distribution, i.e., a weighted summation of all individual MOI distributions", metavar="character"),
  make_option(c("--util"), type="character", default="./s_givenMOI_list",
              help="the list which stores the probability distribution of the number of non-upsA DBLa types given any MOI", metavar="character"),
  make_option(c("--output"), type="character", default="./out.txt",
              help="path and name of the output file [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

a <- 10 # lowest number of non-upsA DBLa types sequenced in a single infection, from the empirical repertoire size distribution.
b <- 45 # highest number of non-upsA DBLa types sequenced in a single infection, from the empirical repertoire size distribution
MOI_max <- 20

inputFile <- read.csv(opt$input, header = T, row.names = NULL)
utilFile <- opt$util
load(utilFile)
p_s_givenMOI <- s_givenMOI_list


simNegBinomPrior <- function(n, prob, MOI_min = 1, MOI_max = 20) {
  baseline <- 0.00001
  MOI_priors_test <- rnbinom(20000, size = n, prob = prob)
  MOI_priors_test <- MOI_priors_test[MOI_priors_test>=MOI_min & MOI_priors_test <= MOI_max]
  MOI_priors_test <- prop.table(table(MOI_priors_test))
  labels <- rownames(MOI_priors_test)
  MOI_priors <- rep(baseline, MOI_max)
  names(MOI_priors) <- as.character(1:MOI_max)
  MOI_priors[labels] <- MOI_priors_test + baseline
  MOI_priors <- MOI_priors/sum(MOI_priors)
  return(MOI_priors)
} 


if (opt$prior == "uniform") {
  MOI_priors <- rep(1/MOI_max, MOI_max)
  names(MOI_priors) <- as.character(1:MOI_max)
} else if (opt$prior == "negBinom") {
  if (opt$params == "low") {
    set.seed(0)
    MOI_priors <- simNegBinomPrior(n = 1, prob = 0.75, MOI_max = MOI_max)
  } else if (opt$params == "medium") {
    set.seed(1)
    MOI_priors <- simNegBinomPrior(n = 5, prob = 0.55, MOI_max = MOI_max)
  } else if (opt$params == "high") {
    set.seed(2)
    MOI_priors <- simNegBinomPrior(n = 4, prob = 0.35, MOI_max = MOI_max)
  }
}

p <- p_s_givenMOI[[1]]
s_single <- max(45, b)

s_temp <- inputFile$NumDBLaTypes
names <- inputFile$HostID

dfAll <- NULL
for (j in 1:length(s_temp)) {
  s <- s_temp[j]
  if (s < 10 | s > 900) {
    warning("The number of non-upsA DBLa type is either fewer than 10 or greater than 900. Preprocess the input file to filter those isolates out.")
  }
  
  numerator <- 0
  for (b in 1:MOI_max) {
    if (!is.na(p_s_givenMOI[[as.character(b)]][as.character(s)])) {
      numerator <- numerator + p_s_givenMOI[[as.character(b)]][as.character(s)]*MOI_priors[as.character(b)]
    }
  }
  
  stopifnot(numerator != 0)
  
  p_c_givens_all <- rep(NA, length(1:MOI_max))
  for (c in 1:MOI_max) {
    if (!is.na(p_s_givenMOI[[as.character(c)]][as.character(s)])) {
      temp <- p_s_givenMOI[[as.character(c)]][as.character(s)]*MOI_priors[as.character(c)]/numerator
      names(temp) <- NULL
    } else {
      temp <- 0
    }
    p_c_givens_all[c] <- temp
  }
  
  if (opt$aggregate == "pool") {
    df <- data.frame("p" = max(p_c_givens_all), "MOI" = c(1:MOI_max)[which.max(p_c_givens_all)], "NumDBLaTypes" = s, "HostID" = names[j])
  } else if (opt$aggregate == "mixtureDist") {
    df <- data.frame("p" = p_c_givens_all, "MOI" = 1:MOI_max, "numDBLaTypes" = s, "HostID" = names[j])
  }
  dfAll <- rbind(dfAll, df)
}

if (opt$aggregate == "pool") {
  dfAllPop <- dfAll %>% group_by(MOI) %>% summarise("n" = n()) %>% mutate("n_total" = sum(n), "prob" = n/n_total) %>% select(MOI, prob)
} else if (opt$aggregate == "mixtureDist") {
  dfAllPop <- dfAll %>% group_by(MOI) %>% summarise("n" = sum(p)) %>% mutate("n_total" = sum(n), "prob" = n/n_total) %>% select(MOI, prob)
} else if (opt$aggregate == "mixtureDist") {
}
if (opt$verbose == TRUE) {
  outputList <- list("prior" = MOI_priors,
                     "indLevel" = dfAll,
                     "popLevel" = dfAllPop)
} else {
  outputList <- list("indLevel" = dfAll,
                     "popLevel" = dfAllPop)
}

save(outputList, file = opt$output)
