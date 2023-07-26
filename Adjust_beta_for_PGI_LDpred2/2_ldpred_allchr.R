#!/usr/bin/env Rscript

## This script runs LDpred2-auto on each chromosome in parallel. 
  # required inputs (see arguments):
      # - matched summary stat file (from 1_match_sumstat.R)
      # - LD estimates from reference panel
  # output: LD-adjusted beta estimates

library(bigsnpr)
library(data.table)
library(future.apply)

# Collect arguments
   # 1: path to mapping file for the reference panel
   # 2: path to LD estimates of the reference panel
   # 3: output file prefix (i.e., prefix for the matched summary stat file)
arguments <- commandArgs(trailingOnly = T)
map_path <- as.character(arguments[1])
LD <- as.character(arguments[2]) 
out <- as.character(arguments[3])

############################################################
# Functions

# Function: prepare the inputs for LDpred2 for each chromosome.
    # input: chromosome number (CHR)
    # output: a list containing:
        # - path to temporary on-disk LD matrix filtered for the target SNPs
        # - LDSC estimates  
prep_chr <-  function(CHR){
    ind <- which(map[chr==CHR, marker.ID] %in% match[chr==CHR, rsid])  # get index of matched snps in genetic file

    LD <- readRDS(paste0(LD, "_chr", CHR, ".Rds")) # read LD matrix 
    LD <- LD[ind, ind]  # leave only the matched snps.

    ldsc <- snp_ldsc2(LD, match[chr==CHR])
    cat("CHR ", CHR,  ", LDSC h2: ", ldsc[2], "\n")

    tmp <- tempfile(tmpdir = "/scratch-local/hkweon/")
    LD_sfbm <- as_SFBM(LD, tmp, compact=TRUE)

    rm(ind, LD)
    gc(verbose=FALSE)

    return(list(SFBM=LD_sfbm, LDSC=ldsc))
}


# Function: run LDpred2-auto function for each chromosome.
    # input: chromosome number (CHR)
    # output: NULL
        # the results are saved in RDS file
run_chr <- function(CHR){
    s <- Sys.time()
    h2_i = ifelse(CHR_LD[[CHR]]$LDSC[2]<0 | is.na(CHR_LD[[CHR]]$LDSC[2]), 1e-100, CHR_LD[[CHR]]$LDSC[2])
    new_beta <- snp_ldpred2_auto(CHR_LD[[CHR]]$SFBM, match[chr==CHR], h2_init = h2_i,  vec_p_init = 0.2,   shrink_corr = 1, allow_jump_sign = TRUE)

    file <- paste0("../TEMP/", out, "_LDpred_chr", CHR, ".Rds")
    saveRDS(new_beta, file)

    e <- Sys.time()
    cat("\n Results for chr ", CHR, " saved to ", file, ", elapsed time: ", difftime(e,s,units="mins"), "mins  \n")

    cat("CHR ", CHR, " estimated p: ", new_beta[[1]]$p_est, "\n")
    cat("CHR ", CHR, " estimated h2: ", new_beta[[1]]$h2_est, "\n")
    gc(verbose=FALSE)

    return(NULL)
}
############################################################333



### MAIN RUN ###

s <- Sys.time()

cat("Loading data \n")
map <- readRDS(map_path)
match <- readRDS(paste0("../TEMP/", out, "_matched.Rds"))   # matched sumstats


cat("Reading LD matrix and computing h2 estimate for each chromosome \n")

options(future.globals.maxSize = 5 * 1024^3)
plan(multicore, workers = 14)
CHR_LD <- future_lapply(1:22, prep_chr)
gc(verbose=FALSE)


cat("\n\n Computing LDpred-adjusted Beta for each chromosome \n")

options(future.globals.maxSize= 3.2*1024^3)
plan(multicore, workers = 18)
RUN <- future_lapply(1:22, run_chr)

e <- Sys.time()

cat("Completed,  elapsed time: ", difftime(e,s,units="mins"), "mins \n")


# Delete the on-disk temp files
for (CHR in 1:22){
system(paste0("rm ", CHR_LD[[CHR]]$SFBM$backingfile))
}