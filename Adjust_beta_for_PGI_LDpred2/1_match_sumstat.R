#!/usr/bin/env Rscript

## This script prepares the summary statistics input file for LDpred2. 
  # required inputs (see arguments):
      # - summary stat file
      # - reference panel
      # - prediction sample plink file 
  # output: the summary stat file matched with the reference panel     


library(bigsnpr)
library(data.table)

# Collect arguments
   # 1: path to summary statistics 
   # 2: plink file prefix of prediction sample
   # 3: path to mapping file for the reference panel
   # 4: input file format   
   # 5: output file prefix
arguments <- commandArgs(trailingOnly = T)
sumstat_path <- as.character(arguments[1])
g_path <- as.character(arguments[2])
map_path <- as.character(arguments[3])
format <- as.character(arguments[4])
out <- as.character(arguments[5])

# Read inputs
map <- readRDS(map_path)

cat("Reading the list of SNPs in the prediction sample \n")
snp_pred <- fread(paste0(g_path, ".bim"))[[2]]
map <- map[marker.ID %in% snp_pred]  # leave only SNPs available in prediction sample 

cat("Reading summary statistics file \n")
sumstat <- fread(sumstat_path)

# Re-structure the summary stat given the input format
cat("Input format = ", format, "\n")
switch(format,
         "easyqc" = {
          sumstat <- sumstat[, .(SNP, Chr, position, EFFECT_ALLELE, OTHER_ALLELE, BETA, SE, PVAL, N)]
          names(sumstat) <-  c("rsid", "chr", "pos", "a0", "a1", "beta", "beta_se", "p", "n_eff")  
          },
         "mtag_meta" = {
          sumstat[, n_eff := 1/(mtag_se^2 * 2* meta_freq * (1-meta_freq))]   
          sumstat <- sumstat[, c(1,2,3,4,5, 7,8,10, 11)] 
          names(sumstat) <-  c("rsid", "chr", "pos", "a0", "a1", "beta", "beta_se", "p", "n_eff") 
          },
         "mtag" = {
          sumstat[, n_eff := 1/(mtag_se^2 * 2* meta_freq * (1-meta_freq))]   
          sumstat <- sumstat[, c(1,2,3,4,5, 9,10,12)] },
         "ldpred" = {
          sumstat <- sumstat
          names(sumstat) <-  c("rsid", "chr", "pos", "a0", "a1", "beta", "beta_se", "p", "n_eff") 
          }    )

# Filter out SNPs not in the reference panel
sumstat <- sumstat[rsid %in% map$marker.ID]
sumstat <- sumstat[n_eff > max(n_eff)*0.5] # drop SNPs with noo small N

cat("Summary statistics look: \n")
print(head(sumstat))

# Match SNPs with prediction sample 
cat("Matching SNPs with prediction sample \n")
match <- data.table(snp_match(sumstat, map[, -2]))

cat("Average N = ", mean(sumstat$n_eff), "\n")

# Save
file <- paste0("../TEMP/", out, "_matched.Rds")

saveRDS(match, file)
cat("Saved to ", file, "\n")