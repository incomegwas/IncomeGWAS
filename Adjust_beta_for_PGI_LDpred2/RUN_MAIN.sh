######################
###LDpred2 pipeline###
######################


##### prepare inputs ######3

# reference panel #
ref_map={1}
LD={2}

# prediction sample plink file prefix #
g_path={3}

# GWAS summary statistics and its format
sumstat={4}
format={5}

# output prefix
output={6}


##################################
## 1. Create matched sumstat file
##################################
Rscript --verbose 1_match_sumstat.R ${sumstat} ${g_path} ${ref_map} ${output} ${format} 

##################################
## 2. Run LDpred main step 
##################################
Rscript --verbose ${ref_map} ${LD} ${output}