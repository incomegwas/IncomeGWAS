library(GenomicSEM)

# Prerequisite: the input files are formatted in a compatible way for genomic SEM.
    # see: https://github.com/GenomicSEM/GenomicSEM/wiki

##################################################
# prepare inputs

##### LDSC #####
list <- c("EA.sumstats.gz", "INC.sumstats.gz")
list_name <- c("EA", "INC")

LD="../eur_w_ld_chr"
LDSC <- ldsc(traits=list, sample.prev=rep(NA, length(list)), population.prev=rep(NA, length(list)), ld=LD, wld=LD, trait.names=list_name)

##### matched sumstats #####
list_name <- c("EA", "INC")
list <- paste0("../INPUT/", list_name, "_input_GSEM.txt")

p_sumstats <-sumstats(files=list, ref="HRC_ref_GSEM.gz", trait.names=list_name,  se.logit=NULL, OLS=rep(TRUE, length(list)), linprob=NULL,prop=NULL, keep.indel=FALSE, parallel=FALSE, cores=length(list))

###############################################

# Run genomic SEM

model<-'
INC~d*SNP
EA~a*SNP
INC~b*EA
i:=a*b
'
#i = indirect effect via EA
#d = direct effect (NonEA-INC)
output <- userGWAS(covstruc=LDSC, SNPs=p_sumstats, estimation="DWLS", model=model, cores=110, parallel=TRUE)

saveRDS(output, "NonEA_INC_OUTPUT.RDS")
