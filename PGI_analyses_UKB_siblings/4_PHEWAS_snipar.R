library(lme4)
library(lmerTest)

# Load data
UKB <- readRDS("UKBsib.Rds")

phecode_list=data.table(PheWAS::pheinfo)$phecode
phecode_list=paste0("phecode_", phecode_list)

IND_IN = UKB[male==1, lapply(.SD, mean, na.rm=TRUE)>0.01, .SDcols=phecode_list] & UKB[male==0, lapply(.SD, mean, na.rm=TRUE)>0.01, .SDcols=phecode_list]
phecode_list = phecode_list[IND_IN]

LABEL = data.table(PheWAS::pheinfo)$description[IND_IN]
GROUP = data.table(PheWAS::pheinfo)$group[IND_IN]

# Merge PGI data
PGI = fread("/INC_SNIPAR.pgs.txt", 
                    col.names=c("FID", "IID", "INC_PGI", "INC_PGI_M", "INC_PGI_F"))

i="INC"
sd = sd(PGI[, get(paste0(i, "_PGI"))])   
COLS = paste0(i, c("_PGI", "_PGI_M", "_PGI_P"))
PGI[, (COLS) := lapply(.SD, function(x) x / sd), .SDcols=COLS ]     

UKB[, IID := paste0(n_eid, "_", n_eid)]
UKB = merge(UKB, PGI, by="IID")


###################
pc <- paste0("pc", 1:20)
pc <- paste0(pc, collapse="+")
FM=paste0("male*(I((yob-1900)) + I((yob-1900)^2) ) + geno+", pc)

reg_func <- function(i){
LM1 = coeftable(fixest::feols(as.formula(paste0(phecode_list[i], "~ INC_PGI +", FM)), data = UKB, se = "cluster", cluster = "familyID"))[2, ]
LM2 <- coeftable(fixest::feols(as.formula(paste0(phecode_list[i], "~ INC_PGI + INC_PGI_P + INC_PGI_M + ", FM)), data = UKB, se = "cluster", cluster = "familyID"))[2, ]
DF = data.table(rbind(LM1, LM2))

names(DF) <- c("EST", "SE", "T", "P")
DF$M = c("POP", "D")

DF$PHECODE = phecode_list[i]
DF$LABEL = LABEL[i]
DF$GROUP = GROUP[i]
DF$N_case = nrow(UKB[get(phecode_list[i])==1])
DF$N_con = nrow(UKB[get(phecode_list[i])==0])
DF$N = DF$N_case + DF$N_con
DF$PREV = DF$N_case / DF$N
return(DF)
}

output = rbindlist(pbmcapply::pbmclapply(1:length(phecode_list), reg_func, mc.cores=32, ignore.interactive=TRUE))
fwrite(output, "PHEWAS_INC_SNIPAR.csv", sep="\t")



