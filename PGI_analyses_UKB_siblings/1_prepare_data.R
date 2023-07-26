UKB <- readRDS("") # UKB data with occupational and household income

#Identify siblings
related <- fread("ukb_kinship.dat")
related_sibs <- related[ (Kinship<=0.25+0.064&Kinship>=0.25-0.064) & IBS0 >= 0.002 ] # drop those with kinship coef too different from 0.25 and IBS0 < 0.002 (parent-child)

sibIDs <- unique( c(related_sibs$ID1, related_sibs$ID2) )
length(sibIDs) #41046
sibIDs <- sibIDs[sibIDs %in% UKB$n_eid] 
length(sibIDs) #14669 

library(pbmcapply)

familyID <- pbmclapply(sibIDs, function(x){
   ID2 <- related_sibs[which(related_sibs$ID1==x), ID2]
   ID1 <- related_sibs[which(related_sibs$ID2==x), ID1]
   res <- sort( unique( c(x, ID1, ID2) ))
   return(res[1])
})

sibs <- data.table(sibIDs, unlist(familyID))
setnames(sibs, "V2", "familyID")
 
UKB_sibs <- merge(UKB, sibs, by.x="n_eid", by.y="sibIDs")
UKB_sibs[, n_sibs:=.N, by=familyID]
UKB <- UKB_sibs[n_sibs>=2]

######################3
# EA

UKB_edu <- fread("UKB_edu.tab")
UKB <- merge(UKB, UKB_edu, by.x="n_eid", by.y="f.eid", all.x=T)

UKB[, edu:=NULL]
edu_fun <- function(x){
  x <- x[!is.na(x)&x!=-3]
  x <- ifelse(x==-7, 7, x)
  paste(sort(unique(x)), collapse="")
}

library(parallel)
cl <- makeCluster(23, type="FORK")
UKB[, edu_entered := parApply(cl, .SD, 1, edu_fun), .SDcols=grep("f.6138.", names(UKB), value=T)] 
stopCluster(cl)
UKB[, edu_entered := ifelse(edu_entered=="", NA, edu_entered)]


edu_func <- function(x){
edu <- x[!is.na(x)&x!=-3]

if (all(is.na(edu))){
    return(rep(NA, 7))
} else {
res <- sapply(c(1:6, -7), function(y) ifelse(any(edu==y),1,0) )
return(res)
}
}
UKB[, paste0("edu_dm", 1:7) := data.table(t(apply(.SD, 1, edu_func))) , .SDcols=grep("6138", names(UKB), value=T) ]

get_EA  <- function(x){
eduage <- x[8]
edu <- x[1:7]
EA <- c(20, 13, 10, 10, eduage-5, 15, 7)
EA <- edu*EA
max(EA, na.rm=T)
}
UKB[!is.na(edu_entered), edu := apply(.SD, 1, get_EA), .SDcols=c(paste0("edu_dm", 1:7), "eduage")]
UKB[, edu := ifelse(edu<7, NA, edu)]
UKB[, edu := ifelse(edu>=23, NA, edu)]

#################################################################3
# Merge EHR data

# 41280 Date of first in-patient diagnosis - ICD10
# 41281 Date of first in-patient diagnosis - ICD9

# 41270 Diagnoses - ICD10
# 41271 Diagnoses - ICD9

# 40005	Date of cancer diagnosis
# 40006	Type of cancer: ICD10
# 40013	Type of cancer: ICD9

# 40000	Date of death
# 40001	Underlying (primary) cause of death: ICD10
# 40002	Contributory (secondary) causes of death: ICD10

varlist <- names(fread("ukb_disease.tab", nrows=0))
varlist <- grep("eid|41280|41281|41270|41271|40005|40006|40013|40000|40001|40002", varlist, value=T)
to_add <- fread("ukb_disease.tab", select=varlist, colClasses=list(character=varlist) )
 
to_add[, f.eid := as.numeric(f.eid)]

UKB <- merge(UKB, to_add, by="f.eid")

#################
# transform cancer data 
cancer_comb <- function(x){
ICD10 <- x[names(x) %in% grep("40006", names(x), value=T)]
ICD9 <- x[names(x) %in% grep("40013", names(x), value=T)]

out <- rep(NA, 16)
out <- ifelse(is.na(ICD10[1:15]), ICD9[1:15], ICD10[1:15] )
out <- c(out, ICD10[16:17])
return(out)
} 
UKB[, paste0("cancer_", 0:16) := data.table(t(apply(.SD, 1, cancer_comb))) , .SDcols=grep("40006|40013", names(UKB), value=T) ]

######################################################
# phecode case/control
phemap <- data.table(PheWAS::phecode_map)
excl <- data.table(PheWAS::phecode_exclude)
phecode_list=data.table(PheWAS::pheinfo)$phecode

phemap[,  code := gsub("\\.", "", code)]

phecode <- function(x, D){
x <- trimws(x[!is.na(x)])
PHE <- phemap[match(x, phemap$code), phecode]
PHE <- ifelse(is.na(PHE), phemap[match(substr(x, 1,4), phemap$code), phecode], PHE)
PHE <- ifelse(is.na(PHE), phemap[match(substr(x, 1,3), phemap$code), phecode], PHE)
PHE <- unique(PHE[!is.na(PHE)])

sapply(D, function(d) ifelse( any(PHE == d, na.rm=T ), 1, 
                            ifelse(PHE %in% excl[code==d, exclusion_criteria], NA, 0)) )
}

library(future.apply)
plan(multicore, workers=23)

PHECODE <- UKB[, data.table((t(future_apply(.SD, 1, phecode, D=phecode_list)))), .SDcols=grep("41270|41271|cancer|40001|40002", names(UKB), value=T)]

UKB[, paste0("phecode_", phecode_list) := PHECODE]

####################################################

saveRDS(UKB, "UKBsib.Rds")
