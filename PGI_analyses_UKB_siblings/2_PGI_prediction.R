library(fixest)
source("PGI_LM_functions.R")

# Load data
UKB <- readRDS("UKBsib.Rds")

# merge PGI data
PGI <- fread("EA_plink_score.profile", select=c(2,6), col.names=c("bgen_ID", "EA_PGI"))
UKB <- merge(UKB, PGI, by="bgen_ID", all.x=T, sort=F)

PGI <- fread("INC_plink_score.profile", select=c(2,6), col.names=c("bgen_ID", "INC_PGI"))
UKB <- merge(UKB, PGI, by="bgen_ID", all.x=T, sort=F)

# Standardize
UKB[, grep("PGI", names(UKB), value=TRUE) := lapply(.SD, scale), .SDcols=grep("PGI", names(UKB), value=TRUE)  ]

# Transform household income
UKB$log_hinc <- UKB[, dplyr::case_when(hinc == 1 ~ log(3/4*18000),
                                       hinc == 2 ~ log((18000+30999) / 2), 
                                       hinc == 3 ~ log((31000+51999) / 2),
                                       hinc == 4 ~ log((52000+100000) / 2),
                                       hinc == 5 ~ log(4/3*100000))]


#####################
# Get prediction

pc <- paste0("pc", 1:20)
pc <- paste0(pc, collapse="+")
fm2 <- as.formula(paste0("~  geno + ", pc))

######################
Y=c("logyh_hourly", "log_hinc", "edu")
Y_out=c("Occupational wage", "Household income", "EA")
FM=c("~   male*(factor(year_hourly) + age_hourly + I(age_hourly^2) + I(age_hourly^3) )",
     "~   male*( age + I(age^2) + I(age^3)+factor(year))",
     "~   male*(factor(yob))")

P = c("INC", "EA")
P_out = c("Income", "EA")

res = data.table()
for (i in 1:length(Y)){
fm <- as.formula(paste0(Y[i], FM[i]))
set.seed(111)
IN_IID = UKB[!is.na(get(Y[i])), n_eid[sample(1:.N, 1)], by=familyID][[2]]
cat(Y, length(IN_IID), "\n")
    for (p in 1:length(P)){
        LM <- PGI_LM_dr2(fm, fm2, paste0(P[p], "_PGI"), data=UKB[n_eid %in% IN_IID])
        res <- rbind(res, data.table(T=Y_out[i], PGI=P_out[p], R2=LM$d_r2))
}

print(Y[i])
}
#######################

fwrite(res, "UKBsib_PGI_prediction.txt")

########################3
# bootstrap

bs_func <- function(s){
source("PGI_LM_functions.R")

res = NULL
for (i in 1:length(Y)){
set.seed(111)
IN_IID = UKB[!is.na(get(Y[i])), n_eid[sample(1:.N, 1)], by=familyID][[2]]

set.seed(s); IN_IID = sample(IN_IID, length(IN_IID), replace=TRUE)

fm <- as.formula(paste0(Y[i], FM[i]))
    for (p in 1:length(P)){
         res <- rbind(res 
         ,c(PGI_LM_dr2(fm, fm2, paste0(P[p], "_PGI"), data=UKB[n_eid %in% IN_IID])$d_r2))
          }
}
return(res)
}

BS <- pbmcapply::pbmclapply(1:1000, bs_func, mc.cores=32, ignore.interactive=TRUE)
saveRDS(BS, "boostrap.Rds")



