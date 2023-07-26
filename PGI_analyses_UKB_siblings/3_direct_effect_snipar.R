### This script estimates the direct effect of income PGI
library(fixest)

# Load data
UKB <- readRDS("UKBsib.Rds")

# Merge PGI data
PGI = fread("/INC_SNIPAR.pgs.txt", 
                    col.names=c("FID", "IID", "INC_PGI", "INC_PGI_M", "INC_PGI_F"))

i="INC"
sd = sd(PGI[, get(paste0(i, "_PGI"))])   
COLS = paste0(i, c("_PGI", "_PGI_M", "_PGI_P"))
PGI[, (COLS) := lapply(.SD, function(x) x / sd), .SDcols=COLS ]     

UKB[, IID := paste0(n_eid, "_", n_eid)]
UKB = merge(UKB, PGI, by="IID")

UKB[, INC_PGI_PAR := INC_PGI_P + INC_PGI_M]

# transform houseld income
UKB$log_hinc <- UKB[, dplyr::case_when(hinc == 1 ~ log(3/4*18000),
                                       hinc == 2 ~ log((18000+30999) / 2), 
                                       hinc == 3 ~ log((31000+51999) / 2),
                                       hinc == 4 ~ log((52000+100000) / 2),
                                       hinc == 5 ~ log(4/3*100000))]


############################################
# assortative mating correction
get_factor = function(x, D){
FID_IN = D[, .N, by=FID][N==2, FID]
r = D[FID %in% FID_IN, .(get(x)[1], get(x)[2]), by = FID][, cor(V1, V2)]
r_am = 2*r -1
return(r_am)
}

r_am = get_factor("INC_PGI", UKB) # 1.1066098
##############################################3

pc <- paste0("pc", 1:20)
pc <- paste0(pc, collapse="+")

Y=c("logyh_hourly", "log_hinc")
Y_out=c("Occupational wage", "Household income")

FM=c(paste0(" ~ male*(factor(year_hourly) + age_hourly + I(age_hourly^2) + I(age_hourly^3)) + geno + ", pc),
     paste0(" ~ male*(age + I(age^2) + I(age^3)+ factor(year)) + geno + ", pc))

output = NULL
for (i in 1:2){

UKB$Ystd = NULL
UKB[!is.na(get(Y[i])) , Ystd := lm(paste0(Y[i], FM[i]), UKB )$residuals  ]
UKB[!is.na(get(Y[i])) & male==1 , Ystd := scale(Ystd)  ]
UKB[!is.na(get(Y[i])) & male==0 , Ystd := scale(Ystd)  ]

temp <- NULL
temp <- rbind(temp, 
          coeftable(fixest::feols(as.formula(paste0("Ystd ~ INC_PGI ")), data=UKB, se="cluster", cluster="familyID"))[2,])

LM = fixest::feols(as.formula(paste0("Ystd ~ INC_PGI + INC_PGI_PAR")), data=UKB, se="cluster", cluster="familyID")

temp <- rbind(temp, 
          coeftable(LM)[2:3,])
temp <- data.table(temp)

names(temp) <- c("EST", "SE", "T", "P")
POP_COR = car::deltaMethod(LM, "INC_PGI + 1.1066098* INC_PGI_PAR")
temp$POP_COR = POP_COR[1]
temp$POP_COR_SE = POP_COR[2]

RATIO = car::deltaMethod(LM, "INC_PGI / (INC_PGI + 1.1066098*INC_PGI_PAR)")
temp$RATIO = RATIO[1]
temp$RATIO_SE = RATIO[2]

temp$Y = Y_out[i]
temp$N = LM$nobs
temp$M = c("POP", "D", "AVG_NTC")
print(temp)
output <- rbind(output, temp)
}

fwrite(output, "./INC_PGS_COEF_SNIPAR.txt", sep="\t")

