# By Hyeokmoon Kweon, h.kweon@vu.nl

######################################################################################
# function PGS_LM: Estimate PGS coefficient with or without family fixed effects

    # Note: the function depends on "fixest", "data.table" packages

# Inputs:
    # fm: formula object to specify a model. Do not include a PGS here, but only covariates. 
    # PGS: character object of PGS name
    # data: data.frame or data.table 
    # FID: character object of family ID column name, Default to NULL. If provided, SE will be clustered at family. Has to be supplied for within-family estimation.  
    # resid: logical, indicating whether the target outcome should be residualized and standardized first. Default to TRUE
    # within: logical, indicating within-family estimation. Default to FALSE

# Output:
    # fixest oject from "fixest" package
     # see https://www.rdocumentation.org/packages/fixest/versions/0.8.4/topics/feols
     
PGS_LM <- function(fm, PGS, data, FID=NULL, resid=TRUE, std=TRUE, within=FALSE, between=FALSE, d_r2=FALSE, d_r2_fm1=NULL, d_r2_fm2=NULL){
fixest::setFixest_nthreads(flexiblas::flexiblas_get_num_threads(), save = FALSE)

library(data.table)        
data <- data.table(data)
data <- na.omit(data[, .SD, .SDcols=unique(c(all.vars(fm), FID, PGS))])
if (within==TRUE ) {
        data <- data[get(FID) %in% data[, .N, by=get(FID)][N>=2][[1]] ]}

if (between==TRUE & !is.null(FID)) {
        set.seed(123)
        data[, IID := paste0("ID_", 1:nrow(data))]
        data <- data[IID %in% data[, sample(IID[1:.N],1), by=FID][[2]] ]
        # print(data)
        }

if (resid==TRUE){
        LM = lm(fm, data)
        data[, Y := LM$residuals]
        cat("Covar R2 = ", summary(LM)$r.squared, "\n")
        if (std == TRUE) data[, Y := scale(Y)]

        if (within==TRUE){
                FM = paste0("Y ~ PGS | ", FID)
        } else {
                FM = "Y ~ PGS"
        }

} else if (resid==FALSE){
        fm = as.character(fm)        
        print(fm)
        if (within==TRUE){
                FM = paste0(fm[2], "~ PGS + ", fm[3], "| ", FID)
        } else {
                FM = paste0(fm[2], "~ PGS + ", fm[3])
                print(FM)
        }
} else if (d_r2==TRUE){
        fm = as.character(fm)        
        
}

# standardize observed PGS
data[, PGS := scale(get(PGS))]

if (within){
        return(fixest::feols(as.formula(FM), data=data, se="cluster", cluster=FID))
} else if (!is.null(FID) & between==FALSE ){
        print(FM)
        return(fixest::feols(as.formula(FM), data=data, se="cluster", cluster=FID))
} else {
        return(fixest::feols(as.formula(FM), data=data, se="hetero"))
}
}




######## Example ########
# UKB <- fread("../TEMP/UKB_temp.csv")

# pc <- paste0("pc", 1:20)
# pc <- paste0(pc, collapse="+")
# fm = as.formula(paste0("edu ~ male*factor(yob) +geno+ ", pc))

# PGS_LM(fm, "EA_PGS", data=UKB, FID="familyID")
# PGS_LM(fm, "EA_PGS", data=UKB, FID="familyID", resid=FALSE)

# PGS_LM(fm, "EA_PGS", data=UKB, FID="familyID", within=TRUE)
# PGS_LM(fm, "EA_PGS", data=UKB, FID="familyID", within=TRUE, resid=FALSE)



PGS_LM_dr2 <- function(d_r2_fm1=NULL, d_r2_fm2=NULL, PGS, data, CI=FALSE, unrelated=FALSE, FID=NULL){
library(data.table)        
data <- data.table(data)
data <- na.omit(data[, .SD, .SDcols=unique(c(all.vars(d_r2_fm1), all.vars(d_r2_fm2), FID, PGS))])

if (unrelated==TRUE & !is.null(FID)) {
        set.seed(123)
        data[, IID := paste0("ID_", 1:nrow(data))]
        data <- data[IID %in% data[, sample(IID[1:.N],1), by=FID][[2]] ]
        # print(data)
        }

if (!is.null(d_r2_fm1) & !is.null(d_r2_fm2)){
        LM = lm(d_r2_fm1, data)
        data[, Y := LM$residuals]
        cat("Covar R2 = ", summary(LM)$r.squared, "\n")

        d_r2_fm2 = as.character(d_r2_fm2)        
        FM1 = paste0("Y ~       ", d_r2_fm2[length(d_r2_fm2)])
        FM2 = paste0("Y ~ ", paste0(PGS, collapse="+"), " + ", d_r2_fm2[length(d_r2_fm2)])
} else if (!is.null(d_r2_fm1) & is.null(d_r2_fm2)){
        d_r2_fm1 = as.character(d_r2_fm1)        

        FM1 = paste0(d_r2_fm1[2], "~       ", d_r2_fm1[3])
        FM2 = paste0(d_r2_fm1[2], "~ ", paste0(PGS, collapse="+"), " + ", d_r2_fm1[3])
}


# standardize observed PGS
# data[, PGS := scale(get(PGS))]

LM_res = s(lm(as.formula(FM2), data))
COEF = coef(LM_res)
print(COEF)
dr2 <- LM_res$r.squared - s(lm(as.formula(FM1), data))$r.squared
return(list(d_r2 = dr2, N = nrow(data), COEF = COEF[2,]))

}

