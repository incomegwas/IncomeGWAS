library(data.table)

LIST <- fread("SUMSTAT_LIST.csv")
output_INC <- data.table(trait = LIST$trait)

for (i in LIST$trait) {
    file <- paste0("../OUTPUT/rg_INC_", i, ".log")
    RES <- readLines(file)
    RES <- grep("Genetic Correlation:", RES, fixed = TRUE, value = TRUE)
    RES <- stringr::str_split_fixed(RES, " ", 4)
    rg <- as.numeric(RES[, 3])
    se <- as.numeric(gsub("\\(|\\)", "", RES[, 4]))

    Z_RES <- grep("Z-score: ", readLines(file), fixed = TRUE, value = TRUE)
    Z_RES <- as.numeric(stringr::str_split_fixed(Z_RES, " ", 2)[1, 2])

    output_INC[trait == i, EST := rg]
    output_INC[trait == i, SE := se]
    output_INC[trait == i, Z := Z_RES]
}

output_INC[, P := pnorm(-abs(Z)) * 2]
output_INC[, P_FDR := p.adjust(P, "fdr")]