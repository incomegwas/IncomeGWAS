# Script to combine the MTAG meta-analysis from the 4 versions. 

D1 = fread("../OUTPUT/INCOME_GWAS_META_MTAG_all.txt")
D2 = fread("../OUTPUT/INCOME_GWAS_META_MTAG_excl_INDIVIDUAL_PARETNAL.txt")
D3 = fread("../OUTPUT/INCOME_GWAS_META_MTAG_excl_INDIVIDUAL.txt")
D4 = fread("../OUTPUT/INCOME_GWAS_META_MTAG_excl_PARETNAL.txt")

# I(ndividual), O(ccupational), H(ousehold), P(parental)
D1$INCL = "IHOP"
D2$INCL = "OH"
D3$INCL = "OHP"
D4$INCL = "IOH"

# For each SNP, find the one with the largest chi-square
LIST = unique(c(D1$SNP, D2$SNP, D3$SNP, D4$SNP))
Z = cbind(D1[match(LIST, SNP), Z], D2[match(LIST, SNP), Z], D3[match(LIST, SNP), Z], D4[match(LIST, SNP), Z])
IND = apply(Z, 1, function(x) which.max(x^2))

output = rbind(D1[match(LIST[IND==1], SNP)], D2[match(LIST[IND==2], SNP)], D3[match(LIST[IND==3], SNP)], D4[match(LIST[IND==4], SNP)])

output = output[, .(SNP, CHR, BP, A1, A2, FRQ, mtag_beta, mtag_se, mtag_z, mtag_pval, INCL)]
setnames(output, c("FRQ", "mtag_beta", "mtag_se", "mtag_z", "mtag_pval"), c("EAF", "BETA", "SE", "Z", "P"))

# effective sample size per SNP
output[, N := 1/(SE^2 * 2 * EAF * (1-EAF))]

# Drop SNPs if effective N < 70% of max N
output = output[N>=0.7*max(N)]

fwrite(output, "../OUTPUT/INCOME_GWAS_META.txt", sep="\t")