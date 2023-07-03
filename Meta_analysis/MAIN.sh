############################################################################
# 1. Run meta-analysis on the cohort results for each sex & each income phenotype
############################################################################

for pheno in INDIVIDUAL HOUSEHOLD OCCUPATIONAL
do
    for G in MEN WOMEN
    do
        bash ./METAL/METAL_${pheno}_${G}.sh
    done
done


############################################################################
# 2. Run MTAG to combine male + female results for each phenotype
############################################################################

bash MTAG_COMBINE_SEX.sh


############################################################################
# 3. Run MTAG to conduct meta-analysis across four income measures 
############################################################################

bash MTAG_COMBINE_SEX.sh




