
project_path=

# directory path for temporary files 
path=   

###########################################################################
# PREPARE INPUTS FOR MTAG
###########################################################################

for pheno in INDIVIDUAL HOUSEHOLD OCCUPATIONAL
do
(
awk ' NR==FNR {a[$3] = $10; next} ($3 in a) {print $3, a[$3]+$10}' OFS="\t"  \
${project_path}/META_ANALYSIS/OUTPUT/${pheno}_MEN/METAANALYSIS1.TBL \
${project_path}/META_ANALYSIS/OUTPUT/${pheno}_WOMEN/METAANALYSIS1.TBL \
 > ${path}/temp_${pheno} 

max_N=$(awk 'NR>1{print $2}' ${path}/temp_${pheno}  | datamash max 1)

# Drop SNPs with N < 0.5 * max N
awk -v N=$max_N ' ($2>(0.5*N)) {print $1}' OFS="\t"  \
 ${path}/temp_${pheno} >  ${path}/temp_${pheno}_SNPlist

    for G in MEN WOMEN
    do
        awk 'NR==FNR{a[$1]=$1; next} NR>1&&FNR==1{print "snpid","chr","bpos","a1","a2","freq","z","pval","n"}; (FNR>1 && $3 in a) {print $3,$1,$2,$4,$5,$6,$11,$12,$10}' OFS="\t" \
        ${path}/temp_${pheno}_SNPlist ${project_path}/META_ANALYSIS/OUTPUT/${pheno}_${G}/METAANALYSIS1.TBL \
        > ${path}/${pheno}_${G}_MTAG_FORMAT
    done

) 
done

### Parental income results are from a single cohort. Assume that the input files for parental income are already prepared in ${path}


###########################################################################
# RUN MATG to combine male + female
###########################################################################

for pheno in INDIVIDUAL HOUSEHOLD OCCUPATIONAL PARENTAL
do
(
python /home/hkweon/tools/mtag/mtag.py \
    --sumstats ${path}/${pheno}_WOMEN_MTAG_FORMAT,${path}/${pheno}_MEN_MTAG_FORMAT \
    --n_min 0 \
    --out ../OUTPUT/${pheno}.MTAG_META \
    --stream_stdout \
    --perfect_gencov \
    --equal_h \
    --force \
    --incld_ambig_snps
) 
done
