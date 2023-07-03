
project_path=

# directory path for temporary files 
path=   

###########################################################################
## prepare inputs for MTAG Meta-analysis (across four income measures)
###########################################################################

for pheno in INDIVIDUAL HOUSEHOLD OCCUPATIONAL PARENTAL
do 
    awk 'NR==1{OFS="\t"; print "snpid","chr","bpos","a1","a2","freq","z","pval","n"};
        NR>1 {OFS="\t"; print $1,$2,$3,$4,$5,$6,$7/$8,$9,$10}' \
        ../OUTPUT/${pheno}.MTAG_META.txt \
        > ${path}/${pheno}.MTAG_FORMAT 
done

###########################################################################
# Run META-ANALYSIS (4 versions)
###########################################################################

# META all
( python /home/hkweon/tools/mtag/mtag.py \
    --sumstats ${path}/INDIVIDUAL.MTAG_FORMAT,${path}/OCCUPATIONAL.MTAG_FORMAT,${path}/HOUSEHOLD.MTAG_FORMAT,${path}/PARENTAL.MTAG_FORMAT \
    --n_min 0 \
    --out ../OUTPUT/ALL_MEASURES \
    --stream_stdout \
    --perfect_gencov \
    --force \
    --incld_ambig_snps

) &

# Excluding PARENTAL
( python /home/hkweon/tools/mtag/mtag.py \
    --sumstats ${path}/INDIVIDUAL.MTAG_FORMAT,${path}/OCCUPATIONAL.MTAG_FORMAT,${path}/HOUSEHOLD.MTAG_FORMAT \
    --n_min 0 \
    --out ../OUTPUT/ALL_MEASURES.excl_PARENTAL \
    --stream_stdout \
    --perfect_gencov \
    --force \
    --incld_ambig_snps

) &

# Excluding Individual
( python /home/hkweon/tools/mtag/mtag.py \
    --sumstats ${path}/PARENTAL.MTAG_FORMAT,${path}/OCCUPATIONAL.MTAG_FORMAT,${path}/HOUSEHOLD.MTAG_FORMAT \
    --n_min 0 \
    --out ../OUTPUT/ALL_MEASURES.excl_INDIVIDUAL \
    --stream_stdout \
    --perfect_gencov \
    --force \
    --incld_ambig_snps

) &

# Excluding Individual & PARENTAL
( python /home/hkweon/tools/mtag/mtag.py \
    --sumstats ${path}/HOUSEHOLD.MTAG_FORMAT,${path}/OCCUPATIONAL.MTAG_FORMAT \
    --n_min 0 \
    --out ../OUTPUT/ALL_MEASURES.excl_INDIVIDUAL_PARENTAL \
    --stream_stdout \
    --perfect_gencov \
    --force \
    --incld_ambig_snps

) &

wait

mv ../OUTPUT/ALL_MEASURES_trait_2.txt ../OUTPUT/INCOME_GWAS_META_MTAG_all.txt
mv ../OUTPUT/ALL_MEASURES.excl_PARENTAL_trait_2.txt ../OUTPUT/INCOME_GWAS_META_MTAG_excl_PARENTAL.txt
mv ../OUTPUT/ALL_MEASURES.excl_INDIVIDUAL_trait_2.txt ../OUTPUT/INCOME_GWAS_META_MTAG_excl_INDIVIDUAL.txt
mv ../OUTPUT/ALL_MEASURES.excl_INDIVIDUAL_PARENTAL_trait_2.txt ../OUTPUT/INCOME_GWAS_META_MTAG_excl_INDIVIDUAL_PARENTAL.txt

Rscript ../CODE/COMBINE_MTAG.R









# awk 'BEGIN{OFS="\t"} NR==1 {print "rsID","chr","position","EFFECT_ALLELE","OTHER_ALLELE","EAF","BETA","SE","PVAL","N_EQUIV"};
#     NR>1 {print $1,$2,$3,$4,$5,$8,$9,$10,$12,int(1/($10*$10*2*$8*(1-$8)));}' ../OUTPUT/ALL_MEASURES_EXCL_UKBsibs.excl_PARENTAL_trait_2.txt > ../OUTPUT/ALL_MEASURES_EXCL_UKBsibs.excl_PARENTAL.MTAG_META.txt


# # check N of HM3 SNPs
# awk ' NR==FNR {a[$1] = $1; next} ($1 in a) {print $1}' OFS="\t"  \
# /home/hkweon/data/w_hm3.snplist ../OUTPUT/INCOME_GWAS_META_MTAG.txt  \
#   | wc -l  

# awk ' NR==FNR {a[$1] = $1; next} ($1 in a) {print $1}' OFS="\t"  \
# /home/hkweon/data/w_hm3.snplist ../OUTPUT/ALL_MEASURES.MTAG_META_RAISS_COMBINED.txt  \
#   | wc -l  

# awk ' NR==FNR {a[$1] = $1; next} ($1 in a) {print $1}' OFS="\t"  \
# /home/hkweon/data/w_hm3.snplist   ../OUTPUT/ALL_MEASURES_EXCL_UKBsibs_COMBINED_trait_2.txt  \
#   | wc -l  

