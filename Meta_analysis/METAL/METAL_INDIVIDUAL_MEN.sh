project_path=
cohort_results_path={project_path}/QC_cohorts/OUTPUT

mkdir {project_path}/META_ANALYSIS/OUTPUT/INDIVIDUAL_MEN
cd {project_path}/META_ANALYSIS/OUTPUT/INDIVIDUAL_MEN

metal_software_path=
{metal_software_path} << EOF
# Options
SCHEME SAMPLESIZE
MINMAXFREQ ON
AVERAGEFREQ ON
GENOMICCONTROL OFF
COLUMNCOUNTING LENIENT
USESTRAND OFF
TRACKPOSITIONS ON
# OVERLAP ON
## Removed overlap for now due to inconsistent results

## file options ##
SEPARATOR TAB
MARKER rsID
ALLELE EFFECT_ALLELE OTHER_ALLELE
EFFECT BETA
PVALUE PVAL
STDERR SE
FREQ EAF
WEIGHT N
CHROMOSOME chr
POSITION position

#files 
PROCESS ${cohort_results_path}/INCOME_GWAS_INDIVIDUAL_MEN_cohort1.txt
PROCESS ${cohort_results_path}/INCOME_GWAS_INDIVIDUAL_MEN_cohort2.txt
PROCESS ${cohort_results_path}/INCOME_GWAS_INDIVIDUAL_MEN_cohort3.txt

## OUTPUT

## Start META ANALYSIS
ANALYZE HETEROGENEITY

##Exit metal
QUIT

EOF


