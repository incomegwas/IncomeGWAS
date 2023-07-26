
awk 'NR>1' SUMSTAT_LIST.csv | while read line
do 
file=$(echo ${line} | awk -F"," '{print $1}')
trait=$(echo ${line} | awk -F"," '{print $2}')

python /home/hkweon/tools/ldsc/ldsc.py \
--rg INC.sumstats.gz,${file} \
--ref-ld-chr /home/hkweon/tools/ldsc/eur_w_ld_chr/ \
--w-ld-chr /home/hkweon/tools/ldsc/eur_w_ld_chr/ \
--return-silly-things \
--print-delete-vals \
--out ../OUTPUT/INC_${trait} 

done