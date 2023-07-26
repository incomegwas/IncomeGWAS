# LDpred2 wrapper

__Usage__:

bash RUN_MAIN.sh ${ref_map} ${LD} ${g_path} ${sumstat} ${format} ${output}


__Inputs__:
  - ref_map: the plink bim file for the reference panel 
  - LD: LD estimates from the reference panel
  - g_path: the plink file prefix for the prediction sample
  - sumstat: GWAS summary statistics file path
  - format: format of the GWAS summary statistics (easyqc, mtag_meta, mtag, ldpred)
  - output: prefix for the outputs
  