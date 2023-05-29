#first identified the snps present in most old (90%) and mertensiana samples
vcftools --gzvcf TotalRawSNPs_rmhet.recode.vcf.gz --minQ 30 --remove-indels --keep samples_old_mer.txt --max-alleles 2 --min-alleles 2 --max-missing 0.9 --recode --stdout | bcftools 
query -f '%CHROM %POS\n' > snps_old_mer_90.txt
#then for that set of snps I calculated the missingness per individuals
vcftools --gzvcf TotalRawSNPs_rmhet.recode.vcf.gz --positions snps_old_mer_90.txt --missing-indv
#I removed individuals with > 15% missing data
vcftools --gzvcf TotalRawSNPs_rmhet.recode.vcf.gz --positions snps_old_mer_90.txt --keep indv_missing_15.txt --recode --stdout | bgzip -c > dataset_abba_baba.vcf.gz
#ABBA-BABA with Suite
Dsuite Dtrios -t tree_MNT.txt -c dataset_abba_baba.vcf.gz sets_MNT.txt


