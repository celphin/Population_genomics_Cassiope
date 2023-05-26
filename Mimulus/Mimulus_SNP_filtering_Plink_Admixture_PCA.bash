#####################################
# Mimulus analysis
# May 2023
# SNP filtering, plink conversion, Admixture
#####################################
# File contains:

# SNP filtering (vcftools)
## quality, minor allele counts, % found accross individuals
## filter out individuals with lots of missing data
## A script to count the number of potential genotyping errors due to low read depth

# Plink
## make a plink file from vcf
## make .bed file from Plink

# Admixture
## run admixture for different k values with loop
## bootstrap admixture results

###########################
# General SNP filtering starts

# Mimulus

#load vcftools
module load StdEnv/2020 
module load vcftools/0.1.16
module load gnuplot/5.4.2
module load gcc/9.3.0
module load bcftools/1.16
module load plink/1.9b_6.21-x86_64

# rename samples in vcf file
#https://www.biostars.org/p/279195/ 

cd /home/celphin/projects/def-rieseber/Dryas_shared_data/Mimulus/Mimulus_timeseries

bcftools query -l Mimulus_timeseries_filtered_variants.vcf > sample_names.txt

sed 's/,/ /g' M_caridnalis_samples_renamed.csv > M_caridnalis_samples_renamed.txt
bcftools reheader -s M_caridnalis_samples_renamed.txt -o Mimulus_timeseries_filtered_variants_rename.vcf Mimulus_timeseries_filtered_variants.vcf

bcftools query -l Mimulus_timeseries_filtered_variants_rename.vcf

#-----------------------------
# all individuals
# filter for quality, indels, biallelic, missing in less than 90%, monomorphic, LD, and minor allele freq of 0.01
vcftools --vcf Mimulus_timeseries_filtered_variants_rename.vcf \
--minGQ 30 \
--remove-indels \
--max-alleles 2 \
--max-missing 0.60 \
--min-alleles 2 \
--thin 1000 \
--mac 4 \
--recode --recode-INFO-all --out Mimulus_filtered
# After filtering, kept 402 out of 402 Individuals
# After filtering, kept 151 out of a possible 4995632 Sites # 0.95
# After filtering, kept 2161 out of a possible 4995632 Sites #0.70
# After filtering, kept 4622 out of a possible 4995632 Sites #0.60


####################################
#change for other file here: 
#/project/6019339/shared_data/Mimulus_Beagle/
#all_vqsr_filtered_variants.vcf.gz
#all_vqsr_filtered_variants.vcf.gz.tbi

cd /home/celphin/projects/def-rieseber/Dryas_shared_data/Mimulus/

gunzip all_var_snp.vcf.gz 

bcftools query -l all_var_snp.vcf > sample_names.txt

# edit to make rename file in Excel

sed 's/,/ /g' M_caridnalis_samples_renamed.csv > M_caridnalis_samples_renamed.txt
bcftools reheader -s M_caridnalis_samples_renamed.txt -o all_var_snp_rename.vcf all_var_snp.vcf

bcftools query -l all_var_snp_rename.vcf

#-----------------------------
# all individuals
# filter for quality, indels, biallelic, missing in less than 90%, monomorphic, LD, and minor allele freq of 0.01
vcftools --vcf all_var_snp_rename.vcf \
--minGQ 30 \
--remove-indels \
--max-alleles 2 \
--max-missing 0.60 \
--min-alleles 2 \
--thin 1000 \
--mac 4 \
--recode --recode-INFO-all --out Mimulus_filtered
# After filtering, kept 665 out of 665 Individuals
# After filtering, kept 14444 out of a possible 29928667 Sites


################################
# same for both files filtering from here

#list the amount of missing data per individual - find indivdiuals with no reads mapped
vcftools --vcf Mimulus_filtered.recode.vcf --missing-indv --out Mimulus_filtered

#filter out the individuals with greater than 45% missing SNPs 
mawk '$5 > 0.45' Mimulus_filtered.imiss | cut -f1 > lowDP.indv
vcftools --vcf Mimulus_filtered.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out Mimulus_filtered_rm20
# After filtering, kept 364 out of 402 Individuals
# After filtering, kept 576 out of 665 Individuals


# get list of all individuals left
vcftools --vcf Mimulus_filtered_rm20.recode.vcf --missing-indv --out Mimulus_filtered_rm20

#############################
# Filter vcf to include only baseline

grep -E '2006|2007|2008|2009|2010' Mimulus_filtered_rm20.imiss | cut -f1  > baseline_indv.txt

vcftools --vcf Mimulus_filtered_rm20.recode.vcf --keep baseline_indv.txt --recode --recode-INFO-all --out Mimulus_filtered_baseline
# After filtering, kept 285 out of 576 Individuals

vcftools --vcf Mimulus_filtered_baseline.recode.vcf --missing-indv --out Mimulus_filtered_baseline

################################
# Mimulus RUN ADMIXTURE

# Baseline

tmux new-session -s Mimulus_base
tmux attach-session -t Mimulus_base

# make Plink file to be able to run Admixture
mkdir Admixture_Mimulus_Baseline
cp Mimulus_filtered_baseline.recode.vcf Admixture_Mimulus_Baseline
cd Admixture_Mimulus_Baseline

# convert the vcf file to Plink
vcftools --vcf Mimulus_filtered_baseline.recode.vcf --plink --out Mimulus_filtered_baseline

#MAKE A BED FILE 
plink --file Mimulus_filtered_baseline --allow-no-sex --allow-extra-chr 0 --make-bed --out Mimulus_filtered_baseline

# Unsupervised analysis K from 1 to 9 and 5 replicates

# get allocation
salloc -c48 --time 02:50:00 --mem 187G --account rpp-rieseber

module load StdEnv/2020
module load admixture/1.3.0

# run loop through various K values
mkdir Take1

cp Mimulus_filtered_baseline* Take1

cd Take1
for K in 1 2 3 4 5 6 7 8 9 10; \
do admixture --cv=10 -s time -j48 -C 0.0000000001  Mimulus_filtered_baseline.bed $K | tee log${K}.out; done

#to get the CV errors and see which K value is the best model
grep -h CV ./Take*/log*out

CV error (K=1): 0.17135
CV error (K=2): 0.13176
CV error (K=3): 0.12380
CV error (K=4): 0.12098
CV error (K=5): 0.12077
CV error (K=6): 0.12008
CV error (K=7): 0.11967
CV error (K=8): 0.11962
CV error (K=9): 0.12030
CV error (K=10): 0.12241

########################
# Total

tmux new-session -s Mimulus
tmux attach-session -t Mimulus

# make Plink file to be able to run Admixture
mkdir Admixture_Mimulus_Total
cp Mimulus_filtered_rm20.recode.vcf Admixture_Mimulus_Total
cd Admixture_Mimulus_Total

# convert the vcf file to Plink
vcftools --vcf Mimulus_filtered_rm20.recode.vcf --plink --out Mimulus_filtered_rm20

#MAKE A BED FILE 
plink --file Mimulus_filtered_rm20 --allow-no-sex --allow-extra-chr 0 --make-bed --out Mimulus_filtered_rm20

# Unsupervised analysis K from 1 to 9 and 5 replicates

# get allocation
salloc -c48 --time 02:50:00 --mem 187G --account rpp-rieseber

module load StdEnv/2020
module load admixture/1.3.0

# run loop through various K values
mkdir Take1

cp Mimulus_filtered_rm20* Take1

cd Take1
for K in 1 2 3 4 5 6 7 8 9 10; \
do admixture --cv=10 -s time -j48 -C 0.0000000001  Mimulus_filtered_rm20.bed $K | tee log${K}.out; done

#to get the CV errors and see which K value is the best model
grep -h CV ./Take*/log*out

CV error (K=1): 0.16769
CV error (K=2): 0.12637 4000
CV error (K=3): 0.11564 1000
CV error (K=4): 0.11086 500
CV error (K=5): 0.10867
CV error (K=6): 0.10696
CV error (K=7): 0.10500
CV error (K=8): 0.10580
CV error (K=9): 0.10580
CV error (K=10): 0.10537

#############################
# Timeseries
#to get the CV errors and see which K value is the best model
grep -h CV ./Take*/log*out

# Timeseries
CV error (K=1): 0.23570 
CV error (K=2): 0.17157 6000
CV error (K=3): 0.15471 2000
*CV error (K=4): 0.14469 1000
CV error (K=5): 0.13991 500
CV error (K=6): 0.13803 100
CV error (K=7): 0.13766
CV error (K=8): 0.13391
CV error (K=9): 0.13296
CV error (K=10): 0.13251


######################################
# Make PCA GDS file

# copy over PCA file
cd ..

module load nixpkgs/16.09
module load gcc/7.3.0
module load r/3.6.0
module load gdal
module load udunits
module load python
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/3.6/

R

library(SNPRelate)

snpgdsVCF2GDS("./Mimulus_filtered_rm20.recode.vcf",
              "./Mimulus_filtered_total.recode.gds",
              method="biallelic.only")

# Number of samples: 576
# import 14444 variants.


snpgdsVCF2GDS("./Mimulus_filtered_baseline.recode.vcf",
              "./Mimulus_filtered_baseline.recode.gds",
              method="biallelic.only")
# Number of samples: 285
# import 14444 variants.
