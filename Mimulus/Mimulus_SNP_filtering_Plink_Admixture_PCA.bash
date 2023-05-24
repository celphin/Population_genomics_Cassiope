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

cd /home/celphin/projects/def-rieseber/celphin/Mimulus/Mimulus_timeseries

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

cd /home/celphin/projects/def-rieseber/celphin/Mimulus/
gunzip all_var_snp.vcf.gz 

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


################################
# same for both files filtering from here

#list the amount of missing data per individual - find indivdiuals with no reads mapped
vcftools --vcf Mimulus_filtered.recode.vcf --missing-indv --out Mimulus_filtered

#filter out the individuals with greater than 45% missing SNPs 
mawk '$5 > 0.45' Mimulus_filtered.imiss | cut -f1 > lowDP.indv
vcftools --vcf Mimulus_filtered.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out Mimulus_filtered_rm20
# After filtering, kept 364 out of 402 Individuals

# get list of all individuals left
vcftools --vcf Mimulus_filtered_rm20.recode.vcf --missing-indv --out Mimulus_filtered_rm20

# make Plink file to be able to run Admixture
mkdir Admixture_Mimulus
cp Mimulus_filtered_rm20.recode.vcf Admixture_Mimulus
cd Admixture_Mimulus

# convert the vcf file to Plink
vcftools --vcf Mimulus_filtered_rm20.recode.vcf --plink --out Mimulus_filtered_rm20

#MAKE A BED FILE 
plink --file Mimulus_filtered_rm20 --allow-no-sex --allow-extra-chr 0 --make-bed --out Mimulus_filtered_rm20


################################
# Mimulus RUN ADMIXTURE
# Unsupervised analysis K from 1 to 9 and 5 replicates

tmux new-session -s Mimulus
tmux attach-session -t Mimulus

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
              "./Mimulus_filtered_rm20.recode.gds",
              method="biallelic.only")
