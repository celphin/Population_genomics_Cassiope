#####################################
# Entire Cassiope analysis - part 7
# March 2023
# TreeMix Plotting in R
#####################################
# https://github.com/carolindahms/TreeMix
# https://speciationgenomics.github.io/Treemix/

######################
# Copy files over to local computer
# plot in R

setwd("~/GitHub/Population_genomics_Cassiope/Figures_data/TreeMix/TreeMix_replicates") # of course this needs to be adjusted
prefix="TreeMix.0.rep"

library(RColorBrewer)
library(R.utils)

# download or find this file
source("../plotting_funcs.R") # here you need to add the path

# make 10 plots in one
png("../../Plots/Treemix_Total_reps.png",width=2400,height=2000,res = 300)
par(mfrow=c(2,3))
for(edge in 0:10){
  plot_tree(cex=0.8,paste0(prefix,edge))
  title(paste(edge,"edges"))
}
dev.off()

# one plot per edge
for(edge in 0:10){
  png(paste0("../../Plots/Treemix_rep", as.character(edge), ".png"),width=2400,height=2000,res = 300)
  plot_tree(cex=0.8,paste0(prefix,edge))
  title(paste(edge,"rep"))
  dev.off()
}

# Residuals
pop.list <-  c("Outgroup", "Saximontana", "Beringia", "NWT", "Nunavut", "Greenland", "Europe", "Russia")

for(edge in 0:5){
 plot_resid(stem=paste0(prefix,"_",edge), pop_order="pop.list")
}


##########################
# Fastsimcoal2 model - tetragona
# only one model with no edges

setwd("~/GitHub/Rieseberg_Lab/Population_genomics_Cassiope/Figures_data/TreeMix/TreeMix_tet/") # of course this needs to be adjusted
prefix="TreeMix_tet"

library(RColorBrewer)
library(R.utils)

# download or find this file
source("../plotting_funcs.R") # here you need to add the path

# make 5 plots in one
edge=0
plot_tree(cex=0.8,paste0(prefix,".",edge))
title(paste(edge,"edges tet"))

# one plot per edge
png(paste0("../../Plots/Treemix_tet_", as.character(edge), ".png"),width=2400,height=2000,res = 300)
plot_tree(cex=0.8,paste0(prefix,".",edge))
title(paste(edge,"edges tet"))
dev.off()



##########################
setwd("~/GitHub/Rieseberg_Lab/Population_genomics_Cassiope/Figures_data/TreeMix/TreeMix_mer/") # of course this needs to be adjusted
prefix="TreeMix_mer"

library(RColorBrewer)
library(R.utils)

# download or find this file
source("../plotting_funcs.R") # here you need to add the path

# make 5 plots in one
par(mfrow=c(2,3))
for(edge in 0:5){
  plot_tree(cex=0.8,paste0(prefix,".",edge))
  title(paste(edge,"edges mer"))
}

# one plot per edge
for(edge in 0:5){
  png(paste0("../../Plots/Treemix_mer_", as.character(edge), ".png"),width=2400,height=2000,res = 300)
  plot_tree(cex=0.8,paste0(prefix,".",edge))
  title(paste(edge,"edges mer"))
  dev.off()
}


##########################




























############################

library(RColorBrewer)
library(R.utils)
source("/scratch/celphin/GBS_Cassiope/Mar2020_SNPFiltering_PopStats/TreeMix/plotting_funcs.R") 

for(edge in 0:5){
png(paste0("../Plots/Treemix_", as.character(edge), ".png"),width=2400,height=2000,res = 300)
  plot_tree(cex=0.8,paste0(prefix,"_",edge))
  title(paste(edge,"edges"))
dev.off()
}

png("Treemix_Resid.png",width=2400,height=2000,res=300)
for(edge in 0:5){
 plot_resid(stem=paste0(prefix,".",edge), pop_order="pop.list")
}
dev.off()



module load nixpkgs/16.09
module load gcc/7.3.0
module load r/3.6.0
module load gdal
module load udunits
module load python
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/3.6/

R

setwd("/scratch/celphin/Mar2020_SNPFiltering_PopStats/Trees/TreeMix/") 
prefix="FiltXXg9mac5minq30r60i_chrom_rLD.noN"

library(RColorBrewer)
library(R.utils)
source("/scratch/celphin/Mar2020_SNPFiltering_PopStats/Trees/TreeMix/plotting_funcs.R") 

for(edge in 0:5){
png(paste0("Treemix_v2_", as.character(edge), "plot.png"),width=2400,height=2000,res = 300)
  plot_tree(cex=0.8,paste0(prefix,".",edge))
  title(paste(edge,"edges"))
dev.off()
}

png("Treemix_Resid.png",width=2400,height=2000,res=300)
for(edge in 0:5){
 plot_resid(stem=paste0(prefix,".",edge), pop_order="pop.list")
}
dev.off()


##################################
#from Armando
#initial analysis were done with BT 4 samples as the outgroup
#I then added LKTRT as the outgroup and considered BT as two groups of 2 samples
#LKTRT has lots of missing data. Also try one third analysis with BRKTRT instead of LKTRT
#BT interior and BT coastal
#in Flex /home/geraldes/fish/AC_DV/AC_DV_SNPs do
mkdir acdvbtlt184
cd acdvbtlt184
#remove BTBR and 3 bad DV
cat > BTBR_lowcov
BR_UN_UN_102.d.bam
DV_AK_AF_274.d.bam
DV_AK_SS_278.d.bam
DV_AK_SS_280.d.bam
/Linux/bin/vcftools --vcf ../../acdv_align/acdv188.vcf \
--remove BTBR_lowcov \
--remove-indels \
--min-alleles 2 \
--max-alleles 2 \
--recode \
--out acdvbtlt184
#After filtering, kept 184 out of 188 Individuals
#After filtering, kept 3,151,729 out of a possible 4,283,884 Sites
#use the new version of Greg Owen’s script ../../extras/vcf2maxhet.pl
#to eliminate SNPs with Observed Heterozygosity 0.6 or above
cat acdvbtlt184.recode.vcf | perl ../../extras/vcf2maxhet.pl 0.6 > acdvbtlt184_het06_snps.vcf
#3,104,923 printed sites.
#21,711 cut due to heterozygosity.
#25,095 cut due to no data
#remove 30%miss, GQ<10 and singletons
/Linux/bin/vcftools --vcf acdvbtlt184_het06_snps.vcf \
--mac 2 \
--minGQ 10 \
--max-missing 0.7 \
--recode \
--out acdvbtlt184_het06_snps_miss07_mac2_GQ10
#After filtering, kept 184 out of 184 Individuals
#After filtering, kept 614,783 out of a possible 3104923 Sites
#########GENERATE PLINK FILES
#generate ped and map files to use in plink
/Linux/bin/vcftools --vcf acdvbtlt184_het06_snps_miss07_mac2_GQ10.recode.vcf \
--plink \
--out acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink
#VCFTOOLS makes the chromosome 0 and for extimating LD, PLINK does not like that.
#copy to my computer at /Users/ageraldes/Documents/UBC/RickTaylor/AC_DV/acdvbtlt184/
scp geraldes@files.zoology.ubc.ca:flex/fish/AC_DV/AC_DV_SNPs/acdvbtlt184/acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink.map .
#Fix the file in excel so that the contig shows up in the chromosome
#copy to flex
scp acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink.map geraldes@files.zoology.ubc.ca:flex/fish/AC_DV/AC_DV_SNPs/acdvbtlt184/acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink.map
#MAKE A BED FILE
/Linux/bin/plink_1.9 \
--file acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink \
--allow-no-sex \
--allow-extra-chr \
--make-bed \
--out acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink

###############################
#614,783 variants loaded from .bim file.
#184 people (0 males, 0 females, 183 ambiguous) loaded from .fam.
#create a data184.clust file where the first two columns are the first two columns of the fam file
#and the third column is the POP

/Linux/bin/plink_1.9 \
--bfile acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink \
--allow-no-sex \
--allow-extra-chr \
--freq \
--missing \
--within data184.clust \
--out acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink_data184clust

#614,783 variants loaded from .bim file.
#184 people (0 males, 0 females, 183 ambiguous) loaded from .fam.
#Total genotyping rate is 0.883256.
#--freq: Cluster-stratified allele frequencies written to
#acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink_data184clust.frq.strat .
#--missing: Individual missing data report written to
#acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink_data184clust.imiss, and
#variant-based cluster-stratified missing data report written to
#acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink_data184clust.lmiss.

gzip acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink_data184clust.frq.strat
python ../../tools/plink2treemix.py acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink_data184clust.frq.strat.gz acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink_data184clust.frq.strat_treemix.gz
#this creates the treemix input file

#Treemix requires that we tell it where the root is.
#take LD into account with option -k
#k is the number of SNPs to group into blocks to take LD into account
#tell treemix that LT is the root

/Linux/bin/treemix -i acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink_data184clust.frq.strat_treemix.gz -root LT -k 50 -o acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink_data184clust.frq.strat_treemix_k50_LT_stem

#Allow migration
/Linux/bin/treemix -i acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink_data184clust.frq.strat_treemix.gz -root LT -k 50 -m 1 -o acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink_data184clust.frq.strat_treemix_k50_LT_m1_stem
/Linux/bin/treemix -i acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink_data184clust.frq.strat_treemix.gz -root LT -k 50 -m 2 -o acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink_data184clust.frq.strat_treemix_k50_LT_m2_stem
/Linux/bin/treemix -i acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink_data184clust.frq.strat_treemix.gz -root LT -k 50 -m 3 -o acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink_data184clust.frq.strat_treemix_k50_LT_m3_stem
/Linux/bin/treemix -i acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink_data184clust.frq.strat_treemix.gz -root LT -k 50 -m 4 -o acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink_data184clust.frq.strat_treemix_k50_LT_m4_stem
/Linux/bin/treemix -i acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink_data184clust.frq.strat_treemix.gz -root LT -k 50 -m 5 -o acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink_data184clust.frq.strat_treemix_k50_LT_m5_stem
/Linux/bin/treemix -i acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink_data184clust.frq.strat_treemix.gz -root LT -k 50 -m 6 -o acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink_data184clust.frq.strat_treemix_k50_LT_m6_stem
/Linux/bin/treemix -i acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink_data184clust.frq.strat_treemix.gz -root LT -k 50 -m 7 -o acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink_data184clust.frq.strat_treemix_k50_LT_m7_stem
/Linux/bin/treemix -i acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink_data184clust.frq.strat_treemix.gz -root LT -k 50 -m 8 -o acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink_data184clust.frq.strat_treemix_k50_LT_m8_stem
/Linux/bin/treemix -i acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink_data184clust.frq.strat_treemix.gz -root LT -k 50 -m 9 -o acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink_data184clust.frq.strat_treemix_k50_LT_m9_stem
/Linux/bin/treemix -i acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink_data184clust.frq.strat_treemix.gz -root LT -k 50 -m 10 -o acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink_data184clust.frq.strat_treemix_k50_LT_m10_stem

#copy to my computer at /Users/ageraldes/Documents/UBC/RickTaylor/AC_DV/acdvbtlt184/
scp geraldes@files.zoology.ubc.ca:flex/fish/AC_DV/AC_DV_SNPs/acdvbtlt184/acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink_data184clust.frq.strat_treemix_k50_LT*stem* .

#Finally, perform the threepop and fourpop tests
/Linux/bin/threepop -i acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink_data184clust.frq.strat_treemix.gz -k 50 2>&1 | tee acdvbtlt184_3pop.log
/Linux/bin/fourpop -i acdvbtlt184_het06_snps_miss07_mac2_GQ10_plink_data184clust.frq.strat_treemix.gz -k 50 2>&1 | tee acdvbtlt184_4pop.log

#############NOTE: I just figured out that we can bootstrap the tree Treemix generates. but to do that we cannot use the Treemix way of accointing for LD
#I’d like to redo all of the analysis by trimming on LD in plink first and then doing all I’ve done so far but by running bootstraps
#See https://speciationgenomics.github.io/Treemix/

https://speciationgenomics.github.io/Treemix/

#################################

############################
#TreeMix
# use Admixture groups

cd /scratch/celphin/GBS_Cassiope/Mar2020_SNPFiltering_PopStats/TreeMix/

tmux new-session -s Treemix
tmux attach-session -t Treemix

salloc -c1 --time 03:00:00 --mem 120000m --account def-rieseber

# move old files to V2 folder
#mv *.* ./TreeMix_V2

# move old files to FiltXX folder
#mv FiltXX* ./FiltXX/
#mv *.png ./FiltXX/

# get files to run on:
# run on Ancient2 

cp  /scratch/celphin/GBS_Cassiope/Mar2020_SNPFiltering_PopStats/SNP_filtering/Ancient_test/Ancient_test2.* .
#cp  /scratch/celphin/GBS_Cassiope/Mar2020_SNPFiltering_PopStats/SNP_filtering/FiltXXg9mac5minq30r60i_chrom_rLD.* .

#FILE= FiltXXg9mac5minq30r60i_chrom_rLD
FILE=Ancient_test2

module load vcftools
vcftools --vcf $FILE.vcf --max-missing 1 --recode --stdout | gzip > $FILE.noN.vcf.gz

module load bcftools
bcftools query -l $FILE.noN.vcf.gz | awk '{split($1,pop,"."); print $1"\t"$1"\t"pop[2]}' > $FILE.noN.clust

#always need to save clust file based on backup???

## download file to the computer
## save population group next to each indiv in Excel
## use Admixture groups

[[1]]Europe - 45
"LatnjajaureFieldStation" 
"LongyearbyenSvalbard"    
"PaddusSweden"            
"PetuniaSvalbard"        
"ZackenbergGreenland"    

[[2]] N Nunavut - 53
"AlexandraFiord"    
"AlexandraFiordOld" 
"BylotIslandPos1"   
"BylotIslandPos2"   
"BylotIslandPos3"  
"GriseFiord"        
"LakeHazenNunavut"  
"SverdrupPass"     

[[3]] Russia -19
"SamoylovRussia" 
"YedonaRussia"  

[[4]] BC -10
"GentianPeak"

[[5]] XX -0
"DiskoLyngmarksfjeldGreenland"

[[6]] Outgroup - 10
"HarrisonHut" 
"PhelixHut"  

[[7]]Alaska-Sagwon
"ImnavaitCreek"     
"SagwonAlaksa"      
"ToolikMoistAcidic"

[[8]] Yukon -alpine -59
"12MileMeadow"     
"AtqasukDryAlaska" 
"BarrowDryAlaska"  
"DenaliAlaska"     
"KluanePlateau"    
"PikaCamp"        

[[9]] S Nunavut 68
"ClydeRiver"                 
"DiskoQeqertarsuaqGreenland" 
"IlulissatGreenland"         
"IqaluitNunavut"            
"ItilleqFiordGreenland"      
"QikiqtarjuaqNunavut"       

[[10]] NWT - 49
"AxelHeibergIsland"       
"EurekaNunavut"           
"FosheimPeninsula"         
"KugluktukNunavut"        
"PearcePoint"              
"QikiqtarukHerschelIsland" 
"YAMBAbeach"              

# edit line 53 of vcf2treemix to /scratch/celphin/GBS_Cassiope/Mar2020_SNPFiltering_PopStats/TreeMix/plink2treemix.py
# done

# rename .bed, .fam, .bim files

#FILE=FiltXXg9mac5minq30r60i_chrom_rLD.noN
FILE=Ancient_test2.noN


#gunzip -k $FILE.vcf.gz

###########################################

#gunzip -k $FILE.treemix.frq.gz

#salloc -c10 --time 01:00:00 --mem 120000m --account def-rieseber

module load nixpkgs/16.09 
module load gcc/7.3.0
module load treemix/1.13
module load python
module load intel/2018.3
module load vcftools/0.1.16
module load boost/1.68.0 
module load gsl/2.5 
module load plink/1.9b_5.2-x86_64

./vcf2treemix.sh $FILE.vcf $FILE.clust

###########################
# upload each grouping as .clust file

./vcf2treemix.sh $FILE.vcf $FILE.AKYk.clust

for i in {0..5}
do
 treemix -i $FILE.treemix.frq.gz -m $i -o $FILE.AKYk.$i -root PopGEN -bootstrap -k 50 > treemix_${i}_log &
done

###########################

./vcf2treemix.sh $FILE.vcf $FILE.allpop.clust

for i in {0..5}
do
 treemix -i $FILE.treemix.frq.gz -m $i -o $FILE.allpop.$i -root PopGEN -bootstrap -k 50 > treemix_${i}_log &
done

###########################
./vcf2treemix.sh $FILE.vcf $FILE.Ellesmere.clust

for i in {0..5}
do
 treemix -i $FILE.treemix.frq.gz -m $i -o $FILE.Ellesmere.$i -root PopGEN -bootstrap -k 50 > treemix_${i}_log &
done

###########################
./vcf2treemix.sh $FILE.vcf $FILE.Baffin.clust


for i in {0..5}
do
 treemix -i $FILE.treemix.frq.gz -m $i -o $FILE.Baffin.$i -root PopGEN -bootstrap -k 50 > treemix_${i}_log &
done

###########################

./vcf2treemix.sh $FILE.vcf $FILE.Europe.clust


for i in {0..5}
do
 treemix -i $FILE.treemix.frq.gz -m $i -o $FILE.Europe.$i -root PopGEN -bootstrap -k 50 > treemix_${i}_log &
done

###########################

./vcf2treemix.sh $FILE.vcf $FILE.NWT.clust

for i in {0..5}
do
 treemix -i $FILE.treemix.frq.gz -m $i -o $FILE.NWT.$i -root PopGEN -bootstrap -k 50 > treemix_${i}_log &
done



###########################################
# TreeMix plots

cd /scratch/celphin/GBS_Cassiope/Mar2020_SNPFiltering_PopStats/TreeMix/
module load nixpkgs/16.09
module load gcc/7.3.0
module load r/3.6.0
module load gdal
module load udunits
module load python
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/3.6/

R

setwd("/scratch/celphin/GBS_Cassiope/Mar2020_SNPFiltering_PopStats/TreeMix/") 
prefix="FiltXXg9mac5minq30r60i_chrom_rLD.noN.NWT"

library(RColorBrewer)
library(R.utils)
source("/scratch/celphin/GBS_Cassiope/Mar2020_SNPFiltering_PopStats/TreeMix/plotting_funcs.R") 

for(edge in 0:5){
png(paste0("Treemix_NWT_", as.character(edge), "plot.png"),width=2400,height=2000,res = 300)
  plot_tree(cex=0.8,paste0(prefix,".",edge))
  title(paste(edge,"edges"))
dev.off()
}

#-----------------------------
png("Treemix_Resid_allpop.png",width=2400,height=2000,res=300)
for(edge in 0:5){
 plot_resid(stem=paste0(prefix,".",edge), pop_order="pop.list")
}
dev.off()


###########################