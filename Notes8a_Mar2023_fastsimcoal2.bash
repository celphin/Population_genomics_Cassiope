###########################################################
# Entire Cassiope analysis - part 8 - Demography
# March 2023
# try fastsimcoal2 to test a bunch of different demographic scenarios
#####################################
# Information

# http://cmpg.unibe.ch/software/fastsimcoal27/index.html
# https://www.pnas.org/doi/epdf/10.1073/pnas.1600405113 
# https://speciationgenomics.github.io/fastsimcoal2/
# file:///C:/Users/Owner/Downloads/fastsimcoal27.pdf
# https://evomics.org/learning/population-and-speciation-genomics/2022-population-and-speciation-genomics/inferring-demography/

###################################
# data to use

cd ~/scratch/Cassiope/
mkdir Fastsimcoal2_Mar2023
cd ~/scratch/Cassiope/Fastsimcoal2_Mar2023

# copy over vcf and obs files to use for models from globus/Fernando
cp -v ~/scratch/Cassiope/SNP_filtering_March2023/Cassiope_r30i.recode.vcf ~/scratch/Cassiope/Fastsimcoal2_Mar2023


cp -v ~/scratch/Cassiope/SNP_filtering_March2023/Cassiope_noMER_r10i.imiss ~/scratch/Cassiope/Fastsimcoal2_Mar2023/pop_lists
cp -v ~/scratch/Cassiope/Admixture_March2023/Take1/Cassiope_noMER_r10i.5.Q ~/scratch/Cassiope/Fastsimcoal2_Mar2023/pop_lists


####################################
# make population files (see bottom of this file)
cd  ~/scratch/Cassiope/Fastsimcoal2_Mar2023/pop_lists

# remove header
tail -n +2 Cassiope_noMER_r10i.imiss > Cassiope_noMER_r10i.txt
paste Cassiope_noMER_r10i.txt  Cassiope_noMER_r10i.5.Q > Admix_data_5.txt

#--------------------
# columns 6-10 of Admix_data_5.txt have Admix data
# select individuals  > 0.95 of that group

# column 6 (V1) is Russia
mawk '$6 > 0.95' Admix_data_5.txt | cut -f1 > Russia.txt

# column 7 (V2) is Europe
mawk '$7 > 0.95' Admix_data_5.txt | cut -f1 > Europe.txt

# column 8 (V3) is Saximontana
mawk '$8 > 0.95' Admix_data_5.txt | cut -f1 > BCtet.txt

# column 9 (V4) is Greenland
mawk '$9 > 0.95' Admix_data_5.txt | cut -f1 > Greenland.txt

# column 10 (V5) is Alaska
mawk '$10 > 0.95' Admix_data_5.txt | cut -f1 > Alaska.txt

# Nunavut 
touch Nunavut0.txt
grep "AlexNew" Admix_data_5.txt | cut -f1 > tmp.txt
cat Nunavut0.txt tmp.txt > Nunavut1.txt
grep "EUR" Admix_data_5.txt | cut -f1 > tmp.txt
cat Nunavut1.txt tmp.txt > Nunavut2.txt
grep "FOS" Admix_data_5.txt | cut -f1 > tmp.txt
cat Nunavut2.txt tmp.txt > Nunavut3.txt
grep "HAZ" Admix_data_5.txt | cut -f1 > tmp.txt
cat Nunavut3.txt tmp.txt > Nunavut4.txt
grep "GF" Admix_data_5.txt | cut -f1 > tmp.txt
cat Nunavut4.txt tmp.txt > Nunavut5.txt
more Nunavut5.txt

# NWT
touch NWT0.txt
grep "QHI" Admix_data_5.txt | cut -f1 > tmp.txt
cat NWT0.txt tmp.txt > NWT1.txt
grep "PEA" Admix_data_5.txt | cut -f1 > tmp.txt
cat NWT1.txt tmp.txt > NWT2.txt
grep "KUQ" Admix_data_5.txt | cut -f1 > tmp.txt
cat NWT2.txt tmp.txt > NWT3.txt
grep "YAM" Admix_data_5.txt | cut -f1 > tmp.txt
cat NWT3.txt tmp.txt > NWT4.txt
grep "AXE" Admix_data_5.txt | cut -f1 > tmp.txt
cat NWT4.txt tmp.txt > NWT5.txt
more NWT5.txt


# Baffin
touch Baffin0.txt
grep "CR" Admix_data_5.txt | cut -f1 > tmp.txt
cat Baffin0.txt tmp.txt > Baffin1.txt
grep "IQ" Admix_data_5.txt | cut -f1 > tmp.txt
cat Baffin1.txt tmp.txt > Baffin2.txt
grep "Kik" Admix_data_5.txt | cut -f1 > tmp.txt
cat Baffin2.txt tmp.txt > Baffin3.txt
grep "025" Admix_data_5.txt | cut -f1 > tmp.txt
cat Baffin3.txt tmp.txt > Baffin4.txt
grep "IG" Admix_data_5.txt | cut -f1 > tmp.txt
cat Baffin4.txt tmp.txt > Baffin5.txt
more Baffin5.txt

#----------------------
# add column with group name
wc -l Alaska.txt
yes "Alaska" | head -n 63 > Reg.txt
paste Alaska.txt Reg.txt > Alaska_list

wc -l Europe.txt
yes "Europe" | head -n 24 > Reg.txt
paste Europe.txt Reg.txt > Europe_list

wc -l Russia.txt
yes "Russia" | head -n 18 > Reg.txt
paste Russia.txt Reg.txt > Russia_list

wc -l Greenland.txt
yes "Greenland" | head -n 16 > Reg.txt
paste Greenland.txt Reg.txt > Greenland_list

wc -l Nunavut5.txt
yes "Nunavut" | head -n 34 > Reg.txt
paste Nunavut5.txt Reg.txt > Nunavut_list

wc -l NWT5.txt
yes "NWT" | head -n 41 > Reg.txt
paste NWT5.txt Reg.txt > NWT_list

wc -l Baffin5.txt
yes "Greenland" | head -n 38 > Reg.txt
paste Baffin5.txt Reg.txt > Baffin_list

wc -l BCtet.txt
yes "BCtet" | head -n 9 > Reg.txt
paste BCtet.txt Reg.txt > BCtet_list

#---------------------
# add outgroup
# BCmer
cat << EOF > BCmer_list
PopHAR_13 BCmer
PopHAR_1 BCmer
PopHAR_21 BCmer
PopHAR_8 BCmer
PopPHE_1 BCmer
PopPHE_3 BCmer
EOF

sed -i 's/ /\t/'  BCmer_list

#---------------------
# join into population lists
cd  ~/scratch/Cassiope/Fastsimcoal2_Mar2023/pop_lists

# no Greenland
cat Russia_list Alaska_list Europe_list BCmer_list BCtet_list  > list_pops_5pop_mer.txt
cat Russia_list Alaska_list Europe_list NWT_list Nunavut_list > list_pops_5pop_tet.txt

# with Greenland
cat Russia_list Alaska_list Europe_list BCmer_list Greenland_list > list_pops_6pop_mer.txt # removes Saximontana and add Greenland
cat Greenland_list Alaska_list Europe_list NWT_list Nunavut_list > list_pops_6pop_tet.txt # does not include hybrid individuals in Baffin grouped with Greenland

mv list_pops_* ..
cd ..

##############################
# https://github.com/isaacovercast/easySFS

#run easySFS to convert vcf to SFS 
git clone https://github.com/isaacovercast/easySFS.git

cd easySFS
chmod 777 easySFS.py

cd  ~/scratch/Cassiope/Fastsimcoal2_Mar2023/

module load StdEnv/2020 
module load python/3.8.10
module load scipy-stack/2022a

python3 ./easySFS/easySFS.py -i Cassiope_r30i.recode.vcf -p list_pops_5pop_mer.txt --proj 20,20,20,12,18 -o 5pop_mer
python3 ./easySFS/easySFS.py -i Cassiope_r30i.recode.vcf -p list_pops_5pop_tet.txt --proj 20,20,20,20,20 -o 5pop_tet
python3 ./easySFS/easySFS.py -i Cassiope_r30i.recode.vcf -p list_pops_6pop_mer.txt --proj 20,20,20,12,20 -o 6pop_mer # note cannot do 6 populations at once so remove Saximontana
python3 ./easySFS/easySFS.py -i Cassiope_r30i.recode.vcf -p list_pops_6pop_tet.txt --proj 20,20,20,20,20 -o 6pop_tet # note cannot do 6 populations so remove Russia to keep Greenland

python3 ./easySFS/easySFS.py -i Cassiope_r30i.recode.vcf -p list_pops_5pop_mer.txt --proj 20,20,20,12,18 -o 5pop_mer_eur
python3 ./easySFS/easySFS.py -i Cassiope_r30i.recode.vcf -p list_pops_5pop_tet.txt --proj 20,20,20,20,20 -o 5pop_tet_eur

python3 ./easySFS/easySFS.py -i Cassiope_r30i.recode.vcf -p list_pops_6pop_mer.txt --proj 20,20,20,12,20 -o 6pop_mer_eur # note cannot do 6 populations at once so remove Saximontana
python3 ./easySFS/easySFS.py -i Cassiope_r30i.recode.vcf -p list_pops_6pop_tet.txt --proj 20,20,20,20,20 -o 6pop_tet_eur # note cannot do 6 populations so remove Russia to keep Greenland

# Processing 5 populations - ['Russia', 'Alaska', 'Europe', 'BCmer', 'BCtet']
# Processing 5 populations - ['Russia', 'Alaska', 'Europe', 'NWT', 'Nunavut']

# Processing 5 populations - ['Russia', 'Alaska', 'Europe', 'BCmer', 'Greenland']
# Processing 5 populations - ['Greenland', 'Alaska', 'Europe', 'NWT', 'Nunavut']


# Continue, excluding samples not in both pops file and VCF? (yes/no)
# yes
# Sampling one snp per locus (CHROM)
# SFS files written to /scratch/celphin/Cassiope/Fastsimcoal2_Mar2023/output

#-------------------
# copy over files
cp -v ./5pop_mer/fastsimcoal2/Cassiope_r30i_MSFS.obs ./5pops_mertensiana_MSFS.obs
cp -v ./5pop_tet/fastsimcoal2/Cassiope_r30i_MSFS.obs ./5pops_tet_MSFS.obs
cp -v ./6pop_mer/fastsimcoal2/Cassiope_r30i_MSFS.obs ./6pops_mertensiana_MSFS.obs
cp -v ./6pop_tet/fastsimcoal2/Cassiope_r30i_MSFS.obs ./6pops_tet_MSFS.obs
cp -v ./5pop_mer_eur/fastsimcoal2/Cassiope_r30i_MSFS.obs ./5pops_mer_eur_MSFS.obs
cp -v ./5pop_tet_eur/fastsimcoal2/Cassiope_r30i_MSFS.obs ./5pops_tet_eur_MSFS.obs

cp -v ./6pop_mer_eur/fastsimcoal2/Cassiope_r30i_MSFS.obs ./6pops_mer_eur_MSFS.obs
cp -v ./6pop_tet_eur/fastsimcoal2/Cassiope_r30i_MSFS.obs ./6pops_tet_eur_MSFS.obs


#----------------
# double check projection
head -2 *obs

==> 5pops_mertensiana_MSFS.obs <==
1 observations. No. of demes and sample sizes are on next line.
5       20 20 20 12 18

==> 5pops_tet_MSFS.obs <==
1 observations. No. of demes and sample sizes are on next line.
5       20 20 20 20 20

==> 6pops_mertensiana_MSFS.obs <==
1 observations. No. of demes and sample sizes are on next line.
5       20 20 20 12 20

==> 6pops_tet_MSFS.obs <==
1 observations. No. of demes and sample sizes are on next line.
5       20 20 20 20 20


############################################
# setup 5 pop mertensiana model

cat << EOF > 5pops_mertensiana.tpl
//Number of population samples (demes)
5 populations to simulate
//Population effective sizes (number of genes)
N_Russia
N_Alaska
N_Europe
N_mertensiana
N_saximontana
//Sample sizes
20
20
20
12
18
//Growth rates
0
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
0
//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index 
4 historical event
TDIVEurAla 2 1 1 1 0 0
TDIVRusAla 0 1 1 1 0 0
TDIVSaxTet 4 1 1 1 0 0
1860000 3 1 1 N_Anc 0 0 absoluteResize
//Number of independent loci [chromosomes]
1 0
//Per chromosome: Number of linkage blocks
1
//per block: Datatype, numm loci, rec rate and mut rate + optional parameters
FREQ 1 0 1e-7 OUTEXP

EOF

cat << EOF >  5pops_mertensiana.est
// Priors and rules file
// *********************

[PARAMETERS]
//#isInt? #name   #dist.#min  #max 
//all Ns are in number of haploid individuals
1  N_Europe       unif     10    1e6   output
1  N_Alaska       unif     10    1e6   output
1  N_Russia       unif     10    1e6   output
1  N_mertensiana       unif     10    1e6   output
1  N_saximontana       unif     10    1e6   output
1  N_Anc       unif     10    1e6   output
1  TDIVEurAla     unif     10    1e6   output 
1  TDIVRusAla     unif     10    1e6   output 
1  TDIVSaxTet     unif     10    1e6   output 

EOF

############################################
# setup 5 pop tetragona model

cat << EOF > 5pops_tet.tpl
//Number of population samples (demes)
5 populations to simulate
//Population effective sizes (number of genes)
N_Russia
N_Alaska
N_Europe
N_NWT
N_Nunavut
//Sample sizes
20
20
20
20
20
//Growth rates
0
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
2
//Migration matrix 0
0 0 0 0 0 
0 0 0 0 0
0 0 0 0 0
0 0 0 0 migNWTNun
0 0 0 migNunNWT 0
//Migration matrix 1
0 0 0 0 0 
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index 
4 historical event
700 3 1 1 1 0 1
900 4 2 1 1 0 1
TDIVEurAla 2 1 1 1 0 1
TDIVRusAla 0 1 1 1 0 1
//Number of independent loci [chromosomes]
1 0
//Per chromosome: Number of linkage blocks
1
//per block: Datatype, numm loci, rec rate and mut rate + optional parameters
FREQ 1 0 1e-7 OUTEXP

EOF

cat << EOF >  5pops_tet.est
// Priors and rules file
// *********************

[PARAMETERS]
//#isInt? #name   #dist.#min  #max 
//all Ns are in number of haploid individuals
1  N_Europe       unif     10    1e6   output
1  N_Alaska       unif     10    1e6   output
1  N_Russia       unif     10    1e6   output
1  N_NWT       unif     10    1e6   output
1  N_Nunavut       unif     10    1e6   output
1  TDIVEurAla     unif     10    1e6   output 
1  TDIVRusAla     unif     10    1e6   output
0  migNWTNun     unif     0    0.5   output 
0  migNunNWT     unif     0    0.5   output

EOF


############################################
# setup 6 pop mertensiana model

cat << EOF > 6pops_mertensiana.tpl
//Number of population samples (demes)
5 populations to simulate
//Population effective sizes (number of genes)
N_Russia
N_Alaska
N_Europe
N_mertensiana
N_Greenland
//Sample sizes
20
20
20
12
20
//Growth rates
0
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
0
//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index 
4 historical event
TDIVGreEur 4 2 1 1 0 0
TDIVEurAla 2 1 1 1 0 0
TDIVRusAla 0 1 1 1 0 0
1860000 3 1 1 N_Anc 0 0 absoluteResize
//Number of independent loci [chromosomes]
1 0
//Per chromosome: Number of linkage blocks
1
//per block: Datatype, numm loci, rec rate and mut rate + optional parameters
FREQ 1 0 1e-7 OUTEXP

EOF

cat << EOF >  6pops_mertensiana.est
// Priors and rules file
// *********************

[PARAMETERS]
//#isInt? #name   #dist.#min  #max 
//all Ns are in number of haploid individuals
1  N_Europe       unif     10    1e6   output
1  N_Alaska       unif     10    1e6   output
1  N_Russia       unif     10    1e6   output
1  N_mertensiana       unif     10    1e6   output
1  N_Greenland       unif     10    1e6   output
1  N_Anc       unif     10    1e6   output
1  TDIVGreEur     unif     10    1e6   output 
1  TDIVEurAla     unif     10    1e6   output 
1  TDIVRusAla     unif     10    1e6   output 

EOF

############################################
# setup 6 pop tetragona model

cat << EOF > 6pops_tet.tpl
//Number of population samples (demes)
5 populations to simulate
//Population effective sizes (number of genes)
N_Greenland
N_Alaska
N_Europe
N_NWT
N_Nunavut
//Sample sizes
20
20
20
20
20
//Growth rates
0
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
2
//Migration matrix 0
0 0 0 0 0 
0 0 0 0 0
0 0 0 0 0
0 0 0 0 migNWTNun
0 0 0 migNunNWT 0
//Migration matrix 1
0 0 0 0 0 
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index 
4 historical event
700 3 1 1 1 0 1
900 4 2 1 1 0 1
TDIVGreEur 0 2 1 1 0 1
TDIVEurAla 2 1 1 1 0 1
//Number of independent loci [chromosomes]
1 0
//Per chromosome: Number of linkage blocks
1
//per block: Datatype, numm loci, rec rate and mut rate + optional parameters
FREQ 1 0 1e-7 OUTEXP

EOF

cat << EOF >  6pops_tet.est
// Priors and rules file
// *********************

[PARAMETERS]
//#isInt? #name   #dist.#min  #max 
//all Ns are in number of haploid individuals
1  N_Europe       unif     10    1e6   output
1  N_Alaska       unif     10    1e6   output
1  N_Greenland       unif     10    1e6   output
1  N_NWT       unif     10    1e6   output
1  N_Nunavut       unif     10    1e6   output
1  TDIVEurAla     unif     10    1e6   output 
1  TDIVGreEur     unif     10    1e6   output
0  migNWTNun     unif     0    0.5   output 
0  migNunNWT     unif     0    0.5   output

EOF



############################################
# setup 5 pop mertensiana Europe model - split Russia from Sax, Eur from Russia, Alaska from Eur

cat << EOF > 5pops_mer_eur.tpl
//Number of population samples (demes)
5 populations to simulate
//Population effective sizes (number of genes)
N_Russia
N_Alaska
N_Europe
N_mertensiana
N_saximontana
//Sample sizes
20
20
20
12
18
//Growth rates
0
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
0
//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index 
4 historical event
TDIVEurAla 1 2 1 1 0 0
TDIVRusEur 2 0 1 1 0 0
TDIVSaxTet 4 0 1 1 0 0
1860000 3 0 1 N_Anc 0 0 absoluteResize
//Number of independent loci [chromosomes]
1 0
//Per chromosome: Number of linkage blocks
1
//per block: Datatype, numm loci, rec rate and mut rate + optional parameters
FREQ 1 0 1e-7 OUTEXP

EOF

cat << EOF >  5pops_mer_eur.est
// Priors and rules file
// *********************

[PARAMETERS]
//#isInt? #name   #dist.#min  #max 
//all Ns are in number of haploid individuals
1  N_Europe       unif     10    1e6   output
1  N_Alaska       unif     10    1e6   output
1  N_Russia       unif     10    1e6   output
1  N_mertensiana       unif     10    1e6   output
1  N_saximontana       unif     10    1e6   output
1  N_Anc       unif     10    1e6   output
1  TDIVEurAla     unif     10    1e6   output 
1  TDIVRusEur     unif     10    1e6   output 
1  TDIVSaxTet     unif     10    1e6   output 

EOF

############################################
# setup 5 pop tetragona  Europe model

cat << EOF > 5pops_tet_eur.tpl
//Number of population samples (demes)
5 populations to simulate
//Population effective sizes (number of genes)
N_Russia
N_Alaska
N_Europe
N_NWT
N_Nunavut
//Sample sizes
20
20
20
20
20
//Growth rates
0
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
2
//Migration matrix 0
0 0 0 0 0 
0 0 0 0 0
0 0 0 0 0
0 0 0 0 migNWTNun
0 0 0 migNunNWT 0
//Migration matrix 1
0 0 0 0 0 
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index 
4 historical event
700 3 1 1 1 0 1
900 4 2 1 1 0 1
TDIVEurAla 1 2 1 1 0 1
TDIVRusEur 2 0 1 1 0 1
//Number of independent loci [chromosomes]
1 0
//Per chromosome: Number of linkage blocks
1
//per block: Datatype, numm loci, rec rate and mut rate + optional parameters
FREQ 1 0 1e-7 OUTEXP

EOF

cat << EOF >  5pops_tet_eur.est
// Priors and rules file
// *********************

[PARAMETERS]
//#isInt? #name   #dist.#min  #max 
//all Ns are in number of haploid individuals
1  N_Europe       unif     10    1e6   output
1  N_Alaska       unif     10    1e6   output
1  N_Russia       unif     10    1e6   output
1  N_NWT       unif     10    1e6   output
1  N_Nunavut       unif     10    1e6   output
1  TDIVEurAla     unif     10    1e6   output 
1  TDIVRusEur     unif     10    1e6   output
0  migNWTNun     unif     0    0.5   output 
0  migNunNWT     unif     0    0.5   output

EOF


############################################
# setup 6 pop mertensiana model _ Europe
cat << EOF > 6pops_mer_eur.tpl
//Number of population samples (demes)
5 populations to simulate
//Population effective sizes (number of genes)
N_Russia
N_Alaska
N_Europe
N_mertensiana
N_Greenland
//Sample sizes
20
20
20
12
20
//Growth rates
0
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
0
//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index 
4 historical event
TDIVGreAla 4 1 1 1 0 0
TDIVEurAla 2 1 1 1 0 0
TDIVRusAla 0 1 1 1 0 0
1860000 3 1 1 N_Anc 0 0 absoluteResize
//Number of independent loci [chromosomes]
1 0
//Per chromosome: Number of linkage blocks
1
//per block: Datatype, numm loci, rec rate and mut rate + optional parameters
FREQ 1 0 1e-7 OUTEXP

EOF

cat << EOF >  6pops_mer_eur.est
// Priors and rules file
// *********************

[PARAMETERS]
//#isInt? #name   #dist.#min  #max 
//all Ns are in number of haploid individuals
1  N_Europe       unif     10    1e6   output
1  N_Alaska       unif     10    1e6   output
1  N_Russia       unif     10    1e6   output
1  N_mertensiana       unif     10    1e6   output
1  N_Greenland       unif     10    1e6   output
1  N_Anc       unif     10    1e6   output
1  TDIVGreAla     unif     10    1e6   output 
1  TDIVEurAla     unif     10    1e6   output 
1  TDIVRusAla     unif     10    1e6   output 

EOF

############################################
# setup 6 pop tetragona model - Europe

cat << EOF > 6pops_tet_eur.tpl
//Number of population samples (demes)
5 populations to simulate
//Population effective sizes (number of genes)
N_Greenland
N_Alaska
N_Europe
N_NWT
N_Nunavut
//Sample sizes
20
20
20
20
20
//Growth rates
0
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
2
//Migration matrix 0
0 0 0 0 0 
0 0 0 0 0
0 0 0 0 0
0 0 0 0 migNWTNun
0 0 0 migNunNWT 0
//Migration matrix 1
0 0 0 0 0 
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index 
4 historical event
700 3 1 1 1 0 1
900 4 2 1 1 0 1
TDIVGreAla 0 2 1 1 0 1
TDIVEurAla 2 1 1 1 0 1
//Number of independent loci [chromosomes]
1 0
//Per chromosome: Number of linkage blocks
1
//per block: Datatype, numm loci, rec rate and mut rate + optional parameters
FREQ 1 0 1e-7 OUTEXP

EOF

cat << EOF >  6pops_tet_eur.est
// Priors and rules file
// *********************

[PARAMETERS]
//#isInt? #name   #dist.#min  #max 
//all Ns are in number of haploid individuals
1  N_Europe       unif     10    1e6   output
1  N_Alaska       unif     10    1e6   output
1  N_Greenland       unif     10    1e6   output
1  N_NWT       unif     10    1e6   output
1  N_Nunavut       unif     10    1e6   output
1  TDIVEurAla     unif     10    1e6   output 
1  TDIVGreAla     unif     10    1e6   output
0  migNWTNun     unif     0    0.5   output 
0  migNunNWT     unif     0    0.5   output

EOF


##################################
# make 50 folders with files

#mkdir runs; cd runs
cd ~/scratch/Cassiope/Fastsimcoal2_Mar2023/runs

# 5 Mertensiana model
for i in {1..50}
do
mkdir 5pops_mertensiana_$i
cp ../5pops_mertensiana_MSFS.obs 5pops_mertensiana_$i
cp ../5pops_mertensiana.tpl 5pops_mertensiana_$i
cp ../5pops_mertensiana.est 5pops_mertensiana_$i
done

# 5 tet model
for i in {1..50}
do
mkdir 5pops_tet_$i
cp ../5pops_tet_MSFS.obs 5pops_tet_$i
cp ../5pops_tet.tpl 5pops_tet_$i
cp ../5pops_tet.est 5pops_tet_$i
done

# 6 Mertensiana model
for i in {1..50}
do
mkdir 6pops_mertensiana_$i
cp ../6pops_mertensiana_MSFS.obs 6pops_mertensiana_$i
cp ../6pops_mertensiana.tpl 6pops_mertensiana_$i
cp ../6pops_mertensiana.est 6pops_mertensiana_$i
done

# 6 tet model
for i in {1..50}
do
mkdir 6pops_tet_$i
cp ../6pops_tet_MSFS.obs 6pops_tet_$i
cp ../6pops_tet.tpl 6pops_tet_$i
cp ../6pops_tet.est 6pops_tet_$i
done

# 5 Mertensiana Europe model
for i in {1..50}
do
mkdir 5pops_mer_eur_$i
cp ../5pops_mer_eur_MSFS.obs 5pops_mer_eur_$i
cp ../5pops_mer_eur.tpl 5pops_mer_eur_$i
cp ../5pops_mer_eur.est 5pops_mer_eur_$i
done

# 5 tet Europe model
for i in {1..50}
do
mkdir 5pops_tet_eur_$i
cp ../5pops_tet_eur_MSFS.obs 5pops_tet_eur_$i
cp ../5pops_tet_eur.tpl 5pops_tet_eur_$i
cp ../5pops_tet_eur.est 5pops_tet_eur_$i
done

# 6 Mer model Europe
for i in {1..50}
do
mkdir 6pops_mer_eur_$i
cp ../6pops_mer_eur_MSFS.obs 6pops_mer_eur_$i
cp ../6pops_mer_eur.tpl 6pops_mer_eur_$i
cp ../6pops_mer_eur.est 6pops_mer_eur_$i
done

# 6 tet model Europe
for i in {1..50}
do
mkdir 6pops_tet_eur_$i
cp ../6pops_tet_eur_MSFS.obs 6pops_tet_eur_$i
cp ../6pops_tet_eur.tpl 6pops_tet_eur_$i
cp ../6pops_tet_eur.est 6pops_tet_eur_$i
done

#########################################
# try to run one iteration of mertensiana model

tmux new-session -s fastsimcoal
tmux attach-session -t fastsimcoal

salloc -c48 --time 2:50:00 --mem-per-cpu=2G --account rpp-rieseber

cd ~/scratch/Cassiope/Fastsimcoal2_Mar2023

module load StdEnv/2020
module load fastsimcoal2/2.7.0.9

fsc27 -t 5pops_mertensiana.tpl -e 5pops_mertensiana.est -n 100000 -0 -m -M -L 48 -B 48 -c 48 --multiSFS -q

https://github.com/aglaszuk/Polygenic_Adaptation_Heliosperma/tree/main/06_demography/input_files

-M do maximum likelihood estimations
-d estimate the expected derived SFS
-n number of simulations
-L number of optimizations iterations to perform
-c number of threads
-q print as less info as possible
-0 Does not take into account monomorphic sites in observed SFS for parameter inference. Parameters are estimated only from the SFS and mutation rate is ignored. This option requires that one parameter be fixed in the .est file


###############################################
# run in loop to submit 50 iterations of the model with migration
# https://docs.alliancecan.ca/wiki/Cedar
# https://docs.alliancecan.ca/wiki/Running_jobs#Use_sbatch_to_submit_jobs

cd ~/scratch/Cassiope/Fastsimcoal2_Mar2023/scripts/

# tetragona 5 pop

for i in `seq 1 50`; do
cat << EOF > tet5_${i}.sh
#!/bin/bash
#SBATCH --account=rpp-rieseber
#SBATCH --time=0-11:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=187G

module load StdEnv/2020
module load fastsimcoal2/2.7.0.9

cd ~/scratch/Cassiope/Fastsimcoal2_Mar2023/runs/5pops_tet_$i
srun fsc27 -t *.tpl -n 100000 -e *.est -0 -m -M -L 48 -B 48 -c 48 --multiSFS -q

EOF

#sbatch tet5_${i}.sh

done


#-----------------------------
# loop for mertensiana model
cd ~/scratch/Cassiope/Fastsimcoal2_Mar2023/scripts/

for i in `seq 1 50`; do
cat << EOF > mer5_${i}.sh
#!/bin/bash
#SBATCH --account=rpp-rieseber
#SBATCH --time=0-02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=187G

module load StdEnv/2020
module load fastsimcoal2/2.7.0.9

cd ~/scratch/Cassiope/Fastsimcoal2_Mar2023/runs/5pops_mertensiana_$i
srun fsc27 -t *.tpl -n 100000 -e *.est -0 -m -M -L 48 -B 48 -c 48 --multiSFS -q

EOF

#sbatch mer5_${i}.sh

done


#-----------------------------
# loop for mertensiana model
cd ~/scratch/Cassiope/Fastsimcoal2_Mar2023/scripts/

for i in `seq 1 50`; do
cat << EOF > mer6_${i}.sh
#!/bin/bash
#SBATCH --account=rpp-rieseber
#SBATCH --time=0-3:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=187G

module load StdEnv/2020
module load fastsimcoal2/2.7.0.9

cd ~/scratch/Cassiope/Fastsimcoal2_Mar2023/runs/6pops_mertensiana_$i
srun fsc27 -t *.tpl -n 100000 -e *.est -0 -m -M -L 48 -B 48 -c 48 --multiSFS -q

EOF

#sbatch mer6_${i}.sh

done

#-------------------------------
# tetragona 6 pop
cd ~/scratch/Cassiope/Fastsimcoal2_Mar2023/scripts/

for i in `seq 1 50`; do
cat << EOF > tet6_${i}.sh
#!/bin/bash
#SBATCH --account=rpp-rieseber
#SBATCH --time=0-11:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=187G

module load StdEnv/2020
module load fastsimcoal2/2.7.0.9

cd ~/scratch/Cassiope/Fastsimcoal2_Mar2023/runs/6pops_tet_$i
srun fsc27 -t *.tpl -n 100000 -e *.est -0 -m -M -L 48 -B 48 -c 48 --multiSFS -q

EOF

#sbatch tet6_${i}.sh

done

########################
cd ~/scratch/Cassiope/Fastsimcoal2_Mar2023/scripts/

# tetragona Europe 5 pop

for i in `seq 1 50`; do
cat << EOF > tet_eur_5_${i}.sh
#!/bin/bash
#SBATCH --account=rpp-rieseber
#SBATCH --time=0-11:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=187G

module load StdEnv/2020
module load fastsimcoal2/2.7.0.9

cd ~/scratch/Cassiope/Fastsimcoal2_Mar2023/runs/5pops_tet_eur_$i
srun fsc27 -t *.tpl -n 100000 -e *.est -0 -m -M -L 48 -B 48 -c 48 --multiSFS -q

EOF

#sbatch tet_eur_5_${i}.sh

done


#-----------------------------
# loop for mertensiana Europe model
cd ~/scratch/Cassiope/Fastsimcoal2_Mar2023/scripts/

for i in `seq 1 50`; do
cat << EOF > mer_eur_5_${i}.sh
#!/bin/bash
#SBATCH --account=rpp-rieseber
#SBATCH --time=0-02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=187G

module load StdEnv/2020
module load fastsimcoal2/2.7.0.9

cd ~/scratch/Cassiope/Fastsimcoal2_Mar2023/runs/5pops_mer_eur_$i
srun fsc27 -t *.tpl -n 100000 -e *.est -0 -m -M -L 48 -B 48 -c 48 --multiSFS -q

EOF

#sbatch mer_eur_5_${i}.sh

done


#-----------------------------
# loop for 6 mertensiana model Europe
cd ~/scratch/Cassiope/Fastsimcoal2_Mar2023/scripts/

for i in `seq 1 50`; do
cat << EOF > mer6_eur_${i}.sh
#!/bin/bash
#SBATCH --account=rpp-rieseber
#SBATCH --time=0-3:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=187G

module load StdEnv/2020
module load fastsimcoal2/2.7.0.9

cd ~/scratch/Cassiope/Fastsimcoal2_Mar2023/runs/6pops_mer_eur_$i
srun fsc27 -t *.tpl -n 100000 -e *.est -0 -m -M -L 48 -B 48 -c 48 --multiSFS -q

EOF

sbatch mer6_eur_${i}.sh

done

#-------------------------------
# tetragona 6 pop Europe
cd ~/scratch/Cassiope/Fastsimcoal2_Mar2023/scripts/

for i in `seq 1 50`; do
cat << EOF > tet6_eur_${i}.sh
#!/bin/bash
#SBATCH --account=rpp-rieseber
#SBATCH --time=0-11:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=187G

module load StdEnv/2020
module load fastsimcoal2/2.7.0.9

cd ~/scratch/Cassiope/Fastsimcoal2_Mar2023/runs/6pops_tet_eur_$i
srun fsc27 -t *.tpl -n 100000 -e *.est -0 -m -M -L 48 -B 48 -c 48 --multiSFS -q

EOF

#sbatch tet6_eur_${i}.sh

done


#######################
# check output

# # Alaska split to Russia and then to Europe - mer model

cat ~/scratch/Cassiope/Fastsimcoal2_Mar2023/runs/5pops_mertensiana_*/5pops_mertensiana/5pops_mertensiana.bestlhoods > \
~/scratch/Cassiope/Fastsimcoal2_Mar2023/Total_mertensiana.bestlhoods

            N_Europe	N_Alaska	N_Russia	N_mertensiana	N_saximontana	N_Anc	        TDIVEurAla	TDIVRusAla	     TDIVSaxTet	         MaxEstLhood	MaxObsLhood
   33.00	234 283.00	930 650.00	527 321.00	    841 133.00	   195 377.00	5 820 155.00	30 556.00	136 859.00	    671 810.00	            -9679.67	-8446.53


#----------------
## Greenland from Europe mer model

cat ~/scratch/Cassiope/Fastsimcoal2_Mar2023/runs/6pops_mertensiana_*/6pops_mertensiana/6pops_mertensiana.bestlhoods > \
~/scratch/Cassiope/Fastsimcoal2_Mar2023/Total_mertensiana_Greenland.bestlhoods

	    N_Europe	N_Alaska	N_Russia	N_mertensiana	N_Greenland	 N_Anc	        TDIVGreEur	TDIVEurAla	TDIVRusAla	MaxEstLhood	MaxObsLhood
19.00	427 060.00	994 418.00	473 763.00	883 091.00	     25 298.00	5 760 086.00	12 216.00	62 031.00	155 779.00	-8406.85	-7410.13


#----------------
# Alaska split to Russia and then to Europe - tetragona model

cat ~/scratch/Cassiope/Fastsimcoal2_Mar2023/runs/5pops_tet_*/5pops_tet/5pops_tet.bestlhoods > \
~/scratch/Cassiope/Fastsimcoal2_Mar2023/Total_tetragona.bestlhoods

	    N_Europe	N_Alaska	    N_Russia	N_NWT	   N_Nunavut	TDIVEurAla	TDIVRusAla	migNWTNun	migNunNWT	MaxEstLhood	MaxObsLhood
21.00	257 770.00	2 021 998.00	554 466.00	476 666.00	2 280.00	    40 646.00	162 624.00	0.01	0.03	-7113.58	-6864.55

#----------------
# Greenland from Europe tetragona model


cat ~/scratch/Cassiope/Fastsimcoal2_Mar2023/runs/6pops_tet_*/6pops_tet/6pops_tet.bestlhoods > \
~/scratch/Cassiope/Fastsimcoal2_Mar2023/Total_tetragona_Greenland.bestlhoods


	    N_Europe	N_Alaska	N_Greenland	     N_NWT	   N_Nunavut	TDIVEurAla	   TDIVGreEur	migNWTNun	migNunNWT	MaxEstLhood	MaxObsLhood
23.00	910 164.00	3 108 542.00	64 392.00	160 462.00	4 466.00	 129 256.00	   10 808.00	 0.02	    0.03	    -5665.70	-5522.22


#----------------
# Russia split to Europe to Alaska - mer model

cat ~/scratch/Cassiope/Fastsimcoal2_Mar2023/runs/5pops_mer_eur_*/5pops_mer_eur/5pops_mer_eur.bestlhoods > \
~/scratch/Cassiope/Fastsimcoal2_Mar2023/Total_mertensiana_Europe.bestlhoods

N_Europe	N_Alaska	N_Russia	N_mertensiana	N_saximontana	N_Anc	TDIVEurAla	TDIVRusEur	TDIVSaxTet	MaxEstLhood	MaxObsLhood
489347.00	468948.00	887918.00	843483.00	    273946.00	5401852.00	43679.00	118354.00	810614.00	-9785.70	-8446.53


#----------------
# Russia split to Europe to Alaska - tetragona model

cat ~/scratch/Cassiope/Fastsimcoal2_Mar2023/runs/5pops_tet_eur_*/5pops_tet_eur/5pops_tet_eur.bestlhoods > \
~/scratch/Cassiope/Fastsimcoal2_Mar2023/Total_tetragona_Europe.bestlhoods


N_Europe        N_Alaska        N_Russia        N_NWT   N_Nunavut       TDIVEurAla      TDIVRusEur      migNWTNun       migNunNWT       MaxEstLhood    MaxObsLhood
1208897         2589205          3119330        745404  7934            209176          325828           0.4367214       0.4664329       -7292.979       -6865.001


#------------------
# Greenland from Alaska mer model

cat ~/scratch/Cassiope/Fastsimcoal2_Mar2023/runs/6pops_mer_eur_*/6pops_mer_eur/6pops_mer_eur.bestlhoods > \
~/scratch/Cassiope/Fastsimcoal2_Mar2023/Total_mertensiana_Europe_Greenland.bestlhoods

N_Europe	N_Alaska	N_Russia	N_mertensiana	N_Greenland	N_Anc	TDIVGreAla	TDIVEurAla	TDIVRusAla	MaxEstLhood	MaxObsLhood
320672.00	906408.00	384508.00	852034.00	60836.00	6133653.00	66333.00	46623.00	157072.00	-8405.72	-7408.16

#--------------------
# greenland from Alaska tetragona model

cat ~/scratch/Cassiope/Fastsimcoal2_Mar2023/runs/6pops_tet_eur_*/6pops_tet_eur/6pops_tet_eur.bestlhoods > \
~/scratch/Cassiope/Fastsimcoal2_Mar2023/Total_tetragona_Europe_Greenland.bestlhoods

N_Europe	N_Alaska	N_Greenland	N_NWT	    N_Nunavut	TDIVEurAla	TDIVGreAla	migNWTNun	migNunNWT	MaxEstLhood	MaxObsLhood
853813.00	3202335.00	61807.00	15073.00	10098.00	115113.00	9769.00	     0.00	     0.00	    -5590.28	-5522.22


#########################
# Bootstrap estimates on top model parameters

cd ~/scratch/Cassiope/Fastsimcoal2_Mar2023/
cp -vr ./runs/5pops_tet_1/* ./bootstrap_tet5/

# change number of independent loci [chromosomes] to 200000
# and datatype to DNA 

cd ~/scratch/Cassiope/Fastsimcoal2_Mar2023/bootstrap_tet5/5pops_tet/
cp 5pops_tet_maxL.par 5pops_tet_maxL_1boot.par
export var=$(wc -l < 5pops_tet_maxL_1boot.par)
var1=`expr $var - 1`
var2=`expr $var - 5`
sed -i "${var1} s/FREQ/DNA/" 5pops_tet_maxL_1boot.par 
sed -i "${var2} s/1/200000/" 5pops_tet_maxL_1boot.par

more 5pops_tet_maxL_1boot.par

fsc27 -i 5pops_tet_maxL.par -n10 -q -j -s0 -d --multisfs


Observed data file ("5pops_tet_maxL_jointDAFpop1_0.obs") does not seem to be present
Unable to compute the likelihood
Unable to read observed SFS, quitting program

#------------------------
# generate 100 SFS in separate subdirectories 
fsc27 -i 5pops_tet_maxL_1boot.par -n2 -j -d -s0 -x –I -q

Random generator seed : 996918

No population growth detected in input file

fastsimcoal2 is building 400000 genealogies ...

Simulating 200000 independent chromosomes using 12 batches and 1 threads.

Iteration 1/1 done in 235.922 secs


Program total execution time: 236.012 seconds

#####################
# making population files

cat << EOF > list_pops_5pop_mer.txt
PopSAM_1a	Russia
PopSAM_1	Russia
PopSAM_2a	Russia
PopSAM_2	Russia
PopSAM_3	Russia
PopSAM_4	Russia
PopSAM_5	Russia
PopYED_1a	Russia
PopYED_1	Russia
PopYED_2a	Russia
PopYED_2	Russia
PopYED_3a	Russia
PopYED_3	Russia
PopYED_4	Russia
PopYED_5a	Russia
PopYED_5	Russia
PopIMN_11	Alaska
PopIMN_13	Alaska
PopIMN_15	Alaska
PopIMN_17	Alaska
PopIMN_20	Alaska
PopIMN_21	Alaska
PopIMN_2	Alaska
PopIMN_5	Alaska
PopIMN_7	Alaska
PopMAT_13	Alaska
PopMAT_18	Alaska
PopMAT_2	Alaska
PopMAT_3	Alaska
PopMAT_5	Alaska
PopMAT_7	Alaska
PopMIL_25	Alaska
PopMIL_2	Alaska
PopMNT_12	Alaska
PopMNT_17	Alaska
PopMNT_18	Alaska
PopMNT_21	Alaska
PopMNT_22	Alaska
PopMNT_23	Alaska
PopMNT_4	Alaska
PopMNT_5	Alaska
PopMNT_7	Alaska
PopMNT_9	Alaska
PopSAG_16	Alaska
PopSAG_21	Alaska
PopSAG_23	Alaska
PopSAG_24	Alaska
PopSAG_5	Alaska
Pop025_24	Europe
PopDLG_10	Europe
PopDLG_11	Europe
PopDLG_12	Europe
PopDLG_13	Europe
PopDLG_14	Europe
PopDLG_16	Europe
PopDLG_17	Europe
PopDLG_7	Europe
PopDLG_8	Europe
PopDLG_9	Europe
PopDQG_10	Europe
PopDQG_11	Europe
PopDQG_13	Europe
PopDQG_14	Europe
PopDQG_1	Europe
PopDQG_2	Europe
PopDQG_3	Europe
PopDQG_4	Europe
PopDQG_6	Europe
PopIG_10	Europe
PopIG_12	Europe
PopIG_14	Europe
PopIG_2	Europe
PopIG_3	Europe
PopIG_4	Europe
PopIG_5	Europe
PopIG_6	Europe
PopLAJ_13	Europe
PopLAJ_17	Europe
PopLAJ_19	Europe
PopLAJ_20	Europe
PopLAJ_21	Europe
PopLAJ_23	Europe
PopLAJ_25	Europe
PopLAJ_26	Europe
PopLAJ_5	Europe
PopLON_11	Europe
PopLON_15	Europe
PopLON_17	Europe
PopLON_23	Europe
PopLON_28	Europe
PopLON_36	Europe
PopLON_4	Europe
PopLON_9	Europe
PopSW_18	Europe
PopSW_22	Europe
PopSW_29	Europe
PopSW_2	Europe
PopSW_31	Europe
PopSW_34	Europe
PopSW_36	Europe
PopSW_3	Europe
PopSW_45	Europe
PopZAC_1	Europe
PopZAC_20	Europe
PopZAC_21	Europe
PopZAC_23	Europe
PopZAC_24	Europe
PopZAC_25	Europe
PopZAC_26	Europe
PopZAC_28	Europe
PopZAC_7	Europe
PopHAR_13	BCmer
PopHAR_1	BCmer
PopHAR_21	BCmer
PopHAR_8	BCmer
PopPHE_1	BCmer
PopPHE_2	BCmer
PopPHE_3	BCmer
PopGEN_10	BCtet
PopGEN_1	BCtet
PopGEN_2	BCtet
PopGEN_3	BCtet
PopGEN_4	BCtet
PopGEN_5	BCtet
PopGEN_6	BCtet
PopGEN_7	BCtet
PopGEN_8	BCtet
PopGEN_9	BCtet

EOF

cat << EOF > list_pops_5pop_tet.txt
PopSAM_1a	Russia
PopSAM_1	Russia
PopSAM_2a	Russia
PopSAM_2	Russia
PopSAM_3	Russia
PopSAM_4	Russia
PopSAM_5	Russia
PopYED_1a	Russia
PopYED_1	Russia
PopYED_2a	Russia
PopYED_2	Russia
PopYED_3a	Russia
PopYED_3	Russia
PopYED_4	Russia
PopYED_5a	Russia
PopYED_5	Russia
PopIMN_11	Alaska
PopIMN_13	Alaska
PopIMN_15	Alaska
PopIMN_17	Alaska
PopIMN_20	Alaska
PopIMN_21	Alaska
PopIMN_2	Alaska
PopIMN_5	Alaska
PopIMN_7	Alaska
PopMAT_13	Alaska
PopMAT_18	Alaska
PopMAT_2	Alaska
PopMAT_3	Alaska
PopMAT_5	Alaska
PopMAT_7	Alaska
PopMIL_25	Alaska
PopMIL_2	Alaska
PopMNT_12	Alaska
PopMNT_17	Alaska
PopMNT_18	Alaska
PopMNT_21	Alaska
PopMNT_22	Alaska
PopMNT_23	Alaska
PopMNT_4	Alaska
PopMNT_5	Alaska
PopMNT_7	Alaska
PopMNT_9	Alaska
PopSAG_16	Alaska
PopSAG_21	Alaska
PopSAG_23	Alaska
PopSAG_24	Alaska
PopSAG_5	Alaska
Pop025_24	Europe
PopDLG_10	Europe
PopDLG_11	Europe
PopDLG_12	Europe
PopDLG_13	Europe
PopDLG_14	Europe
PopDLG_16	Europe
PopDLG_17	Europe
PopDLG_7	Europe
PopDLG_8	Europe
PopDLG_9	Europe
PopDQG_10	Europe
PopDQG_11	Europe
PopDQG_13	Europe
PopDQG_14	Europe
PopDQG_1	Europe
PopDQG_2	Europe
PopDQG_3	Europe
PopDQG_4	Europe
PopDQG_6	Europe
PopIG_10	Europe
PopIG_12	Europe
PopIG_14	Europe
PopIG_2	Europe
PopIG_3	Europe
PopIG_4	Europe
PopIG_5	Europe
PopIG_6	Europe
PopLAJ_13	Europe
PopLAJ_17	Europe
PopLAJ_19	Europe
PopLAJ_20	Europe
PopLAJ_21	Europe
PopLAJ_23	Europe
PopLAJ_25	Europe
PopLAJ_26	Europe
PopLAJ_5	Europe
PopLON_11	Europe
PopLON_15	Europe
PopLON_17	Europe
PopLON_23	Europe
PopLON_28	Europe
PopLON_36	Europe
PopLON_4	Europe
PopLON_9	Europe
PopSW_18	Europe
PopSW_22	Europe
PopSW_29	Europe
PopSW_2	Europe
PopSW_31	Europe
PopSW_34	Europe
PopSW_36	Europe
PopSW_3	Europe
PopSW_45	Europe
PopZAC_1	Europe
PopZAC_20	Europe
PopZAC_21	Europe
PopZAC_23	Europe
PopZAC_24	Europe
PopZAC_25	Europe
PopZAC_26	Europe
PopZAC_28	Europe
PopZAC_7	Europe
PopKUQ_14	NWT
PopKUQ_20	NWT
PopKUQ_26	NWT
PopKUQ_37	NWT
PopKUQ_4	NWT
PopKUQ_9	NWT
PopPEA_11	NWT
PopPEA_18	NWT
PopPEA_25	NWT
PopPEA_26	NWT
PopPEA_27	NWT
PopPEA_28	NWT
PopPEA_30	NWT
PopPEA_34	NWT
PopPEA_36	NWT
PopYAM_10	NWT
PopYAM_1	NWT
PopYAM_3	NWT
PopYAM_4	NWT
PopYAM_5	NWT
PopYAM_6	NWT
PopYAM_7	NWT
PopYAM_8	NWT
PopYAM_9	NWT
PopAlexNew_21	Nunavut
PopAlexNew_24	Nunavut
PopAlexNew_28	Nunavut
PopAlexNew_32	Nunavut
PopAlexNew_33	Nunavut
PopAlexNew_36	Nunavut
PopAlexNew_42	Nunavut
PopAlexNew_46	Nunavut
PopAlexNew_47	Nunavut
PopAlexNew_49	Nunavut
PopAlexNew_51	Nunavut
PopAlexNew_68	Nunavut
PopAlexNew_6	Nunavut
PopAlexNew_73	Nunavut
PopAXE_11	Nunavut
PopAXE_13	Nunavut
PopAXE_14	Nunavut
PopAXE_21	Nunavut
PopAXE_26	Nunavut
PopAXE_27	Nunavut
PopAXE_2	Nunavut
PopAXE_34	Nunavut
PopAXE_6	Nunavut
PopBY1_14	Nunavut
PopBY1_1	Nunavut
PopBY1_5	Nunavut
PopBY1_7	Nunavut
PopBY2_12	Nunavut
PopBY2_3	Nunavut
PopBY2_9	Nunavut
PopBY3_14	Nunavut
PopBY3_18	Nunavut
PopBY3_19	Nunavut
PopCR_10	Nunavut
PopCR_12	Nunavut
PopCR_13	Nunavut
PopCR_17	Nunavut
PopCR_1	Nunavut
PopCR_2	Nunavut
PopCR_4	Nunavut
PopCR_5	Nunavut
PopCR_7	Nunavut
PopCR_9	Nunavut
PopEUR_1	Nunavut
PopEUR_22	Nunavut
PopEUR_34	Nunavut
PopEUR_5	Nunavut
PopFOS_2	Nunavut
PopFOS_3	Nunavut
PopFOS_7	Nunavut
PopFOS_8	Nunavut
PopGF_1	Nunavut
PopGF_4	Nunavut
PopHAZ_18	Nunavut
PopHAZ_20	Nunavut
PopHAZ_22	Nunavut
PopHAZ_24	Nunavut
PopHAZ_26	Nunavut
PopHAZ_29	Nunavut
PopHAZ_33	Nunavut
PopHAZ_34	Nunavut
PopHAZ_6	Nunavut
PopHAZ_9	Nunavut
PopIq_10	Nunavut
PopIq_13	Nunavut
PopIq_15	Nunavut
PopIq_16	Nunavut
PopIq_17	Nunavut
PopIq_20	Nunavut
PopIq_2	Nunavut
PopIq_3	Nunavut
PopIq_9	Nunavut
PopKik_11	Nunavut
PopKik_13	Nunavut
PopKik_16	Nunavut
PopKik_17	Nunavut
PopKik_18	Nunavut
PopKik_1	Nunavut
PopKik_20	Nunavut
PopKik_2	Nunavut
PopKik_3	Nunavut
EOF