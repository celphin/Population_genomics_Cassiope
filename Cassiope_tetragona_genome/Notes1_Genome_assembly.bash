############################
# Cassiope tetragona genome
# Hifiasm phased scaffolding
# Feb 8 2025
#############################

# Quality check 

cd /lustre04/scratch/celphin/Cassiope_genome/raw_data
mkdir fastQC_reports

module load StdEnv/2020
module load fastqc/0.11.9

gzip SRR32205129.fastq # HiFi long reads
gzip SRR31644788_1.fastq  # HiC forward reads
gzip SRR31644788_2.fastq  # HiC reverse reads

fastqc SRR31644788_1.fastq.gz -o ./fastQC_reports
fastqc SRR31644788_2.fastq.gz -o ./fastQC_reports
fastqc SRR32205129.fastq.gz -o ./fastQC_reports


############################
# Hifiasm
# https://hifiasm.readthedocs.io/en/latest/hic-assembly.html 
# https://github.com/chhylp123/hifiasm

tmux new-session -s Cassiope
tmux attach-session -t Cassiope

cd ~/scratch/Cassiope_genome

# Install hifiasm (requiring g++ and zlib)
git clone https://github.com/chhylp123/hifiasm
cd hifiasm && make

# Hi-C phasing with paired-end short reads in two FASTQ files

# from Oxyra: [M::main] Real time: 20581.846 sec; CPU: 571748.638 sec; Peak RSS: 67.930 GB

salloc -c32 --time 7:00:00 --mem 249G --account def-henryg

# Narval
/lustre04/scratch/celphin/Cassiope_genome/hifiasm/hifiasm \
-o /lustre04/scratch/celphin/Cassiope_genome/Cassiope_genome_Feb2025.asm \
-t32 \
--h1 /lustre04/scratch/celphin/Cassiope_genome/raw_data/SRR31644788_1.fastq.gz \
--h2 /lustre04/scratch/celphin/Cassiope_genome/raw_data/SRR31644788_2.fastq.gz \
/lustre04/scratch/celphin/Cassiope_genome/raw_data/SRR32205129.fastq.gz

# Writing /lustre04/scratch/celphin/Cassiope_genome/Cassiope_genome_Feb2025.asm.hic.hap2.p_ctg.gfa to disk...
# Inconsistency threshold for low-quality regions in BED files: 70%
# [M::main] Version: 0.24.0-r703
# [M::main] CMD: /lustre04/scratch/celphin/Cassiope_genome/hifiasm/hifiasm -o /lustre04/scratch/celphin/Cassiope_genome/Cassiope_genome_Feb2025.asm -t32 --h1 /lustre04/scratch/celphin/Cassiope_genome/raw_data/SRR31644788_1.fastq.gz --h2 /lustre04/scratch/celphin/Cassiope_genome/raw_data/SRR31644788_2.fastq.gz /lustre04/scratch/celphin/Cassiope_genome/raw_data/SRR32205129.fastq.gz
# [M::main] Real time: 8515.393 sec; CPU: 201932.536 sec; Peak RSS: 140.020 GB


#---------------------------------
# get fasta from gfa
# https://hifiasm.readthedocs.io/en/latest/interpreting-output.html

cd /lustre04/scratch/celphin/Cassiope_genome/

# Hap 1
awk '/^S/{print ">"$2;print $3}' Cassiope_genome_Feb2025.asm.hic.hap1.p_ctg.gfa > Cassiope_tetragona_haplotigs.Hap1.fa

# Hap 2
awk '/^S/{print ">"$2;print $3}' Cassiope_genome_Feb2025.asm.hic.hap2.p_ctg.gfa > Cassiope_tetragona_haplotigs.Hap2.fa

#--------------------------
# check fasta N50 etc.

# Assembly stats
#https://github.com/sanger-pathogens/assembly-stats

module load StdEnv/2020 
module load gcc/9.3.0
module load r-bundle-bioconductor/3.14
cd ~

git clone https://github.com/sanger-pathogens/assembly-stats.git

mkdir build
cd build
cmake -DINSTALL_DIR:PATH=/home/celphin/ ..
make
make test
make install

cd /home/celphin/assembly-stats/
./assembly-stats /lustre04/scratch/celphin/Cassiope_genome/Cassiope_tetragona_haplotigs.Hap1.fa

# stats for /lustre04/scratch/celphin/Cassiope_genome/Cassiope_tetragona_haplotigs.Hap1.fa
# sum = 1 232 948 422, n = 757, ave = 1628729.75, largest = 42 585 782
# N50 = 17 322 531, n = 23
# N60 = 12 896 457, n = 31
# N70 = 8 829 443, n = 42
# N80 = 5 347 285, n = 60
# N90 = 2 894 383, n = 93
# N100 = 21 549, n = 757
# N_count = 0
# Gaps = 0

./assembly-stats /lustre04/scratch/celphin/Cassiope_genome/Cassiope_tetragona_haplotigs.Hap2.fa

# stats for /lustre04/scratch/celphin/Cassiope_genome/Cassiope_tetragona_haplotigs.Hap2.fa
# sum = 164 196 174, n = 571, ave = 287558.97, largest = 7419039
# N50 = 663227, n = 53
# N60 = 475573, n = 83
# N70 = 319402, n = 125
# N80 = 200126, n = 191
# N90 = 109105, n = 300
# N100 = 22142, n = 571
# N_count = 0
# Gaps = 0

# 2n = 26 (2x)
# 12 Gbp

#################################################################
# Hi-C mapping and scaffolding - run on Cedar
# https://github.com/rieseberglab/haplotype_aware_scaffolding
########################################################

cd /lustre04/scratch/celphin/Cassiope_genome/

##############################
# Kaede and Eric's code for Juicer
# Using https://github.com/ericgonzalezs/ASSEMBLIES/blob/main/Juicer 

## Juicer and 3D-DNA pipelne for genome scaffolding ## 
#1. Run juicer to produce .hic and .assembly
#2. Open the .hic file in Juicebox (can use cloud version for bigger files)
#3. 3d-dna first step keeps all contigs, 2nd and 3rd step polishes - removes contigs to trash (on debris. 
#4. For placing cnntigs in the right place, look for lines of red/white. If there is a thick white line, move contigs. 

#-----------------
#create directories in the working directory
mkdir juicer
cd juicer
git clone https://github.com/rmdickson/juicer.git
cd ..

mkdir juicier
cd juicier
ln -s ../juicer/juicer/CPU/ scripts
cd scripts/common
wget https://hicfiles.tc4ga.com/public/juicer/juicer_tools.1.9.9_jcuda.0.8.jar
ln -s juicer_tools.1.9.9_jcuda.0.8.jar  juicer_tools.jar

#-----------------
cd ../../
mkdir fastq
cd fastq
ln -s /lustre04/scratch/celphin/Cassiope_genome/raw_data/SRR31644788_1.fastq.gz
ln -s /lustre04/scratch/celphin/Cassiope_genome/raw_data/SRR31644788_2.fastq.gz

#-----------------
cd ..
mkdir references
cd references
ln -s /lustre04/scratch/celphin/Cassiope_genome/Cassiope_tetragona_haplotigs.Hap1.fa
module load StdEnv/2020 bwa/0.7.17
bwa index Cassiope_tetragona_haplotigs.Hap1.fa

#-----------------
cd ..
mkdir restriction_sites
cd restriction_sites

#download and modify generate_site_positions.py https://github.com/aidenlab/juicer/blob/main/misc/generate_site_positions.py
wget https://raw.githubusercontent.com/aidenlab/juicer/refs/heads/main/misc/generate_site_positions.py
chmod 755 generate_site_positions.py

#Add to this section the path to your fasta
filenames = {
    'hg19': '/seq/references/Homo_sapiens_assembly19.fasta',
    'Cassiope_tetragona_H1': '../references/Cassiope_tetragona_haplotigs.Hap1.fa', #here I am adding the path to my fasta
  }

#run generate_site_positions.py
python generate_site_positions.py DpnII Cassiope_tetragona_H1

#after we finish running generate_site_positions.py, we run Chromosome_sizes.sh with this:               
for i in $(ls *DpnII.txt)
do
name=$(echo $i | cut -d "." -f 1 )
awk 'BEGIN{OFS="\t"}{print $1, $NF}'  $i > "$name"".chrom.sizes"
done

#-------------------
# tree looks like

cd /lustre04/scratch/celphin/Cassiope_genome/juicier/

tree
.
├── fastq
│   ├── SRR31644788_1.fastq.gz -> /lustre04/scratch/celphin/Cassiope_genome/raw_data/SRR31644788_1.fastq.gz
│   └── SRR31644788_2.fastq.gz -> /lustre04/scratch/celphin/Cassiope_genome/raw_data/SRR31644788_2.fastq.gz
├── references
│   ├── Cassiope_tetragona_haplotigs.Hap1.fa -> /lustre04/scratch/celphin/Cassiope_genome/Cassiope_tetragona_haplotigs.Hap1.fa
│   ├── Cassiope_tetragona_haplotigs.Hap1.fa.amb
│   ├── Cassiope_tetragona_haplotigs.Hap1.fa.ann
│   ├── Cassiope_tetragona_haplotigs.Hap1.fa.bwt
│   ├── Cassiope_tetragona_haplotigs.Hap1.fa.pac
│   └── Cassiope_tetragona_haplotigs.Hap1.fa.sa
├── restriction_sites
│   ├── Cassiope_tetragona_H1_DpnII.chrom.sizes
│   ├── Cassiope_tetragona_H1_DpnII.txt
│   └── generate_site_positions.py
└── scripts -> ../juicer/juicer/CPU/

chmod -R 770 * # give execute permissions

#-----------------
#run juicier without GPUS and more cores
cd /lustre04/scratch/celphin/Cassiope_genome/juicier/

cat << EOF > Juicer_sbatch.sh
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=def-cronk
#SBATCH --time=22:00:00
#SBATCH --ntasks-per-node=32
#SBATCH --mem=149G
module load StdEnv/2020 bwa/0.7.17 java/17.0.2 samtools/1.15.1  

cd /lustre04/scratch/celphin/Cassiope_genome/juicier/

bash scripts/juicer.sh -D $PWD -g Cassiope_tetragona_H1 -s DpnII -p restriction_sites/Cassiope_tetragona_H1_DpnII.chrom.sizes -y restriction_sites/Cassiope_tetragona_H1_DpnII.txt -z references/Cassiope_tetragona_haplotigs.Hap1.fa -t 32 -S early
EOF

sbatch Juicer_sbatch.sh

# [main] Real time: 6222.129 sec; CPU: 176466.670 sec
# (-:  Align of /lustre04/scratch/celphin/Cassiope_genome/juicier/splits/SRR31644788.fastq.gz.sam done successfully

#---------------
#run juicer with the file RUN_JUICER.sh  for compute canada
# with GPUS

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=def-henryg
#SBATCH --time=22:00:00
#SBATCH --gpus-per-node=a100:4
#SBATCH --ntasks-per-node=24
#SBATCH --mem=510000M
module load StdEnv/2020 bwa/0.7.17 java/17.0.2 samtools/1.15.1  cudacore/.11.7.0 ucx-cuda/1.12.1 gcc/11.3.0 openmpi/4.1.4 cuda/11.7

cd /lustre04/scratch/celphin/Cassiope_genome/juicier/

bash scripts/juicer.sh -D $PWD -g Cassiope_tetragona_H1 -s DpnII \
-p restriction_sites/Cassiope_tetragona_H1_DpnII.chrom.sizes \
-y restriction_sites/Cassiope_tetragona_H1_DpnII.txt \
-z references/Cassiope_tetragona_haplotigs.Hap1.fa  -t 24


########################################
# 3D DNA
# https://github.com/ericgonzalezs/ASSEMBLIES/blob/main/3D-DNA.sh

cd /lustre04/scratch/celphin/Cassiope_genome/
git clone https://github.com/aidenlab/3d-dna.git
cd 3d-dna
chmod -R 770 * # give execute permissions

#create python env
module load StdEnv/2020 python/3.11.2

virtualenv 3ddna

source 3ddna/bin/activate

pip install scipy numpy matplotlib #libraries required for 3d-dna 

deactivate

#----------------------
cd /lustre04/scratch/celphin/Cassiope_genome/3d-dna

# copy over files as symbolic links
ln -s /lustre04/scratch/celphin/Cassiope_genome/juicier/references/Cassiope_tetragona_haplotigs.Hap1.fa
ln -s /lustre04/scratch/celphin/Cassiope_genome/juicier/aligned/merged_nodups.txt 

#----------------
# run 3DDNA

cat << EOF >  3DDNA_Cass_Hap1.sh
#!/bin/bash

#SBATCH --account=def-henryg
#SBATCH --time=0-22:00:00
#SBATCH --cpus-per-task=15
#SBATCH --mem=100G

module load StdEnv/2020 python/3.10.2 java/17.0.2 lastz/1.04.03
export PATH="/lustre04/scratch/celphin/Cassiope_genome/3d-dna:$PATH"
source /lustre04/scratch/celphin/Cassiope_genome/3d-dna/3ddna/bin/activate

cd /lustre04/scratch/celphin/Cassiope_genome/3d-dna
/lustre04/scratch/celphin/Cassiope_genome/3d-dna/run-asm-pipeline.sh -r 2 Cassiope_tetragona_haplotigs.Hap1.fa merged_nodups.txt

deactivate

EOF

sbatch 3DDNA_Cass_Hap1.sh


#------------------------------
# taking a long time
# Finished writing norms
# :| Warning: No input for label1 was provided. Default for label1 is ":::fragment_"
# :| Warning: No input for label2 was provided. Default for label2 is ":::debris"

###########################
# Finalize output

cat << EOF >  3DDNA_Cass_Hap1_final.sh
#!/bin/bash

#SBATCH --account=def-henryg
#SBATCH --time=2-23:00:00
#SBATCH --cpus-per-task=3
#SBATCH --mem=100G

module load StdEnv/2020 python/3.10.2 java/17.0.2 lastz/1.04.03
export PATH="/lustre04/scratch/celphin/Cassiope_genome/3d-dna:$PATH"
source /lustre04/scratch/celphin/Cassiope_genome/3d-dna/3ddna/bin/activate

cd /lustre04/scratch/celphin/Cassiope_genome/3d-dna
run-asm-pipeline.sh -r 0 --stage seal Cassiope_tetragona_haplotigs.Hap1.fa merged_nodups.txt 

deactivate
EOF

sbatch 3DDNA_Cass_Hap1_final.sh



#######################################
# check assembly stats

module load StdEnv/2020 
module load gcc/9.3.0
module load r-bundle-bioconductor/3.14
cd ~/assembly-stats/build

./assembly-stats /lustre04/scratch/celphin/Cassiope_genome/3d-dna/Cassiope_tetragona_haplotigs.Hap1.FINAL.fasta


stats for /lustre04/scratch/celphin/Cassiope_genome/3d-dna/Cassiope_tetragona_haplotigs.Hap1.FINAL.fasta
sum = 1 233 198 922, n = 650, ave = 1897229.11, largest = 105 355 053
N50 = 94 033 818, n = 7
N60 = 89 300 599, n = 8
N70 = 88 331 210, n = 9
N80 = 81 564 107, n = 11
N90 = 76 597 457, n = 12

N100 = 1000, n = 650
N_count = 250500
Gaps = 501

# diploid
# 2n = 26
#(13 chromosomes in HiC map)
# each about 100 Mbp, whole genome is 1.2 Gbp
# https://nature.ca/aaflora/data/www/ercate.htm

###############################################
# BUSCO

tmux new-session -s BUSCO
tmux attach-session -t BUSCO
salloc -c10 --time 2:55:00 --mem 120000m --account def-henryg

module load StdEnv/2020
module load gcc/9.3.0 python augustus hmmer blast+ metaeuk prodigal r bbmap
source /lustre02/home/celphin/busco_env/bin/activate
deactivate

cd /lustre04/scratch/celphin/Cassiope_genome/
busco --offline --in  /lustre04/scratch/celphin/Cassiope_genome/3d-dna/Cassiope_tetragona_haplotigs.Hap1.FINAL.fasta \
--out  BUSCO_Cassiope_Assembly_eudicots  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path /lustre02/home/celphin/BUSCO_downloads/

        # --------------------------------------------------
        # |Results from dataset eudicots_odb10              |
        # --------------------------------------------------
        # |C:95.8%[S:90.5%,D:5.3%],F:0.7%,M:3.5%,n:2326     |
        # |2227   Complete BUSCOs (C)                       |
        # |2104   Complete and single-copy BUSCOs (S)       |
        # |123    Complete and duplicated BUSCOs (D)        |
        # |17     Fragmented BUSCOs (F)                     |
        # |82     Missing BUSCOs (M)                        |
        # |2326   Total BUSCO groups searched               |
        # --------------------------------------------------


#################################################
# To view in Juicer
# Cassiope_tetragona_haplotigs.Hap1.FINAL.fasta

# open in Juicebox 
# https://aidenlab.org/juicebox/

# Cassiope_tetragona_haplotigs.Hap1.FINAL.assembly
# Cassiope_tetragona_haplotigs.Hap1.final.hic

# use desktop version to edit manually
# https://www.youtube.com/watch?v=Nj7RhQZHM18
# need .assembly and .hic files

#################################
# after editing 
# convert back to fasta: juicebox_assembly_converter.py
# copy review.assembly into Cedar
# https://github.com/phasegenomics/juicebox_scripts

git clone https://github.com/phasegenomics/juicebox_scripts.git
chmod -R 770 /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/juicebox_scripts/*

tmux new-session -s Juicer 
tmux attach-session -t Juicer

salloc -c1 --time 2:50:00 --mem 120000m --account def-rieseber
module load StdEnv/2020 lastz/1.04.03 java/1.8.0_192 python/3.10.2 scipy-stack/2022a

cd /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/juicebox_scripts/

python /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/juicebox_scripts/juicebox_scripts/juicebox_assembly_converter.py \
-a /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/juicebox_scripts/data/Oxy_dig_1.scaffolds_FINAL.final.review.assembly  \
-f /home/celphin/projects/def-henryg/celphin/Oxyria/Oxy_dig_1.scaffolds_FINAL.fasta 


########################
# RepeatOBserver

cp /lustre04/scratch/celphin/Cassiope_genome/3d-dna/Cassiope_tetragona_haplotigs.Hap1.FINAL.fasta /lustre04/scratch/celphin/repeats/auto_script/

cd  /lustre04/scratch/celphin/repeats/auto_script/
mv Cassiope_tetragona_haplotigs.Hap1.FINAL.fasta Cassiopesax.fasta

cat << EOF > Auto_Cassiopesax.sh
#!/bin/bash
#SBATCH --account=def-henryg
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=19
#SBATCH --mem=191000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.2.2

srun Setup_Run_Repeats.sh -i Cassiopesax -f Cassiopesax.fasta -h H0 -c 19 -m 191000M -g FALSE

EOF

sbatch Auto_Cassiopesax.sh

# copy over spectra
cd /lustre04/scratch/celphin/Cassiope_genome/
mkdir RepeatOBserver
cp /lustre04/scratch/celphin/repeats/auto_script/output_chromosomes/Cassiopesax_H0-AT/Chr*part*/largeimages.png/All_spec1_Cassiopesax_H0-AT_Chr*part*_bp35_2000seq2501_*TRUE.png .



































#############################
# Old Oxyria example


#-------------------------------------
# Hifiasm
# https://hifiasm.readthedocs.io/en/latest/hic-assembly.html 
# https://github.com/chhylp123/hifiasm

tmux new-session -s Oxy
tmux attach-session -t Oxy

cd /home/celphin/projects/def-henryg/celphin/Oxyria/

# Install hifiasm (requiring g++ and zlib)
git clone https://github.com/chhylp123/hifiasm
cd hifiasm && make

# Hi-C phasing with paired-end short reads in two FASTQ files
hifiasm 
-o Oxyria1_Sept4.asm 
-t32 
--h1 NS.1871.003.IDT_i7_217---IDT_i5_217.Oxyria_Hi-C_R1.fastq.gz 
--h2 NS.1871.003.IDT_i7_217---IDT_i5_217.Oxyria_Hi-C_R2.fastq.gz
HiFi-reads.tar.gz


# try again with ccs reads
salloc -c32 --time 23:50:00 --mem 1510G --account def-henryg

/home/celphin/projects/def-henryg/celphin/Oxyria/hifiasm/hifiasm
 -o /home/celphin/projects/def-henryg/celphin/Oxyria/Oxyria1_Sept6.asm
 -t32 --h1 /home/celphin/projects/def-henryg/celphin/Oxyria/Hi-C/NS.1871.003.IDT_i7_217---IDT_i5_217.Oxyria_Hi-C_R1.fastq.gz
 --h2 /home/celphin/projects/def-henryg/celphin/Oxyria/Hi-C/NS.1871.003.IDT_i7_217---IDT_i5_217.Oxyria_Hi-C_R2.fastq.gz
 /home/celphin/projects/def-henryg/celphin/Oxyria/PacBio/PacBio_Aug2022/HiFi-ccs_reads.tar.gz

[M::ha_analyze_count] left: none
[M::ha_analyze_count] right: count[56] = 7429776
[M::ha_ft_gen] peak_hom: 56; peak_het: 28
[M::ha_ct_shrink::979.499*7.87] ==> counted 5093711 distinct minimizer k-mers
[M::ha_ft_gen::982.930*7.85@20.872GB] ==> filtered out 5093711 k-mers occurring 280 or more times
[M::ha_opt_update_cov] updated max_n_chain to 280
[M::yak_count] collected 878997595 minimizers
[M::ha_pt_gen::1514.717*7.90] ==> counted 38661657 distinct minimizer k-mers
[M::ha_pt_gen] count[4095] = 0 (for sanity check)
[M::ha_analyze_count] lowest: count[10] = 26177
[M::ha_analyze_count] highest: count[27] = 630362
[M::ha_hist_line]     1: ****************************************************************************************************> 18966805
[M::ha_hist_line]     2: ****************************************************************************************************> 1124580
[M::ha_hist_line]     3: ************************************************ 304805
[M::ha_hist_line]     4: ********************* 133135
[M::ha_hist_line]     5: ************ 74255
[M::ha_hist_line]     6: ******** 47920
[M::ha_hist_line]     7: ****** 35530
[M::ha_hist_line]     8: ***** 29171
[M::ha_hist_line]     9: **** 26359
[M::ha_hist_line]    10: **** 26177
[M::ha_hist_line]    11: ***** 28906
[M::ha_hist_line]    12: ***** 33250
[M::ha_hist_line]    13: ****** 39507
[M::ha_hist_line]    14: ******** 50732
[M::ha_hist_line]    15: *********** 67151
[M::ha_hist_line]    16: ************** 89142
[M::ha_hist_line]    17: ******************* 121318
[M::ha_hist_line]    18: ************************** 161233
[M::ha_hist_line]    19: ********************************* 209539                                                                                    [156/395]
[M::ha_hist_line]    20: ****************************************** 263922
[M::ha_hist_line]    21: **************************************************** 329234
[M::ha_hist_line]    22: *************************************************************** 395978
[M::ha_hist_line]    23: ************************************************************************** 466064
[M::ha_hist_line]    24: ************************************************************************************ 527763
[M::ha_hist_line]    25: ******************************************************************************************* 576581
[M::ha_hist_line]    26: ************************************************************************************************* 613578
[M::ha_hist_line]    27: **************************************************************************************************** 630362
[M::ha_hist_line]    28: **************************************************************************************************** 629509
[M::ha_hist_line]    29: ************************************************************************************************ 606825
[M::ha_hist_line]    30: ******************************************************************************************* 570854
[M::ha_hist_line]    31: *********************************************************************************** 522042
[M::ha_hist_line]    32: ************************************************************************** 465909
[M::ha_hist_line]    33: **************************************************************** 406468
[M::ha_hist_line]    34: ******************************************************* 346852
[M::ha_hist_line]    35: ********************************************** 291324
[M::ha_hist_line]    36: ************************************** 240437
[M::ha_hist_line]    37: ******************************** 199549
[M::ha_hist_line]    38: *************************** 168181
[M::ha_hist_line]    39: *********************** 144095
[M::ha_hist_line]    40: ********************* 129572
[M::ha_hist_line]    41: ******************** 123532
[M::ha_hist_line]    42: ******************** 125919
[M::ha_hist_line]    43: ********************* 135373
[M::ha_hist_line]    44: *********************** 147966
[M::ha_hist_line]    45: ************************** 165869
[M::ha_hist_line]    46: ****************************** 186248
[M::ha_hist_line]    47: ********************************* 209911
[M::ha_hist_line]    48: ************************************* 232500
[M::ha_hist_line]    49: ***************************************** 256126
[M::ha_hist_line]    50: ******************************************** 276333
[M::ha_hist_line]    51: *********************************************** 295641
[M::ha_hist_line]    52: ************************************************* 311695
[M::ha_hist_line]    53: *************************************************** 323590
[M::ha_hist_line]    54: **************************************************** 329718
[M::ha_hist_line]    55: ***************************************************** 332978
[M::ha_hist_line]    56: ***************************************************** 331651
[M::ha_hist_line]    57: **************************************************** 328722
[M::ha_hist_line]    58: ************************************************** 315036
[M::ha_hist_line]    59: *********************************************** 298332
[M::ha_hist_line]    60: ********************************************* 282808
[M::ha_hist_line]    61: ****************************************** 261853
[M::ha_hist_line]    62: ************************************** 240872
[M::ha_hist_line]    63: *********************************** 219614
[M::ha_hist_line]    64: ******************************* 198130
[M::ha_hist_line]    65: **************************** 175598
[M::ha_hist_line]    66: ************************ 153721
[M::ha_hist_line]    67: ********************* 132674
[M::ha_hist_line]    68: ****************** 114684
[M::ha_hist_line]    69: **************** 99156
[M::ha_hist_line]    70: ************* 84544
[M::ha_hist_line]    71: ************ 72509
[M::ha_hist_line]    72: ********** 62353
[M::ha_hist_line]    73: ********* 54281
[M::ha_hist_line]    74: ******** 47588
[M::ha_hist_line]    75: ******* 42134
[M::ha_hist_line]    76: ****** 37711
[M::ha_hist_line]    77: ***** 34569
[M::ha_hist_line]    78: ***** 32293
[M::ha_hist_line]    79: ***** 30532
[M::ha_hist_line]    80: ***** 29084
[M::ha_hist_line]    81: **** 28160
[M::ha_hist_line]    82: **** 27338
[M::ha_hist_line]    83: **** 26909
[M::ha_hist_line]    84: **** 26319
[M::ha_hist_line]    85: **** 25717
[M::ha_hist_line]    86: **** 25739
[M::ha_hist_line]    87: **** 24374
[M::ha_hist_line]    88: **** 24130
[M::ha_hist_line]    89: **** 23094
[M::ha_hist_line]    90: **** 22365
[M::ha_hist_line]    91: *** 21566
[M::ha_hist_line]    92: *** 20933
[M::ha_hist_line]    93: *** 20154
[M::ha_hist_line]    99: *** 17130
[M::ha_hist_line]   100: *** 16586
[M::ha_hist_line]   101: *** 16288
[M::ha_hist_line]   102: *** 15896
[M::ha_hist_line]   103: ** 15711
[M::ha_hist_line]   104: *** 15771
[M::ha_hist_line]   105: ** 15444
[M::ha_hist_line]   106: ** 15572
[M::ha_hist_line]   107: *** 15797
[M::ha_hist_line]   108: ** 15597
[M::ha_hist_line]   109: ** 15155
[M::ha_hist_line]   110: ** 15383
[M::ha_hist_line]   111: ** 15439
[M::ha_hist_line]   112: ** 15124
[M::ha_hist_line]   113: ** 14983
[M::ha_hist_line]   114: ** 14922
[M::ha_hist_line]   115: ** 14795
[M::ha_hist_line]   116: ** 14731
[M::ha_hist_line]   117: ** 14533
[M::ha_hist_line]   118: ** 13704
[M::ha_hist_line]   119: ** 13550
[M::ha_hist_line]   120: ** 13060
[M::ha_hist_line]   121: ** 12985
[M::ha_hist_line]   122: ** 12542
[M::ha_hist_line]   123: ** 12083
[M::ha_hist_line]   124: ** 11720
[M::ha_hist_line]   125: ** 11500
[M::ha_hist_line]   126: ** 10983
[M::ha_hist_line]   127: ** 10604
[M::ha_hist_line]   128: ** 10132
[M::ha_hist_line]   129: ** 9837
[M::ha_hist_line]   130: ** 9638
[M::ha_hist_line]   131: * 9147
[M::ha_hist_line]   132: * 9131
[M::ha_hist_line]   133: * 8903
[M::ha_hist_line]   134: * 8560
[M::ha_hist_line]   135: * 8440
[M::ha_hist_line]   136: * 8225
[M::ha_hist_line]   137: * 7989



[M::ha_hist_line]   176: * 4685
[M::ha_hist_line]   177: * 4749
[M::ha_hist_line]   178: * 4676
[M::ha_hist_line]   179: * 4629
[M::ha_hist_line]   180: * 4536
[M::ha_hist_line]   181: * 4613
[M::ha_hist_line]   182: * 4435
[M::ha_hist_line]   183: * 4313
[M::ha_hist_line]   184: * 4193
[M::ha_hist_line]   185: * 4206
[M::ha_hist_line]   186: * 4149
[M::ha_hist_line]   187: * 4105
[M::ha_hist_line]   188: * 4018
[M::ha_hist_line]   189: * 3976
[M::ha_hist_line]   190: * 3793
[M::ha_hist_line]   191: * 3825
[M::ha_hist_line]   192: * 3845
[M::ha_hist_line]   193: * 3631
[M::ha_hist_line]   194: * 3637
[M::ha_hist_line]   195: * 3643
[M::ha_hist_line]   196: * 3548
[M::ha_hist_line]   197: * 3386
[M::ha_hist_line]   198: * 3469
[M::ha_hist_line]   199: * 3409
[M::ha_hist_line]   200: * 3353
[M::ha_hist_line]   201: * 3301
[M::ha_hist_line]   202: * 3182
[M::ha_hist_line]   203: * 3177
[M::ha_hist_line]   204: * 3226
[M::ha_hist_line]   205: * 3241
[M::ha_hist_line]   206: * 3174
[M::ha_hist_line]  rest: ************************** 162545
[M::ha_analyze_count] left: none
[M::ha_analyze_count] right: count[55] = 332978
[M::ha_pt_gen] peak_hom: 55; peak_het: 27
[M::ha_ct_shrink::1514.835*7.90] ==> counted 19694852 distinct minimizer k-mers
[M::ha_pt_gen::] counting in normal mode
[M::yak_count] collected 878997595 minimizers
[M::ha_pt_gen::1711.608*9.32] ==> indexed 860030790 positions, counted 19694852 distinct minimizer k-mers



[M::ha_hist_line]  rest: **************************** 175767
[M::ha_analyze_count] left: none
[M::ha_analyze_count] right: count[56] = 331532
[M::ha_pt_gen] peak_hom: 56; peak_het: 28
[M::ha_ct_shrink::7463.273*26.53] ==> counted 18079381 distinct minimizer k-mers
[M::ha_pt_gen::] counting in normal mode
[M::yak_count] collected 871849600 minimizers
[M::ha_pt_gen::7655.420*26.38] ==> indexed 870601322 positions, counted 18079381 distinct minimizer k-mers



Writing /home/celphin/projects/def-henryg/celphin/Oxyria/Oxyria1_Sept6.asm.hic.p_ctg.gfa to disk...
[M::build_unitig_index::44.750] ==> Counting
[M::build_unitig_index::8.382] ==> Memory allocating
[M::build_unitig_index::63.636] ==> Filling pos
[M::build_unitig_index::0.277] ==> Sorting pos
[M::build_unitig_index::117.051] ==> HiC index has been built
[M::write_hc_pt_index] Index has been written.
[M::alignment_worker_pipeline::231.830] ==> Qualification
[M::dedup_hits::1.196] ==> Dedup
[M::mc_solve_core::0.129] ==> Partition
[M::dedup_hits::0.658] ==> Dedup
[M::stat] # misjoined unitigs: 7 (N50: 9858488); # corrected unitigs: 14 (N50: 8339145)
[M::adjust_weight_kv_u_trans_advance::0.697]
[M::mb_solve_core::1.460] ==> Partition
[M::mc_solve_core::1.966] ==> Partition
[M::adjust_weight_kv_u_trans_advance::4.504]
[M::mb_solve_core::1.432] ==> Partition
[M::mc_solve_core::1.929] ==> Partition
[M::adjust_weight_kv_u_trans_advance::4.492]
[M::mb_solve_core::1.421] ==> Partition
[M::mc_solve_core::1.902] ==> Partition
[M::stat] # heterozygous bases: 1064182092; # homozygous bases: 101045007
[M::reduce_hamming_error::0.695] # inserted edges: 3992, # fixed bubbles: 7
[M::adjust_utg_by_trio] primary contig coverage range: [47, infinity]
Writing /home/celphin/projects/def-henryg/celphin/Oxyria/Oxyria1_Sept6.asm.hic.hap1.p_ctg.gfa to disk...
[M::adjust_utg_by_trio] primary contig coverage range: [47, infinity]
Writing /home/celphin/projects/def-henryg/celphin/Oxyria/Oxyria1_Sept6.asm.hic.hap2.p_ctg.gfa to disk...
Inconsistency threshold for low-quality regions in BED files: 70%
[M::main] Version: 0.16.1-r375
[M::main] CMD: /home/celphin/projects/def-henryg/celphin/Oxyria/hifiasm/hifiasm -o /home/celphin/projects/def-henryg/celphin/Oxyria/Oxyria1_Sept6.asm -t32 --h1 /home/celphin/projects/def-henryg/celphin/Oxyria/Hi-C/NS.1871.003.IDT_i7_217---IDT_i5_217.Oxyria_Hi-C_R1.fastq.gz --h2 /home/celphin/projects/def-henryg/celphin/Oxyria/Hi-C/NS.1871.003.IDT_i7_217---IDT_i5_217.Oxyria_Hi-C_R2.fastq.gz /home/celphin/projects/def-henryg/celphin/Oxyria/PacBio/PacBio_Aug2022/HiFi-ccs_reads.tar.gz
[M::main] Real time: 20581.846 sec; CPU: 571748.638 sec; Peak RSS: 67.930 GB

#---------------------------------
# get fasta from gfa
# https://hifiasm.readthedocs.io/en/latest/interpreting-output.html

cd /home/celphin/projects/rpp-rieseber/celphin/Oxyria/Oxyria_haplotigs_Sept6/

awk '/^S/{print ">"$2;print $3}' Oxyria1_Sept6.asm.hic.p_ctg.gfa > Oxyria1_Sept6_haplotigs.p_ctg.fa

# Hap 1
awk '/^S/{print ">"$2;print $3}' Oxyria1_Sept6.asm.hic.hap1.p_ctg.gfa > Oxyria1_Sept6_haplotigs.Hap1.fa

# Hap 2
awk '/^S/{print ">"$2;print $3}' Oxyria1_Sept6.asm.hic.hap2.p_ctg.gfa > Oxyria1_Sept6_haplotigs.Hap2.fa

#--------------------------
# check fasta N50 etc.

# Assembly stats

module load StdEnv/2020 
module load gcc/9.3.0
module load r-bundle-bioconductor/3.14

cd /home/celphin/assembly-stats

./assembly-stats /home/celphin/projects/rpp-rieseber/celphin/Oxyria/Oxyria_haplotigs_Sept6/Oxyria1_Sept6_haplotigs.p_ctg.fa

stats for /home/celphin/projects/rpp-rieseber/celphin/Oxyria/Oxyria_haplotigs_Sept6/Oxyria1_Sept6_haplotigs.p_ctg.fa
sum = 590 373 804, n = 974, ave = 606 133.27, largest = 100 326 276
N50 = 56 351 117, n = 4
N60 = 39 295 833, n = 6
N70 = 34 734 330, n = 7
N80 = 28 174 138, n = 9
N90 = 10 901 797, n = 13
N100 = 14349, n = 974
N_count = 0
Gaps = 0

2n = 14 and 42

Genome size: 800-1000 Mbp
#----------------------------
cd /home/celphin/assembly-stats

./assembly-stats /home/celphin/projects/rpp-rieseber/celphin/Oxyria/Oxyria_haplotigs_Sept6/Oxyria1_Sept6_haplotigs.Hap1.fa

stats for /home/celphin/projects/rpp-rieseber/celphin/Oxyria/Oxyria_haplotigs_Sept6/Oxyria1_Sept6_haplotigs.Hap1.fa
sum = 582 616 305, n = 1101, ave = 529170.12, largest = 76049245
N50 = 38 327 256, n = 6
N60 = 26 598 077, n = 8
N70 = 17086624, n = 11
N80 = 10595219, n = 15
N90 = 2919708, n = 27
N100 = 12798, n = 1101
N_count = 0
Gaps = 0


./assembly-stats /home/celphin/projects/rpp-rieseber/celphin/Oxyria/Oxyria_haplotigs_Sept6/Oxyria1_Sept6_haplotigs.Hap2.fa

stats for /home/celphin/projects/rpp-rieseber/celphin/Oxyria/Oxyria_haplotigs_Sept6/Oxyria1_Sept6_haplotigs.Hap2.fa
sum = 566 594 759, n = 402, ave = 1409439.70, largest = 49230237
N50 = 27 807 003, n = 9
N60 = 26766067, n = 11
N70 = 22934860, n = 13
N80 = 15449834, n = 16
N90 = 5566760, n = 23
N100 = 14349, n = 402
N_count = 0
Gaps = 0


#################################################################
# Hi-C mapping and scaffolding - run on Cedar
# https://github.com/rieseberglab/haplotype_aware_scaffolding
########################################################


##############################
# try Kaede's code for Juicer
# https://github.com/ericgonzalezs/ASSEMBLIES/blob/main/Juicer 

## Juicer and 3D-DNA pipelne for genome scaffolding ## 
#1. Run juicer to produce .hic and .assembly
#2. Open the .hic file in Juicebox (can use cloud version for bigger files)
#3. 3d-dna first step keeps all contigs, 2nd and 3rd step polishes - removes contigs to trash (on debris. 
#4. For placing cnntigs in the right place, look for lines of red/white. If there is a thick white line, move contigs. 

#do this to run it as a pipeline based on Eric's script (https://github.com/ericgonzalezs/ASSEMBLIES/blob/main/Juicer)
#haplotype 1 
ln -s ../Juicer/juicer/CPU/ scripts
cd scripts/common
wget https://hicfiles.tc4ga.com/public/juicer/juicer_tools.1.9.9_jcuda.0.8.jar
ln -s juicer_tools.1.9.9_jcuda.0.8.jar  juicer_tools.jar
cd ../..
mkdir fastq
cd fastq
mv ../../Oxy_HiC_R1.fastq.gz Oxy_HiC_R1.fastq.gz
mv ../../Oxy_HiC_R2.fastq.gz Oxy_HiC_R2.fastq.gz
cd ../
mkdir references
cd references/
mv ../Oxy_draft_assembly_Hap1.fa   Oxy_draft_assembly_Hap1.fa
mv ../Oxy_draft_assembly_Hap2.fa   Oxy_draft_assembly_Hap2.fa

#-----------------------------
tmux new-session -s bwa
tmux attach-session -t bwa

cd /lustre04/scratch/celphin/Cassiope_genome/Juicer_run/references
salloc -c1 --time 2:55:00 --mem 120000m --account def-rieseber

module load bwa
bwa index Oxy_draft_assembly_Hap1.fa #Run this in a job
bwa index Oxy_draft_assembly_Hap2.fa #Run this in a job

#--------------------------------
cd ..
mkdir restriction_sites 
cd restriction_sites
wget https://raw.githubusercontent.com/aidenlab/juicer/main/misc/generate_site_positions.py

#---- here download the juicer/misc/generate_site_positions.py and edit accordingly
  filenames = {
    'hg19': '/seq/references/Homo_sapiens_assembly19.fasta',
    'mm9' : '/seq/references/Mus_musculus_assembly9.fasta',
    'mm10': '/seq/references/Mus_musculus_assembly10.fasta',
    'hg18': '/seq/references/Homo_sapiens_assembly18.fasta',
    'Oxy_Hap1': '../references/Oxy_draft_assembly_Hap1.fa', 
    'Oxy_Hap2': '../references/Oxy_draft_assembly_Hap2.fa',#here you put your contig/assembly and its path
  }

tmux new-session -s python
tmux attach-session -t python

/lustre04/scratch/celphin/Cassiope_genome/Juicer_run/restriction_sites
salloc -c1 --time 2:55:00 --mem 120000m --account def-rieseber

module load python
python generate_site_positions.py DpnII Oxy_Hap1 #Run this in a job
python generate_site_positions.py DpnII Oxy_Hap2 #Run this in a job

#--------------------------------
#generate a file Chromosome_sizes.sh with this:
for i in $(ls *_DpnII.txt)
do
name=$(echo $i | cut -d "." -f 1 )
awk 'BEGIN{OFS="\t"}{print $1, $NF}'  $i > "$name"".chrom.sizes"
done

#-----------------------
# tree looks like

tree
.
├── fastq
│   ├── Oxy_HiC_R1.fastq.gz
│   └── Oxy_HiC_R2.fastq.gz
├── juicer_run_Oxyria_Hap1.sh
├── juicer_run_Oxyria_Hap2.sh
├── references
│   ├── Oxy_draft_assembly_Hap1.fa
│   ├── Oxy_draft_assembly_Hap1.fa.amb
│   ├── Oxy_draft_assembly_Hap1.fa.ann
│   ├── Oxy_draft_assembly_Hap1.fa.bwt
│   ├── Oxy_draft_assembly_Hap1.fa.pac
│   ├── Oxy_draft_assembly_Hap1.fa.sa
│   ├── Oxy_draft_assembly_Hap2.fa
│   ├── Oxy_draft_assembly_Hap2.fa.amb
│   ├── Oxy_draft_assembly_Hap2.fa.ann
│   ├── Oxy_draft_assembly_Hap2.fa.bwt
│   ├── Oxy_draft_assembly_Hap2.fa.pac
│   └── Oxy_draft_assembly_Hap2.fa.sa
├── restriction_sites
│   ├── Oxy_Hap1_DpnII.chrom.sizes
│   ├── Oxy_Hap1_DpnII.txt
│   ├── Oxy_Hap2_DpnII.chrom.sizes
│   ├── Oxy_Hap2_DpnII.txt
│   └── generate_site_positions.py
└── scripts -> ../Juicer/juicer/CPU/

#----------------------------
# run Juicer 

#run juicier with out GPUS and more cores
nano juicer_run_Oxyria_Hap1.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=def-rieseber
#SBATCH --time=3-0
#SBATCH --ntasks-per-node=32
#SBATCH --mem=149G
module load StdEnv/2020 bwa/0.7.17 java/17.0.2 samtools/1.15.1
cd /lustre04/scratch/celphin/Cassiope_genome/Juicer_run
bash scripts/juicer.sh -D $PWD -g Oxy_Hap1 -s DpnII -p restriction_sites/Oxy_Hap1_DpnII.chrom.sizes -y restriction_sites/Oxy_Hap1_DpnII.txt -z references/Oxy_draft_assembly_Hap1.fa -t 32 -S early

sbatch juicer_run_Oxyria_Hap1.sh
(-:  Align of /lustre04/scratch/celphin/Cassiope_genome/Juicer_run/splits/Oxy_HiC.fastq.gz.sam done successfully

#------------------
nano juicer_run_Oxyria_Hap2.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=def-rieseber
#SBATCH --time=3-0
#SBATCH --ntasks-per-node=32
#SBATCH --mem=149G
module load StdEnv/2020 bwa/0.7.17 java/17.0.2 samtools/1.15.1
cd /lustre04/scratch/celphin/Cassiope_genome/Juicer_run
bash scripts/juicer.sh -D $PWD -g Oxy_Hap2 -s DpnII -p restriction_sites/Oxy_Hap2_DpnII.chrom.sizes -y restriction_sites/Oxy_Hap2_DpnII.txt -z references/Oxy_draft_assembly_Hap2.fa -t 32 -S early

sbatch juicer_run_Oxyria_Hap2.sh

# worried that the files overlap between haplotypes

##################################
# try running again in separate folders
cd /lustre04/scratch/celphin/Cassiope_genome/Juicer_run
cp -rv fastq ../Oxy_Hap1_juicer
cp -rv references/ ../Oxy_Hap1_juicer
cp -rv restriction_sites ../Oxy_Hap1_juicer
cp -rv scripts ../Oxy_Hap1_juicer
cd ..
cp -rv Oxy_Hap1_juicer Oxy_Hap2_juicer

#-------------------
cd /lustre04/scratch/celphin/Cassiope_genome/Oxy_Hap1_juicer
nano juicer_run_Oxyria_Hap1.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=def-rieseber
#SBATCH --time=0-22:00:00
#SBATCH --ntasks-per-node=32
#SBATCH --mem=149G
module load StdEnv/2020 bwa/0.7.17 java/17.0.2 samtools/1.15.1
cd /lustre04/scratch/celphin/Cassiope_genome/Oxy_Hap1_juicer
bash scripts/juicer.sh -D $PWD -g Oxy_Hap1 -s DpnII -p restriction_sites/Oxy_Hap1_DpnII.chrom.sizes -y restriction_sites/Oxy_Hap1_DpnII.txt -z references/Oxy_draft_assembly_Hap1.fa -t 32 -S early

sbatch juicer_run_Oxyria_Hap1.sh

#------------------
cd /lustre04/scratch/celphin/Cassiope_genome/Oxy_Hap2_juicer
nano juicer_run_Oxyria_Hap2.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=def-rieseber
#SBATCH --time=0-22:00:00
#SBATCH --ntasks-per-node=32
#SBATCH --mem=149G
module load StdEnv/2020 bwa/0.7.17 java/17.0.2 samtools/1.15.1
cd /lustre04/scratch/celphin/Cassiope_genome/Oxy_Hap2_juicer
bash scripts/juicer.sh -D $PWD -g Oxy_Hap2 -s DpnII -p restriction_sites/Oxy_Hap2_DpnII.chrom.sizes -y restriction_sites/Oxy_Hap2_DpnII.txt -z references/Oxy_draft_assembly_Hap2.fa -t 32 -S early

sbatch juicer_run_Oxyria_Hap2.sh




########################################
# try 3D DNA
# https://github.com/ericgonzalezs/ASSEMBLIES/blob/main/3D-DNA.sh

cd /lustre04/scratch/celphin/Cassiope_genome/
git clone https://github.com/aidenlab/3d-dna.git
cd 3d-dna
chmod -R 770 * # give execute permissions


#create python env
module load StdEnv/2020 python/3.11.2

virtualenv 3ddna

source 3ddna/bin/activate

pip install scipy numpy matplotlib #libraries required for 3d-dna 

deactivate

#----------------------
cd /lustre04/scratch/celphin/Cassiope_genome/3d-dna
# copy over files as symbolic links
ln -s /lustre04/scratch/celphin/Cassiope_genome/Oxy_Hap1_juicer/references/Oxy_draft_assembly_Hap1.fa
ln -s /lustre04/scratch/celphin/Cassiope_genome/Oxy_Hap1_juicer/references/Oxy_draft_assembly_Hap2.fa

ln -s /lustre04/scratch/celphin/Cassiope_genome/Oxy_Hap1_juicer/aligned/merged_nodups_Oxy_Hap1.txt 
ln -s /lustre04/scratch/celphin/Cassiope_genome/Oxy_Hap2_juicer/aligned/merged_nodups_Oxy_Hap2.txt

#----------------
# run 3DDNA

nano 3DDNA_Oxy_Hap1.sh

#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=1-22:00:00
#SBATCH --cpus-per-task=15
#SBATCH --mem=100G

module load StdEnv/2020 python/3.10.2 java/17.0.2 lastz/1.04.03
export PATH="/lustre04/scratch/celphin/Cassiope_genome/3d-dna:$PATH"
source /lustre04/scratch/celphin/Cassiope_genome/3d-dna/3ddna/bin/activate

cd /lustre04/scratch/celphin/Cassiope_genome/3d-dna/Oxy_hap1
/lustre04/scratch/celphin/Cassiope_genome/3d-dna/run-asm-pipeline.sh -r 2 Oxy_draft_assembly_Hap1.fa merged_nodups_Oxy_Hap1.txt

deactivate

sbatch 3DDNA_Oxy_Hap1.sh

#----------------------------
nano 3DDNA_Oxy_Hap2.sh

#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=1-22:00:00
#SBATCH --cpus-per-task=15
#SBATCH --mem=100G

module load StdEnv/2020 python/3.10.2 java/17.0.2 lastz/1.04.03
export PATH="/lustre04/scratch/celphin/Cassiope_genome/3d-dna:$PATH"
source /lustre04/scratch/celphin/Cassiope_genome/3d-dna/3ddna/bin/activate

cd /lustre04/scratch/celphin/Cassiope_genome/3d-dna/Oxy_hap2
/lustre04/scratch/celphin/Cassiope_genome/3d-dna/run-asm-pipeline.sh -r 2 Oxy_draft_assembly_Hap2.fa merged_nodups_Oxy_Hap2.txt

deactivate

sbatch 3DDNA_Oxy_Hap2.sh

#------------------------------
# taking a long time
Finished writing norms
:| Warning: No input for label1 was provided. Default for label1 is ":::fragment_"
:| Warning: No input for label2 was provided. Default for label2 is ":::debris"

###########################
# Finalize output

nano 3DDNA_Oxy_Hap1_final.sh

#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=0-22:00:00
#SBATCH --cpus-per-task=15
#SBATCH --mem=100G

module load StdEnv/2020 python/3.10.2 java/17.0.2 lastz/1.04.03
export PATH="/lustre04/scratch/celphin/Cassiope_genome/3d-dna:$PATH"
source /lustre04/scratch/celphin/Cassiope_genome/3d-dna/3ddna/bin/activate

run-asm-pipeline.sh -r 0 --stage finalize Oxy_draft_assembly_Hap1.fa merged_nodups_Oxy_Hap1.txt 

deactivate

sbatch 3DDNA_Oxy_Hap1_final.sh

#-----------------------------
nano 3DDNA_Oxy_Hap2_final.sh

#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=1-22:00:00
#SBATCH --cpus-per-task=15
#SBATCH --mem=100G

module load StdEnv/2020 python/3.10.2 java/17.0.2 lastz/1.04.03
export PATH="/lustre04/scratch/celphin/Cassiope_genome/3d-dna:$PATH"
source /lustre04/scratch/celphin/Cassiope_genome/3d-dna/3ddna/bin/activate

run-asm-pipeline.sh -r 0 --stage seal Oxy_draft_assembly_Hap2.fa merged_nodups_Oxy_Hap2.txt 

deactivate

sbatch 3DDNA_Oxy_Hap2_final.sh

#######################################
# check assembly stats

module load StdEnv/2020 
module load gcc/9.3.0
module load r-bundle-bioconductor/3.14
cd ~/assembly-stats/build

./assembly-stats /lustre04/scratch/celphin/Cassiope_genome/3d-dna/Oxy_hap1/Oxy_draft_assembly_Hap1.FINAL.fasta

stats for /lustre04/scratch/celphin/Cassiope_genome/3d-dna/Oxy_hap1/Oxy_draft_assembly_Hap1.FINAL.fasta
sum = 582 998 305, n = 1886, ave = 309118.93, largest = 79 385 870
N50 = 74 109 000, n = 4
N60 = 70 648 606, n = 5
N70 = 69 949 313, n = 6
N80 = 67 411 532, n = 7

N90 = 100188, n = 65
N100 = 1000, n = 1886
N_count = 382000
Gaps = 764

#----------------------------
./assembly-stats /lustre04/scratch/celphin/Cassiope_genome/3d-dna/Oxy_hap2/Oxy_draft_assembly_Hap2.FINAL.fasta

stats for /lustre04/scratch/celphin/Cassiope_genome/3d-dna/Oxy_hap2/Oxy_draft_assembly_Hap2.FINAL.fasta
sum = 566 728 759, n = 615, ave = 921510.18, largest = 86 137 797
N50 = 74 731 636, n = 4
N60 = 72 720 770, n = 5
N70 = 72 112 163, n = 6
N80 = 72 112 163, n = 6
N90 = 70 815 000, n = 7

N100 = 1000, n = 615
N_count = 134000
Gaps = 268

#################################################