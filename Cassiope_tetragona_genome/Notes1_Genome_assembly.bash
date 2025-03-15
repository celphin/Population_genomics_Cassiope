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
# https://github.com/celphin/RepeatOBserverV1

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


