############################
# Cassiope tetragona genome
# Data download 
# Feb 8 2025
#############################

# data found here:
# https://www.ncbi.nlm.nih.gov/bioproject/1166291

# PacBio
# https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR32205129&display=download

# HiC
# https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR31644788&display=metadata


#-----------------------
# Download with 
# https://github.com/ncbi/sra-tools/wiki
# https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/

# Narval3
tmux new-session -s Cassiope
tmux attach-session -t Cassiope

module spider sra-toolkit/3.0.9

module load StdEnv/2023  gcc/12.3 sra-toolkit/3.0.9


# on Narval - access to Beluga
cd /lustre04/scratch/celphin
mkdir Cassiope_genome; cd Cassiope_genome
mkdir raw_data

cd /lustre04/scratch/celphin/Cassiope_genome/raw_data

wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos7/sra-pub-run-38/SRR032/32205/SRR32205129/SRR32205129.1
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos7/sra-pub-zq-41/SRR031/31644/SRR31644788/SRR31644788.lite.1


fasterq-dump SRR31644788
fasterq-dump SRR32205129 

# downloaded
SRR32205129.fastq # HiFi long reads
SRR31644788_1.fastq  # HiC forward reads
SRR31644788_2.fastq  # HiC reverse reads

##################################