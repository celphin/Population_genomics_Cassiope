#####################################
# Entire Cassiope analysis - part 1
# March 2023
# Demultiplex and dDocent
#####################################

#######################################################
# Demultiplexing
#######################################################
# see GBS barcodes folder on github for PE barcodes

cp HI*.gz /scratch/celphin/GBS_Cassiope/Dec2019_Demultiplex/
cp CASS_Barcodes_lane* /scratch/celphin/GBS_Cassiope/Dec2019_Demultiplex/
cp GBS2enzymedemultiplex_notrim.pl /scratch/celphin/GBS_Cassiope/Dec2019_Demultiplex/

cd /scratch/celphin/GBS_Cassiope/Dec2019_Demultiplex/
salloc -c2 --time 20:00:00 --mem 120000m --account def-rieseber
start 1:30pm on Cedar 833
finished 8pm

tmux new-session -s lane1
tmux new-session -s lane2

tmux attach-session -t lane1
tmux attach-session -t lane2

perl GBS2enzymedemultiplex_notrim.pl CASS_Barcodes_lane1.txt HI.4977.003.lane1_343A_R1.fastq.gz HI.4977.003.lane1_343A_R2.fastq.gz ./data/L1-  

perl GBS2enzymedemultiplex_notrim.pl CASS_Barcodes_lane2.txt HI.4962.002.lane2_566A_R1.fastq.gz HI.4962.002.lane2_566A_R2.fastq.gz ./data/L2-

# Demultiplexed data found BioProject: PRJNA824830
# https://www.ncbi.nlm.nih.gov/bioproject/?term=Cassiope%20tetragona%20subsp.%20tetragona%5Borgn%5D

###########################
# ddocent

# try run to build the assembly 

# create tmux session (see blog)
tmux new-session -s Cassiope
tmux attach-session -t Cassiope

# get interactive allocation
salloc -c48 --time 23:00:00 --mem 187G --account def-rieseber

#--------------------------
# set up Ddocent environment in session (see blog)
~/miniconda2/bin/conda create --name dDocent  
source ~/miniconda2/bin/activate dDocent

# update conda
conda update -n base -c defaults conda
# Your installed version is: 2.17

#source ~/miniconda2/bin/deactivate dDocent
conda install -c bioconda ddocent
# Update from dDocent 2.6.0 to dDocent 2.9.4

#----------------
# change directory to your fastq files named according to Ddocent:
# https://www.ddocent.com/UserGuide/#quality-filtering

cd /project/6019339/celphin/Cassiope/Feb2023_dDocent/reference_fastq

# activate conda environment
source ~/miniconda2/bin/activate dDocent

#use dDocent to build the assembly (leave k2 parameter blank)
dDocent

55 individuals are detected. Is this correct? Enter yes or no and press [ENTER]
yes

dDocent detects 48 processors available on this system.
Please enter the maximum number of processors to use for this analysis.
48

dDocent detects 187 gigabytes of maximum memory available on this system.
Please enter the maximum memory to use for this analysis in gigabytes
For example, to limit dDocent to ten gigabytes, enter 10
This option does not work with all distributions of Linux.  If runs are hanging at variant calling, enter 0
Then press [ENTER]
180

Do you want to quality trim your reads?
Type yes or no and press [ENTER]?
no

Do you want to perform an assembly?
Type yes or no and press [ENTER].
yes

What type of assembly would you like to perform?  Enter SE for single end, PE for paired-end, RPE for paired-end sequencing for RAD protocols with random shearing, or OL for paired-end sequencing that has substantial overlap.
Then press [ENTER]
PE

Reads will be assembled with Rainbow

CD-HIT will cluster reference sequences by similarity. The -c parameter (% similarity to cluster) may need to be changed for your taxa.
Would you like to enter a new c parameter now? Type yes or no and press [ENTER]
yes
Please enter new value for c. Enter in decimal form (For 90%, enter 0.9)
0.95

Do you want to map reads?  Type yes or no and press [ENTER]
no

Mapping will not be performed
Do you want to use FreeBayes to call SNPs?  Please type yes or no and press [ENTER]
no

Please enter your email address.  dDocent will email you when it is finished running.
Don't worry; dDocent has no financial need to sell your email address to spammers.
cassandra.elphinstone@shaw.ca

dDocent will require input during the assembly stage.  Please wait until prompt says it is safe to move program to the background.


                       Number of Unique Sequences with More than X Coverage (Counted within individuals)

    3e+06 ++----------+-----------+----------+-----------+-----------+-----------+----------+-----------+----------++
          +           +           +          +           +           +           +          +           +           +
          *                                                                                                         |
          |*                                                                                                        |
  2.5e+06 ++*                                                                                                      ++
          | *                                                                                                       |
          |  *                                                                                                      |
          |   *                                                                                                     |
    2e+06 ++   *                                                                                                   ++
          |    *                                                                                                    |
          |     *****                                                                                               |
  1.5e+06 ++         *                                                                                             ++
          |           *****                                                                                         |
          |                *                                                                                        |
          |                 ******                                                                                  |
    1e+06 ++                      *****                                                                            ++
          |                            ************                                                                 |
          |                                        ******************                                               |
          |                                                          ******************                             |
   500000 ++                                                                           *****************************+
          |                                                                                                         *
          |                                                                                                         |
          +           +           +          +           +           +           +          +           +           +
        0 ++----------+-----------+----------+-----------+-----------+-----------+----------+-----------+----------++
          2           4           6          8           10          12          14         16          18          20
                                                           Coverage

Please choose data cutoff.  In essence, you are picking a minimum (within individual) coverage level for a read (allele) to be used in the reference assembly
3



                                 Number of Unique Sequences present in more than X Individuals
 Number of Unique Sequences
   300000 ++----------------+----------------+-----------------+-----------------+----------------+----------------++
          +                 +                +                 +                 +                +                 +
          |      *                                                                                                  |
          |       *                                                                                                 |
   250000 ++      *                                                                                                ++
          |        *                                                                                                |
          |        *                                                                                                |
          |         *                                                                                               |
   200000 ++        *                                                                                              ++
          |          *                                                                                              |
          |           *                                                                                             |
   150000 ++           *                                                                                           ++
          |            *                                                                                            |
          |             ***                                                                                         |
          |                *                                                                                        |
   100000 ++                ***                                                                                    ++
          |                                                                                                         |
          |                    ****                                                                                 |
          |                        ***                                                                              |
    50000 ++                          ****                                                                         ++
          |                               ***                                                                       |
          |                                  ****                                                                   |
          +                 +                +   ******************              +                +                 +
        0 ++----------------+----------------+-----------------+---***************************************---------++
          0                 5                10                15                20               25                30
                                                     Number of Individuals

Please choose data cutoff.  Pick point right before the assymptote. A good starting cutoff might be 10% of the total number of individuals
5


At this point, all configuration information has been entered and dDocent may take several hours to run.
It is recommended that you move this script to a background operation and disable terminal input and output.
All data and logfiles will still be recorded.
To do this:
Press control and Z simultaneously
Type 'bg' without the quotes and press enter
Type 'disown -h' again without the quotes and press enter

Now sit back, relax, and wait for your analysis to finish

dDocent assembled 102787 sequences (after cutoffs) into 23890 contigs

dDocent has finished with an analysis in /project/6019339/celphin/Cassiope/Feb2023_dDocent/fastq/reference_fastq

dDocent started Mon Mar 6 16:58:48 PST 2023

dDocent finished Mon Mar 6 17:04:39 PST 2023

dDocent 2.9.4
The 'd' is silent, hillbilly.

#########################
# move reference to main directory above and copy to tetragona and outgroup directories

source ~/miniconda2/bin/deactivate dDocent

mv reference.fasta ..

cd /project/6019339/celphin/Cassiope/Feb2023_dDocent/

head reference.fasta
wc -l reference.fasta 
47780 reference.fasta

mawk '/>/' reference.fasta | wc -l # old reference had 38k contigs
23890

grep '^NAATT.*N*.*CCGN$' reference.fasta | wc -l
0

#####################################
# try mapping for all the tetragona individuals (including ancient)
# Note: did not include the ancient, outgroup and subspp individuals in the reference

# create tmux session (see blog)
tmux new-session -s Cassiope2
tmux attach-session -t Cassiope2

salloc -c48 --time 23:00:00 --mem 187G --account def-rieseber

cd /project/6019339/celphin/Cassiope/Feb2023_dDocent/tetragona

source ~/miniconda2/bin/activate dDocent 

dDocent

351 individuals are detected. Is this correct? Enter yes or no and press [ENTER]
yes
Proceeding with 351 individuals

dDocent detects 48 processors available on this system.
Please enter the maximum number of processors to use for this analysis.
48

dDocent detects 187 gigabytes of maximum memory available on this system.
Please enter the maximum memory to use for this analysis in gigabytes
For example, to limit dDocent to ten gigabytes, enter 10
This option does not work with all distributions of Linux.  If runs are hanging at variant calling, enter 0
Then press [ENTER]
180

Do you want to quality trim your reads?
Type yes or no and press [ENTER]?
no

Do you want to perform an assembly?
Type yes or no and press [ENTER].
no

Reference contigs need to be in a file named reference.fasta

Do you want to map reads?  Type yes or no and press [ENTER]
yes

BWA will be used to map reads.  You may need to adjust -A -B and -O parameters for your taxa.
Would you like to enter a new parameters now? Type yes or no and press [ENTER]
yes
Please enter new value for A (match score).  It should be an integer.  Default is 1.
1
Please enter new value for B (mismatch score).  It should be an integer.  Default is 4.
4
Please enter new value for O (gap penalty).  It should be an integer.  Default is 6.
6

Do you want to use FreeBayes to call SNPs?  Please type yes or no and press [ENTER]
no

Please enter your email address.  dDocent will email you when it is finished running.
Don't worry; dDocent has no financial need to sell your email address to spammers.
cassandra.elphinstone@shaw.ca

At this point, all configuration information has been entered and dDocent may take several hours to run.
It is recommended that you move this script to a background operation and disable terminal input and output.
All data and logfiles will still be recorded.
To do this:
Press control and Z simultaneously
Type bg and press enter
Type disown -h and press enter
Now sit back, relax, and wait for your analysis to finish

Using BWA to map reads

Creating alignment intervals



dDocent has finished with an analysis in /project/6019339/celphin/Cassiope/Feb2023_dDocent/tetragona

dDocent started Mon Mar 6 17:20:58 PST 2023

dDocent finished Mon Mar 6 22:11:18 PST 2023



dDocent 2.9.4
The 'd' is silent, hillbilly.
#/home/celphin/miniconda2/envs/dDocent/bin/dDocent: line 1360: ee: command not found
#/home/celphin/miniconda2/envs/dDocent/bin/dDocent: line 1361: syntax error near unexpected token `fi'
#/home/celphin/miniconda2/envs/dDocent/bin/dDocent: line 1361: `fi'

#####################################
# try mapping of saximontana and mertensiana with a bit more relaxed scores

# create tmux session (see blog)
tmux new-session -s Cassiope
tmux attach-session -t Cassiope

salloc -c48 --time 23:00:00 --mem 187G --account def-rieseber

cd /project/6019339/celphin/Cassiope/Feb2023_dDocent/outgroup

source ~/miniconda2/bin/activate dDocent 

dDocent

20 individuals are detected. Is this correct? Enter yes or no and press [ENTER]
yes
Proceeding with 20 individuals

dDocent detects 48 processors available on this system.
Please enter the maximum number of processors to use for this analysis.
48

dDocent detects 187 gigabytes of maximum memory available on this system.
Please enter the maximum memory to use for this analysis in gigabytes
For example, to limit dDocent to ten gigabytes, enter 10
This option does not work with all distributions of Linux.  If runs are hanging at variant calling, enter 0
Then press [ENTER]
180

Do you want to quality trim your reads?
Type yes or no and press [ENTER]?
no

Do you want to perform an assembly?
Type yes or no and press [ENTER].
no

Reference contigs need to be in a file named reference.fasta

Do you want to map reads?  Type yes or no and press [ENTER]
yes

BWA will be used to map reads.  You may need to adjust -A -B and -O parameters for your taxa.
Would you like to enter a new parameters now? Type yes or no and press [ENTER]
yes
Please enter new value for A (match score).  It should be an integer.  Default is 1.
1
Please enter new value for B (mismatch score).  It should be an integer.  Default is 4.
3
Please enter new value for O (gap penalty).  It should be an integer.  Default is 6.
5

Do you want to use FreeBayes to call SNPs?  Please type yes or no and press [ENTER]
no

Please enter your email address.  dDocent will email you when it is finished running.
Don't worry; dDocent has no financial need to sell your email address to spammers.
cassandra.elphinstone@shaw.ca

At this point, all configuration information has been entered and dDocent may take several hours to run.
It is recommended that you move this script to a background operation and disable terminal input and output.
All data and logfiles will still be recorded.
To do this:
Press control and Z simultaneously
Type bg and press enter
Type disown -h and press enter
Now sit back, relax, and wait for your analysis to finish

Using BWA to map reads

Creating alignment intervals

dDocent has finished with an analysis in /project/6019339/celphin/Cassiope/Feb2023_dDocent/outgroup

dDocent started Mon Mar 6 17:13:40 PST 2023

dDocent finished Mon Mar 6 17:26:16 PST 2023

dDocent 2.9.4


####################################################
# SNP calling with FreeBayes 
# http://www.ddocent.com/UserGuide/#snp-calling-customization

#-------------------------
# move all the mapping data into one folder

source ~/miniconda2/bin/deactivate dDocent

cd /project/6019339/celphin/Cassiope/Feb2023_dDocent/

cp -v ./tetragona/Pop* ./SNP_calling
cp -v ./tetragona/reference.fasta.* ./SNP_calling
cp -v ./outgroup/Pop* ./SNP_calling

#------------------------
tmux new-session -s Cassiope
tmux attach-session -t Cassiope

salloc -c48 --time 23:00:00 --mem 187G --account def-rieseber

cd /project/6019339/celphin/Cassiope/Feb2023_dDocent/SNP_calling/

source ~/miniconda2/bin/activate dDocent 

dDocent

371 individuals are detected. Is this correct? Enter yes or no and press [ENTER]
yes
Proceeding with 371 individuals

dDocent detects 48 processors available on this system.
Please enter the maximum number of processors to use for this analysis.
48

Do you want to quality trim your reads?
Type yes or no and press [ENTER]?
no

Do you want to perform an assembly?
Type yes or no and press [ENTER].
no

Reference contigs need to be in a file named reference.fasta

Do you want to map reads?  Type yes or no and press [ENTER]
no

Do you want to use FreeBayes to call SNPs?  Please type yes or no and press [ENTER]
yes

Please enter your email address.  dDocent will email you when it is finished running.
Don't worry; dDocent has no financial need to sell your email address to spammers.
cassandra.elphinstone@shaw.ca


At this point, all configuration information has been entered and dDocent may take several hours to run.
It is recommended that you move this script to a background operation and disable terminal input and output.
All data and logfiles will still be recorded.
To do this:
Press control and Z simultaneously
Type 'bg' without the quotes and press enter
Type 'disown -h' again without the quotes and press enter

Now sit back, relax, and wait for your analysis to finish

Genomic interval creation completed successfully.

Using FreeBayes to call SNPs
Using existing popmap file
100% 102:0=0s 7

Using VCFtools to parse TotalRawSNPS.vcf for SNPs that are called in at least 90% of individuals

dDocent has finished with an analysis in /project/6019339/celphin/Cassiope/Feb2023_dDocent/SNP_calling

dDocent started Tue Mar 7 10:36:45 PST 2023

dDocent finished Tue Mar 7 14:03:39 PST 2023

After filtering, kept 96 796 out of a possible 559 815 Sites

dDocent 2.9.4

# 4 hours

###################################################
# https://www.cyberciti.biz/faq/how-to-tar-a-file-in-linux-using-command-line/
# tar older folders

tar -zcvf Jan2023_BLAST_Kluane.tar.gz ./Jan2023_BLAST_Kluane
tar -zcvf Oct2022_dadi.tar.gz ./Oct2022_dadi

mkdir Pre-2023_analysis
# move all tar folders over

mkdir SNP_filtering_March2023
mkdir PopStats_Splitstree_March2023
mkdir ABBA_BABA_March2023
mkdir TreeMix_March2023

##############################################
# try SNP calling and mapping of outgroup again
# try mapping of mertensiana with more relaxed scores
# moved GEN indiv to own folder with logs from mapping run above

# create tmux session (see blog)
tmux new-session -s Cassiope
tmux attach-session -t Cassiope

salloc -c48 --time 23:00:00 --mem 187G --account def-rieseber

cd /project/6019339/celphin/Cassiope/Feb2023_dDocent/outgroup

source ~/miniconda2/bin/activate dDocent 

dDocent

10 individuals are detected. Is this correct? Enter yes or no and press [ENTER]
yes
Proceeding with 10 individuals

dDocent detects 48 processors available on this system.
Please enter the maximum number of processors to use for this analysis.
48

dDocent detects 187 gigabytes of maximum memory available on this system.
Please enter the maximum memory to use for this analysis in gigabytes
For example, to limit dDocent to ten gigabytes, enter 10
This option does not work with all distributions of Linux.  If runs are hanging at variant calling, enter 0
Then press [ENTER]
180

Do you want to quality trim your reads?
Type yes or no and press [ENTER]?
no

Do you want to perform an assembly?
Type yes or no and press [ENTER].
no

Reference contigs need to be in a file named reference.fasta

Do you want to map reads?  Type yes or no and press [ENTER]
yes

BWA will be used to map reads.  You may need to adjust -A -B and -O parameters for your taxa.
Would you like to enter a new parameters now? Type yes or no and press [ENTER]
yes
Please enter new value for A (match score).  It should be an integer.  Default is 1.
3
Please enter new value for B (mismatch score).  It should be an integer.  Default is 4.
3
Please enter new value for O (gap penalty).  It should be an integer.  Default is 6.
4

Do you want to use FreeBayes to call SNPs?  Please type yes or no and press [ENTER]
no

Please enter your email address.  dDocent will email you when it is finished running.
Don't worry; dDocent has no financial need to sell your email address to spammers.
cassandra.elphinstone@shaw.ca

At this point, all configuration information has been entered and dDocent may take several hours to run.
It is recommended that you move this script to a background operation and disable terminal input and output.
All data and logfiles will still be recorded.
To do this:
Press control and Z simultaneously
Type bg and press enter
Type disown -h and press enter
Now sit back, relax, and wait for your analysis to finish

#----------------------
# are the bam files larger?

# old run
-rw-r----- 1 celphin rpp-rieseber    33544397 Mar  6 23:34 PopHAR_21-RG.bam
-rw-r----- 1 celphin rpp-rieseber    30674044 Mar  6 23:34 PopHAR_27-RG.bam
-rw-r----- 1 celphin rpp-rieseber    62374752 Mar  6 23:34 PopHAR_30-RG.bam
-rw-r----- 1 celphin rpp-rieseber    69817099 Mar  6 23:34 PopHAR_6-RG.bam


grep 'PopHAR_*' Cassiope_miss90_mac2_Q30_rmLD.imiss
PopHAR_13       11368   0       1543    0.135732
PopHAR_1        11368   0       3334    0.293279
PopHAR_21       11368   0       2861    0.251671
PopHAR_27       11368   0       893     0.0785538
PopHAR_30       11368   0       96      0.00844476
PopHAR_6        11368   0       77      0.0067734
PopHAR_8        11368   0       3073    0.27032

# if bam files are larger try rerunning SNP calling - not much larger...
# 2,3,4

-rw-r----- 1 celphin rpp-rieseber  34520456 Mar  7 17:19 PopHAR_21-RG.bam
-rw-r----- 1 celphin rpp-rieseber  31423358 Mar  7 17:19 PopHAR_27-RG.bam
-rw-r----- 1 celphin rpp-rieseber  63939654 Mar  7 17:20 PopHAR_30-RG.bam
-rw-r----- 1 celphin rpp-rieseber  71305339 Mar  7 17:20 PopHAR_6-RG.bam

# try again with even less strict values 3,3,4
ls -la *-RG-RG.bam
-rw-r----- 1 celphin rpp-rieseber 35424588 Mar  7 18:02 PopHAR_21-RG-RG.bam
-rw-r----- 1 celphin rpp-rieseber 32076436 Mar  7 18:02 PopHAR_27-RG-RG.bam
-rw-r----- 1 celphin rpp-rieseber 63939654 Mar  7 18:02 PopHAR_30-RG-RG.bam
-rw-r----- 1 celphin rpp-rieseber 71305339 Mar  7 18:02 PopHAR_6-RG-RG.bam

# try again 3,3,3

-rw-r----- 1 celphin rpp-rieseber 35424588 Mar  7 18:32 PopHAR_21-RG-RG-RG.bam
-rw-r----- 1 celphin rpp-rieseber 32076436 Mar  7 18:32 PopHAR_27-RG-RG-RG.bam
-rw-r----- 1 celphin rpp-rieseber 65089545 Mar  7 18:32 PopHAR_30-RG-RG-RG.bam
-rw-r----- 1 celphin rpp-rieseber 72449946 Mar  7 18:32 PopHAR_6-RG-RG-RG.bam

# worth remapping??
# no
####################################################

tmux new-session -s Cassiope3
tmux attach-session -t Cassiope3

# zip/tar files not needed now
tar cfvz Feb2023_dDocent.tar.gz ./Feb2023_dDocent/*
tar cfvz Apr2022_Genbank_upload.tar.gz ./Apr2022_Genbank_upload/*

# remove the full directories
rm -r Feb2023_dDocent
rm -r Apr2022_Genbank_upload

