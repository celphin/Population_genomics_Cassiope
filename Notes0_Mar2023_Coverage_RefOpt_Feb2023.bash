# Feb 2023
# Check coverage of old assembly and optiminze reference assembly parameters

##############################################
# get all bam files to check mapping
mkdir Dec2019_Mapping_bam_files
cd /home/celphin/projects/rpp-rieseber/celphin/Cassiope/Dec2019_Mapping_bam_files

tmux new-session -s fastsimcoal2
tmux attach-session -t fastsimcoal2

cd /home/celphin/projects/rpp-rieseber/celphin/Cassiope/
tar -xzvf Dec2019_demultiplex_dDocent.tgz --wildcards project/6003374/celphin/Cassiope/Dec2019_demultiplex_dDocent/Dec2019_Mapping/*

mv Dec2019_Mapping /home/celphin/projects/rpp-rieseber/celphin/Cassiope/

#-------------------
# check mapping and coverage for each
cd /home/celphin/projects/rpp-rieseber/celphin/Cassiope/Dec2019_Mapping/

# https://www.makeuseof.com/how-to-use-sort-in-linux/

# cov.stats file in .bed format that contains overall levels of coverage across reference contigs
sort -k4 -rn cov.stats > cov.stats_sorted
more cov.stats_sorted
dDocent_Contig_9597     0       272     983658
dDocent_Contig_1391     0       249     977974
dDocent_Contig_23574    0       245     948466
dDocent_Contig_24780    0       201     911186
dDocent_Contig_5208     1       216     811530


cd /home/celphin/projects/rpp-rieseber/celphin/Cassiope/Dec2019_Mapping/cov_stats

sort -k4 -rn PopFOS_3.cov.stats > PopFOS_3.cov.stats_sorted
more PopFOS_3.cov.stats_sorted

dDocent_Contig_24357    0       201     11142
dDocent_Contig_1391     0       249     7466
dDocent_Contig_22484    103     233     7160
dDocent_Contig_22484    1       97      7148


sort -k4 -rn PopGEN_9.cov.stats > PopGEN_9.cov.stats_sorted
more PopGEN_9.cov.stats_sorted

dDocent_Contig_9597     0       272     4737
dDocent_Contig_1391     0       249     3872
dDocent_Contig_5208     1       216     3669
dDocent_Contig_23574    0       245     2319



sort -k4 -rn PopGEN_9.cov.stats > PopGEN_9.cov.stats_sorted
more PopGEN_9.cov.stats_sorted

# sort all cov files by first column
cd /home/celphin/projects/rpp-rieseber/celphin/Cassiope/Dec2019_Mapping/cov_stats
for i in *.cov.stats
do
sort ${i} > ${i}_sorted
done 

sort ../cov.stats > cov.stats_sorted

# join cov files by contig ids
cat cov.stats_sorted > total_cov_stats 
sed -i '1i Contig\tStart\tEnd\tTotal' total_cov_stats 
#i="PopZAC_28.cov.stats_sorted"

for i in *.cov.stats_sorted
do 
echo ${i}
awk '{print $4}' ${i} > file_column
sed "1i ${i}" file_column > add_column
paste -d'\t' total_cov_stats add_column > total_plus
mv total_plus total_cov_stats
done

more total_cov_stats

# sort by last column - total coverage
sort -k4 -rn total_cov_stats > total_cov_stats_sorted
more total_cov_stats_sorted


#------------------------
# plot in R with imagenan
tmux new-session -s R
tmux attach-session -t R

module load StdEnv/2020 r/4.1.2
R

library(RepeatObserver)
library(Hmisc)
library(seqinr)
library(vcfR)
library(ape)
library(ade4)
library(dendextend) #colouring dendogram
library(matrixStats)
library(pracma)
library(stringr)
library(parallel)
library(dplyr)

setwd ("/home/celphin/projects/rpp-rieseber/celphin/Cassiope/Dec2019_Mapping/cov_stats/")

mydata <- read.table("./total_cov_stats", header=TRUE)

mydata[1:5, 1:10]

--
mydata1 <- mydata[, -c(1:4)]

mydata2 <- mydata1
colnames(mydata2) <- sub(".cov.stats_sorted", "", colnames(mydata2))
colnames(mydata2) <- sub("Pop", "", colnames(mydata2))

mydata3 <- arrange(mydata2, desc(QHI_32))
mydata4 <- mydata3[1:10000,]
jpeg("./QHI_Cassiope_Coverage.jpg", width = 7000, height = 2000)
imagenan(mydata4, zlim=c(0,50), yline=1,yma=15,xline=3,xma=20,lnumr=39,lnumc=38,lasval=2,cex.axis=6, col = topo.colors(255),outside.below.color='black',outside.above.color='red',na.color='gray')
dev.off()

mydata3 <- arrange(mydata2, desc(SW_31))
mydata4 <- mydata3[1:10000,]
jpeg("./SW_Cassiope_Coverage.jpg", width = 7000, height = 2000)
imagenan(mydata4, zlim=c(0,50), yline=1,yma=15,xline=3,xma=20,lnumr=39,lnumc=38,lasval=2,cex.axis=6, col = topo.colors(255),outside.below.color='black',outside.above.color='red',na.color='gray')
dev.off()

mydata3 <- arrange(mydata2, desc(ZAC_20))
mydata4 <- mydata3[1:10000,]
jpeg("./ZAC_Cassiope_Coverage.jpg", width = 7000, height = 2000)
imagenan(mydata4, zlim=c(0,50), yline=1,yma=15,xline=3,xma=20,lnumr=39,lnumc=38,lasval=2,cex.axis=6, col = topo.colors(255),outside.below.color='black',outside.above.color='red',na.color='gray')
dev.off()


mydata3 <- arrange(mydata2, desc(GEN_3))
mydata4 <- mydata3[1:10000,]
jpeg("./GEN3_Cassiope_Coverage.jpg", width = 7000, height = 2000)
imagenan(mydata4, zlim=c(0,50), yline=1,yma=15,xline=3,xma=20,lnumr=39,lnumc=38,lasval=2,cex.axis=6, col = topo.colors(255),outside.below.color='black',outside.above.color='red',na.color='gray')
dev.off()

mydata3 <- arrange(mydata2, desc(GEN_4))
mydata4 <- mydata3[1:10000,]
jpeg("./GEN4_Cassiope_Coverage.jpg", width = 7000, height = 2000)
imagenan(mydata4, zlim=c(0,50), yline=1,yma=15,xline=3,xma=20,lnumr=39,lnumc=38,lasval=2,cex.axis=6, col = topo.colors(255),outside.below.color='black',outside.above.color='red',na.color='gray')
dev.off()


mydata3 <- arrange(mydata2, desc(HAR_13))
mydata4 <- mydata3[1:10000,]
jpeg("./HAR_Cassiope_Coverage.jpg", width = 7000, height = 2000)
imagenan(mydata4, zlim=c(0,50), yline=1,yma=15,xline=3,xma=20,lnumr=39,lnumc=38,lasval=2,cex.axis=6, col = topo.colors(255),outside.below.color='black',outside.above.color='red',na.color='gray')
dev.off()

mydata3 <- arrange(mydata2, desc(AlexOld_33))
mydata4 <- mydata3[1:10000,]
jpeg("./Alex_old_Cassiope_Coverage.jpg", width = 7000, height = 2000)
imagenan(mydata4, zlim=c(0,50), yline=1,yma=15,xline=3,xma=20,lnumr=39,lnumc=38,lasval=2,cex.axis=6, col = topo.colors(255),outside.below.color='black',outside.above.color='red',na.color='gray')
dev.off()

#############################
# rerun dDocent with new assembly

# https://www.ddocent.com/assembly/

# first sart dDocent and check old assembly
cd /home/celphin/projects/rpp-rieseber/celphin/Cassiope/
mkdir Feb2023_dDocent; cd Feb2023_dDocent

mv /home/celphin/projects/rpp-rieseber/celphin/Cassiope/Dec2019_Mapping/fastq /home/celphin/projects/rpp-rieseber/celphin/Cassiope/Feb2023_dDocent

cd /project/6019339/celphin/Cassiope/Feb2023_dDocent/fastq

mkdir reference_fastq

cp -v PopLAJ*.fq.gz reference_fastq/
cp -v PopSAM*.fq.gz reference_fastq/
cp -v PopATQ*.fq.gz reference_fastq/

cp -v PopPET*.fq.gz reference_fastq/
cp -v PopYED*.fq.gz reference_fastq/
cp -v PopBARD*.fq.gz reference_fastq/

# should be 60 individuals

###################
# check best similarity to cluster value (used 0.9 previously)
# https://www.ddocent.com/assembly/

# create tmux session (see blog)
tmux new-session -s Cassiope
tmux attach-session –t Cassiope

# get interactive allocation
salloc -c1 --time 2:50:00 --mem 120000m --account def-rieseber

cd /project/6019339/celphin/Cassiope/Feb2023_dDocent/fastq/reference_fastq/

ls *.F.fq.gz > namelist
sed -i'' -e 's/.F.fq.gz//g' namelist
AWK1='BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}'
AWK2='!/>/'
AWK3='!/NNN/'
PERLT='while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}'

cat namelist | parallel --no-notice -j 8 "zcat {}.F.fq.gz | mawk '$AWK1' | mawk '$AWK2' > {}.forward"
cat namelist | parallel --no-notice -j 8 "zcat {}.R.fq.gz | mawk '$AWK1' | mawk '$AWK2' > {}.reverse"
cat namelist | parallel --no-notice -j 8 "paste -d '-' {}.forward {}.reverse | mawk '$AWK3' | sed 's/-/NNNNNNNNNN/' | perl -e '$PERLT' > {}.uniq.seqs"

cat *.uniq.seqs > uniq.seqs
for i in {2..20};
do 
echo $i >> pfile
done
cat pfile | parallel --no-notice "echo -n {}xxx && mawk -v x={} '\$1 >= x' uniq.seqs | wc -l" | mawk  '{gsub("xxx","\t",$0); print;}'| sort -g > uniqseq.data
rm pfile

more uniqseq.data

2       2771955
3       1707760
4       1373274
5       1190459
6       1066595
7       972731
8       896482
9       832507
10      776133
11      726674
12      682364
13      642483
14      605845
15      571841
16      540376
17      511328
18      484238
19      458768
20      434908

module load  StdEnv/2020 gnuplot/5.4.2

gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale
set xrange [2:20] 
unset label
set title "Number of Unique Sequences with More than X Coverage (Counted within individuals)"
set xlabel "Coverage"
set ylabel "Number of Unique Sequences"
plot 'uniqseq.data' with lines notitle
pause -1
EOF



                        Number of Unique Sequences with More than X Coverage (Counted within individuals)
       3e+06 +------------------------------------------------------------------------------------------------------+
             |          +           +          +           +          +           +          +           +          |
             |                                                                                                      |
             |*                                                                                                     |
     2.5e+06 |-*                                                                                                  +-|
             | *                                                                                                    |
             |  *                                                                                                   |
             |   *                                                                                                  |
       2e+06 |-+  *                                                                                               +-|
             |    *                                                                                                 |
             |     **                                                                                               |
     1.5e+06 |-+     **                                                                                           +-|
             |         ***                                                                                          |
             |            ***                                                                                       |
             |               ****                                                                                   |
       1e+06 |-+                 ******                                                                           +-|
             |                         ***********                                                                  |
             |                                    ******************                                                |
             |                                                      *****************                               |
      500000 |-+                                                                     **************************** +-|
             |                                                                                                   ***|
             |                                                                                                      |
             |          +           +          +           +          +           +          +           +          |
           0 +------------------------------------------------------------------------------------------------------+
             2          4           6          8           10         12          14         16          18         20
                                                            Coverage

# K1=4 or 5 seems like a reasonable number of times the sequence needs to appear in an individual to be considered real

parallel --no-notice -j 8 mawk -v x=4 \''$1 >= x'\' ::: *.uniq.seqs | cut -f2 | perl -e 'while (<>) {chomp; $z{$_}++;} while(($k,$v) = each(%z)) {print "$v\t$k\n";}' > uniqCperindv
wc -l uniqCperindv

for ((i = 2; i <= 10; i++));
do
echo $i >> ufile
done

cat ufile | parallel --no-notice "echo -n {}xxx && mawk -v x={} '\$1 >= x' uniqCperindv | wc -l" | mawk  '{gsub("xxx","\t",$0); print;}'| sort -g > uniqseq.peri.data
rm ufile


gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Number of Unique Sequences present in more than X Individuals"
set xlabel "Number of Individuals"
set ylabel "Number of Unique Sequences"
plot 'uniqseq.peri.data' with lines notitle
pause -1
EOF


                                  Number of Unique Sequences present in more than X Individuals
     250000 +-------------------------------------------------------------------------------------------------------+
            |            +            +            +            +            +            +            +            |
            |                                                                                                       |
            |**                                                                                                     |
            |  **                                                                                                   |
     200000 |-+  **                                                                                               +-|
            |      **                                                                                               |
            |        **                                                                                             |
            |          **                                                                                           |
     150000 |-+          ***                                                                                      +-|
            |               ****                                                                                    |
            |                   ****                                                                                |
            |                       *****                                                                           |
            |                            ****                                                                       |
     100000 |-+                              ****                                                                 +-|
            |                                    ******                                                             |
            |                                          ******                                                       |
            |                                                **********                                             |
      50000 |-+                                                        *************                              +-|
            |                                                                       *************                   |
            |                                                                                    *************      |
            |                                                                                                 ******|
            |            +            +            +            +            +            +            +            |
          0 +-------------------------------------------------------------------------------------------------------+
            2            3            4            5            6            7            8            9            10
                                                      Number of Individuals


# K2=4 seems like a reasonable number of individuals for the read to need to appear in


mawk -v x=4 '$1 >= x' uniqCperindv > uniq.k.4.c.4.seqs
wc -l uniq.k.4.c.4.seqs

cut -f2 uniq.k.4.c.4.seqs > totaluniqseq
mawk '{c= c + 1; print ">Contig_" c "\n" $1}' totaluniqseq > uniq.fasta

# Note the real dDocent run should be done because it will check for Illumina adapters

sed -e 's/NNNNNNNNNN/\t/g' uniq.fasta | cut -f1 > uniq.F.fasta

source ~/miniconda2/bin/activate dDocent 
cd-hit-est -i uniq.F.fasta -o xxx -c 0.8 -T 0 -M 0 -g 1


   114749  finished      23094  clusters

Apprixmated maximum memory consumption: 934M
writing new database
writing clustering information
program completed !

Total CPU time 768.42

mawk '{if ($1 ~ /Cl/) clus = clus + 1; else  print $3 "\t" clus}' xxx.clstr | sed 's/[>Contig_,...]//g' | sort -g -k1 > sort.contig.cluster.ids
paste sort.contig.cluster.ids totaluniqseq > contig.cluster.totaluniqseq
sort -k2,2 -g contig.cluster.totaluniqseq | sed -e 's/NNNNNNNNNN/\t/g' > rcluster

# Read_ID	Cluster_ID	Forward_Read	Reverse_Read

cut -f2 rcluster | uniq | wc -l 
23,094

#-------------------
# f value here is important, f= min freq of an allele to make it into its own cluster, make higher for more diverse reads
rainbow div -i rcluster -o rbdiv.out -f 0.5 -K 10

# -r parameter, which is the minimum number of reads to assemble (default 5) - for many individuals to make ref a cutoff of 2 is better
rainbow merge -o rbasm.out -a -i rbdiv.out -r 2

# rbasm output lists optimal and suboptimal contigs
cat rbasm.out <(echo "E") |sed 's/[0-9]*:[0-9]*://g' | mawk ' {
if (NR == 1) e=$2;
else if ($1 ~/E/ && lenp > len1) {c=c+1; print ">dDocent_Contig_" e "\n" seq2 "NNNNNNNNNN" seq1; seq1=0; seq2=0;lenp=0;e=$2;fclus=0;len1=0;freqp=0;lenf=0}
else if ($1 ~/E/ && lenp <= len1) {c=c+1; print ">dDocent_Contig_" e "\n" seq1; seq1=0; seq2=0;lenp=0;e=$2;fclus=0;len1=0;freqp=0;lenf=0}
else if ($1 ~/C/) clus=$2;
else if ($1 ~/L/) len=$2;
else if ($1 ~/S/) seq=$2;
else if ($1 ~/N/) freq=$2;
else if ($1 ~/R/ && $0 ~/0/ && $0 !~/1/ && len > lenf) {seq1 = seq; fclus=clus;lenf=len}
else if ($1 ~/R/ && $0 ~/0/ && $0 ~/1/) {seq1 = seq; fclus=clus; len1=len}
else if ($1 ~/R/ && $0 ~!/0/ && freq > freqp && len >= lenp || $1 ~/R/ && $0 ~!/0/ && freq == freqp && len > lenp) {seq2 = seq; lenp = len; freqp=freq}
}' > rainbow.fasta

# resulting contigs need to be aligned and clustered by sequence similarity
# memory usage (-M) and number of threads (-T)
cd-hit-est -i rainbow.fasta -o referenceRC.fasta -M 0 -T 0 -c 0.9

#----------------
# remake reference - using code above
# varying the uniq sequence copy cutoff and the final clustering similarity have the the largest effect on the number of final contigs.

tmux new-session -s Cassiope
tmux attach-session -t Cassiope
salloc -c32 --time 2:00:00 --mem 120000m --account def-rieseber
source ~/miniconda2/bin/activate dDocent 

curl -L -O https://github.com/jpuritz/dDocent/raw/master/scripts/remake_reference.sh

module load StdEnv/2020 fastp/0.23.1

bash remake_reference.sh 4 4 0.90 PE 2 

dDocent assembled 114,749 sequences (after cutoffs) into 22,940 contigs

#---
# examine the reference
head reference.fasta
wc -l reference.fasta
22676 reference.fasta

mawk '/>/' reference.fasta | wc -l 
11338

#mawk '/^NAATT.*N*.*CCGN$/' reference.fasta | wc -l
grep '^NAATT.*N*.*CCGN$' reference.fasta | wc -l
0


#-----------------
# optimize reference % similarity and cutoffs (k1 and k2)
tmux new-session -s Cassiope
tmux attach-session –t Cassiope

# get interactive allocation
salloc -c32 --time 5:00:00 --mem 120000m --account def-rieseber

module load StdEnv/2020 fastp/0.23.1
source ~/miniconda2/bin/activate dDocent 

curl -L -O https://github.com/jpuritz/dDocent/raw/master/scripts/ReferenceOpt.sh 

bash ReferenceOpt.sh 4 8 4 8 PE 16 # takes a long time 


K1 is 7 K2 is 6 c is 0.88
K1 is 7 K2 is 6 c is 0.90
K1 is 7 K2 is 6 c is 0.92
K1 is 7 K2 is 6 c is 0.94
K1 is 7 K2 is 6 c is 0.96
K1 is 7 K2 is 6 c is 0.98
K1 is 7 K2 is 7 c is 0.80
ReferenceOpt.sh: line 259: / 100 + 1: syntax error: operand expected (error token is "/ 100 + 1")
K1 is 7 K2 is 7 c is 0.82
K1 is 7 K2 is 7 c is 0.84
K1 is 7 K2 is 7 c is 0.86
K1 is 7 K2 is 7 c is 0.88
K1 is 7 K2 is 7 c is 0.90
K1 is 7 K2 is 7 c is 0.92
K1 is 7 K2 is 7 c is 0.94
# ran out of time

#----
# restart with more time
salloc -c32 --time 8:00:00 --mem 120000m --account def-rieseber

module load StdEnv/2020 fastp/0.23.1
source ~/miniconda2/bin/activate dDocent 

bash ReferenceOpt.sh 4 8 4 8 PE 32 # takes a long time 

more kopt.data
4 4 0.80 22940
4 4 0.82 22996
4 4 0.84 23051
4 4 0.86 23091
4 4 0.88 23124
4 4 0.90 23156
4 4 0.92 23316
4 4 0.94 23434
4 4 0.96 23547
4 4 0.98 23708
4 5 0.80 19622
4 5 0.82 19665
4 5 0.84 19699
4 5 0.86 19730
4 5 0.88 19757
4 5 0.90 19781
4 5 0.92 19894
4 5 0.94 19995
4 5 0.96 20072
4 5 0.98 20180
4 6 0.80 17475

5 8 0.80 13583  # maybe best mapping coverage but low #contigs
5 8 0.82 13595
5 8 0.84 13605
5 8 0.86 13618
5 8 0.88 13626
5 8 0.90 13634
5 8 0.92 13682
5 8 0.94 13725
5 8 0.96 13754


source ~/miniconda2/bin/activate dDocent 
module load  StdEnv/2020 gnuplot/5.4.2 fastp/0.23.1

cd /project/6019339/celphin/Cassiope/Feb2023_dDocent/fastq/reference_fastq/

gnuplot << \EOF
set terminal dumb size 120, 30
set autoscale
unset label
set title "Histogram of number of reference contigs"
set ylabel "Number of Occurrences"
set xlabel "Number of reference contigs"
max = `sort -g plot.kopt.data | tail -1`
binwidth = max/250.0
bin(x,width)=width*floor(x/width) + binwidth/2.0
#set xtics 10
plot 'plot.kopt.data' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF

                                         Histogram of number of reference contigs
     7 +------------------------------------------------------------------------------------------------------------+
       |          *    +          *** +               +           *  +               +              +               |
       |          *               ***                       'plot.*opt.data' using (bin($1,binwidth)):(1.0) ******* |
     6 |-+        * ****    **    ***         **                  *                                               +-|
       |          * *  *    **    ***         **                  *                                                 |
       |          * *  *    **    ***         **                  *                                                 |
       |          * *  *    **    ***         **                  *                                                 |
     5 |-+        * *  *  **** ** **********  ***  ***    **      *                                               +-|
       |          * *  *  **** ** ******   *  ***  ***    **      *                                                 |
       |          * *  *  **** ** ******   *  ***  ***    **      *                                                 |
     4 |-+        * *  *  **** ** ******** ******  *****  ****    * **                                            +-|
       |          * *  *  **** ** ******** ******  *****  ****    * **                                              |
       |          * *  *  **** ** ******** ******  *****  ****    * **                                              |
     3 |-+        ***  ** ************************ *****  ************   ***       *   **                    *    +-|
       |          * *  ** ************************ *****  ***** * ****   * *       *   **                    *      |
       |          * *  ** ************************ *****  ***** * ****   * *       *   **                    *      |
     2 |-+        * *  *************************** ****** ***** * *****  * *      **  ***                    *    +-|
       |          * *  *************************** ****** ***** * *****  * *      **  ***                    *      |
       |          * *  *************************** ****** ***** * *****  * *      **  ***                    *      |
       |          * *  *************************** ****** ***** * *****  * *      **  ***                    *      |
     1 |-+        * *  ****************************************** ************************************************+-|
       |          * *  *************************** ************ * ****** * ****   **********       *        ******  |
       |          * *  *************************** ************ * ****** * ****   **********       *+       ******  |
     0 +------------------------------------------------------------------------------------------------------------+
     10000           12000          14000           16000          18000           20000          22000           24000
                                                Number of reference contigs

Average contig number = 15721.7
The top three most common number of contigs
X       Contig number
2       17577
2       17544
1       23708
The top three most common number of contigs (with values rounded)
X       Contig number
9       16100
8       13700
7       15800

#-----
# what is the best -c 

# probably 80%??
# or does this merge too many different contigs creating far more diversity than is real

# maybe 95% is better separating contigs only if they differ by 8-15bp 


5 8 0.80 13583
5 8 0.82 13595
5 8 0.84 13605
5 8 0.86 13618
5 8 0.88 13626
5 8 0.90 13634
5 8 0.92 13682
5 8 0.94 13725
5 8 0.96 13754
5 8 0.98 13805

6 7 0.80 13914
6 7 0.82 13926
6 7 0.84 13939
6 7 0.86 13953
6 7 0.88 13961
6 7 0.90 13970
6 7 0.92 14028
6 7 0.94 14078
6 7 0.96 14111
6 7 0.98 14166

8 6 0.82 13662
8 6 0.84 13674
8 6 0.86 13689
8 6 0.88 13698
8 6 0.90 13707
8 6 0.92 13768
8 6 0.94 13823
8 6 0.96 13858
8 6 0.98 13904



#-----------------
# optimize mapping
tmux new-session -s Cassiope
tmux attach-session –t Cassiope

# get interactive allocation
salloc -c32 --time 5:00:00 --mem 120000m --account def-rieseber

module load StdEnv/2020 fastp/0.23.1
source ~/miniconda2/bin/activate dDocent 
#conda deactivate
curl -L -O https://raw.githubusercontent.com/jpuritz/WinterSchool.2016/master/RefMapOpt.sh

# RefMapOpt minK1 maxK1 minK2 maxK2 cluster_similarity Assembly_Type Num_of_Processors optional_list_of_individuals
RefMapOpt.sh 4 8 4 8 0.9 PE 32 # takes a long time

 more mapping.results
Cov     Non0Cov Contigs MeanContigsMapped       K1      K2      SUM Mapped      SUM Properly    Mean Mapped     Mean Properly   MisMatched
67.1496 76.251  17503   19338.8 4       4       29355565        28366498        1.46778e+06     1.41832e+06     10362.8
81.1426 86.5919 14530   16087.5 4       5       27717282        26816784        1.38586e+06     1.34084e+06     7087.55
81.0009 87.023  17748   15959.5 4       6       27663670        26763892        1.38318e+06     1.33819e+06     7601.55

[bam_sort_core] merging from 0 files and 32 in-memory blocks...
[bam_sort_core] merging from 0 files and 32 in-memory blocks...
/home/celphin/miniconda2/envs/dDocent/bin/RefMapOpt.sh: line 201: / 100 + 1: syntax error: operand expected (error token is "/ 100 + 1")

#---
# restart with more time
tmux new-session -s Cassiope
tmux attach-session –t Cassiope

# get interactive allocation
salloc -c32 --time 23:00:00 --mem 120000m --account def-rieseber

module load StdEnv/2020 fastp/0.23.1
source ~/miniconda2/bin/activate dDocent 

# RefMapOpt minK1 maxK1 minK2 maxK2 cluster_similarity Assembly_Type Num_of_Processors optional_list_of_individuals
RefMapOpt.sh 4 8 4 8 0.9 PE 32 # takes a long time
# did not finish

RefMapOpt.sh 5 8 5 8 0.9 PE 32 # takes a long time

# The output contains the average coverage per contig, the average coverage per contig not counting zero coverage contigs, 
# the number of contigs, the mean number of contigs mapped, the two cutoff values used, the sum of all mapped reads, the sum of 
# all properly mapped reads, the mean number of mapped reads, the mean number of properly mapped reads, and the number of reads 
# that are mapped to mismatching contigs. Here, we are looking to values that maximize properly mapped reads, the mean number 
# of contigs mapped, and the coverage. 

more mapping.results # for 0.9
Cov     Non0Cov Contigs MeanContigsMapped       K1      K2      SUM Mapped      SUM Properly    Mean Mapped     Mean Properly   MisMatched
84.1015 92.2217 17884   16287.1                 5       5       30083088        29201949        1.50415e+06     1.4601e+06      7928.65
91.3291 97.5153 16159   15115.2                 5       6       29517589        28653428        1.47588e+06     1.43267e+06     6556.35
96.472  101.371 14823   14091.5                 5       7       28602006        27794691        1.4301e+06      1.38973e+06     5217.1
101.374 105.224 13635   13124.2                 5       8       27646742        26889954        1.38234e+06     1.3445e+06      3963
88.9378 96.4815 16691   15363.1                 6       5       29691009        28834804        1.48455e+06     1.44174e+06     6583.55
95.3637 101.269 15219   14313.1                 6       6       29028692        28225357        1.45143e+06     1.41127e+06     5562.4
100.115 104.78  13972   13335.8                 6       7       27978072        27237389        1.3989e+06      1.36187e+06     4027.85
105.622 109.329 12785   12340.7                 6       8       27009726        26328125        1.35049e+06     1.31641e+06     3059.5
92.1775 99.4138 15847   14671.2                 7       5       29216570        28414516        1.46083e+06     1.42073e+06     6045.5
98.6326 104.394 14430   13617                   7       6       28467348        27719389        1.42337e+06     1.38597e+06     4556.4
103.934 108.476 13184   12619.3                 7       7       27407510        26710421        1.37038e+06     1.33552e+06     3265.95
109.726 113.329 12026   11633.5                 7       8       26393606        25751268        1.31968e+06     1.28756e+06     2665.15
95.3817 102.444 15070   14010.8                 8       5       28749978        27983014        1.4375e+06      1.39915e+06     5163.75
102.034 107.567 13710   12989.5                 8       6       27979823        27269466        1.39899e+06     1.36347e+06     3454.1
108.002 112.404 12458   11958.5                 8       7       26911932        26247933        1.3456e+06      1.3124e+06      2805.55
113.725 117.234 11338   10989.4                 8       8       25790517        25182737        1.28953e+06     1.25914e+06     2114.45


#---------------------
# what are the best K1 and K2?
# probably 5 and 8

101.374 105.224 13635   13124.2                 5       8       27646742        26889954        1.38234e+06     1.3445e+06      3963
105.622 109.329 12785   12340.7                 6       8       27009726        26328125        1.35049e+06     1.31641e+06     3059.5
113.725 117.234 11338   10989.4                 8       8       25790517        25182737        1.28953e+06     1.25914e+06     2114.45

#----------------------
# try again with guessed values

tmux new-session -s Cassiope
tmux attach-session -t Cassiope

cd /project/6019339/celphin/Cassiope/Feb2023_dDocent/fastq/reference_fastq/

salloc -c48 --time 23:00:00 --mem 187000m --account def-rieseber

module load StdEnv/2020 fastp/0.23.1
source ~/miniconda2/bin/activate dDocent 

RefMapOpt.sh 5 8 6 10 0.8 PE 48 # takes a long time

more mapping.results # for 0.8 may have been lost when ran 0.95, see above for 0.9


#-----

module load StdEnv/2020 fastp/0.23.1
source ~/miniconda2/bin/activate dDocent 

RefMapOpt.sh 5 8 6 10 0.95 PE 48 # takes a long time

more mapping.results # for 0.95

Cov     Non0Cov Contigs MeanContigsMapped       K1      K2      SUM Mapped      SUM Properly    Mean Mapped     Mean Properly   MisMatched
52.5015 69.462  28504   20934.2                 5       6       29931112        29052582        1.49656e+06     1.45263e+06     16525.7
57.3907 72.3444 26074   20107.8                 5       7       29929251        29092419        1.49646e+06     1.45462e+06     14127.1
60.8474 74.7563 24593   19469.5                 5       8       29929645        29106019        1.49648e+06     1.4553e+06      13473.1
63.8164 76.9617 23441   18918.2                 5       9       29919657        29120814        1.49598e+06     1.45604e+06     12095.8
66.2114 78.8148 22526   18428.2                 5       10      29830888        29053621        1.49154e+06     1.45268e+06     11444.9
57.3144 73.6304 25842   19564                   6       6       29623539        28807232        1.48118e+06     1.44036e+06     13630.4
62.2246 76.6562 23903   18882.4                 6       7       29748298        28952092        1.48741e+06     1.4476e+06      12334.8
65.4936 78.9305 22689   18330.1                 6       8       29721006        28945140        1.48605e+06     1447257 11674.1
68.1264 80.8311 21775   17878                   6       9       29670385        28916151        1.48352e+06     1.44581e+06     10467
70.0994 82.1827 21066   17512.2                 6       10      29535651        28780645        1.47678e+06     1.43903e+06     10493.6
61.1564 77.0677 24041   18571                   7       6       29406467        28603519        1.47032e+06     1.43018e+06     13079.5
65.7927 79.8157 22429   18008.4                 7       7       29514574        28741624        1.47573e+06     1.43708e+06     11277.3
68.8082 81.8175 21410   17544.8                 7       8       29465037        28715563        1.47325e+06     1.43578e+06     10079.6
71.1847 83.4812 20695   17203.5                 7       9       29464780        28721929        1473239 1.4361e+06      9977.1
73.2803 85.0234 20027   16835.5                 7       10      29353140        28623804        1467657 1.43119e+06     10013.7
64.1664 79.7226 22770   17854.8                 8       6       29222648        28445091        1.46113e+06     1.42225e+06     11847.6
68.6408 82.3968 21387   17366                   8       7       29361783        28613737        1.46809e+06     1.43069e+06     9811.15
71.6633 84.2968 20490   16984.3                 8       8       29369036        28628668        1.46845e+06     1.43143e+06     9907.7
73.8492 85.8351 19824   16637.7                 8       9       29281208        28545100        1.46406e+06     1427255 10167.9
75.8058 87.2522 19239   16312                   8       10      29170038        28458352        1.4585e+06      1.42292e+06     9370.85


# 90%
101.374 105.224 13635   13124.2                 5       8       27646742        26889954        1.38234e+06     1.3445e+06      3963
105.622 109.329 12785   12340.7                 6       8       27009726        26328125        1.35049e+06     1.31641e+06     3059.5
113.725 117.234 11338   10989.4                 8       8       25790517        25182737        1.28953e+06     1.25914e+06     2114.45


# The output contains the average coverage per contig, the average coverage per contig not counting zero coverage contigs, 
# the number of contigs, the mean number of contigs mapped, the two cutoff values used, the sum of all mapped reads, the sum of 
# all properly mapped reads, the mean number of mapped reads, the mean number of properly mapped reads, and the number of reads 
# that are mapped to mismatching contigs. Here, we are looking to values that maximize properly mapped reads, the mean number 
# of contigs mapped, and the coverage. 

# K1=4 or 5 seems like a reasonable number of times the sequence needs to appear in an individual to be considered real
# K2=5 seems like a reasonable number of individuals for the read to need to appear in to be real

# Average contig number = 15721.7

###################
# final choice
# 95% similarity to merge
# K1=4 and K2=5


tmux new-session -s Cassiope
tmux attach-session -t Cassiope

cd /project/6019339/celphin/Cassiope/Feb2023_dDocent/fastq/reference_fastq/

salloc -c48 --time 23:00:00 --mem 187000m --account def-rieseber

module load StdEnv/2020 fastp/0.23.1
source ~/miniconda2/bin/activate dDocent 

bash remake_reference.sh 4 5 0.95 PE 2 

dDocent assembled 89 723 sequences (after cutoffs) into 20 044 contigs

#-------------
# maybe want 90% similarity since we want to be able to compare and have overlap where there is overlap
# but previously all populations had -Fis= more heterozygous than expected, possibly due to some different reads mapping to same contig in the reference 
# go with 95% only adds a few more contigs 


bash remake_reference.sh 4 5 0.90 PE 2 

dDocent assembled 89723 sequences (after cutoffs) into 19781 contigs

#---
# examine the reference # old reference had 330 000 contigs!!
head reference.fasta
wc -l reference.fasta 

mawk '/>/' reference.fasta | wc -l 

#mawk '/^NAATT.*N*.*CCGN$/' reference.fasta | wc -l
grep '^NAATT.*N*.*CCGN$' reference.fasta | wc -l

#------------
# https://bio-bwa.sourceforge.net/bwa.shtml
# Mapping 
# what are the best -A -B and -O ?
# stick with default or make more flexible?? reduce a bit from 4 to 3 and from 6 to 5

-A INT 	Matching score. [1]
-B INT 	Mismatch penalty. The sequence error rate is approximately: {.75 * exp[-log(4) * B/A]}. [3]
-O INT 	Gap open penalty. [5] 

