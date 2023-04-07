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
14
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

