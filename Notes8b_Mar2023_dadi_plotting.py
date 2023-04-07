########################
# Entire Cassiope analysis - part 8b - Demography
# dadi - general plotting/stats code in python
# March 2023
#########################
# easySFS output data is in:
cd /home/celphin/scratch/Cassiope/Fastsimcoal2_Mar2023/5pop_mer/dadi
cd /home/celphin/scratch/Cassiope/Fastsimcoal2_Mar2023/6pop_mer/dadi
cd /home/celphin/scratch/Cassiope/Fastsimcoal2_Mar2023/5pop_tet/dadi
cd /home/celphin/scratch/Cassiope/Fastsimcoal2_Mar2023/6pop_mer/dadi

##########################
# dadi setup

cd /home/celphin/scratch/Cassiope/Fastsimcoal2_Mar2023

wget https://bitbucket.org/gutenkunstlab/dadi/get/ad9de6f61a9c.zip
unzip ad9de6f61a9c.zip

mv gutenkunstlab-dadi-ad9de6f61a9c gutenkunstlab-dadi

#-------------------
# Python environment
module load StdEnv/2020
module load python/3.10.2
module load scipy-stack/2022a
module load nlopt/2.7.0
module load ipykernel/2022a

python3 -m venv /home/celphin/scratch/Cassiope/dadi_vir_env/
source /home/celphin/scratch/Cassiope/dadi_vir_env/bin/activate
# deactivate # to close

#-----------------
# install

python3 -m pip install dadi


################################
# start python session

tmux new-session -s dadi
tmux attach-session -t dadi

salloc -c1 --time 11:50:00 --mem 120000m --account def-cronk

module load StdEnv/2020
module load python/3.10.2
module load scipy-stack/2022a
module load nlopt/2.7.0
module load ipykernel/2022a
source /home/celphin/projects/rpp-rieseber/celphin/Cassiope/Oct2022_dadi/dadi_vir_env/bin/activate

cd /home/celphin/scratch/Cassiope/Fastsimcoal2_Mar2023/

python3

################################################
import os
import pickle
import nlopt
import matplotlib.pyplot as plt
import numpy
import scipy
import dadi
import csv
import matplotlib.pyplot as pyplot

os.getcwd()

##############################
# Not used - SFS made in easySFS

# make the data dictionary
#dd = dadi.Misc.make_data_dict_vcf("example.vcf.gz", "popfile.txt")
dd_5mer = dadi.Misc.make_data_dict_vcf("Cassiope_r30i.recode.vcf", "list_pops_5pop_mer.txt")
dd_6mer = dadi.Misc.make_data_dict_vcf("Cassiope_r30i.recode.vcf", "list_pops_6pop_mer.txt")
dd_5tet = dadi.Misc.make_data_dict_vcf("Cassiope_r30i.recode.vcf", "list_pops_5pop_tet.txt")
dd_6tet = dadi.Misc.make_data_dict_vcf("Cassiope_r30i.recode.vcf", "list_pops_6pop_tet.txt")

#-------------------------------------
# import the population file

pop_ids_5mer = ["Russia","Alaska","Europe","BCmer","BCtet"]
pop_ids_6mer = ["Russia", "Alaska", "Europe", "BCmer", "Greenland"]
pop_ids_5tet = ["Russia", "Alaska", "Europe", "NWT", "Nunavut"]
pop_ids_6tet = ["Greenland", "Alaska", "Europe", "NWT", "Nunavut"]

#-----------------------------
# make the site frequency spectrum

# fs = dadi.Spectrum.from_data_dict(data_dict, pop_ids, projections, mask_corners=True, polarized=True)
fs_5mer = dadi.Spectrum.from_data_dict(dd_5mer, pop_ids_5mer, projections= [10] * 5, mask_corners=True, polarized=False)
fs_6mer = dadi.Spectrum.from_data_dict(dd_6mer, pop_ids_6mer, projections= [10] * 5, mask_corners=True, polarized=True)
fs_5tet = dadi.Spectrum.from_data_dict(dd_5tet, pop_ids_5tet , projections= [10] * 5, mask_corners=True, polarized=True)
fs_6tet = dadi.Spectrum.from_data_dict(dd_6tet, pop_ids_6tet, projections= [10] * 5, mask_corners=True, polarized=True)

# download
#fs.to_file('dadi_admix_pop_frequency_spectrum.fs')

#------------------------------------------
# make site specific fs
Russia_20 = dadi.Spectrum.from_data_dict(dd_5mer, ['Russia'], projections = [20], mask_corners=True, polarized = False)
Alaska_20 = dadi.Spectrum.from_data_dict(dd_5mer, ['Alaska'], projections = [20], mask_corners=True, polarized = False)
Europe_20 = dadi.Spectrum.from_data_dict(dd_5mer, ['Europe'], projections = [20], mask_corners=True, polarized = False)
BCmer_6 = dadi.Spectrum.from_data_dict(dd_5mer, ['BCmer'], projections = [12], mask_corners=True, polarized = False)
BCtet_9 = dadi.Spectrum.from_data_dict(dd_5mer, ['BCtet'], projections = [18], mask_corners=True, polarized = False)
Greenland_20 = dadi.Spectrum.from_data_dict(dd_6mer, ['Greenland'], projections = [20], mask_corners=True, polarized = False)
NWT_20 = dadi.Spectrum.from_data_dict(dd_5tet, ['NWT'], projections = [20], mask_corners=True, polarized = False)
Nunavut_20 = dadi.Spectrum.from_data_dict(dd_5tet, ['Nunavut'], projections = [20], mask_corners=True, polarized = False)


#------------------------
# multi poppulation fs

Alaska_Russia = dadi.Spectrum.from_data_dict(dd_5mer, ['Alaska', 'Russia'], projections = [20,20], mask_corners=True, polarized = False)
Alaska_Europe = dadi.Spectrum.from_data_dict(dd_5mer, ['Alaska', 'Europe'], projections = [20,20], mask_corners=True, polarized = False)
Alaska_BCmer = dadi.Spectrum.from_data_dict(dd_5mer, ['Alaska', 'BCmer'], projections = [6,6], mask_corners=True, polarized = False)
Alaska_BCtet = dadi.Spectrum.from_data_dict(dd_5mer, ['Alaska', 'BCtet'], projections = [9,9], mask_corners=True, polarized = False)
Alaska_NWT = dadi.Spectrum.from_data_dict(dd_5tet, ['Alaska','NWT'], projections = [20,20], mask_corners=True, polarized = False)
Alaska_Nunavut = dadi.Spectrum.from_data_dict(dd_5tet, ['Alaska', 'Nunavut'], projections = [20,20], mask_corners=True, polarized = False)
Alaska_Greenland = dadi.Spectrum.from_data_dict(dd_6mer, ['Alaska', 'Greenland'], projections = [20,20], mask_corners=True, polarized = False)

Europe_Greenland = dadi.Spectrum.from_data_dict(dd_6mer, ['Europe', 'Greenland'], projections = [20,20], mask_corners=True, polarized = False)
Europe_Nunavut = dadi.Spectrum.from_data_dict(dd_5tet, ['Europe', 'Nunavut'], projections = [20,20], mask_corners=True, polarized = False)
Europe_NWT = dadi.Spectrum.from_data_dict(dd_5tet, ['Europe', 'NWT'], projections = [20,20], mask_corners=True, polarized = False)
Europe_Russia = dadi.Spectrum.from_data_dict(dd_5mer, ['Europe', 'Russia'], projections = [20,20], mask_corners=True, polarized = False)

BCmer_BCtet = dadi.Spectrum.from_data_dict(dd_5mer, ['BCmer', 'BCtet'], projections = [6,6], mask_corners=True, polarized = False)

#################################
# here upload data from easy SFS
#fs_5mer = dadi.Spectrum.from_file('./5pop_mer/dadi/Russia-Alaska-Europe-BCmer-BCtet.sfs')
#fs_6mer = dadi.Spectrum.from_file('./6pop_mer/dadi/Russia-Alaska-Europe-BCmer-Greenland.sfs')
#fs_5tet = dadi.Spectrum.from_file('./5pop_tet/dadi/Russia-Alaska-Europe-NWT-Nunavut.sfs')
#fs_6tet = dadi.Spectrum.from_file('./6pop_tet/dadi/Greenland-Alaska-Europe-NWT-Nunavut.sfs')

# Russia_20 = dadi.Spectrum.from_file('./5pop_mer/dadi/Russia-20.sfs')
# Alaska_20 = dadi.Spectrum.from_file('./5pop_mer/dadi/Alaska-20.sfs')
# Europe_20 = dadi.Spectrum.from_file('./5pop_mer/dadi/Europe-20.sfs')
# BCmer_6 = dadi.Spectrum.from_file('./5pop_mer/dadi/BCmer-12.sfs')
# BCtet_9 = dadi.Spectrum.from_file('./5pop_mer/dadi/BCtet-18.sfs')
# Greenland_20 = dadi.Spectrum.from_file('./6pop_mer/dadi/Greenland-20.sfs')
# NWT_20 = dadi.Spectrum.from_file('./5pop_tet/dadi/NWT-20.sfs')
# Nunavut_20 = dadi.Spectrum.from_file('./5pop_tet/dadi/Nunavut-20.sfs')

# Alaska_Russia = dadi.Spectrum.from_file('./5pop_mer/dadi/Russia-Alaska.sfs')
# Alaska_Europe = dadi.Spectrum.from_file('./5pop_mer/dadi/Alaska-Europe.sfs')
# Alaska_BCmer = dadi.Spectrum.from_file('./5pop_mer/dadi/Alaska-BCmer.sfs')
# Alaska_BCtet = dadi.Spectrum.from_file('./5pop_mer/dadi/Alaska-BCtet.sfs')
# Alaska_NWT = dadi.Spectrum.from_file('./5pop_tet/dadi/Alaska-NWT.sfs')
# Europe_Greenland = dadi.Spectrum.from_file('./6pop_mer/dadi/Europe-Greenland.sfs')
# Europe_Nunavut = dadi.Spectrum.from_file('./5pop_tet/dadi/Europe-Nunavut.sfs')

#################################################
# Pop stats
help(dadi.Spectrum)
# Spectrum(data, mask=False, mask_corners=True, data_folded=None, check_folding=True, dtype=<class 'float'>, copy=True, fill_value=nan, keep_mask=True, shrink=True, pop_ids=None, extrap_x=None)

#------------------------
# single popualtion fs
# https://dadi.readthedocs.io/en/latest/user-guide/plotting/
# https://matplotlib.org/stable/tutorials/colors/colormaps.html 

def single_pop(SFS):
  dadi.Spectrum.mask_corners(SFS)
  Theta = dadi.Spectrum.Watterson_theta(SFS)
  Pi = dadi.Spectrum.pi(SFS)
  TajimaD = dadi.Spectrum.Tajima_D(SFS) 
  Sites = dadi.Spectrum.S(SFS)
  fig = plt.figure(219033)
  fig.clear()
  dadi.Plotting.plot_1d_fs(SFS, fig_num=None, show=True)
  name = SFS.pop_ids
  plt.savefig('%s.jpg' % name)
  return Sites, Theta, Pi, TajimaD
  

def two_pop(SFS):
  dadi.Spectrum.mask_corners(SFS)
  Sites = dadi.Spectrum.S(SFS)
  FST = dadi.Spectrum.Fst(SFS)
  fig = plt.figure(219033)
  fig.clear()
  dadi.Plotting.plot_single_2d_sfs(SFS, vmin = 0.1, cmap='gist_rainbow')
  name = SFS.pop_ids
  plt.savefig('%s.jpg' % name)
  return Sites, FST

#------------------------
# run code for single SFS
Russia_20_data = single_pop(Russia_20)
Alaska_20_data = single_pop(Alaska_20)
Europe_20_data = single_pop(Europe_20)
BCmer_12_data = single_pop(BCmer_6)
BCtet_18_data = single_pop(BCtet_9)
NWT_20_data = single_pop(NWT_20)
Nunavut_20_data = single_pop(Nunavut_20)
Greenland_20_data = single_pop(Greenland_20)

Russia_20_data
Alaska_20_data
Europe_20_data
BCmer_12_data 
BCtet_18_data
NWT_20_data
Nunavut_20_data
Greenland_20_data 

>>> Russia_20_data
(2737.6050065739732, 771.6476605214161, 849.1951862526107, 0.4215392367223979)
>>> Alaska_20_data
(3461.2383011999336, 975.6178963781713, 976.1984314093639, 0.002496380750770671)
>>> Europe_20_data
(2549.7848822398755, 718.706874983248, 761.2085381135562, 0.24803758660780437)
>>> BCmer_12_data
(1606.0, 531.8096785368709, 603.8030303030303, 0.640543616494993)
>>> BCtet_18_data
(592.0, 172.11541213665922, 214.98692810457516, 1.0640592456222735)
>>> NWT_20_data
(3163.2605093303687, 891.627011852707, 912.6994352782726, 0.09914421627458395)
>>> Nunavut_20_data
(2960.19427019545, 834.3888098538577, 858.569248298491, 0.12156570514233783)
>>> Greenland_20_data
(1803.7460208710977, 508.4211907260721, 610.9205815553302, 0.845290430474282)


#--------------------------------
# run code for multiple SFS

Alaska_Russia_data = two_pop(Alaska_Russia)
Alaska_Europe_data = two_pop(Alaska_Europe)
Alaska_BCmer_data = two_pop(Alaska_BCmer)
Alaska_BCtet_data = two_pop(Alaska_BCtet)
Alaska_NWT_data = two_pop(Alaska_NWT)
Alaska_Nunavut_data = two_pop(Alaska_Nunavut)
Alaska_Greenland_data = two_pop(Alaska_Greenland)

Europe_Nunavut_data = two_pop(Europe_Nunavut)
Europe_Greenland_data = two_pop(Europe_Greenland)
Europe_NWT_data = two_pop(Europe_NWT)
Europe_Russia_data = two_pop(Europe_Russia)

BCmer_BCtet_data = two_pop(BCmer_BCtet)

Alaska_Russia_data
Alaska_Europe_data 
Alaska_BCmer_data 
Alaska_BCtet_data
Alaska_NWT_data 
Europe_Nunavut_data
Europe_Greenland_data 
Alaska_Nunavut_data
Europe_NWT_data
Alaska_Greenland_data
Europe_Russia_data

>>> Alaska_Russia_data
(4480.038790023957, 0.2202639095077612)
>>> Alaska_Europe_data
(4029.4609091930633, 0.1383667563974237)
>>> Alaska_BCmer_data
(5833.401499299587, 0.8456753234104071)
>>> Alaska_BCtet_data
(4097.475526848526, 0.7595698709141386)
>>> Alaska_NWT_data
(4011.292204815832, 0.060395899378273554)
>>> Europe_Nunavut_data
(3534.694915071199, 0.10354201257817025)
>>> Europe_Greenland_data
(2929.905676619488, 0.2050600887188706)
>>> Alaska_Nunavut_data
(4021.867495177445, 0.09563054126015441)
>>> Europe_NWT_data
(3767.3529527491887, 0.13925094134964885)
>>> Alaska_Greenland_data
(3756.8057656925175, 0.20333136905156918)
>>> Europe_Russia_data
(4004.5093656273007, 0.3000457267705876)
>>> BCmer_BCtet_data
(5068.3594465991555, 0.8962311592119808)
