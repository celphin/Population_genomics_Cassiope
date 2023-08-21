########################
# Entire Cassiope analysis - part 8b - Demography
# dadi - general plotting/stats code in python
# March 2023
#########################
# easySFS output data is in:
cd ~/projects/def-rieseber/Dryas_shared_data/Fernando_SFS

scp ./pops* celphin@cedar.computecanada.ca:~/projects/def-rieseber/Dryas_shared_data/Fernando_SFS
scp ./vcf.fsc.nomaf.miss0.1.vcf celphin@cedar.computecanada.ca:~/projects/def-rieseber/Dryas_shared_data/Fernando_SFS

for f in popsCA1.txt popsCA4.txt popsCO.txt popsMA.txt; do sed -i "s/$/\t$f/" $f ; done

cat popsCA1.txt popsCA4.txt popsCO.txt popsMA.txt > popFern.txt

awk '{$2=""; print $0}' popFern.txt > popsFern.txt

sed -i 's/pops//g' popsFern.txt
sed -i 's/\.txt//g' popsFern.txt

##########################
# dadi setup

cd ~/projects/def-rieseber/Dryas_shared_data/Fernando_SFS/

wget https://bitbucket.org/gutenkunstlab/dadi/get/ad9de6f61a9c.zip
unzip ad9de6f61a9c.zip

mv gutenkunstlab-dadi-ad9de6f61a9c gutenkunstlab-dadi

#-------------------
# Python environment
module load StdEnv/2020
module load python/3.10.2

python3 -m venv ~/projects/def-rieseber/Dryas_shared_data/Fernando_SFS/dadi_vir_env/
source ~/projects/def-rieseber/Dryas_shared_data/Fernando_SFS/dadi_vir_env/bin/activate

module load scipy-stack/2022a
module load nlopt/2.7.0
module load ipykernel/2022a

python3 -m pip install dadi

# deactivate # to close vir env

#-----------------
# install

#cd ~/projects/def-rieseber/Dryas_shared_data/Fernando_SFS/gutenkunstlab-dadi/
#python setup.py install


################################
# start python session

tmux new-session -s dadi
tmux attach-session -t dadi

salloc -c1 --time 2:50:00 --mem 120000m --account def-cronk

module load StdEnv/2020
module load python/3.10.2
module load scipy-stack/2022a
module load nlopt/2.7.0
module load ipykernel/2022a
source ~/projects/def-rieseber/Dryas_shared_data/Fernando_SFS/dadi_vir_env/bin/activate

cd ~/projects/def-rieseber/Dryas_shared_data/Fernando_SFS/

python3

################################################
import os
import dadi
import pickle
import matplotlib.pyplot as plt
import numpy
import scipy
import csv
import matplotlib.pyplot as pyplot
import nlopt

os.getcwd()

##############################
# Not used - SFS made in easySFS

# make the data dictionary
#dd = dadi.Misc.make_data_dict_vcf("example.vcf.gz", "popfile.txt")
dd_CA1 = dadi.Misc.make_data_dict_vcf("vcf.fsc.nomaf.miss0.1.vcf", "popsCA1.txt")
dd_CA4 = dadi.Misc.make_data_dict_vcf("vcf.fsc.nomaf.miss0.1.vcf", "popsCA4.txt")
dd_MA = dadi.Misc.make_data_dict_vcf("vcf.fsc.nomaf.miss0.1.vcf", "popsMA.txt")
dd_CO = dadi.Misc.make_data_dict_vcf("vcf.fsc.nomaf.miss0.1.vcf", "popsCO.txt")

dd_Fern = dadi.Misc.make_data_dict_vcf("vcf.fsc.nomaf.miss0.1.vcf", "popsFern.txt")


#-------------------------------------
# import the population file

pop_ids_CA1 = ["P1","P2","P3"]
pop_ids_CA4 = ["P1","P2","P3"]
pop_ids_MA = ["P1","P2","P3"]
pop_ids_CO = ["P1","P2","P3"]

pop_ids_Fern = ["CA1","CA4","MA", "CO"]

#-----------------------------
# make the site frequency spectrum

# fs = dadi.Spectrum.from_data_dict(data_dict, pop_ids, projections, mask_corners=True, polarized=True)
fs_CA1 = dadi.Spectrum.from_data_dict(dd_CA1, pop_ids_CA1, projections= [16,16,14], mask_corners=True, polarized=False)
fs_CA4 = dadi.Spectrum.from_data_dict(dd_CA4, pop_ids_CA4, projections= [18,20,14], mask_corners=True, polarized=False)
fs_MA = dadi.Spectrum.from_data_dict(dd_MA, pop_ids_MA, projections= [54,48,14], mask_corners=True, polarized=False)
fs_CO = dadi.Spectrum.from_data_dict(dd_CO, pop_ids_CO, projections= [18,18,14], mask_corners=True, polarized=False)

# download
#fs.to_file('dadi_admix_pop_frequency_spectrum.fs')

#------------------------------------------
# make site specific fs
CA1_P1 = dadi.Spectrum.from_data_dict(dd_CA1, ['P1'], projections = [16], mask_corners=True, polarized = False)
CA1_P2 = dadi.Spectrum.from_data_dict(dd_CA1, ['P2'], projections = [16], mask_corners=True, polarized = False)
CA1_P3 = dadi.Spectrum.from_data_dict(dd_CA1, ['P3'], projections = [14], mask_corners=True, polarized = False)

CA4_P1 = dadi.Spectrum.from_data_dict(dd_CA4, ['P1'], projections = [18], mask_corners=True, polarized = False)
CA4_P2 = dadi.Spectrum.from_data_dict(dd_CA4, ['P2'], projections = [20], mask_corners=True, polarized = False)
CA4_P3 = dadi.Spectrum.from_data_dict(dd_CA4, ['P3'], projections = [14], mask_corners=True, polarized = False)

MA_P1 = dadi.Spectrum.from_data_dict(dd_MA, ['P1'], projections = [54], mask_corners=True, polarized = False)
MA_P2 = dadi.Spectrum.from_data_dict(dd_MA, ['P2'], projections = [48], mask_corners=True, polarized = False)
MA_P3 = dadi.Spectrum.from_data_dict(dd_MA, ['P3'], projections = [14], mask_corners=True, polarized = False)

CO_P1 = dadi.Spectrum.from_data_dict(dd_CO, ['P1'], projections = [18], mask_corners=True, polarized = False)
CO_P2 = dadi.Spectrum.from_data_dict(dd_CO, ['P2'], projections = [18], mask_corners=True, polarized = False)
CO_P3 = dadi.Spectrum.from_data_dict(dd_CO, ['P3'], projections = [14], mask_corners=True, polarized = False)

#-----------
CA1 = dadi.Spectrum.from_data_dict(dd_Fern, ['CA1'], projections = [20], mask_corners=True, polarized = False)
CA4 = dadi.Spectrum.from_data_dict(dd_Fern, ['CA4'], projections = [20], mask_corners=True, polarized = False)
MA = dadi.Spectrum.from_data_dict(dd_Fern, ['MA'], projections = [116], mask_corners=True, polarized = False)
CO = dadi.Spectrum.from_data_dict(dd_Fern, ['CO'], projections = [20], mask_corners=True, polarized = False)


#------------------------
# multi poppulation fs

CA1_CA4 = dadi.Spectrum.from_data_dict(dd_Fern, ['CA1', 'CA4'], projections = [20,20], mask_corners=True, polarized = False)
CA1_MA = dadi.Spectrum.from_data_dict(dd_Fern, ['CA1', 'MA'], projections = [20,20], mask_corners=True, polarized = False)
CA1_CO = dadi.Spectrum.from_data_dict(dd_Fern, ['CA1', 'CO'], projections = [20,20], mask_corners=True, polarized = False)
CA4_MA = dadi.Spectrum.from_data_dict(dd_Fern, ['CA4', 'MA'], projections = [20,20], mask_corners=True, polarized = False)
CA4_CO = dadi.Spectrum.from_data_dict(dd_Fern, ['CA4','CO'], projections = [20,20], mask_corners=True, polarized = False)
MA_CO = dadi.Spectrum.from_data_dict(dd_Fern, ['MA', 'CO'], projections = [20,20], mask_corners=True, polarized = False)


#################################
# here upload data from easy SFS
# fs_5mer = dadi.Spectrum.from_file('./5pop_mer/dadi/Russia-Alaska-Europe-BCmer-BCtet.sfs')
# Russia_20 = dadi.Spectrum.from_file('./5pop_mer/dadi/Russia-20.sfs')
# Alaska_20 = dadi.Spectrum.from_file('./5pop_mer/dadi/Alaska-20.sfs')

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
CA1_P1_data = single_pop(CA1_P1)
CA1_P2_data = single_pop(CA1_P2)
CA1_P3_data = single_pop(CA1_P3)

CA1_P1_data
CA1_P2_data
CA1_P3_data

>>> CA1_P1_data
(4667.862057103559, 1406.7329490003724, 1235.3645252838023, -0.5354954105454023)
>>> CA1_P2_data
(3481.1759201926493, 1049.106594902328, 996.5645768833882, -0.22012121537611698)
>>> CA1_P3_data
(2507.420381836955, 788.4638115579809, 775.7275937770545, -0.07329788390324658)

#-------------
CA4_P1_data = single_pop(CA4_P1)
CA4_P2_data = single_pop(CA4_P2)
CA4_P3_data = single_pop(CA4_P3)

CA4_P1_data
CA4_P2_data
CA4_P3_data

>>> CA4_P1_data
(4897.555160628822, 1423.8931164419796, 1233.9687619226268, -0.5715573958677461)
>>> CA4_P2_data
(3184.2683982683875, 897.5484973528388, 966.0795169742507, 0.3203077014950567)
>>> CA4_P3_data
(2507.420381836955, 788.4638115579809, 775.7275937770545, -0.07329788390324658)

#------------
MA_P1_data = single_pop(MA_P1)
MA_P2_data = single_pop(MA_P2)
MA_P3_data = single_pop(MA_P3)

MA_P1_data
MA_P2_data
MA_P3_data

>>> MA_P1_data
(6893.795299169446, 1512.8217283833092, 1228.6005276768603, -0.6827127485357927)
>>> MA_P2_data
(5380.206206562125, 1212.3141148585541, 1130.8704751962057, -0.24742491609218167)
>>> MA_P3_data
(2507.420381836955, 788.4638115579809, 775.7275937770545, -0.07329788390324658)

#--------
CO_P1_data = single_pop(CO_P1)
CO_P2_data = single_pop(CO_P2)
CO_P3_data = single_pop(CO_P3)

CO_P1_data
CO_P2_data
CO_P3_data

>>> CO_P1_data
(4897.555160628822, 1423.8931164419796, 1233.9687619226268, -0.5715573958677461)
>>> CO_P2_data
(1587.7023923444908, 461.60143852826326, 648.5631891672115, 1.734032752325183)
>>> CO_P3_data
(2507.420381836955, 788.4638115579809, 775.7275937770545, -0.07329788390324658)
#---------------------------
CA1_data = single_pop(CA1)
CA4_data = single_pop(CA4)
MA_data = single_pop(MA)
CO_data = single_pop(CO)

#--------------------------------
# run code for multiple SFS

CA1_CA4_data = two_pop(CA1_CA4)
CA1_MA_data = two_pop(CA1_MA)
CA1_CO_data = two_pop(CA1_CO)
CA4_MA_data = two_pop(CA4_MA)
CA4_CO_data = two_pop(CA4_CO)
MA_CO_data = two_pop(MA_CO)


