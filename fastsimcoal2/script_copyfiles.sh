#! /bin/bash
for i in {1..50}
do
mkdir 5pops_mertensiana_$i
cp 5pops_mertensiana_MSFS.obs 5pops_mertensiana_$i
cp 5pops_mertensiana.tpl 5pops_mertensiana_$i
cp 5pops_mertensiana.est 5pops_mertensiana_$i
done

for i in {1..50}
do cd 5pops_mertensiana_$i
nice ../../../../fastsimcoal/fsc2709 -t *.tpl -n 100000 -e *.est -0 -m -M -L 40 -c 30 --multiSFS -q
cd ..
done

