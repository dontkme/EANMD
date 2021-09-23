#!/bin/bash
for i in $(ls -d Brain*)
do
cd $i
perl ../EANMD_filterPSI.pl -o $i.PSIfilter.out -i 0.1 -d 20 -m 2 -f 0.05 -c1 6 -c2 2 -mf 0.05 SE.MATS.JCEC.txt
echo "Done $i.PSIfilter"
cd ..
done
