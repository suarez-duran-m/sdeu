#!/bin/bash

cat ../../coincident_histos/CCH_list.dat | awk '{print $1}' | sort -n | uniq -c | awk '{print $2}' > stlist.dat 

nl=$(wc -l stlist.dat | awk '{print $1}')

for i in $(seq 1 ${nl});
do
  st_i=$(head -${i} stlist.dat | tail -1)
  cat ../../coincident_histos/CCH_list.dat | awk -v st=${st_i} '{if($1==st)print $i}' > evtsForSt.dat
  for pmt in $(seq 1 3);
  do
    ./readTraces evtsForSt.dat ${pmt} ../reconstructing/forIoana/ADST2022_02_09_histos.root ../reconstructing/forIoana/ADST2022_02_10_histos.root 
  done
done
