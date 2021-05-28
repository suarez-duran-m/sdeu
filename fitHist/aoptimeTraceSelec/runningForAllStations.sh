#!/bin/bash

month=$1

list=(863 1222 1219 1211 1740 1743 1221 1223 1217 1747 1741 1745 1818 1851 1729 1735 1746 1819 1791)

echo 

for st in $(seq 0 0); #18);
do
  echo ${list[${st}]} > list${list[${st}]}.dat
  for pmt in $(seq 1 3);
  do
    ./analyseSDEUcalib list${list[${st}]}.dat ${pmt} ${month} ../../gettingTraces/*${month}*.root > kk${list[${st}]}${month} 

    #python3 root2jsonAoPchpk.py uub ${pmt} ${list[${st}]} ${month}
  done
done
