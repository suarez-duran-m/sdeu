#!/bin/bash

readList=$1
make clean
make

nl=$(wc -l ${readList} | awk '{print $1}')

for line in $(seq 1 ${nl});
do
  stid=$(head -${line} ${readList} | tail -1)
  echo ${stid} > ${stid}.dat
  for pmt in $(seq 1 4);
  do
    ./analyseSDEUcalib ${stid}.dat ${pmt} 35 Dec ../Raid_CDAS_test_2022_02_02/data/Sd/2022/02/sd_2022_02_02_13h47.root > kk
  done
  rm ${stid}.dat
done 
