#!/bin/bash

list=$1

nl=$(wc -l $HOME/2021/sdeu/listStations${list}.txt | awk '{print $1}')

for st in $(seq 1 ${nl});
do
  stId=$(head -${st} $HOME/2021/sdeu/listStations${list}.txt | tail -1)
  for pmt in $(seq 1 3);
  do
    cp ReadSdfiles.xml.in ReadSdfiles.xml
    sed -i "s/XXX/${stId}/g" ReadSdfiles.xml
    sed -i "s/CCC/${pmt}/g" ReadSdfiles.xml
    echo ""
    echo "Doing for ${stId} and PMT${pmt}"
    echo ""
    time ./userAugerOffline 1>kk 2> >(gzip > stderr.gz)
  done
done
