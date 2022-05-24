#!/bin/bash

list=$1

nl=$(wc -l $HOME/2021/sdeu/listStations${list}.txt | awk '{print $1}')

for st in $(seq 1 ${nl});
do
  stId=$(head -${st} $HOME/2021/sdeu/listStations${list}.txt | tail -1)
  for pmt in $(seq 1 3);
  do
    cp bootstrap_base.xml bootstrap${list}${st}.xml
    sed -i "s/XXX/${stId}/g" bootstrap${list}${st}.xml 
    sed -i "s/CCC/${pmt}/g" bootstrap${list}${st}.xml
    echo ""
    echo "Doing for ${stId} and PMT${pmt}"
    echo ""
    time ./userAugerOffline -b bootstrap${list}${st}.xml 1>kk${list}${st} 2> >(gzip > stderr${list}${st}.gz)
  done
done
