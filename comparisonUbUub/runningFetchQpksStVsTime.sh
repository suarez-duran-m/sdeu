#!/bin/bash

#nl=$(wc -l ../fullUubStationsListVert.txt | awk '{print $1}')
nl=$(wc -l listStatOk.dat | awk '{print $1}')

for st in $(seq 1 ${nl});
do
  #stId=$(head -${st} ../fullUubStationsListVert.txt | tail -1)
  stId=$(head -${st} listStatOk.dat | tail -1)
  root -l "fetchingQpksPerStationVsTime.C(true, ${stId}, false)"
done
