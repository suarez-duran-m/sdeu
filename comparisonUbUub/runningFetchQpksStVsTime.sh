#!/bin/bash

listStations="../fullUubStationsListVert.txt"
#listStations="listStatOk.dat"


nl=$(wc -l ${listStations} | awk '{print $1}')

for st in $(seq 1 ${nl});
do
  stId=$(head -${st} ${listStations} | tail -1)
  root -l "fetchingQpksPerStationVsTime.C(false, ${stId}, false)"
done
