#!/bin/bash

listfile=$1
nl=$(wc -l ../listStations${listfile}.txt | awk '{print $1}')

for st in $(seq 1 ${nl});
do
  head -${st} ../listStations${listfile}.txt | tail -1 > list${listfile}${st}.txt
  ./readuubmc list${listfile}${st}.txt 1 ~/Documents/augerData/monit/2021/08/mc_2021_08_* ~/Documents/augerData/monit/2021/09/mc_2021_09_* ~/Documents/augerData/monit/2021/10/mc_2021_10_*
done
