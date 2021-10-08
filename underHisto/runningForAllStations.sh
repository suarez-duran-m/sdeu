#!/bin/bash

listfile=$1
nl=$(wc -l ../listStations${listfile}.txt | awk '{print $1}')

for st in $(seq 1 ${nl});
do
  head -${st} ../listStations${listfile}.txt | tail -1 > list${listfile}${st}.txt
  for pmt in $(seq 1 3);
  do
      ./analyseSDEUcalib list${listfile}${st}.txt ${pmt} 35 Sep ../gettingTraces/sd*Sep2021Traces.root > kk${listfile}${st}
  done
done
