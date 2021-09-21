#!/bin/bash
 
nl=$(wc -l $HOME/2021/sdeu/fullUubStationsListVert.txt | awk '{print $1}')

for st in $(seq 1 ${nl});
do
  stId=$(head -${st} $HOME/2021/sdeu/fullUubStationsListVert.txt | tail -1)
  for pmt in $(seq 1 3);
  do
    root -l "readingChOffCdas.C(${stId}, ${pmt})"
  done
done
