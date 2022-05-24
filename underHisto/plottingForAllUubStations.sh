#!/bin/bash

nl=$(wc -l ../fullUubStationsListVert.txt | awk '{print $1}')

for st in $(seq 1 ${nl});
do
  stId=$(head -${st} ../fullUubStationsListVert.txt | tail -1)
  root -l -x "readingChMonthsFitting.C(${stId})"
done
