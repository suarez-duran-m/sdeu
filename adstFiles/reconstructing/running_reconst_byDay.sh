#!/bin/bash

for day in $(seq 1 31);
do
  cp baseEventFileReader.xml EventFileReader.xml
  cp baseBootstrap.xml bootstrap.xml
  if [ ${day} -le 9 ];
  then   
    sed -i "s/DAY/0${day}/" EventFileReader.xml
    sed -i "s/DAY/0${day}/" bootstrap.xml
  fi
  if [ ${day} -gt 9 ];
  then
    sed -i "s/DAY/${day}/" EventFileReader.xml
    sed -i "s/DAY/${day}/" bootstrap.xml
  fi
  ./userAugerOffline
done
