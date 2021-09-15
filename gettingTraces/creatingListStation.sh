./analyseSDEU ../fullStationsList.txt test.root ~/Documents/augerData/2021/08/sd_2021_08_31_* > kk

cat kk | awk '{if($1=="#")printf( "%s \n", $5)}' > kk1

sort -n kk1 | uniq -c | awk '{print $2}' > kk2

for i in $(seq 9 9 74);do head -$i kk2 | tail -6 > ../listStations${i}.txt;done
