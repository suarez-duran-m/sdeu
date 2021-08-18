#/bin/bash

st=(0863 1222 1219 1211 1740 1743 1221 1223 1217 1747 1741 1745 1818 1851 1729 1735 1746 1819 1791)

make clean
make

for i in $(seq 0 18);do
  echo ${st[$i]} > tmp.txt
  ./analyseSDEUbl tmp.txt bl${i}.root ../???202?Traces.root
done

