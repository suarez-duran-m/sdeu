* How to run analyseSDEUbins to create baseline100bins.root file:
make clean
make
./analyseSDEUbins ../../select/listStations.txt pmtId ../dec2020Traces.root ../jan2021Traces.root

Where pmtId regards to the PMT you want to run:
* 1 For PMT1
* 2 For PMT2
* 3 For PMT3
* 4 For SPMT
* 5 For PMTSSD


* To plot baseline100bins.root:
  - Open a new terminal
  - root -L plottingBl100bins.C 
