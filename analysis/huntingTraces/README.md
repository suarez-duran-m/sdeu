* How to run analyseSDEUhunting:
make clean
make
./analyseSDEUhunting ../../select/listStations.txt pmtId ../dec2020Traces.root ../jan2021Traces.root

# The Aim
Identify mu-Traces and c-Traces. Here, mu-Traces are traces with a muon type signal in the last 100 bins; whereas c-Traces are traces with low frequency signal in th las 100 bins (see slide 2 on slides_results.pdf).

* For mu-Traces: mean cummulative
* For c-Traces: bins below Mean-2sigma.


