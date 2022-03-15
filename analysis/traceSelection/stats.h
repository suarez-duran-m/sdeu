#ifndef STATS_H
#define STATS_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <math.h>

using namespace std;

class stats {
  public: 
    double getMean ( vector<int> *arr, unsigned int nb, bool fok );
    double getRms ( vector<int> *arr, double meanarr, unsigned int nb, bool fok );
};

#endif
