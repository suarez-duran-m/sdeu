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
    double getMean ( const vector<int> &arr, 
				const unsigned int nb, const bool fok );
    double getRms ( const vector<int> &arr, 
				const double meanarr, const unsigned int nb, const bool fok );
};

#endif
