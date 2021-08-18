#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "stats.h"

double stats::getMean( vector<int> *arr, unsigned int nb, bool fok ) {
	double mean = 0.;
  int lb = arr->size() - 1;
    for ( unsigned int i=0; i<nb; i++ ){
      if ( fok )
        mean += (*arr)[i];
      else
        mean += (*arr)[lb-i];
    }
  return mean/nb;
}

double stats::getRms( vector<int> *arr, double meanarr, unsigned int nb, bool fok ) {
  double rms = 0.;
  int lb = arr->size() - 1;
  for ( unsigned int i=0; i<nb; i++ ){
    if ( fok )
      rms += ( (*arr)[i] - meanarr )*( (*arr)[i] - meanarr );
    else
      rms += ( (*arr)[lb-i] - meanarr )*( (*arr)[lb-i] - meanarr );
  }
  return sqrt(rms/nb);
}
