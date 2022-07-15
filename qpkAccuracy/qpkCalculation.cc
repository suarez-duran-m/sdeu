#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "include/readSdFile.h"
//#include "include/plotDiffDist.h"

using namespace std;

/********************/
/* Global variables */
/********************/


int main ( int argc, char *argv[]) {
  if ( argc < 4 ) {
    cout << endl << "=========================" << endl << endl
      << "Usage: " << argv[0] << " file_StIds Date sd_files" << endl;
    cout << endl << "file_StIds: ascii with station Ids" << endl
      << "Date: MonthYear, ex. Sep2021" << endl
      << "sd_files: sd*.root" << endl << endl;
    exit(0);
  } 

  const char *filename = argv[1];
  char *date = argv[2];
  ifstream fileWithIds(filename, ios::in);
  if (!fileWithIds.is_open()) {
    cerr << "Could not open file: " << argv[1] << endl; 
    exit(0);
  }

  vector < double > st2read;
  unsigned int id = 0;

  while ( fileWithIds.good() ) {
    id = 0;
    fileWithIds >> id;
    st2read.push_back( id );
  }

  int binsLeftRigh = 35;

  readSdFile sdFile(argc-3, argv+3, binsLeftRigh);
  sdFile.GetAndWriteChHisto(st2read[0], date);
  //readSdFile sdFile;
  //vector < vector < double > > qpks = sdFile.GetQpk(st2read[0], date);  


  return 0;
}
