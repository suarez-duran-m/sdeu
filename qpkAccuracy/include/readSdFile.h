#ifndef READSDFILE_H
#define READSDFILE_H

#include "IoAuger.h"
#include <vector>
#include <string>
#include <iostream>

#include "TTree.h"
#include "TGraphErrors.h"

using namespace std;

class readSdFile {
  public:
    readSdFile();
    readSdFile(int nfiles, char * fileNames[], int binLR);
    virtual ~readSdFile(){}

    int binsLeftRigh;

    void GetAndWriteChHisto(int StId, char *date);
    const vector < vector < double > > GetQpk(int StId, char *date);

  private:
    int totalEvents;
    AugerIoSd *sdFile;

    int IdAndTime[2];
    double fitVar[10];
    double tmpTestQpk;
    TGraphErrors *chPmt1;
    TGraphErrors *chPmt2;
    TGraphErrors *chPmt3;

    void SetTree(TTree *pmtTree);
    TGraphErrors *FillingTree(TH1F *hist, int Offset, int binsLeftRigh);
};
#endif
