#ifndef READADSTFILE_H
#define READADSTFILE_H

#include <RecEvent.h>
#include <RecEventFile.h>
#include <Traces.h>

#include <vector>

using namespace std;

class readAdstFile {
  public:
    readAdstFile();

    void setAdstFileName(const char* const file);
    void checkNumberEvts();
    void readSignals(int st);

    const vector < vector < double > > getSglPmt();
    const vector < vector < double > > getMuonTauPmt();
    const vector < double > getTotSgl(); 
    const vector < double > getTotSglErr();
    const vector < double > getSglGpsSec();

  private:
    // from ADST/RecEvent/src/
    const char* fileName;
    RecEvent *theRecEvent;
    void readEvents(int st);

    vector < vector < double > > sglPmt;
    vector < vector < double > > muonTauPmt;
    vector < double > totSgl;
    vector < double > totSglErr;
    vector < double > sglGpsSec;
};
#endif
