#ifndef RAWCOINCIHISTODATA_H 
#define RAWCOINCIHISTODATA_H 

#include <iostream>
#include <vector>

using namespace std;

class rawCoinciHistoData {
  public:
		rawCoinciHistoData();
    ~rawCoinciHistoData(){}

    void readData( const char *file );
    vector < int > getUtcEvtWidthChisto() {
      return utcEvtWidthChisto;
    }
    vector < int > getStWidthChisto() {
      return stWidthChisto;
    }
    vector < vector < vector < int > > > getCQhisto() {
      return cQhisto;
    }
    vector < vector < vector < int > > > getCheight() {
      return cHeight;
    }


	private:
    void fetch_coinciHistos( char *buff, int st, int gps );
    // evtWithCQ[utc]
    vector < int > utcEvtWidthChisto;
    // stWidthCQ[utc-st]
    vector < int > stWidthChisto;
    // cQhisto[utc-st][pmt][bin-cnts]
    vector < vector < vector < int > > > cQhisto;
    // cHeight[utc-st][pmt][bin-cnts]
    vector < vector < vector < int > > > cHeight;
};

#endif
