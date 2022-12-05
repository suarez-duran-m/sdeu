#include "msg_unknown_read.h"
#include <cstring>
#include <arpa/inet.h>

#include <iostream>
#include <fstream>
#include "rawCoincHistoData.h"

#include <IoAuger.h>

// Initializing class
rawCoincHistoData::rawCoincHistoData() { 
  cout << "MSD" << endl;
}

struct muon_extra_histo {
  /* charge and peak for the coincidence (WCD+SSD) histograms.
   *  It will include only WCD PMTs (index from 0 to 2)
   *
   */
  uint16_t has_histo;
  uint16_t ssd_th; /* the applied threshold */

  uint16_t Charge[3][600];
  uint16_t Peak[3][150];
};


void rawCoincHistoData::fetch_coinciHistos( char *buff, int st, int gps ) {
  struct muon_extra_histo h;
  memcpy(&h, buff, sizeof(struct muon_extra_histo));
  h.has_histo = ntohs(h.has_histo);
  h.ssd_th = ntohs(h.ssd_th);
  vector < vector < int > > PmtCharge;
  vector < vector < int > > PmtHeight;
  // Storing data into private variables
  // Fetching cCharge and cHeight
  //
  for(int pmt=0; pmt<3; pmt++) {
    vector < int > charge;
    for(int bin=0; bin<600; bin++)
      charge.push_back( ntohs(h.Charge[pmt][bin]) );

    vector < int > height;
    for(int bin=0; bin<150; bin++)
      height.push_back( ntohs(h.Peak[pmt][bin]) );

    PmtCharge.push_back( charge );
    PmtHeight.push_back( height );
  }
  // 
  // Writing histo into ascii file for cross-check  
  /*
  for ( int pmt_i=0; pmt_i<3; pmt_i++ ) {
    fstream outPutHistoCrossCheck;
    TString outputName = Form("%d_%d_%d", gps+2, st, pmt_i);
    outPutHistoCrossCheck.open("histos_crossCheck2/"+outputName+".dat", ios_base::out);
    for ( auto &i : PmtCharge[pmt_i] )
      outPutHistoCrossCheck << i << endl;
    outPutHistoCrossCheck.close();
  }
  */
  cQhisto.push_back( PmtCharge );
  cHeight.push_back( PmtHeight ); 
  PmtCharge.clear();
  PmtCharge.erase(PmtCharge.begin(), PmtCharge.end());
  PmtHeight.clear();
  PmtHeight.erase(PmtHeight.begin(), PmtHeight.end());
}


void rawCoincHistoData::readData( const char *file ) {
  cout << "MSD rawCoincHistoData::readData" << endl;
  char buff[1048576]; 
  int flag;
  const int gps2utc = 315964782;
  struct cl_msg_unknown_pack_h h;
  struct cl_msg_unknown_pack_H H;
  msg_unknown_read ff(file, "!!_T3_!!", 40960);
  
  flag = 0;
  while( flag==0 ){
    flag = ff.search_preamble();
    if(flag != 0){
      break;
    }
    flag = ff.get_header(&H);
    if(flag != 0){
      printf("header error\n");
      break;
    }
    int tmpCnt = 0;
    int entries = 0;
    while(ff.get_pack(&h, buff, 1048570)==0)
      if( h.type==3 && h.version==1 ) {
        utcEvtWidthChisto.push_back( H.timestamp_sec+gps2utc );
        stWidthChisto.push_back( H.LsId );
        fetch_coinciHistos( buff, H.LsId, H.timestamp_sec );
        tmpCnt++;
      }
    ff.skip();
  }
}


void rawCoincHistoData::SetClear() {
  utcEvtWidthChisto.clear();
  stWidthChisto.clear();
  cQhisto.clear();
  cHeight.clear();
}