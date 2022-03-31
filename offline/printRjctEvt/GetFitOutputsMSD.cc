#include "GetFitOutputsMSD.h" 
#include <fwk/CentralConfig.h>

#include <det/Detector.h>

#include <sdet/SDetector.h>
#include <sdet/Station.h>
#include <sdet/PMT.h>
#include <sdet/PMTConstants.h>

#include <sevt/PMTRecData.h>
#include <sevt/Station.h>
#include <sevt/StationGPSData.h>

using namespace std;
using namespace utl;
using namespace fwk;

VModule::ResultFlag GetFitOutputsMSD::Init() {
    INFO("GetFitOutputsMSD::Init()");
 
    Branch topB = 
      CentralConfig::GetInstance()->GetTopBranch("GetFitOutputsMSD");

    topB.GetChild("OutPutFileName").GetData(fOutPutFileName);
    topB.GetChild("SelectStId").GetData(fSelectStId);
    topB.GetChild("SelectPmt").GetData(fSelectPmt);
    
    hfile = new TFile(fOutPutFileName.c_str(), "RECREATE","");

    chEvtId = 0;
    chGpstime = 0;
    qpk = 0.;

    charge = new TTree("charge","");
    charge->Branch("evtId",&chEvtId,"chEvtId/I");
    charge->Branch("GpsTime",&chGpstime,"chGpstime/I");
    charge->Branch("qpkVal",&qpk,"qpk/D");
    return eSuccess;
}


VModule::ResultFlag GetFitOutputsMSD::Run(evt::Event& event) {

  INFO("GetFitOutputsMSD::Run()");
  if ( !event.HasSEvent() )
    return eSuccess;

  const sevt::SEvent& sEvent = event.GetSEvent();

  bool stFound = false;

  for (sevt::SEvent::ConstStationIterator sIt = sEvent.StationsBegin(); sIt != sEvent.StationsEnd(); ++sIt) {
    //const sdet::Station& dStation = det::Detector::GetInstance().GetSDetector().GetStation(*sIt);
    if (sIt->GetId() == fSelectStId)
     stFound = true; 
    
    //if ( !dStation.IsUUB() )
      //continue;
    /*
    if ( sIt->GetId() == fSelectStId ) {
      sevt::PMT& pmt = event.GetSEvent().GetStation(sIt->GetId()).GetPMT( fSelectPmt );
      if ( !pmt.HasRecData() )
        pmt.MakeRecData();
      
      const sevt::PMTRecData& recpmtDat = pmt.GetRecData();
      const sevt::StationGPSData& gpstime = event.GetSEvent().GetStation(sIt->GetId()).GetGPSData();

      chEvtId = sEvent.GetHeader().GetId();
      chGpstime = gpstime.GetSecond();
 
      if ( recpmtDat.IsVEMChargeFromHistogram() )
        qpk = recpmtDat.GetVEMCharge(); 
      else
        qpk = 0.;
      charge->Fill();
    }
    */
  }
  stFound = true;
  if ( stFound )
    return eSuccess;
  else
    return eContinueLoop;
}


VModule::ResultFlag GetFitOutputsMSD::Finish() {  
  INFO("GetFitOutputsMSD::Finish()");
  hfile->Write();
  hfile->Close();
  return eSuccess;
}
