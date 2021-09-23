#include "GetFitOutputs.h"

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

VModule::ResultFlag GetFitOutputs::Init() {
    INFO("GetFitOutputs::Init()");
 
    Branch topB = 
      CentralConfig::GetInstance()->GetTopBranch("GetFitOutputs");

    topB.GetChild("OutPutFileName").GetData(fOutPutFileName);
    topB.GetChild("SelectStId").GetData(fSelectStId);
    topB.GetChild("SelectPMT").GetData(fSelectPMT);
    
    hfile = new TFile(fOutPutFileName.c_str(), "RECREATE","");

    pkEvtId = 0;
    pkGpstime = 0;
    pk = 0.;
    pkChi2 = 0.;
    pkNdof = 0;
    pkLow = 0.;
    pkHigh = 0.;
    pkP0 = 0.;
    pkP1 = 0.;
    pkP2 = 0.;

    chEvtId = 0;
    chGpstime = 0;
    ch = 0.;
    chChi2 = 0.;
    chNdof = 0;
    chLow = 0.;
    chHigh = 0.;
    chP0 = 0.;
    chP1 = 0.;
    chP2 = 0.;

    peak = new TTree("peak","");
    peak->Branch("evtId",&pkEvtId,"pkEvtId/I");
    peak->Branch("GpsTime",&pkGpstime,"pkGpstime/I");
    peak->Branch("peakVal",&pk,"pk/D");
    peak->Branch("peakChi2",&pkChi2,"pkChi2/D");
    peak->Branch("peakNdofl",&pkNdof,"pkNdof/I");
    peak->Branch("peakLow",&pkLow,"pkLow/D");
    peak->Branch("peakHigh",&pkHigh,"pkHigh/D");
    peak->Branch("peakP0",&pkP0,"pkP0/D");
    peak->Branch("peakP1",&pkP1,"pkP1/D");
    peak->Branch("peakP2",&pkP2,"pkP2/D");

    charge = new TTree("charge","");
    charge->Branch("evtId",&chEvtId,"chEvtId/I");
    charge->Branch("GpsTime",&chGpstime,"chGpstime/I");
    charge->Branch("chargeVal",&ch,"ch/D");
    charge->Branch("chargeChi2",&chChi2,"chChi2/D");
    charge->Branch("chargeNdofl",&chNdof,"chNdof/I");
    charge->Branch("chargeLow",&chLow,"chLow/D");
    charge->Branch("chargeHigh",&chHigh,"chHigh/D");
    charge->Branch("chargeP0",&chP0,"chP0/D");
    charge->Branch("chargeP1",&chP1,"chP1/D");
    charge->Branch("chargeP2",&chP2,"chP2/D");

    return eSuccess;
}


VModule::ResultFlag GetFitOutputs::Run(evt::Event& event) {

  INFO("GetFitOutputs::Run()");
  if ( !event.HasSEvent() )
    return eSuccess;

  const sevt::SEvent& sEvent = event.GetSEvent();

  for (sevt::SEvent::ConstStationIterator sIt = sEvent.StationsBegin(); sIt != sEvent.StationsEnd(); ++sIt) {
    const sdet::Station& dStation = det::Detector::GetInstance().GetSDetector().GetStation(*sIt);
    if ( !dStation.IsUUB() )
      continue;
    
    if ( sIt->GetId() == fSelectStId ) {
      sevt::PMT& pmt = event.GetSEvent().GetStation(sIt->GetId()).GetPMT( fSelectPMT );
      if ( !pmt.HasRecData() )
        pmt.MakeRecData();
      
      const sevt::PMTRecData& recpmtDat = pmt.GetRecData();
      const sevt::StationGPSData& gpstime = event.GetSEvent().GetStation(sIt->GetId()).GetGPSData();

      pkEvtId = sEvent.GetHeader().GetId();
      pkGpstime = gpstime.GetSecond();
      
      if ( recpmtDat.IsVEMPeakFromHistogram() ) {
        pk = recpmtDat.GetVEMPeak();
        pkChi2 = recpmtDat.GetVEMPeakChi2();
        pkNdof = recpmtDat.GetVEMPeakNdof();
        pkLow = recpmtDat.GetVEMPeakLow();
        pkHigh = recpmtDat.GetVEMPeakHigh();
        pkP0 = recpmtDat.GetVEMPeakP0();
        pkP1 = recpmtDat.GetVEMPeakP1();
        pkP2 = recpmtDat.GetVEMPeakP2();
      }
      else {
        pk = 0.;
        pkChi2 = 0.;
        pkNdof = 0;
        pkLow = 0.;
        pkHigh = 0.;
        pkP0 = 0.;
        pkP1 = 0.;
        pkP2 = 0.;
      }

      chEvtId = sEvent.GetHeader().GetId();
      chGpstime = gpstime.GetSecond();
 
      //if ( recpmtDat.IsVEMChargeFromHistogram() ) {
        ch = recpmtDat.GetVEMCharge();
        chChi2 = recpmtDat.GetVEMChargeChi2();
        chNdof = recpmtDat.GetVEMChargeNdof();
        chLow = recpmtDat.GetVEMChargeLow();
        chHigh= recpmtDat.GetVEMChargeHigh();
        chP0 = recpmtDat.GetVEMChargeP0();
        chP1 = recpmtDat.GetVEMChargeP1();
        chP2 = recpmtDat.GetVEMChargeP2();
      /*
      }      
      else {
        ch = 0.;
        chChi2 = 0.;
        chNdof = 0;
        chLow = 0.;
        chHigh = 0.;
        chP0 = 0.;
        chP1 = 0.;
        chP2 = 0.;
      }
      */
      peak->Fill();
      charge->Fill();
    } 
  }
  
  return eSuccess;  
}


VModule::ResultFlag GetFitOutputs::Finish() {  
  INFO("GetFitOutputs::Finish()");
  hfile->Write();
  hfile->Close();
  return eSuccess;
}
