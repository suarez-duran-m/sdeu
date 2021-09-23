#include "SelectingStation.h"

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

VModule::ResultFlag SelectingStation::Init() {
    INFO("SelectingStation::Init()");

    Branch topB = 
      CentralConfig::GetInstance()->GetTopBranch("GetFitOutputs");
    topB.GetChild("SelectStId").GetData(fChoseStId);

    return eSuccess;
}


VModule::ResultFlag SelectingStation::Run(evt::Event& event) {
  INFO("SelectingStation::Run()");

  const sevt::SEvent& sEvent = event.GetSEvent();

  for (sevt::SEvent::ConstStationIterator sIt = sEvent.StationsBegin(); sIt != sEvent.StationsEnd(); ++sIt)
    if ( fChoseStId == sIt->GetId() )
      return eSuccess;

  return eContinueLoop;  
}


VModule::ResultFlag SelectingStation::Finish() {  
  INFO("SelectingStation::Finish()");
  return eSuccess;
}
