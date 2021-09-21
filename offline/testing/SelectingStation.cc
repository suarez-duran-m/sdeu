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
      CentralConfig::GetInstance()->GetTopBranch("ReadSdfiles");
    topB.GetChild("OutNameRootStId").GetData(fChoseStId);

    return eSuccess;
}


VModule::ResultFlag SelectingStation::Run(evt::Event& event) {
  INFO("SelectingStation::Run()");

  const sevt::SEvent& sEvent = event.GetSEvent();

  stOk = false;
  for (sevt::SEvent::ConstStationIterator sIt = sEvent.StationsBegin(); sIt != sEvent.StationsEnd(); ++sIt)
    if ( fChoseStId == sIt->GetId() )
      stOk = true;

  if ( !stOk )
    return eContinueLoop;

  return eSuccess;
}


VModule::ResultFlag SelectingStation::Finish() {  
  INFO("SelectingStation::Finish()");
  return eSuccess;
}
