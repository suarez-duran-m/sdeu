#include <fwk/VModule.h>

#include <utl/ErrorLogger.h>

#include <evt/Event.h>

#include <string>
#include <time.h>

#include <sevt/SEvent.h>
#include <sevt/Header.h>
#include <sevt/Station.h>
#include <sevt/StationRecData.h>
#include <sevt/EventTrigger.h>
#include <sevt/StationTriggerData.h>
#include <sevt/StationGPSData.h>

#include <det/Detector.h>

#include <sdet/SDetector.h>
#include <sdet/Station.h>
#include <sdet/StationTriggerAlgorithm.h>
#include <sdet/PMTConstants.h>

#include <TFile.h>
#include <TTree.h>

using namespace std;

class SelectingStation : public fwk::VModule {
  public:
    SelectingStation() { }
    virtual ~SelectingStation() { }

    fwk::VModule::ResultFlag Init();
    fwk::VModule::ResultFlag Run(evt::Event& event);
    fwk::VModule::ResultFlag Finish();

    //Functions

    private:
    REGISTER_MODULE("SelectingStation",SelectingStation);
    bool stOk;

    // XML config
    int fGETstId = 0;
};
