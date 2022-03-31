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

class GetFitOutputsMSD : public fwk::VModule {
  public:
    GetFitOutputsMSD() { }
    virtual ~GetFitOutputsMSD() { }

    fwk::VModule::ResultFlag Init();
    fwk::VModule::ResultFlag Run(evt::Event& event);
    fwk::VModule::ResultFlag Finish();

    //Functions

    private:
    REGISTER_MODULE("GetFitOutputsMSD",GetFitOutputsMSD);

    TFile *hfile; 
    TTree *charge;

    unsigned int chEvtId;
    unsigned int chGpstime;
    double qpk;

    // XML config
    string fOutPutFileName;
    int fSelectPmt;
    int fSelectStId;
};
