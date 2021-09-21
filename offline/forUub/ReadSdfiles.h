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

class ReadSdfiles : public fwk::VModule {
  public:
    ReadSdfiles() { }
    virtual ~ReadSdfiles() { }

    fwk::VModule::ResultFlag Init();
    fwk::VModule::ResultFlag Run(evt::Event& event);
    fwk::VModule::ResultFlag Finish();

    //Functions

    private:
    REGISTER_MODULE("ReadSdfiles",ReadSdfiles);

    TFile *hfile; 
    TTree *peak;
    TTree *charge;
    
    unsigned int pkEvtId;
    unsigned int pkGpstime;
    double pk;
    double pkChi2;
    unsigned int pkNdof;
    double pkLow;
    double pkHigh;
    double pkP0;
    double pkP1;
    double pkP2;

    unsigned int chEvtId;
    unsigned int chGpstime;
    double ch;
    double chChi2;
    unsigned int chNdof;
    double chLow;
    double chHigh;
    double chP0;
    double chP1;
    double chP2;

    // XML config
    int fOutNameRootYear;
    int fOutNameRootMonth;
    int fOutNameRootPmt;
    int fOutNameRootStId;

    string monthUub[8] = {"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug"};
};
