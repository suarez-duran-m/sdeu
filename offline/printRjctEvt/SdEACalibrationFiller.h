/** 
   \file
   declaration of SdEACalibrationFiller

   \author Alvaro Taboada Nunez
   \date 05 Mar 2018
  */

#ifndef _SdEACalibrationFillerKG_SdEACalibrationFiller_h_
#define _SdEACalibrationFillerKG_SdEACalibrationFiller_h_

#include <fwk/VModule.h>
#include <sevt/Station.h>
#include <utl/Histogram.h>

namespace evt {
  class Event;
}

namespace SdEACalibrationFillerKG {

  class SdEACalibrationFiller : public fwk::VModule {

  public:
    fwk::VModule::ResultFlag Init();
    fwk::VModule::ResultFlag Run(evt::Event& event);
    fwk::VModule::ResultFlag Finish();

  private:
    void FillCalibrationInfo(sevt::Station& station);
    void FillCalibrationInfo_SmallPMT(sevt::Station& station);
    
    std::vector<double> fRawChargeUUB;
    std::vector<double> fRawPeakUUB;
    std::vector<double> fHgBaselineUUB;
    std::vector<double> fLgBaselineUUB;
    std::vector<double> fDARatioUUB;

    std::vector<double> fRawChargePPA;
    std::vector<double> fRawPeakPPA;

    double rawPeak;
    double rawCharge;
    double hgBaseline;
    double lgBaseline;
    double DARatio;
    int vmax;

    double fBeta;
    double fBetaErr;
    double fBetaChi2;
    double fBetaCorrectionFact;

    REGISTER_MODULE("SdEACalibrationFillerKG", SdEACalibrationFiller);
  };

}


#endif
