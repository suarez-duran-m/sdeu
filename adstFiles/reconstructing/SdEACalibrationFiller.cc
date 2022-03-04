/**
   \file 
   Implementation of SdEACalibrationFiller

   \author Alvaro Taboada Nunez

   \date 05 Mar 2018
*/

#include "SdEACalibrationFiller.h"

#include <fwk/CentralConfig.h>
#include <det/VManager.h>
#include <utl/Branch.h>

#include <evt/Event.h>
#include <sevt/SEvent.h>
#include <sevt/Station.h>
#include <sevt/SmallPMT.h>
#include <sevt/SmallPMTCalibData.h>
#include <sdet/PMTConstants.h>
#include <sdet/SDetector.h>

#include <utl/ErrorLogger.h>

using namespace std;

using namespace fwk;
using namespace det;
using namespace utl;
using namespace evt;
using namespace sevt;


namespace SdEACalibrationFillerKG {

  VModule::ResultFlag
  SdEACalibrationFiller::Init()
  {
    const Branch topB = CentralConfig::GetInstance()->GetTopBranch("SdEACalibrationFiller");
    
    Branch onlineValuesUUB = topB.GetChild("onlineValuesUUB");
    onlineValuesUUB.GetChild("rawCharge").GetData(fRawChargeUUB);
    onlineValuesUUB.GetChild("rawPeak").GetData(fRawPeakUUB);
    onlineValuesUUB.GetChild("hgBaseline").GetData(fHgBaselineUUB);
    onlineValuesUUB.GetChild("lgBaseline").GetData(fLgBaselineUUB);
    onlineValuesUUB.GetChild("DARatio").GetData(fDARatioUUB);

    Branch onlineValuesKitBox = topB.GetChild("onlineValuesPPA");
    onlineValuesKitBox.GetChild("rawCharge").GetData(fRawChargePPA);
    onlineValuesKitBox.GetChild("rawPeak").GetData(fRawPeakPPA);

    Branch offlineValuesUUB = topB.GetChild("offlineValuesUUB");
    offlineValuesUUB.GetChild("beta").GetData(fBeta);
    offlineValuesUUB.GetChild("betaErr").GetData(fBetaErr);
    offlineValuesUUB.GetChild("betaChi2").GetData(fBetaChi2);
    offlineValuesUUB.GetChild("betaCorrectionFact").GetData(fBetaCorrectionFact);
    
    return eSuccess;
  }


  VModule::ResultFlag
  SdEACalibrationFiller::Run(Event& event)
  {
    INFO("Filling PMT Calibration values for Upgraded Stations according to information from EA data analysis!");

    if (!event.HasSEvent())
      return eSuccess;

    SEvent& sEvent = event.GetSEvent();

    for (SEvent::StationIterator sIt = sEvent.StationsBegin(), end = sEvent.StationsEnd(); sIt != end; ++sIt) {

      const sdet::Station& dStation =
        det::Detector::GetInstance().GetSDetector().GetStation(*sIt);

      const bool isUub = dStation.IsUUB();
      const bool hasSmallPMT = dStation.HasSmallPMT();
      const bool hasScintillator = dStation.HasScintillator(); // for stations that have SSD conected to UB through kit-box
                                                               // up to now (Mar 2018) UUB stations are EA stations
      if (isUub || hasScintillator){
        FillCalibrationInfo(*sIt);
        if(hasSmallPMT)
          FillCalibrationInfo_SmallPMT(*sIt);
      }
      else
        continue;

    }

    return eSuccess;
  }

  void
  SdEACalibrationFiller::FillCalibrationInfo(Station& station) {

    /*
    Values needed for Pre Production Array (PPA) stations:
    Estimate MIP Charge (for SSD)
    Estimate of MIP Peak (in CDAS is estimated as Charge / 1.8)
    SetIsTubeOk & SetIsLowGainOk for SSD

    Values needed for Engineering Array stations (with UUB):
    Estimate of VEM and MIP Charge
    Estimate of VEM and MIP Peak
    Estimate of HG and LG baseline
    Set values of DA ratio
    SetIsLowGainOk
    SetIsTubeOk
    */

    const sdet::Station& dStation = 
      det::Detector::GetInstance().GetSDetector().GetStation(station);

    const bool hasScintillator = dStation.HasScintillator();
    const bool isUub = dStation.IsUUB();

    typedef VariableBinHistogramWrap<short, int> CalibHistogram;

    for (Station::PMTIterator pmtIt = station.PMTsBegin(sdet::PMTConstants::eAnyType);
      pmtIt != station.PMTsEnd(sdet::PMTConstants::eAnyType); ++pmtIt) {

      if (!pmtIt->HasCalibData())
        continue;

      PMTCalibData& pmtCalibData = pmtIt->GetCalibData();
      StationCalibData& stationCalibData = station.GetCalibData();
      const auto& chargeHistoBinning = dStation.GetMuonChargeHistogramBinning<short>(pmtIt->GetType(), stationCalibData.GetVersion());

      // load raw values from XML
      if (isUub) {
        
        const unsigned int pmtId = pmtIt->GetId();
        rawPeak = fRawPeakUUB[pmtId-1];
        rawCharge = fRawChargeUUB[pmtId-1];
        hgBaseline = fHgBaselineUUB[pmtId-1];
        lgBaseline = fLgBaselineUUB[pmtId-1];
        DARatio = fDARatioUUB[pmtId-1];

      } else if (hasScintillator) {
        
        const unsigned int pmtId = (pmtIt->GetType() == sdet::PMTConstants::eScintillator) ? 2 : 1;
        rawPeak = fRawPeakPPA[pmtId-1];
        rawCharge = fRawChargePPA[pmtId-1];

        // change raw values by online estimates if they exist (only for WCD)
        if (pmtIt->GetType() == sdet::PMTConstants::eWaterCherenkovLarge) {
          double vemPeak = pmtCalibData.GetVEMPeak();
          double vemCharge = pmtCalibData.GetVEMCharge();

          if (vemPeak > 0)
            rawPeak = vemPeak;
          if (vemCharge > 0)
            rawCharge = vemCharge;

          // for WCD we trust the online estimates (from LS)
          pmtCalibData.SetVEMPeak(rawPeak);
          pmtCalibData.SetVEMCharge(rawCharge);
          continue;

        } else if (pmtIt->GetType() == sdet::PMTConstants::eScintillator) {
          pmtCalibData.SetIsTubeOk(true);
          pmtCalibData.SetIsLowGainOk(true);
        }
      }

      // estimate HG and LG baseline for EA stations
      if (isUub && pmtIt->HasFADCTrace()) {

        PMT& pmt = *pmtIt;

        const TraceI& lgTrace = pmt.GetFADCTrace(sdet::PMTConstants::eLowGain, sevt::StationConstants::SignalComponent::eTotal);
        const TraceI& hgTrace = pmt.GetFADCTrace(sdet::PMTConstants::eHighGain, sevt::StationConstants::SignalComponent::eTotal);
      
        // compute baselines of HG and LG traces as median of the first 500 bins
        vector<double> hgPiece;
        vector<double> lgPiece;

        hgPiece.reserve(500);
        lgPiece.reserve(500);

        for (int i = 0; i < 500; ++i) {
          hgPiece.push_back(hgTrace[i]);
          lgPiece.push_back(lgTrace[i]);
        }

        sort(hgPiece.begin(), hgPiece.end());
        sort(lgPiece.begin(), lgPiece.end());

        if (hgPiece.size() % 2 == 0) 
          hgBaseline = (hgPiece[hgPiece.size() / 2 - 1] + hgPiece[hgPiece.size() / 2]) / 2;
        else
          hgBaseline = hgPiece[hgPiece.size() / 2];

        if (lgPiece.size() % 2 == 0) 
          lgBaseline = (lgPiece[lgPiece.size() / 2 - 1] + lgPiece[lgPiece.size() / 2]) / 2;
        else
          lgBaseline = lgPiece[lgPiece.size() / 2];
        

        pmtCalibData.SetBaseline(hgBaseline, 2);
        pmtCalibData.SetBaseline(lgBaseline, 2, sdet::PMTConstants::eLowGain);
        pmtCalibData.SetDynodeAnodeRatio(DARatio, 0);
        pmtCalibData.SetIsTubeOk(true);
        pmtCalibData.SetIsLowGainOk(true);
      } 

      // skip small PMT in VEM / MIP estimations
      if (pmtIt->GetType() == sdet::PMTConstants::eWaterCherenkovSmall)
        continue;

      // find estimate of VEM / MIP charge using calibration histograms
      double estimate = 0;
      double res = 0;
      const int muonChargeHistoSize = pmtCalibData.GetMuonChargeHisto().size();

      if (!muonChargeHistoSize || int(chargeHistoBinning.size()) - 1 != muonChargeHistoSize) {
        WARNING("There should be a muon charge histogram!");
      } else {

        const CalibHistogram chargeHisto(chargeHistoBinning, pmtCalibData.GetMuonChargeHisto());
        const double baseEstimate = 
          pmtCalibData.GetBaseline() * 20 - pmtCalibData.GetMuonChargeHistoOffset();
        const double base = (fabs(baseEstimate) < 20) ? baseEstimate : 0;

        if (chargeHisto.GetMaximum() > 500) {

          int max = 0;
          unsigned int bin = chargeHisto.GetNBins() - 1;
          unsigned int nbbin = bin;
          int v = 0, v1 = 0, v2 = 0;
          unsigned int binmin = 0, binmax = 0;
          int dir = -1;
          // range for the MIP / VEM charge taken from EA Analysis (only UUB stations)
          const double humpvMin = (pmtIt->GetType() == sdet::PMTConstants::eScintillator) ? 100 : 800;
          const double humpvMax = (pmtIt->GetType() == sdet::PMTConstants::eScintillator) ? 600 : 2500;
          // max value "vmax" different between UB and UUB (due to different histogram binning)
          // also different between SSD and WCD
          if (isUub)
            vmax = (pmtIt->GetType() == sdet::PMTConstants::eScintillator) ? 20 : 40;
          else if (hasScintillator)
            vmax = (pmtIt->GetType() == sdet::PMTConstants::eScintillator) ? 100 : 300;
          
          while (!binmax) {
            v  = int(chargeHisto.GetBinAverage(bin));
            v1 = int(chargeHisto.GetBinAverage(bin + dir));
            v2 = int(chargeHisto.GetBinAverage(bin + 2*dir));

            if (v > vmax && v > max) {
              max = v;
            }

            if (binmin && double(max) / v > 1.3 && double(max) / v1 > 1.3 && double(max) / v2 > 1.3)
              binmax = bin;
            if (!binmin && max && double(max) / v > 1.3 && double(max) / v1 > 1.3 && double(max) / v2 > 1.3) {
              binmin = bin;
              dir = 1;
            }
            bin += dir;
            if (bin == nbbin)
              binmax = -1;
            if (bin == 2)
              binmax = -1;
          }

          estimate = (binmin + binmax) / 2;

          if (!isUub) {
            res = estimate;
          } else if (estimate < binmax && estimate > binmin && binmin > 0 && binmax < chargeHisto.GetNBins() && binmin < binmax) {
              double x = 0, x2 = 0, x3 = 0, x4 = 0, y = 0, xy = 0, x2y = 0;
              int nb = 0;
              for (unsigned int i = binmin; i < binmax; ++i) {
                const double v = chargeHisto.GetBinAverage(i);
                const double a = chargeHisto.GetBinCenter(i);
                x += a;
                x2 += a * a;
                x3 += a * a * a;
                x4 += a * a * a * a;
                y += v;
                xy += a * v;
                x2y += a * a * v;
                ++nb;
              }
              
              double B =
                (y * (x4 * x - x2 * x3) + xy * (x2 * x2 - nb * x4) +
                 x2y * (nb * x3 - x * x2));
              double A =
                (y * (x2 * x2 - x * x3) + xy * (nb * x3 - x * x2) +
                 x2y * (x * x - nb * x2));
              res = -B / (2 * A);
              
              if (res > chargeHisto.GetBinCenter(binmin) && res < chargeHisto.GetBinCenter(binmax)) 
                res -= base; 
            }
            
            // check that value is within range for UUB
            if (isUub && !(res > humpvMin && res < humpvMax))
              res = 0; 

        } // if max value > 500
      } // if there is histogram
      
      
      if (res) {
        rawCharge = res / 1.045;
        /*
          A factor is multiplied to the "online" estimate
          of the VEM charge in the SdCalibrator (fOnlineChargeVEMFactor).
          This factor is intended to correct the estimate given by the LS
          (which is off from the value of the "hump" by ~4.5%),
          but it should not be applied to the estimates obtained here.
        */
        if (!isUub)
          rawPeak = rawCharge / 1.8; // hardcoded for SSD-UB (also in CDAS...)
      }
      else {
        rawCharge = -1; // no calibration available
        pmtCalibData.SetIsTubeOk(false);
      }

      pmtCalibData.SetVEMPeak(rawPeak); 
      pmtCalibData.SetVEMCharge(rawCharge);

    }
  }

  void
  SdEACalibrationFiller::FillCalibrationInfo_SmallPMT(Station& station) {

    /*
    Filling SmallPMT calibration (different from standard PMTCalibData)
    -) Beta factor ( Q[VEM] = beta * Q[FADC counts] )
    -) Beta factor error
    -) Chi2 of the spectra comparison calibration procedure
    -) Correction factor ( != 1 only if needed)
    Also SetIsTubeOk for the SmallPMTCalibData class
    The standard PMTCalibData for the SmallPMT (pmt id 4) is left empty 
    */ 
    
    sevt::SmallPMT& spmt = station.GetSmallPMT();
    if(spmt.HasCalibData()) return;
    
    spmt.MakeCalibData();
    sevt::SmallPMTCalibData& spmtCalib = spmt.GetCalibData();
    spmtCalib.SetIsTubeOk(true);
    spmtCalib.SetVersion(0);
    spmtCalib.SetBeta(fBeta);
    spmtCalib.SetBetaError(fBetaErr);
    spmtCalib.SetChi2(fBetaChi2);
    spmtCalib.SetCorrectionFactor(fBetaCorrectionFact);

    //The values w.r.t. each single LPMT are put equal to the overall ones
    for (Station::PMTIterator pmtIt = station.PMTsBegin(sdet::PMTConstants::eWaterCherenkovLarge);
      pmtIt != station.PMTsEnd(sdet::PMTConstants::eWaterCherenkovLarge); ++pmtIt) {
      spmtCalib.SetBeta(fBeta, pmtIt->GetId());
      spmtCalib.SetBetaError(fBetaErr, pmtIt->GetId());
      spmtCalib.SetChi2(fBetaChi2, pmtIt->GetId());
      spmtCalib.SetCorrectionFactor(fBetaCorrectionFact, pmtIt->GetId());
    }

  }

  VModule::ResultFlag
  SdEACalibrationFiller::Finish()
  {
    return eSuccess;
  }

}
