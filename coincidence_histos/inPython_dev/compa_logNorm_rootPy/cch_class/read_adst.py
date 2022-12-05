"""
Created on September 13th, 2022
@author: Mauricio Suarez-Duran
"""

import sys
import numpy as np
from cch_class import adst
import ROOT

class READ_ADST:

    def __init__(self):
        return

    def load_data(self, adst_file):
        dict_cq = {}
        dict_q = {}
        ch_split = 403 # number of channel :403 from where the tail is extended
        isHistosRead = False
        """
        Read calibration histograms and coincidence histograms
        from ADST file        
        """
        ROOT.gROOT.SetBatch(1)  # batch mode, no ROOT graphic
        # read only shower-level observables like energy, direction, ...
        #
        # use adst.FULL to ld
        #everything/Calibration histograms
        degreeOfDetail = adst.FULL 
        print ("Reading ADST data")
        iRecEvent = 0
        for recEvent in adst.RecEventProvider(adst_file, degreeOfDetail):
            sdEvent = recEvent.GetSDEvent()
            stations = sdEvent.GetStationVector()
            #
            # Getting SdRecShower
            ldf = sdEvent.GetSdRecShower().GetLDF()
            energyEvt = sdEvent.GetSdRecShower().GetEnergy()
            energyErrEvt = sdEvent.GetSdRecShower().GetEnergyError()
            zenithEvt = sdEvent.GetSdRecShower().GetZenith()
            for st_i in stations:
                if not st_i.IsUUB():
                    continue
                #
                # Filter for engineer array
                if st_i.GetId() < 100:
                    continue;
                #
                # Reading per PMT from traces
                traces = st_i.GetPMTTraces()
                cCharge = np.zeros(shape=[600, 4], dtype='float')
                charge = np.zeros(shape=[600, 4], dtype='float')
                isHistosRead = False
                qpk = [0, 0, 0]
                qpkErr = [0, 0, 0]
                cqpk = [0, 0, 0]
                cqpkErr = [0, 0, 0]
                isVEMfromHisto = [0, 0, 0]
                vemSignal = [0., 0., 0.]
                ldfVal = 0.
                spDistance = 0.
                enerVal = 0.
                enerErrVal = 0.
                angleVal = 0.
                #
                for tr_i in traces:
                    if tr_i.GetType() != 0: # 0: eTotalTrace
                        continue
                    if tr_i.GetPMTId() > 3:
                        continue
                    histoCQ = tr_i.GetChargeHistogramCoinc().GetValues()
                    histoQ = tr_i.GetChargeHistogram().GetValues()
                    #
                    # Filter wrong histograms
                    if len(histoCQ) < 600 or len(histoQ) < 600:
                        continue
                    if np.sum(histoCQ) < 1e4:
                        continue
                    tmpArgMax = np.argmax(histoCQ)
                    if tmpArgMax > 400 or tmpArgMax < 120:
                        continue
                    #
                    # Kata's filter
                    if tmpArgMax < int(0.06*len(histoCQ)):
                        continue
                    #
                    pmt_i = tr_i.GetPMTId()
                    histoFADC = tr_i.GetChargeHistogramCoinc().GetBinning()
                    cCharge[:, pmt_i] = histoCQ
                    cCharge[:, 0] = histoFADC
                    charge[:, pmt_i] = histoQ
                    charge[:, 0] = histoFADC
                    #
                    # dividing the bins by bin width
                    cCharge[:ch_split, pmt_i] /= 8.
                    charge[:ch_split-3, pmt_i] /= 8.
                    cCharge[ch_split:, pmt_i] /= 32. 
                    charge[ch_split-3:, pmt_i] /= 32.
                    #
                    # 1.01 from chargeConversionFactor
                    qpk[pmt_i-1] = tr_i.GetCharge()*1.01
                    qpkErr[pmt_i-1] = tr_i.GetChargeError()*1.01
                    cqpk[pmt_i-1] = tr_i.GetChargeCoinc()
                    cqpkErr[pmt_i-1] = tr_i.GetChargeCoincError()
                    isVEMfromHisto[pmt_i-1] = tr_i.IsVEMChargeFromHistogram()
                    isHistosRead = True
                    vemSignal[pmt_i-1] = tr_i.GetVEMSignal()
                    spDistance = st_i.GetSPDistance()
                    ldfVal = ldf.Evaluate(spDistance, 0)
                    enerVal = energyEvt
                    enerErrVal = energyErrEvt
                    angleVal = zenithEvt
                #
                # Writing into dictionary
                if isHistosRead:
                    keyName = str(st_i.GetId())+' '+str(sdEvent.GetGPSSecond())
                    dict_cq[keyName] = {}
                    dict_cq[keyName]['cqpk'] = cqpk
                    dict_cq[keyName]['cqpkErr'] = cqpkErr
                    dict_cq[keyName]['fromHistogram'] = isVEMfromHisto
                    dict_cq[keyName]['hist'] = cCharge.tolist()
                    dict_cq[keyName]['signal'] = vemSignal

                    dict_q[keyName] = {}
                    dict_q[keyName]['qpk'] = qpk
                    dict_q[keyName]['qpkErr'] = qpkErr
                    dict_q[keyName]['fromHistogram'] = isVEMfromHisto
                    dict_q[keyName]['hist'] = charge.tolist()
                    dict_q[keyName]['signal'] = vemSignal
                    dict_q[keyName]['ldfR'] = ldfVal
                    dict_q[keyName]['spDist'] = spDistance
                    dict_q[keyName]['Energy'] = enerVal
                    dict_q[keyName]['EnergyErr'] = enerErrVal
                    dict_q[keyName]['Zenith'] = angleVal
        return dict_cq, dict_q