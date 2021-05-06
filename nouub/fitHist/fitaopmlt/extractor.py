import numpy as np

class extractor:
    def __init__(self, hist):
        self.hist = hist # Receive data
        self.lstentry = 32000 #hist.GetEntries() # Set the total of entries

        self.distaopHb = [] # It stores the dist. AoP for HBase
        self.aveaopHb = [] # It stores Average AoP for HBase
        self.rmsaopHb = [] # It stores RMS AoP for HBase

        self.distaopCa = [] # It stores the dist. AoP for Calib
        self.aveaopCa = [] # It stores Average AoP for Calib
        self.rmsaopCa = [] # It stores RMS AoP for Calib

        self.totEvts = [] # It stores total events per Station
        
        self.distXiPkHb = [] # It stores the dist. Xi for HBase
        self.aveXiPkHb = [] # It stores Average Xi for HBase
        self.rmsXiPkHb = [] # It stores RMS Xi for HBase

        self.distXiPkCa = [] # It stores the dist. Xi for Calib
        self.aveXiPkCa = [] # It stores Average Xi for Calib
        self.rmsXiPkCa = [] # It stores RMS AoP Xi Calib

        self.distXiChHb = [] # It stores the dist. Xi for HBase
        self.aveXiChHb = [] # It stores Average Xi for HBase
        self.rmsXiChHb = [] # It stores RMS Xi for HBase

        self.distXiChCa = [] # It stores the dist. Xi for Calib
        self.aveXiChCa = [] # It stores Average Xi for Calib
        self.rmsXiChCa = [] # It stores RMS AoP Xi Calib
        
    # =========================
    # *** For AoP for HBase ***
    def getAoPDistHb(self, listStat):
        print("Doing AoP Distribution for HBase")
        tmpdist = []
        tmpentry = 0
        for st in range(0, len(listStat)): 
            tmpdist = []
            for evt in range(0, self.lstentry):
                self.hist.GetEntry(evt)
                if tmpentry != self.hist.eventStat[st+1]:
                    tmpentry = self.hist.eventStat[st+1]
                    if self.hist.apHbase.GetBinContent(st) > 0:
                        tmpdist.append(self.hist.apHbase.GetBinContent(st))
            if len(tmpdist) > 0:
                self.distaopHb.append(tmpdist)
            else:
                 self.distaopHb.append(0)

    def getAoPAveRmsHb(self):
        print("Doing AoP Average and RMS for HBase")
        tmpmean = 0
        tmprms = 0
        for st in range(0, len(self.distaopHb)):
            tmpmean = np.average(np.array(self.distaopHb[st]))
            tmprms = 0
            if type (self.distaopHb[st]) is not int:
                for bn in self.distaopHb[st]:
                    tmprms += (bn - tmpmean)*(bn - tmpmean)
                self.rmsaopHb.append( np.sqrt(tmprms/len(self.distaopHb[st])) )
                self.aveaopHb.append(tmpmean)
            else:
                self.rmsaopHb.append(0)
                self.aveaopHb.append(0)

    # =========================
    # *** For AoP for Calib ***
    def getAoPDistCa(self, listStat):
        print("Doing AoP Distribution for Calib")
        tmpdist = []
        tmpentry = 0
        for st in range(0, len(listStat)):
            tmpdist = []
            for evt in range(0, self.lstentry):
                self.hist.GetEntry(evt)
                if tmpentry != self.hist.eventStat[st+1]:
                    tmpentry = self.hist.eventStat[st+1]
                    if self.hist.apCalib.GetBinContent(st) > 0:
                        tmpdist.append(self.hist.apCalib.GetBinContent(st))
            if len(tmpdist) > 0:
                self.distaopCa.append(tmpdist)
            else:
                self.distaopCa.append(0)

    def getAoPAveRmsCa(self):
        print("Doing AoP Average and RMS for Calib")
        tmpmean = 0
        tmprms = 0
        for st in range(0, len(self.distaopCa)):
            tmpmean = np.average(np.array(self.distaopCa[st]))
            tmprms = 0
            if type(self.distaopCa[st]) is not int :
                for bn in self.distaopCa[st]:
                    tmprms += (bn - tmpmean)*(bn - tmpmean)
                self.rmsaopCa.append( np.sqrt(tmprms/len(self.distaopCa[st])) )
                self.aveaopCa.append(tmpmean)
            else:
                self.rmsaopCa.append(0)
                self.aveaopCa.append(0)

    # ============================
    # *** For Total of Events ***
    def getTotEvents(self, listStat):
        tmpEntries = 0
        for st in range(0, len(listStat)):
            tmpEntries = 0
            for evt in range(0, self.lstentry):
                self.hist.GetEntry(evt)
                tmpEntries = self.hist.eventStat.GetBinContent(st+1)
            self.totEvts.append( tmpEntries )

    # ===================
    # *** For Xi-Peak ***

    # *** For HBase ***
    def getXiPkDistHb(self, listStat):
        print("Doing Xi-Peak Distribution for HBase")
        tmpdist = []
        tmpentry = 0
        for st in range(0, len(listStat)):
            tmpdist = []
            for evt in range(0, self.lstentry):
                self.hist.GetEntry(evt)
                if tmpentry != self.hist.eventStat[st+1]:
                    tmpentry = self.hist.eventStat[st+1]
                    if self.hist.apHbase.GetBinContent(st) > 0:
                        tmpdist.append(self.hist.chisHbasePk.GetBinContent(st))
            if len(tmpdist) > 0:
                self.distXiPkHb.append(tmpdist)
            else:
                self.distXiPkHb.append(0)

    def getXiPkveRmsHb(self):
        print("Doing Xi-Peak Average and RMS for HBase")
        tmpmean = 0
        tmprms = 0
        for st in range(0, len(self.distXiPkHb)):
            tmpmean = np.average(np.array(self.distXiPkHb[st]))
            tmprms = 0
            if type(self.distXiPkHb[st]) is not int :
                for bn in self.distXiPkHb[st]:
                    tmprms += (bn - tmpmean)*(bn - tmpmean)
                self.rmsXiPkHb.append( np.sqrt(tmprms/len(self.distXiPkHb[st])) )
                self.aveXiPkHb.append(tmpmean)
            else:
                self.rmsXiPkHb.append(0)
                self.aveXiPkHb.append(0)

    # *** For Calib ***
    def getXiPkDistCa(self, listStat):
        print("Doing Xi-Peak Distribution for Calib")
        tmpdist = []
        tmpentry = 0
        for st in range(0, len(listStat)):
            tmpdist = []
            for evt in range(0, self.lstentry):
                self.hist.GetEntry(evt)
                if tmpentry != self.hist.eventStat[st+1]:
                    tmpentry = self.hist.eventStat[st+1]
                    if self.hist.apCalib.GetBinContent(st) > 0:
                        tmpdist.append(self.hist.chisCalibPk.GetBinContent(st))
            if len(tmpdist) > 0:
                self.distXiPkCa.append(tmpdist)
            else:
                self.distXiPkCa.append(0)

    def getXiPkveRmsCa(self):
        print("Doing Xi-Peak Average and RMS for Calib")
        tmpmean = 0
        tmprms = 0
        for st in range(0, len(self.distXiPkCa)):
            tmpmean = np.average(np.array(self.distXiPkCa[st]))
            tmprms = 0
            if type(self.distXiPkCa[st]) is not int :
                for bn in self.distXiPkCa[st]:
                    tmprms += (bn - tmpmean)*(bn - tmpmean)
                self.rmsXiPkCa.append( np.sqrt(tmprms/len(self.distXiPkCa[st])) )
                self.aveXiPkCa.append(tmpmean)
            else:
                self.rmsXiPkCa.append(0)
                self.aveXiPkCa.append(0)

    # ===================
    # *** For Xi-Charge ***

    # *** For HBase ***
    def getXiChDistHb(self, listStat):
        print("Doing Xi-Charge Distribution for HBase")
        tmpdist = []
        tmpentry = 0
        for st in range(0, len(listStat)):
            tmpdist = []
            for evt in range(0, self.lstentry):
                self.hist.GetEntry(evt)
                if tmpentry != self.hist.eventStat[st+1]:
                    tmpentry = self.hist.eventStat[st+1]
                    if self.hist.chisHbaseCh.GetBinContent(st) > 0:
                        tmpdist.append(self.hist.chisHbaseCh.GetBinContent(st))
            if len(tmpdist) > 0:
                self.distXiChHb.append(tmpdist)
            else:
                self.distXiChHb.append(0)

    def getXiChveRmsHb(self):
        print("Doing Xi-Charge Average and RMS for HBase")
        tmpmean = 0
        tmprms = 0
        for st in range(0, len(self.distXiChHb)):
            tmpmean = np.average(np.array(self.distXiChHb[st]))
            tmprms = 0
            if type(self.distXiChHb[st]) is not int :
                for bn in self.distXiChHb[st]:
                    tmprms += (bn - tmpmean)*(bn - tmpmean)
                self.rmsXiChHb.append( np.sqrt(tmprms/len(self.distXiChHb[st])) )
                self.aveXiChHb.append(tmpmean)
            else:
                self.rmsXiChHb.append(0)
                self.aveXiChHb.append(0)

    # *** For Calib ***
    def getXiChDistCa(self, listStat):
        print("Doing Xi-Charge Distribution for Calib")
        tmpdist = []
        tmpentry = 0
        for st in range(0, len(listStat)):
            tmpdist = []
            for evt in range(0, self.lstentry):
                self.hist.GetEntry(evt)
                if tmpentry != self.hist.eventStat[st+1]:
                    tmpentry = self.hist.eventStat[st+1]
                    if self.hist.chisCalibCh.GetBinContent(st) > 0:
                        tmpdist.append(self.hist.chisCalibCh.GetBinContent(st))
            if len(tmpdist) > 0:
                self.distXiChCa.append(tmpdist)
            else:
                self.distXiChCa.append(0)

    def getXiChveRmsCa(self):
        print("Doing Xi-Charge Average and RMS for Calib")
        tmpmean = 0
        tmprms = 0
        for st in range(0, len(self.distXiChCa)):
            tmpmean = np.average(np.array(self.distXiChCa[st]))
            tmprms = 0
            if type(self.distXiChCa[st]) is not int :
                for bn in self.distXiChCa[st]:
                    tmprms += (bn - tmpmean)*(bn - tmpmean)
                self.rmsXiChCa.append( np.sqrt(tmprms/len(self.distXiChCa[st])) )
                self.aveXiChCa.append(tmpmean)
            else:
                self.rmsXiChCa.append(0)
                self.aveXiChCa.append(0)
