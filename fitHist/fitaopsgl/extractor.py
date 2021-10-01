import numpy as np

class extractor:
    def __init__(self, hist):
        self.hist = hist # Receive data
        self.lstentry = hist.GetEntries() # Set the total of entries

        self.distaopHb = [] # It stores the dist. AoP for HBase
        self.aveaopHb = [] # It stores Average AoP for HBase
        self.rmsaopHb = [] # It stores RMS AoP for HBase

        self.distaopCa = [] # It stores the dist. AoP for Calib
        self.aveaopCa = [] # It stores Average AoP for Calib
        self.rmsaopCa = [] # It stores RMS AoP for Calib

        self.totEvts = 0 # It stores total events per Station
        
        self.nfitokAoPHb = 0 # It stores the dist. Xi for HBase
        self.nfitokAoPCa = 0 # It stores the dist. Xi for HBase
        self.nfitokCh = 0 # It stores the dist. Xi for HBase

    # =========================
    # *** For AoP for HBase ***
    def getAoPDistHb(self):
        print("Doing AoP Distribution for HBase")
        tmpentry = 0
        for evt in range(0, self.lstentry):
            self.hist.GetEntry(evt)
            if tmpentry != self.hist.eventStat[1]:
                tmpentry = self.hist.eventStat[1]
                if self.hist.apHbase.GetBinContent(0) > 0:
                    self.distaopHb.append(self.hist.apHbase.GetBinContent(0))

    def getAoPAveRmsHb(self):
        print("Doing AoP Average and RMS for HBase")
        tmprms = 0
        if len(self.distaopHb) > 0:
            tmpmean = np.average(np.array(self.distaopHb))
            for bn in self.distaopHb:
                tmprms += (bn - tmpmean)*(bn - tmpmean)
            self.rmsaopHb.append( np.sqrt(tmprms/len(self.distaopHb)) )
            self.aveaopHb.append(tmpmean)
        else:
            self.rmsaopHb.append(0)
            self.aveaopHb.append(0)


    # =========================
    # *** For AoP for Calib ***
    def getAoPDistCa(self):
        print("Doing AoP Distribution for Calib")
        tmpentry = 0
        for evt in range(0, self.lstentry):
            self.hist.GetEntry(evt)
            if tmpentry != self.hist.eventStat[1]:
                tmpentry = self.hist.eventStat[1]
                if self.hist.apCalib.GetBinContent(0) > 0:
                    self.distaopCa.append(self.hist.apCalib.GetBinContent(0))

    def getAoPAveRmsCa(self):
        print("Doing AoP Average and RMS for Calib")
        tmprms = 0
        if len(self.distaopCa) > 0:
            tmpmean = np.average(np.array(self.distaopCa))
            tmprms = 0
            for bn in self.distaopCa:
                tmprms += (bn - tmpmean)*(bn - tmpmean)
            self.rmsaopCa.append( np.sqrt(tmprms/len(self.distaopCa)) )
            self.aveaopCa.append(tmpmean)
        else:
            self.rmsaopCa.append(0)
            self.aveaopCa.append(0)

    # ============================
    # *** For Total of Events ***
    def getTotEvents(self):
        self.totEvts = 0
        for evt in range(0, self.lstentry):
            self.hist.GetEntry(evt)
            self.totEvts = self.hist.eventStat.GetBinContent(1)

        
        return self.totEvts

    # ===================
    # *** For Fit-Oks ***

    # *** For HBase ***
    def getFitokHb(self):
        print("Doing Get Fits Ok AoP for HBase")
        tmpentry = 0
        self.nfitokAoPHb = 0
        for evt in range(0, self.lstentry):
            self.hist.GetEntry(evt)
            if tmpentry != self.hist.eventStat[1]:
                tmpentry = self.hist.eventStat[1]
                if self.hist.apHbase.GetBinContent(0) > 0:
                    self.nfitokAoPHb += 1
        return self.nfitokAoPHb

    # *** For Calib ***
    def getFitokCa(self):
        print("Doing Get Fits Ok AoP for Calib")
        tmpentry = 0
        self.nfitokAoPCa = 0
        for evt in range(0, self.lstentry):
            self.hist.GetEntry(evt)
            if tmpentry != self.hist.eventStat[1]:
                tmpentry = self.hist.eventStat[1]
                if self.hist.apCalib.GetBinContent(0) > 0:
                    self.nfitokAoPCa += 1
        return self.nfitokAoPCa



