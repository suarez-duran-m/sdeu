import numpy as np

class extractor:
    def __init__(self, hist):
        self.hist = hist # Receive data
        self.lstentry = 32000 #hist.GetEntries() # Set the total of entries

        self.difstoffset = [] # it stores offset distribution per station for Pk
        self.avestoffset = [] # it stores offset average per station for Pk
        self.rmsstoffset = [] # it stores offset rms values per station for Pk

        self.distoffsetCh = [] # it stores offset distribution per station for Ch
        self.avestoffsetCh = [] # it stores offset average per station for Ch
        self.rmsstoffsetCh = [] # it stores offset rms values per station for Ch

        self.distFrstBin = [] # It stores first bin distribution per Station
        self.aveFrstBin = [] # It stores Average first bin for Peak
        self.rmsFrstBin = [] # It stores RMS first bin for Peak

        self.distHbase = [] # It stores HBase per Station
        self.aveHbase = [] # It stores Average HBase
        self.rmsHbase = [] # It stores RMS HBase

        self.distCalib = [] # It stores Calib per Station
        self.aveCalib = [] # It stores Average Calib
        self.rmsCalib = [] # It stores RMS Calib

        self.distBinCntrHbase = [] # It stores Calib per Station
        self.aveBinCntrHbase = [] # It stores Average Calib
        self.rmsBinCntrHbase = [] # It stores RMS Calib

        self.distBinCntrCalib = [] # It stores Calib per Station
        self.aveBinCntrCalib = [] # It stores Average Calib
        self.rmsBinCntrCalib = [] # It stores RMS Calib

    # ==================
    # *** For Offset ***
    def getOffsetDist(self, listStat):
        print("Doing Distribution for Offset")
        tmpdist = []
        tmpentry = 0
        for st in range(0, len(listStat)):
            tmpdist = []
            for evt in range(0, self.lstentry):
                self.hist.GetEntry(evt)
                if tmpentry != self.hist.eventStat[st+1]:
                    tmpentry = self.hist.eventStat[st+1]
                    tmpdist.append(self.hist.offSetPk.GetBinContent(st))
            self.difstoffset.append(tmpdist)

    def getAveOffsetAveRms(self):
        print("Doing Average and RMS for Offset")
        tmpmean = 0
        tmprms = 0
        for st in range(0, len(self.difstoffset )):
            tmpmean = np.average(np.array(self.difstoffset [st]))
            tmprms = 0
            for bn in self.difstoffset[st]:
                tmprms += (bn - tmpmean)*(bn - tmpmean)
            self.rmsstoffset.append( np.sqrt(tmprms/len(self.difstoffset[st])) )
            if (tmpmean != 0):
                self.avestoffset.append(tmpmean)
            else:
                self.avestoffset.append(0)

    # *** For Charge ***
    def getOffsetDistCh(self, listStat):
        print("Doing Distribution for Offset for Charge")
        tmpdist = []
        tmpentry = 0
        for st in range(0, len(listStat)):
            tmpdist = []
            for evt in range(0, self.lstentry):
                self.hist.GetEntry(evt)
                if tmpentry != self.hist.eventStat[st+1]:
                    tmpentry = self.hist.eventStat[st+1]
                    tmpdist.append(self.hist.offSetCh.GetBinContent(st))
            self.distoffsetCh.append(tmpdist)

    def getAveOffsetAveRmsCh(self):
        print("Doing Average and RMS for Offset for Charge")
        tmpmean = 0
        tmprms = 0
        for st in range(0, len(self.distoffsetCh)):
            tmpmean = np.average(np.array(self.distoffsetCh[st]))
            tmprms = 0
            for bn in self.distoffsetCh[st]:
                tmprms += (bn - tmpmean)*(bn - tmpmean)
            self.rmsstoffsetCh.append( np.sqrt(tmprms/len(self.distoffsetCh[st])) )
            if (tmpmean != 0):
                self.avestoffsetCh.append(tmpmean)
            else:
                self.avestoffsetCh.append(0)

    # =====================
    # *** For First Bin ***
    def getFrstBinDist(self, listStat):
        print("Doing Distribution for first bin of Peak")
        tmpdist = []
        tmpentry = 0
        for st in range(0, len(listStat)):
            tmpdist = []
            totetrySt = 0
            #if st == 16:
                #averageStat.append(0)
                #continue
            for evt in range(0, self.lstentry):
                self.hist.GetEntry(evt)
                if tmpentry != self.hist.eventStat[st+1]:
                    tmpentry = self.hist.eventStat[st+1]
                    tmpdist.append( self.hist.firstBinCntPk.GetBinContent(st))
            self.distFrstBin.append( tmpdist )

    def getFrstBinAveRms(self):
        print("Doing Average and RMS for first bin") 
        tmpmean = 0
        tmprms = 0
        for st in range(0, len(self.distFrstBin)):
            tmpmean = np.average(np.array(self.distFrstBin[st]))
            tmprms = 0
            for bn in self.distFrstBin[st]:
                tmprms += (bn - tmpmean)*(bn - tmpmean)
            self.rmsFrstBin.append( np.sqrt(tmprms/len(self.distFrstBin[st])) )
            if (tmpmean != 0):
                self.aveFrstBin.append(tmpmean)
            else:
                self.aveFrstBin.append(0)


    # ====================
    # *** For Baseline ***

    # *** For HBase ***
    def getHbaseDist(self, listStat):
        print("Doing Distribution for HBase")
        tmpdist = []
        tmpentry = 0
        for st in range(0, len(listStat)):
            tmpdist = []
            totetrySt = 0
            for evt in range(0, self.lstentry):
                self.hist.GetEntry(evt)
                if tmpentry != self.hist.eventStat[st+1]:
                    tmpentry = self.hist.eventStat[st+1]
                    tmpdist.append( self.hist.baselineHbase.GetBinContent(st))
            self.distHbase.append( tmpdist )

    def getHbaseAveRms(self):
        print("Doing Average and RMS for HBase") 
        tmpmean = 0
        tmprms = 0
        for st in range(0, len(self.distHbase)):
            tmpmean = np.average(np.array(self.distHbase[st]))
            tmprms = 0
            for bn in self.distHbase[st]:
                tmprms += (bn - tmpmean)*(bn - tmpmean)
            self.rmsHbase.append( np.sqrt(tmprms/len(self.distHbase[st])) )
            if (tmpmean != 0):
                self.aveHbase.append(tmpmean)
            else:
                self.aveHbase.append(0)

    # *** For Calib ***
    def getCalibDist(self, listStat):
        print("Doing Distribution for Calib")
        tmpdist = []
        tmpentry = 0
        for st in range(0, len(listStat)):
            tmpdist = []
            totetrySt = 0
            for evt in range(0, self.lstentry):
                self.hist.GetEntry(evt)
                if tmpentry != self.hist.eventStat[st+1]:
                    tmpentry = self.hist.eventStat[st+1]
                    tmpdist.append( self.hist.baselineCalib.GetBinContent(st))
            self.distCalib.append( tmpdist )

    def getCalibAveRms(self):
        print("Doing Average and RMS for Calib") 
        tmpmean = 0
        tmprms = 0
        for st in range(0, len(self.distCalib)):
            tmpmean = np.average(np.array(self.distCalib[st]))
            tmprms = 0
            for bn in self.distCalib[st]:
                tmprms += (bn - tmpmean)*(bn - tmpmean)
            self.rmsCalib.append( np.sqrt(tmprms/len(self.distCalib[st])) )
            if (tmpmean != 0):
                self.aveCalib.append(tmpmean)
            else:
                self.aveCalib.append(0)

    # ============================
    # *** For First bin center ***

    # *** For HBase ***
    def getBinCenterDistHbase(self, listStat):
        print("Doing Distribution for first bin center HBase")
        tmpdist = []
        tmpentry = 0
        for st in range(0, len(listStat)):
            tmpdist = []
            for evt in range(0, self.lstentry):
                self.hist.GetEntry(evt)
                if tmpentry != self.hist.eventStat[st+1]:
                    tmpentry = self.hist.eventStat[st+1]
                    tmpdist.append(self.hist.bincentHbasePk.GetBinContent(st))
            self.distBinCntrHbase.append(tmpdist)

    def getBinCenterAveRmsHbase(self):
        print("Doing Average and RMS for first bin center HBase")
        tmpmean = 0
        tmprms = 0
        for st in range(0, len(self.distBinCntrHbase)):
            tmpmean = np.average(np.array(self.distBinCntrHbase[st]))
            tmprms = 0
            for bn in self.distBinCntrHbase[st]:
                tmprms += (bn - tmpmean)*(bn - tmpmean)
            self.rmsBinCntrHbase.append( np.sqrt(tmprms/len(self.distBinCntrHbase[st])) )
            if (tmpmean != 0):
                self.aveBinCntrHbase.append(tmpmean)
            else:
                self.aveBinCntrHbase.append(0)

    # *** For Calib ***
    def getBinCenterDistCalib(self, listStat):
        print("Doing Distribution for first bin center Calib")
        tmpdist = []
        tmpentry = 0
        for st in range(0, len(listStat)):
            tmpdist = []
            for evt in range(0, self.lstentry):
                self.hist.GetEntry(evt)
                if tmpentry != self.hist.eventStat[st+1]:
                    tmpentry = self.hist.eventStat[st+1]
                    tmpdist.append(self.hist.bincentCalibPk.GetBinContent(st))
            self.distBinCntrCalib.append(tmpdist)

    def getBinCenterAveRmsCalib(self):
        print("Doing Average and RMS for first bin center Calib")
        tmpmean = 0
        tmprms = 0
        for st in range(0, len(self.distBinCntrCalib)):
            tmpmean = np.average(np.array(self.distBinCntrCalib[st]))
            tmprms = 0
            for bn in self.distBinCntrCalib[st]:
                tmprms += (bn - tmpmean)*(bn - tmpmean)
            self.rmsBinCntrCalib.append( np.sqrt(tmprms/len(self.distBinCntrCalib[st])) )
            if (tmpmean != 0):
                self.aveBinCntrCalib.append(tmpmean)
            else:
                self.aveBinCntrCalib.append(0)
