void makeForIoana(bool ifIsUub, int pmt, int st_id) {
    
  TString bnCdas = (ifIsUub) ?
    "~/2021/sdeu/underHisto/results/uubChPkPMT" :
    "~/2021/sdeu/nouub/underHistos/results/ubChPkPMT"; 

  TString strIfUub = (ifIsUub) ? "Uub" : "Ub";
  TString strPmt;
  TString strSt;

  strPmt.Form("%d", pmt);
  strSt.Form("%d", st_id);

  TFile *outRoot = new TFile("qpks"+strIfUub+"Pmt"+strPmt+"St"+strSt+".root","RECREATE");
  double retQpk = 0.;
  int retTim = 0;
  TTree *treeQpkVals = new TTree("QpkValues","");
  treeQpkVals->Branch("qpk",&retQpk,"retQpk/D");
  treeQpkVals->Branch("timeGPS",&retTim,"retTim/I");

  TString pmtId;
  TString strStId;
  TString fname;
  pmtId.Form("%d", pmt);
  strStId.Form("St%d", (int)st_id);

  TString monthUub[4] = {"Aug", "Sep", "Oct", "Nov"};
  int nMonths = (ifIsUub) ?
    sizeof(monthUub)/sizeof(*monthUub) : 
    sizeof(monthUub)/sizeof(*monthUub)-1;
  
  int nYears = 5;
  TString strYear[5] = {"2016", "2018", "2019", "2020", "2021"};
  int stYear = (ifIsUub) ? nYears-1 : 1;
  int lstYear = (ifIsUub) ? nYears-1 : 1;
  
  TString strChargeData = "ChargeData";
  TFile *f;
  TTree *chargeInfo;
  double fetchQpkVals = 0.;
  int fetchTime = 0;

  for ( int year=stYear; year<=lstYear; year++ ) {
    for ( int month_i=0; month_i<nMonths; month_i++ ) {
      fname = bnCdas + pmtId + strStId + "lrb35" + monthUub[month_i] + strYear[year];
      f = TFile::Open(fname+".root");
      chargeInfo = (TTree*)f->Get(strChargeData);
      chargeInfo->SetBranchAddress("chargeVal", &fetchQpkVals);
      chargeInfo->SetBranchAddress("timeEvnt", &fetchTime);
      fetchQpkVals = 0.;
      fetchTime = 0;
      for( int etry=0; etry<chargeInfo->GetEntries(); etry++) {
        chargeInfo->GetEntry(etry);
        if ( fetchQpkVals <= 0 )
          continue;
        retQpk = fetchQpkVals ;
        retTim = fetchTime;
        treeQpkVals->Fill();
      }
      f->Clear();
      f->Close();
    }
  }
  outRoot->Write();
  outRoot->Close();
  exit(0);
}
