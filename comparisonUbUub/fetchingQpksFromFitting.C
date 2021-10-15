TH1D *getQpkValues( TString bname, double StId, int pmt, bool ifIsUub,
   int nbins, double frstBin, double lstBin ) {

  TString monthUub[2] = {"Aug", "Sep"};
  TString pmtId;
  pmtId.Form("%d", pmt);
  TString strStId;
  strStId.Form("St%d", (int)StId);
  TString strYear[3] = {"2019", "2020", "2021"};
  TString fname;
  int stYear = (ifIsUub) ? 2 : 0;
  int lstYear = (ifIsUub) ? 2 : 1;

  TString strChargeData = "ChargeData";
  TFile *f;
  TTree *chargeInfo;
  double fetchQpkVals = 0.;
  TH1D *retQpkDist = new TH1D("retQpkDist"+strStId+"Pmt"+pmtId,"", 
      nbins, frstBin, lstBin);

  for ( int year=stYear; year<=lstYear; year++ ) {
    for ( int month_i=0; month_i<sizeof(monthUub)/sizeof(*monthUub); month_i++ ) {
      fname = bname + pmtId + strStId + "lrb35" + monthUub[month_i] + strYear[year];
      f = TFile::Open(fname+".root");
      chargeInfo = (TTree*)f->Get(strChargeData);
      chargeInfo->SetBranchAddress("chargeVal", &fetchQpkVals);
      
      fetchQpkVals = 0.;
      for( int etry=0; etry<chargeInfo->GetEntries(); etry++) {
        chargeInfo->GetEntry(etry);
        if ( fetchQpkVals == 0 || fetchQpkVals < 0 )
          continue;
        retQpkDist->Fill( fetchQpkVals );
      }
      f->Clear();
      f->Close();
    }
  }
  chargeInfo->Delete();
  delete f;
  return retQpkDist;
}

void fetchingQpksFromFitting(bool ifUub) {

  //string stListPath = "/home/msd/2021/sdeu/listStationsShort.txt";
  string stListPath = "/home/msd/2021/sdeu/fullUubStationsListVert.txt";
  TString bnCdas = (ifUub) ?
    "~/2021/sdeu/underHisto/results/uubChPkPMT" :
    "~/2021/sdeu/nouub/underHistos/results/ubChPkPMT";

  ifstream fileStList;
  double St_i;
  vector < double > stListId;
  fileStList.open(stListPath);
  while( fileStList.good() ) {
    fileStList >> St_i;
    stListId.push_back(St_i);
  }
  stListId.pop_back();
  fileStList.close();

  int nPmts = 3;
  TTree *stations;
  TString typeStation = (ifUub) ? "StationsUub" : "StationsUb";
  TFile *f = new TFile("qpkVal"+typeStation+".root","recreate");
  stations = new TTree(typeStation,"");
  int nbins = (ifUub) ? 3000 : 300;
  double frstBin = 0.;
  double lstBin = (ifUub) ? 3000. : 300.; 
  int stid;
  TH1D *qpkValuesPmt1 = new TH1D();
  TH1D *qpkValuesPmt2 = new TH1D();
  TH1D *qpkValuesPmt3 = new TH1D();
  stations->Branch("qpkValuesPmt1","TH1D", &qpkValuesPmt1, 32000, 0);
  stations->Branch("qpkValuesPmt2","TH1D", &qpkValuesPmt2, 32000, 0);
  stations->Branch("qpkValuesPmt3","TH1D", &qpkValuesPmt3, 32000, 0);
  stations->Branch("stId",&stid, "stid/I"); 

  for ( auto & st_i : stListId ) {
    qpkValuesPmt1 = getQpkValues(bnCdas, st_i, 1, ifUub, nbins, frstBin, lstBin);
    qpkValuesPmt2 = getQpkValues(bnCdas, st_i, 2, ifUub, nbins, frstBin, lstBin);
    qpkValuesPmt3 = getQpkValues(bnCdas, st_i, 3, ifUub, nbins, frstBin, lstBin);
    stid = (int)st_i;
    stations->Fill();
  }

  f->Write();
  f->Close();
  exit(0);
} 
