TCanvas *canvasStyle(TString name) {
  TCanvas *canvas = new TCanvas(name, name, 1600, 900);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(2);
  canvas->SetLeftMargin(0.11); 
  canvas->SetRightMargin(0.03);
  canvas->SetTopMargin(0.02); 
  canvas->SetBottomMargin(0.15);
  canvas->SetFrameBorderMode(0);
  return canvas;
}

void histoStyle(TH1F *hist) {
  hist->GetXaxis()->SetTitleOffset(1.3);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleOffset(0.9);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
}

void histoStyle(TGraphErrors *hist) {
  hist->GetXaxis()->SetTitleOffset(1.3);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleOffset(1.1);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
}

double getmean( vector<double> arr ) {
  double mean = 0.;
  int nb = arr.size();
  int ngoodb = 0;
    for (auto & elem : arr)
      if ( elem > 0 ) { 
        mean += elem;
        ngoodb++;
      }
  return mean/ngoodb;
}

double getrms( vector<double> arr, double meanarr ) {
  double rms = 0.;
  int nb = arr.size();
  int ngoodb = 0;
  for (auto & elem : arr)
    if ( elem > 0 ) {
      rms += (elem - meanarr)*(elem - meanarr);
      ngoodb++;
    }
  return sqrt(rms/ngoodb)/meanarr;
}

void getEvtIdQpk( TString bname, TString st, int pmt, bool ifIsUub, bool ifoff, 
   vector < int > &retId, vector < double > &retQpk ) {

  TString monthUub[] = {"Aug"};
  TString pmtId;
  TString strYear[] = {"2019", "2020", "2021"};
  TString fname;
  int stYear = (ifIsUub) ? 2 : 0;
  int lstYear = (ifIsUub) ? 2 : 1;
  // Sept. 1st. for 2019, 2020, and 2021
  int beforeSept [3] = {1251331218, 1282953618, 1314489618};
  TString strChargeData = (ifoff) ? "charge" : "ChargeData";
  pmtId.Form("%d", pmt);

  TFile *f;
  TTree *chargeInfo;
  double getQpkVals = 0.;
  int month = 0;

  int evttime = 0;
  int prevTime = 0;
  int EvtId = 0;

  for ( int year=stYear; year<=lstYear; year++ ) {
    fname = (ifoff) ? 
      (bname + monthUub[month] + strYear[year] + st + "Pmt" + pmtId)
      : (bname + pmtId + st + "lrb35" + monthUub[month] + strYear[year]);
   
    f = TFile::Open(fname+".root");
    chargeInfo = (TTree*)f->Get(strChargeData);
    chargeInfo->SetBranchAddress("chargeVal", &getQpkVals); 
    if ( ifoff ) {
      chargeInfo->SetBranchAddress("GpsTime", &evttime);
      chargeInfo->SetBranchAddress("evtId", &EvtId);
    }
    else {
      chargeInfo->SetBranchAddress("timeEvnt", &evttime);
      chargeInfo->SetBranchAddress("eventId", &EvtId);
    }

    prevTime = 0;
    getQpkVals = 0.;
    for( int etry=0; etry<chargeInfo->GetEntries(); etry++) {
      chargeInfo->GetEntry(etry);
      if ( evttime < beforeSept[year] ) { // Events before Sept. 1st, 2021
        // To avoid read the same event twice from Offline
        if ( ifoff ) {
          if ( prevTime != evttime ) {
            retId.push_back( EvtId );
            retQpk.push_back( getQpkVals );
            prevTime = evttime;
          }
        }
        else {
          retId.push_back( EvtId );
          retQpk.push_back( getQpkVals );
        }
      }
    }
  }

  chargeInfo->Delete();
  f->Delete();
}

void doMatch(vector<int> evtIdOff, vector<int> evtIdCdas, 
    vector<double> QpkOff, vector<double> QpkCdas, 
    vector<double> &matchOff, vector<double> &matchCdas) {

  for ( int idOff=0; idOff<evtIdOff.size(); idOff++ )
    for ( int idCdas=0; idCdas<evtIdCdas.size(); idCdas++ )
      if ( evtIdOff[idOff] == evtIdCdas[idCdas] ) {
        matchOff.push_back( QpkOff[idOff] );
        matchCdas.push_back( QpkCdas[idCdas] );
      }
}

void MakeQpkCdasVsQpkOff(int st, int pmt) {

  TString bnOffl = "~/2021/sdeu/offline/forUb/Aug/offlineUb";
  TString bnUubOffl = "~/2021/sdeu/offline/forUub/AugResults/offlineUub";
  TString bnCdas = "~/2021/sdeu/nouub/underHistos/AugResults/ubChPkPMT";
  TString bnUubCdas = "~/2021/sdeu/underHisto/AugResults/uubChPkPMT";

  TPaveStats *ptstats;
  TLegend *leg;
  TString strFittedEvts;
  TString strSt;
  TString strPmt;
  strSt.Form("St%d", st);
  strPmt.Form("%d", pmt);
  bool ifoff = false;
  bool ifIsUub = false;
  TString whinfo = "chargeVal";

  vector < int > evtIdUbOff; 
  vector < int > evtIdUubOff; 
  vector < int > evtIdUbCdas; 
  vector < int > evtIdUubCdas; 

  vector < double > QpkUbOff; 
  vector < double > QpkUubOff; 
  vector < double > QpkUbCdas; 
  vector < double > QpkUubCdas; 
  
  ifoff = true;
  ifIsUub = false; 
  getEvtIdQpk(bnOffl, strSt, pmt, ifIsUub, ifoff, evtIdUbOff, QpkUbOff);
  ifIsUub = true;
  getEvtIdQpk(bnUubOffl, strSt, pmt, ifIsUub, ifoff, evtIdUubOff, QpkUubOff);

  ifoff = false;
  ifIsUub = false;
  getEvtIdQpk(bnCdas, strSt, pmt, ifIsUub, ifoff, evtIdUbCdas, QpkUbCdas);
  ifIsUub = true;
  getEvtIdQpk(bnUubCdas, strSt, pmt, ifIsUub, ifoff, evtIdUubCdas, QpkUubCdas);

  vector < double > QpkUbMatchIdOff;
  vector < double > QpkUbMatchIdCdas;
  vector < double > QpkUubMatchIdOff;
  vector < double > QpkUubMatchIdCdas;

  doMatch(evtIdUbOff, evtIdUbCdas, QpkUbOff, QpkUbCdas, 
      QpkUbMatchIdOff, QpkUbMatchIdCdas);

  doMatch(evtIdUubOff, evtIdUubCdas, QpkUubOff, QpkUubCdas,
      QpkUubMatchIdOff, QpkUubMatchIdCdas);

  TGraph *grpForUb = new TGraph(QpkUbMatchIdOff.size(), &QpkUbMatchIdOff[0], &QpkUbMatchIdCdas[0]);
  TGraph *grpForUub = new TGraph(QpkUubMatchIdOff.size(), &QpkUubMatchIdOff[0], &QpkUubMatchIdCdas[0]);
/*
  TCanvas *c1 = canvasStyle("c1");
  c1->cd();
  grpForUb->SetTitle("");
  grpForUb->SetMarkerStyle(8);
  grpForUb->GetXaxis()->SetRangeUser(0, 200);
  grpForUb->GetYaxis()->SetRangeUser(0, 200);
  grpForUb->GetYaxis()->SetTitle("Q^{pk}_{VEM} From Local Method [FADC]");
  grpForUb->GetXaxis()->SetTitle("Q^{pk}_{VEM} From SdCalibrator [FADC]");
  grpForUb->Draw("AP");

  leg = new TLegend(0.7, 0.8, 0.9, 0.95);
  leg->SetHeader("UB: "+strSt+" PMT"+strPmt);
  leg->SetTextSize(0.05);
  leg->SetFillStyle(1001);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();
*/
  TCanvas *c2 = canvasStyle("c2");
  c2->cd();
  grpForUub->SetTitle("");
  grpForUub->SetMarkerStyle(8);
  grpForUub->GetXaxis()->SetRangeUser(0, 2000);
  grpForUub->GetYaxis()->SetRangeUser(0, 2000);
  grpForUub->GetYaxis()->SetTitle("Q^{pk}_{VEM} From Local Method [FADC]");
  grpForUub->GetXaxis()->SetTitle("Q^{pk}_{VEM} From SdCalibrator [FADC]");
  grpForUub->Draw("AP");

  leg = new TLegend(0.7, 0.45, 0.8, 0.8);
  leg->SetHeader("UUB: "+strSt+" PMT"+strPmt);
  leg->AddEntry(grpForUub, "","");
  strFittedEvts.Form("%zu", evtIdUubCdas.size());
  leg->AddEntry(grpForUub, "Events fitted by","");
  leg->AddEntry(grpForUub, "Local method: "+strFittedEvts,"");
  strFittedEvts.Form("%zu", evtIdUubOff.size());
  leg->AddEntry(grpForUub, "SdCalibrator: "+strFittedEvts,"");
  leg->AddEntry(grpForUub, "","");
  strFittedEvts.Form("%zu", QpkUubMatchIdOff.size());
  leg->AddEntry(grpForUub, "Events Matched: "+strFittedEvts,"");
  leg->SetTextSize(0.05);
  leg->SetFillStyle(1001);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();
}
