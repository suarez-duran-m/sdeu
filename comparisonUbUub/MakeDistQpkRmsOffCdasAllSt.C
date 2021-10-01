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

TH1F *getFitVals( TString bname, vector < int > stList, int pmt, bool ifIsUub, bool ifoff ) {
  TString monthUub[] = {"Aug"};
  TString pmtId;
  TString strYear[] = {"2019", "2020", "2021"};
  TString fname;
  TString st;
  int stYear = (ifIsUub) ? 2 : 0;
  int lstYear = (ifIsUub) ? 2 : 1;
  // Sept. 1st. for 2019, 2020, and 2021
  int beforeSept [3] = {1251331218, 1282953618, 1314489618};  
  TString strChargeData = (ifoff) ? "charge" : "ChargeData";
  pmtId.Form("%d", pmt);

  TFile *f;
  TTree *chargeInfo;
  double fetchQpkVals = 0.;
  int nBins = (ifIsUub) ? 2000 : 300;
  double stBin = 0.; //(ifIsUub) ? 1000. : 0.;
  double endBin = (ifIsUub) ? 2000. : 300.;
  TString strUborUub = (ifIsUub) ? "Uub" : "Ub";
  TString strOfforCdas = (ifoff) ? "OffLine" : "CDAS";
  cout << "QpkDist"+strUborUub+strOfforCdas << endl;
  TH1F *retQpkDist = new TH1F("QpkDist"+strUborUub+strOfforCdas, "", nBins, stBin, endBin);
  int month = 0;

  int evttime = 0;
  int prevTime = 0;
  int EvtId = 0;

  for ( int year=stYear; year<=lstYear; year++ ) {
    for ( auto & st_i : stList ) {
      st.Form("%d", st_i);
      fname = (ifoff) ? 
        (bname + monthUub[month] + strYear[year] + "St" + st + "Pmt" + pmtId)
        : (bname + pmtId + "St" + st + "lrb35" + monthUub[month] + strYear[year]);
 
      f = TFile::Open(fname+".root");
      chargeInfo = (TTree*)f->Get(strChargeData);
      chargeInfo->SetBranchAddress("chargeVal", &fetchQpkVals); 
      if ( ifoff ) {
        chargeInfo->SetBranchAddress("GpsTime", &evttime);
        chargeInfo->SetBranchAddress("evtId", &EvtId);
      }
      else {
        chargeInfo->SetBranchAddress("timeEvnt", &evttime);
        chargeInfo->SetBranchAddress("eventId", &EvtId);
      }

      prevTime = 0;
      for( int etry=0; etry<chargeInfo->GetEntries(); etry++) {
        chargeInfo->GetEntry(etry);

        if ( fetchQpkVals > 0 && evttime < beforeSept[year] ) { // Events before Sept. 1st, 2021
          // To avoid read the same event twice from Offline
          if ( ifoff ) {
            if ( prevTime != evttime ) {
              retQpkDist->Fill( fetchQpkVals );
              prevTime = evttime;     
            }
          }
          else
            retQpkDist->Fill( fetchQpkVals );
        }
      }
    }
  }

  chargeInfo->Delete();
  f->Close();
  f->Delete();
  TH1F *tmp = (TH1F*)retQpkDist->Clone();
  retQpkDist->Delete();
  return tmp;
}

vector < double > getFailsVals( TString bname, TString st, int pmt, int year, bool ifoff) {
  //TString monthUub[] = {"dec", "jan", "feb", "mar", "apr", "may", "jun", "jul", "aug"};
  // Still without program properly.
  TString monthUub[] = {"Aug"};
  TString pmtId;
  TString fname;
  TString strYear;
  TString strChargeData;
  pmtId.Form("%d", pmt);
  strYear.Form("%d", year);

  TFile *f;
  TTree *chargeInfo;
  double getQpk = 0.;
  vector < double > retQpk;
  int nMonths = 0;

  int evttime = 0;
  int prevTime = 0;

  if ( !ifoff )
    strChargeData = "ChargeData";
  else
    strChargeData = "charge";

  for ( int month=nMonths; month<nMonths+1; month++ ) {
    if ( !ifoff )
      fname = bname + pmtId + st + "lrb35" + monthUub[month] + strYear;
    else if ( ifoff )
      fname = bname + monthUub[month] + strYear + st + "Pmt" + pmtId;

    f = TFile::Open(fname+".root");
    chargeInfo = (TTree*)f->Get(strChargeData);

    getQpk = 0.;
    chargeInfo->SetBranchAddress("chargeVal", &getQpk);
    if ( ifoff && year==2021 )
      chargeInfo->SetBranchAddress("GpsTime", &evttime);

    for( int etry=0; etry<chargeInfo->GetEntries(); etry++) {
      chargeInfo->GetEntry(etry);

      if ( evttime <= 0 || int(getQpk) <= 0 )  
        continue;

      if ( ifoff && year==2021 ) {
        if ( prevTime != evttime ) {
          if ( getQpk == 0 )
            retQpk.push_back( getQpk );
          prevTime = evttime;
        }
      }
      else {
       if ( getQpk == 0 )
         retQpk.push_back( getQpk );
      }
    }
  }
  chargeInfo->Delete();
  f->Delete();
  return retQpk;
}

TGraph *getGraph( vector<int> time, vector<double> values) {
  int nPoints = time.size();
  double xtime[ nPoints ];
  double yvalues[ nPoints ];
  for ( int i=0; i<nPoints; i++ )
  {
    xtime[i] = time[i];
    yvalues[i] = values[i];
  }

  TGraph *grp = new TGraph(nPoints,xtime,yvalues);
  return grp;
}


// =================================
// *** *** *** MAIN CODE *** *** ***

void MakeDistQpkRmsOffCdasAllSt(int pmt) {

  TString stListLoc = "/home/msd/2021/sdeu/fullUubStationsList.txt";
  //TString stListLoc = "/home/msd/2021/sdeu/fullUubStationsListVert.txt";
  TString bnOffl = "~/2021/sdeu/offline/forUb/Aug/offlineUb";
  TString bnUubOffl = "~/2021/sdeu/offline/forUub/AugResults/offlineUub";
  TString bnCdas = "~/2021/sdeu/nouub/underHistos/AugResults/ubChPkPMT";
  TString bnUubCdas = "~/2021/sdeu/underHisto/AugResults/uubChPkPMT";

  TPaveStats *ptstats;
  TLegend *leg;
  TString strMean;
  TString strRms;
  TString strFails;
  TString strGoods;
  TString strPmt;
  strPmt.Form("%d", pmt);
  bool ifoff = false;
  bool ifIsUub = false;
  TString whinfo = "chargeVal";
  double uub2ub = 11.4;

  TH1F *qpkUbOff;  
  TH1F *qpkUubOff; 
  TH1F *qpkUbCdas; 
  TH1F *qpkUubCdas;

  ifstream fileStList;
  double St_i;
  vector < int > stList;
  fileStList.open(stListLoc);
  while( fileStList.good() ) {
    fileStList >> St_i;
    stList.push_back(St_i);
  }
  stList.pop_back();
  fileStList.close();

  TFile *file;
  TH1F *distOffQpkUb;
  TH1F *distOffQpkUub;
  TH1F *distCdasQpkUb;
  TH1F *distCdasQpkUub;

  bool writeFile = true;
  if ( writeFile ) {
    file = new TFile("DistQpkRmsOffCdasAllStPmt"+strPmt+".root", "RECREATE");
    ifIsUub = false;
    ifoff = true;
    distOffQpkUb = getFitVals( bnOffl, stList, pmt, ifIsUub, ifoff);
    file->cd();
    distOffQpkUb->Write();

    ifIsUub = true;
    distOffQpkUub = getFitVals( bnUubOffl, stList, pmt, ifIsUub, ifoff);
    file->cd();
    distOffQpkUub->Write();

    ifIsUub = false;
    ifoff = false;
    distCdasQpkUb = getFitVals( bnCdas, stList, pmt, ifIsUub, ifoff);
    file->cd();
    distCdasQpkUb->Write();
    ifIsUub = true;
    distCdasQpkUub = getFitVals( bnUubCdas, stList, pmt, ifIsUub, ifoff);
    file->cd();
    distCdasQpkUub->Write();

    file->Write();
    file->Close();
  }
  else {    
    file = new TFile("DistQpkRmsOffCdasAllStPmt"+strPmt+".root");
    distOffQpkUb = (TH1F*)file->Get("QpkDistUbOffLine");
    distOffQpkUub = (TH1F*)file->Get("QpkDistUubOffLine");
    distCdasQpkUb = (TH1F*)file->Get("QpkDistUbCDAS");
    distCdasQpkUub = (TH1F*)file->Get("QpkDistUubCDAS");
  }

  // Doing plots
  TCanvas *c1 = canvasStyle("c1");
  c1->cd(0);
  distOffQpkUb->SetStats(kFALSE);
  distOffQpkUb->SetLineColor(kGreen+3);
  distOffQpkUb->SetLineWidth(2);
  distOffQpkUb->SetFillStyle(3001);
  distOffQpkUb->SetFillColor(kGreen+3);
  distOffQpkUb->GetXaxis()->SetTitle("Q^{pk}_{VEM} [FADC]");
  distOffQpkUb->GetYaxis()->SetTitle("Counts [au]");
  distOffQpkUb->GetXaxis()->SetRangeUser(2, 300);
  histoStyle(distOffQpkUb);
  distOffQpkUb->Draw();
  distCdasQpkUb->SetLineColor(kMagenta+2);
  distCdasQpkUb->SetLineWidth(2);
  distCdasQpkUb->SetFillStyle(3001);
  distCdasQpkUb->SetFillColor(kMagenta+2);
  distCdasQpkUb->Draw("same");

  leg = new TLegend(0.55, 0.6, 0.9, 0.95);
  leg->SetHeader("UB: PMT"+strPmt);
  leg->AddEntry(distOffQpkUb, "From SdCalibrator (Q^{pk}_{VEM}/1.01)");
  strMean.Form("%.2f", distOffQpkUb->GetMean());
  strRms.Form("%.4f", distOffQpkUb->GetRMS()/distOffQpkUb->GetMean());
  leg->AddEntry(distOffQpkUb, "#mu="+strMean+"; RMS/#mu="+strRms,"");
  leg->AddEntry(distOffQpkUb, "", "");
  leg->AddEntry(distCdasQpkUb, "Own Implementation");
  strMean.Form("%.2f", distCdasQpkUb->GetMean());
  strRms.Form("%.4f", distCdasQpkUb->GetRMS()/distCdasQpkUb->GetMean()); 
  leg->AddEntry(distCdasQpkUb, "#mu="+strMean+"; RMS/#mu="+strRms,"");
  leg->SetTextSize(0.045);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();
  c1->Print("../plots/QpkDistributionUbPmt"+strPmt+".pdf");

  TCanvas *c2 = canvasStyle("c2");
  c2->cd(0);
  distOffQpkUub->SetStats(kFALSE);
  distOffQpkUub->SetLineColor(kGreen+3);
  distOffQpkUub->SetLineWidth(2);
  distOffQpkUub->SetFillStyle(3001);
  distOffQpkUub->SetFillColor(kGreen+3);
  distOffQpkUub->GetXaxis()->SetTitle("Q^{pk}_{VEM} [FADC]");
  distOffQpkUub->GetYaxis()->SetTitle("Counts [au]");
  distOffQpkUub->GetXaxis()->SetRangeUser(5, 2000);
  histoStyle(distOffQpkUub);
  distOffQpkUub->Draw();
  distCdasQpkUub->SetLineColor(kMagenta+2);
  distCdasQpkUub->SetLineWidth(2);
  distCdasQpkUub->SetFillStyle(3001);
  distCdasQpkUub->SetFillColor(kMagenta+2);
  distCdasQpkUub->Draw("same");

  leg = new TLegend(0.15, 0.6, 0.5, 0.95);
  leg->SetHeader("UUB: PMT"+strPmt);
  leg->AddEntry(distOffQpkUub, "From SdCalibrator (Q^{pk}_{VEM}/1.01)");
  strMean.Form("%.2f", distOffQpkUub->GetMean());
  strRms.Form("%.4f", distOffQpkUub->GetRMS()/distOffQpkUub->GetMean());
  leg->AddEntry(distOffQpkUub, "#mu="+strMean+"; RMS/#mu="+strRms,"");
  leg->AddEntry(distOffQpkUub, "", "");
  leg->AddEntry(distCdasQpkUub, "Own Implementation");
  strMean.Form("%.2f", distCdasQpkUub->GetMean());
  strRms.Form("%.4f", distCdasQpkUub->GetRMS()/distCdasQpkUub->GetMean());
  leg->AddEntry(distCdasQpkUub, "#mu="+strMean+"; RMS/#mu="+strRms,"");
  leg->SetTextSize(0.045);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();
  c2->Print("../plots/QpkDistributionUubPmt"+strPmt+".pdf");
 
  exit(0);
}
