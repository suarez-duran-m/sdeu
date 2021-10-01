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

double getmean( TH1F *arr ) {
  double mean = 0.;
  int nb = arr->GetXaxis()->GetNbins();
  for (int bin_i=0; bin_i<nb; bin_i++)
    mean += arr->GetBinCenter(bin_i)*arr->GetBinContent(bin_i);

  return mean/arr->GetEntries();
}

double getrms( TH1F *arr, double meanarr ) {
  double rms = 0.;
  int nb = arr->GetXaxis()->GetNbins();
  double tmpVal = 0.;
  for (int bin_i=0; bin_i<nb; bin_i++) {
    tmpVal = arr->GetBinLowEdge(bin_i);
    rms += (tmpVal-meanarr)*(tmpVal-meanarr)*arr->GetBinContent(bin_i);
  }

  return sqrt(rms/arr->GetEntries());
}

void getEvtIdQpk( TString bname, vector<int> listSt, int pmt, 
    bool ifIsUub, bool ifoff, vector<int> &retId, 
    vector<double> &retQpk, vector<int> &retStId, ofstream &outFile) {

  TString monthUub[] = {"Aug"};
  TString pmtId;
  TString strStId;
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
    for ( auto & stId_i : listSt ) {
      strStId.Form("St%d", stId_i);
      fname = (ifoff) ? 
        (bname + monthUub[month] + strYear[year] + strStId + "Pmt" + pmtId)
        : (bname + pmtId + strStId + "lrb35" + monthUub[month] + strYear[year]);
  
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
              retStId.push_back( stId_i );
              prevTime = evttime;
              //outFile << "Offline " << ifIsUub << " " << EvtId << " " << stId_i << " " << getQpkVals << "\n";
            }
          }
          else {
            retId.push_back( EvtId );
            retQpk.push_back( getQpkVals );
            retStId.push_back( stId_i );
            //outFile << "CDAS " << ifIsUub << " " << EvtId << " " << stId_i << " " << getQpkVals << "\n";
          }
        }
      }
    }
  }

  chargeInfo->Delete();
  f->Delete();
}

void doMatch(vector<int> evtIdOff, vector<int> evtIdCdas, 
    vector<double> QpkOff, vector<double> QpkCdas,
    vector<int> stIdOff, vector<int> stIdCdas,
    vector<double> &matchOff, vector<double> &matchCdas, 
    vector<int> &matchId, vector<int> &matchSt) {

  for ( int idOff=0; idOff<evtIdOff.size(); idOff++ )
    for ( int idCdas=0; idCdas<evtIdCdas.size(); idCdas++ )
      if ( evtIdOff[idOff] == evtIdCdas[idCdas] )
        if ( stIdOff[idOff] == stIdCdas[idCdas] ) {
          matchOff.push_back( QpkOff[idOff] );
          matchCdas.push_back( QpkCdas[idCdas] );
          matchId.push_back( evtIdOff[idOff] );
          matchSt.push_back( stIdOff[idOff] );
          if ( QpkCdas[idCdas] < 500 && QpkCdas[idCdas] > 0 && QpkOff[idOff] == 0 )
            cout << "EvtId " << evtIdOff[idOff] 
              << " StId " << stIdOff[idOff] 
              << " QpkCdas " << QpkCdas[idCdas] 
              << " QpkOff " << QpkOff[idOff] << endl;
        }
}

void QpkCdasVsQpkOffAllSt(int pmt) {

  TString stListLoc = "/home/msd/2021/sdeu/fullUubStationsListVert.txt";
  TString bnOffl = "~/2021/sdeu/offline/forUb/Aug/offlineUb";
  TString bnUubOffl = "~/2021/sdeu/offline/forUub/AugResults/offlineUub";
  TString bnCdas = "~/2021/sdeu/nouub/underHistos/AugResults/ubChPkPMT";
  TString bnUubCdas = "~/2021/sdeu/underHisto/AugResults/uubChPkPMT";

  TPaveStats *ptstats;
  TLegend *leg;
  TString strEntr;
  TString strMean;
  TString strRms;
  TString strPmt;
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

  vector < int > stIdUbOff;
  vector < int > stIdUubOff; 
  vector < int > stIdUbCdas;
  vector < int > stIdUubCdas;

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

  ofstream outFile;
  outFile.open ("tableQpkValues.dat");
  
  ifoff = true;
  ifIsUub = false;
  getEvtIdQpk(bnOffl, stList, pmt, ifIsUub, ifoff, evtIdUbOff, QpkUbOff, 
      stIdUbOff, outFile);
  ifIsUub = true;
  getEvtIdQpk(bnUubOffl, stList, pmt, ifIsUub, ifoff, evtIdUubOff, QpkUubOff, 
      stIdUubOff, outFile);

  ifoff = false;
  ifIsUub = false;
  getEvtIdQpk(bnCdas, stList, pmt, ifIsUub, ifoff, evtIdUbCdas, QpkUbCdas, 
      stIdUbCdas, outFile);
  ifIsUub = true;
  getEvtIdQpk(bnUubCdas, stList, pmt, ifIsUub, ifoff, evtIdUubCdas, QpkUubCdas, 
      stIdUubCdas, outFile);
  outFile.close();

  vector < double > QpkUbMatchIdOff;
  vector < double > QpkUbMatchIdCdas;
  vector < double > QpkUubMatchIdOff;
  vector < double > QpkUubMatchIdCdas;

  vector < int > matchUb_evtId;
  vector < int > matchUub_evtId;
  vector < int > matchUb_stId;
  vector < int > matchUub_stId;

  doMatch(evtIdUbOff, evtIdUbCdas, QpkUbOff, QpkUbCdas,
      stIdUbOff, stIdUbCdas, QpkUbMatchIdOff, QpkUbMatchIdCdas, 
      matchUb_evtId, matchUb_stId);
  doMatch(evtIdUubOff, evtIdUubCdas, QpkUubOff, QpkUubCdas,
      stIdUubOff, stIdUubCdas, QpkUubMatchIdOff, QpkUubMatchIdCdas, 
      matchUub_evtId, matchUub_stId);

  TGraph *grpForUb = new TGraph(
      QpkUbMatchIdOff.size(), &QpkUbMatchIdOff[0], &QpkUbMatchIdCdas[0]);
  TGraph *grpForUub = new TGraph(
      QpkUubMatchIdOff.size(), &QpkUubMatchIdOff[0], &QpkUubMatchIdCdas[0]);

  TCanvas *c1 = canvasStyle("c1");
  c1->cd();
  grpForUb->SetTitle("");
  grpForUb->SetMarkerStyle(8);
  grpForUb->GetXaxis()->SetRangeUser(0, 300);
  grpForUb->GetYaxis()->SetRangeUser(0, 300);
  grpForUb->GetYaxis()->SetTitle("Q^{pk}_{VEM} From Local Method [FADC]");
  grpForUb->GetXaxis()->SetTitle("Q^{pk}_{VEM}/1.01 From SdCalibrator [FADC]");
  grpForUb->SetMarkerColor(kGray+3);
  grpForUb->Draw("AP");

  leg = new TLegend(0.7, 0.8, 0.9, 0.95);
  leg->SetHeader("UB: PMT"+strPmt);
  leg->SetTextSize(0.05);
  leg->SetFillStyle(1001);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();
  c1->Print("../plots/ubQpkCdasVsQpkOffAllStPMT"+strPmt+".pdf");

  TCanvas *c2 = canvasStyle("c2");
  c2->cd();
  grpForUub->SetTitle("");
  grpForUub->SetMarkerStyle(8);
  grpForUub->GetXaxis()->SetRangeUser(0, 3000);
  grpForUub->GetYaxis()->SetRangeUser(0, 3000);
  grpForUub->GetYaxis()->SetTitle("Q^{pk}_{VEM} From Local Method [FADC]");
  grpForUub->GetXaxis()->SetTitle("Q^{pk}_{VEM}/1.01 From SdCalibrator [FADC]");
  grpForUub->SetMarkerColor(kGray+3);
  grpForUub->Draw("AP");

  leg = new TLegend(0.7, 0.8, 0.9, 0.95);
  leg->SetHeader("UUB: PMT"+strPmt);
  leg->SetTextSize(0.05);
  leg->SetFillStyle(1001);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();
  c2->Print("../plots/uubQpkCdasVsQpkOffAllStPMT"+strPmt+".png");//.pdf");


  // Doing Relative differences
  TH1F *relDiffQpkUbMatch = new TH1F("relDiffQpkUbMatch","", 40000, -2., 2.);
  for (int id_match=0; id_match<QpkUbMatchIdOff.size(); id_match++)
    if ( QpkUbMatchIdCdas[id_match] > 100 && QpkUbMatchIdOff[id_match] > 100 )
      relDiffQpkUbMatch->Fill(QpkUbMatchIdCdas[id_match]/(1.01*QpkUbMatchIdOff[id_match]) -1.); 
  
  TH1F *relDiffQpkUubMatch = new TH1F("relDiffQpkUubMatch","", 40000, -2., 2.); 
  for (int id_match=0; id_match<QpkUubMatchIdOff.size(); id_match++)
    if ( QpkUubMatchIdCdas[id_match] > 1000 && QpkUubMatchIdOff[id_match] > 1000 )
      relDiffQpkUubMatch->Fill(QpkUubMatchIdCdas[id_match]/(1.01*QpkUubMatchIdOff[id_match]) -1.);

  double xlim[2] = {-0.2, 0.2};

  TCanvas *c3 = canvasStyle("c3");
  c3->cd();

  //gPad->SetLogy();
  relDiffQpkUbMatch->SetStats(kFALSE);
  relDiffQpkUbMatch->GetXaxis()->SetTitle("Relative Difference Q^{pk}_{VEM} (Local/SdCalibrator - 1) [au]");
  relDiffQpkUbMatch->GetYaxis()->SetTitle("Counts [au]");
  relDiffQpkUbMatch->GetXaxis()->SetRangeUser(xlim[0], xlim[1]);
  relDiffQpkUbMatch->GetXaxis()->SetTitleOffset(1.3);
  relDiffQpkUbMatch->GetYaxis()->SetTitleOffset(0.8);
  relDiffQpkUbMatch->Draw();

  leg = new TLegend(0.2, 0.7, 0.4, 0.95);
  leg->SetHeader("UB: PMT"+strPmt);
  strEntr.Form("%d", (int)relDiffQpkUbMatch->GetEntries());
  double tmpMean = 0.;
  tmpMean = getmean(relDiffQpkUbMatch);
  strMean.Form("%.5f", tmpMean);
  strRms.Form("%.5f", getrms(relDiffQpkUbMatch, tmpMean));
  leg->AddEntry(relDiffQpkUbMatch, "Entries: "+strEntr,"");
  leg->AddEntry(relDiffQpkUbMatch, "Mean: "+strMean,"");
  leg->AddEntry(relDiffQpkUbMatch, "RMS: "+strRms,"");
  leg->SetTextSize(0.04);
  leg->SetFillStyle(1001);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();
  c3->Print("../plots/ubRelDiffQpkUbMatchAllStPmt"+strPmt+".pdf");

  TCanvas *c4 = canvasStyle("c4");
  c4->cd();

  //gPad->SetLogy();
  relDiffQpkUubMatch->SetStats(kFALSE);
  relDiffQpkUubMatch->GetXaxis()->SetTitle("Relative Difference Q^{pk}_{VEM} (Local/SdCalibrator - 1) [au]");
  relDiffQpkUubMatch->GetYaxis()->SetTitle("Counts [au]");
  relDiffQpkUubMatch->GetXaxis()->SetRangeUser(xlim[0], xlim[1]);
  relDiffQpkUubMatch->GetXaxis()->SetTitleOffset(1.3);
  relDiffQpkUubMatch->GetYaxis()->SetTitleOffset(1.1);
  relDiffQpkUubMatch->Draw();

  leg = new TLegend(0.2, 0.7, 0.4, 0.95);
  leg->SetHeader("UUB: PMT"+strPmt);
  strEntr.Form("%d", (int)relDiffQpkUubMatch->GetEntries());
  tmpMean = getmean(relDiffQpkUubMatch);
  strMean.Form("%.5f", tmpMean);
  strRms.Form("%.5f", getrms(relDiffQpkUubMatch, tmpMean));
  leg->AddEntry(relDiffQpkUubMatch, "Entries: "+strEntr,"");
  leg->AddEntry(relDiffQpkUubMatch, "Mean: "+strMean,"");
  leg->AddEntry(relDiffQpkUubMatch, "RMS: "+strRms,"");
  leg->SetTextSize(0.04);
  leg->SetFillStyle(1001);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();
  c4->Print("../plots/uubRelDiffQpkUbMatchAllStPmt"+strPmt+".pdf");
  
  exit(0);
} 
