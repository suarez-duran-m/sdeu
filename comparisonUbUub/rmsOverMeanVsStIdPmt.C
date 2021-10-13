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
  int nb = arr.size();
  int goodVals = 0;
  double mean = 0.;
  for (int bin_i=0; bin_i<nb; bin_i++)
    if ( arr[bin_i] > 0 ) {
      mean += arr[bin_i];
      goodVals++;
    }
  return mean/goodVals;
}

double getrms( vector<double> arr, double meanarr ) {
  int nb = arr.size();
  double rms = 0.;
  int goodVals = 0;
  for (int bin_i=0; bin_i<nb; bin_i++)
    if ( arr[bin_i]> 0 ) {
      rms += (arr[bin_i]-meanarr)*(arr[bin_i]-meanarr);
      goodVals++;
    }

  return sqrt(rms/goodVals)/meanarr;
}

vector < vector < double > > doRmsRelVsStId( TString bname, vector<double> listSt, int pmt, 
    bool ifIsUub ) {

  TString monthUub[2] = {"Aug", "Sep"};
  TString pmtId;
  TString strStId;
  TString strYear[3] = {"2019", "2020", "2021"};
  TString fname;
  int stYear = (ifIsUub) ? 2 : 0;
  int lstYear = (ifIsUub) ? 2 : 1;

  TString strChargeData = "ChargeData";
  pmtId.Form("%d", pmt);

  TFile *f;
  TTree *chargeInfo;
  double fetchQpkVals = 0.;

  vector < vector < double > > qpkStDist;
  qpkStDist.resize( listSt.size() );

  for ( int stId_i=0; stId_i<listSt.size(); stId_i++ ) {
    for ( int year=stYear; year<=lstYear; year++ )
      for ( int month_i=0; month_i<sizeof(monthUub)/sizeof(*monthUub); month_i++ ) {
        strStId.Form("St%d", (int)listSt[stId_i]);
        fname = bname + pmtId + strStId + "lrb35" + monthUub[month_i] + strYear[year];
  
        f = TFile::Open(fname+".root");
        chargeInfo = (TTree*)f->Get(strChargeData);
        chargeInfo->SetBranchAddress("chargeVal", &fetchQpkVals);

        fetchQpkVals = 0.;
        for( int etry=0; etry<chargeInfo->GetEntries(); etry++) {
          chargeInfo->GetEntry(etry);
          if ( fetchQpkVals == 0 || fetchQpkVals < 0 )
            continue;
          qpkStDist[stId_i].push_back( fetchQpkVals );
        }
      }
  }
  chargeInfo->Delete();
  f->Delete();
  return qpkStDist;
}

void rmsOverMeanVsStIdPmt(int pmt) {

  string stListPath = "/home/msd/2021/sdeu/listStationsShort.txt";
  //string stListPath = "/home/msd/2021/sdeu/fullUubStationsListVert.txt";
  TString bnOffl = "~/2021/sdeu/offline/forUb/results/offlineUb";
  TString bnUubOffl = "~/2021/sdeu/offline/forUub/results/offlineUub";
  TString bnCdas = "~/2021/sdeu/nouub/underHistos/results/ubChPkPMT";
  TString bnUubCdas = "~/2021/sdeu/underHisto/results/uubChPkPMT";

  TPaveStats *ptstats;
  TLegend *leg;
  TString strEntr;
  TString strMean;
  TString strRms;
  TString strPmt;
  strPmt.Form("%d", pmt);
  bool ifIsUub = false;

  ifstream fileStList;
  double St_i;
  vector < double > stListId;
  vector< double > stList;
  fileStList.open(stListPath);
  int cnt = 1;
  while( fileStList.good() ) {
    fileStList >> St_i;
    stListId.push_back(St_i);
    stList.push_back(cnt);
    cnt++;
  }
  stListId.pop_back();
  stList.pop_back();
  fileStList.close();

  int nStations = stList.size();
  vector < vector < double > > distQpkStUb;
  vector < vector < double > > distQpkStUub;

  ifIsUub = false;
  distQpkStUb = doRmsRelVsStId(bnCdas, stListId, pmt, ifIsUub);

  ifIsUub = true;
  distQpkStUub = doRmsRelVsStId(bnUubCdas, stListId, pmt, ifIsUub);

  vector < double > rmsMeanUb;
  vector < double > rmsRmsMeanUb;
  vector < double > rmsMeanUub;
  vector < double > rmsRmsMeanUub;
  rmsMeanUb.resize(stList.size());
  rmsRmsMeanUb.resize(stList.size());
  rmsMeanUub.resize(stList.size());
  rmsRmsMeanUub.resize(stList.size());

  TH1F *tmpDistUb;
  TH1F *tmpDistUub;
  vector < double > meanRmsRelUb;
  vector < double > meanRmsRelUub;
  for ( int st_i=0; st_i<stList.size(); st_i++ ) {
    tmpDistUb = new TH1F("tmpDistUb", "", 250, 50, 300);
    tmpDistUub = new TH1F("tmpDistUub", "", 2500, 500, 3000);
    for ( auto & qpk_i : distQpkStUb[st_i] )
      tmpDistUb->Fill( qpk_i );
    for ( auto & qpk_i : distQpkStUub[st_i] )
      tmpDistUub->Fill( qpk_i );
    rmsMeanUb[st_i] = 100.*tmpDistUb->GetRMS()/tmpDistUb->GetMean();
    rmsRmsMeanUb[st_i] = tmpDistUb->GetRMSError(1);
    rmsMeanUub[st_i] = 100.*tmpDistUub->GetRMS()/tmpDistUub->GetMean();
    rmsRmsMeanUub[st_i] = tmpDistUub->GetRMSError(1);

    meanRmsRelUb.push_back( rmsMeanUb[st_i] );
    meanRmsRelUub.push_back( rmsMeanUub[st_i] );    
    tmpDistUb->Clear();
    tmpDistUub->Clear();
  }

  TString strStName;
  Int_t binForSt;

  TH1F *grpRmsRelUbCdas = new TH1F("grpRmsRelUbCdas", "", nStations, &stList[0]);
  TH1F *grpRmsRelUubCdas = new TH1F("grpRmsRelUubCdas", "", nStations, &stList[0]);

  for ( int i=0; i<stList.size(); i++ ) {
    strStName.Form("%d", (int)stListId[i]);
    if ( rmsMeanUb[i] > 0 ) {
      grpRmsRelUbCdas->SetBinContent( stList[i], rmsMeanUb[i] );
      grpRmsRelUbCdas->SetBinError( stList[i], rmsRmsMeanUb[i] );
      grpRmsRelUbCdas->GetXaxis()->SetBinLabel( stList[i], strStName);
    }
    if ( rmsMeanUub[i] > 0 ) {
      grpRmsRelUubCdas->SetBinContent( stList[i], rmsMeanUub[i] );
      grpRmsRelUubCdas->SetBinError( stList[i], rmsRmsMeanUub[i] );
      grpRmsRelUubCdas->GetXaxis()->SetBinLabel( stList[i], strStName);
    }
  }

  /*
  TF1 *poly0 = new TF1("poly0","[0]", stList[0], stList[nStations]);
  grpRmsRelUbCdas->Fit("poly0","QR");
  grpRmsRelUubCdas->Fit("poly0","QR");
  grpRmsRelUbCdas->GetFunction("poly0")->SetLineColor(kBlack);
  grpRmsRelUubCdas->GetFunction("poly0")->SetLineColor(kRed);
  */

  TLine *meanUb = new TLine(stList[1], getmean(meanRmsRelUb), stList.size(), getmean(meanRmsRelUb));
  TLine *meanUub = new TLine(stList[1], getmean(meanRmsRelUub), stList.size(), getmean(meanRmsRelUub));
  
  TCanvas *c1 = canvasStyle("c1");
  c1->cd();
  
  grpRmsRelUubCdas->SetStats(kFALSE);
  grpRmsRelUubCdas->SetTitle("");
  grpRmsRelUubCdas->GetYaxis()->SetTitle("#LT RMS/#LT Q^{pk}_{VEM} #GT #GT [%]");
  grpRmsRelUubCdas->GetXaxis()->SetTitle("Station Id.");
  grpRmsRelUubCdas->GetXaxis()->SetTitleOffset(1.5);
  grpRmsRelUubCdas->SetMarkerStyle(8);
  grpRmsRelUubCdas->SetMarkerColor(kRed);
  grpRmsRelUubCdas->SetMarkerSize(1.5);
  grpRmsRelUubCdas->SetLineColor(kRed);
  grpRmsRelUubCdas->Draw("E1");

  grpRmsRelUbCdas->SetMarkerColor(kBlack);
  grpRmsRelUbCdas->SetMarkerStyle(22);
  grpRmsRelUbCdas->SetMarkerSize(1.5);
  grpRmsRelUbCdas->SetLineColor(kBlack);
  grpRmsRelUbCdas->Draw("E1 same");

  meanUb->SetLineColor(kBlack);
  meanUb->SetLineWidth(2);
  meanUb->Draw();
  meanUub->SetLineColor(kRed);
  meanUub->SetLineWidth(2);
  meanUub->Draw();

  leg = new TLegend(0.15,0.65,0.4,0.95);
  leg->SetHeader("PMT"+strPmt);
  leg->AddEntry(grpRmsRelUbCdas, "UB stations", "p");
  TString strAve;
  strAve.Form("%.3f", getmean(meanRmsRelUb));
  leg->AddEntry(meanUb, "Average: "+strAve, "l" );
  strAve.Form("%.3f", getrms(meanRmsRelUb, getmean(meanRmsRelUb)));
  leg->AddEntry(meanUb, "RMS: "+strAve, "");
  leg->AddEntry(grpRmsRelUubCdas, "UUB stations", "p");
  strAve.Form("%.3f", getmean(meanRmsRelUub));
  leg->AddEntry(meanUub, "Average: "+strAve, "l" );
  strAve.Form("%.3f", getrms(meanRmsRelUub, getmean(meanRmsRelUb)));
  leg->AddEntry(meanUub, "RMS: "+strAve, "");
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();
  c1->Print("../plots/accuracyQpksFitsCdasPmt"+strPmt+".pdf");


  tmpDistUub = new TH1F("tmpDistUub", "", 250, 50, 300);
  //tmpDistUub = new TH1F("tmpDistUub", "", 2500, 500, 3000);
  double rmsErrHand = 0.;
  double tmpMean =0.;
  double tmpRms =0.;
  double rmsRel = 0.;
  for ( int st_i=0; st_i<stList.size(); st_i++ ) {
    if ( stListId[st_i] != 862 )
      continue;
    for ( auto & qpk_i : distQpkStUb[st_i] )
      tmpDistUub->Fill( qpk_i );
    tmpMean = getmean(distQpkStUb[st_i]);
    rmsRel = getrms(distQpkStUb[st_i], tmpMean);
    tmpRms = tmpMean*rmsRel;
    rmsErrHand = tmpRms/sqrt(distQpkStUb[st_i].size());
  }

  TCanvas *c2 = canvasStyle("c2");
  c2->cd();
  tmpDistUub->SetStats(kFALSE);
  tmpDistUub->Fit("gaus");
  tmpDistUub->GetXaxis()->SetRangeUser(1450, 1850);
  tmpDistUub->Draw();

  leg = new TLegend(0.13,0.4,0.4,0.95);
  leg->SetHeader("UUB: Station 1211, PMT"+strPmt);
  TString strRMS;
  strRMS.Form("%.3f", tmpDistUub->GetMean());
  leg->AddEntry(tmpDistUub, "Mean-ROOT: "+strRMS,"");
  strRMS.Form("%.3f", tmpDistUub->GetRMS());
  leg->AddEntry(tmpDistUub, "RMS-ROOT: "+strRMS,"");
  strRMS.Form("%.3f", 100.*tmpDistUub->GetRMS()/tmpDistUub->GetMean());
  leg->AddEntry(tmpDistUub, "100*(RMS/Mean): "+strRMS, "");
  strRMS.Form("%.3f", 100.*rmsRel);
  leg->AddEntry(tmpDistUub, "100*(RMS/Mean), by hand: "+strRMS, "");
  strRMS.Form("%.3f", tmpDistUub->GetRMSError(1));
  leg->AddEntry(tmpDistUub, "ErrRMS-ROOT: "+strRMS, "");
  strRMS.Form("%.3f", rmsErrHand);
  leg->AddEntry(tmpDistUub, "ErrRMS-Hand: "+strRMS, "");
  leg->AddEntry(tmpDistUub, "", "");
  strRMS.Form("%.3f", tmpDistUub->GetFunction("gaus")->GetParameter(1));
  leg->AddEntry(tmpDistUub, "Guas-Mean: "+strRMS,"");
  strRMS.Form("%.3f", tmpDistUub->GetFunction("gaus")->GetParameter(2));
  leg->AddEntry(tmpDistUub, "Guas-Sigma: "+strRMS,"");
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();

  leg = new TLegend(0.65,0.75,0.9,0.95);
  leg->AddEntry(tmpDistUub, "ErrRMS-Hand:", "");
  leg->AddEntry(tmpDistUub, "#sqrt{ VAR/n }", "");
  leg->AddEntry(tmpDistUub, "ErrRMS-ROOT:", "");
  leg->AddEntry(tmpDistUub, "#sqrt{ VAR/2n }", "");
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();
  //c2->Print("qpksDistriSt1211Pmt1.pdf");
  
  //exit(0);
} 
