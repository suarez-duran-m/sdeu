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

void doRmsRelVsStId( TString bname, vector<double> listSt, int pmt, 
    bool ifIsUub, bool ifoff, vector<double> &retMeanRmsRel5days, 
    vector<double> &retRmsQpksTot ) {

  TString monthUub[2] = {"Aug", "Sep"};
  TString pmtId;
  TString strStId;
  TString strYear[3] = {"2019", "2020", "2021"};
  TString fname;
  int stYear = (ifIsUub) ? 2 : 0;
  int lstYear = (ifIsUub) ? 2 : 1;
  // Oct. 1st. for 2019, 2020, and 2021
  int beforeOct [3] = {1253923218, 1285545618, 1317081618};
  // Sep. 1st. for 2019, 2020, and 2021
  int aug1st [3] = {1248652818, 1280275218, 1311811218};
  TString strChargeData = (ifoff) ? "charge" : "ChargeData";
  pmtId.Form("%d", pmt);

  TFile *f;
  TTree *chargeInfo;
  double fetchQpkVals = 0.;

  int prevTime = 0;
  int fetchEvttime = 0;
  int fetchEvtId = 0;

  vector< vector < double > > qpksIn5days;
  int nWeeks = 13;//7; //13; // 5-days Weeks in August-Sept.
  qpksIn5days.resize(nWeeks); 
  int currentWeek = 0;
  int nsec5days = 432000;//864000; //432000;
  vector < double > meanQpksIn5days;
  vector < double > rmsQpkIn5days;

  for ( int stId_i=0; stId_i<listSt.size(); stId_i++ ) {
    for ( int year=stYear; year<=lstYear; year++ ) {
      currentWeek = 0;
      meanQpksIn5days.resize(nWeeks);
      rmsQpkIn5days.resize(nWeeks);
      for ( int month_i=0; month_i<sizeof(monthUub)/sizeof(*monthUub); month_i++ ) {
        strStId.Form("St%d", (int)listSt[stId_i]);
        fname = (ifoff) ?
          (bname + monthUub[month_i] + strYear[year] + strStId + "Pmt" + pmtId)
          : (bname + pmtId + strStId + "lrb35" + monthUub[month_i] + strYear[year]);
  
        f = TFile::Open(fname+".root");
        chargeInfo = (TTree*)f->Get(strChargeData);
        chargeInfo->SetBranchAddress("chargeVal", &fetchQpkVals);

        if ( ifoff )
          chargeInfo->SetBranchAddress("GpsTime", &fetchEvttime);
        else 
          chargeInfo->SetBranchAddress("timeEvnt", &fetchEvttime);

        prevTime = 0;
        fetchQpkVals = 0.;
        for( int etry=0; etry<chargeInfo->GetEntries(); etry++) {
          chargeInfo->GetEntry(etry);
          if ( fetchQpkVals == 0 || fetchQpkVals < 0 )
            continue;
          currentWeek = (fetchEvttime-aug1st[year])/nsec5days;
          if ( fetchEvttime < beforeOct[year] ) { // Events before Sept. 1st, 2021
            // To avoid read the same event twice from Offline
            if ( ifoff ) {
              if ( prevTime != fetchEvttime ) {
                qpksIn5days[currentWeek].push_back( fetchQpkVals );
                prevTime = fetchEvttime;
              }
            }
            else
              qpksIn5days[currentWeek].push_back( fetchQpkVals );
          }
        }
      }
      double tmpMean = 0.;
      for ( int week=0; week<nWeeks; week++ ) {
        tmpMean = getmean(qpksIn5days[week]);
        meanQpksIn5days[week] = (tmpMean > 0) ? tmpMean : 0.;
        // 100.*getrms to express RMS/mu in % units
        rmsQpkIn5days[week] = (tmpMean > 0) ? 100.*getrms(qpksIn5days[week], tmpMean) : 0;  
        qpksIn5days[week].clear();
      }
      tmpMean = getmean( rmsQpkIn5days );
      if ( strYear[year]=="2020" ) {
        retMeanRmsRel5days[stId_i] += tmpMean;
        retMeanRmsRel5days[stId_i] *= 0.5;
        retRmsQpksTot[stId_i] += getrms(rmsQpkIn5days, tmpMean)*tmpMean;
        retRmsQpksTot[stId_i] *= 0.5;       
      }
      else {
        retMeanRmsRel5days[stId_i] = tmpMean;
        retRmsQpksTot[stId_i] = getrms(rmsQpkIn5days, tmpMean)*tmpMean;
      }

      meanQpksIn5days.clear();
      rmsQpkIn5days.clear();
    }
  }
  chargeInfo->Delete();
  f->Delete();
}

void rmsOverMeanVsStId(int pmt) {

  //string stListPath = "/home/msd/2021/sdeu/listStations9.txt";
  string stListPath = "/home/msd/2021/sdeu/fullUubStationsListVert.txt";
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
  bool ifoff = false;
  bool ifIsUub = false;
  TString whinfo = "chargeVal";

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
  vector < double > meanRmsRelUbCdas;
  vector < double > rmsRmsRelUbCdas;
  vector < double > meanRmsRelUubCdas;
  vector < double > rmsRmsRelUubCdas;

  meanRmsRelUbCdas.resize( nStations );
  rmsRmsRelUbCdas.resize( nStations );
  meanRmsRelUubCdas.resize( nStations );
  rmsRmsRelUubCdas.resize( nStations );
/*
  vector < double > meanRmsRelUbOff;
  vector < double > rmsRmsRelUbOff;
  vector < double > meanRmsRelUubOff;
  vector < double > rmsRmsRelUubOff;

  meanRmsRelUbOff.resize( nStations );
  rmsRmsRelUbOff.resize( nStations );
  meanRmsRelUubOff.resize( nStations );
  rmsRmsRelUubOff.resize( nStations );
  
  ifoff = true;
  ifIsUub = false;
  doRmsRelVsStId(bnOffl, stListId, pmt, ifIsUub, ifoff, 
      meanRmsRelUbOff, rmsRmsRelUbOff);

  ifIsUub = true;
  doRmsRelVsStId(bnUubOffl, stListId, pmt, ifIsUub, ifoff, 
      meanRmsRelUubOff, rmsRmsRelUubOff);  
  */

  ifoff = false;
  ifIsUub = false;
  doRmsRelVsStId(bnCdas, stListId, pmt, ifIsUub, ifoff,
      meanRmsRelUbCdas, rmsRmsRelUbCdas);

  ifoff = false;
  ifIsUub = true;
  doRmsRelVsStId(bnUubCdas, stListId, pmt, ifIsUub, ifoff
      ,meanRmsRelUubCdas, rmsRmsRelUubCdas);

  TString strStName;
  Int_t binForSt;

  //TH1F *grpRmsRelUbOff = new TH1F("grpRmsRelUbOff", "", nStations, &stList[0]);
  //TH1F *grpRmsRelUubOff = new TH1F("grpRmsRelUubOff", "", nStations, &stList[0]);
  TH1F *grpRmsRelUbCdas = new TH1F("grpRmsRelUbCdas", "", nStations, &stList[0]);
  TH1F *grpRmsRelUubCdas = new TH1F("grpRmsRelUubCdas", "", nStations, &stList[0]);

  for ( int i=0; i<stList.size(); i++ ) {
/*
    if ( meanRmsRelUbOff[i] > 0 )
      grpRmsRelUbOff->SetBinContent( stList[i], meanRmsRelUbOff[i] );
    else
      grpRmsRelUbOff->SetBinContent( stList[i], -0.1 );
*/
    /*
    if ( meanRmsRelUubOff[i] > 0 )
      grpRmsRelUubOff->SetBinContent( stList[i], meanRmsRelUubOff[i] ); 
    else
      grpRmsRelUubOff->SetBinContent( stList[i], -0.1 );
*/
    strStName.Form("%d", (int)stListId[i]);
    if ( meanRmsRelUbCdas[i] > 0 ) {
      grpRmsRelUbCdas->SetBinContent( stList[i], meanRmsRelUbCdas[i] );
      grpRmsRelUbCdas->SetBinError( stList[i], rmsRmsRelUbCdas[i] );
      grpRmsRelUbCdas->GetXaxis()->SetBinLabel( stList[i], strStName);
    }
    if ( meanRmsRelUubCdas[i] > 0 ) {
      grpRmsRelUubCdas->SetBinContent( stList[i], meanRmsRelUubCdas[i] ); 
      grpRmsRelUubCdas->SetBinError( stList[i], rmsRmsRelUubCdas[i] );
      grpRmsRelUubCdas->GetXaxis()->SetBinLabel( stList[i], strStName);
    }

    //grpRmsRelUubOff->SetBinError( stList[i], rmsRmsRelUubOff[i] );
    //grpRmsRelUubOff->GetXaxis()->SetBinLabel( stList[i], strStName);

    //grpRmsRelUbOff->SetBinError( stList[i], rmsRmsRelUbOff[i] );
    //grpRmsRelUbOff->GetXaxis()->SetBinLabel( stList[i], strStName);
  }
  
  TF1 *poly1 = new TF1("poly1","[0]+[1]*x", stList[0], stList[nStations]);

  grpRmsRelUubCdas->Fit("poly1","QR");
  grpRmsRelUbCdas->Fit("poly1","QR");
  grpRmsRelUbCdas->GetFunction("poly1")->SetLineColor(kBlack);

  TCanvas *c1 = canvasStyle("c1");
  c1->cd();
  
  grpRmsRelUubCdas->SetStats(kFALSE);
  grpRmsRelUubCdas->SetTitle("");
  grpRmsRelUubCdas->GetYaxis()->SetTitle("#LT (RMS)_{5-days}/#LT Q^{pk}_{VEM} #GT_{5-days} #GT [%]");
  grpRmsRelUubCdas->GetXaxis()->SetTitle("Station Id.");
  grpRmsRelUubCdas->GetXaxis()->SetTitleOffset(1.5);
  grpRmsRelUubCdas->GetYaxis()->SetRangeUser(0, 5.);
  grpRmsRelUubCdas->SetMarkerStyle(8);
  grpRmsRelUubCdas->SetMarkerColor(kRed);
  grpRmsRelUubCdas->SetMarkerSize(1.5);
  grpRmsRelUubCdas->SetLineColor(kRed);
  grpRmsRelUubCdas->Draw("E1");
  c1->Update();

  grpRmsRelUbCdas->SetMarkerStyle(22);
  grpRmsRelUbCdas->SetMarkerSize(1.5);
  grpRmsRelUbCdas->SetLineColor(kBlack);
  grpRmsRelUbCdas->Draw("E1 same");
  c1->Update();

  leg = new TLegend(0.15,0.65,0.4,0.95);
  leg->SetHeader("PMT"+strPmt);
  leg->AddEntry(grpRmsRelUubCdas, "UUB stations", "p");
  TString strSlope;
  strSlope.Form("%5.3g", grpRmsRelUubCdas->GetFunction("poly1")->GetParameter(1));
  TString strInter;
  strInter.Form("%5.3g", grpRmsRelUubCdas->GetFunction("poly1")->GetParameter(0));
  leg->AddEntry(grpRmsRelUubCdas, "Linear fit", "l");
  leg->AddEntry(grpRmsRelUubCdas, "Slope: "+strSlope+" Inter.: "+strInter,"");

  leg->AddEntry(grpRmsRelUbCdas, "UB stations", "p");
  leg->AddEntry(grpRmsRelUbCdas, "Linear fit", "l");
  strSlope.Form("%5.3g", grpRmsRelUbCdas->GetFunction("poly1")->GetParameter(1));
  strInter.Form("%5.3g", grpRmsRelUbCdas->GetFunction("poly1")->GetParameter(0));
  leg->AddEntry(grpRmsRelUbCdas, "Slope: "+strSlope+" Inter.: "+strInter,"");
  leg->SetTextSize(0.035);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();

  c1->Print("../plots/accuracyQpksFitsCdasPmt"+strPmt+".pdf");
/*
  TCanvas *c2 = canvasStyle("c2");
  c2->cd();
  
  grpRmsRelUubOff->SetStats(kFALSE);
  grpRmsRelUubOff->SetTitle("");
  grpRmsRelUubOff->GetYaxis()->SetTitle("RMS/#LT Q^{pk}_{VEM} #GT [au]");
  grpRmsRelUubOff->GetXaxis()->SetTitle("Station Id.");
  grpRmsRelUubOff->GetXaxis()->SetTitleOffset(1.5);
  grpRmsRelUubOff->GetYaxis()->SetRangeUser(0, 0.05);
  grpRmsRelUubOff->SetMarkerStyle(8);
  grpRmsRelUubOff->SetMarkerSize(1.5);
  grpRmsRelUubOff->SetLineColor(kBlack);
  grpRmsRelUubOff->Draw();
  */
/*
  grpRmsRelUbOff->SetMarkerStyle(26);
  grpRmsRelUbOff->SetMarkerSize(1.5);
  grpRmsRelUbOff->SetLineColor(kGray+1);
  grpRmsRelUbOff->Draw("same");
*/
  /*
  leg = new TLegend(0.15,0.75,0.4,0.95);
  leg->AddEntry(grpRmsRelUubOff, "UUB stations", "p");
  //leg->AddEntry(grpRmsRelUbOff, "UB stations", "p");
  leg->SetTextSize(0.035);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();
*/
  //c2->Print("temp2.pdf");
  exit(0);
} 
