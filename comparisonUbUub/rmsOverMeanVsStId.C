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

TH1D *getRelQpk( TH1D *qpkVals, vector<double> &relQpkAllPmts, int pmt,//) {
    bool isCtionSt, vector<double> ctionSt ) {
  double mean = qpkVals->GetMean();
  TString strMean;
  strMean.Form("%.4f", mean);
  int nbins = 200;
  double frtBin = 0.;
  double lstBin = 2.;
  double qpk = 0.;
  double relQpk = 0.;
  double relQpk2Bin = nbins/2.;
  int bin = 0;
  int nCnts = 0;
  bool doCnts = false;
  TH1D *retRelQpk = new TH1D ("retRelQpk"+strMean, "", nbins, frtBin, lstBin);
  for ( int qpk_i=1; qpk_i<qpkVals->GetNbinsX(); qpk_i++ ) {
    if ( qpkVals->GetBinContent(qpk_i) <= 0 )
      continue;
    qpk = qpkVals->GetBinLowEdge(qpk_i);
    if ( isCtionSt ) {
      switch( pmt ) {
      case 3:
        if ( qpk > ctionSt[1]-ctionSt[2] && qpk < ctionSt[1]+ctionSt[2] )
          relQpk = qpk/ctionSt[1];
        else if ( qpk > ctionSt[3]-ctionSt[4] && qpk < ctionSt[3]+ctionSt[4] )
          relQpk = qpk/ctionSt[3];
        else if ( qpk > ctionSt[5]-ctionSt[6] && qpk < ctionSt[5]+ctionSt[6] )
          relQpk = qpk/ctionSt[5];
        else if ( qpk > ctionSt[7]-ctionSt[8] && qpk < ctionSt[7]+ctionSt[8] )
          relQpk = qpk/ctionSt[7];
        else
          relQpk = 0.;
        break;
      default:
        if ( qpk > ctionSt[1]-ctionSt[2] && qpk < ctionSt[1]+ctionSt[2] )
          relQpk = qpk/ctionSt[1];
        else if ( qpk > ctionSt[3]-ctionSt[4] && qpk < ctionSt[3]+ctionSt[4] )
          relQpk = qpk/ctionSt[3];
        else if ( qpk > ctionSt[5]-ctionSt[6] && qpk < ctionSt[5]+ctionSt[6] )
          relQpk = qpk/ctionSt[5];
        else
          relQpk = 0.;
      }
      if ( relQpk == 0 )
        continue;
      bin = relQpk2Bin*(relQpk+(lstBin/nbins));
      nCnts = qpkVals->GetBinContent(qpk_i);
      nCnts += retRelQpk->GetBinContent( bin-1 );
      relQpkAllPmts[bin] += nCnts;
      // Add the entries for same qpk/mean bins, before the
      // new SetBinContent
      retRelQpk->SetBinContent( bin-1, nCnts );
    }
    else {
      relQpk = qpk/mean;
      bin = relQpk2Bin*(relQpk+(lstBin/nbins));
      nCnts = qpkVals->GetBinContent(qpk_i);
      relQpkAllPmts[bin] += nCnts;
      // Add the entries for same qpk/mean bins, before the
      // new SetBinContent
      nCnts += retRelQpk->GetBinContent( bin-1 );
      retRelQpk->SetBinContent( bin-1, nCnts );
    }
  }
  return retRelQpk;
}

void fillRelQpkOk( TH1D *relQpkVals, TH1D *modifingRelQpkOk ) {
  int nCnts = 0;
  for ( int qpk_i=1; qpk_i<relQpkVals->GetNbinsX(); qpk_i++ )
    if ( relQpkVals->GetBinContent(qpk_i) > 0 ) {
      nCnts = relQpkVals->GetBinContent(qpk_i);
      // Add entries for same relQpk bin for previous PMTs
      nCnts += modifingRelQpkOk->GetBinContent( qpk_i );
      modifingRelQpkOk->SetBinContent( qpk_i, nCnts );
    }
}

vector<vector<double>> fecthCtionData(int pmt) {
  vector<vector< double >> retData;
  
  ifstream fData;
  int stId = 0.;
  double mean1 = 0.;
  double sgm1 = 0.;
  double mean2 = 0.;
  double sgm2 = 0.;
  double mean3 = 0.;
  double sgm3 = 0.;
  double mean4 = 0.;
  double sgm4 = 0.;
  string fileName = "cautionStationsDataPmt";
  int maxLines = 0;
  switch( pmt ) {
    case 1:
      fileName += "1.dat";
      maxLines = 21;
      retData.resize(maxLines);
      break;
    case 2:
      fileName += "2.dat";
      maxLines = 16;
      retData.resize(maxLines);
      break;
    default:
      fileName += "3.dat";
      maxLines = 22;
      retData.resize(maxLines);
  }
  int cntLines = 0;

  fData.open(fileName);
  while( fData.good() ) {
    if ( cntLines == maxLines )
      break;
    switch( pmt ) {
      case 3:
        fData >> stId 
          >> mean1 >> sgm1
          >> mean2 >> sgm2
          >> mean3 >> sgm3
          >> mean4 >> sgm4;
        retData[cntLines].push_back(stId);
        retData[cntLines].push_back(mean1);
        retData[cntLines].push_back(sgm1);
        retData[cntLines].push_back(mean2);
        retData[cntLines].push_back(sgm2);
        retData[cntLines].push_back(mean3);
        retData[cntLines].push_back(sgm3);
        retData[cntLines].push_back(mean4);
        retData[cntLines].push_back(sgm4);
        break;
      default:
        fData >> stId 
          >> mean1 >> sgm1 >> mean2 >> sgm2 >> mean3 >> sgm3;
        retData[cntLines].push_back(stId);
        retData[cntLines].push_back(mean1);
        retData[cntLines].push_back(sgm1);
        retData[cntLines].push_back(mean2);
        retData[cntLines].push_back(sgm2);
        retData[cntLines].push_back(mean3);
        retData[cntLines].push_back(sgm3);
    }
    cntLines++;
  }
  retData.pop_back();
  fData.close();
  return retData;
}

struct IfCtionAndLnForSt {
  bool IfCtion;
  int lnForSt;
};

IfCtionAndLnForSt ifCationSt(vector<vector<double>> listCtionSts, 
    int currSt, bool isUub) {
  bool retIfCtion = false;
  int lnForSt = 0;
  if ( isUub ) {
    for ( int ln=0; ln<listCtionSts.size(); ln++ )
      if ( listCtionSts[ln][0] == currSt ) {
        retIfCtion = true;
        lnForSt = ln;
      }
  }
  else {
    retIfCtion = false;
    lnForSt = 0;
  }
  IfCtionAndLnForSt tmp{retIfCtion, lnForSt};
  return tmp;
}


double getmean( TH1D *arr ) {
  int nb = arr->GetNbinsX();
  int goodVals = 0;
  double mean = 0.;
  for (int bin_i=0; bin_i<nb; bin_i++)
    if ( arr->GetBinContent(bin_i) > 0 && arr->GetBinCenter(bin_i)<3.5) {
      mean += arr->GetBinCenter(bin_i)*arr->GetBinContent(bin_i);
      goodVals += arr->GetBinContent(bin_i);
    }
  return mean/goodVals;
}

double getrms( TH1D *arr, double meanarr ) {
  int nb = arr->GetNbinsX();
  double rms = 0.;
  int goodVals = 0;
  for (int bin_i=0; bin_i<nb; bin_i++)
    if ( arr->GetBinContent(bin_i) > 0 && arr->GetBinCenter(bin_i)<3.5 ) {
      rms += arr->GetBinContent(bin_i)*arr->GetBinCenter(bin_i)*arr->GetBinCenter(bin_i);
      goodVals += arr->GetBinContent(bin_i);
    }
  return sqrt(rms/goodVals-(meanarr*meanarr));
}

TH1D *getDistRelRMSqpk(bool ifUub) {
  TString typeSt = (ifUub) ? "StationsUub" : "StationsUb";
  TString strFname = (ifUub) ? "qpkValStationsUub.root" : "qpkValStationsUb.root"; 
  TFile *f = TFile::Open(strFname);
  TTree *stIds = (TTree*)f->Get(typeSt);
  TH1D *fetchQpkValsPmt1 = new TH1D();
  TH1D *fetchQpkValsPmt2 = new TH1D();
  TH1D *fetchQpkValsPmt3 = new TH1D();

  int stid = 0;
  
  stIds->SetBranchAddress("qpkValuesPmt1", &fetchQpkValsPmt1);
  stIds->SetBranchAddress("qpkValuesPmt2", &fetchQpkValsPmt2);
  stIds->SetBranchAddress("qpkValuesPmt3", &fetchQpkValsPmt3);
  stIds->SetBranchAddress("stId", &stid);

  bool isCtionSt = false;
  int lnCationSt = 0;

  vector < vector < double > > cautionStPmt1;
  cautionStPmt1 = fecthCtionData(1);
  vector < vector < double > > cautionStPmt2;
  cautionStPmt2 = fecthCtionData(2);
  vector < vector < double > > cautionStPmt3;
  cautionStPmt3 = fecthCtionData(3);

  int nbins = 200;
  double frtBin = 0.;
  double lstBin = 2.;
  TH1D *distRelQpkPmt1 = new TH1D();
  TH1D *distRelQpkPmt2 = new TH1D(); 
  TH1D *distRelQpkPmt3 = new TH1D(); 
  vector < double > relQpkPmts;

  for ( int i=0; i<nbins; i++ )
    relQpkPmts.push_back(0.);
  TH1D *distRelQpkPmts = new TH1D("distRelQpkPmts", "", 
      nbins, frtBin, lstBin);
  TH1D *distRelQpkOk = new TH1D("distRelQpkOk", "",
      nbins, frtBin, lstBin);

  TH1D *retDistRelRMSqpk = new TH1D("distRelRMSqpk", "", 100, 0., 20.);

  double rmsRelQpkAllPmt = 0.;
  double rmsRelQpkPmt = 0.;
  double farFromSigma = 1.3;
  int cutForMeanByPmt = 0;

  //TH1D *diffMuRMS = new TH1D ("diffMuRMS","", 2000, -1, 1);

  for ( int etry=0; etry<stIds->GetEntries(); etry++ ) {
    // 8 for st827; 12 for st833; 67 for 1746; 14 for 843
    stIds->GetEntry( etry );
    //if ( stid != 1878 )
      //continue;
            
    for ( int i=0; i<nbins; i++ )
      relQpkPmts.push_back(0.);

    IfCtionAndLnForSt tmp = ifCationSt( cautionStPmt1, stid, ifUub );
    isCtionSt = tmp.IfCtion;
    lnCationSt = tmp.lnForSt;
    distRelQpkPmt1 = getRelQpk( fetchQpkValsPmt1, relQpkPmts, 1, 
        isCtionSt, cautionStPmt1[lnCationSt]);

    tmp = ifCationSt( cautionStPmt2, stid, ifUub);
    isCtionSt = tmp.IfCtion;
    lnCationSt = tmp.lnForSt;
    distRelQpkPmt2 = getRelQpk( fetchQpkValsPmt2, relQpkPmts, 2, 
        isCtionSt, cautionStPmt2[lnCationSt]);

    tmp = ifCationSt( cautionStPmt3, stid, ifUub);
    isCtionSt = tmp.IfCtion;
    lnCationSt = tmp.lnForSt;
    distRelQpkPmt3 = getRelQpk( fetchQpkValsPmt3, relQpkPmts, 3, 
        isCtionSt, cautionStPmt3[lnCationSt]);

    for ( int bin=1; bin<relQpkPmts.size(); bin++ )
      if ( relQpkPmts[bin] > 0 )
        distRelQpkPmts->SetBinContent( bin-1, relQpkPmts[bin] );

    rmsRelQpkAllPmt = farFromSigma*distRelQpkPmts->GetRMS();
    
    /*
    TCanvas *c000 = canvasStyle("c000");
    c000->cd();
    distRelQpkPmts->Draw();
    gPad->WaitPrimitive();
    */
    cutForMeanByPmt = (ifUub) ? 500 : 50;

    rmsRelQpkPmt = distRelQpkPmt1->GetRMS();
    if ( rmsRelQpkPmt < rmsRelQpkAllPmt 
        && fetchQpkValsPmt1->GetMean() > cutForMeanByPmt )
      fillRelQpkOk(distRelQpkPmt1, distRelQpkOk);

    rmsRelQpkPmt = distRelQpkPmt2->GetRMS();
    if ( rmsRelQpkPmt < rmsRelQpkAllPmt 
        && fetchQpkValsPmt2->GetMean() > cutForMeanByPmt )
      fillRelQpkOk(distRelQpkPmt2, distRelQpkOk);

    rmsRelQpkPmt = distRelQpkPmt3->GetRMS();
    if ( rmsRelQpkPmt < rmsRelQpkAllPmt 
        && fetchQpkValsPmt3->GetMean() > cutForMeanByPmt )
      fillRelQpkOk(distRelQpkPmt3, distRelQpkOk);
    
    Int_t status = distRelQpkOk->Fit("gaus","","",distRelQpkOk->GetMean()-2.*distRelQpkOk->GetRMS(),
        distRelQpkOk->GetMean()+2.*distRelQpkOk->GetRMS());
   
    if ( status == 0 )
      rmsRelQpkAllPmt = 100.*distRelQpkOk->GetFunction("gaus")->GetParameter(2)/distRelQpkOk->GetFunction("gaus")->GetParameter(1);

    if ( rmsRelQpkAllPmt > 0 )
      retDistRelRMSqpk->Fill( rmsRelQpkAllPmt );
    if ( !ifUub && rmsRelQpkAllPmt > 4.5 ) // && rmsRelQpkAllPmt < 2.9 )
      cout << endl << endl
        << "=============================" << endl
        << "MSD stid " << stid 
        << " ifUub " << ifUub
        << " etry " << etry 
        << " " << distRelQpkOk->GetRMS()
        << " " << distRelQpkOk->GetMean() << endl << endl;
    
    if ( !ifUub ) {
      TLegend *leg;
      TString strTitle;
      strTitle.Form("Station %d UUB",stid);
      //strTitle.Form("Station %d UB",stid);

      /*
      TCanvas *c0 = canvasStyle("c0");
      c0->cd();
        
      //distRelQpkOk->Fit("gaus","","",distRelQpkOk->GetMean()-distRelQpkOk->GetRMS(),
        //  distRelQpkOk->GetMean()+distRelQpkOk->GetRMS());
      
      //if ( distRelQpkOk->GetMean() > 0 )
        //diffMuRMS->Fill( distRelQpkOk->GetFunction("gaus")->GetParameter(2) - distRelQpkOk->GetRMS() );
      //cout << "MSD gaus " << distRelQpkOk->GetFunction("gaus")->GetParameter(2) << endl;
      //cout << "MSD dist " << rmsRelQpkAllPmt << endl;
            
      distRelQpkOk->SetStats(kFALSE);
      distRelQpkOk->GetXaxis()->SetRangeUser(0.78, 1.26);
      distRelQpkOk->GetXaxis()->SetRangeUser(0.9, 1.1); 
      distRelQpkOk->GetXaxis()->SetTitle("Q^{pk}_{i}/#LT Q^{pk}_{VEM}#GT_{PMT}");
      distRelQpkOk->GetYaxis()->SetTitle("Counts [au]");
      distRelQpkOk->SetLineColor(kMagenta-3);
      distRelQpkOk->SetTitle(strTitle);
      distRelQpkOk->Draw();
      
      distRelQpkPmt1->SetLineColor(kBlack);
      distRelQpkPmt1->Draw("same");
      
      distRelQpkPmt2->SetLineColor(kBlue);
      distRelQpkPmt2->Draw("same");
      distRelQpkPmt3->SetLineColor(kGreen+2);
      distRelQpkPmt3->Draw("same");
      
      leg = new TLegend(0.8,0.5,0.9,0.9);
      //leg->SetHeader("Filtered");
      leg->AddEntry(distRelQpkOk->GetFunction("gaus"),"Gaus fit", "l");
      strTitle.Form("%.5f", 
          distRelQpkOk->GetFunction("gaus")->GetParameter(1));
      leg->AddEntry(distRelQpkOk->GetFunction("gaus"),"#mu = "+strTitle, "");
      strTitle.Form("%.5f",
          distRelQpkOk->GetFunction("gaus")->GetParameter(2));
      leg->AddEntry(distRelQpkOk->GetFunction("gaus"),"#sigma = "+strTitle, "");
      leg->AddEntry(distRelQpkOk,"All PMTs", "l");
      leg->AddEntry(distRelQpkPmt1, "PMT1", "l");
      leg->AddEntry(distRelQpkPmt2, "PMT2", "l");
      leg->AddEntry(distRelQpkPmt3, "PMT3", "l");
      leg->SetTextSize(0.03);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->Draw();
      strTitle.Form("filteredSt%d",stid);
      //strTitle.Form("filteredUbSt%d",stid);
      //strTitle.Form("nofilterSt%d",stid);
      //c0->Print("../plots/"+strTitle+".pdf");
      gPad->WaitPrimitive();
      */
      /*
      TCanvas *c00 = canvasStyle("c00");
      c00->cd();
      strTitle.Form("Station %d UUB",stid);
      fetchQpkValsPmt1->SetStats(kFALSE);
      fetchQpkValsPmt1->SetTitle(strTitle);
      fetchQpkValsPmt1->SetLineColor(kRed);
      fetchQpkValsPmt1->GetXaxis()->SetTitle("Q^{pk}_{VEM}");
      fetchQpkValsPmt1->GetYaxis()->SetTitle("Counts [au]");
      fetchQpkValsPmt1->GetXaxis()->SetRangeUser(1450, 1850);
      //fetchQpkValsPmt1->GetXaxis()->SetRangeUser(145, 195);
      fetchQpkValsPmt1->GetYaxis()->SetRangeUser(0, 45);
      fetchQpkValsPmt1->Draw();
      fetchQpkValsPmt2->SetLineColor(kBlue);
      fetchQpkValsPmt2->Draw("same");
      fetchQpkValsPmt3->SetLineColor(kGreen+2);
      fetchQpkValsPmt3->Draw("same");
      strTitle.Form("filteredPMTsSt%d",stid);
      //strTitle.Form("filteredUbPMTsSt%d",stid);
      //c00->Print("../plots/"+strTitle+".pdf");
      gPad->WaitPrimitive();
      */
    }
    distRelQpkPmt1->Reset();
    distRelQpkPmt2->Reset();
    distRelQpkPmt3->Reset();
    distRelQpkPmts->Reset();
    distRelQpkOk->Reset();

    relQpkPmts.clear();
  }
  /*
  if( ifUub ) {
    TCanvas *c0 = canvasStyle("c0");
    c0->cd();
    diffMuRMS->GetYaxis()->SetTitle("Counts [au]");
    diffMuRMS->GetXaxis()->SetTitle("(#mu - RMS) [au]");
    diffMuRMS->GetXaxis()->SetRangeUser(-0.03, 0.04);
    diffMuRMS->Draw();
    c0->Print("../plots/diffMuRMS.pdf");
  }
  */
  return retDistRelRMSqpk; 
  stIds->Clear();
  f->Close();
}

void rmsOverMeanVsStId() {

  TPaveStats *ptstats;
  TLegend *leg;
  TString strMean;
  TString strRms;
 
  bool ifIsUub = true;
  TH1D *distRelRMSqpkUub;
  distRelRMSqpkUub = getDistRelRMSqpk(ifIsUub);

  ifIsUub = false;
  TH1D *distRelRMSqpkUb;
  distRelRMSqpkUb = getDistRelRMSqpk(ifIsUub);


  TCanvas *c1 = canvasStyle("c1");
  c1->cd();

  distRelRMSqpkUub->SetStats(kFALSE);
  distRelRMSqpkUub->SetLineColor(kRed);
  distRelRMSqpkUub->SetLineWidth(2);
  distRelRMSqpkUub->GetYaxis()->SetTitle("Counts [au]");
  distRelRMSqpkUub->GetXaxis()->SetTitle("#sigma/#mu [%]");
  distRelRMSqpkUub->GetXaxis()->SetRangeUser(0, 10);
  distRelRMSqpkUub->GetYaxis()->SetRangeUser(0, 13.5);
  distRelRMSqpkUub->SetFillStyle(3345);
  distRelRMSqpkUub->SetFillColor(kRed);
  distRelRMSqpkUub->Draw();

  distRelRMSqpkUb->SetLineColor(kBlack);
  distRelRMSqpkUb->SetLineWidth(2);
  distRelRMSqpkUb->SetFillStyle(3354);
  distRelRMSqpkUb->SetFillColor(kBlack);
  distRelRMSqpkUb->Draw("same");  

  leg = new TLegend(0.7,0.65,0.9,0.95);
  
  strMean.Form("%.3f", getmean(distRelRMSqpkUb));
  strRms.Form("%.3f", getrms(distRelRMSqpkUb, getmean(distRelRMSqpkUb)));
  leg->AddEntry(distRelRMSqpkUb, "UB", "l" );
  leg->AddEntry(distRelRMSqpkUb, "MEAN (<3.5 %): "+strMean, "");
  leg->AddEntry(distRelRMSqpkUb, "RMS (<3.5 %): "+strRms, "");
  
  strMean.Form("%.3f", getmean(distRelRMSqpkUub));
  strRms.Form("%.3f", getrms(distRelRMSqpkUub, getmean(distRelRMSqpkUub)));
  leg->AddEntry(distRelRMSqpkUub, "UUB", "l" );
  leg->AddEntry(distRelRMSqpkUub, "MEAN (<3.5 %): "+strMean, "");
  leg->AddEntry(distRelRMSqpkUub, "RMS (<3.5 %): "+strRms, "");
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();
  //c1->Print("kk.pdf");
  c1->Print("../plots/accuracyQpksFitsUbUubAllStAllPmt_v2.pdf");

  //exit(0);
} 
