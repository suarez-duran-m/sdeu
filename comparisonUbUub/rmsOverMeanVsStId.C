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

void histoStyle(TH1D *hist) {
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
     bool isCtionSt, vector<double> ctionSt, double r1Random ) {
  double mean = qpkVals->GetMean();
  TString strMean;
  TString strRandom;
  strMean.Form("%.4f", mean);
  strRandom.Form("%.6f", r1Random);
  int nbins = 200;
  double frtBin = 0.;
  double lstBin = 2.;
  double qpk = 0.;
  double relQpk = 0.;
  double relQpk2Bin = nbins/2.;
  int bin = 0;
  int nCnts = 0;
  bool doCnts = false;
  TH1D *retRelQpk = new TH1D ("retRelQpk"+strRandom, "", nbins, frtBin, lstBin);

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

vector<vector<double>> fecthCtionData(int pmt, bool isUub) {
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
  string fileName = (isUub) ? 
    "cautionStationsDataPmt" : 
    "cautionStationsDataUbPmt";
  int maxLines = 0;
  switch( pmt ) {
    case 1:
      fileName += "1.dat";
      maxLines = (isUub) ? 23 : 7;
      retData.resize(maxLines);
      break;
    case 2:
      fileName += "2.dat";
      maxLines = (isUub) ? 17 : 7;
      retData.resize(maxLines);
      break;
    default:
      fileName += "3.dat";
      maxLines = (isUub) ? 24 : 5;
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
  fData.close();
  return retData;
}

struct IfCtionAndLnForSt {
  bool IfCtion;
  int lnForSt;
};

IfCtionAndLnForSt ifCationSt(vector<vector<double>> listCtionSts, 
    int currSt) {
  bool retIfCtion = false;
  int lnForSt = 0;
  for ( int ln=0; ln<listCtionSts.size(); ln++ )
    if ( listCtionSts[ln][0] == currSt ) {
      retIfCtion = true;
      lnForSt = ln;
    }

  IfCtionAndLnForSt tmp{retIfCtion, lnForSt};
  return tmp;
}

double getmean( TH1D *arr ) {
  int nb = arr->GetNbinsX();
  int goodVals = 0;
  double mean = 0.;
  for (int bin_i=0; bin_i<nb; bin_i++)
    if ( arr->GetBinContent(bin_i) > 0 
        && arr->GetBinCenter(bin_i) < 4.
        && arr->GetBinCenter(bin_i) > 0.2 ) {
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
    if ( arr->GetBinContent(bin_i) > 0 
        && arr->GetBinCenter(bin_i) < 4.
        && arr->GetBinCenter(bin_i) > 0.2 ) {
      rms += arr->GetBinContent(bin_i)*arr->GetBinCenter(bin_i)*arr->GetBinCenter(bin_i);
      goodVals += arr->GetBinContent(bin_i);
    }
  return sqrt(rms/goodVals-(meanarr*meanarr));
}

TH1D *getDistRelRMSqpk(bool ifUub, vector<double> &relQpkVect, vector<double> &relSgmVect, vector<double> &stidVect,
    vector<vector<int>> &discartedPMTs) {
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
  cautionStPmt1 = fecthCtionData(1, ifUub);
  vector < vector < double > > cautionStPmt2;
  cautionStPmt2 = fecthCtionData(2, ifUub);
  vector < vector < double > > cautionStPmt3;
  cautionStPmt3 = fecthCtionData(3, ifUub);

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
  double farFromSigma = 1.3;
  int cutForMeanByPmt = 0;
  TRandom *r1Random =new TRandom();
  //TH1D *diffMuRMS = new TH1D ("diffMuRMS","", 2000, -1, 1);

  for ( int etry=0; etry<stIds->GetEntries(); etry++ ) {
    
    stIds->GetEntry( etry );
    //if ( stid != 545 )
      //continue;
            
    for ( int i=0; i<nbins; i++ )
      relQpkPmts.push_back(0.);

    //IfCtionAndLnForSt tmp = ifCationSt( cautionStPmt1, stid );
    //isCtionSt = tmp.IfCtion;
    //lnCationSt = tmp.lnForSt;
    distRelQpkPmt1 = getRelQpk( fetchQpkValsPmt1, relQpkPmts, 1,
        isCtionSt, cautionStPmt1[lnCationSt], r1Random->Rndm());
    
    //tmp = ifCationSt( cautionStPmt2, stid );
    //isCtionSt = tmp.IfCtion;
    //lnCationSt = tmp.lnForSt;
    distRelQpkPmt2 = getRelQpk( fetchQpkValsPmt2, relQpkPmts, 2, 
        isCtionSt, cautionStPmt2[lnCationSt], r1Random->Rndm());

    //isCtionSt = false;
    //lnCationSt = 0;
    //tmp = ifCationSt( cautionStPmt3, stid );
    //isCtionSt = tmp.IfCtion;
    //lnCationSt = tmp.lnForSt;
    distRelQpkPmt3 = getRelQpk( fetchQpkValsPmt3, relQpkPmts, 3, 
        isCtionSt, cautionStPmt3[lnCationSt], r1Random->Rndm());

    for ( int bin=1; bin<relQpkPmts.size(); bin++ )
      if ( relQpkPmts[bin] > 0 )
        distRelQpkPmts->SetBinContent( bin-1, relQpkPmts[bin] );

    rmsRelQpkAllPmt = farFromSigma*distRelQpkPmts->GetRMS();    
    cutForMeanByPmt = (ifUub) ? 500 : 50;
  
    if ( distRelQpkPmt1->GetRMS() < rmsRelQpkAllPmt
        && fetchQpkValsPmt1->GetMean() > cutForMeanByPmt )
      fillRelQpkOk(distRelQpkPmt1, distRelQpkOk);

    if ( distRelQpkPmt2->GetRMS() < rmsRelQpkAllPmt 
        && fetchQpkValsPmt2->GetMean() > cutForMeanByPmt )
      fillRelQpkOk(distRelQpkPmt2, distRelQpkOk);
    
    if ( distRelQpkPmt3->GetRMS() < rmsRelQpkAllPmt 
        && fetchQpkValsPmt3->GetMean() > cutForMeanByPmt )
      fillRelQpkOk(distRelQpkPmt3, distRelQpkOk);
    
    Int_t status = distRelQpkOk->Fit("gaus","0Q","",0,2);
    
   
    /*
    TCanvas *c000 = canvasStyle("c000");
    c000->cd();
    distRelQpkPmts->Draw();
    cout << distRelQpkPmts->GetRMS() << endl;
    //distRelQpkPmt2->Draw();
    //distRelQpkOk->Draw();
    gPad->WaitPrimitive();
    */
 
    double errSgmMu = 0.;
    double tmpMuDelSgm = 0.;
    double tmpSgmDelMu = 0.;
    double mu = 0.;
    double sgm = 0.;
    rmsRelQpkAllPmt = 0.;
    if ( status == 0 ) {
      mu = distRelQpkOk->GetFunction("gaus")->GetParameter(1);
      sgm = distRelQpkOk->GetFunction("gaus")->GetParameter(2);
      rmsRelQpkAllPmt = 100.*sgm/mu;
    }

    if ( rmsRelQpkAllPmt > 0 ) {
      retDistRelRMSqpk->Fill( rmsRelQpkAllPmt );
      relQpkVect.push_back( rmsRelQpkAllPmt ); 
     
      errSgmMu = 0.;
      tmpMuDelSgm = (1./mu)*distRelQpkOk->GetRMSError(1);
      tmpSgmDelMu = (sgm/(mu*mu))*distRelQpkOk->GetMeanError(1);       
      errSgmMu = sqrt( tmpMuDelSgm*tmpMuDelSgm + tmpSgmDelMu*tmpSgmDelMu ); 
      relSgmVect.push_back( 100.*errSgmMu);
      stidVect.push_back( stid );
    }
    else {
      relQpkVect.push_back( 0. );
      relSgmVect.push_back( 0. );
      stidVect.push_back( stid );
    }

    if ( ifUub && rmsRelQpkAllPmt > 4. )
      cout << endl << endl
        << "=============================" << endl
        << "MSD stid " << stid 
        << " ifUub " << ifUub
        << " etry " << etry 
        << " " << distRelQpkOk->GetRMS()
        << " " << distRelQpkOk->GetMean() 
        << " " << rmsRelQpkAllPmt << endl << endl;
    
    rmsRelQpkAllPmt = 0.;

    bool ifPlot = false;
    if ( ifUub && ifPlot ) {
      TLegend *leg;
      TString strTitle;
      strTitle.Form("Station %d UUB",stid);
      //strTitle.Form("Station %d UB",stid);
          
      TCanvas *c0 = canvasStyle("c0");
      c0->cd();
                    
      //if ( distRelQpkOk->GetMean() > 0 )
        //diffMuRMS->Fill( distRelQpkOk->GetFunction("gaus")->GetParameter(2) - distRelQpkOk->GetRMS() );

      distRelQpkOk->SetStats(kFALSE);
      distRelQpkOk->GetXaxis()->SetRangeUser(0.92, 1.07);
      //distRelQpkOk->GetXaxis()->SetRangeUser(0.78, 1.26);
      //distRelQpkOk->GetXaxis()->SetRangeUser(0.9, 1.1); 
      distRelQpkOk->GetXaxis()->SetTitle("Q^{pk}_{i}/#LT Q^{pk}_{VEM}#GT_{PMT}");
      distRelQpkOk->GetYaxis()->SetTitle("Counts [au]");
      distRelQpkOk->SetLineColor(kMagenta-3);
      distRelQpkOk->SetLineWidth(2);
      distRelQpkOk->SetTitle(strTitle);
      histoStyle(distRelQpkOk);
      distRelQpkOk->Draw();
      
      distRelQpkPmt1->SetLineColor(kRed);
      distRelQpkPmt1->Draw("same");
      
      distRelQpkPmt2->SetLineColor(kBlue);
      distRelQpkPmt2->Draw("same");
      distRelQpkPmt3->SetLineColor(kGreen+2);
      distRelQpkPmt3->Draw("same");
      
      leg = new TLegend(0.8,0.5,0.9,0.9);
      leg->AddEntry(distRelQpkOk,"All PMTs", "l");
      leg->AddEntry(distRelQpkPmt1, "PMT1", "l");
      leg->AddEntry(distRelQpkPmt2, "PMT2", "l");
      leg->AddEntry(distRelQpkPmt3, "PMT3", "l");
      leg->SetTextSize(0.04);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->Draw();
      strTitle.Form("filteredSt%d",stid);
      //strTitle.Form("filteredUbSt%d",stid);
      //c0->Print("../plots2/"+strTitle+".pdf");
      gPad->WaitPrimitive();
      

      TCanvas *c00 = canvasStyle("c00");
      c00->cd();
      strTitle.Form("Station %d UUB",stid);
      //strTitle.Form("Station %d UB",stid);
      fetchQpkValsPmt1->SetStats(kFALSE);
      fetchQpkValsPmt1->SetTitle(strTitle);
      fetchQpkValsPmt1->SetLineColor(kRed);
      fetchQpkValsPmt1->GetXaxis()->SetTitle("Q^{pk}_{VEM}");
      fetchQpkValsPmt1->GetYaxis()->SetTitle("Counts [au]");
      fetchQpkValsPmt1->GetXaxis()->SetRangeUser(1.3e3, 1.9e3);
      //fetchQpkValsPmt1->GetXaxis()->SetRangeUser(1450, 1850);
      //fetchQpkValsPmt1->GetXaxis()->SetRangeUser(140, 225);
      fetchQpkValsPmt1->GetYaxis()->SetRangeUser(0,17);
      histoStyle(fetchQpkValsPmt1);
      fetchQpkValsPmt1->Draw();
      fetchQpkValsPmt2->SetLineColor(kBlue);
      fetchQpkValsPmt2->Draw("same");
      fetchQpkValsPmt3->SetLineColor(kGreen+2);
      fetchQpkValsPmt3->Draw("same");
      strTitle.Form("filteredPMTsSt%d",stid);
      //strTitle.Form("filteredUbPMTsSt%d",stid);

      leg = new TLegend(0.8,0.5,0.9,0.9);
      leg->AddEntry(fetchQpkValsPmt1, "PMT1", "l");
      leg->AddEntry(fetchQpkValsPmt2, "PMT2", "l");
      leg->AddEntry(fetchQpkValsPmt3, "PMT3", "l");
      leg->SetTextSize(0.04);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->Draw();
      //c00->Print("../plots2/"+strTitle+".pdf");
      gPad->WaitPrimitive();
      
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
    //diffMuRMS->GetXaxis()->SetRangeUser(-0.03, 0.04);
    diffMuRMS->Draw();
    gPad->WaitPrimitive();
    //c0->Print("../plots2/diffMuRMS.pdf");
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
 
  vector<double> relQpkVectUub;
  vector<double> relSgmVectUub;
  vector<double> stidVectUub;

  vector<double> relQpkVectUb;
  vector<double> relSgmVectUb;
  vector<double> stidVectUb;

  vector < vector < int > > discartedPMTs;
  discartedPMTs.resize(1900);
  for ( int i=0; i<1900; i++ ) {
    discartedPMTs[i].resize(3);
    for ( int j=0; j<3; j++ )
      discartedPMTs[i][j] = 0;
  }

  bool ifIsUub = true;
  TH1D *distRelRMSqpkUub;
  distRelRMSqpkUub = getDistRelRMSqpk(ifIsUub, relQpkVectUub, relSgmVectUub, stidVectUub, discartedPMTs);

  ifIsUub = false;
  TH1D *distRelRMSqpkUb;
  distRelRMSqpkUb = getDistRelRMSqpk(ifIsUub, relQpkVectUb, relSgmVectUb, stidVectUb, discartedPMTs);

  gStyle->SetErrorX(0);
  TCanvas *c1 = canvasStyle("c1");
  c1->cd();

  distRelRMSqpkUub->SetStats(kFALSE);
  distRelRMSqpkUub->SetLineColor(kRed);
  distRelRMSqpkUub->SetLineWidth(2);
  distRelRMSqpkUub->GetYaxis()->SetTitle("Counts [au]");
  distRelRMSqpkUub->GetXaxis()->SetTitle("Relative uncertainty [%]"); //#sigma/#mu [%]");
  distRelRMSqpkUub->GetXaxis()->SetRangeUser(0.4, 2.0);
  distRelRMSqpkUub->GetYaxis()->SetRangeUser(0, 12.);
  distRelRMSqpkUub->SetFillStyle(3345);
  distRelRMSqpkUub->SetFillColor(kRed);
  histoStyle(distRelRMSqpkUub);
  distRelRMSqpkUub->Draw();

  distRelRMSqpkUb->SetLineColor(kBlue);
  distRelRMSqpkUb->SetLineWidth(2);
  distRelRMSqpkUb->SetFillStyle(3354);
  distRelRMSqpkUb->SetFillColor(kBlue);
  distRelRMSqpkUb->Draw("same");  

  leg = new TLegend(0.6,0.6,0.9,0.95);
  strMean.Form("%.2f", getmean(distRelRMSqpkUb));
  strRms.Form("%.2f", getrms(distRelRMSqpkUb, getmean(distRelRMSqpkUb)));
  TString strMeanErr;
  strMeanErr.Form("%.2f", distRelRMSqpkUb->GetMeanError());
  leg->AddEntry(distRelRMSqpkUb, "UB", "l" );
  leg->AddEntry(distRelRMSqpkUb, "Mean: "+strMean+" % #pm "+strMeanErr+" %", "");
  
  strMean.Form("%.2f", getmean(distRelRMSqpkUub));
  strRms.Form("%.2f", getrms(distRelRMSqpkUub, getmean(distRelRMSqpkUub)));
  leg->AddEntry(distRelRMSqpkUub, "UUB", "l" );
  strMeanErr.Form("%.2f", distRelRMSqpkUub->GetMeanError());
  leg->AddEntry(distRelRMSqpkUub, "Mean: "+strMean+" % #pm "+strMeanErr+" %", "");
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();
  //c1->Print("kk.pdf");
  c1->Print("../plots2/accuracyQpksFitsUbUubAllStAllPmt_v2.pdf");

  int nStations = 23;
  TH1D *RelQpkUub = new TH1D("RelQpkUub", "", nStations, 0, nStations);
  TH1D *RelQpkUb = new TH1D("RelQpkUb", "", nStations, 0, nStations);
  TH1D *RelDiffQpk = new TH1D("RelDiffQpk", "", nStations, 0, nStations);
  TH1D *DistRelDiffQpk = new TH1D("DistRelDiffQpk", "", 100, -5.0, 5.0);

  double aveUub = 0.;
  double aveUb = 0.;
  double tmpRelQpkUub = 0.;
  double tmpRelQpkUb = 0.;
  TString strStName;
  for ( int bin_i=1; bin_i<nStations+1; bin_i++ ) {
    tmpRelQpkUub = relQpkVectUub[bin_i-1];
    tmpRelQpkUb = relQpkVectUb[bin_i-1];
    strStName.Form("%d", (int)stidVectUub[bin_i-1]);
 
    RelQpkUub->SetBinContent(bin_i, tmpRelQpkUub);
    RelQpkUub->SetBinError( bin_i, relSgmVectUub[bin_i-1] );
    RelQpkUub->GetXaxis()->SetBinLabel(bin_i, strStName);
    aveUub += tmpRelQpkUub;
 
    RelQpkUb->SetBinContent(bin_i, tmpRelQpkUb);
    RelQpkUb->SetBinError( bin_i, relSgmVectUb[bin_i-1] );
    RelQpkUb->GetXaxis()->SetBinLabel(bin_i, strStName); 
    aveUb += tmpRelQpkUb;

    //if ( tmpRelQpkUb > 0 && tmpRelQpkUb < 4. && tmpRelQpkUub < 4. ) {
      RelDiffQpk->SetBinContent(bin_i, (tmpRelQpkUb - tmpRelQpkUub));
      RelDiffQpk->SetBinError(bin_i, sqrt( pow(relSgmVectUb[bin_i-1], 2) + pow(relSgmVectUub[bin_i-1],2 )) );
      DistRelDiffQpk->Fill( tmpRelQpkUb - tmpRelQpkUub );
    //}
    /*
    else {
      RelDiffQpk->SetBinContent(bin_i, -100.);
      DistRelDiffQpk->Fill(-100.);
      cout << strStName << endl;
    }
    */
    RelDiffQpk->GetXaxis()->SetBinLabel(bin_i, strStName);
  }
  aveUub /= relQpkVectUub.size();
  aveUb /= relQpkVectUb.size();
 
  TCanvas *c2 = canvasStyle("c2"); 
  c2->cd();

  RelQpkUub->SetTitle("");
  RelQpkUub->SetStats(kFALSE);
  RelQpkUub->GetXaxis()->SetTitle("Station ID.");
  RelQpkUub->GetXaxis()->SetTitleOffset(1.5);
  //RelQpkUub->GetYaxis()->SetRangeUser(0.2, 2.4);
  RelQpkUub->GetYaxis()->SetTitle("#sigma/#mu [%]");
  RelQpkUub->SetMarkerStyle(71);
  RelQpkUub->SetMarkerColor(kRed);
  RelQpkUub->SetLineColor(kRed);
  RelQpkUub->SetMarkerSize(1.2);
  histoStyle(RelQpkUub);
  RelQpkUub->Draw("E1");

  TLine *lineUub = new TLine(0.5,aveUub,22.5,aveUub);
  lineUub->SetLineWidth(2);
  lineUub->SetLineColor(kRed);
  lineUub->Draw();

  RelQpkUb->SetMarkerStyle(73);
  RelQpkUb->SetMarkerColor(kBlue);
  RelQpkUb->SetLineColor(kBlue);
  RelQpkUb->SetMarkerSize(1.2);
  RelQpkUb->Draw("E1 same");

  TLine *lineUb = new TLine(0.5,aveUb,22.5,aveUb);
  lineUb->SetLineWidth(2);
  lineUb->SetLineColor(kBlue);
  lineUb->Draw();
 
  leg = new TLegend(0.88,0.8,0.98,0.95);
  //leg->SetHeader("Station 827, UUB");
  leg->AddEntry(RelQpkUub, "UUB", "p");
  leg->AddEntry(RelQpkUb, "UB", "p");
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw(); 
  c2->Print("../plots2/relSigmaQpkVsStationsId.pdf");

  TCanvas *c3 = canvasStyle("c3");
  c3->cd();

  RelDiffQpk->SetTitle("");
  RelDiffQpk->SetStats(kFALSE);                       
  RelDiffQpk->GetXaxis()->SetTitle("Station ID.");
  RelDiffQpk->GetXaxis()->SetTitleOffset(1.5);
  //RelDiffQpk->GetYaxis()->SetRangeUser(-0.9, 0.9);
  RelDiffQpk->GetYaxis()->SetTitle("#left(#sigma/#mu#right)_{UB} - #left(#sigma/#mu#right)_{UUB} [%]");
  RelDiffQpk->SetMarkerStyle(71);
  RelDiffQpk->SetMarkerColor(kGreen+3);
  RelDiffQpk->SetLineColor(kGreen+3);
  RelDiffQpk->SetMarkerSize(1.2);
  histoStyle(RelDiffQpk);
  RelDiffQpk->Draw("P");

  lineUb = new TLine(0.5,0.,22.5,0.);
  lineUb->SetLineWidth(2);
  lineUb->SetLineColor(kGray);
  lineUb->SetLineStyle(2);
  lineUb->Draw();
  c3->Print("../plots2/diffRelSigmaQpkVsStationsId.pdf");

  TCanvas *c4 = canvasStyle("c4");
  c4->cd();

  DistRelDiffQpk->SetTitle("");
  DistRelDiffQpk->SetStats(kFALSE); 
  DistRelDiffQpk->GetXaxis()->SetTitle("#left(#sigma/#mu#right)_{UB} - #left(#sigma/#mu#right)_{UUB} [%]");
  DistRelDiffQpk->GetXaxis()->SetRangeUser(-1., 1.);
  //DistRelDiffQpk->GetXaxis()->SetTitleOffset(1.5);
  //DistRelDiffQpk->GetYaxis()->SetRangeUser(0.9, 1.07);
  DistRelDiffQpk->GetYaxis()->SetTitle("Counts [au]"); 
  DistRelDiffQpk->SetLineColor(kGreen+3);
  histoStyle(DistRelDiffQpk);
  DistRelDiffQpk->Draw();

  leg = new TLegend(0.72,0.78,0.98,0.9);
  strMean.Form("%.3f", DistRelDiffQpk->GetMean());
  strRms.Form("%.3f", DistRelDiffQpk->GetRMS());
  leg->AddEntry(DistRelDiffQpk, "MEAN: "+strMean, "");
  leg->AddEntry(DistRelDiffQpk, "RMS: "+strRms, "");
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();

  lineUb = new TLine(0.,0.2,0.,4.2);
  lineUb->SetLineWidth(2);
  lineUb->SetLineColor(kGray);
  lineUb->SetLineStyle(2);
  lineUb->Draw();
  
  c4->Print("../plots2/distDiffSigmaQpk.pdf");

  //exit(0);
} 
