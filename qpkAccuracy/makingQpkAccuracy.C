/* ===================================
 * Code to calculate the Qpk Accuracy 
 * applying a moving window algorithm
 * Author: Mauricio Suárez Durán
 * ULB, IIHE
 * ===================================
*/

TCanvas *canvasStyle(TString name) {
  TCanvas *canvas = new TCanvas(name, name, 102, 76, 1600, 900);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(2);
  canvas->SetRightMargin(0.04);
  canvas->SetLeftMargin(0.13);
  canvas->SetTopMargin(0.014);
  canvas->SetBottomMargin(0.15);
  canvas->SetFrameBorderMode(0);
  canvas->SetFrameBorderMode(0);
  return canvas;
}

bool cutForChi2SlpLogPval(double chi2ndf, double slpNorm, int pmt, bool ifUub) {
  double cutChi2 = 0.;
  switch ( pmt ) {
    case 1 :
      cutChi2 = (ifUub) ? (3.78+3.11) : (5.61+4.21);
      break;
    case 2 :
      cutChi2 = (ifUub) ? (4.29+3.55) : (5.56+4.35);
      break;
    case 3 :
      cutChi2 = (ifUub) ? (3.92+3.45) : (6.01+4.53);
        break;
  }
  if ( chi2ndf > cutChi2 )
    return true;
  double cutSlpMin = 0.;
  double cutSlpMax = 0.;
  switch ( pmt ) {
    case 1 :
      cutSlpMin = (ifUub) ? (-9.42e-04-2.70e-03) : (-6.23e-04-2.65e-03);
      cutSlpMax = (ifUub) ? (-9.42e-04+2.70e-03) : (-6.23e-04+2.65e-03);
      break;
    case 2 :
      cutSlpMin = (ifUub) ? (-1.10e-03-4.13e-03) : (-6.25e-04-3.16e-03);
      cutSlpMax = (ifUub) ? (-1.10e-03+4.13e-03) : (-6.25e-04+3.16e-03);
      break;
    case 3 :
      cutSlpMin = (ifUub) ? (-9.70e-04-4.29e-03) : (-6.42e-04-6.02e-03);
      cutSlpMax = (ifUub) ? (-9.70e-04+4.29e-03) : (-6.42e-04+6.02e-03);
      break;
  }
  if ( slpNorm < cutSlpMin || cutSlpMax < slpNorm)
    return true;
  if ( log10(TMath::Prob(chi2ndf*5., 5)) < -6. )
    return true;
  return false;
}

void fillQpkTimeVals(bool ifIsUub, int pmt, int st_id, 
    vector<double> &retQpkVect, vector<double> &retTimeVect) {

  // String variable to read root-files with fitted Qpk values
  TString bnCdas = (ifIsUub) ?
    "~/postdoc/sdeu/underHisto/results/uubChPkPMT" :
    "~/postdoc/sdeu/nouub/underHistos/results/ubChPkPMT"; 

  TString pmtId;
  TString strStId;
  TString fname;
  pmtId.Form("%d", pmt);
  strStId.Form("St%d", (int)st_id);

  // Interval time for Qpk accuracy calculation
  TString monthUub[4] = {"Aug", "Sep", "Oct", "Nov"};
  int nMonths = sizeof(monthUub)/sizeof(*monthUub);

  int nYears = 5;
  TString strYear[5] = {"2016", "2018", "2019", "2020", "2021"};
  int stYear = (ifIsUub) ? nYears-1 : 1;
  int lstYear = (ifIsUub) ? nYears-1 : 1;

  // Reading TTree with Qpk values 
  TString strChargeData = "ChargeData";
  TFile *f;
  TTree *chargeInfo;
  double fetchQpkVals = 0.;
  int fetchTime = 0;

  // Fetching Qpk as function of time
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
        // Storing and returning Qpk and time values
        retQpkVect.push_back( fetchQpkVals );
        retTimeVect.push_back( fetchTime );
      }
      f->Clear();
      f->Close();
    }
  }
}

void doAvePerTime(bool ifUub, vector<double> qpkVect, vector<double> timeVect, 
    vector<int> &retNday, vector<double> &retAveQpkDay, 
    vector<double> &retQpkErr, vector<int> &retQpkPerDay, 
    vector<vector<double>> &retqpksInDay) {
  // Starting date (1st Aug.) for average calculation
  const int firstDay = (ifUub) ? 1311811218 : 1217116818;
  int currentDay = firstDay;
  const int oneDay = 86400;
  double qpkAveDay = 0.;
  double qpk2 = 0.;
  double rms = 0.;
  int qpkInDay = 0;
  int timeDiff = 0;
  // Counting the number of days from 1st Aug.
  int nDayOfMonth = 0;
  // Doing Qpk average per day 
  // Moving through Qpks fitted
  for ( int qpk_i=0; qpk_i < qpkVect.size(); qpk_i++ ) {
    // Time difference to identify if current day has gone
    timeDiff = timeVect[qpk_i] - currentDay;
    // Checking if current day has gone
    if ( timeDiff > oneDay ) {
      // Calculating current number day from 1st Aug.
      nDayOfMonth = (int)(currentDay-firstDay)/86400;
      // If <Qpk> > 0, so storing the respective variables
      if ( qpkAveDay > 0 ) {
        qpkAveDay /= qpkInDay;
        // Using the error of the mean
        // rms/sqrt(qpkInDay)
        rms = sqrt(qpk2/qpkInDay - qpkAveDay*qpkAveDay);
        // Filling vectors corresponding with current Qpk-average-day
        retNday[nDayOfMonth] = nDayOfMonth;
        retAveQpkDay[nDayOfMonth] = qpkAveDay;
        retQpkErr[nDayOfMonth] = rms/sqrt(qpkInDay);
        retQpkPerDay[nDayOfMonth] = qpkInDay;
      }
      else {
        retNday[nDayOfMonth] = nDayOfMonth;
        retAveQpkDay[nDayOfMonth] = 0.;
        retQpkErr[nDayOfMonth] = 0.;
      }
      qpkAveDay = 0.;
      qpk2 = 0.;
      qpkInDay = 0;
      currentDay += oneDay;
    }
    // Doing Qpk day average
    retqpksInDay[nDayOfMonth].push_back(qpkVect[qpk_i]);
    qpkAveDay += qpkVect[qpk_i];
    qpk2 += qpkVect[qpk_i]*qpkVect[qpk_i];
    qpkInDay++;
    // Checking for time discontinuity in Qpk series.
    if ( timeDiff > 2*oneDay )
      currentDay += oneDay*(timeDiff/oneDay)-oneDay;
  }
}

TH1D *doingQpkAveDayDist(vector<vector<vector<int>>> qpkPerDay, TString a) {
  TH1D *dist = new TH1D(a, "", 50, 1, 51);
  for ( int st_i=0; st_i<qpkPerDay.size(); st_i++ )
    for ( int pmt_i=0; pmt_i<3; pmt_i++ )
      for ( auto & qpk_i : qpkPerDay[st_i][pmt_i] )
        dist->Fill(qpk_i);
  return dist;
}

vector< double > doPlotDistQpkPerDay(vector<vector<vector<int>>> qpkPerDayUb,
    vector<vector<vector<int>>> qpkPerDayUub, bool ifPlot) {
  // Vector to return fitting parameters
  vector < double > retMusgm(4);  
  TString ifPlotGaus = (ifPlot) ? "Q" : "Q0";
  // Building the Qpk per day distributions
  TH1D *qpkDayDistUb = doingQpkAveDayDist(qpkPerDayUb, "UB");
  qpkDayDistUb->Fit("gaus", ifPlotGaus, "");
  TH1D *qpkDayDistUub = doingQpkAveDayDist(qpkPerDayUub, "UUB");
  qpkDayDistUub->Fit("gaus", ifPlotGaus, "");

  if ( ifPlot ) {
    TCanvas *c0 = canvasStyle("QpkPerDay");
    c0->SetLeftMargin(0.1);
    c0->cd();
    qpkDayDistUb->SetTitle("");
    qpkDayDistUb->SetStats(kFALSE);
    qpkDayDistUb->GetXaxis()->SetTitle("Number of Q^{pk} per day");
    qpkDayDistUb->GetXaxis()->SetTitleOffset(1.);
    qpkDayDistUb->GetXaxis()->SetTitleSize(0.06);
    qpkDayDistUb->GetXaxis()->SetLabelSize(0.05);
    qpkDayDistUb->GetYaxis()->SetTitle("Counts [au]");
    qpkDayDistUb->GetYaxis()->SetRangeUser(0., 1900.);
    qpkDayDistUb->GetYaxis()->SetTitleSize(0.06);
    qpkDayDistUb->GetYaxis()->SetLabelSize(0.05);
    qpkDayDistUb->GetYaxis()->SetTitleOffset(0.8);
    qpkDayDistUb->SetLineColor(kBlue);
    qpkDayDistUb->SetLineWidth(2);
    qpkDayDistUb->GetFunction("gaus")->SetLineColor(kBlack);
    qpkDayDistUb->Draw();
    // Plotting for UUB
    qpkDayDistUub->SetStats(kFALSE);
    qpkDayDistUub->SetLineColor(kRed);
    qpkDayDistUub->SetLineWidth(2);
    qpkDayDistUub->GetFunction("gaus")->SetLineColor(kBlack);
    qpkDayDistUub->Draw("same");
    // Drawing Legend
    TLegend *lege= new TLegend(0.58,0.45,0.9,0.95);  
    lege->AddEntry(qpkDayDistUb, "UB", "l");
    lege->AddEntry(qpkDayDistUb, Form("Entries: %.f", 
          qpkDayDistUb->GetEntries()), "");
    lege->AddEntry(qpkDayDistUb, Form("Mean: %.2e #pm %.2e", 
          qpkDayDistUb->GetMean(), qpkDayDistUb->GetMeanError()), "");
    lege->AddEntry(qpkDayDistUb, Form("Std. Dev.: %.2e #pm %.2e", 
          qpkDayDistUb->GetRMS(), qpkDayDistUb->GetRMSError()), "");
    lege->AddEntry(qpkDayDistUb, Form("#mu: %.2e", 
          qpkDayDistUb->GetFunction("gaus")->GetParameter(1)), "");
    lege->AddEntry(qpkDayDistUb, Form("#sigma: %.2e", 
          qpkDayDistUb->GetFunction("gaus")->GetParameter(2)), "");
    // Drawing legend for UUB
    lege->AddEntry(qpkDayDistUub, "UUB", "l");
    lege->AddEntry(qpkDayDistUub, Form("Entries: %.f", 
          qpkDayDistUub->GetEntries()), "");
    lege->AddEntry(qpkDayDistUub, Form("Mean: %.2e #pm %.2e", 
          qpkDayDistUub->GetMean(), qpkDayDistUub->GetMeanError()), "");
    lege->AddEntry(qpkDayDistUub, Form("Std. Dev.: %.2e #pm %.2e", 
          qpkDayDistUub->GetRMS(), qpkDayDistUub->GetRMSError()), "");
    lege->AddEntry(qpkDayDistUub, Form("#mu: %.2e", 
          qpkDayDistUub->GetFunction("gaus")->GetParameter(1)), "");
    lege->AddEntry(qpkDayDistUub, Form("#sigma: %.2e", 
          qpkDayDistUub->GetFunction("gaus")->GetParameter(2)), "");
    lege->SetTextSize(0.04);
    lege->SetBorderSize(0);
    lege->SetFillStyle(0);
    lege->Draw();
    c0->Print("../plots2/qpkDistPerDayUbUub.pdf");
  }

  retMusgm[0] = qpkDayDistUb->GetFunction("gaus")->GetParameter(1);
  retMusgm[1] = qpkDayDistUb->GetFunction("gaus")->GetParameter(2);
  retMusgm[2] = qpkDayDistUub->GetFunction("gaus")->GetParameter(1);
  retMusgm[3] = qpkDayDistUub->GetFunction("gaus")->GetParameter(2);

  return retMusgm;
}

void doMovingWindow(vector<vector<vector<int>>>nDayMonth,
    vector<vector<vector<double>>>qpkAveDay,
    vector<vector<vector<double>>>errAveDay,
    vector<vector<vector<int>>>qpkPerDay,
    vector<vector<vector<double>>> &retNormSlp, 
    vector<vector<vector<double>>> &retChi2,
    vector<vector<vector<int>>> &retFitFrsDay ) {
  const int nDaysForWindow = 7;
  // To count consecutive days for MW
  int nConsDays = 0;
  // Vectors to storage <Qpk> for MW
  vector < double > qpkDayN(7);
  vector < double > errQpkDayN(7);
  // TGraph to apply fit and get Slope and Chi2
  vector < double > xAxis{1., 2., 3., 4., 5., 6., 7.};
  TGraphErrors *distQpks;
  bool ifFitOk = false;
  double slp = 0.;
  double chi2 = 0.;
  // <Qpk> in the time window to normalize the slope
  double meanTimeWindow = 0.;
  int forCanvas = 0;
  // Consecutive days from 1st Aug.
  vector < double > dayNumbers;
  for ( double i=0.; i<122.; i++ ) 
    dayNumbers.push_back(i);

  // Moving through stations
  for ( int st_i=0; st_i<nDayMonth.size(); st_i++ )
    // Moving through PMTs
    for ( int pmt_i=0; pmt_i<3; pmt_i++ ) {
      // Moving through Days
      nConsDays = 0;
      distQpks = new TGraphErrors(nDayMonth[0][0].size(), &dayNumbers[0], 
          &qpkAveDay[st_i][pmt_i][0], 0, &errAveDay[st_i][pmt_i][0]);
      for ( int day_i=0; day_i<nDayMonth[st_i][pmt_i].size(); day_i++ ) {
        // Cutting for n-Qpk in day <Qpk> > 0
        if ( qpkPerDay[st_i][pmt_i][day_i] < 10
            || qpkAveDay[st_i][pmt_i][day_i] < 1 ) {
          nConsDays = 0;
          continue;
        }
        // Checking if there are 7 days in a row
        if ( nConsDays == nDaysForWindow ) {
          ifFitOk = distQpks->Fit("pol1","Q", "", day_i - nDaysForWindow, day_i-1);
          // Double Checking for 7 days in a row
          if ( distQpks->GetFunction("pol1")->GetNDF() != 5 ) {
            cerr << endl << " ########################## " << endl;
            cerr << "     Wrong NDF     " << endl;
            cerr << distQpks->GetFunction("pol1")->GetNDF() << endl;
            cerr << endl << " ========================== " << endl;
            exit(0);
          }
          slp = distQpks->GetFunction("pol1")->GetParameter(1);
          chi2 = distQpks->GetFunction("pol1")->GetChisquare();
          // If the movie is wanted to be printed.
          /*
          if ( st_i == 45 && pmt_i == 2 ) {
            TCanvas *tmpCanvas = canvasStyle("tmpCanvas");
            TLegend *leg = new TLegend(0.,0.2,0.6,0.5);
            tmpCanvas->cd();
            tmpCanvas->SetTopMargin(0.06);
            tmpCanvas->SetGridy();
            distQpks->SetTitle("UUB Station 1223, PMT 3");
            distQpks->GetYaxis()->SetTitle("#LTQ^{Pk}_{VEM}#GT_{7days} [FADC]");
            //distQpks->GetYaxis()->SetRangeUser(.8e2, 3.5e2);
            distQpks->GetYaxis()->SetRangeUser(.8e3, 2.4e3);
            distQpks->GetYaxis()->SetTitleSize(0.06);
            distQpks->GetYaxis()->SetLabelSize(0.05);
            distQpks->GetYaxis()->SetTitleOffset(1.);
            distQpks->GetXaxis()->SetTitle("Day");
            distQpks->GetXaxis()->SetTitleSize(0.06);
            distQpks->GetXaxis()->SetLabelSize(0.05);
            distQpks->SetLineWidth(3);
            distQpks->SetMarkerSize(1);
            distQpks->SetMarkerStyle(92);
            distQpks->SetMarkerColor(kBlue);
            distQpks->Draw("AP");
            //gPad->WaitPrimitive();
            tmpCanvas->Print(Form("../plots2/pltsForGif/weekUubSt1223pmt3%d.png",
                  forCanvas));
            forCanvas++;             
            cout << day_i << " " 
              << qpkPerDay[st_i][pmt_i][day_i-nDaysForWindow-1] << " "
              << qpkPerDay[st_i][pmt_i][day_i] << " " 
              << qpkPerDay[st_i][pmt_i][day_i+1] << endl;
            //gPad->WaitPrimitive();
          }
          */  
          // Moving to next day
          for ( int i=0; i<nDaysForWindow-1; i++ ) {
            qpkDayN[i] = qpkDayN[i+1];
            errQpkDayN[i] = errQpkDayN[i+1];
          }
          nConsDays = nDaysForWindow-1;
          // Checking if the fit was Ok to store variables
          if ( ifFitOk==0 ) {
            // Normalizing the slope
            meanTimeWindow = 0.;
            for ( int i=0; i<nDaysForWindow; i++ )
              meanTimeWindow += qpkDayN[i];
            meanTimeWindow /= nDaysForWindow;
            // Storing fitted variables
            retNormSlp[st_i][pmt_i].push_back( slp/meanTimeWindow );
            retChi2[st_i][pmt_i].push_back( chi2 );
            retFitFrsDay[st_i][pmt_i].push_back(day_i - nDaysForWindow);
          }
        }
        // Filling time window
        qpkDayN[nConsDays] = qpkAveDay[st_i][pmt_i][day_i];
        errQpkDayN[nConsDays] = errAveDay[st_i][pmt_i][day_i];
        // Counting consecutive days
        nConsDays++;
      }
    }
}

void fillingChi2VsSlop(vector<vector<vector<double>>> chi2, 
    vector<vector<vector<double>>> normslp, TH2D *hist ) {
  double pval = 0.;
  for ( int st_i=0; st_i<chi2.size(); st_i++ )
    for ( int pmt_i=0; pmt_i<3; pmt_i++ )
      for ( int val_i=0; val_i<chi2[st_i][pmt_i].size(); val_i++ ) {
        pval = TMath::Prob(chi2[st_i][pmt_i][val_i], 5);
        hist->Fill( normslp[st_i][pmt_i][val_i], log10(pval) );
      }
}

void plottingAndSaving(TString canvasName, TString outputName, TH2D *histo) {
  TCanvas *canvas = canvasStyle(canvasName); 
  gStyle->SetStatX(0.85);
  gStyle->SetStatY(0.95);
  gStyle->SetOptStat("neMR");
  gStyle->SetPalette(56);
  TString strTmp;
 
  TPaveStats *ptstats = new TPaveStats(0.6,0.55,0.95,0.98,"brNDC");
  ptstats->SetName("stats");
  ptstats->SetBorderSize(1);
  ptstats->SetFillColor(0);
  ptstats->SetTextAlign(12);
  ptstats->SetTextFont(42);
  ptstats->SetStatFormat(".2e");
  TText *ptstats_LaTex = ptstats->AddText(Form("%s", histo->GetName()));
  ptstats_LaTex->SetTextSize(0.06);
  ptstats->SetOptFit(0);
  ptstats->Draw();
  histo->GetListOfFunctions()->Add(ptstats);
  ptstats->SetParent(histo); 
  
  histo->GetYaxis()->SetTitle("Log10(Pval)");
  histo->GetXaxis()->SetTitle("Slope [FADC/day/#LTFADC#GT_{7days}]");
  histo->GetYaxis()->SetTitleSize(0.06);
  histo->GetYaxis()->SetLabelSize(0.05);
  histo->GetXaxis()->SetTitleSize(0.06);
  histo->GetXaxis()->SetLabelSize(0.05);
  //histo->GetYaxis()->SetRangeUser(-20., 1);
  //histo->GetXaxis()->SetRangeUser(-0.02, 0.03);
  histo->GetZaxis()->SetTitle("Counts [au]");
  histo->Draw("COLZ");  

  canvas->Print("../plots2/"+outputName+"2.pdf");
  //canvas->Print("../plots2/"+outputName+"AfterCuts2.pdf");
  //canvas->Close();
}

TH1D *cutByPval(bool ifUub, vector<vector<vector<int>>> fitFrsDay, 
    vector<vector<vector<double>>> chi2, 
    vector<vector<vector<double>>> normSlp, 
    vector<vector<vector<vector<double>>>> qpksInDay,
    vector<vector<vector<double>>> &retFinalQpkDist){
  TString histName = (ifUub) ? "UUB" : "UB";
  // TH1D to return the new slope distribution
  TH1D *retNormSlp = new TH1D("NormSlp"+histName, histName, 80, -0.04, 0.04);
  double pval = 0.;
  double cutNormSlpMin = (ifUub) ? -7.92e-4-3.44e-3 : -6.32e-4-3.18e-3;
  double cutNormSlpMax = (ifUub) ? -7.92e-4+3.44e-3 : -6.32e-4+3.18e-3;
  // Tmp Normalized slope to select the one closest to zero
  double tmpNormSlp = 10.;
  // Set the first day for a chosen the 7-day window
  int bestFrsDay = 0;
  // Flag for pval cut
  bool filterOk = false;
  for ( int st_i=0; st_i<fitFrsDay.size(); st_i++ ) {
    for ( int pmt_i=0; pmt_i<3; pmt_i++ ) {
      filterOk = false;
      tmpNormSlp = 10.;
      bestFrsDay = 0;
      for ( int val_i=0; val_i<fitFrsDay[st_i][pmt_i].size(); val_i++ ) {
        // First cut by Chi2 < 10        
        if ( chi2[st_i][pmt_i][val_i] > 10. )
          continue;
        // Second cut by normSlp 
        if ( normSlp[st_i][pmt_i][val_i] < cutNormSlpMin 
            || normSlp[st_i][pmt_i][val_i] > cutNormSlpMax )
          continue;          
        pval = TMath::Prob(chi2[st_i][pmt_i][val_i], 5);
        // Third cut by Pval
        if ( log10(pval) < -5. )
          continue;
        retNormSlp->Fill( normSlp[st_i][pmt_i][val_i] );
        filterOk = true;
        // Choosing the Norm. Slope closest to zero
        if ( tmpNormSlp > abs(normSlp[st_i][pmt_i][val_i]) ) {
          tmpNormSlp = abs(normSlp[st_i][pmt_i][val_i]);
          // Setting the first day for this 7-day window          
          bestFrsDay = fitFrsDay[st_i][pmt_i][val_i];
        }  
      }
      // Returning the all Qpk during the Chosen 7-days
      if ( filterOk ) {
        for ( int i=0; i<7; i++ )
          for ( auto & qpk_i : qpksInDay[st_i][pmt_i][bestFrsDay+i] )
            retFinalQpkDist[st_i][pmt_i].push_back( qpk_i );
      }
      else 
        retFinalQpkDist[st_i][pmt_i].push_back( 0 );
    }
  }
  return retNormSlp;
}

void plottingFinalSlp(TH1D *histo, bool ifUub) {
  TString printName = (ifUub) ? "UUB" : "UB"; 
  TCanvas *finalSlpCanvas = canvasStyle("Slp after Pval Cut, "+printName);
  finalSlpCanvas->cd();
  finalSlpCanvas->SetLogy(kTRUE);
 
  gStyle->SetStatX(0.85);
  gStyle->SetStatY(0.95);
  gStyle->SetOptStat("eMR");
  gStyle->SetPalette(56);

  TPaveStats *ptstats = new TPaveStats(0.6,0.55,0.95,0.98,"brNDC");      
  ptstats->SetName("stats");
  ptstats->SetBorderSize(1);
  ptstats->SetFillColor(0);
  ptstats->SetTextAlign(12);
  ptstats->SetTextFont(42);
  ptstats->SetStatFormat(".2e");
  TText *ptstats_LaTex = ptstats->AddText(Form("%s", histo->GetName()));
  ptstats_LaTex->SetTextSize(0.06);
  ptstats->SetOptFit(0);
  ptstats->Draw();
  histo->GetListOfFunctions()->Add(ptstats);
  ptstats->SetParent(histo);
  histo->GetYaxis()->SetTitle("Counts [au]");
  histo->GetXaxis()->SetTitle("Slope [FADC/day/#LTFADC#GT_{7days}]");
  histo->GetXaxis()->SetRangeUser(-0.01, 0.01);
  histo->GetYaxis()->SetTitleSize(0.06); 
  histo->GetYaxis()->SetLabelSize(0.05);
  histo->GetXaxis()->SetTitleSize(0.06);
  histo->GetXaxis()->SetLabelSize(0.05); 
  histo->Draw();
  printName = (ifUub) ? "../plots2/chi2VsSlopUubProjSlop2.pdf" 
    : "../plots2/chi2VsSlopUbProjSlop2.pdf";
  finalSlpCanvas->Print(printName);
}

void plotPmtQpksDist(TString ubUub, int st, TH1D *qpkPmt1, TH1D *qpkPmt2,
    TH1D *qpkPmt3) {
  double ymax = qpkPmt1->GetBinContent( qpkPmt1->GetMaximumBin() );
  int tmp = qpkPmt2->GetBinContent( qpkPmt2->GetMaximumBin() );
  ymax = (ymax < tmp) ? tmp : ymax;
  tmp = qpkPmt3->GetBinContent( qpkPmt3->GetMaximumBin() );
  ymax = (ymax < tmp) ? tmp : ymax;
  double xmin = (ubUub == "UB") ? 1.e2 : 1.e3;
  double xmax = (ubUub == "UB") ? 2.7e2 : 2.5e3;
  if ( st==1746 ) {
    xmin = (ubUub == "UB") ? .5e2 : 1.e3;
    xmax = (ubUub == "UB") ? 2.e2 : 2.5e3;
  }
  TCanvas *qpksDist = canvasStyle(Form("qpksDist%d",st));
  qpksDist->cd();
  qpksDist->SetLogy(kTRUE);
  qpkPmt1->SetStats(kFALSE);
  qpkPmt1->GetXaxis()->SetTitle("Q^{pk} [FADC]");
  qpkPmt1->GetXaxis()->SetRangeUser(xmin, xmax);
  qpkPmt1->GetYaxis()->SetTitle("Counts [au]");
  qpkPmt1->GetYaxis()->SetRangeUser(0.5, ymax*1.1);
  qpkPmt1->SetLineColor(kRed);
  qpkPmt1->SetLineWidth(2);
  qpkPmt1->Draw();
  qpkPmt2->SetLineColor(kBlue);
  qpkPmt2->SetLineWidth(2);
  qpkPmt2->Draw("same");
  qpkPmt3->SetLineColor(kGreen+3);
  qpkPmt3->SetLineWidth(2);
  qpkPmt3->Draw("same");
  TLegend *lgnd = new TLegend(0.7, 0.4, 0.95, 0.95);
  lgnd->AddEntry(qpkPmt1, Form(ubUub+" St. %d",st),"");
  lgnd->AddEntry(qpkPmt1, "PMT 1", "l");
  lgnd->AddEntry(qpkPmt1, Form("Entries: %.f", 
        qpkPmt1->GetEntries()), "-");
  lgnd->AddEntry(qpkPmt1, Form("Mean: %.2f #pm %.2f",
        qpkPmt1->GetMean(), qpkPmt1->GetMeanError()), "");
  lgnd->AddEntry(qpkPmt1, Form("RMS: %.2f #pm %.2f",
        qpkPmt1->GetRMS(), qpkPmt1->GetRMSError()), "");
  lgnd->AddEntry(qpkPmt2, "PMT 2", "l");
  lgnd->AddEntry(qpkPmt2, Form("Entries: %.f", 
        qpkPmt2->GetEntries()), "-");
  lgnd->AddEntry(qpkPmt2, Form("Mean: %.2f #pm %.2f",
        qpkPmt2->GetMean(), qpkPmt2->GetMeanError()), "");
  lgnd->AddEntry(qpkPmt2, Form("RMS: %.2f #pm %.2f",
        qpkPmt2->GetRMS(), qpkPmt2->GetRMSError()), "");
  lgnd->AddEntry(qpkPmt3, "PMT 3", "l");
  lgnd->AddEntry(qpkPmt3, Form("Entries: %.f", 
        qpkPmt3->GetEntries()), "-");
  lgnd->AddEntry(qpkPmt3, Form("Mean: %.2f #pm %.2f",
        qpkPmt3->GetMean(), qpkPmt3->GetMeanError()), "");
  lgnd->AddEntry(qpkPmt3, Form("RMS: %.2f #pm %.2f",
        qpkPmt3->GetRMS(), qpkPmt3->GetRMSError()), "");
  lgnd->SetTextSize(0.03);
  lgnd->SetBorderSize(0);
  lgnd->SetFillStyle(0);
  lgnd->Draw();
  qpksDist->Print(Form("../plots2/qpkdDist"+ubUub+"St%d.pdf",st));
  qpksDist->Close();
}

void plotPmtNormQpksDist(TString ubUub, int st, TH1D *qpkNorm) {
  double ymax = qpkNorm->GetBinContent( qpkNorm->GetMaximumBin() );
  double xmin = 0.7;
  double xmax = 1.4;
  TCanvas *qpksNormDist = canvasStyle(Form("qpksNormDist%d",st));
  qpksNormDist->cd();
  qpksNormDist->SetLogy(kTRUE);
  qpkNorm->SetStats(kFALSE);
  qpkNorm->GetXaxis()->SetTitle("Q^{pk}_{i, PMT}/#LTQ^{Pk}#GT_{PMT}");
  qpkNorm->GetXaxis()->SetRangeUser(xmin, xmax);
  qpkNorm->GetYaxis()->SetTitle("Counts [au]");
  qpkNorm->GetYaxis()->SetRangeUser(0.5, ymax*1.1);
  qpkNorm->SetLineColor(kBlack);
  qpkNorm->SetLineWidth(2);
  qpkNorm->Draw();
  TLegend *lgnd = new TLegend(0.7, 0.6, 0.95, 0.95);
  lgnd->AddEntry(qpkNorm, Form(ubUub+" St. %d",st),"l");
  lgnd->AddEntry(qpkNorm, Form("Entries: %.f", 
        qpkNorm->GetEntries()), "");
  lgnd->AddEntry(qpkNorm->GetFunction("gaus"), "Gauss Fit", "l");
  lgnd->AddEntry(qpkNorm, Form("#mu = %.2e #pm %.2e", 
        qpkNorm->GetFunction("gaus")->GetParameter(1), 
        qpkNorm->GetFunction("gaus")->GetParError(1)), ""); 
  lgnd->AddEntry(qpkNorm, Form("#sigma = %.2e #pm %.2e",
        qpkNorm->GetFunction("gaus")->GetParameter(2), 
        qpkNorm->GetFunction("gaus")->GetParError(2)), "");
  lgnd->SetTextSize(0.03);
  lgnd->SetBorderSize(0);
  lgnd->SetFillStyle(0);
  lgnd->Draw();
  qpksNormDist->Print(Form("../plots2/qpkdNormDist"+ubUub+"St%d.pdf",st));
  qpksNormDist->Close();
}


void doingAcc(vector<double> stIds,
    vector<vector<vector<double>>> finalQpkDistUb,
    vector<vector<vector<double>>> finalQpkDistUub, 
    vector<double> &retAccPerIdUb, vector<double> &retAccPerIdUub, 
    vector<double> &retErrAccPerIdUb, vector<double> &retErrAccPerIdUub,
    TH1D *retAccuUb, TH1D *retAccuUub, TH1D *retDistDiff, 
    vector<double> &retAccDiffPerId, vector<double> &retErrAccDiffPerId) {
  // Vector to store the <Qpk> per PMT for UB and UUB
  vector < double > avePmt(2);
  double mu = 0.;
  double errMu = 0.;
  double sgm = 0.;
  double errSgm = 0.;
  double errAcc = 0.;
  double diff = 0.;
  double diffErr = 0.;
  TString histName;
  // Flag to check if the Gaus fit was Ok
  vector < bool > ifFitOk(2);
  // Vectors to store the Qpk distributions
  vector < double > qpkDistUb;
  vector < double > qpkDistUub;
  // TH1D to plot and fit Qpk distribution
  TH1D *qpkNormUb;
  TH1D *qpkNormUub;
  TH1D *qpkPmt1Ub;
  TH1D *qpkPmt2Ub;
  TH1D *qpkPmt3Ub;
  TH1D *qpkPmt1Uub;
  TH1D *qpkPmt2Uub;
  TH1D *qpkPmt3Uub;
  // Flag to check if PMT_i was ok for UB and UUB 
  bool pmtCoinci = false;
  for ( int st_i=0; st_i<finalQpkDistUb.size(); st_i++ ) {
    pmtCoinci = false;
    qpkNormUb = new TH1D(Form("UB%d", st_i),"UB", 200, 0., 2.);
    qpkNormUub = new TH1D(Form("UUB%d",st_i), "UUB", 200, 0., 2.);
    qpkPmt1Ub = new TH1D(Form("qpkPmt1Ub%d",st_i),"", 300, 0, 300);
    qpkPmt2Ub = new TH1D(Form("qpkPmt2Ub%d",st_i),"", 300, 0, 300);;
    qpkPmt3Ub = new TH1D(Form("qpkPmt3Ub%d",st_i),"", 300, 0, 300);;
    qpkPmt1Uub = new TH1D(Form("qpkPmt1Uub%d",st_i),"", 2000, 1e3, 3e3);
    qpkPmt2Uub = new TH1D(Form("qpkPmt2Uub%d",st_i),"", 2000, 1e3, 3e3);
    qpkPmt3Uub = new TH1D(Form("qpkPmt3Uub%d",st_i),"", 2000, 1e3, 3e3);
    for ( int pmt_i=0; pmt_i<3; pmt_i++ ) {
      // Checking if pmt_i has valid Qpk dist. (see cutByPval)
      if ( finalQpkDistUb[st_i][pmt_i].size() < 2 ||
          finalQpkDistUub[st_i][pmt_i].size() < 2) {
        continue;
      }
      avePmt[0] = 0.;
      avePmt[1] = 0.;
      // Getting Qpk distribution and calculating <Qpk>
      for ( auto & qpk_i : finalQpkDistUb[st_i][pmt_i] ) {
        // Filling for individual PMTs distribution
        switch( pmt_i )  {
          case 0 :
            qpkPmt1Ub->Fill( qpk_i );
            break;
          case 1 :
            qpkPmt2Ub->Fill( qpk_i );
            break;
          case 2 :
            qpkPmt3Ub->Fill( qpk_i );
            break;
        }
        qpkDistUb.push_back( qpk_i );
        avePmt[0] += qpk_i;
      }
      avePmt[0] /= finalQpkDistUb[st_i][pmt_i].size();
      for ( auto & qpk_i : finalQpkDistUub[st_i][pmt_i] ) {
        // Filling for individual PMTs distribution
        switch( pmt_i )  {
          case 0 :
            qpkPmt1Uub->Fill( qpk_i );
            break;
          case 1 :
            qpkPmt2Uub->Fill( qpk_i );
            break;
          case 2 :
            qpkPmt3Uub->Fill( qpk_i );
            break; 
        }
        qpkDistUub.push_back( qpk_i );
        avePmt[1] += qpk_i;
      }
      avePmt[1] /= finalQpkDistUub[st_i][pmt_i].size();
      // Filling TH1D 
      for ( auto & qpk_i : qpkDistUb )
        qpkNormUb->Fill( qpk_i/avePmt[0] );
      for ( auto & qpk_i : qpkDistUub )
        qpkNormUub->Fill( qpk_i/avePmt[1] );
      pmtCoinci = true;
      qpkDistUb.clear();
      qpkDistUub.clear();
    }
    if ( !pmtCoinci )
      continue;
    // Doing Gauss fitting to calculate accuracy
    ifFitOk[0] = qpkNormUb->Fit("gaus", "Q", "", 0.4, 1.6);
    ifFitOk[1] = qpkNormUub->Fit("gaus", "Q", "", 0.4, 1.6);
    /*
    if ( st_i == 69 ) {
    TCanvas *tmpCanvas2 = canvasStyle("tmpCanvas2");
    tmpCanvas2->cd();
    qpkNormUb->Draw();
    gPad->WaitPrimitive();
    }
    */
    plotPmtQpksDist("UB", stIds[st_i], qpkPmt1Ub, qpkPmt2Ub, qpkPmt3Ub);    
    plotPmtQpksDist("UUB", stIds[st_i], qpkPmt1Uub, qpkPmt2Uub, qpkPmt3Uub);
    if ( !ifFitOk[0] && !ifFitOk[1] ) {
      plotPmtNormQpksDist("UB", stIds[st_i], qpkNormUb);
      mu = qpkNormUb->GetFunction("gaus")->GetParameter(1);
      sgm = qpkNormUb->GetFunction("gaus")->GetParameter(2);
      diff = 100.*sgm/mu;
      retAccuUb->Fill( diff );
      retAccPerIdUb[st_i] = diff;
      errMu = (1./mu)*qpkNormUb->GetFunction("gaus")->GetParError(2);
      errSgm = (sgm/(mu*mu))*qpkNormUb->GetFunction("gaus")->GetParError(1);
      errAcc = 100.*( errSgm*errSgm + errMu*errMu );
      retErrAccPerIdUb[st_i] = errAcc;
      diffErr = errAcc*errAcc;
      if ( diff > 10 ) cout << "UB: Acc > 10: " << st_i << endl;
      // Doing for UUB
      plotPmtNormQpksDist("UUB", stIds[st_i], qpkNormUub);
      mu = qpkNormUub->GetFunction("gaus")->GetParameter(1);
      sgm = qpkNormUub->GetFunction("gaus")->GetParameter(2);
      retAccuUub->Fill( 100.*sgm/mu );
      retAccPerIdUub[st_i] = 100.*sgm/mu;
      errMu = (1./mu)*qpkNormUub->GetFunction("gaus")->GetParError(2);
      errSgm = (sgm/(mu*mu))*qpkNormUub->GetFunction("gaus")->GetParError(1);
      errAcc = 100.*( errSgm*errSgm + errMu*errMu );       
      retErrAccPerIdUub[st_i] = errAcc;
      retDistDiff->Fill( diff - 100.*sgm/mu );
      retAccDiffPerId[st_i] = diff - 100.*sgm/mu;
      diffErr += errAcc*errAcc;
      retErrAccDiffPerId[st_i] = sqrt(diffErr);
      if ( (100.*sgm/mu) > 10. ) cout << "UUB: Acc > 10: " << st_i << " " << 
        stIds[st_i] << " " << 100.*sgm/mu << endl;
    }
    else {
      retAccuUb->Fill(-1.);
      retAccuUub->Fill(-1);
      retErrAccPerIdUub[st_i] = -1.;
    }
    qpkNormUb->Reset();
    qpkNormUb->Delete();
    qpkNormUub->Reset();
    qpkNormUub->Delete();
    qpkPmt1Ub->Reset();
    qpkPmt1Ub->Delete();
    qpkPmt2Ub->Reset();
    qpkPmt2Ub->Delete();
    qpkPmt3Ub->Reset();
    qpkPmt3Ub->Delete();
  }
}

void fillingPerId(vector<double> accVec, vector<double> errAccVec,
    vector<double> stIdVec, TH1D *retHist) {
  for ( int acc_i=0; acc_i<accVec.size(); acc_i++ ) {
    retHist->SetBinContent(acc_i+1, accVec[acc_i]);
    retHist->SetBinError(acc_i+1, errAccVec[acc_i]);
    retHist->GetXaxis()->SetBinLabel(acc_i+1, Form("%.f", stIdVec[acc_i]));
  }
}

void makingQpkAccuracy() {

  // output file for storing results
  TFile *outFile = new TFile("makingQpkAccuracy.root","RECREATE");
 
  // Reading IDs stations from list
  ifstream stationSelected;
  //TString strFileSelecSt = "../halfList.txt";
  TString strFileSelecSt = "../fullUubStationsListVert.txt";
  int readStId = 0;
  vector < int > stId;
  stationSelected.open(strFileSelecSt);
  while( stationSelected.good() ) {
    stationSelected >> readStId;
    stId.push_back( readStId );
  }
  stId.pop_back();
  stationSelected.close();

  // Vector for store Qpk and time data
  vector < vector < double > > qpkUb(3);
  vector < vector < double > > timeUb(3);
  vector < vector < double > > qpkUub(3);
  vector < vector < double > > timeUub(3);

  // Vectors for store qpkAveDay and its errors
  //
  // Day of the month, starting from 1st of August
  vector < vector < vector < int > > > nDayMonthUb(stId.size());
  vector < vector < vector < int > > > nDayMonthUub(stId.size());
  // Number of Qpk in a single day
  vector < vector < vector < int > > > qpkPerDayUb(stId.size());
  vector < vector < vector < int > > > qpkPerDayUub(stId.size());
  // Average Qpk per day
  vector < vector < vector < double > > > qpkAveDayUb(stId.size());
  vector < vector < vector < double > > > qpkAveDayUub(stId.size());
  vector < vector < vector < double > > > errAveDayUb(stId.size());
  vector < vector < vector < double > > > errAveDayUub(stId.size());
  // All the Qpk's fitted in a single day
  vector < vector < vector < vector < double > > > > qpksInDayUb(stId.size());
  vector < vector < vector < vector < double > > > > qpksInDayUub(stId.size());

  // Vector to store Station IDs.
  vector < double > vecStID;
   
  // Counter for the stations readied
  int cntSts = 0;
  // Days from 1st of Aug to 30th of Nov
  const int daysAugNov = 122;

  // Reading the Qpk fitted per station, per PMT
  // and returning <Qpk> per day, etc
  for ( auto & st_i : stId ) {
    cout << "Doing Qpk temporal series for station " << st_i << endl;
    // Resizing by PMTs
    nDayMonthUb[cntSts].resize(3);
    qpkAveDayUb[cntSts].resize(3);
    errAveDayUb[cntSts].resize(3);
    qpkPerDayUb[cntSts].resize(3);
    qpksInDayUb[cntSts].resize(3);
    nDayMonthUub[cntSts].resize(3);
    qpkAveDayUub[cntSts].resize(3);
    errAveDayUub[cntSts].resize(3);
    qpkPerDayUub[cntSts].resize(3);
    qpksInDayUub[cntSts].resize(3);
    vecStID.push_back(st_i);
    // Fetching Qpk and time values
    for ( int pmt_i=1; pmt_i<4; pmt_i++ ) {
      // Resizing by n days from Aug to Nov: 122
      nDayMonthUb[cntSts][pmt_i-1].resize(daysAugNov);
      qpkAveDayUb[cntSts][pmt_i-1].resize(daysAugNov);
      errAveDayUb[cntSts][pmt_i-1].resize(daysAugNov);
      qpkPerDayUb[cntSts][pmt_i-1].resize(daysAugNov);
      qpksInDayUb[cntSts][pmt_i-1].resize(daysAugNov);
      nDayMonthUub[cntSts][pmt_i-1].resize(daysAugNov);
      qpkAveDayUub[cntSts][pmt_i-1].resize(daysAugNov);
      errAveDayUub[cntSts][pmt_i-1].resize(daysAugNov);
      qpkPerDayUub[cntSts][pmt_i-1].resize(daysAugNov);
      qpksInDayUub[cntSts][pmt_i-1].resize(daysAugNov);
      //if ( st_i != 1746 )
        //continue;
      // Filling vectors with the Qpk fitted 
      fillQpkTimeVals(false, pmt_i, st_i, qpkUb[pmt_i-1], timeUb[pmt_i-1]);
      fillQpkTimeVals(true, pmt_i, st_i, qpkUub[pmt_i-1], timeUub[pmt_i-1]);
      // Building vectors for <Qpk>day, time, etc
      doAvePerTime(false, qpkUb[pmt_i-1], timeUb[pmt_i-1], 
          nDayMonthUb[cntSts][pmt_i-1], qpkAveDayUb[cntSts][pmt_i-1], 
          errAveDayUb[cntSts][pmt_i-1], qpkPerDayUb[cntSts][pmt_i-1], 
          qpksInDayUb[cntSts][pmt_i-1]);
      doAvePerTime(true, qpkUub[pmt_i-1], timeUub[pmt_i-1], 
          nDayMonthUub[cntSts][pmt_i-1], qpkAveDayUub[cntSts][pmt_i-1], 
          errAveDayUub[cntSts][pmt_i-1], qpkPerDayUub[cntSts][pmt_i-1],
          qpksInDayUub[cntSts][pmt_i-1]);
    }
    cntSts++;
    qpkUb.clear();
    qpkUub.clear();
    timeUb.clear();
    timeUub.clear();
    qpkUb.resize(3);
    qpkUub.resize(3);
    timeUb.resize(3);
    timeUub.resize(3);
  }

  // Doing Distribution for <Qpk>_day 
  // and returning fitting parameters
  vector < double > muSgmForCutAveDay 
    = doPlotDistQpkPerDay(qpkPerDayUb, qpkPerDayUub, false);
  
  // Applying Moving-Window
  // Vectors to store fitted information per Station per PMT
  vector < vector < vector < double > > > normSlpUb(stId.size());
  vector < vector < vector < double > > > normSlpUub(stId.size());
  vector < vector < vector < double > > > chi2Ub(stId.size());
  vector < vector < vector < double > > > chi2Uub(stId.size());
  // Vector to store the first day of the window
  vector < vector < vector < int > > > fitFrsDayUb(stId.size());
  vector < vector < vector < int > > > fitFrsDayUub(stId.size());
  // Vector to store the all the Qpk in during the 7 days
  vector < vector < vector < double > > > finalQpkDistUb(stId.size());
  vector < vector < vector < double > > > finalQpkDistUub(stId.size());
  for (int i=0; i<stId.size(); i++) {
    normSlpUb[i].resize(3);
    normSlpUub[i].resize(3);
    chi2Ub[i].resize(3);
    chi2Uub[i].resize(3);
    fitFrsDayUb[i].resize(3);
    fitFrsDayUub[i].resize(3);
    finalQpkDistUb[i].resize(3);
    finalQpkDistUub[i].resize(3);
  }
  doMovingWindow(nDayMonthUb, qpkAveDayUb, errAveDayUb, qpkPerDayUb,
      normSlpUb, chi2Ub, fitFrsDayUb);
  doMovingWindow(nDayMonthUub, qpkAveDayUub, errAveDayUub, qpkPerDayUub, 
      normSlpUub, chi2Uub, fitFrsDayUub);
  
  // TH2 to plot Log10(Pval) chi2 and NormSlope 
  int nBinsX = 70; // bins for Norm. Slope's X axis
  double xLow = -0.03;
  double xUp = 0.04;
  int nBinsY = 62; // bins for Log10(Pval) Y axis
  double yLow = -30.;
  double yUp = 1;
  TH2D *chi2VsSlopUb = new TH2D("chi2VsSlopUb", "", 
      nBinsX, xLow, xUp, nBinsY, yLow, yUp);
  TH2D *chi2VsSlopUub = new TH2D("chi2VsSlopUub", "", 
      nBinsX, xLow, xUp, nBinsY, yLow, yUp);
  // Filling and plotting TH2 plots
  fillingChi2VsSlop( chi2Ub, normSlpUb, chi2VsSlopUb );
  fillingChi2VsSlop( chi2Uub, normSlpUub, chi2VsSlopUub );   
  plottingAndSaving("UB", "chi2VsSlopUb", chi2VsSlopUb);
  plottingAndSaving("UUB", "chi2VsSlopUub", chi2VsSlopUub);

  // Applying cut by Log10(Pval) and chosen the best slope.
  // The Qpks values in the chosen 7 days is returned
  TH1D *normSlpDistUb = cutByPval(false, fitFrsDayUb, chi2Ub, normSlpUb, 
      qpksInDayUb, finalQpkDistUb); 
  TH1D *normSlpDistUub = cutByPval(true, fitFrsDayUub, chi2Uub, normSlpUub, 
      qpksInDayUub, finalQpkDistUub);
  // Plotting the normalized slope dist after cut
  plottingFinalSlp(normSlpDistUb, false);
  plottingFinalSlp(normSlpDistUub, true);  

  // Doing Accuracy calculations  
  // TH1 histos to store the respective results
  TH1D *accDistUb = new TH1D("accDistUb", "", 200, 0, 20);
  vector < double > accPerIdUb(stId.size());
  vector < double > errAccPerIdUb(stId.size());
  TH1D *accDistUub = new TH1D("accDistUub", "", 200, 0, 20);
  vector < double > accPerIdUub(stId.size());
  vector < double > errAccPerIdUub(stId.size()); 

  TH1D *distDiff = new TH1D("distDiff", "", 80, -4, 4);
  vector < double > accDiffPerId(stId.size());
  vector < double > errAccDiffPerId(stId.size());

  doingAcc(vecStID, finalQpkDistUb, finalQpkDistUub, accPerIdUb, 
      accPerIdUub, errAccPerIdUb, errAccPerIdUub, accDistUb, accDistUub, 
      distDiff, accDiffPerId, errAccDiffPerId); 

  // Plotting Accuracy distribution
  TCanvas *accCanvas = canvasStyle("accCanvas");
  accCanvas->cd();
  accDistUb->SetStats(kFALSE);
  accDistUb->SetLineColor(kBlue);
  accDistUb->GetXaxis()->SetRangeUser(0., 5.);
  accDistUb->GetXaxis()->SetTitle("#sigma/#mu [%]");  
  accDistUb->GetXaxis()->SetTitleSize(0.06);
  accDistUb->GetXaxis()->SetLabelSize(0.05);
  //accDistUb->GetYaxis()->SetRangeUser(0., 9.5);
  accDistUb->GetYaxis()->SetTitle("Counts [au]");
  accDistUb->GetYaxis()->SetTitleSize(0.06);
  accDistUb->GetYaxis()->SetLabelSize(0.05);
  accDistUb->SetLineWidth(2);
  accDistUb->SetFillStyle(3345);
  accDistUb->SetFillColor(kBlue);
  accDistUb->Draw();
  accDistUub->SetLineColor(kRed);
  accDistUub->SetLineWidth(2);
  accDistUub->SetFillStyle(3354);
  accDistUub->SetFillColor(kRed);
  accDistUub->Draw("same");
  TLegend *leg = new TLegend(0.59,0.6,0.9,0.95);  
  leg->AddEntry(accDistUb, "UB", "l" );
  leg->AddEntry(accDistUb, Form("Mean: %.2f %% #pm %.2f %%",
        accDistUb->GetMean(), accDistUb->GetMeanError()), "");
  leg->AddEntry(accDistUb, Form("RMS: %.2f %% #pm %.2f %%",accDistUb->GetRMS(),
        accDistUb->GetRMSError()),"");
  leg->AddEntry(accDistUub, "UUB", "l" );
  leg->AddEntry(accDistUub, Form("Mean: %.2f %% #pm %.2f %%",
        accDistUub->GetMean(), accDistUub->GetMeanError()), "");
  leg->AddEntry(accDistUub, Form("RMS: %.2f %% #pm %.2f %%",accDistUub->GetRMS(),
        accDistUub->GetRMSError()) ,"");
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();
  accCanvas->Print("../plots2/accQpkFitUbUubAllStAllPmt_Stats2.pdf");
   
  // Plotting accuracy per Station ID 
  TCanvas *accPerIdCanvas = canvasStyle("accPerIdCanvas");
  accPerIdCanvas->cd();
  TH1D *histAccPerIdUb = new TH1D("histAcccPerIdUb","UB", vecStID.size(), 0, 
      vecStID.size());
  TH1D *histAccPerIdUub = new TH1D("histAccPerIdUub","UUB", vecStID.size(), 0, 
      vecStID.size());
  fillingPerId(accPerIdUb, errAccPerIdUb, vecStID, histAccPerIdUb);
  fillingPerId(accPerIdUub, errAccPerIdUub, vecStID, histAccPerIdUub);
  
  TLine *lineAveAccUb;
  TLine *lineAveAccUub;
  histAccPerIdUb->SetTitle("");
  histAccPerIdUb->SetStats(kFALSE);
  histAccPerIdUb->GetXaxis()->SetTitle("Station ID.");
  histAccPerIdUb->GetXaxis()->SetTitleOffset(1.5); 
  histAccPerIdUb->GetXaxis()->SetTitleSize(0.05);
  histAccPerIdUb->GetXaxis()->SetLabelSize(0.04);
  histAccPerIdUb->GetXaxis()->SetTitleOffset(1.1);
  histAccPerIdUb->GetYaxis()->SetRangeUser(0., 5.);
  histAccPerIdUb->GetYaxis()->SetTitle("#sigma/#mu [%]");
  histAccPerIdUb->GetYaxis()->SetTitleSize(0.06);
  histAccPerIdUb->GetYaxis()->SetLabelSize(0.05);
  histAccPerIdUb->SetMarkerStyle(71);
  histAccPerIdUb->SetMarkerColor(kBlue);
  histAccPerIdUb->SetLineColor(kBlue);
  histAccPerIdUb->SetMarkerSize(1.2);
  histAccPerIdUb->Draw("E1");
  lineAveAccUb = new TLine(0.2, accDistUb->GetMean(), 74.5, accDistUb->GetMean());
  lineAveAccUb->SetLineWidth(2);
  lineAveAccUb->SetLineColor(kBlue);
  lineAveAccUb->Draw();

  histAccPerIdUub->SetMarkerStyle(73);
  histAccPerIdUub->SetMarkerColor(kRed);
  histAccPerIdUub->SetLineColor(kRed);
  histAccPerIdUub->SetMarkerSize(1.2);
  histAccPerIdUub->Draw("E1 same");
  lineAveAccUb = new TLine(0.2, accDistUub->GetMean(), 74.5, accDistUub->GetMean());
  lineAveAccUb->SetLineWidth(2);
  lineAveAccUb->SetLineColor(kRed);
  lineAveAccUb->Draw(); 
  
  leg = new TLegend(0.4,0.7,0.95,0.95);
  leg->AddEntry(histAccPerIdUb,Form("UB Ave.: %.2f %% #pm %.2f %%",
        accDistUb->GetMean(), accDistUb->GetMeanError()), "l");
  leg->AddEntry(histAccPerIdUub,Form("UUB Ave.: %.2f %% #pm %.2f %%",
        accDistUub->GetMean(), accDistUub->GetMeanError()), "l");
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();
  accPerIdCanvas->Print("../plots2/accQpkFitUbUubPerSt_Stats2.pdf");

  // Plotting Diff
  TCanvas *diffCanvas = canvasStyle("diffCanvas");
  diffCanvas->cd();

  distDiff->SetTitle("");
  distDiff->SetStats(kFALSE);
  distDiff->GetXaxis()->SetTitle("#left(#sigma/#mu#right)_{UB} - #left(#sigma/#mu#right)_{UUB} [%]");
  distDiff->GetXaxis()->SetTitleOffset(1.3);
  distDiff->GetXaxis()->SetRangeUser(-3., 3.);
  distDiff->GetXaxis()->SetTitleSize(0.06);
  distDiff->GetXaxis()->SetLabelSize(0.05);
  distDiff->GetXaxis()->SetTitleOffset(1.);
  distDiff->GetYaxis()->SetTitle("Counts [au]");
  distDiff->GetYaxis()->SetTitleSize(0.06); 
  distDiff->GetYaxis()->SetLabelSize(0.05);
  distDiff->GetYaxis()->SetRangeUser(0., 12.);
  distDiff->SetLineColor(kGreen+3);
  distDiff->Draw();

  leg = new TLegend(0.5,0.78,0.98,0.9);
  leg->AddEntry(distDiff, Form("MEAN: %.3f %% #pm %.3f %%", distDiff->GetMean(), 
        distDiff->GetMeanError()), "");
  leg->AddEntry(distDiff, Form("RMS:     %.3f %% #pm %.3f %%", distDiff->GetRMS(), 
      distDiff->GetRMSError()), ""); 
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();

  lineAveAccUb = new TLine(0.,0.,0.,17.);
  lineAveAccUb->SetLineWidth(2);
  lineAveAccUb->SetLineColor(kGray);
  lineAveAccUb->SetLineStyle(2);
  lineAveAccUb->Draw();
  diffCanvas->Print("../plots2/accQpkFitUbUubDistDiff_Stats2.pdf");

  TCanvas *diffPerIdCanvas = canvasStyle("diffPerIdCanvas");
  diffPerIdCanvas->cd();
  TH1D *histDiffPerId = new TH1D("histDiffPerId","", vecStID.size(), 0,
      vecStID.size());
  fillingPerId(accDiffPerId, errAccDiffPerId, vecStID, histDiffPerId);

  histDiffPerId->SetTitle("");
  histDiffPerId->SetStats(kFALSE);                       
  histDiffPerId->GetXaxis()->SetTitle("Station ID.");
  histDiffPerId->GetXaxis()->SetTitleOffset(1.1);
  histDiffPerId->GetXaxis()->SetTitleSize(0.05);
  histDiffPerId->GetXaxis()->SetLabelSize(0.04);
  histDiffPerId->GetYaxis()->SetRangeUser(-2., 2.);
  histDiffPerId->GetYaxis()->SetTitleSize(0.06);
  histDiffPerId->GetYaxis()->SetLabelSize(0.05);
  histDiffPerId->GetYaxis()->SetTitleOffset(1.0);
  histDiffPerId->GetYaxis()->SetTitle("#left(#sigma/#mu#right)_{UB} - #left(#sigma/#mu#right)_{UUB} [%]");
  histDiffPerId->SetMarkerStyle(71);
  histDiffPerId->SetMarkerColor(kGreen+3);
  histDiffPerId->SetLineColor(kGreen+3);
  histDiffPerId->SetMarkerSize(1.2);
  histDiffPerId->Draw("P");
 
  lineAveAccUub = new TLine(0.5,0.,74.5,0.);
  lineAveAccUub->SetLineWidth(2);
  lineAveAccUub->SetLineColor(kGray);
  lineAveAccUub->SetLineStyle(2);
  lineAveAccUub->Draw();

  diffPerIdCanvas->Print("../plots2/accQpkFitUbUubDiffPerSt_Stats2.pdf");
  //exit(0);
}
