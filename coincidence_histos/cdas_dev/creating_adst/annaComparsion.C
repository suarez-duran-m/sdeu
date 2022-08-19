double getMean(vector<double> vect);
double getRMS(vector<double> vect, double mean);
void setCanvasStyle(TCanvas& canvas);
void plotRelQpk(vector<double> grp, vector<double> grpRms,
    vector<double> labels, int sizeGrp, int pmt, int grpN, TString name, 
    int isCQpk);
void plotCoinciFactor(TH1D &kat, TH1D &mao, TString pmt);

void annaComparsion() {
  // Fetching data files
  ifstream dataMao("qpksForDeltas.dat");  
  ifstream dataKat("cch_fits.csv");

  // Defining vectors to store Qpk, CQpk and their errors
  // Mao's  
  vector < vector < vector < double > > > qpkStMao(1900);
  vector < vector < vector < double > > > qpkErrStMao(1900);
  vector < vector < vector < double > > > cqpkStMao(1900);
  vector < vector < vector < double > > > cqpkErrStMao(1900);
  vector < vector < vector < int > > > timeStMao(1900);
  // Katarina's    
  vector < vector < vector < double > > > qpkStKat(1900);
  vector < vector < vector < double > > > qpkErrStKat(1900);
  vector < vector < vector < double > > > cqpkStKat(1900);
  vector < vector < vector < double > > > cqpkErrStKat(1900);
  vector < vector < vector < int > > > timeStKat(1900);  
  // For plotting
  vector < vector < vector < double > > > qpkRel(1900);
  vector < vector < vector < double > > > cqpkRel(1900);
  // Re-sizing per PMT
  for ( int i=0; i<1900; i++ ) {
    qpkStMao[i].resize(3);
    qpkErrStMao[i].resize(3);
    cqpkStMao[i].resize(3);
    cqpkErrStMao[i].resize(3);
    timeStMao[i].resize(3);    
    qpkStKat[i].resize(3);
    qpkErrStKat[i].resize(3);
    cqpkStKat[i].resize(3);
    cqpkErrStKat[i].resize(3);
    timeStKat[i].resize(3);
    qpkRel[i].resize(3);
    cqpkRel[i].resize(3);
  }  
  
  // Temporal variables to fetch data
  double tmpStId = 0.;
  double tmpPmtId = 0.;
  double tmpQpk = 0.;
  double tmpQpkErr = 0.;
  double tmpCQpk = 0.;
  double tmpCQpkErr = 0.;
  double tmpTime = 0.;   
  
  // Reading and storing data from Mao's  
  while ( dataMao >> tmpStId >> tmpPmtId >> tmpTime
      >> tmpQpk >> tmpQpkErr >> tmpCQpk >> tmpCQpkErr ) {
    qpkStMao[int(tmpStId)][int(tmpPmtId-1)].push_back( tmpQpk );
    qpkErrStMao[int(tmpStId)][int(tmpPmtId-1)].push_back( tmpQpkErr );
    cqpkStMao[int(tmpStId)][int(tmpPmtId-1)].push_back( tmpCQpk );
    cqpkErrStMao[int(tmpStId)][int(tmpPmtId-1)].push_back( tmpCQpkErr );
    timeStMao[int(tmpStId)][int(tmpPmtId-1)].push_back( tmpTime );
  }  
  dataMao.close();
  // Reading and storing data from Kat's  
  while ( dataKat >> tmpStId >> tmpPmtId >> tmpTime 
      >> tmpCQpk >> tmpCQpkErr >> tmpQpk >> tmpQpkErr ) {
    qpkStKat[int(tmpStId)][int(tmpPmtId-1)].push_back( tmpQpk );
    qpkErrStKat[int(tmpStId)][int(tmpPmtId-1)].push_back( tmpQpkErr );
    cqpkStKat[int(tmpStId)][int(tmpPmtId-1)].push_back( tmpCQpk );
    cqpkErrStKat[int(tmpStId)][int(tmpPmtId-1)].push_back( tmpCQpkErr );
    timeStKat[int(tmpStId)][int(tmpPmtId-1)].push_back( tmpTime );
  }  
  dataKat.close();

  vector < vector < double > > qpkRelMean(3);
  vector < vector < double > > qpkRelRms(3);
  vector < vector < double > > cqpkRelMean(3);
  vector < vector < double > > cqpkRelRms(3);
  vector < vector < double > > stIds(3);
  vector < vector < double > > stIdsLabel(3);
  TH1D coinciFactorKatPmt1("coinciFactorKatPmt1", "", 300, 0, 3);
  TH1D coinciFactorKatPmt2("coinciFactorKatPmt2", "", 300, 0, 3);
  TH1D coinciFactorKatPmt3("coinciFactorKatPmt3", "", 300, 0, 3);
  TH1D coinciFactorMaoPmt1("coinciFactorMaoPmt1", "", 300, 0, 3);
  TH1D coinciFactorMaoPmt2("coinciFactorMaoPmt2", "", 300, 0, 3);
  TH1D coinciFactorMaoPmt3("coinciFactorMaoPmt3", "", 300, 0, 3);
  vector < int > cntStRelQpk(3);
  double tmpMean = 0.;
  // Building RelErr vs Station
  for ( int st_i=0; st_i<1900; st_i++ ) {
    for ( int pmt_i=0; pmt_i<3; pmt_i++ ) {      
      if ( timeStMao[st_i][pmt_i].size() < 1 || timeStKat[st_i][pmt_i].size() < 1 )
        continue;
      for ( int utcMao = 0; utcMao < timeStMao[st_i][pmt_i].size(); utcMao++ ) {
        for ( int utcKat = 0; utcKat < timeStKat[st_i][pmt_i].size(); utcKat++ ) {
          if ( timeStMao[st_i][pmt_i][utcMao] == timeStKat[st_i][pmt_i][utcKat] ) {
            qpkRel[st_i][pmt_i].push_back( 
                100.*(qpkStKat[st_i][pmt_i][utcKat] - qpkStMao[st_i][pmt_i][utcMao])
                / qpkStMao[st_i][pmt_i][utcMao] );
            cqpkRel[st_i][pmt_i].push_back( 
                100.*(cqpkStKat[st_i][pmt_i][utcKat] - cqpkStMao[st_i][pmt_i][utcMao])
                / cqpkStMao[st_i][pmt_i][utcMao] );
            switch( pmt_i ) {
              case 0:
                coinciFactorKatPmt1.Fill( cqpkStKat[st_i][pmt_i][utcKat] / qpkStKat[st_i][pmt_i][utcKat] );
                coinciFactorMaoPmt1.Fill( cqpkStMao[st_i][pmt_i][utcMao] / qpkStMao[st_i][pmt_i][utcMao] );
                break;
              case 1:
                coinciFactorKatPmt2.Fill( cqpkStKat[st_i][pmt_i][utcKat] / qpkStKat[st_i][pmt_i][utcKat] );
                coinciFactorMaoPmt2.Fill( cqpkStMao[st_i][pmt_i][utcMao] / qpkStMao[st_i][pmt_i][utcMao] );
                break;
              case 2:
                coinciFactorKatPmt3.Fill( cqpkStKat[st_i][pmt_i][utcKat] / qpkStKat[st_i][pmt_i][utcKat] );
                coinciFactorMaoPmt3.Fill( cqpkStMao[st_i][pmt_i][utcMao] / qpkStMao[st_i][pmt_i][utcMao] );
                break;
            
            }         }
        }
      }
      if ( qpkRel[st_i][pmt_i].size() > 0 ) {
        cntStRelQpk[pmt_i]++;
        // Doing for Qpk
        tmpMean = getMean( qpkRel[st_i][pmt_i] );
        qpkRelMean[pmt_i].push_back( tmpMean );
        qpkRelRms[pmt_i].push_back( getRMS( qpkRel[st_i][pmt_i], tmpMean ) );
        // Doing for CQpk
        tmpMean = getMean( cqpkRel[st_i][pmt_i] );
        cqpkRelMean[pmt_i].push_back( tmpMean );
        cqpkRelRms[pmt_i].push_back( getRMS( cqpkRel[st_i][pmt_i], tmpMean ) );

        stIds[pmt_i].push_back( cntStRelQpk[pmt_i]-1 );
        stIdsLabel[pmt_i].push_back( st_i );
      }
    }
  }

  // Plotting results for RelQpk
  // Creating groups for plotting
  // relQpkPmtGrp[grp_i][st_i]
  vector < vector < double > > relQpkPmt1Grp(3);
  vector < vector < double > > relQpkPmt2Grp(3);
  vector < vector < double > > relQpkPmt3Grp(3);
  vector < vector < double > > relQpkRmsPmt1Grp(3);
  vector < vector < double > > relQpkRmsPmt2Grp(3);
  vector < vector < double > > relQpkRmsPmt3Grp(3);

  vector < vector < double > > relCQpkPmt1Grp(3);
  vector < vector < double > > relCQpkPmt2Grp(3);
  vector < vector < double > > relCQpkPmt3Grp(3);
  vector < vector < double > > relCQpkRmsPmt1Grp(3);
  vector < vector < double > > relCQpkRmsPmt2Grp(3);
  vector < vector < double > > relCQpkRmsPmt3Grp(3);
 
  for ( int grp_i=0; grp_i<3; grp_i++ ) {
    if ( grp_i < 2 ) {
      relQpkPmt1Grp[grp_i].resize(50);
      relQpkPmt2Grp[grp_i].resize(50);
      relQpkPmt3Grp[grp_i].resize(50);
      relQpkRmsPmt1Grp[grp_i].resize(50);
      relQpkRmsPmt2Grp[grp_i].resize(50);
      relQpkRmsPmt3Grp[grp_i].resize(50);

      relCQpkPmt1Grp[grp_i].resize(50);    
      relCQpkPmt2Grp[grp_i].resize(50);
      relCQpkPmt3Grp[grp_i].resize(50);    
      relCQpkRmsPmt1Grp[grp_i].resize(50);
      relCQpkRmsPmt2Grp[grp_i].resize(50);
      relCQpkRmsPmt3Grp[grp_i].resize(50);
    }
    
    if ( grp_i == 2 ) {
      relQpkPmt1Grp[grp_i].resize(stIds[0].size()-100);
      relQpkPmt2Grp[grp_i].resize(stIds[1].size()-100);
      relQpkPmt3Grp[grp_i].resize(stIds[2].size()-100);
      relQpkRmsPmt1Grp[grp_i].resize(stIds[0].size()-100);
      relQpkRmsPmt2Grp[grp_i].resize(stIds[1].size()-100);
      relQpkRmsPmt3Grp[grp_i].resize(stIds[2].size()-100);

      relCQpkPmt1Grp[grp_i].resize(stIds[0].size()-100);
      relCQpkPmt2Grp[grp_i].resize(stIds[1].size()-100);
      relCQpkPmt3Grp[grp_i].resize(stIds[2].size()-100);
      relCQpkRmsPmt1Grp[grp_i].resize(stIds[0].size()-100);
      relCQpkRmsPmt2Grp[grp_i].resize(stIds[1].size()-100);
      relCQpkRmsPmt3Grp[grp_i].resize(stIds[2].size()-100);      
    }
  }

  for ( int pmt_i=0; pmt_i<3; pmt_i++ ) {
    for ( int st_i=0; st_i<stIds[pmt_i].size(); st_i++ ) {
      switch ( pmt_i ) {
        case 0:
          if ( st_i < 50 ) {
            relQpkPmt1Grp[0][st_i] = qpkRelMean[pmt_i][st_i];
            relQpkRmsPmt1Grp[0][st_i] = qpkRelRms[pmt_i][st_i];
            relCQpkPmt1Grp[0][st_i] = cqpkRelMean[pmt_i][st_i];
            relCQpkRmsPmt1Grp[0][st_i] = cqpkRelRms[pmt_i][st_i];
          }
          if ( st_i >= 50 && st_i < 100 ) {
            relQpkPmt1Grp[1][st_i-50] = qpkRelMean[pmt_i][st_i];
            relQpkRmsPmt1Grp[1][st_i-50] = qpkRelRms[pmt_i][st_i];
            relCQpkPmt1Grp[1][st_i-50] = cqpkRelMean[pmt_i][st_i];
            relCQpkRmsPmt1Grp[1][st_i-50] = cqpkRelRms[pmt_i][st_i];
          }
          if ( st_i >= 100 ) {
            relQpkPmt1Grp[2][st_i-100] = qpkRelMean[pmt_i][st_i];
            relQpkRmsPmt1Grp[2][st_i-100] = qpkRelRms[pmt_i][st_i];
            relCQpkPmt1Grp[2][st_i-100] = cqpkRelMean[pmt_i][st_i];
            relCQpkRmsPmt1Grp[2][st_i-100] = cqpkRelRms[pmt_i][st_i];
          }
          break;
        case 1:
          if ( st_i < 50 ) {
            relQpkPmt2Grp[0][st_i] = qpkRelMean[pmt_i][st_i];
            relQpkRmsPmt2Grp[0][st_i] = qpkRelRms[pmt_i][st_i];
            relCQpkPmt2Grp[0][st_i] = cqpkRelMean[pmt_i][st_i];
            relCQpkRmsPmt2Grp[0][st_i] = cqpkRelRms[pmt_i][st_i];
          }
          if ( st_i >= 50 && st_i < 100 ) {
            relQpkPmt2Grp[1][st_i-50] = qpkRelMean[pmt_i][st_i];
            relQpkRmsPmt2Grp[1][st_i-50] = qpkRelRms[pmt_i][st_i];
            relCQpkPmt2Grp[1][st_i-50] = cqpkRelMean[pmt_i][st_i];
            relCQpkRmsPmt2Grp[1][st_i-50] = cqpkRelRms[pmt_i][st_i];
          }
          
          if ( st_i >= 100 ) {
            relQpkPmt2Grp[2][st_i-100] = qpkRelMean[pmt_i][st_i];
            relQpkRmsPmt2Grp[2][st_i-100] = qpkRelRms[pmt_i][st_i];
            relCQpkPmt2Grp[2][st_i-100] = cqpkRelMean[pmt_i][st_i];
            relCQpkRmsPmt2Grp[2][st_i-100] = cqpkRelRms[pmt_i][st_i];
          }
          break;
          
        case 2:
          if ( st_i < 50 ) {
            relQpkPmt3Grp[0][st_i] = qpkRelMean[pmt_i][st_i];
            relQpkRmsPmt3Grp[0][st_i] = qpkRelRms[pmt_i][st_i];
            relCQpkPmt3Grp[0][st_i] = cqpkRelMean[pmt_i][st_i];
            relCQpkRmsPmt3Grp[0][st_i] = cqpkRelRms[pmt_i][st_i];
          }
          if ( st_i >= 50 && st_i < 100 ) {
            relQpkPmt3Grp[1][st_i-50] = qpkRelMean[pmt_i][st_i];
            relQpkRmsPmt3Grp[1][st_i-50] = qpkRelRms[pmt_i][st_i];
            relCQpkPmt3Grp[1][st_i-50] = cqpkRelMean[pmt_i][st_i];
            relCQpkRmsPmt3Grp[1][st_i-50] = cqpkRelRms[pmt_i][st_i];
          }
          
          if ( st_i >= 100 ) {
            relQpkPmt3Grp[2][st_i-100] = qpkRelMean[pmt_i][st_i];
            relQpkRmsPmt3Grp[2][st_i-100] = qpkRelRms[pmt_i][st_i];
            relCQpkPmt3Grp[2][st_i-100] = cqpkRelMean[pmt_i][st_i];
            relCQpkRmsPmt3Grp[2][st_i-100] = cqpkRelRms[pmt_i][st_i];
          }
          break;
      }     
    }
  }

  // For PMT1
  vector < vector < double > > labelsPmt1Grp(3);
  for ( int i=0; i<50; i++ ) {
    labelsPmt1Grp[0].push_back( stIdsLabel[0][i] );
    labelsPmt1Grp[1].push_back( stIdsLabel[0][i+50] );
  }
  for ( int i=0; i<stIds[0].size()-100; i++ )
    labelsPmt1Grp[2].push_back( stIdsLabel[0][i+100] );

  plotRelQpk(relQpkPmt1Grp[0], relQpkRmsPmt1Grp[0], 
      labelsPmt1Grp[0], 50, 1, 1, "qpkRelPmt", 0);
  plotRelQpk(relQpkPmt1Grp[1], relQpkRmsPmt1Grp[1],
      labelsPmt1Grp[1], 50, 1, 2, "qpkRelPmt", 0);
  plotRelQpk(relQpkPmt1Grp[2], relQpkRmsPmt1Grp[2],
      labelsPmt1Grp[2], stIds[0].size()-100, 1, 3, "qpkRelPmt", 0);

  plotRelQpk(relCQpkPmt1Grp[0], relCQpkRmsPmt1Grp[0], 
      labelsPmt1Grp[0], 50, 1, 1, "cqpkRelPmt", 1);
  plotRelQpk(relCQpkPmt1Grp[1], relCQpkRmsPmt1Grp[1],
      labelsPmt1Grp[1], 50, 1, 2, "cqpkRelPmt", 1);
  plotRelQpk(relCQpkPmt1Grp[2], relCQpkRmsPmt1Grp[2],
      labelsPmt1Grp[2], stIds[0].size()-100, 1, 3, "cqpkRelPmt", 1);

  // For PMT2
  vector < vector < double > > labelsPmt2Grp(3);   
  for ( int i=0; i<50; i++ ) {
    labelsPmt2Grp[0].push_back( stIdsLabel[1][i] );
    labelsPmt2Grp[1].push_back( stIdsLabel[1][i+50] );
  }
  for ( int i=0; i<stIds[1].size()-100; i++ )
    labelsPmt2Grp[2].push_back( stIdsLabel[1][i+100] );

  plotRelQpk(relQpkPmt2Grp[0], relQpkRmsPmt2Grp[0],
      labelsPmt2Grp[0], 50, 2, 1, "qpkRelPmt", 0);
  plotRelQpk(relQpkPmt2Grp[1], relQpkRmsPmt2Grp[1],
      labelsPmt2Grp[1], 50, 2, 2, "qpkRelPmt", 0);
  plotRelQpk(relQpkPmt2Grp[2], relQpkRmsPmt2Grp[2],
      labelsPmt2Grp[2], stIds[1].size()-100, 2, 3, "qpkRelPmt", 0);

  plotRelQpk(relCQpkPmt2Grp[0], relCQpkRmsPmt2Grp[0], 
      labelsPmt2Grp[0], 50, 2, 1, "cqpkRelPmt", 1);
  plotRelQpk(relCQpkPmt1Grp[1], relCQpkRmsPmt2Grp[1],
      labelsPmt2Grp[1], 50, 2, 2, "cqpkRelPmt", 1);
  plotRelQpk(relCQpkPmt2Grp[2], relCQpkRmsPmt2Grp[2],
      labelsPmt2Grp[2], stIds[1].size()-100, 2, 3, "cqpkRelPmt", 1);

  // For PMT3
  vector < vector < double > > labelsPmt3Grp(3);   
  for ( int i=0; i<50; i++ ) {
    labelsPmt3Grp[0].push_back( stIdsLabel[2][i] );
    labelsPmt3Grp[1].push_back( stIdsLabel[2][i+50] );
  }
  for ( int i=0; i<stIds[2].size()-100; i++ )
    labelsPmt3Grp[2].push_back( stIdsLabel[2][i+100] );

  plotRelQpk(relQpkPmt3Grp[0], relQpkRmsPmt3Grp[0],
      labelsPmt3Grp[0], 50, 3, 1, "qpkRelPmt", 0);
  plotRelQpk(relQpkPmt3Grp[1], relQpkRmsPmt3Grp[1],
      labelsPmt3Grp[1], 50, 3, 2, "qpkRelPmt", 0);
  plotRelQpk(relQpkPmt3Grp[2], relQpkRmsPmt3Grp[2],
      labelsPmt3Grp[2], stIds[2].size()-100, 3, 3, "qpkRelPmt", 0);

  plotRelQpk(relCQpkPmt3Grp[0], relCQpkRmsPmt3Grp[0], 
      labelsPmt3Grp[0], 50, 3, 1, "cqpkRelPmt", 1);
  plotRelQpk(relCQpkPmt3Grp[1], relCQpkRmsPmt3Grp[1],
      labelsPmt3Grp[1], 50, 3, 2, "cqpkRelPmt", 1);
  plotRelQpk(relCQpkPmt3Grp[2], relCQpkRmsPmt3Grp[2],
      labelsPmt3Grp[2], stIds[2].size()-100, 3, 3, "cqpkRelPmt", 1);

  // Plotting for relQpk distributions
  TH1D distRelQpkPmt1("distRelQpkPmt1", "", 20, -1, 2);
  TH1D distRelQpkPmt2("distRelQpkPmt2", "", 20, -1, 2);
  TH1D distRelQpkPmt3("distRelQpkPmt3", "", 20, -1, 2);

  TH1D distRelCQpkPmt1("distRelCQpkPmt1", "", 50, -5, 5);
  TH1D distRelCQpkPmt2("distRelCQpkPmt2", "", 50, -5, 5);
  TH1D distRelCQpkPmt3("distRelCQpkPmt3", "", 50, -5, 5);


  for ( int grp_i=0; grp_i<3; grp_i++ ) {
    for ( auto &i : relQpkPmt1Grp[grp_i] )
      distRelQpkPmt1.Fill( i );
    for ( auto &i : relQpkPmt2Grp[grp_i] )
      distRelQpkPmt2.Fill( i );
    for ( auto &i : relQpkPmt3Grp[grp_i] )
      distRelQpkPmt3.Fill( i );

    for ( auto &i : relCQpkPmt1Grp[grp_i] )
      distRelCQpkPmt1.Fill( i );
    for ( auto &i : relCQpkPmt2Grp[grp_i] )
      distRelCQpkPmt2.Fill( i );
    for ( auto &i : relCQpkPmt3Grp[grp_i] )
      distRelCQpkPmt3.Fill( i );
  }

  TCanvas cvnsDistQpkRelPmt1("cvnsDistQpkRelPmt1","", 1.6e3, 9e2);
  setCanvasStyle(cvnsDistQpkRelPmt1);

  distRelQpkPmt1.SetStats(kFALSE);
  distRelQpkPmt1.GetYaxis()->SetTitle("Counts [au]");
  distRelQpkPmt1.GetXaxis()->SetTitle("Q^{Kat}_{pk}/Q^{Offline}_{pk} - 1 [%]");
  distRelQpkPmt1.GetXaxis()->SetTitleOffset(1.4);
  distRelQpkPmt1.Fit("gaus","Q");
  distRelQpkPmt1.Draw();

  TLegend *lgnd = new TLegend(0.7, 0.7, 0.95, 0.95);
  lgnd->AddEntry(&distRelQpkPmt1, "PMT 1", "h");
  lgnd->AddEntry(&distRelQpkPmt1, Form( "Mean: %.2f [%] #pm %.2f [%]",
        distRelQpkPmt1.GetMean(), distRelQpkPmt1.GetMeanError()), "h");
  lgnd->AddEntry(distRelQpkPmt1.GetFunction("gaus"), "Gaussian fit", "l");
  lgnd->AddEntry(distRelQpkPmt1.GetFunction("gaus"), Form("#mu: %.2f #pm %.2f",
        distRelQpkPmt1.GetFunction("gaus")->GetParameter(1),
        distRelQpkPmt1.GetFunction("gaus")->GetParError(1)), "");
  lgnd->AddEntry(distRelQpkPmt1.GetFunction("gaus"), Form("#sigma: %.2f #pm %.2f",
        distRelQpkPmt1.GetFunction("gaus")->GetParameter(2),
        distRelQpkPmt1.GetFunction("gaus")->GetParError(2)), "");
  lgnd->SetBorderSize(0);
  lgnd->SetTextSize(0.045);
  lgnd->Draw();

  cvnsDistQpkRelPmt1.Print("results/distQpkRelPmt1.pdf");

  TCanvas cvnsDistQpkRelPmt2("cvnsDistQpkRelPmt2","", 1.6e3, 9e2);
  setCanvasStyle(cvnsDistQpkRelPmt2);

  distRelQpkPmt2.SetStats(kFALSE);
  distRelQpkPmt2.GetYaxis()->SetTitle("Counts [au]");
  distRelQpkPmt2.GetXaxis()->SetTitle("Q^{Kat}_{pk}/Q^{Offline}_{pk} - 1 [%]");
  distRelQpkPmt2.GetXaxis()->SetTitleOffset(1.4);
  distRelQpkPmt2.Fit("gaus","Q");
  distRelQpkPmt2.Draw();

  lgnd = new TLegend(0.7, 0.7, 0.95, 0.95);
  lgnd->AddEntry(&distRelQpkPmt2, "PMT 2", "h");
  lgnd->AddEntry(&distRelQpkPmt2, Form( "Mean: %.2f [%] #pm %.2f [%]",
        distRelQpkPmt2.GetMean(), distRelQpkPmt2.GetMeanError()), "h");
  lgnd->AddEntry(distRelQpkPmt2.GetFunction("gaus"), "Gaussian fit", "l");
  lgnd->AddEntry(distRelQpkPmt2.GetFunction("gaus"), Form("#mu: %.2f #pm %.2f",
        distRelQpkPmt2.GetFunction("gaus")->GetParameter(1),
        distRelQpkPmt2.GetFunction("gaus")->GetParError(1)), "");
  lgnd->AddEntry(distRelQpkPmt2.GetFunction("gaus"), Form("#sigma: %.2f #pm %.2f",
        distRelQpkPmt2.GetFunction("gaus")->GetParameter(2),
        distRelQpkPmt2.GetFunction("gaus")->GetParError(2)), "");
  lgnd->SetBorderSize(0);
  lgnd->SetTextSize(0.045);
  lgnd->Draw();

  cvnsDistQpkRelPmt2.Print("results/distQpkRelPmt2.pdf");

  TCanvas cvnsDistQpkRelPmt3("cvnsDistQpkRelPmt3","", 1.6e3, 9e2);
  setCanvasStyle(cvnsDistQpkRelPmt3);

  distRelQpkPmt3.SetStats(kFALSE);
  distRelQpkPmt3.GetYaxis()->SetTitle("Counts [au]");
  distRelQpkPmt3.GetXaxis()->SetTitle("Q^{Kat}_{pk}/Q^{Offline}_{pk} - 1 [%]");
  distRelQpkPmt3.GetXaxis()->SetTitleOffset(1.4);
  distRelQpkPmt3.Fit("gaus","Q");
  distRelQpkPmt3.Draw();

  lgnd = new TLegend(0.7, 0.7, 0.95, 0.95);
  lgnd->AddEntry(&distRelQpkPmt3, "PMT 3", "h");
  lgnd->AddEntry(&distRelQpkPmt3, Form( "Mean: %.2f [%] #pm %.2f [%]",
        distRelQpkPmt3.GetMean(), distRelQpkPmt3.GetMeanError()), "h");
  lgnd->AddEntry(distRelQpkPmt3.GetFunction("gaus"), "Gaussian fit", "l");
  lgnd->AddEntry(distRelQpkPmt3.GetFunction("gaus"), Form("#mu: %.2f #pm %.2f",
        distRelQpkPmt3.GetFunction("gaus")->GetParameter(1),
        distRelQpkPmt3.GetFunction("gaus")->GetParError(1)), "");
  lgnd->AddEntry(distRelQpkPmt1.GetFunction("gaus"), Form("#sigma: %.2f #pm %.2f",
        distRelQpkPmt3.GetFunction("gaus")->GetParameter(2),
        distRelQpkPmt3.GetFunction("gaus")->GetParError(2)), "");
  lgnd->SetBorderSize(0);
  lgnd->SetTextSize(0.045);
  lgnd->Draw();

  cvnsDistQpkRelPmt3.Print("results/distQpkRelPmt3.pdf");


  TCanvas cvnsDistCQpkRelPmt1("cvnsDistCQpkRelPmt1","", 1.6e3, 9e2);
  setCanvasStyle(cvnsDistCQpkRelPmt1);

  distRelCQpkPmt1.SetStats(kFALSE);
  distRelCQpkPmt1.GetYaxis()->SetTitle("Counts [au]");
  distRelCQpkPmt1.GetXaxis()->SetTitle("CQ^{Kat}_{pk}/CQ^{Offline}_{pk} - 1 [%]");
  distRelCQpkPmt1.GetXaxis()->SetTitleOffset(1.4);
  distRelCQpkPmt1.Fit("gaus","Q");
  distRelCQpkPmt1.Draw();

  lgnd = new TLegend(0.7, 0.7, 0.95, 0.95);
  lgnd->AddEntry(&distRelCQpkPmt1, "PMT 1", "h");
  lgnd->AddEntry(&distRelCQpkPmt1, Form( "Mean: %.2f [%] #pm %.2f [%]",
        distRelCQpkPmt1.GetMean(), distRelCQpkPmt1.GetMeanError()), "h");
  lgnd->AddEntry(distRelCQpkPmt1.GetFunction("gaus"), "Gaussian fit", "l");
  lgnd->AddEntry(distRelCQpkPmt1.GetFunction("gaus"), Form("#mu: %.2f #pm %.2f",
        distRelCQpkPmt1.GetFunction("gaus")->GetParameter(1),
        distRelCQpkPmt1.GetFunction("gaus")->GetParError(1)), "");
  lgnd->AddEntry(distRelCQpkPmt1.GetFunction("gaus"), Form("#sigma: %.2f #pm %.2f",
        distRelCQpkPmt1.GetFunction("gaus")->GetParameter(2),
        distRelCQpkPmt1.GetFunction("gaus")->GetParError(2)), "");
  lgnd->SetBorderSize(0);
  lgnd->SetTextSize(0.045);
  lgnd->Draw();

  cvnsDistCQpkRelPmt1.Print("results/distCQpkRelPmt1.pdf");

  TCanvas cvnsDistCQpkRelPmt2("cvnsDistCQpkRelPmt2","", 1.6e3, 9e2);
  setCanvasStyle(cvnsDistCQpkRelPmt2);

  distRelCQpkPmt2.SetStats(kFALSE);
  distRelCQpkPmt2.GetYaxis()->SetTitle("Counts [au]");
  distRelCQpkPmt2.GetXaxis()->SetTitle("CQ^{Kat}_{pk}/CQ^{Offline}_{pk} - 1 [%]");
  distRelCQpkPmt2.GetXaxis()->SetTitleOffset(1.4);
  distRelCQpkPmt2.Fit("gaus","Q");
  distRelCQpkPmt2.Draw();

  lgnd = new TLegend(0.7, 0.7, 0.95, 0.95);
  lgnd->AddEntry(&distRelCQpkPmt2, "PMT 2", "h");
  lgnd->AddEntry(&distRelCQpkPmt2, Form( "Mean: %.2f [%] #pm %.2f [%]",
        distRelCQpkPmt2.GetMean(), distRelCQpkPmt2.GetMeanError()), "h");
  lgnd->AddEntry(distRelCQpkPmt2.GetFunction("gaus"), "Gaussian fit", "l");
  lgnd->AddEntry(distRelCQpkPmt2.GetFunction("gaus"), Form("#mu: %.2f #pm %.2f",
        distRelCQpkPmt2.GetFunction("gaus")->GetParameter(1),
        distRelCQpkPmt2.GetFunction("gaus")->GetParError(1)), "");
  lgnd->AddEntry(distRelCQpkPmt2.GetFunction("gaus"), Form("#sigma: %.2f #pm %.2f",
        distRelCQpkPmt2.GetFunction("gaus")->GetParameter(2),
        distRelCQpkPmt2.GetFunction("gaus")->GetParError(2)), "");
  lgnd->SetBorderSize(0);
  lgnd->SetTextSize(0.045);
  lgnd->Draw();

  cvnsDistCQpkRelPmt2.Print("results/distCQpkRelPmt2.pdf");

  TCanvas cvnsDistCQpkRelPmt3("cvnsDistCQpkRelPmt3","", 1.6e3, 9e2);
  setCanvasStyle(cvnsDistCQpkRelPmt3);

  distRelCQpkPmt3.SetStats(kFALSE);
  distRelCQpkPmt3.GetYaxis()->SetTitle("Counts [au]");
  distRelCQpkPmt3.GetXaxis()->SetTitle("CQ^{Kat}_{pk}/CQ^{Offline}_{pk} - 1 [%]");
  distRelCQpkPmt3.GetXaxis()->SetTitleOffset(1.4);
  distRelCQpkPmt3.Fit("gaus","Q");
  distRelCQpkPmt3.Draw();

  lgnd = new TLegend(0.7, 0.7, 0.95, 0.95);
  lgnd->AddEntry(&distRelCQpkPmt3, "PMT 3", "h");
  lgnd->AddEntry(&distRelCQpkPmt3, Form( "Mean: %.2f [%] #pm %.2f [%]",
        distRelCQpkPmt3.GetMean(), distRelCQpkPmt3.GetMeanError()), "h");
  lgnd->AddEntry(distRelCQpkPmt3.GetFunction("gaus"), "Gaussian fit", "l");
  lgnd->AddEntry(distRelCQpkPmt3.GetFunction("gaus"), Form("#mu: %.2f #pm %.2f",
        distRelCQpkPmt3.GetFunction("gaus")->GetParameter(1),
        distRelCQpkPmt3.GetFunction("gaus")->GetParError(1)), "");
  lgnd->AddEntry(distRelCQpkPmt3.GetFunction("gaus"), Form("#sigma: %.2f #pm %.2f",
        distRelCQpkPmt3.GetFunction("gaus")->GetParameter(2),
        distRelCQpkPmt3.GetFunction("gaus")->GetParError(2)), "");
  lgnd->SetBorderSize(0);
  lgnd->SetTextSize(0.045);
  lgnd->Draw();

  cvnsDistCQpkRelPmt3.Print("results/distCQpkRelPmt3.pdf");

  plotCoinciFactor(coinciFactorKatPmt1, coinciFactorMaoPmt1, "1");
  plotCoinciFactor(coinciFactorKatPmt2, coinciFactorMaoPmt2, "2");
  plotCoinciFactor(coinciFactorKatPmt3, coinciFactorMaoPmt3, "3");

  exit(0);
} 

double getMean(vector<double> vect) {
  double mean = 0.;
  for ( auto &i : vect )
    mean += i;
  return mean/vect.size();
}

double getRMS(vector<double> vect, double mean) {
  double rms = 0.;
  for ( auto &i : vect )
    rms += (i-mean)*(i-mean); 

  return sqrt( rms/( vect.size() ));
}

void setCanvasStyle(TCanvas& canvas) {
  canvas.SetTopMargin(0.03);
  canvas.SetLeftMargin(0.08);
  canvas.SetRightMargin(0.02);
}

void plotRelQpk(vector<double> grp, vector<double> grpRms, 
    vector<double> labels, int sizeGrp, int pmt, int grpN, TString name,
    int isCQpk) {

  TCanvas cvnsQpkRel("cvnsQpkRel","", 1.6e3, 9e2);
  setCanvasStyle(cvnsQpkRel);

  TH1D *dist = new TH1D(Form(name+"dist%d%d", pmt, grpN), 
      "", sizeGrp, 0, sizeGrp);
  for ( int i=0; i<sizeGrp; i++ ) {
    dist->Fill( i, grp[i] );
    dist->SetBinError( i, grpRms[i] );
  }
  
  dist->SetStats(kFALSE);
  dist->SetMarkerColor(kBlue);
  dist->SetLineColor(kBlue);
  dist->SetMarkerStyle(21);
  dist->GetXaxis()->SetRangeUser(0, 50);
  dist->GetXaxis()->SetTitle("St. ID.");
  dist->GetXaxis()->SetTitleOffset(1.4);
  if ( isCQpk )
    dist->GetYaxis()->SetTitle("CQ^{Kat}_{pk}/CQ^{Offline}_{pk} - 1 [%]");
  else 
    dist->GetYaxis()->SetTitle("Q^{Kat}_{pk}/Q^{Offline}_{pk} - 1 [%]");
  dist->GetYaxis()->SetTitleOffset(1.);
  for (int i=0; i<labels.size(); i++ )
    dist->GetXaxis()->SetBinLabel(i+1, Form( "%.f", labels[i] ) );
  dist->LabelsOption("v");
  dist->Draw();

  TLegend *lgnd = new TLegend(0.7, 0.85, 0.95, 0.95);
  lgnd->AddEntry(dist, Form( "PMT %d", pmt), "h");
  lgnd->AddEntry(dist, Form( "Mean: %.2f [%] #pm %.2f [%]", 
        getMean(grp), getRMS(grp, getMean(grp))/sqrt(grp.size()) ), "h");

  lgnd->SetBorderSize(0);
  lgnd->SetTextSize(0.045);
  lgnd->Draw();
   
  cvnsQpkRel.Print(Form("results/"+name+"%dGrp%d.pdf", pmt, grpN));
}


void plotCoinciFactor(TH1D &kat, TH1D &mao, TString pmt) {
  TCanvas cvnsFactorPmt("cvnsFactorPmt"+pmt,"", 1.6e3, 9e2);
  setCanvasStyle(cvnsFactorPmt);
  cvnsFactorPmt.SetLogy();

  mao.Fit("gaus", "Q");
  mao.GetFunction("gaus")->SetLineStyle(9);
  mao.GetFunction("gaus")->SetLineColor(kRed);
  kat.Fit("gaus", "Q");
  kat.GetFunction("gaus")->SetLineStyle(9);
  kat.GetFunction("gaus")->SetLineColor(kBlue);

  kat.SetStats(kFALSE);
  kat.SetLineColor(kBlue);
  kat.SetLineWidth(2);
  kat.GetXaxis()->SetRangeUser(0.7, 1.5);
  kat.GetXaxis()->SetTitle("Q^{pk}_{CH} / Q^{pk}");
  kat.GetYaxis()->SetTitle("Counts [au]");
  kat.SetFillColorAlpha(kBlue, 0.15);
  kat.Draw();

  mao.SetLineColor(kRed);
  mao.SetLineWidth(2);
  mao.SetFillColorAlpha(kRed, 0.15);
  mao.Draw("same");

  TLegend *lgnd = new TLegend(0.7, 0.3, 0.9, 0.96);
  lgnd->AddEntry(&kat, "PMT "+pmt, "h");
  lgnd->AddEntry(&kat, "Katarina's", "l");
  lgnd->AddEntry(&kat, Form( "Mean: %.2f [%] #pm %.2f [%]",
        kat.GetMean(), kat.GetMeanError()), "h");
  lgnd->AddEntry(&kat, Form( "RMS: %.2f [%] #pm %.2f [%]",
        kat.GetRMS(), kat.GetRMSError()), "h");
  lgnd->AddEntry(kat.GetFunction("gaus"), "Gaussian fit", "h");
  lgnd->AddEntry(kat.GetFunction("gaus"), Form(
        "#mu: %.2f [%] #pm %.2f", 
        kat.GetFunction("gaus")->GetParameter(1),
        kat.GetFunction("gaus")->GetParError(1)),"");
  lgnd->AddEntry(kat.GetFunction("gaus"), Form( 
        "#sigma: %.2f [%] #pm %.2f", 
        kat.GetFunction("gaus")->GetParameter(2),
        kat.GetFunction("gaus")->GetParError(2)),"");
  lgnd->AddEntry(&kat, "", "h");

  lgnd->AddEntry(&mao, "OffLine's", "l");
  lgnd->AddEntry(&mao, Form( "Mean: %.2f [%] #pm %.2f [%]",
        mao.GetMean(), mao.GetMeanError()), "h");
  lgnd->AddEntry(&mao, Form( "RMS: %.2f [%] #pm %.2f [%]",
        mao.GetRMS(), mao.GetRMSError()), "h");

  lgnd->AddEntry(mao.GetFunction("gaus"), "Gaussian fit", "h");
  lgnd->AddEntry(mao.GetFunction("gaus"), Form( 
        "#mu: %.2f [%] #pm %.2f", 
        mao.GetFunction("gaus")->GetParameter(1),
        mao.GetFunction("gaus")->GetParError(1)),"");
  lgnd->AddEntry(mao.GetFunction("gaus"), Form( 
        "#sigma: %.2f [%] #pm %.2f", 
        mao.GetFunction("gaus")->GetParameter(2),
        mao.GetFunction("gaus")->GetParError(2)),"");
  lgnd->AddEntry(&mao, "", "h");

  lgnd->SetBorderSize(0);
  lgnd->SetTextSize(0.04);
  lgnd->Draw();


  cvnsFactorPmt.Print("results/coinciFactorPmt"+pmt+".pdf");
}
