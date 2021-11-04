#include <iostream>
#include <fstream>

#include <TSDMonCal.h>
#include <TChain.h>
#include <iostream>
#include <cstdlib>

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>

using namespace std;


int main(int argc, char* argv[]) {
 
  if ( argc < 4 ) {
    cout << endl
      << "Usage: " << argv[0] << " <stationsFile>  <monitFiles> <isUUbSt>" << endl
      << " <stationsFile>: file with a list of stations" << endl
      << "  <monitFiles>: monit files to be read" << endl
      << " <isUUbSt> (0 or 1) Bool value to identify if you are reading UUB stations." 
      << endl
      << " " << endl;
    exit(0);
  }
  // Reading station list
  const char* stationsFileName = argv[1];
  ifstream stationsFile(stationsFileName, ios::in);

  if (!stationsFile.is_open()) {
    cout << "Could not open file: " << stationsFileName << endl;
    exit(0);
  }
  vector<unsigned int> stationsIds;
  unsigned int st = 0;
  while (stationsFile.good()) {
    st = 0;
    stationsFile >> st;
    if (st)
      stationsIds.push_back(st);
  }
  if (stationsIds.empty()) {
    cout << "Please specify the stations ids in the file " << endl;
    exit(0);
  }
  
  bool ifItUub = (char)atoi(argv[2]);

  TChain chain("SDMonCal");
  chain.SetBranchStatus("fRawMonitoring.fListOfMembers", 0);
  chain.SetBranchStatus("fCalibration.fListOfMembers", 0);

  for (int i = 3; i < argc; ++i) {
    cerr << "Opening " << argv[i] << endl;
    chain.Add(argv[i]);
  }

  TSDMonCal sdMonCal;
  TSDMonCal* mc = &sdMonCal;
  chain.SetBranchAddress("SDMonCalBranch", &mc);
  const int numberOfEntries = chain.GetEntries();
  bool found = false;
  unsigned int stId = 0;

  TString strStId;
  strStId.Form("stId%d", stId);
  TString strIfUub = (ifItUub) ? "Uub" : "Ub";
  strStId.Form("monitVolt"+strIfUub+"St%d", stationsIds[0]);
  TFile hfile(strStId+".root", "RECREATE","");
  TTree *stationsTree = new TTree(strStId,"");
  
  int time = 0;
  double pmv1 = 0.;
  double pmv2 = 0.;
  double pmv3 = 0.;
  double upmv1 = 0.;
  double upmv2 = 0.;
  double upmv3 = 0.;
  int area1 = 0.;
  int area2 = 0.;
  int area3 = 0.;

  stationsTree->Branch("fTime",&time,"time/I");
  stationsTree->Branch("fPMV1",&pmv1,"pmv1/D");
  stationsTree->Branch("fPMV2",&pmv2,"pmv2/D");
  stationsTree->Branch("fPMV3",&pmv3,"pmv3/D");
  stationsTree->Branch("fUPMV1",&upmv1,"upmv1/D");
  stationsTree->Branch("fUPMV2",&upmv2,"upmv2/D");
  stationsTree->Branch("fUPMV3",&upmv3,"upmv3/D");
  stationsTree->Branch("fArea1",&area1,"area1/I");
  stationsTree->Branch("fArea2",&area2,"area2/I");
  stationsTree->Branch("fArea3",&area3,"area3/I");

  for (int i = 0; i < numberOfEntries; ++i) {
    if (!(i % int(1e5)))
      cerr << i << '/' << numberOfEntries << endl;
    chain.GetEvent(i);
  
    if ( ifItUub && !mc->fRawMonitoring.fIsUUB ) 
      continue;

    for (  vector<unsigned int>::const_iterator iter= stationsIds.begin();
        iter!= stationsIds.end(); ++iter) {
      stId = mc->fLsId;
      if (stId == *iter )
        found = true;
    }
    if ( !found )
      continue;
    
    time = mc->fTime;
    pmv1 = mc->fMonitoring.fPMV[0];
    pmv2 = mc->fMonitoring.fPMV[1];
    pmv3 = mc->fMonitoring.fPMV[2];
    upmv1 = mc->fMonitoring.fUPMV[0];
    upmv2 = mc->fMonitoring.fUPMV[1];
    upmv3 = mc->fMonitoring.fUPMV[2];
    area1 = mc->fCalibration.fArea[0];
    area2 = mc->fCalibration.fArea[1];
    area3 = mc->fCalibration.fArea[2];
    stationsTree->Fill(); 
  }
  hfile.Write();
  hfile.Close();
  return EXIT_SUCCESS;
}
