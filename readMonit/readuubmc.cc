#include <TSDMonCal.h>
#include <TChain.h>
#include <iostream>
#include <cstdlib>

using namespace std;


int
main(int argc, char* argv[]) {
  TChain chain("SDMonCal");
  chain.SetBranchStatus("fRawMonitoring.fListOfMembers", 0);
  chain.SetBranchStatus("fCalibration.fListOfMembers", 0);
  for (int i = 1; i < argc; ++i) {
    cerr << "Opening " << argv[i] << endl;
    chain.Add(argv[i]);
  }
  TSDMonCal sdMonCal;
  TSDMonCal* mc = &sdMonCal;
  chain.SetBranchAddress("SDMonCalBranch", &mc);
  const int numberOfEntries = chain.GetEntries();
  for (int i = 0; i < numberOfEntries; ++i) {
    if (!(i % int(1e5)))
      cerr << i << '/' << numberOfEntries << endl;
    chain.GetEvent(i);
    if ( !mc->fRawMonitoring.fIsUUB )
      continue;
    //const TSDMonitoring& m = mc->fMonitoring;
    //const TSDMonitoring& m = mc->fRawMonitoring;
		if ( mc->fLsId==545 ) {
      cout << mc->fTime << " ";
      for( int pmt=0; pmt<3; pmt++ )
        cout << mc->fRawMonitoring.fPMV[pmt] << " " << mc->fRawMonitoring.fPMT[pmt] << " " << mc->fCalibration.fArea[pmt] << " ";
      cout << endl;
			//cout << "fPMV " << mc->fTime << " " << mc->fRawMonitoring.fPMV[0] << " " << mc->fRawMonitoring.fPMV[1] << " " << mc->fRawMonitoring.fPMV[2] << endl;
			//cout << "fUPMV " << mc->fTime << " " << mc->fRawMonitoring.fUPMV[0] << " " << mc->fRawMonitoring.fUPMV[1] << " " << mc->fRawMonitoring.fUPMV[2] << endl;
			/*for (int i=0;i<8;i++) 
				 cout << mc->fRawMonitoring.fExtra[i] << ' ';
			 cout << endl;
			 */
		}

		/*
    cout << mc->fTime << " " 
			<< m.fUPMV[0] << " " 
			<< mc->fRawMonitoring.fDACPM[0]
			<< endl;
			*/
  }
  return EXIT_SUCCESS;
}
