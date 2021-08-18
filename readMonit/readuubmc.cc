#include <TSDMonCal.h>
#include <TChain.h>
#include <iostream>
#include <cstdlib>

using namespace std;


int
main(int argc, char* argv[])
{
  TChain chain("SDMonCal");
  chain.SetBranchStatus("fRawMonitoring.fListOfMembers", 0);
  chain.SetBranchStatus("fCalibration.fListOfMembers", 0);
  for (int i = 1; i < argc; ++i) 
	{
    cerr << "Opening " << argv[i] << endl;
    chain.Add(argv[i]);
  }
  TSDMonCal sdMonCal;
  TSDMonCal* mc = &sdMonCal;
  chain.SetBranchAddress("SDMonCalBranch", &mc);
  const int numberOfEntries = chain.GetEntries();
  for (int i = 0; i < numberOfEntries; ++i) 
	{
    if (!(i % 10000))
      cerr << i << '/' << numberOfEntries << endl;
    chain.GetEvent(i);
    //if (mc->fRawMonitoring.fIsUUB != 1)
      //continue;
    const TSDMonitoring& m = mc->fMonitoring;
		if ( mc->fLsId==1740 ) 
		{
			//cout << mc->fTime << " " << m.fPMV[0] << " " << m.fPMV[1] << " " << m.fPMV[2] << endl;
			cout << mc->fTime << " " << m.fUPMV[0] << " " << m.fUPMV[1] << " " << m.fUPMV[2] << endl;
			//cout << "PMV: " << m.fPMV[0] << " " << m.fPMV[1] << " " << m.fPMV[2] << endl;
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
