#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include <IoAuger.h>
#include <Ec.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TGraphErrors.h>

using namespace std;

// ========================== 
// ******** The MAIN ********
// ==========================
int main (int argc, char *argv[]) {
   if ( argc < 3 ) {
		 cout << endl
         << "Usage: " << argv[0] << " <stationsFile>  <files>" << endl
         << "  <stationsFile>: file with a list of stations" << endl
         << "  <files>: IoSd or IoAuger files to be read" << endl
				 << " " << endl
				 << "In case you want the distribution of all events for a specific Station, " << endl
				 << "just make sure the stationsFile conteins a single station." << endl
				 << endl;
 		 exit(0);
	 }
	
  const char* stationsFileName = argv[1];
  AugerIoSd input(argc-2, argv+2);
  
  ifstream stationsFile(stationsFileName, ios::in);
  if (!stationsFile.is_open()){
    cout << "Could not open file: " << stationsFileName << endl;
    exit(0);
  }
  vector<unsigned int> stationsIds;
  vector< bool > stationsIdsFlag;
  while (stationsFile.good()) {
    unsigned int st = 0;
    stationsFile >> st;
    if (st) {
      stationsIds.push_back(st);
			stationsIdsFlag.push_back(true);
		}
  }
	unsigned int allStread = 0;
 
  if (stationsIds.empty()){
    cout << "Please specify the stations ids in the file " << endl;
    exit(0);
  }  

  EventPos pos;

  for (pos=input.FirstEvent(); pos<input.LastEvent(); pos=input.NextEvent()) {
		bool found = false;
		IoSdEvent event(pos);

    for (unsigned int i = 0 ; i < event.Stations.size(); ++i) {
      found = false;
      for ( unsigned int cSt = 0; cSt<stationsIds.size(); ++cSt )
        if (event.Stations[i].Id == stationsIds[cSt] ) {
					IoSdEvent event(pos);
					if ( stationsIds[cSt] == 1223 && event.Stations[i].IsUUB )
						cerr << event.Id << endl;
					if ( event.Stations[i].IsUUB && stationsIdsFlag[cSt] ) {
						cout << stationsIds[cSt] << " " 
							<< event.Stations[i].easting() << " "
							<< event.Stations[i].northing() << " "
							<< endl;
						stationsIdsFlag[cSt] = false;
						allStread++;
					}
          found = true;
				}
      if ( !found )
        continue;
		}
		if ( allStread == stationsIdsFlag.size() ) {
			cerr << allStread << " Stations have been read" << endl;
			break;
		}
  }
	return 0;
}
