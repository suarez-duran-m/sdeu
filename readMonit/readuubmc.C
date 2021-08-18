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
  for (int i = 1; i < argc; ++i) {
    cerr << "Opening " << argv[i] << endl;
    chain.Add(argv[i]);
  }
  TSDMonCal sdMonCal;
  TSDMonCal* mc = &sdMonCal;
  chain.SetBranchAddress("SDMonCalBranch", &mc);
  const int numberOfEntries = chain.GetEntries();
  for (int i = 0; i < numberOfEntries; ++i) {
    if (!(i % 10000))
      cerr << i << '/' << numberOfEntries << endl;
    chain.GetEvent(i);
    if (mc->fRawMonitoring.fIsUUB != 1)
      continue;
    const TSDMonitoring& m = mc->fMonitoring;
    cout
      /*01*/ << mc->fLsId << ' '     // station id
      /*02*/ << mc->fTime << ' '     // unix second
      /*03*/ << mc->fCDASTime << ' ' // unix second
      /*04*/ << m.f3V3Analog << ' '
      /*05*/ << m.f_3V3Analog << ' '
      /*06*/ << m.f1V0 << ' '
      /*07*/ << m.f1V2 << ' '
      /*08*/ << m.f1V8 << ' '
      /*09*/ << m.fPMV[0] << ' '
      /*10*/ << m.fPMV[1] << ' '
      /*11*/ << m.fPMV[2] << ' '
      /*12*/ << m.fPMI[0] << ' '
      /*13*/ << m.fPMI[1] << ' '
      /*14*/ << m.fPMI[2] << ' '
      /*15*/ << m.fPMT[0] << ' '
      /*16*/ << m.fPMT[1] << ' '
      /*17*/ << m.fPMT[2] << ' '
      /*18*/ << m.fElectT << ' '
      /*19*/ << m.fBatteryT[0] << ' '
      /*20*/ << m.fBatteryT[1] << ' '
      /*21*/ << m.fBatteryV[0] << ' '
      /*22*/ << m.fBatteryV[1] << ' '
      /*23*/ << m.fSolarPanelV << ' '
      /*24*/ << m.fSolarPanelI << ' '
      /*25*/ << m.fWaterLevel << ' '
      /*26*/ << m.fWaterT << ' '
      /*27*/ << m.fCurrentLoad << ' '
      /*28*/ << m.f12V_WT << ' '
      /*29*/ << m.f12VPM << ' '
      /*30*/ << m.f3V3 << ' '
      /*31*/ << m.f12VRadioI << ' '
      /*32*/ << m.f12VRadio << ' '
      /*33*/ << m.f5VGPS << ' '
      /*34*/ << m.f5VUSB << ' '
      /*35*/ << m.f12VPMT[0] << ' '
      /*36*/ << m.f12VPMT[1] << ' '
      /*37*/ << m.f12VPMT[2] << ' '
      /*38*/ << m.f24VExt[0] << ' '
      /*39*/ << m.f24VExt[1] << ' '
      /*40*/ << m.fUPMV[0] << ' '
      /*41*/ << m.fUPMV[1] << ' '
      /*42*/ << m.fUPMV[2] << ' '
      /*43*/ << m.fUPMI[0] << ' '
      /*44*/ << m.fUPMI[1] << ' '
      /*45*/ << m.fUPMI[2] << ' '
      /*46*/ << m.fUPMT[0] << ' '
      /*47*/ << m.fUPMT[1] << ' '
      /*48*/ << m.fUPMT[2] << ' '
      /*49*/ << m.fExtT << ' '
      /*50*/ << m.fAirT << ' '
      /*51*/ << m.fAirP << ' '
      /*52*/ << m.f12VPMI << ' '
      /*53*/ << m.f3V3AnalogI << ' '
      /*54*/ << m.f_3V3AnalogI << ' '
      /*55*/ << m.f5VGPSI << ' '
      /*56*/ << m.f1V0I << ' '
      /*57*/ << m.f1V2I << ' '
      /*58*/ << m.f1V8I << ' '
      /*59*/ << m.f3V3I << ' '
      /*60*/ << m.fVInI << ' '
      /*61*/ << m.f3V3SCI << ' ';
      /*62+69*/ for (int i=0;i<8;i++) cout << mc->fRawMonitoring.fExtra[i] << ' ';

      cout << endl;
  }
  return EXIT_SUCCESS;
}
