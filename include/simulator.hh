///////////////////////////////////////////////////
//   High field simulation for MuSEUM Experiment 
//
//      Author:  Hideharu Yamauchi 2021/09/18 
///////////////////////////////////////////////////

#if !defined(___header_simulator_)
#define ___header_simulator_ 1

#include "TMath.h"
#include <vector>
#include "TString.h"
#include "TTree.h"
#include "TF1.h"
#include "HFgeometry.hh"
#include "TCanvas.h"

class simulator{
private:
  std::vector<double> position;
  double angle;
  double solid_angle;
  double signal;
  int Mode;
  double state_amp[2];
  int scan_points = 40;
  const double gamma = 1/muon_life;
  double b; // the coefficient to change RF field to frequency
  Double_t decaytime, decaypositionx, decaypositiony, decaypositionz, RF, Effective_RF, magnet_field;
  int entries;
  
public:
  double Non;
  double Noff;
  double L; // microwave term describe increase of positron counting due to microwave irradiation
  simulator(TTree* decaytree, int mode);
  ~simulator(void){std::cout << "finish simulator..." << std::endl;};
  TF1* f1;
  void GetXYZ(int x, int y, int z);
  void timedev(double t, double delta, double gamma, double position[3]);
  void Vis_state_amp(void);
  double ConventionalSignal(double position[3]);
  double OldMuoniumSignal(double position[3], Double_t windowopen, Double_t timewindow);
};
#endif
