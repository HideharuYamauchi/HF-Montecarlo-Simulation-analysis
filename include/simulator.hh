////////////////////////////////////////////////////////
//    High field simulation for MuSEUM Collaboration
//
//        Author:  Hideharu Yamauchi 2021/09/18 
////////////////////////////////////////////////////////
#if !defined(___header_simulator_)
#define ___header_simulator_ 1

#include "TString.h"
#include "TTree.h"
#include "TF1.h"
#include "HFgeometry.hh"
#include "TCanvas.h"
#include <vector>

class simulator{
private:
  std::vector<Double_t> position;
  double angle;
  double solid_angle;
  double signal;
  int Mode;
  double state_amp[2];
  int scan_range = 200; // kHz
  int scan_points = 40; // scan points, same with liu
  const double gamma = 1/muon_life; // muon's natural width
  int entries;
  double Calculate_g(Double_t Gamma, Double_t t);
  
public:
  double Non;
  double Noff;
  std::vector<Int_t> a;
  double L; // microwave term describe increase of positron counting due to microwave irradiation
  double K; // solid_angle integral and cavity volume integral
  simulator(TTree* decaytree, int mode);
  ~simulator(void){std::cout << "finish simulator..." << std::endl;};
  TF1* f1;
  //void GetXYZ(int x, int y, int z);
  void timedev(double t, double b, double delta, double gamma, double position[3]);
  //void Vis_state_amp(void);
  //double ConventionalSignal(double power, double detuning, double position[3]);
  //double OldMuoniumSignal(double power, double detuning, double position[3], Double_t windowopen, Double_t windowclose);
  //double Calculate_Signal(void);
};
#endif
