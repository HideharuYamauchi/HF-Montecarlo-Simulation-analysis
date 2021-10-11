////////////////////////////////////////////////////////
//    High field simulation for MuSEUM Collaboration
//
//        Author:  Hideharu Yamauchi 2021/09/18 
////////////////////////////////////////////////////////
#if !defined(___header_simulator_)
#define ___header_simulator_ 1

#include <vector>
#include "TString.h"
#include "TTree.h"
#include "TF1.h"
#include "HFgeometry.hh"
#include "magnet.hh"
#include "TCanvas.h"

class simulator{
private:
  std::vector<double> position;
  Double_t magnet_field, RF, b; // std::vector<Double_t> field;
  Double_t decaytime, positron_energy;
  std::string* decayvolume;
  TBranch* decayvolume_branch;
  std::vector<Double_t>* muon_position;
  TBranch* muon_position_branch;
  std::vector<Double_t>* muon_momentum;
  TBranch* muon_momentum_branch;
  std::vector<Double_t>* positron_position;
  TBranch* positron_position_branch;
  std::vector<Double_t>* positron_momentum;
  TBranch* positron_momentum_branch;
  
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
  double L; // microwave term describe increase of positron counting due to microwave irradiation
  double K; // solid_angle integral and cavity volume integral
  simulator(TTree* decaytree, int mode, magfield* magnet);
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
