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
#include "RF.hh" // for omega
#include "magnet.hh"
#include "TCanvas.h"
#include <vector>

class simulator{
private:
  TFile* treefile;
  TTree* tree;
  double angle;
  double solid_angle;
  double signal;
  int scan_range = 200; // kHz
  int scan_points = 40; // scan points, same with liu
  const double gamma = 1/muon_life; // muon's natural width
  int entries;
  //double Calculate_g(Double_t Gamma, Double_t t);
  
public:
  Double_t Calculate_EnergySplit();
  Double_t Non;
  Double_t Noff;
  Double_t Amplitude[2];
  double L; // microwave term describe increase of positron counting due to microwave irradiation
  double K; // solid_angle integral and cavity volume integral
  simulator(const char* rootfile);
  ~simulator(void);
  Double_t* timedev(Double_t delta);
  void Vis_State_Amp(Int_t entry);
  //double ConventionalSignal(double power, double detuning, double position[3]);
  //double OldMuoniumSignal(double power, double detuning, double position[3], Double_t windowopen, Double_t windowclose);
  //double Calculate_Signal(void);
};
#endif
