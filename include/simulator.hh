////////////////////////////////////////////////////////
//    High field simulation for MuSEUM Collaboration
//
//        Author:  Hideharu Yamauchi 2021/09/18 
////////////////////////////////////////////////////////
#ifndef ___header_simulator_
#define ___header_simulator_ 1

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TF1.h"
#include "HFgeometry.hh"
#include <vector>

class SIMULATOR{
private:
  TFile* myFile;
  TTree* myTree;
  std::vector<std::string>* myStringVec;
  std::vector<Double_t>* myMuonVec;
  std::vector<Double_t>* myMuonDispersion;
  std::vector<Double_t>* myPositronVec;
  std::vector<Double_t>* myPositronDispersion;
  std::vector<Double_t>* myField;
  std::vector<Double_t>* myAmp;
  std::vector<Double_t>* myAngleVec;
  Int_t scan_range; // kHz
  Int_t scan_points; // same with liu
  Double_t scan_step; //10 kHz
  Int_t minutes;
  const Double_t gamma = (1/muon_life)*1.0e-3; // muon's natural width, 1.0e-3 is for convert Hz to kHz
  const Double_t threshold;
  Int_t entries;
  std::string run_num;
  char method; // conventional method or oldmuonium method
  Int_t plot_method;
  char fit_gamma;
    
public:
  Double_t Non;
  Double_t Noff;
  Double_t Pon;
  Double_t Poff;
  Double_t signal;
  Double_t power_mean;
  Double_t cos_solid_angle_mean;
  Double_t solid_angle_mean;
  Double_t Amplitude[4] = {0.};
  Double_t A[2] = {0.};
  Double_t The_detected;
  Double_t Sim_detected;
  Double_t y; // positron_energy/positron_max_energy
  std::string myTreeTitle;
  std::string tree_TMmode;
  std::string tree_Pressure;
  std::string tree_Temperature;
  Double_t Sim_FWHM;
  Double_t The_FWHM;
  Double_t Sim_Height;
  Double_t The_Height;
  Double_t center_error;
  SIMULATOR(const char* rootfile, Int_t run_time);
  ~SIMULATOR(void);
  Double_t TimeDev(Double_t t, Double_t detune);
  //void Vis_StateAmp(Double_t detune);
  Double_t ConventionalSignal(Double_t power, Double_t detuning, Double_t windowopen, Double_t cos_solid_angle, Double_t solid_angle, bool flag);
  Double_t OldMuoniumSignal(Double_t power, Double_t detuning, Double_t windowopen, Double_t windowclose, Double_t cos_solid_angle, Double_t solid_angle, bool flag);
  void CalculateSignal(bool FWHM_falg, Double_t bp, Double_t start, Double_t end);
  void CalculateFWHM(void);
};
#endif
