////////////////////////////////////////////////////////
//    High field simulation for MuSEUM Collaboration
//
//        Author:  Hideharu Yamauchi 2021/09/18 
////////////////////////////////////////////////////////
#if !defined(___header_simulator_)
#define ___header_simulator_ 1

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "HFgeometry.hh"
#include "RF.hh" // for omega
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
  Double_t signal;
  Int_t scan_range; // kHz
  Int_t scan_points; // same with liu
  Double_t scan_step; //10 kHz
  const Double_t gamma = 1/muon_life; // muon's natural width
  Int_t entries;
  std::string run_num;
  std::vector<Double_t> position;
  
public:
  Double_t Non;
  Double_t Noff;
  Double_t signal;
  Double_t Amplitude[2];
  Double_t L;
  SIMULATOR(const char* rootfile);
  ~SIMULATOR(void);
  //Double_t* timedev(Double_t delta);
  //void Vis_State_Amp(Int_t entry);
  //Double_t Calculate_EnergySplit(void);
  Double_t CalculateAngle(void);
  Double_t Calculate_g(Double_t Gamma, Double_t t);
  Double_t ConventionalSignal(Double_t power, Double_t detuning, Double_t positron_energy);
  Double_t OldMuoniumSignal(Double_t power, Double_t detuning, Double_t windowopen, Double_t windowclose, Double_t positron_energy);
  void CalculateSignal(Int_t minutes);
};
#endif
