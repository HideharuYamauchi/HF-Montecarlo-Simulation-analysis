//////////////////////////////////////////////////////////
//     High field simulation for MuSEUM Collaboration
//
//         Author : Hideharu Yamauchi 2021/10/13
/////////////////////////////////////////////////////////
#ifndef ___header_maketree2_
#define ___header_maketree2_ 1

#include "TTree.h"
#include "RF.hh"
#include "magnet.hh"
#include "HFgeometry.hh"
#include <vector>
#include "TFile.h"
#include "TSystem.h"

class MAKETREE{
private:
  MAGNETFIELD* magnet;
  RFFIELD* RF;
  TFile* file;
  TTree* DecayTree;
  FileStat_t info;
  bool flag;
  Double_t decaytime;
  TBranch* decaytime_branch;
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
  Double_t positron_energy;
  TBranch* positron_energy_branch;
  Double_t muon_energy;
  TBranch* muon_energy_branch;
  
public:
  std::vector<std::string> str_vec;
  TBranch* str_branch;
  std::vector<Double_t> muon_vec;
  TBranch* muon_vec_branch; 
  std::vector<Double_t> muon_dispersion;
  TBranch* muon_dispersion_branch;
  std::vector<Double_t> positron_vec;
  TBranch* positron_vec_branch;
  std::vector<Double_t> positron_dispersion;                    
  TBranch* positron_dispersion_branch;
  std::vector<Double_t> field; // magnet field, RF field, effective RF field
  TBranch* field_branch;
  std::vector<Double_t> state_amp;
  TBranch* state_amp_branch;
  std::vector<Double_t> angle_vec;
  TBranch* angle_branch;
  std::vector<Double_t> position;
  MAKETREE(TTree* decaytree, int mode, std::string run_num);
  ~MAKETREE(void);
  void CalculateAngle(void);
  Double_t cos_solidangle;
  Double_t solidangle;
};
#endif
