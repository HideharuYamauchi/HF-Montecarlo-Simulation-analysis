//////////////////////////////////////////////////////////
//     High field simulation for MuSEUM Collaboration
//
//         Author : Hideharu Yamauchi 2021/10/12
/////////////////////////////////////////////////////////
#ifndef ___header_maketree_
#define ___header_maketree_ 1

#include "TTree.h"
#include "RF.hh"
#include "magnet.hh"
#include "HFgeometry.hh"
#include <vector>
#include "TFile.h"
#include "TSystem.h"

class maketree{
private:
  magfield* magnet;
  RFfield* RF;
  TFile* file;
  FileStat_t info;
  bool flag;
  
public:
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
  std::vector<Double_t> field; // magnet field, RF field, effective RF field
  TBranch* field_branch;
  std::vector<Double_t> state_amp;
  TBranch* state_amp_branch;
  maketree(TTree* decaytree, int mode, std::string run_num);
  ~maketree(void);
};
#endif
