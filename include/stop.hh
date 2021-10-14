//////////////////////////////////////////////////////////               
//     High field simulation for MuSEUM Collaboration  
//                                                  
//         Author : Hideharu Yamauchi 2021/09/19
/////////////////////////////////////////////////////////                     
#ifndef ___header_muonstopping_
#define ___header_muonstopping_ 1

#include <string>
#include <cstring>
#include "TString.h"
#include <vector>
#include "HFgeometry.hh"
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TCanvas.h"
#include "RF.hh"
#include "magnet.hh"

class STOP{
private:
  TFile* file;
  TTree* tree;
  TBranch* branchtime; // need to initialize(nullptr) because of pointer
  Double_t time;
  TBranch* branchparticle; // need to initialize(nullptr) because of pointer
  char particle[3];
  TBranch* branchvolume; // need to initialize(nullptr) because of pointer
  char volume[27]; // av_1_impr_3_FD_S08_L1_pv_16
  TBranch* branchprocess;  
  char process[15]; 
  TBranch* branchX;
  Double_t X;
  TBranch* branchY;
  Double_t Y;
  TBranch* branchZ;
  Double_t Z;
  TBranch* branchPx;
  Double_t Px;
  TBranch* branchPy;
  Double_t Py;
  TBranch* branchPz;
  Double_t Pz;
  TBranch* branchkE;
  Double_t kE;
  TBranch* branchdepE; 
  Double_t depE;
  TBranch* branchtrack;
  Int_t ntrack;
  TBranch* branchstep;
  Int_t nstep;
  TBranch** branchcopyno;
  Int_t copyno;
  TString run_num;
  RFFIELD* RF;
  MAGNETFIELD* magnet;
  
public:
  Int_t entries;
  Int_t nbranches;
  STOP(std::string runfile, const char* envfile, Int_t mode);
  ~STOP(void);
  void CreateRootFile(void);
  TH2D* Vis_Stopping_DistZ(void);
  TH2D* Vis_Stopping_DistXY(Double_t zpoint1, Double_t zpoint2, bool saveflag);
  TTree* GetDecayTree(bool scanflag);
  void Vis_RFPowerHist(void);
  void Vis_FieldHist(void);
  void Vis_PositronEnergyHist(void);
  void Vis_PositronAngleHist(void);
};
#endif
