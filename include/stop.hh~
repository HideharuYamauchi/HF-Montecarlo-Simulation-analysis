///////////////////////////////////////////                
//  High field simulation(MuSEUM project)                               
//                                                  
//      Hideharu Yamauchi 2021/09/19
///////////////////////////////////////////                                    
#ifndef ___header_muonstopping_
#define ___header_muonstopping_ 1

#include <string>
#include <cstring>
#include "TString.h"
#include <vector>
#include "./HFgeometry.hh"
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TGraph.h"
#include "TCanvas.h"

class muonstopping{
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
  
public:
  Int_t entries;
  Int_t nbranches;
  muonstopping(std::string runfile, const char* envfile);
  ~muonstopping(void){delete tree;}
  void CreateRootFile(void);
  void Vis_stopping_distZ();
  void Vis_stopping_distXY(Double_t posZ);
  int GetNumber(Double_t* pos);
  TCanvas* c;
  TCanvas* c2;
  TH2D* dtxy;
  TH2D* dtz;
};
#endif
