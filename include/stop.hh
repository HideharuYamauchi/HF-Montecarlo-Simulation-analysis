///////////////////////////////////////////                
//  High field simulation(MuSEUM project)                               
//                                                  
//  Hideharu Yamauchi 2021/09/19
///////////////////////////////////////////                                    
#ifndef ___header_muonstopping_
#define ___header_muonstopping_

#include <string>
#include "TString.h"
#include <vector>
#include "./HFgeometry.hh"
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TObjArray.h"

class muonstopping{
private:
  std::vector<double> position;
  TFile* file;
  TTree* tree;
  TObjArray* branches;
  TEntryList* entrylist;
  Int_t entries;
  Int_t branch_num;
  Int_t size;
  Double_t time;
  std::string sparticle;
  std::string svolume;
  std::string sprocess;
  char* particle;//=sparticle.c_str();
  char* volume;//=svolume.c_str();
  char* process;//=sprocess.c_str();
  Double_t X;
  Double_t Y;
  Double_t Z;
  Double_t Px;
  Double_t Py;
  Double_t Pz;
  Double_t kE;
  Double_t depE;
  Int_t track;
  Int_t step;
  Int_t copyno;
  TString run_num;
  
public:
  //muonstopping(TString runfile){tree = new TTree("tree",""); Long64_t nlines = tree->ReadFile(runfile.Data(),"time/D:particle/C:volume:process:X/D:Y:Z:Px:Py:Pz:kE:depE:track/I:step:copyno");};
  muonstopping(std::string runfile, const char* envfile);
  ~muonstopping(void){delete tree; delete[] particle; delete[]volume; delete[] process;};
  void GetScanTree(const char* varexp, const char* selection)const{tree->Scan(varexp, selection);};
  void GetBranchInfo(void)const{tree->Print();};
  //void GetBranchName(Int_t entry)const{branches->At(entry)->GetName();};
  void MakeRootFile(void){file = TFile::Open(run_num+=".root","RECREATE"); tree->Write(); file->Close(); std::cout << run_num << " is made." << std::endl;} // tree->SaveAs("hoge.root")
  double GetXYZ(double x, double y, double z);
  void Vis_stopping_distZ(void);
  void GetListOfEntry(Int_t entry)const{tree->Show(entry);};
  void Vis_stopping_distXY(double Z);
  TCanvas* c;
  TCanvas* c2;
  TH2D* dtxy;
  TH2D* dtz;
};
#endif
