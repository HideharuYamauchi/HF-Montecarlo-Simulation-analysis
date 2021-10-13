///////////////////////////////////////////////////////
//   High field simulation for MuSEUM Collaboration 
//
//        Author : Hideharu Yamauchi 2021/09/18                
///////////////////////////////////////////////////////
#if !defined(___header_magfield_)
#define ___header_magfield_ 1

#include <cmath>
#include <vector>
#include "HFgeometry.hh"
#include "TH2.h"
#include "TCanvas.h"
#include "TGraph2D.h"
#include "TStyle.h"

class magfield{
private:
  const int moment_num = 3602;
  int mode;
  std::vector<std::vector<double>> moment_coordinate;
  std::vector<std::vector<double>> moment;
  std::vector<std::vector<double>> distance;
  std::vector<double> interval;
  std::vector<double> position;
  const double DSV = 300.; // mm
  
public:
  magfield(const char* magnetfile, int Mode);
  ~magfield(void){;};
  double GetDistance(double x, double y, double z);
  double GetBfieldValue(void);
  void Vis_magfield(double Z);
  //TTree* AddMagnetBranch(TTree* decaytree);
  TCanvas* c;
  TPad* center_pad;
  TPad* top_pad;
  TGraph2D* dt;
  TH2D* dt2;
  const Double_t B_ave = 1.199671277270803; // average field
  const Double_t scaling_factor = B_cons/B_ave;
};
#endif
