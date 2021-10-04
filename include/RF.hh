/////////////////////////////////////////////////////////
//   Header file for calculating the RF field
//      
//         Hideharu Yamauchi 2021/09/16
/////////////////////////////////////////////////////////
#ifndef ___header_RFfield_
#define ___header_RFfield_ 1

#include "TMath.h"
#include <cmath>         
#include <math.h>
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph2D.h"
#include "TStyle.h"
#include "HFgeometry.hh"

class RFfield{  
public:
  RFfield(int Mode);
  ~RFfield(void){;};
  double GetXY(int x, int y);
  double GetXY(double x, double y);
  double TM_mode(void);
  void Vis_RF(void);
  Int_t Effective(TH2D* xy_dist);
  TTree* AddRFBranch(TTree* decaytree);
  
private:
  TString title;
  TString title2;
  const double j_11 = 3.831705970207512315614; 
  const double j_21 = 5.135622301840682556301;
  const double e = 8.854187817e-12; // dielectric constant for krypton     
  const double C = 2.99792458e+8; // speed of light    
  const double permeability = pi*4e-7; // permeability for krypton
  int mode;
  double b; // the coefficient to change RF field to frequency
  double kc; 
  double Kr_freq;
  double Bfield;
  double H_coefficient;
  double angle;
  double distance;
  TCanvas* c;
  TPad* center_pad;
  TPad* top_pad;
  TGraph2D* dt;
  TH2D* dt2;
  TH1D* hist;
};
#endif
