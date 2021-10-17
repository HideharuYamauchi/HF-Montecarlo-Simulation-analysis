/////////////////////////////////////////////////////////
//      Header file for calculating the RF field
//      
//        Author: Hideharu Yamauchi 2021/09/16
/////////////////////////////////////////////////////////
#ifndef ___header_RFfield_
#define ___header_RFfield_ 1

#include <cmath>
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TGraph2D.h"
#include "TStyle.h"
#include "HFgeometry.hh"

class RFFIELD{
public:
  RFFIELD(Int_t Mode);
  ~RFFIELD(void){;};
  Double_t GetXY(int x, int y);
  Double_t TM_mode(void);
  void Vis_RF(void);
  Double_t GetEffectivePower(TH2D* xy_dist);
  
private:
  TString title;
  TString title2;
  Int_t mode;
  const Double_t e = 8.854187817e-12; // dielectric constant for krypton
  const Double_t C = 2.99792458e+8; // speed of light
  const Double_t permeability = pi*4e-7; // permeability for krypton
  const Double_t j_11 = 3.831705970207512315614;
  const Double_t j_21 = 5.135622301840682556301;
  Double_t b; // the coefficient to change RF field to frequency
  Double_t kc; 
  Double_t Kr_freq;
  Double_t Bfield;
  Double_t H_coefficient;
  Double_t angle;
  Double_t distance;
};
#endif
