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
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TPaveLabel.h"
#include "TGraph2D.h"
#include "TStyle.h"
#include <string>

class RFfield{  
public:
  RFfield(int Mode);
  ~RFfield(void){;};
  double GetXYZ(int x, int y, int z);
  double TM_mode(void);
  void Vis_RF();
  
private:
  const double pi = TMath::Pi();
  const double j_11 = 3.831705970207512315614; 
  const double j_21 = 5.135622301840682556301;
  const double e = 8.854187817e-12; // dielectric constant for krypton                                                    
  const double C = 2.99792458e+8; // speed of light                                                                            
  const double permeability = pi*4e-7; // permeability for krypton
  int mode;
  int position[3];
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
};
#endif
