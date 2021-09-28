///////////////////////////////////////////
//          High field simulation                      
//
//  Author: Hideharu Yamauchi 2021/09/19                     
///////////////////////////////////////////

#ifndef ___class_simulator_
#define ___class_simulator_ 1

#include <stdio.h>
#include <string>
#include <math.h>
#include "TString.h"
#include "../include/simulator.hh"
#include "./RF.cc"
#include "TCanvas.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TPaveLabel.h"
#include "TStyle.h"
#endif

simulator::simulator(TTree* runtree, int mode){
  Mode = mode;
  CalculateBtob();
}

void simulator::CalculateBtob(void){
  double X = -B_cons*(gfactor_j*magnetic_moment_j + gfactor_mu_prime*magnetic_moment_mu)/(plank_const*v_exp);
  coefficient_s = sqrt(0.5)*sqrt(1-(X/(sqrt(1+X*X))));
  coefficient_c = sqrt(0.5)*sqrt(1+(X/(sqrt(1+X*X))));
  if(Mode==110){Btob = 0.001*0.25*(coefficient_s*gfactor_j*magnetic_moment_j + coefficient_c*gfactor_mu_prime*magnetic_moment_mu)/plank_const_divided;} // 0.001 is for converting Hz to kHz
  else if(Mode==210){Btob = 0.001*0.25*(-coefficient_s*gfactor_j*magnetic_moment_j + coefficient_c*gfactor_mu_prime*magnetic_moment_mu)/plank_const_divided;}
  std::cout << Mode << ":" << Btob << "(kHz/T)" <<std::endl;
}

int simulator::CalculateEffectiveb(int* position, RFfield* GetRFField){;
  int Effectiveb;
  GetRFfiled->GetXYZ(*position, *(position+1), *(position+2));
  GetRFfiled->TMmode();
  delete GetRFfield;
  return nEffectiveb;
}

void simulator::timedev(double t, double b, double delta, double gamma){
  if(Mode==110){
    state_amp[0]=0.25*(1+polarization);
    state_amp[1]=0.25*(1-(pow(coefficient_c,2.)-pow(coefficient_s,2.))*polarization);
  }else if(Mode==210){
    state_amp[0]=0.25*(1-polarization);
    state_amp[1]=0.25*(1+(pow(coefficient_s,2.)-pow(coefficient_s,2.))*polarization);
  }
  std::cout << "state amplitue(t=" << t << ") of 12transition=" << state_amp[0] << "\t"
            << "state amplitue(t=" << t << ") of 34transition=" << state_amp[1] << std::endl;
  state_amp[0]=(pow(state_amp[0],2.)*(pow(std::cos(0.5*gamma*t),2.)+pow((delta/gamma)*std::sin(0.5*gamma*t),2.))+pow(state_amp[1]*2*b*std::sin(0.5*gamma*t)/gamma,2.)+state_amp[0]*state_amp[1]*2*delta*b*pow(std::sin(0.5*gamma*t)/gamma,2.))*std::exp(-muon_life*t);
  state_amp[1]=(pow(state_amp[0]*2*b*std::sin(0.5*gamma*t)/gamma,2.)+pow(state_amp[1],2.)*(pow(std::cos(0.5*gamma*t),2.)+pow(delta*std::sin(0.5*gamma*t)/gamma,2.))-state_amp[0]*state_amp[1]*2*delta*b*pow(std::sin(0.5*gamma*t)/gamma,2.))*std::exp(-0.5*muon_life*t);
  std::cout << "state amplitue(t=" << t << ") of 12transition=" << state_amp[0] << "\t"
            << "state amplitue(t=" << t << ") of 34transition=" << state_amp[1] << std::endl;
}

void simulator::Vis_state_amp(void){
  f1 = new TF1("amplitude","myfunction",0,10,2);    
  //f1->SetParameters(2,1);                                                                                       
  //f1->SetParNames("constant","coefficient");                                                      
  //f1->Draw();
  delete f1;
}