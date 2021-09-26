///////////////////////////////////////////
//  High field simulation(MuSEUM project)                      
//
//  Hideharu Yamauchi 2021/09/19                     
///////////////////////////////////////////

#ifndef ___class_simulator_
#define ___class_simulator_ 1

#include <stdio.h>
#include <string>
#include <cmath>
#include <vector>
#include "TString.h"
#include "TMath.h"
#include "TEntryList.h"
#include "../include/stop.hh"
#include "../include/HFgeometry.hh"
#include "../include/simulator.hh"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TPaveLabel.h"
#include "TStyle.h"
#endif

/*
simulator::simulator(TTree* runtree, int mode){
  ;
}
*/

double* simulator::timedev(double t, double b, double delta, double gamma){
  if(Mode==110){
    state_amp[0]=0.;
    state_amp[1]=0.5;
  }else if(Mode==210){
    state_amp[0]=0.5;
    state_amp[1]=0.;
  }
  state_amp[0]=(pow(state_amp[0],2.)*(pow(std::cos(0.5*gamma*t),2.)+pow((delta/gamma)*std::sin(0.5*gamma*t),2.))+pow(state_amp[1]*2*b*std::sin(0.5*gamma*t)/gamma,2.)+state_amp[0]*state_amp[1]*2*delta*b*pow(std::sin(0.5*gamma*t)/gamma,2.))*std::exp(-muon_life*t);
  state_amp[1]=(pow(state_amp[0]*2*b*std::sin(0.5*gamma*t)/gamma,2.)+pow(state_amp[1],2.)*(pow(std::cos(0.5*gamma*t),2.)+pow(delta*std::sin(0.5*gamma*t)/gamma,2.))-state_amp[0]*state_amp[1]*2*delta*b*pow(std::sin(0.5*gamma*t)/gamma,2.))*std::exp(-0.5*muon_life*t);
  std::cout << "state amplitue(1,3)=" << state_amp[0] << "\t"
            << "state amplitue(2,4)=" << state_amp[1] << std::endl;
  return state_amp;
}

void simulator::Vis_state_amp(void){
  f1 = new TF1("amplitude","myfunction",0,10,2);    
  //f1->SetParameters(2,1);                                                                                                                                                                               
  //f1->SetParNames("constant","coefficient");                                                                                                                                                           
  //f1->Draw();                                                                                                                                                                                          
}
