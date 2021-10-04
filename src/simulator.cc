/////////////////////////////////////////////////////////////
//          High field simulation for MuSEUM Collaboration
//
//               Author: Hideharu Yamauchi 2021/09/19
/////////////////////////////////////////////////////////////
#ifndef ___class_simulator_
#define ___class_simulator_ 1

#include <stdio.h>
#include <string>
#include <math.h>
#include "../include/simulator.hh"
#include "TStyle.h"
#endif

simulator::simulator(TTree* decaytree, int mode):Non(0.),Noff(0.){
  Mode = mode;
  // initial state amplitude from MuSEUM technical note (2.12)
  if(mode==110){
    state_amp[0]=0.25*(1+polarization);
    state_amp[1]=0.25*(1-(pow(coefficient_c,2.)-pow(coefficient_s,2.))*polarization);
  }else if(mode==210){
    state_amp[0]=0.25*(1-polarization);
    state_amp[1]=0.25*(1+(pow(coefficient_s,2.)-pow(coefficient_s,2.))*polarization);
  }
  std::cout << "state amplitue(t=0) of 12transition=" << state_amp[0] << "\n"
            << "state amplitue(t=0) of 34transition=" << state_amp[1] << std::endl;

  decaytree->SetBranchAddress("decaytime",&decaytime);
  decaytree->SetBranchAddress("decaypositionx",&decaypositionx);
  decaytree->SetBranchAddress("decaypositiony",&decaypositiony);
  decaytree->SetBranchAddress("decaypositionz",&decaypositionz);
  decaytree->SetBranchAddress("magnet_field",&magnet_field);
  decaytree->SetBranchAddress("b",&b);
  decaytree->SetBranchAddress("RF",&RF);
  decaytree->SetBranchAddress("Effective_RF",&Effective_RF);
  decaytree->SetBranchStatus("*",1);
  entries = decaytree->GetEntries();
}

void simulator::timedev(double t, double delta, double gamma, double position[3]){
  state_amp[0]=(pow(state_amp[0],2.)*(pow(std::cos(0.5*gamma*t),2.)+pow((delta/gamma)*std::sin(0.5*gamma*t),2.))+pow(state_amp[1]*2*b*std::sin(0.5*gamma*t)/gamma,2.)+state_amp[0]*state_amp[1]*2*delta*b*pow(std::sin(0.5*gamma*t)/gamma,2.))*std::exp(-muon_life*t);
  state_amp[1]=(pow(state_amp[0]*2*b*std::sin(0.5*gamma*t)/gamma,2.)+pow(state_amp[1],2.)*(pow(std::cos(0.5*gamma*t),2.)+pow(delta*std::sin(0.5*gamma*t)/gamma,2.))-state_amp[0]*state_amp[1]*2*delta*b*pow(std::sin(0.5*gamma*t)/gamma,2.))*std::exp(-0.5*muon_life*t);
  std::cout << "state amplitue(t=" << t << ") of 12transition=" << state_amp[0] << "\n"
            << "state amplitue(t=" << t << ") of 34transition=" << state_amp[1] << std::endl;
}

void simulator::Vis_state_amp(void){
  f1 = new TF1("amplitude","",0,10,2);
  //f1->SetParameters(2,1);             
  //f1->SetParNames("constant","coefficient");                                                      
  //f1->Draw();
  delete f1;
}

double simulator::ConventionalSignal(double position[3]){
  //L = ;
  //return Non/Noff-1;
  return 0;
}

double simulator::OldMuoniumSignal(double position[3], Double_t windowopen, Double_t timewindow){
  //L =;
  //return Non/Noff-1;
  return 0;
}
