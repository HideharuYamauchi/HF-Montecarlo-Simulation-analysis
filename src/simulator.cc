/////////////////////////////////////////////////////////////
//      High field simulation for MuSEUM Collaboration
//
//            Author: Hideharu Yamauchi 2021/09/19
/////////////////////////////////////////////////////////////
#ifndef ___class_simulator_
#define ___class_simulator_ 1

#include <stdio.h>
#include <string>
#include "../include/simulator.hh"
#include "RF.cc"
#include "TStyle.h"

SIMULATOR::SIMULATOR(const char* rootfile)
  :Non(0.),Noff(0.)
{
  treefile = new TFile(rootfile,"");
  if(!treefile){
    std::cout << "Error opening " << rootfile << std::endl;
    exit(-1);
  }
  tree = (TTree*)treefile->Get("DecayTree");
  tree->SetBranchStatus("*", 1);
  tree->Scan("*");
}

SIMULATOR::~SIMULATOR(void){
  treefile->Close();
  //delete tree; // give error!!
  delete treefile;
  std::cout << "finish simulator..." << std::endl;
}
/*
Double_t SIMULATOR::timedev(Double_t delta){
  Amplitude[0]=(pow(Amplitude[0],2.)*(pow(std::cos(0.5*gamma*t),2.)+pow((delta/gamma)*std::sin(0.5*gamma*t),2.))+pow(Amplitude[1]*2*b*std::sin(0.5*gamma*t)/gamma,2.)+Amplitude[0]*Amplitude[1]*2*delta*b*pow(std::sin(0.5*gamma*t)/gamma,2.))*std::exp(-muon_life*t);
  Amplitude[1]=(pow(Amplitude[0]*2*b*std::sin(0.5*gamma*t)/gamma,2.)+pow(Amplitude[1],2.)*(pow(std::cos(0.5*gamma*t),2.)+pow(delta*std::sin(0.5*gamma*t)/gamma,2.))-Amplitude[0]*Amplitude[1]*2*delta*b*pow(std::sin(0.5*gamma*t)/gamma,2.))*std::exp(-0.5*muon_life*t);
  return Amplitude[1];
}

void SIMULATOR::Vis_State_Amp(Int_t entry){
  tree->GetEntry(entry);
  TCanvas* c = new TCanvas("c","c",900,900);
  TF1* f1 = new TF1("amplitude",timedev(),0,30000);
  f1->Draw();
  c->SaveAs("../figure/amplitude.png");
  delete f1;
}
*/
/*
double simulator::Calculate_g(Double_t Gamma, Double_t t){
  double g = std::cos(Gamma*t)-Gamma*std::sin(Gamma*t)/muon_life;
  return g;
}

double simulator::ConventionalSignal(double power, double detuning, double position[3]){
  double Gamma_square = pow(detuning, 2.)+4*pow(power, 2.);
  L = 2*pow(power, 2.)/(Gamma_square + pow(gamma, 2.));
  //K = ;
  return 0;
}

double simulator::OldMuoniumSignal(double power, double detuning, double position[3], Double_t windowopen, Double_t windowclose){
  double Gamma_square = pow(detuning, 2.)+4*pow(power, 2.);
  L = 2*pow(power, 2.)*(
			std::exp(-muon_life*windowopen)*(1-Calculate_g(sqrt(Gamma_square)*windowopen)*pow(muon_life,2.)/(Gamma_square+pow(muon_life,2.)))
			-std::exp(-muon_life*windowclose)*(1-Calculate_g(sqrt(Gamma_square)*windowopen)*pow(muon_life,2.)/(Gamma_square+pow(muon_life,2.)))
			)/Gamma_square;
  //K = ;
  return 0;
}

double simulator::Calculate_Signal(void){
  for(int i=0; i<entries; i++){
    ;
  }
  return 0;
}
*/
#endif
