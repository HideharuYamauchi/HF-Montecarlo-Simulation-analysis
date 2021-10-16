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
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TGraphErrors.h"

SIMULATOR::SIMULATOR(const char* rootfile)
  : myStringVec(0),myMuonVec(0),myMuonDispersion(0),myPositronVec(0),myPositronDispersion(0),myField(0),myAmp(0),
    Non(0), scan_range(400), scan_points(40), scan_step(scan_range/scan_points), signal(0.)
{
  std::string myString(rootfile);
  run_num = myString.substr(myString.find("run"), myString.find(".root")-myString.find("run"));
  
  myFile = new TFile(rootfile,"");
  if(!myFile){
    std::cout << "Error opening " << rootfile << std::endl;
    exit(-1);
  }
  myTree = (TTree*)myFile->Get("DecayTree");

  myTree->SetBranchAddress("str_vec",&myStringVec);
  myTree->SetBranchAddress("muon_vec",&myMuonVec);
  myTree->SetBranchAddress("muon_dispersion",&myMuonDispersion);
  myTree->SetBranchAddress("positron_vec",&myPositronVec);
  myTree->SetBranchAddress("positron_dispersion",&myPositronDispersion);
  myTree->SetBranchAddress("field",&myField);
  myTree->SetBranchAddress("state_amp",&myAmp);
  entries = myTree->GetEntries();
  Noff = entries;
}

SIMULATOR::~SIMULATOR(void){
  myFile->Close();
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
Double_t SIMULATOR::Calculate_g(Double_t Gamma, Double_t t){
  Double_t g = std::cos(Gamma*t)-Gamma*std::sin(Gamma*t)/gamma;
  return g;
}

Double_t SIMULATOR::ConventionalSignal(Double_t power, Double_t detuning){
  Double_t Gamma_square = pow(detuning, 2.)+4*pow(power, 2.); // MuSEUM technical note (2.50)
  L = 2*pow(power, 2.)/(Gamma_square + pow(gamma, 2.));
  K = 1;
  std::cout << "probability:" << L*K << std::endl;
  return L*K;
}

Double_t SIMULATOR::OldMuoniumSignal(Double_t power, Double_t detuning, Double_t windowopen, Double_t windowclose){
  Double_t Gamma_square = pow(detuning, 2.)+4*pow(power, 2.); // MuSEUM technical note (2.50)
  L = 2*pow(power, 2.)*(
			std::exp(-muon_life*windowopen)*(1-Calculate_g(sqrt(Gamma_square), windowopen)*pow(muon_life,2.)/(Gamma_square+pow(muon_life,2.)))
			-std::exp(-muon_life*windowclose)*(1-Calculate_g(sqrt(Gamma_square), windowopen)*pow(muon_life,2.)/(Gamma_square+pow(muon_life,2.)))
			)/Gamma_square;
  K = 1;
  std::cout << "probability:" << L*K << std::endl;
  return L*K;
}

void SIMULATOR::CalculateSignal(Int_t minutes=20){
  TCanvas* c = new TCanvas("c","c", 900, 900);
  TGraphErrors* curve = new TGraphErrors();
  
  Double_t detuning;
  Double_t error;
  Double_t ppm = 1.0e-6;
  char method;
  Int_t gates;
  do{
    // choose method 
    std::cout << "Conventional[c/0] or OldMuonium[o/1]:" << std::endl;
    std::cin >> method;

    // determine time for one point of resonance curve, default is 20 mins
    if(minutes!=20){
      std::cout << "How many minutes for one point:" << std::endl;
      std::cin >> minutes;
    }
  }while((method!='c'&&method!='C'&&method!='o'&&method!='O'&&method!='0'&&method!='1')||minutes<0);
  
  gates = 1500*minutes*0.5; //60(sec)*25(gates/sec) = 1500 gates for one minute, 0.5 is for beam off
  std::cout << gates << " is for BeamOn." << "\n"
	    << 1500*minutes-gates << " is for BeamOff." << std::endl;
  
  for(int w=0; w<scan_points; w++){
    detuning = -1*scan_range*0.5 + scan_step*w;
    std::cout << "START detuning " << detuning << "[/kHz]..." << std::endl;
    for(int p=0; p<10/*gates*/; p++){ // BEAM ON
      if((p+1)%7500==0) std::cout << "Elapsed Time since detuning "<< detuning << "[/kHz] starts: " << minutes*0.5+5*(p+1)/7500 << "[mins]" << std::endl;
      for(int i=0; i<10/*entries*/; i++){
	myTree->GetEntry(i);
	if(method=='c'||method=='C'||method=='0') signal += ConventionalSignal((*myField)[2], detuning);
	else if(method=='o'||method=='O'||method=='1') signal += OldMuoniumSignal((*myField)[2], // power
										  detuning,
										  3.12*ppm, // windowopen 3.12 or 6.92, same with liu
										  4.07*ppm); // windowclose 4.07 or 7.87, same with liu
	std::cout << "SIGNAL INTENSITY:" << signal << std::endl; 
      }
    }
    std::cout << "TOTAL ELAPSED TIME SINCE RUN START: " << (w+1)*minutes/60  << "[hours]" << std::endl;
    error = (Non/Noff)*TMath::Sqrt(1.0/Non+1.0/Noff);
    //curve->SetPoint(w, detuning, signal);
    //curve->SetPointError(w, 0, error);
  }
  curve->GetXaxis()->SetTitle("Frequency Detuning [/kHz]");
  curve->GetYaxis()->SetTitle("Signal");
  //curve->Draw("AP");
  //c->SaveAs(("../figure/"+run_num+".png").c_str());
  delete curve;
  delete c;
}
/*
Double_t* SIMULATOR::CalculateSolid(std::vector<Double_t> pos){
  Double solid_angle[2] = {0}; // second is for cos term
  
  return solid_angle;
}
*/
#endif
