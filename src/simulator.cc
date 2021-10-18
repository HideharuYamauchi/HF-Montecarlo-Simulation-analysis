/////////////////////////////////////////////////////////////
//      High field simulation for MuSEUM Collaboration
//
//            Author: Hideharu Yamauchi 2021/09/19
/////////////////////////////////////////////////////////////
#ifndef ___class_simulator_
#define ___class_simulator_ 1

#include <stdio.h>
#include <string>
#include <time.h>
#include "../include/simulator.hh"
#include "RF.cc"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TGraphErrors.h"

SIMULATOR::SIMULATOR(const char* rootfile)
  : myStringVec(0),myMuonVec(0),myMuonDispersion(0),myPositronVec(0),myPositronDispersion(0),myField(0),myAmp(0),myAngleVec(0),
    Non(0), scan_range(400), scan_points(40), scan_step(scan_range/scan_points), signal(0.), position(3), cos_solidangle(0.), solidangle(0.)
#ifndef ___header_simulator_
#define ___header_simulator_ 1
  ,myStringVec_branch(0),myMuonVec_branch(0),myMuonDispersionVec_branch(0),myPositronVec_branch(0),myPositronDispersionVec_branch(0),myField_branch(0),myAmp_branch(0),AngleBranch(0)
#endif
{
  std::string myString(rootfile);
  run_num = myString.substr(myString.find("run"), myString.find(".root")-myString.find("run"));
  
  myFile = new TFile(rootfile,"");
  if(!myFile){
    std::cout << "Error opening " << rootfile << std::endl;
    exit(-1);
  }
  
  myTree = (TTree*)myFile->Get("DecayTree");
  //myTree->SetBranchAddress("str_vec",&myStringVec);
  myTree->SetBranchAddress("muon_vec",&myMuonVec);
  myTree->SetBranchAddress("muon_dispersion",&myMuonDispersion);
  myTree->SetBranchAddress("positron_vec",&myPositronVec);
  myTree->SetBranchAddress("positron_dispersion",&myPositronDispersion);
  myTree->SetBranchAddress("field",&myField);
  myTree->SetBranchAddress("state_amp",&myAmp);
  myTree->SetBranchAddress("Angle",&myAngleVec);
  entries = myTree->GetEntries();
  Noff = entries;
#ifndef ___header_simulator_
#define ___header_simulator_ 1
  AngleBranch = myTree->Branch("Angle", &angle_vec);
  
  for(int i=0; i<entries; i++){
    myTree->GetEntry(i);
    position[0] = (*myPositronVec)[1];                                                              
    position[1] = (*myPositronVec)[2];
    position[2] = (*myPositronVec)[3];
    CalculateAngle();
    angle_vec[0] = cos_solidangle;
    angle_vec[1] = solidangle;
    AngleBranch->Fill();
  }
  myTree->Write("", TObject::kOverwrite);
  myTree->Print();
  myTree->Scan("*");
#endif
}

SIMULATOR::~SIMULATOR(void){
  myFile->Close();
  std::cout << "Finish simulator of "<< run_num << "..." << "\n"
	    << "See you next time."<< std::endl;
}
/*
void SIMULATOR::timedev(Double_t delta){
  Amplitude[0]=(pow(Amplitude[0],2.)*(pow(std::cos(0.5*gamma*t),2.)+pow((delta/gamma)*std::sin(0.5*gamma*t),2.))+pow(Amplitude[1]*2*b*std::sin(0.5*gamma*t)/gamma,2.)+Amplitude[0]*Amplitude[1]*2*delta*b*pow(std::sin(0.5*gamma*t)/gamma,2.))*std::exp(-muon_life*t);
  Amplitude[1]=(pow(Amplitude[0]*2*b*std::sin(0.5*gamma*t)/gamma,2.)+pow(Amplitude[1],2.)*(pow(std::cos(0.5*gamma*t),2.)+pow(delta*std::sin(0.5*gamma*t)/gamma,2.))-Amplitude[0]*Amplitude[1]*2*delta*b*pow(std::sin(0.5*gamma*t)/gamma,2.))*std::exp(-0.5*muon_life*t);
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
void SIMULATOR::CalculateAngle(){
  Double_t x = -119.5, y = 119.5;
  Double_t ratio = (DetectorD_center-position[2])/(cavity_downfoil_center-position[2]);
  std::vector<Double_t> positron_at_foil(2);
  Double_t r, R;
  Double_t norm;
  std::vector<Double_t> distance(3);
  std::vector<Double_t> basic_vector(3);
  cos_solidangle = 0.;
  solidangle = 0.;
  
  for(int i=0; i<240; i++){ // get detector's x position
    x = -119.5+i;  
    for(int n=0; n<240; n++){ // get detector's y position
      y = 119.5-n;
      
      distance[0] = x-position[0];
      distance[1] = y-position[1];
      distance[2] = DetectorD_center-position[2];
      norm = sqrt(pow(distance[0],2.) + pow(distance[1],2.) + pow(distance[2],2.));
      basic_vector[0] = distance[0]/norm;
      basic_vector[1] = distance[1]/norm;
      basic_vector[2] = distance[2]/norm;
      positron_at_foil[0] = position[0]+basic_vector[0]*ratio;
      positron_at_foil[1] = position[1]+basic_vector[1]*ratio;
      r = sqrt(pow(positron_at_foil[0],2.)+pow(positron_at_foil[1],2.));

      if(0<=r&&r<=cavity_radius*1.e+3){ // discriminate projection position is on cavity foil or not
	R = sqrt(pow(distance[0],2.)+pow(distance[1],2.)+pow(distance[2],2.));
	cos_solidangle += (distance[2]/R)*DetectorD_center*pow(pow(R,2.),-1.5);
	solidangle += DetectorD_center*pow(pow(R,2.),-1.5);
      }
      else if(cavity_radius*1.e+3<r){
	cos_solidangle += 0.;
	solidangle += 0.;
      }
    }
  }
}

Double_t SIMULATOR::Calculate_g(Double_t Gamma, Double_t t){
  Double_t g = std::cos(Gamma*t)-Gamma*std::sin(Gamma*t)/gamma;
  return g;
}

Double_t SIMULATOR::ConventionalSignal(Double_t power, Double_t detuning, Double_t windowopen, Double_t cos_solid_angle, Double_t solid_angle){
  Double_t Gamma_square = pow(detuning, 2.)+4*pow(power, 2.); // MuSEUM technical note (2.50)
  Double_t L = 2*pow(power, 2.)/(Gamma_square + pow(gamma, 2.));
  
  Double_t probability
    = A[0]*(cos_solid_angle)*(polarization*L)/((A[0]*(cos_solid_angle)*polarization+A[1]*(solid_angle))*(0-std::exp(-1*gamma*windowopen)));
  if(false) std::cout << "probability:" << probability << std::endl;
  return probability;
}

Double_t SIMULATOR::OldMuoniumSignal(Double_t power, Double_t detuning, Double_t windowopen, Double_t windowclose, Double_t cos_solid_angle, Double_t solid_angle, bool flag){
  Double_t Gamma_square = pow(detuning, 2.)+4*pow(power, 2.); // MuSEUM technical note (2.50)
  Double_t L = 2*pow(power, 2.)*(
				 std::exp(-1*gamma*windowopen)*(1-Calculate_g(sqrt(Gamma_square), windowopen)*pow(gamma,2.)/(Gamma_square+pow(gamma,2.)))
				 -std::exp(-1*gamma*windowclose)*(1-Calculate_g(sqrt(Gamma_square), windowclose)*pow(gamma,2.)/(Gamma_square+pow(gamma,2.)))
				 )/Gamma_square;

  Double_t probability
    = A[0]*(cos_solid_angle)*(polarization*L)/((A[0]*(cos_solid_angle)*polarization+A[1]*(solid_angle))*(std::exp(-1*gamma*windowclose)-std::exp(-1*gamma*windowopen)));
  if(flag) std::cout << "probability:" << probability << std::endl;
  return probability;
}

void SIMULATOR::CalculateSignal(Int_t minutes=20){
  TCanvas* c = new TCanvas("c","c", 900, 900);
  c->SetGrid();
  gPad->SetRightMargin(0.03); // let the figure more right 
  TGraphErrors* curve = new TGraphErrors();
  
  char date[64];
  time_t t = time(NULL);
  Double_t y; // positron_energy/positron_max_energy
  Double_t detuning;
  Double_t error;
  Double_t ppm = 1.0e-6;
  char method;
  Int_t gates;
  
  do{
    std::cout << "Conventional[c/0] or OldMuonium[o/1]:" << std::endl;
    std::cin >> method;
    if(minutes!=20){ // determine time for one point of resonance curve, default is 20 mins
      std::cout << "How many minutes for one point:" << std::endl;
      std::cin >> minutes;
    }
  }while((method!='c'&&method!='C'&&method!='o'&&method!='O'&&method!='0'&&method!='1')||minutes<0);
  
  gates = 1500*minutes*0.5; // 60(sec)*25(gates/sec) = 1500 gates for one minute, 0.5 is for beam off
  std::cout << gates << " is for BeamOn." << "\n"
	    << 1500*minutes-gates << " is for BeamOff." << "\n"
	    << std::string(40, '*') << "\n"
	    << "RUN START: " << ctime(&t) << std::endl;
  
  for(int w=0; w<scan_points; w++){
    Non = 0;
    Noff = entries;
    detuning = -1*scan_range*0.5 + scan_step*w;
    std::cout << "START detuning " << detuning << "[/kHz]..." << "\n"
	      << "Elapsed Time since detuning "<< detuning << "[/kHz] starts...." << std::endl;
    for(int p=0; p<gates; p++){
      if((p+1)%7500==0) std::cout << minutes*0.5+5*(p+1)/7500 << "[mins]" << std::endl;
      for(int i=0; i<entries; i++){
	myTree->GetEntry(i);
	if(35<=(*myPositronDispersion)[0]*1.0e-3){ // cut positrons below threshold energy, 35 MeV for liu
	  y = (*myPositronDispersion)[0]*1.0e-3/positron_max_energy;
	  A[0] = 2*y-1;
	  A[1] = 3-2*y;
	  if(method=='c'||method=='C'||method=='0') Non += ConventionalSignal((*myField)[2], // power
									      detuning,
									      0., // windowopen 
									      (*myAngleVec)[0], // cos_solid_angle
									      (*myAngleVec)[1]); // solid_angle
	  else if(method=='o'||method=='O'||method=='1') Non += OldMuoniumSignal((*myField)[2], // power
										 detuning,
										 3.12*ppm, // windowopen 3.12 or 6.92, same with liu
										 4.07*ppm, // windowclose 4.07 or 7.87, same with liu
										 (*myAngleVec)[0], // cos_solid_angle
										 (*myAngleVec)[1], // solid_angle
										 false); // show the value
	}
      }
    }
    Noff = gates*Noff;
    signal = Non/Noff-1;
    std::cout << "Non: " << Non << "\t" << "Noff: " << Noff << "\n" 
	      << "SIGNAL INTENSITY(=Non/Noff-1): " << signal << "\n"
	      << "TOTAL ELAPSED TIME SINCE RUN START: " << std::fixed << std::setprecision(6) << (w+1)*minutes/60  << "[hours] \n" << std::endl;
    error = (Non/Noff)*TMath::Sqrt(1.0/Non+1.0/Noff);
    curve->SetPoint(w, detuning, signal);
    curve->SetPointError(w, 0, error);
  }
  std::cout << "RUN FINISH: " << ctime(&t) << std::string(40, '*') << std::endl;
  curve->GetXaxis()->SetTitle("Frequency Detuning [/kHz]");
  curve->GetYaxis()->SetTitle("Signal");
  curve->GetYaxis()->SetTitleOffset(1.4);
  curve->Draw("AP");
  c->SaveAs(("../figure/"+run_num+".png").c_str());
  delete curve;
  delete c;
}
#endif
