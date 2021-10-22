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
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TRandom3.h"

SIMULATOR::SIMULATOR(const char* rootfile)
  : myMuonVec(0),myMuonDispersion(0),myPositronVec(0),myPositronDispersion(0),myField(0),myAmp(0),myAngleVec(0),
    Non(0), scan_range(800), scan_step(10), scan_points(81), signal(0.), power_mean(0.)
#ifndef ___header_simulator_
  ,myMuonVec_branch(0),myMuonDispersionVec_branch(0),myPositronVec_branch(0),myPositronDispersionVec_branch(0),myField_branch(0),myAmp_branch(0),AngleBranch(0)
#endif
{
  gStyle->SetOptFit(1111);
  std::string myString(rootfile);
  run_num = myString.substr(myString.find("run"), myString.find(".root")-myString.find("run"));
  
  myFile = new TFile(rootfile,"");
  if(!myFile){
    std::cout << "Error opening " << rootfile << std::endl;
    exit(-1);
  }
  
  myTree = (TTree*)myFile->Get("DecayTree");
  myTree->SetBranchAddress("muon_vec",&myMuonVec);
  myTree->SetBranchAddress("muon_dispersion",&myMuonDispersion);
  myTree->SetBranchAddress("positron_vec",&myPositronVec);
  myTree->SetBranchAddress("positron_dispersion",&myPositronDispersion);
  myTree->SetBranchAddress("field",&myField);
  myTree->SetBranchAddress("state_amp",&myAmp);
  myTree->SetBranchAddress("Angle",&myAngleVec);
  entries = myTree->GetEntries();
  
  std::string myTreeName = myTree->GetName();
  myTreeTitle = myTree->GetTitle();
  tree_TMmode = myTreeTitle.substr(myTreeTitle.find("TM"), 5);
  tree_Pressure = myTreeTitle.substr(myTreeTitle.find("atmosphere")-4, 14);
  tree_Temperature = myTreeTitle.substr(myTreeTitle.find("kelvin")-4, 10);

  int l=0;
  for(int i=0; i<entries; i++){
    myTree->GetEntry(i);
    if(35<=(*myPositronDispersion)[0]*1.0e-3){
      power_mean += (*myField)[2];
      l++;
    }
    Amplitude[0] += (*myAmp)[0];
    Amplitude[1] += (*myAmp)[1];
    Amplitude[2] += (*myAmp)[2];
    Amplitude[3] += (*myAmp)[3];
  }
  
  power_mean = power_mean/l;
  Amplitude[0] = Amplitude[0]/entries;
  Amplitude[1] = Amplitude[1]/entries;
  Amplitude[2] = Amplitude[2]/entries;
  Amplitude[3] = Amplitude[3]/entries;
  std::cout << "microwave power average is " << power_mean << "[/kHz]" << "\n"
	    << "Amplitude of |1> = " << Amplitude[0] << "\n"
	    << "Amplitude of |2> = " <<	Amplitude[1] << "\n"
	    << "Amplitude of |3> = " <<	Amplitude[2] << "\n"
	    << "Amplitude of |4> = " <<	Amplitude[3] << std::endl;
}

SIMULATOR::~SIMULATOR(void){
  myFile->Close();
  std::cout << "Finish simulator of "<< run_num << std::endl;
}

Double_t SIMULATOR::TimeDev(Double_t t, Double_t detune){
  Double_t GAMMA = sqrt(4*pow(pi,2.)*pow(detune, 2.)+4*pow(power_mean, 2.));

  Amplitude[0]=(pow(Amplitude[0],2.)*(pow(std::cos(0.5*GAMMA*t),2.)+pow((detune/GAMMA)*std::sin(0.5*GAMMA*t),2.))+pow(Amplitude[1]*2*power_mean*std::sin(0.5*GAMMA*t)/GAMMA,2.)+Amplitude[0]*Amplitude[1]*2*detune*power_mean*pow(std::sin(0.5*GAMMA*t)/GAMMA,2.))*std::exp(-gamma*t);
  Amplitude[1]=(pow(Amplitude[0]*2*power_mean*std::sin(0.5*GAMMA*t)/GAMMA,2.)+pow(Amplitude[1],2.)*(pow(std::cos(0.5*GAMMA*t),2.)+pow(detune*std::sin(0.5*GAMMA*t)/GAMMA,2.))-Amplitude[0]*Amplitude[1]*2*detune*power_mean*pow(std::sin(0.5*GAMMA*t)/GAMMA,2.))*std::exp(-gamma*t);
 
  return Amplitude[0];
}
/*
void SIMULATOR::Vis_StateAmp(Double_t detune){
  TCanvas* c = new TCanvas("c","c",900,900);
  TF1* f1 = new TF1("amplitude", TimeDev, 0, 12000, 1);
  f1->SetParNames("detune");
  f1->SetParameter("detune", detune);  
  f1->Draw();
  c->SaveAs("../figure/amplitude.png");
  delete f1;
  delete c;
}
*/
Double_t SIMULATOR::ConventionalSignal(Double_t power, Double_t detuning, Double_t windowopen, Double_t cos_solid_angle, Double_t solid_angle, bool flag){
  Double_t GAMMA = sqrt(4*pow(pi,2.)*pow(detuning, 2.)+4*pow(power, 2.)); // MuSEUM technical note (2.50)
  Double_t L = 2*pow(power, 2.)/(pow(GAMMA, 2.) + pow(gamma, 2.));

  Double_t sig;
  sig = A[0]*(cos_solid_angle)*(polarization*L)/((A[0]*(cos_solid_angle)*polarization+A[1]*(solid_angle))*(0-std::exp(-1*gamma*windowopen)));
  
  if(flag) std::cout << "L:" << L << "\t"
		     << "signal:" << sig << "\t"
		     << "energy:" << y << std::endl;
  return sig;
}

Double_t SIMULATOR::OldMuoniumSignal(Double_t power, Double_t detuning, Double_t windowopen, Double_t windowclose, Double_t cos_solid_angle, Double_t solid_angle, bool flag){
  Double_t GAMMA = sqrt(4*pow(pi,2.)*pow(detuning, 2.)+4*pow(power, 2.)); // MuSEUM technical note (2.50)
  Double_t g_open = std::cos(GAMMA*windowopen)-GAMMA*std::sin(GAMMA*windowopen)/gamma;
  Double_t g_close = std::cos(GAMMA*windowclose)-GAMMA*std::sin(GAMMA*windowclose)/gamma;
  
  Double_t L = 2*pow(power, 2.)*(
				 std::exp(-1*gamma*windowopen)*(1-g_open*pow(gamma,2.)/(pow(GAMMA,2.)+pow(gamma,2.)))
				 -std::exp(-1*gamma*windowclose)*(1-g_close*pow(gamma,2.)/(pow(GAMMA,2.)+pow(gamma,2.)))
				 )/pow(GAMMA,2.);

  Double_t sig;
  sig  = A[0]*(cos_solid_angle)*(polarization*L)/((A[0]*(cos_solid_angle)*polarization+A[1]*(solid_angle))*(std::exp(-1*gamma*windowclose)-std::exp(-1*gamma*windowopen)));
  
  if(flag) std::cout << "L:" << L << "\t"
		     << "signal:" << sig << "\t"
		     << "energy:" << y << std::endl;
  return sig;
}

void SIMULATOR::CalculateSignal(Int_t minutes=20){
  TCanvas* c = new TCanvas("c","c", 900, 900);
  c->SetGrid();
  gPad->SetRightMargin(0.03); // let the figure more right
  TGraphErrors* curve = new TGraphErrors();
  TF1* f1;
  
  time_t t = time(NULL);
  Double_t detuning;
  Double_t error;
  Double_t ppm = 1.0e-6;
  char method;
  Int_t gates;
  Double_t counton;
  Double_t countoff;
  //gRandom = new TRandom3(time(NULL));
  
  do{
    std::cout << "Conventional[c/C] or OldMuonium[o/O]:" << std::endl;
    std::cin >> method;
  }while((method!='c'&&method!='C'&&method!='o'&&method!='O')||minutes<0);
  
  gates = 1500*minutes*0.5; // 60(sec)*25(gates/sec) = 1500 gates for one minute, 0.5 is for beam off
  std::cout << gates << "[/pulse] is for BeamOn." << "\n"
	    << 1500*minutes-gates << "[/pulse] is for BeamOff." << "\n"
	    << std::string(60, '*') << "\n"
	    << "RUN START: " << ctime(&t) << std::endl;
  
  for(int w=0; w<scan_points; w++){
    Non = 0;
    Noff = 0;
    detuning = -1*scan_range*0.5 + scan_step*w;
    std::cout << "START detuning " << detuning << "[/kHz]..." << std::endl;
    for(int i=0; i<entries; i++){
      myTree->GetEntry(i);
      if(35<=(*myPositronDispersion)[0]*1.0e-3){ // cut positrons below threshold energy, 35 MeV for liu
	y = (*myPositronDispersion)[0]*1.0e-3/positron_max_energy;
	A[0] = 2*y-1;
	A[1] = 3-2*y;
	if(method=='c'||method=='C'){
	  Non += ConventionalSignal((*myField)[2], // power
				    detuning, // frequency detune
				    0., // windowopen 
				    (*myAngleVec)[0], // cos_solid_angle
				    (*myAngleVec)[1], // solid_angle
				    false); // show the value
	  Noff++;
	}
	else if(method=='o'||method=='O'){
	  Non += OldMuoniumSignal((*myField)[2], // power
				  detuning, // frequency detune
				  3.12*ppm*1.0e+3, // windowopen 3.12 or 6.92, 1.0e+3 is for convert s to ms
				  4.07*ppm*1.0e+3, // windowclose 4.07 or 7.87, 
				  (*myAngleVec)[0], // cos_solid_angle
				  (*myAngleVec)[1], // solid_angle
				  false); // show the value
	  Noff++;
	}
      }
    }
    Non = Non*gates;
    Noff = Noff*gates;
    signal = Non/Noff;
    /*
    std::cout << "Non: " << Non << "\t" << "Noff: " << Noff << "\t" << "SIGNAL INTENSITY(=Non/Noff-1): " << signal << "\n"
	      << "TOTAL ELAPSED TIME SINCE RUN START: " << std::fixed << std::setprecision(6) << (w+1)*minutes/60  << "[hours] \n" << std::endl;
    */
    error = (Non/Noff)*TMath::Sqrt(1.0/Non+1.0/Noff);
    curve->SetPoint(w, detuning, signal);
    curve->SetPointError(w, 0, error);
  }
  
  t = time(NULL);
  std::cout << "RUN FINISH: " << ctime(&t) << std::string(60, '*') << std::endl;
  
  curve->GetXaxis()->SetTitle("Frequency Detuning [/kHz]");
  curve->GetYaxis()->SetTitle("Signal");
  curve->GetYaxis()->SetTitleOffset(1.4);
  curve->Draw("AP");
  
  if(method=='c'||method=='C')
    f1 = new TF1("f1"," [0]+[4]*2*[1]*[1]/(4*TMath::Pi()*TMath::Pi()*(x-[3])*(x-[3])+4*[1]*[1]+[2]*[2]) ", -scan_range*0.5, scan_range*0.5);
  else if(method=='o'||method=='O')
    f1 = new TF1("f1",
"[0]+[4]*2*[1]*[1]*(TMath::Exp(-[2]*[5])*(1-(TMath::Cos((TMath::Sqrt(4*TMath::Pi()*TMath::Pi()*(x-[3])*(x-[3])+4*[1]*[1]))*[5])-(TMath::Sqrt(4*TMath::Pi()*TMath::Pi()*(x-[3])*(x-[3])+4*[1]*[1]))*TMath::Sin((TMath::Sqrt(4*TMath::Pi()*TMath::Pi()*(x-[3])*(x-[3])+4*[1]*[1]))*[5])/[2])*[2]*[2]/(4*TMath::Pi()*TMath::Pi()*(x-[3])*(x-[3])+4*[1]*[1]+[2]*[2])) - TMath::Exp(-[2]*[6])*(1-(TMath::Cos(TMath::Sqrt(4*TMath::Pi()*TMath::Pi()*(x-[3])*(x-[3])+4*[1]*[1])*[6])-TMath::Sqrt(4*TMath::Pi()*TMath::Pi()*(x-[3])*(x-[3])+4*[1]*[1])*TMath::Sin(TMath::Sqrt(4*TMath::Pi()*TMath::Pi()*(x-[3])*(x-[3])+4*[1]*[1])*[6])/[2])*[2]*[2]/(4*TMath::Pi()*TMath::Pi()*(x-[3])*(x-[3])+4*[1]*[1]+[2]*[2])))/(4*TMath::Pi()*TMath::Pi()*(x-[3])*(x-[3])+4*[1]*[1])",
		 -scan_range*0.5, scan_range*0.5);
  
  f1->FixParameter(0, 0); // offset
  if(tree_TMmode=="TM110") f1->SetParameter(1, power_mean); // b12 for TM110
  else if(tree_TMmode=="TM210") f1->SetParameter(1, power_mean); // b34 for TM210
  f1->SetParameter(2, gamma); // gamma
  f1->SetParameter(3, 0); // center                                                                                                                                                                     
  f1->SetParameter(4, 1); // scaling
  f1->SetParNames("Offset", "b [/kHz]", "#gamma [/kHz]", "Center", "Scaling");
  if(method=='o'||method=='O'){
    f1->SetParameter(5, 3.12e-3); // t1, ms
    f1->SetParameter(6, 4.07e-3); // t2, ms
  }
  
  curve->Fit("f1","EM", "", -scan_range*0.5, scan_range*0.5);
  
  if(method=='c'||method=='C') c->SaveAs(("../figure/"+run_num+":"+tree_TMmode+":"+tree_Pressure+":"+"conventional"+".png").c_str());
  else if(method=='o'||method=='O') c->SaveAs(("../figure/"+run_num+":"+tree_TMmode+":"+tree_Pressure+":"+"oldmuonium"+".png").c_str());
  
  delete f1;
  delete curve;
  delete c;
}
#endif
