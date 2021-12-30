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
#include <Fit/FitResult.h>
#include "TFitResult.h"
#include "RtypesCore.h" // for tmatrix
#include "TStyle.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TRandom3.h"

//ROOT::Fit::FitResult* FitResult = new ROOT::Fit::FitResult();

SIMULATOR::SIMULATOR(const char* rootfile, Int_t run_time=20)
  : myMuonVec(0),myMuonDispersion(0),myPositronVec(0),myPositronDispersion(0),myField(0),myAmp(0),myAngleVec(0),
    Non(0), scan_range(400), scan_step(10), scan_points(41), signal(0), power_mean(0), threshold(35.), solid_angle_mean(0), cos_solid_angle_mean(0)
{
  pola_hist->SetXTitle("P_{ON}");
  pola_hist->SetYTitle("");
  gStyle->SetOptFit(1111);
  std::string myString(rootfile);
  run_num = myString.substr(myString.find("run"), myString.find(".root")-myString.find("run"));
  minutes = run_time;
  
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

  std::cout << "total:" << entries << std::endl;
  
  std::string myTreeName = myTree->GetName();
  myTreeTitle = myTree->GetTitle();
  tree_TMmode = myTreeTitle.substr(myTreeTitle.find("TM"), 5);
  tree_Pressure = myTreeTitle.substr(myTreeTitle.find("atmosphere")-4, 14);
  tree_Temperature = myTreeTitle.substr(myTreeTitle.find("kelvin")-4, 10);

  std::cout << "Our mode is " << tree_TMmode << std::endl;

  for(int i=0; i<entries; i++){
    myTree->GetEntry(i);
    if(threshold<=(*myPositronDispersion)[0]*1.0e-3){
      power_mean += (*myField)[2];
      cos_solid_angle_mean += (*myAngleVec)[0];
      solid_angle_mean += (*myAngleVec)[1];
      l++;
    }
    Amplitude[0] += (*myAmp)[0];
    Amplitude[1] += (*myAmp)[1];
    Amplitude[2] += (*myAmp)[2];
    Amplitude[3] += (*myAmp)[3];
  }
  
  power_mean = power_mean/l;
  cos_solid_angle_mean = cos_solid_angle_mean/l;
  solid_angle_mean = solid_angle_mean/l;
  Amplitude[0] = Amplitude[0]/entries;
  Amplitude[1] = Amplitude[1]/entries;
  Amplitude[2] = Amplitude[2]/entries;
  Amplitude[3] = Amplitude[3]/entries;
  std::cout << "microwave power average is " << power_mean << "[/kHz]" << "\n"
	    << "solid angle average is " << solid_angle_mean << "\n"
	    << "cosin of solid angle average is " << cos_solid_angle_mean << "\n"
	    << "Amplitude of |1> = " << Amplitude[0] << "\n"
	    << "Amplitude of |2> = " <<	Amplitude[1] << "\n"
	    << "Amplitude of |3> = " <<	Amplitude[2] << "\n"
	    << "Amplitude of |4> = " <<	Amplitude[3] << std::endl;

  do{
    std::cout << "Conventional[c/C] or OldMuonium[o/O]:" << std::endl;
    std::cin >> method;
  }while((method!='c'&&method!='C'&&method!='o'&&method!='O')||minutes<0);

  fit_gamma='n';
  /*
  do{ // use gamma as free parameter or not
    std::cout << "Fitting using gamma? [y/Y]or[n/N]:";
    std::cin >> fit_gamma;
  }while(fit_gamma!='y'&&fit_gamma!='Y'&&fit_gamma!='n'&&fit_gamma!='N');
  */
}

SIMULATOR::~SIMULATOR(void){
  myFile->Close();
  std::cout << "Finish simulator of "<< run_num << "..." <<std::endl;
}

Double_t SIMULATOR::TimeDev(Double_t t, Double_t detune){
  Double_t GAMMA = sqrt(4*pow(pi,2.)*pow(detune, 2.)+4*pow(power_mean, 2.));

  Amplitude[0]=(pow(Amplitude[0],2.)*(pow(std::cos(0.5*GAMMA*t),2.)+pow((detune/GAMMA)*std::sin(0.5*GAMMA*t),2.))+pow(Amplitude[1]*2*power_mean*std::sin(0.5*GAMMA*t)/GAMMA,2.)+Amplitude[0]*Amplitude[1]*2*detune*power_mean*pow(std::sin(0.5*GAMMA*t)/GAMMA,2.))*std::exp(-gamma*t);
  Amplitude[1]=(pow(Amplitude[0]*2*power_mean*std::sin(0.5*GAMMA*t)/GAMMA,2.)+pow(Amplitude[1],2.)*(pow(std::cos(0.5*GAMMA*t),2.)+pow(detune*std::sin(0.5*GAMMA*t)/GAMMA,2.))-Amplitude[0]*Amplitude[1]*2*detune*power_mean*pow(std::sin(0.5*GAMMA*t)/GAMMA,2.))*std::exp(-gamma*t);
 
  return Amplitude[0];
}

Double_t SIMULATOR::ConventionalSignal(Double_t power, Double_t detuning, Double_t windowopen, Double_t cos_solid_angle, Double_t solid_angle, bool flag){
  if(flag) std::cout << std::string(60, '*') << std::endl;
  Double_t GAMMA = sqrt(4*pow(pi,2.)*pow(detuning, 2.)+4*pow(power, 2.)); // MuSEUM technical note (2.50)
  Double_t L = 2*pow(power, 2.)/(pow(GAMMA, 2.) + pow(gamma, 2.));

  Pon += 0.25*pow(y,2.)*(A[1]*solid_angle+polarization*A[0]*cos_solid_angle*( 1-L ))/pi;
  Poff += 0.25*pow(y,2.)*(A[1]*solid_angle+polarization*A[0]*cos_solid_angle)/pi;
  
  Double_t sig = A[0]*(cos_solid_angle)*(polarization*L)/((A[0]*cos_solid_angle*polarization+A[1]*solid_angle)*(0-std::exp(-1*gamma*windowopen)));
  
  if(flag) std::cout << "L:" << L << "\n"
		     << "Pon:" << Pon << "\t"
		     << "Poff:" << Poff << "\n"
		     << "signal(Pon/Poff-1):" << Pon/Poff-1 << "\t"
		     << "signal(Non/Noff-1):" << sig << std::endl;
  
  if(flag) std::cout << std::string(60, '*') << std::endl;
  
  return sig;
}

Double_t SIMULATOR::Signal2(Double_t power, Double_t detuning, Double_t Decay, Double_t cos_solid_angle, Double_t solid_angle, Double_t state13, Double_t state24, bool flag){
  if(flag) std::cout << std::string(60, '~') << std::endl;
  
  Double_t Gamma = sqrt(4*pow(pi,2.)*pow(detuning, 2.)+4*pow(power, 2.)); // MuSEUM technical note (2.50)
  Double_t L = 2*pow(power, 2.)/(pow(Gamma, 2.) + pow(gamma, 2.));
  
  Double_t pola13_on = (state13*(pow(TMath::Cos(0.5*Gamma*Decay),2.)+pow(2*TMath::Pi()*detuning*TMath::Sin(0.5*Gamma*Decay)/Gamma,2.))+state24*pow(2*power*TMath::Sin(0.5*Gamma*Decay)/Gamma,2.)+4*state13*state24*power*2*TMath::Pi()*detuning*pow(TMath::Sin(0.5*Gamma*Decay)/Gamma,2.))*TMath::Exp(-gamma*Decay);
  Double_t pola24_on = (state13*pow(2*power*TMath::Sin(0.5*Gamma*Decay)/Gamma,2.)+state24*(pow(TMath::Cos(0.5*Gamma*Decay),2.)+pow(2*TMath::Pi()*detuning*TMath::Sin(0.5*Gamma*Decay)/Gamma,2.))-4*state13*state24*power*2*TMath::Pi()*detuning*pow(TMath::Sin(0.5*Gamma*Decay)/Gamma,2.))*TMath::Exp(-gamma*Decay);
  
  Double_t pola_on = 2*(pola13_on-pola24_on);
  Double_t pola_off = 2*(state13-state24)*TMath::Exp(-gamma*Decay);

  if(flag) std::cout << "decay time: " << Decay*1.0e+3 << " [us]"<< "\n"
		     << "RF:ON, Sz=" << pola_on << "\n"
		     << "RF:OFF, Sz=" << pola_off << std::endl;
  /*
  Pz_on += pola_on;
  Pz_off += pola_off;
  */
  Pon += 0.25*y*y*(A[1]*solid_angle + A[0]*(pola_on*TMath::Exp(gamma*Decay))*cos_solid_angle)/TMath::Pi();
  Poff += 0.25*y*y*(A[1]*solid_angle + A[0]*(pola_off*TMath::Exp(gamma*Decay))*cos_solid_angle)/TMath::Pi();

  if(flag) std::cout << "Pon:" << Pon << "\t"
		     << "Poff:" << Poff << "\n"
		     << "signal(Pon/Poff-1):" << Pon/Poff-1 << std::endl;
  if(flag) std::cout << std::string(60, '~') << std::endl;

  //pola_hist->Fill(pola_on);
  return pola_on;
}

Double_t SIMULATOR::OldMuoniumSignal(Double_t power, Double_t detuning, Double_t windowopen, Double_t windowclose, Double_t cos_solid_angle, Double_t solid_angle, bool flag){
  Double_t GAMMA = sqrt(4*pow(pi,2.)*pow(detuning, 2.)+4*pow(power, 2.)); // MuSEUM technical note (2.50)
  Double_t g_open = std::cos(GAMMA*windowopen)-GAMMA*std::sin(GAMMA*windowopen)/gamma;
  Double_t g_close = std::cos(GAMMA*windowclose)-GAMMA*std::sin(GAMMA*windowclose)/gamma;
  
  Double_t L = 2*pow(power, 2.)*(
				 std::exp(-1*gamma*windowopen)*(1-g_open*pow(gamma,2.)/(pow(GAMMA,2.)+pow(gamma,2.)))
				 -std::exp(-1*gamma*windowclose)*(1-g_close*pow(gamma,2.)/(pow(GAMMA,2.)+pow(gamma,2.)))
				 )/pow(GAMMA,2.);
 
  Pon += 0.5*pow(y,2.)*(A[1]*(std::exp(-1*gamma*windowopen)-std::exp(-1*gamma*windowclose))*solid_angle
			+polarization*A[0]*cos_solid_angle*((std::exp(-1*gamma*windowopen)-std::exp(-1*gamma*windowclose))-L))/pi;
  Poff += 0.5*pow(y,2.)*(A[1]*(std::exp(-1*gamma*windowopen)-std::exp(-1*gamma*windowclose))*solid_angle
			 +polarization*A[0]*(std::exp(-1*gamma*windowopen)-std::exp(-1*gamma*windowclose))*cos_solid_angle)/pi;

  Double_t sig = A[0]*(cos_solid_angle)*(polarization*L)/((A[0]*(cos_solid_angle)*polarization+A[1]*(solid_angle))*(std::exp(-1*gamma*windowclose)-std::exp(-1*gamma*windowopen)));
  
  if(flag) std::cout << "L:" << L << "\n"
		     << "Pon:" << Pon << "\t"
                     << "Poff:" << Poff << "\n"
                     << "signal(Pon/Poff-1):" << Pon/Poff-1 << "\t"
		     << "signal(Non/Noff-1):" << sig << y << std::endl;
  return sig;
}

void SIMULATOR::CalculateSignal(bool FWHM_falg=false, Double_t bp=0., Double_t start=6.5, Double_t end=7.5){
  std::string power_str = std::to_string(bp);
  TGraphErrors* curve = new TGraphErrors();
  TF1* f1;
  
  Double_t Totalon;
  Double_t Totaloff;
  time_t t = time(NULL);
  Double_t detuning;
  Double_t error;
  Int_t gates;
  Double_t normalization;  
  TRandom3* gRandom = new TRandom3(time(NULL));
  Double_t Amp13, Amp24;

  Double_t time[l]={0}, p[l]={0};
  
  std::string start_str;
  std::string end_str;
  if(method=='o'||method=='O'){
    start_str = std::to_string(start) + "e-3";
    end_str = std::to_string(end) + "e-3";
    std::cout << "window open time;" << start_str << "[/ms]\n"
	      << "window close time:" << end_str << "[/ms]" << std::endl;
  }else if(method=='c'||method=='C') start = 0.;
  
  gates = 1500*minutes*0.5; // 60(sec)*25(gates/sec) = 1500 gates for one minute, 0.5 is for beam off
  std::cout << gates << "[/pulse] is for BeamOn." << "\n"
	    << 1500*minutes-gates << "[/pulse] is for BeamOff." << "\n"
	    << std::string(60, '*') << std::endl;
  
  Sim_detected = 0;
  for(int w=1; w<scan_points; w++){
    Non = 0;
    Noff = 0;
    Totalon = 0;
    Totaloff = 0;
    Pon = 0;
    Poff = 0;
    signal = 0;
    normalization = 0;
    detuning = -1*scan_range*0.5 + scan_step*w;
    //std::cout << "START detuning " << detuning << "[/kHz]..." << std::endl;
    for(int i=0; i<entries; i++){
      myTree->GetEntry(i);
      if(threshold<=(*myPositronDispersion)[0]*1.0e-3){ // cut positrons below threshold energy, 35 MeV for liu
	y = (*myPositronDispersion)[0]*1.0e-3/positron_max_energy;
	A[0] = 2*y-1;
	A[1] = 3-2*y;
	if(tree_TMmode=="TM110"){
	  Amp13 = (*myAmp)[0]; // state 1
	  Amp24 = (*myAmp)[1]; // state 2
	}else if(tree_TMmode=="TM210"){
	  Amp13 = (*myAmp)[2]; // state 3
	  Amp24 = (*myAmp)[3]; // state 4
	}
	
	if((method=='c'||method=='C')&&FWHM_falg==false){
	  /*
	  signal += ConventionalSignal((*myField)[2], // power
				       detuning, // frequency detune
				       start, // windowopen 
				       (*myAngleVec)[0], // cos_solid_angle
				       (*myAngleVec)[1], // solid_angle
				       false); // show the value
	  */
	  
	  Signal2((*myField)[2],
                  detuning,
                  (*myMuonVec)[0]*1.0e-6,
                  (*myAngleVec)[0],
                  (*myAngleVec)[1],
                  Amp13,
                  Amp24,
                  false);
	  
	  normalization++;
	}else if((method=='c'||method=='C')&&FWHM_falg==true){
	  /*
	  signal += ConventionalSignal(bp, // power
                                       detuning, // frequency detune
                                       start, // windowopen
                                       (*myAngleVec)[0], // cos_solid_angle
                                       (*myAngleVec)[1], // solid_angle                      
                                       false); // show the value
	  
	  
	  signal += OldMuoniumSignal(bp, // power                                                                                                                                              
                                     detuning, // frequency detune                                                                                                                                        
                                     (*myMuonVec)[0]*1.0e-6, // windowopen, 1.0e-6 is for convert ns to ms
                                     (*myMuonVec)[0]*1.0e-6+1, // decaytime + 0.001[ns] , 0~15*1.0e-3
                                     (*myAngleVec)[0], // cos_solid_angle                                                                                                                                 
                                     (*myAngleVec)[1], // solid_angle                                                                                                                                     
                                     false); // show the value
	  */
	  //time[i]= (*myMuonVec)[0]*1.0e-3;
	  Signal2(bp,
		  detuning,
		  (*myMuonVec)[0]*1.0e-6,
		  (*myAngleVec)[0],
		  (*myAngleVec)[1],
		  Amp13,
		  Amp24,
		  false);
	  
          normalization++;
	}else if((method=='o'||method=='O')&&FWHM_falg==false&&
		 ( start*1.0e+3<=(*myMuonVec)[0] && (*myMuonVec)[0]<=end*1.0e+3 )){ // cut positrons whose decaytime is not in the timewindow
	  /*
	  signal += OldMuoniumSignal((*myField)[2], // power
				     detuning, // frequency detune
				     start*1.0e-3, // windowopen, 1.0e-3 is for convert us to ms
				     end*1.0e-3, // windowclose 
				     (*myAngleVec)[0], // cos_solid_angle
				     (*myAngleVec)[1], // solid_angle
				     false); // show the value
	  */
	  
	  Signal2((*myField)[2],
                  detuning,
                  (*myMuonVec)[0]*1.0e-6,
                  (*myAngleVec)[0],
                  (*myAngleVec)[1],
                  Amp13,
                  Amp24,
                  false);
	  
	  normalization++;
	}else if((method=='o'||method=='O')&&FWHM_falg==true&&
		 ( start*1.0e+3<=(*myMuonVec)[0] && (*myMuonVec)[0]<=end*1.0e+3) ){ // cut positrons whose decaytime is not in the timewindow
	  /*
	  signal += OldMuoniumSignal(bp, // power
                                     detuning, // frequency detune
                                     start*1.0e-3, // windowopen, 1.0e-3 is for convert us to ms
                                     end*1.0e-3, // windowclose , end*1.0e-3
                                     (*myAngleVec)[0], // cos_solid_angle
                                     (*myAngleVec)[1], // solid_angle
                                     false); // show the value
	  */
	  
	  Signal2(bp,
                  detuning,
                  (*myMuonVec)[0]*1.0e-6,
                  (*myAngleVec)[0],
                  (*myAngleVec)[1],
                  Amp13,
                  Amp24,
                  false);
	  
          normalization++;
	}
      }
    }
   
    for(int p=0; p<gates; p++){
      Non += gRandom->Gaus(Pon, TMath::Sqrt(Pon));
      Noff += gRandom->Gaus(Poff, TMath::Sqrt(Poff));
    }
    signal = Non/Noff-1;
    error = (Non/Noff)*TMath::Sqrt(1.0/Non+1.0/Noff);

    //std::cout << "Non:" << Non << "\t" << "Noff:" << Noff << "\t" << "error" << error << std::endl;

    Sim_detected += Non;
    curve->SetPoint(w, detuning, signal);
    curve->SetPointError(w, 0, error);
  }
  /*
  TCanvas* pola = new TCanvas("pola", "pola", 900, 900);
  pola_hist->Draw();
  pola->SaveAs("../figure/pola_hist.png");

  TCanvas* c3 = new TCanvas("c3","c3", 1200, 1200);
  c3->SetGrid();
  TGraph* gr = new TGraph(l, time, p);
  gr->SetTitle("#Delta#omega=200 [kHz]"); //SetTextColor(kRed);                                                                                                                                         
  gr->SetLineColor(2);
  gr->SetLineWidth(2);
  gr->SetFillStyle(0);  
  gr->Draw("A*");
  gr->GetXaxis()->SetTitle("Time [/#muSec]");
  gr->GetYaxis()->SetTitle("P_{z}"); //#LTP_{z}#GT
  c3->SaveAs("../figure/Polarization_off_200.png");
  */
  TCanvas* c = new TCanvas("c","c", 1000, 900);
  c->SetGrid();
  gPad->SetRightMargin(0.03); // let the figure more right
  curve->GetXaxis()->SetTitle("Frequency Detuning [/kHz]");
  curve->GetYaxis()->SetTitle("Signal");
  curve->GetYaxis()->SetTitleOffset(1.4);
  curve->Draw("AP");
  
  if((fit_gamma=='y'||fit_gamma=='Y') && (method=='c'||method=='c'))
    f1 = new TF1("f1","[0]+[3]*2*[1]*[1]/(4*TMath::Pi()*TMath::Pi()*(x-[2])*(x-[2])+4*[1]*[1]+[4]*[4])", -scan_range*0.5, scan_range*0.5);
  else if((fit_gamma=='n'||fit_gamma=='N') && (method=='c'||method=='c')){
    std::string fit_no_con = "[0]+[3]*2*[1]*[1]/(4*TMath::Pi()*TMath::Pi()*(x-[2])*(x-[2])+4*[1]*[1]+"+std::to_string(gamma)+'*'+std::to_string(gamma)+')';
    f1 = new TF1("f1", fit_no_con.c_str(), -scan_range*0.5, scan_range*0.5);
  }else if((method=='o'||method=='O') && (fit_gamma=='n'||fit_gamma=='N')){
    std::string fit_no_old = "[0]+([3]*2*[1]*[1]*(TMath::Exp(-" + std::to_string(gamma) + '*';
    std::string fit_no_old2 = ")*(1-(TMath::Cos((TMath::Sqrt(4*TMath::Pi()*TMath::Pi()*(x-[2])*(x-[2])+4*[1]*[1]))*";
    std::string fit_no_old3 = ")-(TMath::Sqrt(4*TMath::Pi()*TMath::Pi()*(x-[2])*(x-[2])+4*[1]*[1]))*TMath::Sin((TMath::Sqrt(4*TMath::Pi()*TMath::Pi()*(x-[2])*(x-[2])+4*[1]*[1]))*";
    std::string fit_no_old4 = ")/"+std::to_string(gamma)+")*"+std::to_string(gamma)+'*'+std::to_string(gamma)+"/(4*TMath::Pi()*TMath::Pi()*(x-[2])*(x-[2])+4*[1]*[1]+"+std::to_string(gamma)+'*'+std::to_string(gamma)+"))-TMath::Exp(-"+std::to_string(gamma)+'*';
    std::string fit_no_old5 = ")*(1-(TMath::Cos(TMath::Sqrt(4*TMath::Pi()*TMath::Pi()*(x-[2])*(x-[2])+4*[1]*[1])*";
    std::string fit_no_old6 = ")-TMath::Sqrt(4*TMath::Pi()*TMath::Pi()*(x-[2])*(x-[2])+4*[1]*[1])*TMath::Sin(TMath::Sqrt(4*TMath::Pi()*TMath::Pi()*(x-[2])*(x-[2])+4*[1]*[1])*";
    std::string fit_no_old7 = ")/"+std::to_string(gamma)+")*"+std::to_string(gamma)+'*'+std::to_string(gamma)+"/(4*TMath::Pi()*TMath::Pi()*(x-[2])*(x-[2])+4*[1]*[1]+"+std::to_string(gamma)+'*'+std::to_string(gamma)+")))/(4*TMath::Pi()*TMath::Pi()*(x-[2])*(x-[2])+4*[1]*[1])";
    std::string fit_no_old8 = ")/(TMath::Exp(-"+std::to_string(gamma)+'*'+start_str+")-TMath::Exp(-"+std::to_string(gamma)+'*'+end_str+"))";
    f1 = new TF1("f1",(fit_no_old+start_str+fit_no_old2+start_str+fit_no_old3+start_str+fit_no_old4
		       +end_str+fit_no_old5+end_str+fit_no_old6+end_str+fit_no_old7+fit_no_old8).c_str(),-scan_range*0.5,scan_range*0.5);
    
  }else if((method=='o'||method=='O') && (fit_gamma=='y'||fit_gamma=='Y')){
    std::string fit_yes_old = "([0]+[3]*2*[1]*[1]*(TMath::Exp(-[4]*"+start_str;
    std::string fit_yes_old2 = ")*(1-(TMath::Cos((TMath::Sqrt(4*TMath::Pi()*TMath::Pi()*(x-[2])*(x-[2])+4*[1]*[1]))*"+start_str;
    std::string fit_yes_old3 = ")-(TMath::Sqrt(4*TMath::Pi()*TMath::Pi()*(x-[2])*(x-[2])+4*[1]*[1]))*TMath::Sin((TMath::Sqrt(4*TMath::Pi()*TMath::Pi()*(x-[2])*(x-[2])+4*[1]*[1]))*"+start_str;
    std::string fit_yes_old4 = ")/[4])*[4]*[4]/(4*TMath::Pi()*TMath::Pi()*(x-[2])*(x-[2])+4*[1]*[1]+[4]*[4])) - TMath::Exp(-[4]*"+end_str;
    std::string fit_yes_old5 = ")*(1-(TMath::Cos(TMath::Sqrt(4*TMath::Pi()*TMath::Pi()*(x-[2])*(x-[2])+4*[1]*[1])*"+end_str;
    std::string fit_yes_old6 = ")-TMath::Sqrt(4*TMath::Pi()*TMath::Pi()*(x-[2])*(x-[2])+4*[1]*[1])*TMath::Sin(TMath::Sqrt(4*TMath::Pi()*TMath::Pi()*(x-[2])*(x-[2])+4*[1]*[1])*"+end_str;
    std::string fit_yes_old7 = ")/[4])*[4]*[4]/(4*TMath::Pi()*TMath::Pi()*(x-[2])*(x-[2])+4*[1]*[1]+[4]*[4])))/(4*TMath::Pi()*TMath::Pi()*(x-[2])*(x-[2])+4*[1]*[1])";
    std::string fit_yes_old8 = ")/(TMath::Exp(-[4]*"+start_str+")-TMath::Exp(-[4]*"+end_str+"))";
    f1 = new TF1("f1", (fit_yes_old+fit_yes_old2+fit_yes_old3+fit_yes_old4+fit_yes_old5+fit_yes_old6+fit_yes_old7+fit_yes_old8).c_str(), -scan_range*0.5, scan_range*0.5);
  }
   
  f1->FixParameter(0, 0); // offset
  //f1->SetParameter(0, 0); 
  if(FWHM_falg) f1->SetParameter(1, bp); // microwave power,bp
  else  f1->SetParameter(1, power_mean);
  f1->SetParameter(2, 0); // center
  if((method=='c'||method=='C')||((method=='o'||method=='O')&&(start>=2.)))
    f1->SetParameter(3, 2); // scaling
  else if((method=='o'||method=='O')&&(start<=1.))
    f1->SetParameter(3, 1); // scaling
  
  if(fit_gamma=='y'||fit_gamma=='Y'){
    f1->SetParameter(4, gamma); // gamma
    f1->SetParNames("Offset", "b [/kHz]", "Center", "Scaling", "#gamma [/kHz]");
  }else if(fit_gamma=='n'||fit_gamma=='N')
    f1->SetParNames("Offset", "b [/kHz]", "Center", "Scaling");

  
  //curve->Fit("f1","EM", "", -scan_range*0.5, scan_range*0.5);
  TFitResultPtr FitResult = curve->Fit("f1","S", "", -scan_range*0.5, scan_range*0.5);
  TMatrixDSym GetCorrelation_fit = FitResult->GetCorrelationMatrix();
  TMatrixDSym Covariance_fit = FitResult->GetCovarianceMatrix();

  // Serach FWHM
  Double_t fit_power = f1->GetParameter("b [/kHz]");
  Double_t fit_offset = f1->GetParameter("Offset");
  Double_t fit_Scaling = f1->GetParameter("Scaling");
  fit_center = f1->GetParameter(2);
  fit_center_error = f1->GetParError(2);
  
  Double_t fit_gamma_para;
  if(fit_gamma==true) fit_gamma_para = f1->GetParameter("#gamma [/kHz]");
  center_error = f1->GetParError(2);
  
  if(FWHM_falg) std::cout << "Input Microwave power:" << bp << "[/kHz]" << std::endl;
  else std::cout << "Fitting Microwave power:" << fit_power << "[/kHz]" << std::endl;
  
  std::cout << std::string(60, '*') << std::endl;
  Sim_Height = f1->Eval(fit_center)/fit_Scaling;
  Double_t HM_left = f1->GetX(0.5*f1->Eval(fit_center), -200, 0.);
  Double_t HM_right = f1->GetX(0.5*f1->Eval(fit_center), 0., 200);
  Sim_FWHM = HM_right - HM_left;

  std::cout << "Simulation FWHM:" << Sim_FWHM << "[/kHz]" << "\n"
	    << "Simulation Signal Height:" << Sim_Height << std::endl;

  std::cout << "FOM:" << Sim_detected*pow(Sim_Height+1,2.)/pow(Sim_FWHM,2.) << std::endl;

  if(method=='c'||method=='C'){
    The_FWHM = sqrt(4*pow(bp,2.)+pow(gamma,2.))/pi;
    The_Height = 2*pow(bp, 2.)/(4*pow(bp,2.)+pow(gamma,2.));
    std::cout << "Theoritical FWHM:" << The_FWHM << "[/kHz]" << "\n"
	      << "Theoritiacl Singal Height:" << The_Height << std::endl;
  }else if(method=='o'||method=='O'){
    std::string HM_power = std::to_string(bp)+'*'+std::to_string(bp);
    //std::string HM_formula = HM_power+"*2*(1-cos(TMath::Sqrt(4*TMath::Pi()*TMath::Pi()*x*x+4*"+HM_power+")*"+start_str+"))/(4*TMath::Pi()*TMath::Pi()*x*x+4*"+HM_power+')';
    std::string HM_formula =
      HM_power+"*2/(((2*TMath::Pi()*x)**2)+4*"+HM_power+")*(TMath::Exp(-"+std::to_string(gamma)+'*'+start_str+")*(1-(TMath::Cos((TMath::Sqrt(((2*TMath::Pi()*x)**2)+4*"+HM_power+"))*"+start_str+")-(TMath::Sqrt(((2*TMath::Pi()*x)**2)+4*"+HM_power+"))/"+std::to_string(gamma)+"*TMath::Sin((TMath::Sqrt(((2*TMath::Pi()*x)**2)+4*"+HM_power+"))*"+start_str+"))*"+std::to_string(gamma)+'*'+std::to_string(gamma)+"/((((2*TMath::Pi()*x)**2)+4*"+HM_power+")+"+std::to_string(gamma)+'*'+std::to_string(gamma)+"))-TMath::Exp(-"+std::to_string(gamma)+'*'+end_str+")*(1-(TMath::Cos((TMath::Sqrt(((2*TMath::Pi()*x)**2)+4*"+HM_power+"))*"+end_str+")-(TMath::Sqrt(((2*TMath::Pi()*x)**2)+4*"+HM_power+"))/"+std::to_string(gamma)+"*TMath::Sin((TMath::Sqrt(((2*TMath::Pi()*x)**2)+4*"+HM_power+"))*"+end_str+"))*"+std::to_string(gamma)+'*'+std::to_string(gamma)+"/((((2*TMath::Pi()*x)**2)+4*"+HM_power+")+"+std::to_string(gamma)+'*'+std::to_string(gamma)+")))/(TMath::Exp(-"+std::to_string(gamma)+'*'+start_str+")-TMath::Exp(-"+std::to_string(gamma)+'*'+end_str+"))";
    TF1 *fa1 = new TF1("fa1",
                       HM_formula.c_str(),
                       -scan_range, scan_range);
    The_Height = fa1->Eval(0.);
    Double_t HM_left_the = fa1->GetX(0.5*The_Height, -200, 0.);
    Double_t HM_right_the = fa1->GetX(0.5*The_Height, 0., 200);
    The_FWHM = HM_right_the - HM_left_the;
    std::cout << "Theoritical FWHM:" << The_FWHM << "[/kHz]" << "\n"
	      << "Theoritiacl Singal Height:" << The_Height << std::endl;
  }
  if(method=='c'||method=='C') curve->SetTitle("Conventional");
  else if(method=='o'||method=='O') curve->SetTitle(("t1="+std::to_string(start)+" [#muSec]").c_str());

  //if((method=='c'||method=='C')&&FWHM_falg==false) c->SaveAs(("../figure/"+run_num+":"+tree_TMmode+":"+tree_Pressure+":"+"conventional"+".png").c_str());
  if(method=='c'||method=='C') c->SaveAs(("../figure/"+run_num+":"+tree_TMmode+":"+tree_Pressure+":"+"conventional"+".png").c_str());
  //else if((method=='o'||method=='O')&&FWHM_falg==false) c->SaveAs(("../figure/"+run_num+":"+tree_TMmode+":"+tree_Pressure+":"+"oldmuonium"+".png").c_str());
  else if(method=='o'||method=='O') c->SaveAs(("../figure/"+run_num+":"+tree_TMmode+":"+tree_Pressure+":"+start_str+".png").c_str());
  
  delete gRandom;
  delete f1;
  delete curve;
  delete c;
}
#endif
