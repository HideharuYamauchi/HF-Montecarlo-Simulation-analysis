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
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TRandom3.h"

//ROOT::Fit::FitResult* FitResult = new ROOT::Fit::FitResult();

SIMULATOR::SIMULATOR(const char* rootfile, Int_t run_time=20)
  : myMuonVec(0),myMuonDispersion(0),myPositronVec(0),myPositronDispersion(0),myField(0),myAmp(0),myAngleVec(0),
    Non(0), scan_range(800), scan_step(10), scan_points(81), signal(0), power_mean(0), threshold(35.), solid_angle_mean(0), cos_solid_angle_mean(0)
{
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
  
  std::string myTreeName = myTree->GetName();
  myTreeTitle = myTree->GetTitle();
  tree_TMmode = myTreeTitle.substr(myTreeTitle.find("TM"), 5);
  tree_Pressure = myTreeTitle.substr(myTreeTitle.find("atmosphere")-4, 14);
  tree_Temperature = myTreeTitle.substr(myTreeTitle.find("kelvin")-4, 10);

  int l=0;
  for(int i=0; i<entries; i++){
    myTree->GetEntry(i);
    if(threshold<=(*myPositronDispersion)[0]*1.0e-3){
      power_mean += (*myField)[2];
      solid_angle_mean += (*myAngleVec)[0];
      cos_solid_angle_mean += (*myAngleVec)[1];
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
  
  do{
    std::cout << "Choose the plot method:" << std::endl;
    std::cin >> plot_method;
  }while(plot_method!=1 && plot_method!=2);
  
  do{ // use gamma as free parameter or not
    std::cout << "Fitting using gamma? [y/Y]or[n/N]:";
    std::cin >> fit_gamma;
  }while(fit_gamma!='y'&&fit_gamma!='Y'&&fit_gamma!='n'&&fit_gamma!='N');

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
  
  Pon += 0.5*pow(y,2.)*(A[1]*solid_angle+polarization*A[0]*cos_solid_angle*( 1-L ))/pi;
  Poff += 0.5*pow(y,2.)*(A[1]*solid_angle+polarization*A[0]*cos_solid_angle)/pi;

  Double_t sig = A[0]*(cos_solid_angle)*(polarization*L)/((A[0]*cos_solid_angle*polarization+A[1]*solid_angle)*(0-std::exp(-1*gamma*windowopen)));
  
  if(flag) std::cout << "L:" << L << "\n"
		     << "Pon:" << Pon << "\t"
		     << "Poff:" << Poff << "\n"
		     << "signal(Pon/Poff-1):" << Pon/Poff-1 << "\t"
		     << "signal(Non/Noff-1):" << sig << "\n"
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
 
  Pon += 0.5*pow(y,2.)*(A[1]*(std::exp(-1*gamma*windowopen)-std::exp(-1*gamma*windowclose))*solid_angle
			+polarization*A[0]*cos_solid_angle*((std::exp(-1*gamma*windowopen)-std::exp(-1*gamma*windowclose))-L))/pi;
  Poff += 0.5*pow(y,2.)*(A[1]*(std::exp(-1*gamma*windowopen)-std::exp(-1*gamma*windowclose))*solid_angle
			 +polarization*A[0]*(std::exp(-1*gamma*windowopen)-std::exp(-1*gamma*windowclose))*cos_solid_angle)/pi;

  Double_t sig = A[0]*(cos_solid_angle)*(polarization*L)/((A[0]*(cos_solid_angle)*polarization+A[1]*(solid_angle))*(std::exp(-1*gamma*windowclose)-std::exp(-1*gamma*windowopen)));
  
  if(flag) std::cout << "L:" << L << "\n"
		     << "Pon:" << Pon << "\t"
                     << "Poff:" << Poff << "\n"
                     << "signal(Pon/Poff-1):" << Pon/Poff-1 << "\t"
		     << "signal(Non/Noff-1):" << sig << "\n"
		     << "energy:" << y << std::endl;
  return sig;
}

void SIMULATOR::CalculateSignal(bool FWHM_falg=false, Double_t bp=0., Double_t start=6.5, Double_t end=7.5){
  //
  // start = 3.12, 6.92 second                                                                                                                                                                     
  // end = 4.07, 7.87 second
  //
  
  TCanvas* c = new TCanvas("c","c", 1000, 900);
  c->SetGrid();
  gPad->SetRightMargin(0.03); // let the figure more right
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
	    << std::string(60, '*') << "\n"
	    << "RUN START: " << ctime(&t) << std::endl;

  //std::ofstream ofs_con("../figure/reso_conv.dat", std::ios::app);
  //std::ofstream ofs_old("../figure/reso_old.dat", std::ios::app);

  Sim_detected = 0;
  for(int w=0; w<scan_points; w++){
    Non = 0;
    Noff = 0;
    Totalon = 0;
    Totaloff = 0;
    Pon = 0;
    Poff = 0;
    signal = 0;
    normalization = 0;
    detuning = -1*scan_range*0.5 + scan_step*w;
    if(FWHM_falg==false) std::cout << "START detuning " << detuning << "[/kHz]..." << std::endl;
    for(int i=0; i<entries; i++){
      myTree->GetEntry(i);
      if(threshold<=(*myPositronDispersion)[0]*1.0e-3/* && // cut positrons below threshold energy, 35 MeV for liu
							( start*1.0e+3<=(*myMuonVec)[0] && (*myMuonVec)[0]<=end*1.0e+3 )*/){ // cut positrons whose decaytime is not in the timewindow
	y = (*myPositronDispersion)[0]*1.0e-3/positron_max_energy;
	A[0] = 2*y-1;
	A[1] = 3-2*y;
	if((method=='c'||method=='C')&&FWHM_falg==false){
	  signal += ConventionalSignal((*myField)[2], // power
				       detuning, // frequency detune
				       start, // windowopen 
				       (*myAngleVec)[0], // cos_solid_angle
				       (*myAngleVec)[1], // solid_angle
				       false); // show the value
	  normalization++;
	}else if((method=='c'||method=='C')&&FWHM_falg==true){
	  signal += ConventionalSignal(bp, // power
                                       detuning, // frequency detune
                                       start, // windowopen
                                       (*myAngleVec)[0], // cos_solid_angle
                                       (*myAngleVec)[1], // solid_angle                                                                             
                                       false); // show the value
          normalization++;
	}else if((method=='o'||method=='O')&&FWHM_falg==false){
	  signal += OldMuoniumSignal((*myField)[2], // power
				     detuning, // frequency detune
				     start*1.0e-3, // windowopen, 1.0e-3 is for convert us to ms
				     end*1.0e-3, // windowclose 
				     (*myAngleVec)[0], // cos_solid_angle
				     (*myAngleVec)[1], // solid_angle
				     false); // show the value
	  normalization++;
	}else if((method=='o'||method=='O')&&FWHM_falg==true){
	  signal += OldMuoniumSignal(bp, // power
                                     detuning, // frequency detune
                                     start*1.0e-3, // windowopen, 1.0e-3 is for convert us to ms
                                     end*1.0e-3, // windowclose
                                     (*myAngleVec)[0], // cos_solid_angle
                                     (*myAngleVec)[1], // solid_angle
                                     false); // show the value
          normalization++;
	}
      }
    }

    //std::cout << "total e+:" << normalization++ << std::endl;
    
    if(plot_method==1){
      signal = signal/normalization;
      y = threshold/positron_max_energy;
      Noff = normalization*0.25*( (1-(2*pow(y,3.)-pow(y,4.)))*solid_angle_mean + polarization*(1-(3*pow(y,4.)-2*pow(y,3.)))*cos_solid_angle_mean/3)/pi;
      Non = Noff*(1+signal);
      for(int p=0; p<gates; p++) Totalon += gRandom->Gaus(Non, TMath::Sqrt(Non));
      Totaloff = gRandom->Gaus(Noff*gates, TMath::Sqrt(Noff*gates));
      signal = Totalon/Totaloff-1;
      error = (Totalon/Totaloff)*TMath::Sqrt(1.0/Totalon+1.0/Totaloff);
    }
    
    if(plot_method==2){
      for(int p=0; p<gates; p++){
	Non += gRandom->Gaus(Pon, TMath::Sqrt(Pon));
	Noff += gRandom->Gaus(Poff, TMath::Sqrt(Poff));
      }
      signal = Non/Noff-1;
      error = (Non/Noff)*TMath::Sqrt(1.0/Non+1.0/Noff);
    }
    /*
    if(method=='c'||method=='C')
      ofs_con << detuning << "\t" << Non  << "\t" << Noff << "\t" << error << "\n";
    else if(method=='o'||method=='O')
      ofs_old << detuning << "\t" << Non  << "\t" << Noff << "\t" << error << "\n";
    */

    Sim_detected += Non;
    curve->SetPoint(w, detuning, signal);
    curve->SetPointError(w, 0, error);
  }
  std::cout << "Detected e+:" << Sim_detected << std::endl;

  // calculate the detected e+ of theory
  /*
  Double_t interval = std::exp(-gamma*end*1.0e-3)-std::exp(-gamma*start*1.0e-3);
  y = threshold/positron_max_energy;
  The_detected = (0.25*normalization/pi)*( interval*(1-(2*pow(y,3.)-pow(y,4.)))*solid_angle_mean + 0.5*polarization*(1-(3*pow(y,4.)-2*pow(y,3.)))*()*cos_solid_angle_mean/3 );
  */

  
  //ofs_con.close();
  //ofs_old.close();
  t = time(NULL);
  std::cout << "RUN FINISH: " << ctime(&t) << std::string(60, '*') << std::endl;
  
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
  if(FWHM_falg) f1->SetParameter(1, bp); // microwave power
  else  f1->SetParameter(1, power_mean);
  f1->SetParameter(2, 0); // center                                                                                                                                                               
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

  
  std::cout << "Correlation Matrix:";
  GetCorrelation_fit.Print("v");

  
  std::cout << "Covariance Matrix:";
  Covariance_fit.Print("v");

  // Serach FWHM
  Double_t fit_power = f1->GetParameter("b [/kHz]");
  Double_t fit_offset = f1->GetParameter("Offset");
  Double_t fit_Scaling = f1->GetParameter("Scaling");
  Double_t fit_center = f1->GetParameter("Center");
  Double_t fit_gamma_para;
  if(fit_gamma==true) fit_gamma_para = f1->GetParameter("#gamma [/kHz]");
  center_error = f1->GetParError(2);
  
  if(FWHM_falg) std::cout << "Input Microwave power:" << bp << "[/kHz]" << std::endl;
  else std::cout << "Fitting Microwave power:" << fit_power << "[/kHz]" << std::endl;
  
  std::cout << std::string(60, '*') << std::endl;
  Sim_Height = f1->Eval(fit_center)/fit_Scaling;
  Double_t HM_left = f1->GetX(0.5*f1->Eval(fit_center), -400, 0.);
  Double_t HM_right = f1->GetX(0.5*f1->Eval(fit_center), 0., 400);
  Sim_FWHM = HM_right - HM_left;

  std::cout << "Simulation FWHM:" << Sim_FWHM << "[/kHz]" << "\n"
	    << "Simulation Signal Height:" << Sim_Height << std::endl;

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
    Double_t HM_left_the = fa1->GetX(0.5*The_Height, -400, 0.);
    Double_t HM_right_the = fa1->GetX(0.5*The_Height, 0., 400);
    The_FWHM = HM_right_the - HM_left_the;
    std::cout << "Theoritical FWHM:" << The_FWHM << "[/kHz]" << "\n"
	      << "Theoritiacl Singal Height:" << The_Height << std::endl;
  }
  
  if((method=='c'||method=='C')&&FWHM_falg==false) c->SaveAs(("../figure/"+run_num+":"+tree_TMmode+":"+tree_Pressure+":"+"conventional"+".png").c_str());
  //if(method=='c'||method=='C') c->SaveAs(("../figure/"+run_num+":"+tree_TMmode+":"+tree_Pressure+":"+"conventional"+".png").c_str());
  //else if((method=='o'||method=='O')&&FWHM_falg==false) c->SaveAs(("../figure/"+run_num+":"+tree_TMmode+":"+tree_Pressure+":"+"oldmuonium"+".png").c_str());
  else if(method=='o'||method=='O') c->SaveAs(("../figure/"+run_num+":"+tree_TMmode+":"+tree_Pressure+":"+"oldmuonium"+".png").c_str());

  delete gRandom;
  delete f1;
  delete curve;
  delete c;
}
#endif
