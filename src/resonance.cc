////////////////////////////////////////////////////////
//   High Field Simulation for MuSUEUM Collaboration
//
//       Author : Hideharu Yamauchi 2021/10/15
////////////////////////////////////////////////////////
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "simulator.cc"
#include "TH1D.h"

int main(int argc, const char** argv){
  SIMULATOR* sim = new SIMULATOR(argv[1]);

  // true is for const RF field, false is for non-const RF field
  //sim->CalculateSignal(true, 230., 0., 12.);
  //sim->Vis_StateAmp(0.);

  //for(int i=0; i<12; i++) sim->CalculateSignal(true, 230., double(i), double(i+1));

  int power_min = 100;
  int power_max = 240;
  int power;

  /*
  for(int i=0; i<1+(power_max-power_min)/10; i++){                                                                                                                                                       
      power = power_min + 10*i;                                                                                                                                
      sim->CalculateSignal(true, power, 9., 10.);
  }
  */
  
  std::ofstream ofs("../figure/FWHM_conv.dat", std::ios::out);
  std::ofstream ofs2("../figure/FWHM_conv_hight.dat", std::ios::out);
  //std::ofstream ofs("../figure/FWHM_old.dat", std::ios::out);
  //std::ofstream ofs2("../figure/FWHM_old_hight.dat", std::ios::out);
  
  
  if(ofs&&ofs2){
    std::cout << "Successful to open file." << std::endl;
    for(int i=0; i<1+(power_max-power_min)/10; i++){
      power = power_min + 10*i;
      sim->CalculateSignal(true, power, 9., 10.);
      ofs << sim->Sim_FWHM << "\t"
	  << sim->The_FWHM << std::endl;
      ofs2 << sim->Sim_Height << "\t"
	   << sim->The_Height << std::endl;
    }
    ofs.close();
    ofs2.close();
  }
  /*
  std::ofstream ofs("../figure/t1.dat", std::ios::out);
  if(ofs){
    std::cout << "Successful to open file." << std::endl;
    for(int t=3; t<13; t++){
      sim->CalculateSignal(true, 125., t, t+1);
      ofs << t << "\t"
	  << sim->center_error << std::endl;
    }
    ofs.close();
  }
  Double_t error;
  for(int i=0; i<10; i++){
    sim->CalculateSignal(true, 125., 9., 10.);
    error += sim->center_error;
  }
  std::cout << "average of center error by conventional:" << error/10 << std::endl;
  */

  /*
  std::ofstream ofs("../figure/FWHM_windowopen.dat", std::ios::out);
  if(ofs){ 
    std::cout << "Successful to open file." << std::endl;                                                                                                                                                 
    for(int t=2; t<12; t++){
      sim->CalculateSignal(true, 230., t, t+1);
      ofs << t << "\t"
          << sim->Sim_FWHM << "\t"
	  << sim->The_FWHM << std::endl;
    }                                                                                                                                              
    ofs.close();
  }
  */
  /*
  std::ofstream ofs("../figure/merit_230_tm110_2.dat", std::ios::out);
  Double_t det;
  if(ofs){
    std::cout << "Successful to open file." << std::endl;                                                                                                                                                 
    for(int t=2; t<12; t++){
      sim->CalculateSignal(true, 230., t, t+1);
      ofs << t << "\t"
          << sim->Sim_FWHM << "\t"
	  << sim->Sim_Height << "\t"
	  << sim->Sim_detected << std::endl;
    }
    ofs.close();
  }
  */
  /*
  std::ofstream ofs("../figure/t1_height.dat", std::ios::out);
  if(ofs){
    std::cout << "Successful to open file." << std::endl;                                                                                                                                                 
    for(int t=0; t<12; t++){
      sim->CalculateSignal(true, 230., t, t+1);   
      ofs << t << "\t"   
          << sim->Sim_Height << "\t"                                                                                                                                                                      
          << sim->The_Height << std::endl;
    }
    ofs.close();    
  }
  */

  // true is for const RF power, false is for non-const RF power
  /*    
  TCanvas* cc = new TCanvas("cc","cc", 1200, 1200);
  gPad->SetRightMargin(0.03);
  TGraphErrors* ge = new TGraphErrors();
  for(int tt=1; tt<12; tt++){
    sim->CalculateSignal(true, 230, tt, tt+1);
    ge->SetPoint(tt, double(tt), sim->fit_center); //Int_t i, Double_t x, Double_y y
    ge->SetPointError(tt, 0, sim->fit_center_error); //Int_t i, Double_t ex, Double_y ey
  }
  ge->SetTitle("TM_{210}, b=230 [kHz]");
  ge->GetXaxis()->SetTitle("t1 [/#muSec]");
  ge->GetYaxis()->SetTitle("Frequency Center [kHz]");
  ge->GetYaxis()->SetTitleOffset(1.2);
  ge->SetMarkerColor(4);
  ge->SetMarkerStyle(21);
  ge->Draw("AP");
  TF1 *f2 = new TF1("f2","[0]",0, 12);
  f2->SetParameter(0, 0.);
  f2->SetParNames("Center [kHz]");
  ge->Fit("f2");
  cc->SaveAs("../figure/virtual_trial_tm210_old.png");
  */
  
  /*
  TCanvas* cc = new TCanvas("cc","cc", 1200, 1200);
  gPad->SetRightMargin(0.03);
  TGraphErrors* ge2 = new TGraphErrors();
  TF1* f2 = new TF1("f2","[0]", 0, 12);
  f2->SetParNames("Center [kHz]");
  f2->SetParameter(0, 0.);
  for(int i=0; i<100; i++){
    TGraphErrors* ge = new TGraphErrors();
    for(int tt=0; tt<12; tt++){
      sim->CalculateSignal(false, 230, tt, tt+1);    
      ge->SetPoint(tt, double(tt), sim->fit_center); //Int_t i, Double_t x, Double_y y
      ge->SetPointError(tt, 0, sim->fit_center_error); //Int_t i, Double_t ex, Double_y ey
    }
    ge->Fit("f2");
    ge2->SetPoint(i, double(i+1), f2->GetParameter(0)); //Int_t i, Double_t x, Double_y y
    ge2->SetPointError(i, 0, f2->GetParError(0)); //Int_t i, Double_t ex, Double_y ey
    delete ge;
  }
  TF1* f2a = new TF1("f2a","[0]",0, 100);
  f2a->SetParNames("Center [kHz]");
  f2a->SetParameter(0, 0.);
  ge2->SetTitle("TM_{210}, b=230 [kHz]");
  ge2->GetXaxis()->SetTitle("Number of Virtual Trial");
  ge2->GetYaxis()->SetTitle("Frequency Center [kHz]");
  ge2->GetYaxis()->SetTitleOffset(1.4);
  ge2->SetMarkerColor(4);
  ge2->SetMarkerStyle(21);
  ge2->Draw("AP");
  ge2->Fit("f2a");
  cc->SaveAs("../figure/virtual_trial_tm210_dist_old2.png"); 
  */
  /*
  TCanvas* cc = new TCanvas("cc","cc", 1200, 1200);
  cc->Divide(1, 2);
  TGraphErrors* ge = new TGraphErrors();
  TH1D* hist = new TH1D("hist","",100 ,-0.1 ,0.1);
  double start = 0.;
  for(int i=0; i<100; i++){
    sim->CalculateSignal(false, 230, start, start+1);
    ge->SetPoint(i, double(i+1), sim->fit_center); //Int_t i, Double_t x, Double_y y
    ge->SetPointError(i, 0, sim->fit_center_error); //Int_t i, Double_t ex, Double_y ey
    hist->Fill(sim->fit_center);
  }
  
  cc->cd(1);
  gPad->SetGridy(1); 
  ge->GetXaxis()->SetTitle("Number of Virtual Trial");
  ge->GetYaxis()->SetTitle("Frequency Center [kHz]");
  ge->SetMarkerColor(4);
  ge->SetMarkerStyle(21);
  ge->Draw("AP");

  cc->cd(2);
  hist->GetXaxis()->SetTitle("Frequency Center [kHz]");
  hist->GetYaxis()->SetTitle("");
  hist->Draw();

  TF1 *f2 = new TF1("f2","[0]",0, 100);
  f2->SetParameter(0, 0.);
  f2->SetParNames("Center [kHz]");

  ge->Fit("f2");
  cc->SaveAs("../figure/virtual_trial_con_power_tm210.png");
  */

  delete sim;
  return 0;
}
