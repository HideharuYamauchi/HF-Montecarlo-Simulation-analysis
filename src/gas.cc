////////////////////////////////////////////////////////    
//   High Field Simulation for MuSUEUM Collaboration      
//                              
//       Author : Hideharu Yamauchi 2022/01/04
////////////////////////////////////////////////////////
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "simulator.cc"

int main(){
  int mode;
  
  do{
    std::cout << "TM110[1] or TM210[2]:" << std::endl;
    std::cin >> mode;
  }while((mode!=1)&&(mode!=2));

  
  double Pressure = 1.; // atm
  double Pressure2 = 0.5; // atm
  double Pressure3 = 0.3; //atm
  Pressure = Pressure*101325; // change atm to Pascal
  Pressure2 = Pressure2*101325; // change atm to Pascal
  Pressure3 = Pressure3*101325; // change atm to Pascal
  const double temperature = 293.15; // kelvin
  const double R = 8.31446261815324; // J K^-1 mol^-1
  double Pressure_cal;
  double Pressure_cal2;
  double Pressure_cal3;
  double mol_density = 2.39953279e-02; // molar volume(1 atm, 293.15 kelvin), m^3 mol^-1
  double mol_density2 = 4.80541163e-02; // molar volume(0.5 atm, 293.15 kelvin), m^3 mol^-1
  double mol_density3 = 8.01279653e-02; // molar volume(0.5 atm, 293.15 kelvin), m^3 mol^-1
  double av = 0.2325; // Pa m^6 mol^-2
  double bv = 3.96*1.0e-5; // m^3 mol^-1

  // calculate the pressure at 0 273.15 Kelvin
  Pressure_cal = ((273.15-temperature)*av/(mol_density*mol_density)+273.15*Pressure)/temperature;
  Pressure_cal = Pressure_cal/101325; // change Pa to atm
  std::cout << "calibration pressure is " << Pressure_cal << " [atm]" << std::endl; // 0.931504 Pa

  Pressure_cal2 = ((273.15-temperature)*av/(mol_density2*mol_density2)+273.15*Pressure2)/temperature;
  Pressure_cal2 = Pressure_cal2/101325; // change Pa to atm
  std::cout << "calibration pressure is " << Pressure_cal2 << " [atm]" << std::endl; // Pa

  Pressure_cal3 = ((273.15-temperature)*av/(mol_density3*mol_density3)+273.15*Pressure3)/temperature;
  Pressure_cal3 = Pressure_cal3/101325; // change Pa to atm                                                                                                                                               
  std::cout << "calibration pressure is " << Pressure_cal3 << " [atm]" << std::endl; // Pa
  
  gStyle->SetOptFit();
  
  TCanvas* c = new TCanvas("c","c", 1200, 1200);
  gStyle->SetOptFit(1111);
  gStyle->SetTitleXOffset(1.5);
  gStyle->SetTitleYOffset(0.1);
  gStyle->SetPadRightMargin(-0.2);

  double a = -8.2*1.0e-6; // atm^-1                                                                                                                                                                   
  double b = -6.6*1.0e-9; // atm^-2
  
  double v12 = 1897431.594; // kHz
  double v34 = 2565871.274; // kHz
  double v = 4463302.776; // LAMPF
  double v0, v1, v2;
  if(mode==1){
    v0 = v12*(1+a*Pressure_cal+b*Pressure_cal*Pressure_cal);
    v1 = v12*(1+a*Pressure_cal2+b*Pressure_cal2*Pressure_cal2);
    v2 = v12*(1+a*Pressure_cal3+b*Pressure_cal3*Pressure_cal3);
    v0 = v0 - v12;
    v1 = v1 - v12;
    v2 = v2 - v12;
  }
  else if(mode==2){
    v0 = v34*(1+a*Pressure_cal+b*Pressure_cal*Pressure_cal);
    v1 = v34*(1+a*Pressure_cal2+b*Pressure_cal2*Pressure_cal2);
    v2 = v34*(1+a*Pressure_cal3+b*Pressure_cal3*Pressure_cal3);
    v0 = v0 - v34;
    v1 = v1 - v34;
    v2 = v2 - v34;
  }
  /*
  std::cout << "MuHFS at 1 atm: " << v0 << " [kHz]" << std::endl;
  std::cout << "MuHFS Simulated at 1 atm: " << v1 << " [kHz]" << std::endl;
  std::cout << "MuHFS at 0.5 atm: " << v1 << " [kHz]" << std::endl;
  std::cout << "MuHFS Simulated at 1 atm: " << v2 << " [kHz]" << std::endl;
  std::cout << "MuHFS at 0.3 atm: " << v2 << " [kHz]" << std::endl;
  */
  double x[3] = {0}, y[3] = {0}, ex[3] = {0,0,0}, ey[3] = {0,0,0};
  if(mode==1){
    ey[0] = 0.0036;
    ey[1] = 0.0036;
    ey[2] = 0.01;
  }
  else if(mode==2){
    ey[0] = 0.0049;
    ey[1] = 0.0049;
    ey[2] = 0.0144;
  }
  double shift = 0. ; // Statistic error only
  //double shift = 6.71*1.0e-4; // Kr density fluctuations, atm
  //double shift = 3*1.0e-5;// Calibration of Kr density
  //double shift = 0.2; // Drift of Kr density calibration
  double b_unce = 1.3*1.0e-9; // Quadratic pressure shift
  
  x[0] = Pressure_cal;
  //ex[0] = shift; // Kr density fluctuations
  //ex[0] = Pressure_cal*shift; // Calibration of Kr density
  //ex[0] = (shift*R/(mol_density-bv))/101325; // Drift of Kr density calibration
  y[0] = v0; 
  //ey[0] += pow(v*(a+2*b*Pressure_cal)*ex[0], 2.);
  ey[0] += pow(v*Pressure_cal*Pressure_cal*b_unce, 2.); // Quadratic_pressure_shift
  ey[0] = sqrt(ey[0]);
  std::cout << "x at 1 atm : " << x[0] << " [atm]"<< std::endl;
  std::cout << "y at 1 atm : " << y[0] << " [kHz]"<< std::endl;
  std::cout << "error of x at 1 atm : " << ex[0] << " [atm]"<< std::endl;
  std::cout << "error of y at 1 atm : " << ey[0] << " [kHz]"<< std::endl;

  x[1] = Pressure_cal2;
  //ex[1] = shift; // Kr density fluctuations
  //ex[1] = Pressure_cal2*shift; // Calibration of Kr density
  //ex[1] = (shift*R/(mol_density2-bv))/101325; // Drift of Kr density calibration
  y[1] = v1; 
  //ey[1] += pow(v*(a+2*b*Pressure_cal2)*ex[1], 2.);
  ey[1] += v*Pressure_cal2*Pressure_cal2*b_unce; // Quadratic_pressure_shift
  ey[1] = sqrt(ey[1]);
  std::cout << "x at 0.5 atm : " << x[1] << " [atm]"<< std::endl;
  std::cout << "y at 0.5 atm : " << y[1] << " [kHz]"<< std::endl;
  std::cout << "error of x at 0.5 atm : " << ex[1] << " [atm]"<< std::endl;
  std::cout << "error of y at 0.5 atm : " << ey[1] << " [kHz]"<< std::endl;

  x[2] = Pressure_cal3;
  //ex[2] = shift; // Kr density fluctuations
  //ex[2] = Pressure_cal3*shift; // Calibration of Kr density
  //ex[2] = (shift*R/(mol_density3-bv))/101325; // Drift of Kr density calibration
  y[2] = v2;
  //ey[2] += pow(v*(a+2*b*Pressure_cal3)*ex[2], 2.);
  ey[2] += v*Pressure_cal3*Pressure_cal3*b_unce; // Quadratic_pressure_shift
  ey[2] = sqrt(ey[2]);
  std::cout << "x at 0.3 atm : " << x[2] << " [atm]"<< std::endl;
  std::cout << "y at 0.3 atm : " << y[2] << " [kHz]"<< std::endl;
  std::cout << "error of x at 0.3 atm : " << ex[2] << " [atm]"<< std::endl;
  std::cout << "error of y at 0.3 atm : " << ey[2] << " [kHz]"<< std::endl;

  double b_upper = b + b_unce;
  double b_lower = b - b_unce;

  std::string fo;
  if(mode==1){
    //fo = "[0]*(1+[1]*x+"+std::to_string(b)+"*x*x)-1897431";
    fo = "[0]*(1+[1]*x+[2]*x*x)-1897431";
  }
  else if(mode==2){
    //fo = "[0]*(1+[1]*x+"+std::to_string(b)+"*x*x)-2565871";
    fo = "[0]*(1+[1]*x+[2]*x*x)-2565871";
  }
  TF1* f2 = new TF1("f2", fo.c_str(), 0., 1.);
  //f2->SetParNames("#nu_{0} [kHz]", "a [atm^{-1}]"); 
  f2->SetParNames("#nu_{0} [kHz]", "a [atm^{-1}]", "b [atm^{-2}]"); // Quadratic_pressure_shift
    
  std::string fo_upper = "[0]*(1+[1]*x+"+std::to_string(b_upper)+"*x*x)-4463302";
  TF1* f2_upper = new TF1("f2_upper", fo_upper.c_str(), 0., 1.);
  f2_upper->SetParNames("#nu_{0} [kHz]", "a [atm^{-1}]");

  std::string fo_lower = "[0]*(1+[1]*x+"+std::to_string(b_lower)+"*x*x)-4463302";
  TF1* f2_lower = new TF1("f2_lower", fo_lower.c_str(), 0., 1.);
  f2_lower->SetParNames("#nu_{0} [kHz]", "a [atm^{-1}]");
  
  
  //TF1* f2 = new TF1("f2","[0]*(1+[1]*x)-4463302", 0., 1.);
  //f2->SetParNames("#nu_{0}", "a [atm^{-1}]");
  if(mode==1)
    f2->SetParameter(0, v12);
  else if(mode==2)
    f2->SetParameter(0, v34);
  f2->SetParameter(1, a);//+0.61*1.0e-7); // Calibration of Kr density
  //f2->SetParameter(1, a); 
  //f2->SetParameter(2, b+0.6*1.0e-10); // Calibration of Kr density
  f2->SetParameter(2, b); // Quadratic_pressure_shift
  /*
  f2_upper->SetParameter(0, v);
  f2_upper->SetParameter(1, a);
  f2_lower->SetParameter(0, v);
  f2_lower->SetParameter(1, a);
  */
  TGraphErrors* gr = new TGraphErrors(3,x,y,ex,ey);
  if(mode==1){
    //gr->SetTitle("Drift of Kr density calibration(TM110)");
    gr->SetTitle("Quadratic pressure shift(TM110)");
    //gr->SetTitle("Calibration of Kr density(TM110)");
    //gr->SetTitle("Kr density fluctuations(TM110)");
  }
  else if(mode==2){
    //gr->SetTitle("Drift of Kr density calibration(TM210)");
    gr->SetTitle("Quadratic pressure shift(TM210)");
    //gr->SetTitle("Kr density fluctuations(TM110)");
    //gr->SetTitle("Calibration of Kr density(TM210)");
  }
  gr->GetXaxis()->SetTitle("Pressure [atm]");
  if(mode==1)
    gr->GetYaxis()->SetTitle("Frequency - 1897431 [kHz]");
  else if(mode==2)
    gr->GetYaxis()->SetTitle("Frequency - 2565871 [kHz]");
  gr->GetYaxis()->SetTitleOffset(1.5);
  gr->SetMarkerColor(4);
  gr->SetMarkerStyle(21);
  f2->SetRange(0., 1.);
  gr->Draw("AP");
  gr->Fit("f2");

  if(mode==1)
    std::cout << "frequency shift : " << 1000*(f2->GetParameter(0)-1897431) << " [Hz]"<< std::endl;
  else if(mode==2)
    std::cout << "frequency shift : " << 1000*(f2->GetParameter(0)-2565871) << " [Hz]"<< std::endl;
  /*
  double high = x[1]*(y[1]+std::abs(ey[1])-y[0]+std::abs(ey[0]))/(x[0]-x[1]) + y[1]+std::abs(ey[1]);
  double low = x[1]*(y[1]-std::abs(ey[1])-y[0]-std::abs(ey[0]))/(x[0]-x[1]) + y[1]+std::abs(ey[1]);
  std::cout << "y+切片" << high << std::endl;
  std::cout << "y-切片" << low << std::endl;
  std::cout << "systematic error : " << (high -low)*1000 << " [Hz]" << std::endl;
  */
  
  if(mode==1){
    //c->SaveAs("../figure/Drift_of_Kr_density_calibration_tm110.png");
    //c->SaveAs("../figure/Calibration_of_Kr_density_tm110.png");
    //c->SaveAs("../figure/density_fluctuations_tm110.png");
    c->SaveAs("../figure/Quadratic_pressure_shift_tm110.png");
  }
  else if(mode==2){
    //c->SaveAs("../figure/Drift_of_Kr_density_calibration_210.png");
    //c->SaveAs("../figure/density_fluctuations_tm210.png");
    // c->SaveAs("../figure/Calibration_of_Kr_density_tm210.png");
    c->SaveAs("../figure/Quadratic_pressure_shift_tm210.png");
  }
  return 0;
}
