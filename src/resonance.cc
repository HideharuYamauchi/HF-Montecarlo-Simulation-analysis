////////////////////////////////////////////////////////
//   High Field Simulation for MuSUEUM Collaboration
//
//       Author : Hideharu Yamauchi 2021/10/15
////////////////////////////////////////////////////////
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "simulator.cc"

int main(int argc, const char** argv){
  SIMULATOR* sim = new SIMULATOR(argv[1]);
  sim->CalculateSignal(true, 230., 4., 5.);
  //sim->Vis_StateAmp(0.);

  int power_min = 100;
  int power_max = 240;
  int power;

  /*
  for(int i=0; i<1+(power_max-power_min)/10; i++){                                                                                                                                                       
      power = power_min + 10*i;                                                                                                                                
      sim->CalculateSignal(true, power, 9., 10.);
  }
  */
  
  //std::ofstream ofs("../figure/FWHM_conv.dat", std::ios::out);
  //std::ofstream ofs2("../figure/FWHM_conv_hight.dat", std::ios::out);
  //std::ofstream ofs("../figure/FWHM_old.dat", std::ios::out);
  //std::ofstream ofs2("../figure/FWHM_old_hight.dat", std::ios::out);
  
  /*
  if(ofs){
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
    for(int t=3; t<13; t++){
      sim->CalculateSignal(true, 125., t, t+1);
      ofs << t << "\t"
          << sim->Sim_FWHM << "\t"
	  << sim->The_FWHM << std::endl;
    }                                                                                                                                              
    ofs.close();
  }
  
  std::ofstream ofs("../figure/merit_290.dat", std::ios::out);
  Double_t det;
  if(ofs){
    std::cout << "Successful to open file." << std::endl;                                                                                                                                                 
    for(int t=2; t<12; t++){
      sim->CalculateSignal(true, 290., t, t+1);
      ofs << t << "\t"
          << sim->Sim_FWHM << "\t"
	  << sim->Sim_Height << "\t"
	  << sim->Sim_detected << std::endl;
    }
    ofs.close();
  }
  
  std::ofstream ofs("../figure/t1_height.dat", std::ios::out);
  if(ofs){
    std::cout << "Successful to open file." << std::endl;                                                                                                                                                 
    for(int t=0; t<13; t++){
      sim->CalculateSignal(true, 240., t, t+1);   
      ofs << t << "\t"   
          << sim->Sim_Height << "\t"                                                                                                                                                                      
          << sim->The_Height << std::endl;
    }
    ofs.close();                                                                                                                                                                               
  }
  */  

  delete sim;
  return 0;
}
