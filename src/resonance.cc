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
  sim->CalculateSignal();
  //sim->Vis_StateAmp(0.);

  int power_min = 100;
  int power_max = 250;
  int power;
  /*
  std::ofstream ofs("../figure/FWHM_conv.dat", std::ios::out);
  std::ofstream ofs2("../figure/FWHM_conv_hight.dat", std::ios::out);
  //std::ofstream ofs("../figure/FWHM_old.dat", std::ios::out);
  //std::ofstream ofs2("../figure/FWHM_old_hight.dat", std::ios::out);
    
  if(ofs){
    std::cout << "Successful to open file." << std::endl;
    for(int i=0; i<1+(power_max-power_min)/10; i++){
      power = power_min + 10*i;
      sim->CalculateSignal(true, power);
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
  
  std::ofstream ofs("../figure/merit.dat", std::ios::out);
  Double_t det;
  
  if(ofs){
    std::cout << "Successful to open file." << std::endl;                                                                                                                                                 
    for(int t=1; t<13; t++){
      sim->CalculateSignal(true, 125., t, t+1);
      ofs << t << "\t"
          << sim->Sim_FWHM << "\t"
          << sim->The_FWHM << "\t"
	  << sim->Sim_Height << "\t"
	  << sim->The_Height << "\t"
	  << sim->detected << std::endl;
    }
    ofs.close();
  }
  */
  delete sim;
  return 0;
}
