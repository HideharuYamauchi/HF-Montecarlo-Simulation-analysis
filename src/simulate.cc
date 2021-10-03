////////////////////////////////////////////////////
//  High Field Simulation for MuSUEUM Collaboration
//
//      Author: Hideharu Yamauchi 2021/9/23
////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include "RF.cc"
#include "magnet.cc"
#include "stop.cc"
#include "simulator.cc"

int main(int argc, const char** argv){
  /*
  if(argc!=4){
    std::cout << "argument lists:[TMmode] [PATH/to/magnetic_moment_dist] [PATH/to/Geant4_simulation_result] [PATH/to/Geant4_simulation_enviroment]" << std::endl;
    return 0;
  }
  */
  //-------------------------muon stopping distribution info 
  muonstopping* run = new muonstopping(argv[3], argv[4]);
  
  // if you want to visualize muon stopping distribution       
  //run->Vis_stopping_distXY(1050.);        
  //run->Vis_stopping_distZ();
  //----------------------------------------------------------------------

  
  //--------------------------radio frequency info(TMmode)
  RFfield* RF = new RFfield(atol(argv[1]));
  //RF->AddRFBranch(run->GetDecayTree());
  
  // if you want to visualize RF field
  //RF->Vis_RF();

  // if you want to get the effective RF field
  //RF->Effective(run->Vis_stopping_distXY(1050.));
  //----------------------------------------------------------------------

  
  //-------------------------superconductive magnet info
  // using copy constructor to initialize default constructor
  //magfield magnet = magfield();
  //magfield* magnet = new magfield(argv[2]);
  //magnet->AddMagnetBranch(RF->AddRFBranch(run->GetDecayTree()));
  
  // if you want to visualize magnet field at Z
  //magnet.Vis_magfield(-25.);
  //----------------------------------------------------------------------

  
  //-------------------------calculation
  //simulator sim(atol(argv[1]), &field_RF);
  //simulator sim(atol(argv[1]));

  // if you want to get state amplitude
  //sim.timedev(double t, RF.Effective(run->Vis_stopping_distXY(1050.)), double delta, double gamma);
  
  //----------------------------------------------------------------------

  delete run;
  return 0;
}
