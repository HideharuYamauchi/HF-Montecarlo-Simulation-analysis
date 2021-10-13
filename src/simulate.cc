////////////////////////////////////////////////////////
//   High Field Simulation for MuSUEUM Collaboration
//
//       Author : Hideharu Yamauchi 2021/9/23
////////////////////////////////////////////////////////
#ifndef ___simulation_
#define ___simulation_ 1

#include <stdio.h>
#include <stdlib.h>
#include "stop.cc"
#include "make_tree2.cc"
//#include "simulator.cc"
#endif

int main(int argc, const char** argv){
  if(argc!=5){
    std::cout << "argument lists:[./execute.out] [TMmode] [PATH/to/magnetic_moment_dist] [PATH/to/Geant4_simulation_result] [PATH/to/Geant4_simulation_enviroment]" << std::endl;
    return 0;
  }
  //-------------------------muon stopping distribution info 
  muonstopping* run = new muonstopping(argv[3], argv[4]);
  //TTree* tree = run->GetDecayTree(true);
  
  // if you want to visualize muon stopping distribution       
  //run->Vis_stopping_distXY(1030., 1350.);  
  //run->Vis_stopping_distZ();
  //----------------------------------------------------------------------

  
  //-------------------------superconductive magnet info
  // using copy constructor to initialize default constructor
  //magfield magnet = magfield();
  //magfield* magnet = new magfield(argv[2], atol(argv[1]));
  // magnet->AddMagnetBranch(tree);
  
  // if you want to visualize magnet field at Z       
  //magnet.Vis_magfield(-25.); 
  //----------------------------------------------------------------------

  
  //--------------------------radio frequency info(TMmode)
  //RFfield* RF = new RFfield(atol(argv[1]));
  
  // if you want to visualize RF field
  //RF->Vis_RF();

  // if you want to get the effective RF field
  //RF->Effective(run->Vis_stopping_distXY(1050.));
  //----------------------------------------------------------------------

  TTree* tree = run->GetDecayTree(false);
  maketree2* create = new maketree2(tree, atol(argv[1]), "run01"); 

  //-------------------------calculation
  //simulator* sim = new simulator("../data/run01.root");

  // if you want to get state amplitude
  //sim.timedev(double t, RF.Effective(run->Vis_stopping_distXY(1050.)), double delta, double gamma);
  
  //----------------------------------------------------------------------

  //delete RF;
  //delete magnet;
  //delete sim;
  delete run;
  return 0;
}
