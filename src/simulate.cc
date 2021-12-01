////////////////////////////////////////////////////////
//   High Field Simulation for MuSUEUM Collaboration
//
//       Author : Hideharu Yamauchi 2021/9/23
////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include "stop.cc"
//#include "make_tree.cc"
//#include "make_tree_power.cc"
//#include "make_tree_dist.cc"
#include "make_tree_dist2.cc"
//#include "make_tree_angle.cc"

int main(int argc, const char** argv){
  if(argc!=4){
    std::cout << "argument lists:[./execute.out] [TMmode] [PATH/to/Geant4_simulation_result] [PATH/to/Geant4_simulation_enviroment]" << std::endl;
    return 0;
  }
  //-------------------------muon stopping distribution info 
  STOP* run = new STOP(argv[2], argv[3], atol(argv[1]));
  //TTree* tree = run->GetDecayTree();
  
  // if you want to visualize the RFpower histgram
  //run->Vis_RFPowerHist();

  // if you want to visualize the magnet field histgram 
  //run->Vis_FieldHist();

  // if you want to visualize the positron energy histgram
  //run->Vis_PositronEnergyHist();
  //run->Vis_PositronEnergyAtDetector();

  // if you want to visualize the positron angle histgram 
  //run->Vis_PositronAngleHist();
  
  // if you want to visualize muon stopping distribution
  //run->Vis_Stopping_DistXY(1030., 1350., true);  
  //run->Vis_Stopping_DistZ();

    
  //-------------------------superconductive magnet info
  //MAGNETFIELD* magnet = new MAGNETFIELD(atol(argv[1]));
  
  // if you want to visualize magnet field at Z       
  //magnet->Vis_MagField(-75.);
  
  
  //--------------------------radio frequency info(TMmode)
  //RFFIELD* RF = new RFFIELD(atol(argv[1]));
  
  // if you want to visualize RF field
  //RFFIELD->Vis_RF();

  TTree* tree = run->GetDecayTree(false);
  MAKETREE* create = new MAKETREE(tree, atol(argv[1]), "run01");
  
  delete run;
  return 0;
}
