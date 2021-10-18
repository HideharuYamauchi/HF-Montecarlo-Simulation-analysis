////////////////////////////////////////////////////////
//   High Field Simulation for MuSUEUM Collaboration
//
//       Author : Hideharu Yamauchi 2021/10/15
////////////////////////////////////////////////////////
#include <stdio.h>
#include "simulator.cc"

int main(int argc, const char** argv){
  SIMULATOR* sim = new SIMULATOR(argv[1]);
  //sim->CalculateSignal();

  delete sim;
  return 0;
}
