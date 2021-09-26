///////////////////////////////////////////////////
//   High field simulation for MuSEUM Experiment 
//
//      Author:  Hideharu Yamauchi 2021/09/18 
///////////////////////////////////////////////////

#if !defined(___header_simulator_)
#define ___header_simulator_ 1

#include "TMath.h"
#include "TTree.h"
#include "TF1.h"
#include "./magnet.hh" // using B_cons
#include "./HFgeometry.hh"
#include "./RF.hh"
#include "./stop.hh"

class simulator{
private:
  std::vector<double> position;
  double angle;
  double solid_angle;
  int Mode;
  double state_amp[2];
  double Btob;
  
public:
  simulator(int mode){Mode=mode;};
  simulator(TTree* runtree, int mode);
  ~simulator(void){;};
  TF1* f1;
  void GetXYZ(int x, int y, int z);
  double* timedev(double t, double b, double delta, double gamma);
  void Vis_state_amp(void);
  void CalculateBtob(void);
  void CalculateEffectiveb(void);
};
#endif
