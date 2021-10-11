/////////////////////////////////////////////////////////////
//      High field simulation for MuSEUM Collaboration
//
//            Author: Hideharu Yamauchi 2021/09/19
/////////////////////////////////////////////////////////////
#ifndef ___class_simulator_
#define ___class_simulator_ 1

#include <stdio.h>
#include <string>
#include "../include/simulator.hh"
#include "magnet.cc"
#include "TStyle.h"
#endif

simulator::simulator(TTree* decaytree, int mode, magfield* magnet)
  :Non(0.),Noff(0.),decayvolume(0), decayvolume_branch(0),muon_position(0),muon_position_branch(0),muon_momentum(0),muon_momentum_branch(0),positron_position(0),positron_position_branch(0)
{
  std::cout << magnet->B_ave << std::endl;
  Mode = mode;
  Double_t X_temp, coefficientS, coefficientC;
  decaytree->SetBranchAddress("decaytime",&decaytime);
  decaytree->SetBranchAddress("decayvolume",&decayvolume,&decayvolume_branch);
  decaytree->SetBranchAddress("muon_position",&muon_position,&muon_position_branch);
  decaytree->SetBranchAddress("muon_momentum",&muon_momentum,&muon_momentum_branch);
  decaytree->SetBranchAddress("positron_position",&positron_position,&positron_position_branch);
  decaytree->SetBranchAddress("positron_momentum",&positron_momentum,&positron_momentum_branch);
  decaytree->SetBranchAddress("positron_energy",&positron_energy);
  decaytree->SetBranchStatus("*",1);
  entries = decaytree->GetEntries();
  /*
  // add new branches
  auto magnet_field_Branch = decaytree->Branch("magnet_field",&magnet_field,"magnet_field/D");
  auto coefficientS_Branch = decaytree->Branch("coefficientS",&coefficientS,"coefficientS/D");
  auto coefficientC_Branch = decaytree->Branch("coefficientC",&coefficientC,"coefficientC/D");
  auto b_Branch = decaytree->Branch("b",&b,"b/D");
  */
  /*
  Long64_t tentry;
  for(int n=0; n<10; n++){
    tentry = decaytree->LoadTree(n);
    muon_position_branch->GetEntry(tentry);
    //magnet->GetDistance((*muon_position)[0], (*muon_position)[1], (*muon_position)[2]-cavity_center);
    //magnet_field = (magnet->B_ave+magnet->GetBfieldValue())*magnet->scaling_factor; // scaling magnet field to ~1.7
    std::cout << magnet->GetDistance((*muon_position)[0], (*muon_position)[1], (*muon_position)[2]-cavity_center) << std::endl;
    //magnet_field_Branch->Fill();
      
    X_temp = magnet_field*(gfactor_j*magnetic_moment_j + gfactor_mu_prime*magnetic_moment_mu)/(plank_const*v_exp);
    coefficientS = sqrt(0.5)*sqrt(1-X_temp/sqrt(1+X_temp*X_temp));                                                                                         
    coefficientS_Branch->Fill();                                                                                                                                         
    coefficientC = sqrt(0.5)*sqrt(1+X_temp/sqrt(1+X_temp*X_temp));
    coefficientC_Branch->Fill();
    if(mode==110){                                                                                                                                                                           
      b = 0.001*0.25*(coefficientS*gfactor_j*magnetic_moment_j + coefficientC*gfactor_mu_prime*magnetic_moment_mu)/plank_const_divided;
      b_Branch->Fill();                                                                    
    }else if(mode==210){
      b = 0.001*0.25*(coefficientS*gfactor_j*magnetic_moment_j - coefficientC*gfactor_mu_prime*magnetic_moment_mu)/plank_const_divided;
      b_Branch->Fill();                                                                                                                                                                           
    }
  }
  */
  decaytree->Scan("*");
  // initial state amplitude from MuSEUM technical note (2.12)
  /*
  if(mode==110){
    state_amp[0]=0.25*(1+polarization);
    state_amp[1]=0.25*(1-(pow(coefficientC,2.)-pow(coefficientS,2.))*polarization);
  }else if(mode==210){
    state_amp[0]=0.25*(1-polarization);
    state_amp[1]=0.25*(1+(pow(coefficientC,2.)-pow(coefficientS,2.))*polarization);
  }
  std::cout << "state amplitue(t=0) of 12transition=" << state_amp[0] << "\n"
            << "state amplitue(t=0) of 34transition=" << state_amp[1] << std::endl;
  */
  delete magnet;
}

void simulator::timedev(double t, double b, double delta, double gamma, double position[3]){
  state_amp[0]=(pow(state_amp[0],2.)*(pow(std::cos(0.5*gamma*t),2.)+pow((delta/gamma)*std::sin(0.5*gamma*t),2.))+pow(state_amp[1]*2*b*std::sin(0.5*gamma*t)/gamma,2.)+state_amp[0]*state_amp[1]*2*delta*b*pow(std::sin(0.5*gamma*t)/gamma,2.))*std::exp(-muon_life*t);
  state_amp[1]=(pow(state_amp[0]*2*b*std::sin(0.5*gamma*t)/gamma,2.)+pow(state_amp[1],2.)*(pow(std::cos(0.5*gamma*t),2.)+pow(delta*std::sin(0.5*gamma*t)/gamma,2.))-state_amp[0]*state_amp[1]*2*delta*b*pow(std::sin(0.5*gamma*t)/gamma,2.))*std::exp(-0.5*muon_life*t);
  std::cout << "state amplitue(t=" << t << ") of 12transition=" << state_amp[0] << "\n"
            << "state amplitue(t=" << t << ") of 34transition=" << state_amp[1] << std::endl;
}
/*
void simulator::Vis_state_amp(void){
  f1 = new TF1("amplitude","time development",0,10,2);
  //f1->SetParameters(2,1);             
  //f1->SetParNames("constant","coefficient");                                                      
  //f1->Draw();
  delete f1;
}

double simulator::Calculate_g(Double_t Gamma, Double_t t){
  double g = std::cos(Gamma*t)-Gamma*std::sin(Gamma*t)/muon_life;
  return g;
}

double simulator::ConventionalSignal(double power, double detuning, double position[3]){
  double Gamma_square = pow(detuning, 2.)+4*pow(power, 2.);
  L = 2*pow(power, 2.)/(Gamma_square + pow(gamma, 2.));
  //K = ;
  return 0;
}

double simulator::OldMuoniumSignal(double power, double detuning, double position[3], Double_t windowopen, Double_t windowclose){
  double Gamma_square = pow(detuning, 2.)+4*pow(power, 2.);
  L = 2*pow(power, 2.)*(
			std::exp(-muon_life*windowopen)*(1-Calculate_g(sqrt(Gamma_square)*windowopen)*pow(muon_life,2.)/(Gamma_square+pow(muon_life,2.)))
			-std::exp(-muon_life*windowclose)*(1-Calculate_g(sqrt(Gamma_square)*windowopen)*pow(muon_life,2.)/(Gamma_square+pow(muon_life,2.)))
			)/Gamma_square;
  //K = ;
  return 0;
}

double simulator::Calculate_Signal(void){
  for(int i=0; i<entries; i++){
    ;
  }
  return 0;
}
*/
