/////////////////////////////////////////////////////
//   High field simulation for MuSEUM Collaboration
//
//        Author: Hideharu Yamauchi 2021/10/13
/////////////////////////////////////////////////////
#include <stdio.h>
#include <random>
#include "../include/make_tree.hh"

#ifndef ___class_muonstopping_
#define ___class_muonstopping_ 1
#include "magnet.cc"
#include "RF.cc"
#endif

#include "TString.h"

MAKETREE::MAKETREE(TTree* decaytree, int mode, std::string run_num)
  : decayvolume(0),decayvolume_branch(0),muon_position(0),muon_position_branch(0),muon_momentum(0),muon_momentum_branch(0),muon_energy_branch(0),positron_position(0),positron_position_branch(0),
    positron_momentum(0),positron_momentum_branch(0),positron_energy_branch(0),decaytime_branch(0),
    str_vec(4),muon_vec(4),muon_dispersion(4),positron_vec(4),positron_dispersion(4),field(3),state_amp(4),angle_vec(2),position(3),
    str_branch(0),muon_vec_branch(0),muon_dispersion_branch(0),positron_vec_branch(0),positron_dispersion_branch(0),field_branch(0),state_amp_branch(0),angle_branch(0)
{
  TString path = "../data/" + run_num + "_power.root";
  TString MODE;
  TString ATM;
  TString KEL;
  double power;

  const Double_t instability = 0.02;
  Double_t shift;
  std::random_device rnd;
  std::mt19937 mt(rnd());
  std::uniform_real_distribution<> randpower(1-instability/100, 1+instability/100);
  
  if(gSystem->GetPathInfo(path, info)==0) std::cout << run_num + "_power is already exist" << std::endl;
  else{    
    flag = true;
    MODE = "TM"+std::to_string(mode)+"mode";
    ATM = "1.0*atmosphere";
    KEL = "300*kelvin";
    DecayTree = new TTree("DecayTree",MODE+":"+ATM+":"+KEL);
    RF = new RFFIELD(mode);
    magnet = new MAGNETFIELD(mode);
      
    decaytree->SetBranchAddress("decaytime",&decaytime,&decaytime_branch);
    decaytree->SetBranchAddress("decayvolume",&decayvolume,&decayvolume_branch);
    decaytree->SetBranchAddress("muon_position",&muon_position,&muon_position_branch);
    decaytree->SetBranchAddress("muon_momentum",&muon_momentum,&muon_momentum_branch);
    decaytree->SetBranchAddress("muon_energy",&muon_energy,&muon_energy_branch);
    decaytree->SetBranchAddress("positron_position",&positron_position,&positron_position_branch);
    decaytree->SetBranchAddress("positron_momentum",&positron_momentum,&positron_momentum_branch);
    decaytree->SetBranchAddress("positron_energy",&positron_energy,&positron_energy_branch);
    decaytree->SetBranchStatus("*",1);
    
    // add new branch
    DecayTree->Branch("muon_vec",&muon_vec);
    DecayTree->Branch("muon_dispersion",&muon_dispersion);
    DecayTree->Branch("positron_vec",&positron_vec);
    DecayTree->Branch("positron_dispersion",&positron_dispersion);
    DecayTree->Branch("field",&field);
    DecayTree->Branch("state_amp",&state_amp);
    DecayTree->Branch("Angle", &angle_vec);
  
    Int_t entries = decaytree->GetEntries();
    Double_t X_temp, coefficientS, coefficientC, b;
    for(int n=0; n<entries; n++){
      decaytree->GetEntry(n);
      for(int i=0;i<4;i++){
	if(i==0){
	  muon_vec[i] = decaytime;
	  muon_dispersion[i] = muon_energy;
	  positron_vec[i] = decaytime;
	  positron_dispersion[i] = positron_energy;
	}else if(i!=0){
	  muon_vec[i] = (*muon_position)[i-1];
	  muon_dispersion[i] = (*muon_momentum)[i-1];
	  positron_vec[i] = (*positron_position)[i-1];
	  positron_dispersion[i] = (*positron_momentum)[i-1];
	}
      }
      magnet->GetDistance((*muon_position)[0], (*muon_position)[1], (*muon_position)[2]-cavity_center);
      field[0] = B_cons; // scaling magnet field to ~1.7
      X_temp = field[0]*(gfactor_j*magnetic_moment_j + gfactor_mu_prime*magnetic_moment_mu)/(plank_const*v_exp);
      coefficientS = sqrt(0.5)*sqrt(1-X_temp/sqrt(1+X_temp*X_temp));
      coefficientC = sqrt(0.5)*sqrt(1+X_temp/sqrt(1+X_temp*X_temp));
      if(mode==110)
	b = 0.001*0.25*(coefficientS*gfactor_j*magnetic_moment_j + coefficientC*gfactor_mu_prime*magnetic_moment_mu)/plank_const_divided;
      else if(mode==210)
	b = 0.001*0.25*(coefficientS*gfactor_j*magnetic_moment_j - coefficientC*gfactor_mu_prime*magnetic_moment_mu)/plank_const_divided;
      // initial state amplitude from MuSEUM technical note (2.12)
      state_amp[0]=0.25*(1+polarization);
      state_amp[1]=0.25*(1-(pow(coefficientC,2.)-pow(coefficientS,2.))*polarization);
      state_amp[2]=0.25*(1-polarization);
      state_amp[3]=0.25*(1+(pow(coefficientC,2.)-pow(coefficientS,2.))*polarization);
      RF->GetXY((*muon_position)[0], (*muon_position)[1]);
      field[1] = RF->TM_mode();
      shift = randpower(mt); // add 0.02% power instability
      field[2] = b*field[1]*(1+shift); // kHz
      position[0] = (*muon_position)[0];
      position[1] = (*muon_position)[1];
      position[2] = (*muon_position)[2];
      CalculateAngle();
      angle_vec[0] = cos_solidangle;
      angle_vec[1] = solidangle;
      if(n%1000==0) std::cout << "Present Entry: " << n << std::endl;
      DecayTree->Fill();
    }
    //decaytree->Scan("*");
  
    file = new TFile(("../data/"+run_num+"_power.root").c_str(),"RECREATE");
    
    if(DecayTree->Write()) std::cout  << run_num << "_power.root is made." << std::endl;
    file->Close();
  }
}

MAKETREE::~MAKETREE(void){
  if(flag){
    delete RF;
    delete magnet;
    delete file;
  }
}

void MAKETREE::CalculateAngle(){
  Double_t x = -119.5, y = 119.5;
  Double_t ratio = (DetectorD_center-position[2])/(cavity_downfoil_center-position[2]);
  std::vector<Double_t> positron_at_foil(2);
  Double_t r, R;
  Double_t norm;
  std::vector<Double_t> distance(3);
  std::vector<Double_t> basic_vector(3);
  cos_solidangle = 0.;
  solidangle = 0.;

  for(int i=0; i<240; i++){ // get detector's x position
    x = -119.5+i;
    for(int n=0; n<240; n++){ // get detector's y position
      y = 119.5-n;
      
      distance[0] = x-position[0];
      distance[1] = y-position[1];
      distance[2] = DetectorD_center-position[2];
      norm = sqrt(pow(distance[0],2.) + pow(distance[1],2.) + pow(distance[2],2.));
      basic_vector[0] = distance[0]/norm;
      basic_vector[1] = distance[1]/norm;
      basic_vector[2] = distance[2]/norm;
      positron_at_foil[0] = position[0]+basic_vector[0]*ratio;
      positron_at_foil[1] = position[1]+basic_vector[1]*ratio;
      r = sqrt(pow(positron_at_foil[0],2.)+pow(positron_at_foil[1],2.));

      if(0<=r&&r<=cavity_radius*1.e+3){ // discriminate projection position is on cavity foil or not
        R = sqrt(pow(distance[0],2.)+pow(distance[1],2.)+pow(distance[2],2.));
        cos_solidangle += (distance[2]/R)*distance[2]*pow(pow(R,2.),-1.5);
        solidangle += distance[2]*pow(pow(R,2.),-1.5);
      }
      else if(cavity_radius*1.e+3<r){
        cos_solidangle += 0.;
        solidangle += 0.;
      }
    }
  }
}
