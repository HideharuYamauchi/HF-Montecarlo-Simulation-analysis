/////////////////////////////////////////////////////
//   High field simulation for MuSEUM Collaboration
//
//        Author: Hideharu Yamauchi 2021/10/13
/////////////////////////////////////////////////////
#include <stdio.h>
#include "../include/make_tree2.hh"

#ifndef ___class_muonstopping_
#define ___class_muonstopping_ 1
#include "magnet.cc"
#include "RF.cc"
#endif

#include "TString.h"

MAKETREE::MAKETREE(TTree* decaytree, int mode, std::string run_num)
  : decayvolume(0),decayvolume_branch(0),muon_position(0),muon_position_branch(0),muon_momentum(0),muon_momentum_branch(0),muon_energy_branch(0),positron_position(0),positron_position_branch(0),
    positron_momentum(0),positron_momentum_branch(0),positron_energy_branch(0),decaytime_branch(0),
    str_vec(4),muon_vec(4),muon_dispersion(4),positron_vec(4),positron_dispersion(4),field(3),state_amp(4),
    str_branch(0),muon_vec_branch(0),muon_dispersion_branch(0),positron_vec_branch(0),positron_dispersion_branch(0),field_branch(0),state_amp_branch(0)
{
  TString path = "../data/" + run_num + ".root";
  if(gSystem->GetPathInfo(path, info)==0) std::cout << run_num + " is already exist" << std::endl;
  else{
    flag = true;
    DecayTree = new TTree("DecayTree","tree of decay muons");
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
    DecayTree->Branch("str_vec",&str_vec);
    DecayTree->Branch("muon_vec",&muon_vec);
    DecayTree->Branch("muon_dispersion",&muon_dispersion);
    DecayTree->Branch("positron_vec",&positron_vec);
    DecayTree->Branch("positron_dispersion",&positron_dispersion);
    DecayTree->Branch("field",&field);
    DecayTree->Branch("state_amp",&state_amp);
  
    Int_t entries = decaytree->GetEntries();
    Double_t X_temp, coefficientS, coefficientC, b;
    for(int n=0; n<entries; n++){
      decaytree->GetEntry(n);
      str_vec[0] = *decayvolume;
      if(n==0){
	str_vec[1] = "TM"+std::to_string(mode)+"mode";
	str_vec[2] = "1*atmosphere";
	str_vec[3] = "300*kelvin";
      }else if(n!=0){
	str_vec[1] = "";
	str_vec[2] = "";
	str_vec[3] = "";
      }
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
      field[0] = (magnet->B_ave + magnet->GetBfieldValue())*magnet->scaling_factor; // scaling magnet field to ~1.7
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
      field[2] = b*field[1]; // kHz
      DecayTree->Fill();
    }
    //decaytree->Scan("*");
  
    file = new TFile(("../data/"+run_num+".root").c_str(),"RECREATE"); //std::string
    if(DecayTree->Write()) std::cout  << run_num << ".root is made." << std::endl;
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
