/////////////////////////////////////////////////////
//   High field simulation for MuSEUM Collaboration
//
//        Author: Hideharu Yamauchi 2021/10/12
/////////////////////////////////////////////////////
#ifndef ___class_maketree_
#define ___class_maketree_ 1

#include "../include/make_tree.hh"
#include "magnet.cc"
#include "RF.cc"
#endif

maketree::maketree(TTree* decaytree, int mode, std::string run_num)
  : decayvolume(0),decayvolume_branch(0),muon_position(0),muon_position_branch(0),muon_momentum(0),muon_momentum_branch(0),positron_position(0),positron_position_branch(0),field(3),field_branch(0),
    positron_energy_branch(0),decaytime_branch(0)
{
  filename = run_num;
  RF = new RFfield(mode);
  magnet = new magfield("../data/BRECON_MOM_20200716_6.txt", mode);
  decaytree->SetBranchAddress("decaytime",&decaytime,&decaytime_branch);
  decaytree->SetBranchAddress("decayvolume",&decayvolume,&decayvolume_branch);
  decaytree->SetBranchAddress("muon_position",&muon_position,&muon_position_branch);
  decaytree->SetBranchAddress("muon_momentum",&muon_momentum,&muon_momentum_branch);
  decaytree->SetBranchAddress("positron_position",&positron_position,&positron_position_branch);
  decaytree->SetBranchAddress("positron_momentum",&positron_momentum,&positron_momentum_branch);
  decaytree->SetBranchAddress("positron_energy",&positron_energy,&positron_energy_branch);
  decaytree->SetBranchStatus("*",1);
  
  // add new branch
  decaytree->Branch("field",&field);
  
  Long64_t tentry;
  Double_t X_temp, coefficientS, coefficientC, b;
  for(int n=0; n<decaytree->GetEntries(); n++){
    tentry = decaytree->LoadTree(n);
    muon_position_branch->GetEntry(tentry);
    magnet->GetDistance((*muon_position)[0], (*muon_position)[1], (*muon_position)[2]-cavity_center);
    field[0] = (magnet->B_ave + magnet->GetBfieldValue())*magnet->scaling_factor; // scaling magnet field to ~1.7
    X_temp = field[0]*(gfactor_j*magnetic_moment_j + gfactor_mu_prime*magnetic_moment_mu)/(plank_const*v_exp);
    coefficientS = sqrt(0.5)*sqrt(1-X_temp/sqrt(1+X_temp*X_temp));
    coefficientC = sqrt(0.5)*sqrt(1+X_temp/sqrt(1+X_temp*X_temp));
    if(mode==110) b = 0.001*0.25*(coefficientS*gfactor_j*magnetic_moment_j + coefficientC*gfactor_mu_prime*magnetic_moment_mu)/plank_const_divided;
    else if(mode==210) b = 0.001*0.25*(coefficientS*gfactor_j*magnetic_moment_j - coefficientC*gfactor_mu_prime*magnetic_moment_mu)/plank_const_divided;                         
    RF->GetXY((*muon_position)[0], (*muon_position)[1]);
    field[1] = RF->TM_mode();
    field[2] = b*field[1]; // kHz
    decaytree->Fill();
  }
  //decaytree->Scan("*");
  
  file = new TFile((run_num+".root").c_str(),"RECREATE");
  if(decaytree->Write()) std::cout  << run_num << ".root is made." << std::endl;
  file->Close();
}

maketree::~maketree(void){
  delete RF;
  delete magnet;
  delete file;
}
