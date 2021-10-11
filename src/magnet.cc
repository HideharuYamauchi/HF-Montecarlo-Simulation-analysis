/////////////////////////////////////////////////////////          
//         Class for calculation of Magnet field     
//                                                  
//        Author:   Hideharu Yamauchi 2021/09/16        
/////////////////////////////////////////////////////////
#ifndef ___class_magfield_
#define ___class_magfield_ 1

#include <stdio.h>
#include <string>
#include "../include/magnet.hh"
#include <cstdio>
#include <fstream>
#include <iostream>
#include <numeric>
#endif

magfield::magfield(const char* magnetfile, int Mode):moment_coordinate(moment_num,std::vector<double>(3)),moment(moment_num,std::vector<double>(3)),distance(moment_num,std::vector<double>(3)),interval(moment_num),position(3){
  mode = Mode;
  std::ifstream ifs(magnetfile);
  if(ifs.fail()){
    std::cout << "Failed to open the file..." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  int i=0;
  int dummy;
  while(!ifs.eof()){
    if(i==moment_num) break;
    ifs >> dummy >> moment_coordinate[i][0] >> moment_coordinate[i][1] >> moment_coordinate[i][2] >> moment[i][0] >> moment[i][1] >> moment[i][2];
    i++;
  }
}

double magfield::GetDistance(double x, double y, double z){
  if(sqrt(pow(x, 2.)+pow(y, 2.)+pow(z, 2.)) < 300.){ 
    position[0] = x*1.0e-3; // convert mm to m
    position[1] = y*1.0e-3;
    position[2] = z*1.0e-3;
    for(int i=0;i<moment_num;i++){
      for(int j=0;j<3;j++){
	distance[i][j] = position[j]-moment_coordinate[i][j];
      }
      interval[i] = sqrt(pow(distance[i][0], 2.)+pow(distance[i][1], 2.)+pow(distance[i][2], 2.));
    }
  }else{std::cout << "distance is over the DSV" << DSV << "(mm)..." << std::endl;}
  return sqrt(pow(x, 2.)+pow(y, 2.)); // mm
}

double magfield::GetBfieldValue(void){
  double BvalueZ = 0, total_BvalueZ = 0;
  for(int i=0;i<moment_num;i++){
    BvalueZ = (1.0e-7)*(1/pow(interval[i], 3.))*(3*(moment[i][0]*distance[i][0] + moment[i][1]*distance[i][1] + moment[i][2]*distance[i][2])*distance[i][2]/pow(interval[i], 2.)-moment[i][2]);
    total_BvalueZ += BvalueZ;
  }
  return total_BvalueZ; // T
}

void magfield::Vis_magfield(double Z){ // unit of Z:mm, range: from -152 to +152
  char title_dt[48];
  sprintf(title_dt,"z=%.1f [/mm]; X (mm); Y (mm); B^{ER} [/Gauss]", Z);
  char title_dt2[16];
  sprintf(title_dt2,"z=%.1f [/mm]", Z);
  char title_dt3[26];
  sprintf(title_dt3,"magnet_field_z=%.1f.png", Z);
  
  c = new TCanvas("c", "c",1600,600);
  gStyle->SetOptStat(0);
  gStyle->SetTitleXOffset(1.5);
  gStyle->SetTitleYOffset(2); 
  center_pad = new TPad("center_pad", "center_pad",0.5,0,1.0,1.0);
  center_pad->Draw();
  top_pad = new TPad("top_pad", "top_pad",0,0,0.5,1.0);
  top_pad->Draw();
  
  dt = new TGraph2D();
  dt->SetTitle(title_dt);
  dt2 = new TH2D("dt", title_dt2, 201,-100.0 , 100.0 ,201, -100.0, 100.0);
  dt2->SetXTitle("X [/mm]");
  dt2->SetYTitle("Y [/mm]");
  dt2->SetZTitle("B^{ER} [/Gauss]");
  int i=0;
  for(int X=-100; X<101; X++){
    for(int Y=-100; Y<101; Y++){      
      if(GetDistance(double(X),double(Y),Z) < cavity_radius*1.0e+3){
	dt->SetPoint(i, double(X), double(Y), GetBfieldValue()*1.0e+4); // convert T to Gauss
        dt2->Fill(X, Y, GetBfieldValue()*1.0e+4); // convert T to Gauss
	i++;
      }
    }
  }
  top_pad->cd();
  gStyle->SetPalette(1);
  dt->Draw("surf1");
  dt->GetHistogram()->GetZaxis()->SetTitleOffset(1.4);
  center_pad->cd();
  center_pad->SetGridx(1);
  center_pad->SetGridy(1);
  center_pad->SetLeftMargin(0.2);
  center_pad->SetRightMargin(0.2);
  dt2->Draw("colz");
  dt2->GetZaxis()->SetTitleOffset(1.3);
  c->SaveAs(title_dt3); 
  delete center_pad;
  delete top_pad;
  delete dt;
  delete dt2;
  delete c;
}

TTree* magfield::AddMagnetBranch(TTree* decaytree){
  Double_t magnet_field, b;
  Double_t decaytime, positron_energy;
  std::string* decayvolume = 0;
  std::vector<Double_t>* muon_position = 0;
  TBranch* muon_position_branch = 0;
  std::vector<Double_t>* muon_momentum = 0;
  TBranch* muon_momentum_branch = 0;
  std::vector<Double_t>* positron_position = 0;
  TBranch* positron_position_branch = 0;
  std::vector<Double_t>* positron_momentum = 0;
  TBranch* positron_momentum_branch = 0;
  Double_t X_temp, coefficientS, coefficientC;
  decaytree->SetBranchAddress("decaytime",&decaytime);
  decaytree->SetBranchAddress("decayvolume",&decayvolume);
  decaytree->SetBranchAddress("muon_position",&muon_position,&muon_position_branch);
  decaytree->SetBranchAddress("muon_momentum",&muon_momentum,&muon_momentum_branch);
  decaytree->SetBranchAddress("positron_position",&positron_position,&positron_position_branch);
  decaytree->SetBranchAddress("positron_momentum",&positron_momentum,&positron_momentum_branch);
  decaytree->SetBranchAddress("positron_energy",&positron_energy);
  decaytree->SetBranchStatus("*",1);
  auto magnet_field_Branch = decaytree->Branch("magnet_field",&magnet_field,"magnet_field/D");
  auto coefficientS_Branch = decaytree->Branch("coefficientS",&coefficientS,"coefficientS/D");
  auto coefficientC_Branch = decaytree->Branch("coefficientC",&coefficientC,"coefficientC/D");
  auto b_Branch = decaytree->Branch("b",&b,"b/D");
 
  for(int n=0; n<10; n++){
    Long64_t tentry = decaytree->LoadTree(n);
    muon_position_branch->GetEntry(tentry);
    //std::cout << (*muon_position)[2] << std::endl;
    std::cout << (*muon_position)[0]  <<"\t"<< (*muon_position)[1] <<"\t"<< (*muon_position)[2] <<std::endl;
    /*   
    GetDistance((*muon_position)[0], (*muon_position)[1], (*muon_position)[2]);
    magnet_field = (B_ave+GetBfieldValue())*scaling_factor; // scaling magnet field to ~1.7
    magnet_field_Branch->Fill();
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
    */
  }
  //decaytree->Scan("*");
  return decaytree;  
}
