/////////////////////////////////////////////////////        
//   High field simulation for MuSEUM Collaboration             
//                                               
//         Author: Hideharu Yamauchi 2021/09/19                          
/////////////////////////////////////////////////////               
#ifndef ___class_muonstopping_
#define ___class_muonstopping_ 1

#include <stdio.h>
#include <ostream>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include "TMath.h"
#include "../include/stop.hh"
#include "TStyle.h"
#include "TPad.h"
#endif

muonstopping::muonstopping(std::string runfile, const char* envfile):tree(nullptr),branchtime(nullptr),branchX(nullptr),branchY(nullptr),branchZ(nullptr),branchPx(nullptr),branchPy(nullptr),branchPz(nullptr),branchkE(nullptr),branchdepE(nullptr),branchtrack(nullptr),branchstep(nullptr),branchcopyno(nullptr){
  TString runfile2 = runfile;
  run_num = runfile.substr(runfile.find("run"), 7); // get the present run number
  
  tree = new TTree("tree",run_num);
  tree->ReadFile(runfile2.Data(),"time/D:particle/C:volume:process:X/D:Y:Z:Px:Py:Pz:kE:depE:track/I:step:copyno");
  tree->SetBranchAddress("time",&time, &branchtime);
  tree->SetBranchAddress("particle",particle, &branchparticle);
  tree->SetBranchAddress("volume",volume, &branchvolume);
  tree->SetBranchAddress("process",process, &branchprocess);
  tree->SetBranchAddress("X",&X, &branchX);
  tree->SetBranchAddress("Y",&Y, &branchY);
  tree->SetBranchAddress("Z",&Z, &branchZ);
  tree->SetBranchAddress("Px",&Px, &branchPx);
  tree->SetBranchAddress("Py",&Py, &branchPy);
  tree->SetBranchAddress("Pz",&Pz, &branchPz);
  tree->SetBranchAddress("kE",&kE, &branchkE);
  tree->SetBranchAddress("depE",&depE, &branchdepE);
  tree->SetBranchAddress("track",&ntrack, &branchtrack);
  tree->SetBranchAddress("step",&nstep, &branchstep);
  tree->SetBranchAddress("copyno",&copyno, branchcopyno);    
  entries = tree->GetEntries();
  nbranches = tree->GetListOfBranches()->GetEntriesFast();
}

void muonstopping::CreateRootFile(void){
  file = TFile::Open(run_num+=".root","RECREATE");
  if(tree->Write()) std::cout << run_num << " is made." << std::endl;             
  file->Close();
  delete file;
}

void muonstopping::Vis_stopping_distXY(Double_t posZ){
  TString title_dt = "z=" + std::to_string(int(posZ)) + " [/mm]";
  TString title_dt2 = "XY-Distribution_z:" + std::to_string(int(posZ));
  c = new TCanvas("c", "c",900,900);
  TPad* center_pad = new TPad("center_pad", "center_pad",0.0,0.0,0.5,0.5);
  center_pad->Draw();
  TPad* right_pad = new TPad("right_pad", "",0.5,0.0,1.0,0.5);
  right_pad->Draw();
  TPad* top_pad = new TPad("top_pad", "",0.0,0.5,0.5,1.0);
  top_pad->Draw();
  dtxy = new TH2D("XY-Dist", title_dt, 240, -120, 120, 240, -120, 120);
  center_pad->cd();
  dtxy->SetXTitle("X [/mm]");
  dtxy->SetYTitle("Y [/mm]");                         
  center_pad->SetLeftMargin(0.15);
  center_pad->SetRightMargin(-0.03);
  
  for(int n=0;n<entries;n++){
    tree->GetEntry(n);
    if((Z<=posZ)&&(posZ<=Z+1.)&&(std::string(particle)=="mu+")&&(std::string(process)=="DecayWithSpin")) dtxy->Fill(X,Y);
  }
  TH1D* projdtx = dtxy->ProjectionX();
  projdtx->SetTitle("X Projection");
  projdtx->SetXTitle("X [/mm]");
  TH1D* projdty = dtxy->ProjectionY();
  projdty->SetTitle("Y Projection");
  projdty->SetXTitle("Y [/mm]");

  center_pad->cd();
  gStyle->SetPalette(1);
  dtxy->Draw("Colz");
  dtxy->GetZaxis()->SetTitleOffset(1.3);
  dtxy->SetStats(1); // set the stats table
  gStyle->SetStatX(0.9); // set the xposition of stats table
  gStyle->SetStatY(0.95); 
  top_pad->cd();
  projdtx->SetFillColor(kBlue+1);
  top_pad->SetLeftMargin(0.15);
  top_pad->SetRightMargin(-0.03);
  projdtx->Draw("bar");
  projdtx->SetStats(0); // do not set stats table
  right_pad->cd();
  projdty->SetFillColor(kBlue-2);
  right_pad->SetLeftMargin(0.15);
  right_pad->SetRightMargin(0.03);
  projdty->Draw("hbar");
  projdty->SetStats(0);

  c->SaveAs(title_dt2+=".png");
  delete dtxy;
  delete c;
}

void muonstopping::Vis_stopping_distZ(){
  TString title_dt = "Z-Distribution";
  c2 = new TCanvas("c2", "c2",1400,900);
  gStyle->SetOptStat(0); // do not set the stat tabel 
  dtz = new TH2D("Z-Dist", "", 400, 1000, 1400, 240, -120, 120);
  //dtz = new TH2D("Z-Dist", "", 4000, -2000, 2000, 240, -120, 120);
  dtz->SetXTitle("Position on the beam axis [/mm]");
  dtz->SetYTitle("Position on the vertical axis [/mm]");
  for(Int_t n=0;n<entries;n++){
    tree->GetEntry(n);
    if(std::string(process)=="DecayWithSpin") dtz->Fill(Z,X);
    //if(std::string(particle)=="mu+") dtz->Fill(Z,X);
  }
  dtz->Draw("Colz");
  dtz->GetXaxis()->SetTitleOffset(1.3);
  dtz->GetYaxis()->SetTitleOffset(1.2);
  c2->SaveAs(title_dt+=".png");
  delete dtz;
  delete c2;
}

int* muonstopping::GetMuonDist(void){
  int* dist = nullptr;
  int number=0,mu_number=0,e_number=0;
  /*
  for(int n=0;n<entries;n++){
    tree->GetEntry(n);
    //if((X<*pos)&&(*pos<X+1)&&(Y<*(pos+1))&&(*(pos)<Y+1)&&(Z<*(pos+2))&&(*(pos+2)<Z+1)&&((std::string(particle)=="mu+")&&(std::string(process)=="DecayWithSpin")) number++;
    //if(number<=nstep) number=nstep;
    if(std::string(process)=="DecayWithSpin") number++;
  }
  */
  for(int i=0;i<entries;i++){
    tree->GetEntry(i);
    if((std::string(particle)=="e+")&&(std::string(process)=="initStep")) {
      //std::cout << "step:" << nstep <<"\t" <<"particle:"<< particle<< "\t" << "volume:" << volume <<"\t"<< "X:" << X << "\t" << "Y:" << Y << "\t" << "Z:" << Z << "\t" << "depE:" << kE << std::endl;
      e_number++;
    }else if(std::string(process)=="DecayWithSpin") {
      //std::cout << "step:" << nstep <<"\t" <<"particle:"<< particle<< "\t" << "volume:" << volume <<"\t"<< "X:" << X << "\t" << "Y:" << Y << "\t" << "Z:" << Z << "\t" << "depE:" << kE << std::endl;
      mu_number++;
    }
  }
  /*
  for(int i=0;i<29093;i++){ // 29093=max step number
    for(int n=0;n<entries;n++){
      tree->GetEntry(n);
      if((ntrack==1)&&(nstep==i)&&(std::string(particle)=="mu+")) number++;
    }
    if(number!=0) std::cout << "step:" << i << "\t" << number << std::endl;
    number=0;
  }
  */
  std::cout << "mu+ number:" << mu_number << "\t" << "e+ number:" << e_number << std::endl;
  return dist;
}
