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
#include "../include/stop.hh"
#include "TStyle.h"
#endif

muonstopping::muonstopping(std::string runfile, const char* envfile):tree(nullptr),branchtime(nullptr),branchX(nullptr),branchY(nullptr),branchZ(nullptr),branchPx(nullptr),branchPy(nullptr),branchPz(nullptr),branchkE(nullptr),branchdepE(nullptr),branchtrack(nullptr),branchstep(nullptr),branchcopyno(nullptr){
  TString runfile2 = runfile;
  gStyle->SetPalette(1); // set the color plot
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

TH2D* muonstopping::Vis_stopping_distXY(Double_t posZ){
  TString title_dt = "z=" + std::to_string(int(posZ)) + " [/mm]";
  TString title_dt2 = "XY-Distribution_z:" + std::to_string(int(posZ));
  c = new TCanvas("c","c",900,900);
  TPad* center_pad = new TPad("center_pad","",0.0,0.0,0.5,0.5);
  center_pad->Draw();
  TPad* right_pad = new TPad("right_pad","",0.5,0.0,1.0,0.5);
  right_pad->Draw();
  TPad* top_pad = new TPad("top_pad","",0.0,0.5,0.5,1.0);
  top_pad->Draw();
  TH2D* dtxy = new TH2D("XY-Dist", title_dt, 240, -120, 120, 240, -120, 120);
  center_pad->cd();
  dtxy->SetXTitle("X [/mm]");
  dtxy->SetYTitle("Y [/mm]");                         
  center_pad->SetLeftMargin(0.15);
  center_pad->SetRightMargin(-0.03);
  for(int n=0;n<entries;n++){
    tree->GetEntry(n);
    if(std::string(process)=="DecayWithSpin") dtxy->Fill(X,Y); // vis all z components
    //if(std::string(process)=="DecayWithSpin"&&int(posZ)<=Z&&Z<=int(posZ)+1) dtxy->Fill(X,Y); // vis posZ<z<posZ+1
  }
  TH1D* projdtx = dtxy->ProjectionX();
  projdtx->SetTitle("Position on the horizontal axis");
  projdtx->SetXTitle("X [/mm]");
  TH1D* projdty = dtxy->ProjectionY();
  projdty->SetTitle("Position on the vertical axis");
  projdty->SetXTitle("Y [/mm]");

  center_pad->cd();
  dtxy->Draw("Colz");
  dtxy->GetZaxis()->SetTitleOffset(1.3);
  dtxy->SetStats(000001111); // set the stats table
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
  delete c;
  return dtxy;
}

TH2D* muonstopping::Vis_stopping_distZ(void){
  TString title_dt = "Z-Distribution";
  c2 = new TCanvas("c2", "c2",1400,900);
  gStyle->SetOptStat(0); // do not set the stat tabel 
  //TH2D* dtz = new TH2D("Z-Dist", "", 500, 900, 1400, 240, -120, 120);
  TH2D* dtz = new TH2D("Z-Dist", "", 400, -200, 200, 240, -120, 120);
  dtz->SetXTitle("Position on the beam axis [/mm]");
  dtz->SetYTitle("Position on the vertical axis [/mm]");
  double radius = cavity_radius*1.0e+3; // convert 0.0935 m to 93.5 mm
  double foil = cavity_foil_position*0.5*1.0e+3; // convert 0.304 m to 152 mm
  for(Int_t n=0;n<entries;n++){
    tree->GetEntry(n);
    Z = Z -1050.; // change 1050 to cavity_center
    if(std::string(process)=="DecayWithSpin"&&(-radius<=X&&X<=radius)&&(-radius<=Y&&Y<=radius)&&(-70<=Z&&Z<=foil)) dtz->Fill(Z*2.,X*1.5); // 2. and 1.5 is for scaling
    //if(std::string(process)=="DecayWithSpin") dtz->Fill(Z,Y);
  }
  dtz->Draw("Colz");
  dtz->GetXaxis()->SetTitleOffset(1.3);
  dtz->GetYaxis()->SetTitleOffset(1.2);
  c2->SaveAs(title_dt+=".png");
  delete c2;
  return dtz;
}

TTree* muonstopping::GetDecayTree(void){
  TTree* decaytree = new TTree("decaytree","decay muons");
  Double_t decaytime;
  Double_t decaypositionx;
  Double_t decaypositiony;
  Double_t decaypositionz;
  decaytree->Branch("decaytime",&decaytime,"decaytime/D");
  decaytree->Branch("decaypositionx",&decaypositionx,"decaypositionx/D");
  decaytree->Branch("decaypositiony",&decaypositiony,"decaypositiony/D");
  decaytree->Branch("decaypositionz",&decaypositionz,"decaypositionz/D");
  double radius = cavity_radius*1.0e+3; // convert 0.0935 m to 93.5 mm
  double foil = cavity_foil_position*0.5*1.0e+3; // convert 0.304 m to 152 mm
  for(int n=0;n<entries;n++){
    tree->GetEntry(n);
    //if(std::string(process)=="DecayWithSpin"&&(std::string(volume)=="Cavity"||std::string(volume)=="CavityFoil")){
    Z = Z -1050.;
    if(std::string(process)=="DecayWithSpin"&&(-radius<=X&&X<=radius)&&(-radius<=Y&&Y<=radius)&&(-70<=Z&&Z<=foil)){
      decaytime = time;
      decaypositionx = X*1.5; // 1.5 for scaling
      decaypositiony = Y*1.5;
      //decaypositionz = Z-cavity_center;
      decaypositionz = Z*2.; // 2. for scaling
      decaytree->Fill();
    }
  }
  //decaytree->Scan("*");
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
  return decaytree;
}
