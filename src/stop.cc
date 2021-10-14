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
#include "RF.cc"
#include "magnet.cc"
#endif

muonstopping::muonstopping(std::string runfile, const char* envfile, Int_t mode)
  : tree(nullptr),branchtime(nullptr),branchX(nullptr),branchY(nullptr),branchZ(nullptr),branchPx(nullptr),branchPy(nullptr),branchPz(nullptr),branchkE(nullptr),branchdepE(nullptr),branchtrack(nullptr),
    branchstep(nullptr),branchcopyno(nullptr)
{
  RF = new RFfield(mode);
  magnet = new magfield("../data/BRECON_MOM_20200716_6.txt", mode);
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
  tree->SetBranchStatus("*",1);
  entries = tree->GetEntries();
  nbranches = tree->GetListOfBranches()->GetEntriesFast();
}

muonstopping::~muonstopping(void){
  delete RF;
  delete magnet;
  delete tree;
}

void muonstopping::CreateRootFile(void){
  file = TFile::Open(run_num+=".root","RECREATE");
  if(tree->Write()) std::cout << run_num << " is made." << std::endl;   
  file->Close();
  delete file;
}

TH2D* muonstopping::Vis_Stopping_DistXY(Double_t zpoint1, Double_t zpoint2, bool saveflag=true){
  TString title = std::to_string(int(zpoint1)) + "~" + std::to_string(int(zpoint2));
  TCanvas* c = new TCanvas("c","c",900,900);
  TPad* center_pad = new TPad("center_pad","",0.0,0.0,0.5,0.5);
  center_pad->Draw();
  TPad* right_pad = new TPad("right_pad","",0.5,0.0,1.0,0.5);
  right_pad->Draw();
  TPad* top_pad = new TPad("top_pad","",0.0,0.5,0.5,1.0);
  top_pad->Draw();
  TH2D* dtxy = new TH2D("XY-Dist", "", 240, -120, 120, 240, -120, 120);
  center_pad->cd();
  dtxy->SetXTitle("X [/mm]");
  dtxy->SetYTitle("Y [/mm]");                         
  center_pad->SetLeftMargin(0.15);
  center_pad->SetRightMargin(-0.03);
  for(int n=0;n<entries;n++){
    tree->GetEntry(n);
    if(std::string(process)=="DecayWithSpin"&&(zpoint1<=Z&&Z<=zpoint2)) dtxy->Fill(X,Y); // vis all z components
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
  if(saveflag){
    c->SaveAs("../figure/"+title+".png");
    delete dtxy;
  }
  delete c;
  return dtxy;
}

TH2D* muonstopping::Vis_Stopping_DistZ(void){
  TCanvas* c2 = new TCanvas("c2", "c2",1400,900);
  gStyle->SetOptStat(0); // do not set the stat tabel 
  TH2D* dtz = new TH2D("Z-Dist", "", 500, 1000, 1400, 240, -120, 120);
  dtz->SetXTitle("Position on the beam axis [/mm]");
  dtz->SetYTitle("Position on the vertical axis [/mm]");
  for(Int_t n=0;n<entries;n++){
    tree->GetEntry(n);
    if(std::string(process)=="DecayWithSpin") dtz->Fill(Z,Y);
  }
  dtz->Draw("Colz");
  dtz->GetXaxis()->SetTitleOffset(1.3);
  dtz->GetYaxis()->SetTitleOffset(1.2);
  c2->SaveAs("../figure/Z-Distribution.png");
  delete c2;
  return dtz;
}

TTree* muonstopping::GetDecayTree(bool scanflag=false){
  TTree* decaytree = new TTree("decaytree","tree of decay muons");
  Double_t decaytime;
  std::string decayvolume;
  std::vector<Double_t> muon_position(3);
  std::vector<Double_t> muon_momentum(3);
  std::vector<Double_t> positron_position(3);
  std::vector<Double_t> positron_momentum(3);
  Double_t positron_energy, muon_energy;
  decaytree->Branch("decaytime",&decaytime,"decaytime/D");
  decaytree->Branch("decayvolume",&decayvolume);
  decaytree->Branch("muon_position",&muon_position);
  decaytree->Branch("muon_momentum",&muon_momentum);
  decaytree->Branch("muon_energy",&muon_energy,"muon_energy/D");
  decaytree->Branch("positron_position",&positron_position);
  decaytree->Branch("positron_momentum",&positron_momentum);
  decaytree->Branch("positron_energy",&positron_energy,"positron_energy/D");
  for(int n=0;n<entries;n++){
    tree->GetEntry(n);
    if((std::string(particle)=="mu+"&&std::string(process)=="DecayWithSpin")
       ||(std::string(particle)=="e+"&&std::string(process)=="initStep")){
      if(std::string(process)=="DecayWithSpin"
	 &&(std::string(volume)=="Cavity"||std::string(volume)=="CavityFoil"||std::string(volume)=="CavityFlange"||std::string(volume)=="TargetGas")){
	decaytime = time;
	decayvolume = std::string(volume);
	muon_position[0] = X;
	muon_position[1] = Y;
	muon_position[2] = Z;
	muon_momentum[0] = Px;
	muon_momentum[1] = Py;
	muon_momentum[2] = Pz;
	muon_energy = kE;
	for(int l=1 ;l<entries-n;l++){
	  tree->GetEntry(n+l);
	  if(std::string(particle)=="e+"&&std::string(process)=="initStep"&&std::string(volume)==decayvolume
	     &&(X==muon_position[0]&&Y==muon_position[1]&&Z==muon_position[2]&&time==decaytime)){
	    positron_position[0] = X;
	    positron_position[1] = Y;
	    positron_position[2] = Z;
	    positron_momentum[0] = Px;
	    positron_momentum[1] = Py;
	    positron_momentum[2] = Pz;
	    positron_energy = kE; // keV
	    break;
	  }
	}
	decaytree->Fill();
      }
    }
  }
  if(scanflag) decaytree->Scan("*");
  return decaytree;
}

void muonstopping::Vis_RFPowerDist(void){
  RF->GetEffectivePower(Vis_Stopping_DistXY(cavity_upfoil_center,
  					    cavity_downfoil_center,
					    false)
			);
}

void muonstopping::Vis_FieldDist(void){
  TCanvas* c = new TCanvas("c", "c",900,900);
  TH1D* hist = new TH1D("hist","",1000,-0.001,0.001);
  hist->SetXTitle("B^{ER} [/Gauss]");
  hist->SetYTitle("");
  for(int k=0; k<entries; k++){
    tree->GetEntry(k);
    if(std::string(process)=="DecayWithSpin"
       &&(std::string(volume)=="Cavity"||std::string(volume)=="CavityFoil"||std::string(volume)=="CavityFlange"||std::string(volume)=="TargetGas")){
      magnet->GetDistance(X, Y, Z-cavity_center);
      hist->Fill(magnet->GetBfieldValue()*1.0e+4);
    }
  }
  hist->Draw();
  c->SaveAs("../figure/magnetdist.png");
  Int_t mean = hist->GetMean();
  Int_t stddev = hist->GetStdDev();                                                                                                                                              
  Int_t RMS = hist->GetRMS();                                                                        
  Int_t entries = hist->GetEntries();
  delete hist;
  delete c;
}
