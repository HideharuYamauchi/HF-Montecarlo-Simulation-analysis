///////////////////////////////////////////        
//  High field simulation(MuSEUM project)             
//                                               
//  Hideharu Yamauchi 2021/09/19                          
///////////////////////////////////////////                         
#ifndef ___class_muonstopping_
#define ___class_muonstopping_ 1

#include <stdio.h>
#include <ostream>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include "TString.h"
#include "TMath.h"
#include "TEntryList.h"
#include "../include/stop.hh"
#include "../include/HFgeometry.hh"
#include "TCanvas.h"
#include "TStyle.h"
#include "TObject.h"
#include <cstring>
#endif

muonstopping::muonstopping(std::string runfile, const char* envfile):tree(nullptr),volume(nullptr),process(nullptr),branchtime(nullptr),branchX(nullptr),branchY(nullptr),branchZ(nullptr),branchPx(nullptr),branchPy(nullptr),branchPz(nullptr),branchkE(nullptr),branchdepE(nullptr),branchtrack(nullptr),branchstep(nullptr),branchcopyno(nullptr){
    TString runfile2 = runfile;
    run_num = runfile.substr(runfile.find("run"), 7); // get the present run number
    tree = new TTree("tree",run_num);
    tree->ReadFile(runfile2.Data(),"time/D:particle/C:volume:process:X/D:Y:Z:Px:Py:Pz:kE:depE:track/I:step:copyno");
    tree->SetBranchAddress("time",&time, branchtime);
    //tree->SetBranchAddress("particle",&particle, branchparticle);
    tree->SetBranchAddress("X",&X, branchX);
    tree->SetBranchAddress("Y",&Y, branchY);
    tree->SetBranchAddress("Z",&Z, branchZ);
    tree->SetBranchAddress("Px",&Px, branchPx);
    tree->SetBranchAddress("Py",&Py, branchPy);
    tree->SetBranchAddress("Pz",&Pz, branchPz);
    tree->SetBranchAddress("kE",&kE, branchkE);
    tree->SetBranchAddress("depE",&depE, branchdepE);
    tree->SetBranchAddress("track",&ntrack, branchtrack);
    tree->SetBranchAddress("step",&step, branchstep);
    tree->SetBranchAddress("copyno",&copyno, branchcopyno);
    
    entries = tree->GetEntries();
    nbranches = tree->GetListOfBranches()->GetEntriesFast();                              
    //for(int i=0;i<4;i++){tree->GetEntry(i); std::cout << *particle << std::endl;}
}

void muonstopping::MakeRootFile(void){
  file = TFile::Open(run_num+=".root","RECREATE");
  if(tree->Write()) std::cout << run_num << " is made." << std::endl;             
  file->Close();
  delete file;
}

void muonstopping::Vis_particle_track(Int_t track){
  char title_dt[65];
  sprintf(title_dt,"Tracking number%d; X axis [/m]; Y axis [/m]; Z axis [/m]", track);
  char title_dt2[35];
  sprintf(title_dt2,"tracking number%d_3D.png", track);
  TCanvas* c3 = new TCanvas("c2","c2",900,900); 
  TGraph2D *dt = new TGraph2D();
  dt->SetTitle(title_dt);
  for(Int_t n=0;n<entries;n++){
    tree->GetEntry(n);
    if(ntrack==track) dt->SetPoint(n,X*1.0e-3,Y*1.0e-3,Z*1.0e-3);
  }
  gStyle->SetPalette(1);
  dt->Draw("PO");                                             
  dt->GetXaxis()->SetTitleOffset(2.0);
  dt->GetYaxis()->SetTitleOffset(2.0);
  dt->GetZaxis()->SetTitleOffset(1.0);
  c3->SaveAs(title_dt2);
}

void muonstopping::Vis_stopping_distXY(double Z){
  char title_dt[35];
  sprintf(title_dt,"XY-distribution(z=%.1f [/mm])", Z);
  char title_dt2[35];
  sprintf(title_dt2,"XY_distribution(z=%.1f [/mm]).png", Z);
  
  c = new TCanvas("c", "c",1600,600);
  gStyle->SetOptStat(0);
  gStyle->SetTitleXOffset(1.5);
  gStyle->SetTitleYOffset(2);
  
  dtxy = new TH2D("dtxy", title_dt, 201,-100.0 , 100.0 ,201, -100.0, 100.0);

  //c->SaveAs(title_dt2);
  delete dtxy;
  delete c;
}

void muonstopping::Vis_stopping_distZ(void){
  c2 = new TCanvas("c2", "c2",1600,600);
  gStyle->SetOptStat(0);
  gStyle->SetTitleXOffset(1.5);
  gStyle->SetTitleYOffset(2);
  
  dtz = new TH2D("dtz", "Z-distribution", 201, -100.0, 100.0, 201, -100.0, 100.0);

  //c2->SaveAs("Z-distribution.png");
  delete dtz;
  delete c2;
}
