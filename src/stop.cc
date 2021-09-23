///////////////////////////////////////////        
//  High field simulation(MuSEUM project)             
//                                               
//  Hideharu Yamauchi 2021/09/19                          
///////////////////////////////////////////                         
#ifndef ___class_muonstopping_
#define ___class_muonstopping_

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

muonstopping::muonstopping(std::string runfile, const char* envfile):position(3),time(0.),tree(nullptr){
  run_num = runfile.substr(runfile.find("run"), 7); // get the present run number
  std::ifstream ifs(runfile);
  if(ifs.fail()){
    std::cout << "Fail to open" << runfile << std::endl;
    std::exit(EXIT_FAILURE);
  }
  tree = new TTree("tree","");
  tree->Branch("time",&time,"time/D");
  tree->Branch("particle",particle,"*particle/C");
  tree->Branch("volume",volume,"*volume/C");
  tree->Branch("process",process,"*process/C");
  tree->Branch("X",&X,"X/D");
  tree->Branch("Y",&Y,"Y/D");
  tree->Branch("Z",&Z,"Z/D");
  tree->Branch("Px",&Px,"Px/D");
  tree->Branch("Py",&Py,"Py/D");
  tree->Branch("Pz",&Pz,"Pz/D");
  tree->Branch("kE",&kE,"kE/D");
  tree->Branch("depE",&depE,"depE/D");
  tree->Branch("track",&track,"track/I");
  tree->Branch("step",&step,"step/I");
  tree->Branch("copyno",&copyno,"copyno/I");
  while(!ifs.eof()){
    ifs >> time >> sparticle >> svolume >> sprocess >> X >> Y >> Z >> Px >> Py >> Pz >> kE >> depE >> track >> step >> copyno;
    particle = new char[sparticle.size() + 1];
    std::strcpy(particle, sparticle.c_str());
    volume = new char[svolume.size() + 1];
    std::strcpy(volume, svolume.c_str());
    process = new char[sprocess.size() + 1];
    std::strcpy(process, sprocess.c_str());
    //std::cout << time<<"\t"<< particle<<"\t" << volume<<"\t" << process<<"\t" << X<<"\t" << Y<<"\t" << Z<<"\t" << Px<<"\t" << Py<<"\t" << Pz<<"\t" << kE<<"\t" << depE<<"\t" << track<<"\t" << step<<"\t" << copyno << std::endl;
    tree->Fill();
  }
  ifs.close();
  //branches = (TObjArray*)tree->GetListOfBranches();
  //entrylist = (TEntryList*)tree->GetEntryList();
  // get the number of entries and branches
  //entries = tree->GetEntries(); // 1251294
  //branch_num = branches->GetEntries(); // 15
  //size = branches->GetSize(); // 16=15+row
}

double muonstopping::GetXYZ(double x, double y, double z){ // mm
  position[0] = x*1.0e-3; // convert mm to m
  position[1] = y*1.0e-3;
  position[2] = z*1.0e-3;
  return sqrt(pow(x, 2.)+pow(y, 2.)+pow(z, 2.)); // mm
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
