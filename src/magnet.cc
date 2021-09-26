/////////////////////////////////////////////////////////          
//      Class for calculation of Magnet field     
//                                                  
//  Author:   Hideharu Yamauchi 2021/09/16        
/////////////////////////////////////////////////////////
#ifndef ___class_magfield_
#define ___class_magfield_ 1

#include <stdio.h>
#include <string>
#include "../include/magnet.hh"
#include "../include/HFgeometry.hh"
#include <gsl/gsl_sf_bessel.h>
#include "Math/IFunction.h"
#include <Riostream.h>
#include "TMath.h"
#include "TF1.h"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <math.h>
#include <vector>
#include <numeric>
#endif

magfield::magfield(const char* magnetfile):moment_coordinate(moment_num,std::vector<double>(3)),moment(moment_num,std::vector<double>(3)),distance(moment_num,std::vector<double>(3)),interval(moment_num),position(3){
  std::ifstream ifs(magnetfile);
  if(ifs.fail()) {
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
  }else{std::cout << "distance is over the DSV(30 mm)..." << std::endl;}
  return sqrt(pow(x, 2.)+pow(y, 2.)); // mm
}

double magfield::Bfield_value(void){
  double BvalueZ = 0, total_BvalueZ = 0;
  for(int i=0;i<moment_num;i++){
    BvalueZ = (1.0e-7)*(1/pow(interval[i], 3.))*(3*(moment[i][0]*distance[i][0] + moment[i][1]*distance[i][1] + moment[i][2]*distance[i][2])*distance[i][2]/pow(interval[i], 2.)-moment[i][2]);
    total_BvalueZ += BvalueZ;
  }
  return total_BvalueZ;
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
  //dt->SetTitle("z=? [/mm]; X (mm); Y (mm); B^{ER} [/Gauss]");
  dt->SetTitle(title_dt);

  //dt2 = new TH2D("dt", "z=? [/mm]",201,-100.0 , 100.0 ,201, -100.0, 100.0);
  dt2 = new TH2D("dt", title_dt2, 201,-100.0 , 100.0 ,201, -100.0, 100.0);
  dt2->SetXTitle("X [/mm]");
  dt2->SetYTitle("Y [/mm]");
  dt2->SetZTitle("B^{ER} [/Gauss]");

  int i=0;
  for(int X=-100; X<101; X++){
    for(int Y=-100; Y<101; Y++){      
      if(GetDistance(double(X),double(Y),Z) < cavity_radius*1.0e+3){
	dt->SetPoint(i, double(X), double(Y), Bfield_value()*1.0e+4); // convert T to Gauss
        dt2->Fill(X, Y, Bfield_value()*1.0e+4); // convert T to Gauss
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
