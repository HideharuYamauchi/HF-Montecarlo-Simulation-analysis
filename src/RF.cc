////////////////////////////////////////////////////
//      Class for calculation of the RF field     
//             
//      Author:   Hideharu Yamauchi 2021/09/16                            
////////////////////////////////////////////////////
#ifndef ___class_RFfield_
#define ___class_RFfield_ 1

#include <stdio.h>
#include "../include/RF.hh"
#include "../include/HFgeometry.hh"
#include <gsl/gsl_sf_bessel.h> 
#include "Math/IFunction.h"
#include <Riostream.h>
#include "TMath.h"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <math.h>

RFfield::RFfield(int Mode){
  R__LOAD_LIBRARY(libMathMore);
  mode = Mode;
}

double RFfield::GetXYZ(int x, int y, int z){
  position[0]=x;
  position[1]=y;
  position[2]=z;
  angle = std::atan2(position[1], position[0]);
  distance = sqrt(pow(position[0], 2.0) + pow(position[1], 2.0))*1.0e-3; // convert mm to m
  return distance*1.0e+3; // convert m to mm
}

double RFfield::TM_mode(void){
  if(mode==110){
    kc = j_11/cavity_radius;
    Kr_freq = sqrt(kc*kc/(e*permeability)); // angluar frequency of Kr
    H_coefficient = cavity_power[0]*Q_value[0]/(2*Kr_freq*permeability*cavity_volume*pow((gsl_sf_bessel_Jn(2,j_11)), 2.0));
    Bfield=H_coefficient*(pow(gsl_sf_bessel_Jn(2,kc*distance), 2.0)+pow(gsl_sf_bessel_J0(kc*distance), 2.0)-2*(gsl_sf_bessel_Jn(2,kc*distance))*(gsl_sf_bessel_J0(kc*distance))*std::cos(2*angle));
  }else if(mode==210){
    kc = j_21/cavity_radius; 
    Kr_freq =sqrt(kc*kc/(e*permeability)); 
    H_coefficient = cavity_power[1]*Q_value[1]/(2*Kr_freq*permeability*cavity_volume*pow((gsl_sf_bessel_Jn(3,j_21)), 2.0));
    Bfield=H_coefficient*(pow(gsl_sf_bessel_Jn(3,kc*distance), 2.0)+pow(gsl_sf_bessel_J1(kc*distance), 2.0)-2*(gsl_sf_bessel_Jn(3,kc*distance))*(gsl_sf_bessel_J1(kc*distance))*std::cos(4*angle));
    //std::cout << Bfield << std::endl;
  }
  return permeability*sqrt(Bfield);
}

void RFfield::Vis_RF(){
  c = new TCanvas("c", "c",1600,600);
  gStyle->SetOptStat(0);
  gStyle->SetTitleXOffset(1.5);
  gStyle->SetTitleYOffset(2);
  center_pad = new TPad("center_pad", "center_pad",0.5,0,1.0,1.0);  
  center_pad->Draw();
  top_pad = new TPad("top_pad", "top_pad",0,0,0.5,1.0);
  top_pad->Draw();
  dt = new TGraph2D();
  if(mode==110) {
    dt->SetTitle("TM110_1; X [/mm]; Y [/mm]; B [/T]");
    dt2 = new TH2D("dt", "TM110_2",201,-100.0 , 100.0 ,201, -100.0, 100.0);
  }else if(mode==210) {
    dt->SetTitle("TM210_1; X [/mm]; Y [/mm]; B [/T]");
    dt2 = new TH2D("dt", "TM210_2",201,-100.0 , 100.0 ,201, -100.0, 100.0);
  }
  dt->GetXaxis()->SetTitleOffset(2.0);
  dt2->SetXTitle("X [/mm]");
  dt2->SetYTitle("Y [/mm]");
  dt2->SetZTitle("B [/T]");

  int i=0;                                                                                                                                                                                                
  for(int yy=-100; yy<101; yy++){                                                                                                                                                                         
    for(int xx=-100; xx<101; xx++){
      if(GetXYZ(xx,yy,0)<cavity_radius*1.0e+3){
        dt->SetPoint(i,double(xx),double(yy),TM_mode());	
        dt2->Fill(xx,yy,TM_mode());                                                                                                                                                               
        i++;                                                                                                                                                                                              
      }                                                                                                                                                                                                   
    }                                                                                                                                                                                                     
  }
  //std::cout << H_coefficient << std::endl;
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
  if(mode==110) {c->SaveAs("TM110.png");}
  else if(mode==210) {c->SaveAs("TM210.png");}

  delete dt;
  delete dt2;
  delete center_pad;
  delete top_pad;
  delete c;
}
#endif