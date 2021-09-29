////////////////////////////////////////////////////////////////////////////
//                 High field simulation for MuSEUM Collaboration 
//
//                    Author: Hideharu Yamauchi 2021/09/18
//
//     All physics constants are referenced to the following nist website
//         https://physics.nist.gov/cuu/Constants/Table/allascii.txt
///////////////////////////////////////////////////////////////////////////

#if !defined(___header_geometry_)
#define ___header_geometry_ 1 

#include "TMath.h"
#include <math.h>

//------------Physics constants--------------
const double pi = TMath::Pi();
const double magnet_center = 1.2; // m
const double muon_life = 2.1969811*1.0e-6; // s
const double v_exp = 4463302776; // from liu's experiment
const double plank_const = 6.626070040e-34;
const double plank_const_divided = 1.054571800e-34; // h-bar
const double a = 7.2973525693e-3; // finestructure
const double mass_ratio = 4.83633169e-3;
const double magnetic_moment_mu = -4.49044830e-26;
const double magnetic_moment_j = -9.2847647043e-24;
const double gfactor_mu = -2.0023318418;
const double gfactor_mu_prime = gfactor_mu*(1-pow(a,2.0)/3+pow(a,2.0)*mass_ratio*0.5);
const double gfactor_e = -2.00231930436256;
const double gfactor_j = gfactor_e*(1-pow(a,2.0)/3+pow(a,2.0)*mass_ratio*0.5+pow(a,3.0)*0.25/pi);

const double B_cons = 1.7;
const double polarization = -1.;
const double X = B_cons*(gfactor_j*magnetic_moment_j + gfactor_mu_prime*magnetic_moment_mu)/(plank_const*v_exp);
const double coefficient_s = sqrt(0.5)*sqrt(1-X/sqrt(1+X*X));
const double coefficient_c = sqrt(0.5)*sqrt(1+X/sqrt(1+X*X));

// the coefficient to change RF magnetic field to frequency, 0.001 is for change Hz to kHz
const double b_12 = 0.001*0.25*(coefficient_s*gfactor_j*magnetic_moment_j + coefficient_c*gfactor_mu_prime*magnetic_moment_mu)/plank_const_divided; // kHz/T 
const double b_34 = 0.001*0.25*(coefficient_s*gfactor_j*magnetic_moment_j - coefficient_c*gfactor_mu_prime*magnetic_moment_mu)/plank_const_divided;

//------------HF cavity--------------
const double cavity_radius = 0.0935; // m
double cavity_power[2] = {8., 10.}; // same with liu exp         
double Q_value[2] = {14000., 20000.};  // same with liu exp
const double cavity_volume = 0.304*pi*0.0935*0.0935; 

//------------Positron counter---------------
const double counter_sizeXY = 240.; // mm
double positron_counterU = 300.; 
double positron_counterD = 340.; 

//------------Gas Chamber-------------------
const double chamber_length = 390.; // mm
const double chamber_diameter = 400.; // inner diameter
const double chamber_thickness = 30.; 
  
#endif