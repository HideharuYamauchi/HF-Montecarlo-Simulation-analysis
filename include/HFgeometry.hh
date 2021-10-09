///////////////////////////////////////////////////////////////////////////////
//                 High field simulation for MuSEUM Collaboration 
//
//                    Author : Hideharu Yamauchi 2021/09/18
//
//     All physics constants are referenced to the following nist website
//         https://physics.nist.gov/cuu/Constants/Table/allascii.txt
//////////////////////////////////////////////////////////////////////////////

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
const double a = 7.2973525693e-3; // fine structure
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

// the coefficient to change RF field to frequency(Hz), 0.001 is for coverting Hz/T to kHz/T
const double b_12 = 0.001*0.25*(coefficient_s*gfactor_j*magnetic_moment_j + coefficient_c*gfactor_mu_prime*magnetic_moment_mu)/plank_const_divided; // kHz/T 
const double b_34 = 0.001*0.25*(coefficient_s*gfactor_j*magnetic_moment_j - coefficient_c*gfactor_mu_prime*magnetic_moment_mu)/plank_const_divided; // kHz/T

// beam setting
const double beam = 27.4; // MeV/c
const double beam_center = 0.; // mm
const double beam_x_sigma = 34.; // mm
const double beam_y_sigma = 30.; // mm

const double muon_mass = 105.6583755; // MeV/c^2
const double positron_max_momentum = 52.83; // MeV/c

// Geant4 setting 
const double Magnet_center = 1200.; // mm
const double G4WorldSizeXY = 1200.; // mm
const double G4WorldSizeZ = 3600.; // mm

//------------Gas Chamber-------------------
const double chamber_length = 390.; // mm
const double chamber_diameter = 400.; // inner diameter
const double chamber_thickness = 30.;
const double chamber_window_diameter = 100.; // mm, front window
const double chamber_foil_thickness = 0.1; // mm
const double chamber_flange_diameter = 430.; // mm
const double chamber_flange_thickness_u = 30.; // mm

const double chamber_center = Magnet_center; // 1200 mm
const double chamber_upflange_center = Magnet_center-(chamber_length+chamber_flange_thickness_u)*0.5; // 1200-(390+30)*0.5= 990 mm
const double chamber_downflange_center = Magnet_center+(chamber_length+chamber_flange_thickness_d)*0.5; // 1200+(390+30)*0.5= 1410 mm
const double chamber_foil_center = Magnet_center-(chamber_length+chamber_foil_thickness)*0.5-chamber_flange_thickness_u; // 1200-(390+0.1)*0.5-30=1034.95 mm

const double overhang_thickness = 25.; // mm
  
//------------HF cavity--------------
const double cavity_radius = 0.0935; // m
const double cavity_radial_thickness = 0.015; //m
double cavity_power[2] = {8., 10.}; // same with liu exp
double Q_value[2] = {14000., 20000.};  // same with liu exp
const double cavity_length = 0.27; // m
const double cavity_flange_length = 0.031; // m, up=down, 0.013+0.009+0.009
const double cavity_foil_position = 0.304; // =304mm
const double cavity_volume = cavity_foil_position*pow(cavity_radius, 2.)*pi; 
const double cavity_foil_thickness = 25*1.0e-6; // m

// G4 setting
const double cavity_center = Magnet_center-overhang_thickness*0.5; // mm, 1200-25*0.5 = 1187.5 mm
const double cavity_upflange_center = cavity_center-(cavity_flange_length+cavity_length)*0.5; // 1187.5-(44+244)*0.5= 1043.5
const double cavity_downflange_center = cavity_center+(cavity_flange_length+cavity_length)*0.5; // 1187.5+(44+244)*0.5= 1331.5
const double cavity_upfoil_center = cavity_center-26.-cavity_length*0.5; // 1187.5-26-244*0.5= 1039.5
const double cavity_downfoil_center = cavity_center+26.+cavity_length*0.5; // 1187.5+26+244*0.5= 1335.5
  
//------------Positron counter---------------
const double counter_sizeXY = 240.; // mm
const double positron_counterU = 300.; 
const double positron_counterD = 340.;
const double counter_couver_thickness = 2.; 

//------Beam window kapton foil;G4 setting--------
const double kapton_diameter = 300.;
const double kapton_thickness = 0.075; // 0.075*0.5 = 0.0375
const double kapton_center = -500.0;
const double kapton_position = Magnet_center+kapton_center; // 700 mm

const double beam_vacuum_region = Magnet_center+kapton_center+(G4WorldSizeZ-kapton_thickness)*0.5; // 1200-500+(3600-0.075)*0.5=700+1800-0.0325=2499.9675 mm
const double beam_vacuum_region_center = (beam_vacuum_region-G4WorldSizeZ)*0.5; // 1249.98375-1800 = -550.01625 mm
const double beam_vacuum_region_end = beam_vacuum_region_center+beam_vacuum_region*0.5; // -550.01625 + 1249.98375 = 699.9675 mm

//------Beam Profile Monitor;G4 setting---------
const double bpm_sizeXY = 100.; // mm
const double bpm_thickness = 0.15;
const double bpm_positionU = kapton_center+5.; // -495 mm
const double bpm_positionD = kapton_center+10.; // -490 mm

const double upbpm_center = Magnet_center+bpm_positionU; // 1200-495 = 705 mm
const double downbpm_center = Magnet_center+bpm_positionD; // 1200-490 = 710 mm

#endif
