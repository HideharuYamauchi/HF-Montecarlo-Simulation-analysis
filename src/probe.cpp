//山括弧(<>)を使用するとプリプロセッサはシステム内の標準の場所、通常は/usr/includeにあるヘッダーファイルを検索。
#include <stdio.h>
#include <math.h>
#include "TMath.h"
#include <ostream>
#include <iostream>
#include <string>
#include <iomanip>
#define B 1.7
#define reference_temp 25.
double cal_x_error(double h, double delta_B ,double v_exp ,double delta_v_exp ,double gfactor_mu_prime ,double delta_gfactor_mu_prime ,double gfactor_j ,double delta_gfactor_j ,double mag_mom_mu ,double delta_mag_mom_mu ,double mag_mom_j ,double delta_mag_mom_j ,double x){
  double delta_x_H = delta_B*(gfactor_j*mag_mom_j+gfactor_mu_prime*mag_mom_mu)/(h*v_exp);
  double delta_x_v_exp = -delta_v_exp*(x/v_exp);
  double delta_x_gfactor_mu_prime = delta_gfactor_mu_prime*mag_mom_mu*B/(h*v_exp);
  double delta_x_gfactor_j = delta_gfactor_j*mag_mom_j*B/(h*v_exp);
  double delta_x_mag_mom_mu = delta_mag_mom_mu*gfactor_mu_prime*B/(h*v_exp);
  double delta_x_mag_mom_j = delta_mag_mom_j*gfactor_j*B/(h*v_exp); 
  //double total_x_error = sqrt(pow(delta_x_H,2.0)+pow(delta_x_v_exp,2.0)+pow(delta_x_gfactor_mu_prime,2.0)+pow(delta_x_gfactor_j,2.0)+pow(delta_x_mag_mom_mu,2.0)+pow(delta_x_mag_mom_j,2.0));
  //std::cout << "delta_x_H="  << delta_x_H << "\n" << std::endl;
  //std::cout << "total_x_error="  << total_x_error << "\n" << std::endl;
  return delta_x_H;
}

double cal_v12_error(double h, double mag_mom_mu, double gfactor_mu_prime, double v_exp, double x, double delta_B, double delta_v_exp, double delta_x, double delta_mag_mom_mu, double delta_gfactor_mu_prime){
  double v12 = (-mag_mom_mu*gfactor_mu_prime*B)/h+(0.5*v_exp*((1+x)-sqrt(1+std::pow(x,2.0)))); // unit GHz
  double v12_H = -(delta_B*mag_mom_mu*gfactor_mu_prime)/h;
  double v12_mag_mom_mu = -(delta_mag_mom_mu*delta_gfactor_mu_prime*B)/h;
  double v12_gfactor_mu_prime = -(delta_gfactor_mu_prime*delta_mag_mom_mu*B)/h;
  double v12_v_exp = delta_v_exp*0.5*((1+x)-sqrt(1+std::pow(x,2.0)));
  double v12_x = delta_x*0.5*v_exp*(1 - x*std::pow(1+x*x,-0.5));
  double total_v12_error = sqrt(pow(v12_H,2.0)+pow(v12_x,2.0));
  //double total_v12_error = sqrt(pow(v12_H,2.0)+pow(v12_mag_mom_mu,2.0)+pow(v12_gfactor_mu_prime,2.0)+pow(v12_v_exp,2.0)+pow(v12_x,2.0));
  //std::cout << std::setprecision(10) << v12 << std::endl;
  //std::cout << "v12_H=" << v12_H << "\n" << "v12_mag_mom_mu=" << v12_mag_mom_mu << "\n" << "v12_gfactor_mu_prime=" << v12_gfactor_mu_prime << "\n" << "12_v_exp" << v12_v_exp << "\n" << "v12_x=" << v12_x << "\n" << std::endl;
  std::cout << "total_v12_error=" << total_v12_error << "[Hz]" << std::endl;
  return total_v12_error;
}

double cal_v34_error(double h, double mag_mom_mu, double gfactor_mu_prime, double v_exp, double x, double delta_B, double delta_v_exp, double delta_x, double delta_mag_mom_mu, double delta_gfactor_mu_prime){
  double v34 = (mag_mom_mu*gfactor_mu_prime*B)/h+(0.5*v_exp*((1-x)+sqrt(1+std::pow(x,2.0)))); // unit GHz
  double v34_H = (delta_B*mag_mom_mu*gfactor_mu_prime)/h;;
  double v34_mag_mom_mu = delta_gfactor_mu_prime*delta_mag_mom_mu*B/h;
  double v34_gfactor_mu_prime = delta_gfactor_mu_prime*delta_mag_mom_mu*B/h;
  double v34_v_exp = delta_v_exp*0.5*((1-x)+sqrt(1+std::pow(x,2.0)));
  double v34_x = delta_x*0.5*v_exp*(x*std::pow(1+x*x,-0.5) - 1);
  double total_v34_error = sqrt(pow(v34_H,2.0)+pow(v34_x,2.0));
  //double total_v34_error = sqrt(pow(v34_H,2.0)+pow(v34_mag_mom_mu,2.0)+pow(v34_gfactor_mu_prime,2.0)+pow(v34_v_exp,2.0)+pow(v34_x,2.0));
  //std::cout << std::setprecision(10) << v34 << std::endl;
  //std::cout << "v34_H=" << v34_H << "\n" << "v34_mag_mom_mu=" << v34_mag_mom_mu << "\n" << "v34_gfactor_mu_prime=" << v34_gfactor_mu_prime << "\n" << "v34_v_exp=" << v34_v_exp << "\n" << "v34_x=" <<v34_x << std::endl;
  std::cout << "total_v34_error=" << total_v34_error << "[Hz]" << std::endl;
  return total_v34_error;
}

/*
double cal_B_error(double T, double b, double t, double proton_gyromag_ratio, double delta_proton_gyromag_ratio){
  const double pi = TMath::Pi();
  double Delta = ;
  double freq = 0.5*B*delta_proton_gyromag_ratio/pi;
  double freq_proton_gyromag_ratio = -delta_proton_gyromag_ratio*2*pi*freq/pow(proton_gyromag_ratio,2.0);
  double freq_correct = Delta*2*pi/proton_gyromag_ratio;
  double total_B_error = sqrt(pow(freq_proton_gyromag_ratio,2.0)+pow(freq_correct,2.0));
  std::cout << "freq=" << freq << std::endl;
  return total_B_error;
}
*/

int main(){
  ///////////////////////////////////////////////////////////////////////////////////////////////
  //
  //         All physics constants are referenced to the following nist website
  //             (https://physics.nist.gov/cuu/Constants/Table/allascii.txt)
  //
  /////////////////////////////////////////////////////////////////////////////////////////////// 
  const double h = 6.62607015e-34; //プランク定数, exact
  const double h_bar = 1.054571800e-34;
  const double pi = TMath::Pi();
  const double a = 7.2973525693e-3; //0.000 000 0011 e-3
  const double mass_ratio = 4.83633169e-3; //0.000 000 11 e-3
  const double g_factor_mu = -2.0023318418; //0.000 000 0013
  const double g_factor_mu_prime = g_factor_mu*(1-pow(a,2.0)/3+pow(a,2.0)*mass_ratio*0.5);
  const double g_factor_e = -2.00231930436256; //0.000 000 000 000 35
  const double g_factor_j = g_factor_e*(1-pow(a,2.0)/3+pow(a,2.0)*mass_ratio*0.5+pow(a,3.0)*0.25/pi);
  //std::cout << g_factor_mu_prime << "\t" << g_factor_j << std::endl;
  const double magnetic_moment_mu = -4.49044830e-26; // 0.000 000 10 e-26
  const double magnetic_moment_j = -9.2847647043e-24; // 0.000 000 0028 e-24
  const double v_exp = 4463302776; // from liu's experiment, 11 ppb
  const double proton_gyromagnetic_ratio = 42.577478518; // gyromagnetic ratio of the proton, 18 ppb
  const double x = B*(g_factor_j*magnetic_moment_j + g_factor_mu_prime*magnetic_moment_mu)/(h*v_exp);
  const double E1 = 0.25*v_exp+0.5*B*(g_factor_j*magnetic_moment_j + g_factor_mu_prime*magnetic_moment_mu)/h;
  const double E2 = -0.25*v_exp+0.5*v_exp*sqrt(1+x*x);
  const double E3 = 0.25*v_exp-0.5*B*(g_factor_j*magnetic_moment_j + g_factor_mu_prime*magnetic_moment_mu)/h;
  const double E4 = -0.25*v_exp-0.5*v_exp*sqrt(1+x*x);
  std::cout << "E1:" << E1*1.e-9 << " [GHz]" << "\n"
	    << "E2:" <<	E2*1.e-9 << " [GHz]" << "\n"
	    << "E3:" <<	E3*1.e-9 << " [GHz]" << "\n"
	    << "E4:" <<	E4*1.e-9 << " [GHz]" << "\n"
	    << "E1-E2:" << (E1-E2)*1.e-9 << " [GHz]" << "\n"
	    << "E3-E4:" << (E3-E4)*1.e-9 << " [GHz]" << "\n"
	    << "E2-E3:" << (E2-E3)*1.e-9 << " [GHz]" << "\n"
	    << "E4-E1:" << (E1-E4)*1.e-9 << " [GHz]" << "\n"
	    << "MuHFS:" << (E1-E2+E3-E4)*1.e-9 << " [GHz]" << std::endl;
  
  ///////////////////////////////////
  //    set errors for probe
  ///////////////////////////////////
  const double ppm = 1.0e-6;
  const double ppb = 1.0e-9;
  const double ppt = 1.0e-12;

  /*
  const double sigma_water = 25680*ppb; // internal diamagnetic shielding H2O at 25, from Metrologia 51 (2014) 54–60
  const double delta_sigma_water = 2.5*ppb; // from Metrologia 51 (2014) 54–60
  const double delta_T = -10.36*ppb*(reference_temp-25); // from PHYSICAL REVIEW A 103, 042208 (2021)
  //const double delta_delta_T = 0.3*ppb*(reference_temp-25); // from PHYSICAL REVIEW A 103, 042208 (2021)
  const double delta_b = (0.5-1/3)*-9.032*ppb; // from PHYSICAL REVIEW A 103, 042208 (2021) which correspond to seo san's master thesis ,20℃
  const double delta_p = 2*ppb;
  const double delta_t = ;
  const double Delta = delta_sigma_water + delta_T + delta_t + delta_b;
  */
  const double Delta = 15*ppb;
  std::cout << "Delta:" << Delta << std::endl;
  
  //////////////////////////////////////////////////////////////////////
  //         calculate errors of all physics constants
  //////////////////////////////////////////////////////////////////////
  double delta_B = Delta*B;
  double delta_a = a*1.1*ppt;
  double delta_mass_ratio = mass_ratio*0.11*ppb;
  double delta_g_factor_mu = g_factor_mu*1.3*ppb;
  double delta_g_factor_e = g_factor_e*0.35*ppt;
  double delta_magnetic_moment_mu = magnetic_moment_mu*1.0e-33;
  double delta_magnetic_moment_j = magnetic_moment_j*2.8e-33;
  double delta_v_exp = v_exp*11*ppb;
  double delta_g_factor_mu_prime = delta_g_factor_mu*(1-pow(a,2.0)/3+pow(a,2.0)*mass_ratio*0.5) + a*g_factor_mu*(mass_ratio-2/3)*delta_a + 0.5*g_factor_mu*pow(a,2.0)*delta_mass_ratio;
  double delta_g_factor_j = delta_g_factor_e*(1-pow(a,2.0)/3+pow(a,2.0)*mass_ratio*0.5+pow(a,3.0)*0.25/pi) + delta_a*a*g_factor_e*(mass_ratio+0.75*a/pi-2/3) + delta_mass_ratio*0.5*pow(a,2.0)*g_factor_e;
  double delta_proton_gyromagnetic_ratio = 42.577478518*18*ppb;
  
  //double delta_B = cal_B_error(delta_T, delta_b, delta_t, proton_gyromagnetic_ratio, delta_proton_gyromagnetic_ratio);
  double delta_x = cal_x_error(h, delta_B ,v_exp ,delta_v_exp ,g_factor_mu_prime ,delta_g_factor_mu_prime ,g_factor_j ,delta_g_factor_j ,magnetic_moment_mu ,delta_magnetic_moment_mu ,magnetic_moment_j, delta_magnetic_moment_j, x);
  std::cout << delta_x << std::endl;
  double delta_v12 = cal_v12_error(h, magnetic_moment_mu, g_factor_mu_prime, v_exp, x, delta_B, delta_v_exp, delta_x, delta_magnetic_moment_mu, delta_g_factor_mu_prime);
  double delta_v34 = cal_v34_error(h, magnetic_moment_mu, g_factor_mu_prime, v_exp, x, delta_B, delta_v_exp, delta_x, delta_magnetic_moment_mu, delta_g_factor_mu_prime);
  //double delta_v = sqrt(pow(delta_v12,2.0) + pow(delta_v34,2.0));
  //std::cout << "total_vHF_error="  << delta_v << "[Hz]" <<std::endl;
  
  return 0;
}
