//#include <stdio.h>
//#include <fstream>
{
  
  std::string filename = "./merit_130.dat";
  std::ifstream ifs130(filename, std::ios::in);
  std::ifstream ifs150("./merit_150.dat", std::ios::in);
  std::ifstream ifs170("./merit_170.dat", std::ios::in);
  std::ifstream ifs190("./merit_190.dat", std::ios::in);
  std::ifstream ifs210("./merit_210.dat", std::ios::in);
  std::ifstream ifs230("./merit_230.dat", std::ios::in);
  std::ifstream ifs250("./merit_250.dat", std::ios::in);
  std::ifstream ifs270("./merit_270.dat", std::ios::in);
  std::ifstream ifs290("./merit_290.dat", std::ios::in);

  TCanvas* c3 = new TCanvas("c3","c3", 1200, 1200);
  c3->SetGrid();
  TMultiGraph* mg = new TMultiGraph();
  
  TGraph* gr130 = new TGraph();
  //gr130->SetTitle("b=130 [kHz],");
  gr130->SetMarkerStyle(21);
  gr130->SetMarkerColor(2);
  gr130->SetDrawOption("AP");
  gr130->SetLineColor(3);
  gr130->SetLineWidth(3);
  
  TGraph* gr150 = new TGraph();
  //gr150->SetTitle("b=150 [kHz]");
  gr150->SetMarkerStyle(22);
  gr150->SetMarkerColor(4);
  gr150->SetDrawOption("AP");
  gr150->SetLineColor(5);
  gr150->SetLineWidth(3);

  TGraph* gr170 = new TGraph();
  //gr170->SetTitle("b=170 [kHz]");
  gr170->SetMarkerStyle(23);
  gr170->SetMarkerColor(6);
  gr170->SetDrawOption("AP");
  gr170->SetLineColor(7);
  gr170->SetLineWidth(3);

  TGraph* gr190 = new TGraph();
  //gr190->SetTitle("b=190 [kHz] ");
  gr190->SetMarkerStyle(20);
  gr190->SetMarkerColor(2);
  gr190->SetDrawOption("AP");
  gr190->SetLineColor(9);
  gr190->SetLineWidth(3);

  TGraph* gr210 = new TGraph();
  //gr210->SetTitle("b=210 [kHz]");
  gr210->SetMarkerStyle(21);
  gr210->SetMarkerColor(4);
  gr210->SetDrawOption("AP");
  gr210->SetLineColor(11);
  gr210->SetLineWidth(3);

  TGraph* gr230 = new TGraph();
  //gr230->SetTitle("b=230 [kHz]");
  gr230->SetMarkerStyle(22);
  gr230->SetMarkerColor(7);
  gr230->SetDrawOption("AP");
  gr230->SetLineColor(13);
  gr230->SetLineWidth(3);

  TGraph* gr250 = new TGraph();
  //gr250->SetTitle("b=250 [kHz]");
  gr250->SetMarkerStyle(23);
  gr250->SetMarkerColor(3);
  gr250->SetDrawOption("AP");
  gr250->SetLineColor(1);
  gr250->SetLineWidth(3);

  TGraph* gr270 = new TGraph();
  //gr270->SetTitle("b=270 [kHz]");
  gr270->SetMarkerStyle(24);
  gr270->SetMarkerColor(1);
  gr270->SetDrawOption("AP");
  gr270->SetLineColor(5);
  gr270->SetLineWidth(3);

  TGraph* gr290 = new TGraph();
  //gr290->SetTitle("b=290 [kHz]");
  gr290->SetMarkerStyle(25);
  gr290->SetMarkerColor(2);
  gr290->SetDrawOption("AP");
  gr290->SetLineColor(6);
  gr290->SetLineWidth(3);
  
  Double_t t;
  Double_t sim_wid_old;
  Double_t sim_hei_old;
  Double_t num_old;
  Double_t mer_sim_old;
  Double_t sim_wid_con = 165.475;
  Double_t sim_hei_con = 0.116687;
  Double_t num_con = 7.87392e+08;
  Double_t mer_sim_con;
   
  for (int i=0; i<10; i++){
    ifs130 >> t >> sim_wid_old >> sim_hei_old >> num_old;
    mer_sim_old = num_old*pow((1+sim_hei_old)/sim_wid_old,2.);
    gr130->SetPoint(i, t, mer_sim_old);
    
    ifs150 >> t >> sim_wid_old >> sim_hei_old >> num_old;
    mer_sim_old = num_old*pow((1+sim_hei_old)/sim_wid_old,2.);
    gr150->SetPoint(i, t, mer_sim_old);
    
    ifs170 >> t >> sim_wid_old >> sim_hei_old >> num_old;
    mer_sim_old = num_old*pow((1+sim_hei_old)/sim_wid_old,2.);
    gr170->SetPoint(i, t, mer_sim_old);
    
    ifs190 >> t >> sim_wid_old >> sim_hei_old >> num_old;
    mer_sim_old = num_old*pow((1+sim_hei_old)/sim_wid_old,2.);
    gr190->SetPoint(i, t, mer_sim_old);

    ifs210 >> t >> sim_wid_old >> sim_hei_old >> num_old;
    mer_sim_old = num_old*pow((1+sim_hei_old)/sim_wid_old,2.);
    gr210->SetPoint(i, t, mer_sim_old);

    ifs230 >> t >> sim_wid_old >> sim_hei_old >> num_old;
    mer_sim_old = num_old*pow((1+sim_hei_old)/sim_wid_old,2.);
    gr230->SetPoint(i, t, mer_sim_old);

    ifs250 >> t >> sim_wid_old >> sim_hei_old >> num_old;
    mer_sim_old = num_old*pow((1+sim_hei_old)/sim_wid_old,2.);
    gr250->SetPoint(i, t, mer_sim_old);
   }

  for (int i=0; i<9; i++){
    ifs270 >> t >> sim_wid_old >> sim_hei_old >> num_old;
    mer_sim_old = num_old*pow((1+sim_hei_old)/sim_wid_old,2.);
    gr270->SetPoint(i, t, mer_sim_old);

    ifs290 >> t >> sim_wid_old >> sim_hei_old >> num_old;
    mer_sim_old = num_old*pow((1+sim_hei_old)/sim_wid_old,2.);
    gr290->SetPoint(i, t, mer_sim_old);
  }
  /*
  mer_sim_con = num_con*pow((1+sim_hei_old)/sim_wid_old,2.);
  gr2->SetPoint(0, 1, mer_sim_con);
  gr2->SetPoint(1, 12, mer_sim_con);
  */

  //std::cout << "gr130 Integral:" << gr130->Integral() << std::endl;

  mg->Add(gr130);
  std::string gr130str = "b=130 [kHz], Integral="+std::to_string(gr130->Integral());
  gr130->SetTitle(gr130str.c_str());
  mg->Add(gr150);
  std::string gr150str = "b=150 [kHz], Integral="+std::to_string(gr150->Integral());
  gr150->SetTitle(gr150str.c_str());
  mg->Add(gr170);
  std::string gr170str = "b=170 [kHz], Integral="+std::to_string(gr170->Integral());
  gr170->SetTitle(gr170str.c_str());
  mg->Add(gr190);
  std::string gr190str = "b=190 [kHz], Integral="+std::to_string(gr190->Integral());
  gr190->SetTitle(gr190str.c_str());
  mg->Add(gr210);
  std::string gr210str = "b=210 [kHz], Integral="+std::to_string(gr210->Integral());
  gr210->SetTitle(gr210str.c_str());
  mg->Add(gr230);
  std::string gr230str = "b=230 [kHz], Integral="+std::to_string(gr230->Integral());
  gr230->SetTitle(gr230str.c_str());
  mg->Add(gr250);
  std::string gr250str = "b=250 [kHz], Integral="+std::to_string(gr250->Integral());
  gr250->SetTitle(gr250str.c_str());
  mg->Add(gr270);
  std::string gr270str = "b=270 [kHz], Integral="+std::to_string(gr270->Integral());
  gr270->SetTitle(gr270str.c_str());
  mg->Add(gr290);
  std::string gr290str = "b=290 [kHz], Integral="+std::to_string(gr290->Integral());
  gr290->SetTitle(gr290str.c_str());
  mg->Draw("ALP");
  mg->GetYaxis()->SetTitle("Merit [arbitrary unit]");
  mg->GetXaxis()->SetTitle("t1 [/#muSec]");
  c3->BuildLegend(0.64,0.67,0.88,0.89);
  
  //gr1->Draw("ALP");
  //gr1->GetYaxis()->SetTitle("Merit [arbitrary unit]");
  //gr1->GetXaxis()->SetTitle("t1 [/#muSec]");
  c3->SaveAs("merit.png");
}
