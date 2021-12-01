//#include <stdio.h>
//#include <fstream>
{
  
  std::string filename = "./merit.dat";
  std::ifstream ifs(filename, std::ios::in);
  
  if(ifs.fail()) {
    std::cout << "Failed to open file." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  TCanvas* c3 = new TCanvas("c3","c3", 1200, 1200);
  c3->SetGrid();
  TMultiGraph* mg  = new TMultiGraph();
  
  TGraph* gr1 = new TGraph();
  //gr1->SetTitle("Old Muonium");
  gr1->SetMarkerStyle(21);
  gr1->SetDrawOption("AP");
  gr1->SetLineColor(2);
  
  TGraph* gr2 = new TGraph();
  gr2->SetTitle("Conventional");
  gr2->SetMarkerStyle(0);
  gr2->SetMarkerColor(5);
  gr2->SetDrawOption("AP");
  gr2->SetLineColor(5);
  gr2->SetLineWidth(3);
  
  Double_t t;
  Double_t sim_wid_old;
  Double_t sim_hei_old;
  Double_t num_old;
  Double_t mer_sim_old;
  Double_t sim_wid_con = 165.475;
  Double_t sim_hei_con = 0.116687;
  Double_t num_con = 7.87392e+08;
  Double_t mer_sim_con;
   
  for (int i=0; i<16; i++){
    ifs >> t >> sim_wid_old >> sim_hei_old >> num_old;
    mer_sim_old = num_old*pow((1+sim_hei_old)/sim_wid_old,2.);
    gr1->SetPoint(i, t, mer_sim_old);
   }
  
  mer_sim_con = num_con*pow((1+sim_hei_old)/sim_wid_old,2.);
  gr2->SetPoint(0, 1, mer_sim_con);
  gr2->SetPoint(1, 12, mer_sim_con);
  /*
  mg->Add(gr1);
  mg->Add(gr2);
  mg->Draw("ALP");
  mg->GetYaxis()->SetTitle("Merit [arbitrary unit]");
  mg->GetXaxis()->SetTitle("t1 [/#muSec]");
  c3->BuildLegend(0.45,0.2,0.69,0.35);
  */
  gr1->Draw("ALP");
  gr1->GetYaxis()->SetTitle("Merit [arbitrary unit]");
  gr1->GetXaxis()->SetTitle("t1 [/#muSec]");
  c3->SaveAs("merit.png");
}
