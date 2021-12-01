//#include <stdio.h>
//#include <fstream>
{
  
  std::string filename = "./t1_height.dat";
  std::ifstream ifs(filename, std::ios::in);
  
  if(ifs.fail()) {
    std::cout << "Failed to open file." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  TCanvas* c3 = new TCanvas("c3","c3", 1200, 1200);
  c3->SetGrid();
  TMultiGraph* mg  = new TMultiGraph();
  
  TGraph* gr1 = new TGraph();
  gr1->SetTitle("Simulation");
  gr1->SetMarkerStyle(21);
  gr1->SetDrawOption("AP");
  gr1->SetLineColor(2);
  
  TGraph* gr2 = new TGraph();
  gr2->SetTitle("Theoritical");
  gr2->SetMarkerStyle(22);
  gr2->SetMarkerColor(2);
  gr2->SetDrawOption("AP");
  gr2->SetLineColor(3);

  Double_t t;
  Double_t sim_hei;
  Double_t the_hei;
   
  for (int i=0; i<16; i++){
    ifs >> t >> sim_hei >> the_hei;
    gr1->SetPoint(i, t, sim_hei);
    gr2->SetPoint(i, t, the_hei);
   }
  
  mg->Add(gr1);
  mg->Add(gr2);
  mg->Draw("ALP");
  mg->GetYaxis()->SetTitle("Signal_{MAX}");
  mg->GetXaxis()->SetTitle("t1 [/#muSec]");
  
  c3->BuildLegend(0.2,0.65,0.44,0.80);
  c3->SaveAs("t1_Height.png");
}
