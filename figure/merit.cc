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
  Double_t sim_wid;
  Double_t the_wid;
  Double_t sim_hei;
  Double_t the_hei;
  Double_t num;
  Double_t mer_sim;
  Double_t mer_the;
   
  for (int i=0; i<16; i++){
    ifs >> t >> sim_wid >> the_wid
	>> sim_hei >> the_hei >> num;
    mer_sim = num*pow((1+sim_hei)/sim_wid,2.);
    mer_the = num*pow((1+the_hei)/the_wid,2.);
    gr1->SetPoint(i, t, mer_sim);
    gr2->SetPoint(i, t, mer_the);
   }
  
  mg->Add(gr1);
  mg->Add(gr2);
  mg->Draw("ALP");
  mg->GetYaxis()->SetTitle("Merit [arbitrary unit]");
  mg->GetXaxis()->SetTitle("t1 [/#muSec]");
  
  c3->BuildLegend(0.45,0.2,0.69,0.35);
  c3->SaveAs("merit.png");
}
