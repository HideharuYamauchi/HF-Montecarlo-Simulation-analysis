//#include <stdio.h>
//#include <fstream>
{
  
  std::string filename_conv;
  std::string filename_old;
  filename_conv = "FWHM_conv.dat";
  filename_old = "FWHM_old.dat";
  std::ifstream ifs_conv(filename_conv, std::ios::in);
  std::ifstream ifs_old(filename_old, std::ios::in);  
  
  if(ifs_conv.fail()||ifs_old.fail()) {
    std::cout << "Failed to open file." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  TCanvas* c3 = new TCanvas("c3","c3", 1200, 1200);
  c3->SetGrid();
  TMultiGraph* mg  = new TMultiGraph();
  
  TGraph* gr1 = new TGraph();
  gr1->SetTitle("Conventional Simulation");
  gr1->SetMarkerStyle(21);
  gr1->SetDrawOption("AP");
  gr1->SetLineColor(2);
  
  TGraph* gr2 = new TGraph();
  gr2->SetTitle("Conventional Theoritical");
  gr2->SetMarkerStyle(22);
  gr2->SetMarkerColor(2);
  gr2->SetDrawOption("AP");
  gr2->SetLineColor(3);

  TGraph* gr3 = new TGraph();
  gr3->SetTitle("OldMuionium Simulation");
  gr3->SetMarkerStyle(23);
  gr3->SetDrawOption("AP");
  gr3->SetLineColor(4);

  TGraph* gr4 = new TGraph();
  gr4->SetTitle("OldMuonium Theoritical");
  gr4->SetMarkerStyle(24);
  gr4->SetMarkerColor(2);
  gr4->SetDrawOption("AP");
  gr4->SetLineColor(5);

  Double_t conv_sim;
  Double_t conv_the;
  Double_t old_sim;
  Double_t old_the;
  Int_t power;
  
  for (int i=0; i<16; i++){
    power = 100+10*i;
    ifs_conv >> conv_sim >> conv_the;
    ifs_old >> old_sim >> old_the;
    gr1->SetPoint(i, power, conv_sim);
    gr2->SetPoint(i, power, conv_the);
    gr3->SetPoint(i, power, old_sim);
    gr4->SetPoint(i, power, old_the);
   }
  
  mg->Add(gr1);
  mg->Add(gr2);
  mg->Add(gr3);
  mg->Add(gr4);
  mg->Draw("ALP");
  mg->GetXaxis()->SetTitle("Microwave Power [/kHz]");
  mg->GetYaxis()->SetTitle("FWHM [/kHz]");
  
  c3->BuildLegend(0.15,0.7,0.39,0.85);
  c3->SaveAs("FWHM.png");
}
