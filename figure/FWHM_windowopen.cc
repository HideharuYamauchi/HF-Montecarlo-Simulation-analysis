//#include <stdio.h>
//#include <fstream>
{
  const double plank_const = 6.626070040e-34;
  const double plank_const_divided = 1.054571800e-34; // h-bar
  
  std::string filename;  
  filename = "FWHM_windowopen.dat";
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
  gr2->SetDrawOption("P");
  gr2->SetLineColor(3);

  TGraph* gr3 = new TGraph();
  gr3->SetTitle("Uncertainty Principle: #DeltaT#DeltaE=#hbar/2");
  gr3->SetMarkerStyle(23);
  gr3->SetMarkerColor(2);
  gr3->SetDrawOption("P");
  gr3->SetLineColor(4);

  TGraph* gr4 = new TGraph();
  gr4->SetTitle("Natural Width: 145 [/kHz]");
  gr4->SetMarkerStyle(0);
  gr4->SetMarkerColor(5);
  gr4->SetDrawOption("AP");
  gr4->SetLineColor(5);
  gr4->SetLineWidth(4);
  gr4->SetPoint(0, 3, 145);
  gr4->SetPoint(1, 12, 145);

  Double_t sim;
  Double_t the;
  Int_t windowopen;
  Double_t energy;
  Double_t freq;
  
  for (int i=0; i<16; i++){
    ifs >> windowopen >> sim >> the;
    energy = 0.5*plank_const_divided*1.0e+6/(windowopen+1); // 0.5*h_bar(Js)/t
    freq = energy/plank_const_divided; // Hz
    /*
    std::cout << "sim:" << sim << "\t"
    << "the:" << the << std::endl;
    */
    gr1->SetPoint(i, windowopen, sim);
    gr2->SetPoint(i, windowopen, the);
    gr3->SetPoint(i, windowopen, 0.001*freq);
   }
  
  mg->Add(gr1);
  mg->Add(gr2);
  mg->Add(gr3);
  mg->Add(gr4);
  mg->Draw("ALP");
  mg->GetXaxis()->SetTitle("t1 [/#muSec]");
  mg->GetYaxis()->SetTitle("FWHM [/kHz]");
  
  c3->BuildLegend(0.55,0.65,0.79,0.8);
  c3->SaveAs("FWHM_windowopen.png");
}
