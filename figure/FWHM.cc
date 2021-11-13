//#include <stdio.h>
//#include <fstream>
{
  char method;
  do{
    std::cout << "Conventional[c/C] or OldMuonium[o/O]:" << std::endl;
    std::cin >> method;
  }while(method!='c'&&method!='C'&&method!='o'&&method!='O');

  std::string filename;  
  if(method=='c'||method=='C') filename = "FWHM_conv.dat";
  else if (method=='o'||method=='O') filename = "FWHM_old.dat";
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

  Double_t sim;
  Double_t the;
  Int_t power;
  
  for (int i=0; i<16; i++){
    power = 100+10*i;
    ifs >> sim >> the;
    /*
    std::cout << "sim:" << sim << "\t"
    << "the:" << the << std::endl;
    */
    gr1->SetPoint(i, power, sim);
    gr2->SetPoint(i, power, the);
   }
  
  mg->Add(gr1);
  mg->Add(gr2);
  mg->Draw("ALP");
  mg->GetXaxis()->SetTitle("Microwave Power [/kHz]");
  mg->GetYaxis()->SetTitle("FWHM [/kHz]");
  
  if(method=='c'||method=='C') c3->BuildLegend(0.15,0.65,0.39,0.8);
  else if(method=='o'||method=='O') c3->BuildLegend(0.6,0.7,0.84,0.85);
  
  if(method=='c'||method=='C') c3->SaveAs("FWHM_conv.png");
  else if(method=='o'||method=='O') c3->SaveAs("FWHM_old.png");
}
