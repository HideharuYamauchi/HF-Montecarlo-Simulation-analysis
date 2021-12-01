{
  Double_t conv_center_error = 0.113649;
  TCanvas* c = new TCanvas("c3","c3", 1200, 1200);
  c->SetGrid();
  std::ifstream ifs("./t1.dat", std::ios::in);

  if(ifs.fail()) {
    std::cout << "Failed to open file." << std::endl;
    std::exit(EXIT_FAILURE);
  }  

  TMultiGraph* mg  = new TMultiGraph();
  TGraph* gr = new TGraph();
  gr->SetTitle("Old Muonium");
  //gr->Draw("ALP");
  //double xmax = gr->GetUxmax();
  //double xmin = gr->GetUxmin();
  /*
  std::cout << "GetUxmin:" << c->GetUxmin() << "\n"
	    << "GetUxmax:" << c->GetUxmax() << std::endl;
  */
  int t1;
  double error;
  for(int i=0; i<11; i++){
    ifs >> t1 >> error;
    gr->SetPoint(i, t1, error);
    if(t1==12) break;
  }
  
  gr->SetMarkerColor(4);
  gr->SetMarkerStyle(21);
  c->Update();

  
  TGraph* gr4 = new TGraph();
  gr4->SetTitle("Conventional");
  gr4->SetMarkerStyle(0);
  gr4->SetMarkerColor(5);
  gr4->SetDrawOption("AP");
  gr4->SetLineColor(5);
  gr4->SetLineWidth(3);
  gr4->SetPoint(0, 3, conv_center_error); // (number, x, y)
  gr4->SetPoint(1, 12, conv_center_error);

  mg->Add(gr);
  mg->Add(gr4);
  mg->Draw("ALP");
  mg->GetXaxis()->SetTitle("t1 [/#muSec]");
  mg->GetYaxis()->SetTitle("Center Error [/kHz]");
   
  /*
  TLine* l = new TLine(xmin, , xmax, );
  l->SetLineColor(2);
  l->SetLineWidth(3);
  l->SetLineStyle(1);
  l->Draw();
  */
  c->BuildLegend(0.45,0.65,0.69,0.8);
  c->SaveAs("t1.png"); 
}
