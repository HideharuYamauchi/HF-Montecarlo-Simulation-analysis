{
  TCanvas* c = new TCanvas("c3","c3", 1200, 1200);
  c->SetGrid();
  std::ifstream ifs("./t1.dat", std::ios::in);

  if(ifs.fail()) {
    std::cout << "Failed to open file." << std::endl;
    std::exit(EXIT_FAILURE);
  }  
  
  TGraph* gr = new TGraph();
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
  gr->Draw("ALP");
  gr->GetXaxis()->SetTitle("t1 [/#muSec]");
  gr->GetYaxis()->SetTitle("Center Error [/kHz]");
  c->Update();
  
  TLine* l = new TLine(xmin, , xmax, );
  l->SetLineColor(2);
  l->SetLineWidth(3);
  l->SetLineStyle(1);
  l->Draw();
  
  c->SaveAs("t1.png"); 
}
