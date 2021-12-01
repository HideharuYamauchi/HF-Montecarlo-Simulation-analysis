{
  TCanvas* c3 = new TCanvas("c3","c3", 1200, 1200);
  c3->SetGrid();
  TMultiGraph* mg  = new TMultiGraph();
  const Int_t n = 101;
  Double_t x[n]={0}, y[n];
  Double_t y2[n];
  Double_t time = 12*1.0e-3;
  Double_t power =125;
  Double_t detune = 100;
  Double_t Gamma = sqrt(pow(2*TMath::Pi()*detune,2.)+pow(2*power,2.));
  const double muon_life = 2.1969811*1.0e-6;
  const Double_t gamma = (1/muon_life)*1.0e-3;
  for (Int_t i=0;i<n;i++){
    x[i] = time*i/100;
    y[i] = -pow(power*TMath::Sin(0.5*Gamma*x[i])/Gamma,2.)*(pow(TMath::Sin(0.5*detune*x[i]),2.)-pow(TMath::Cos(0.5*detune*x[i]),2.))*TMath::Exp(-gamma*x[i]);
    y2[i] = 2*y[i];
  }
  TGraph* gr1 = new TGraph(n,x,y);
  gr1->SetTitle("#Gammat_{1} =#pi");
  gr1->SetLineColor(2);
  gr1->SetLineWidth(2);
  gr1->SetFillStyle(0);
  TGraph* gr2 = new TGraph(n,x,y2);
  gr2->SetTitle("#Gammat_{1} =3#pi");
  gr2->SetLineColor(3);
  gr2->SetLineWidth(2);
  gr2->SetFillStyle(0);

  //mg->Add(gr1);
  //mg->Add(gr2);
  gr1->Draw("AC");
  mg->GetXaxis()->SetTitle("");
  mg->GetYaxis()->SetTitle("Probability");
  c3->BuildLegend(0.65,0.70,0.80,0.85);
  c3->SaveAs("probability.png");
}
