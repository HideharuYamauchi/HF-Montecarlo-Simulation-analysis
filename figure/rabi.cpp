{
  TCanvas* c3 = new TCanvas("c3","c3", 1200, 1200);
  c3->SetGrid();
  TMultiGraph* mg  = new TMultiGraph("mg","b=125 [kHz]");
  const Int_t n = 101;
  Double_t x[n]={0}, y[n];
  Double_t x_axis[n]={0};
  Double_t y2[n], y3[n], y4[n], y5[n], y6[n];
  Double_t time = 12*1.0e-3;
  Double_t time_for_axis = 12.;
  Double_t power = 125;
  Double_t detune = 0; //400
  Double_t detune34 = 200; //0
  Double_t detune56 = 400; //200
  Double_t Gamma = sqrt(pow(2*TMath::Pi()*detune,2.)+pow(2*power,2.));
  Double_t Gamma34 = sqrt(pow(2*TMath::Pi()*detune34,2.)+pow(2*power,2.));
  Double_t Gamma56 = sqrt(pow(2*TMath::Pi()*detune56,2.)+pow(2*power,2.));
  const double muon_life = 2.1969811*1.0e-6;
  const Double_t gamma = (1/muon_life)*1.0e-3;
  for (Int_t i=0;i<n;i++){
    x[i] = time*i/100;
    x_axis[i] = time_for_axis*i/100;
    y[i] = 2*pow(power*TMath::Sin(0.5*Gamma*x[i])/Gamma,2.)*TMath::Exp(-gamma*x[i]);
    y2[i] = 0.5*(pow(TMath::Cos(0.5*Gamma*x[i]),2.)+pow(detune*TMath::Sin(0.5*Gamma*x[i])/Gamma,2.))*TMath::Exp(-gamma*x[i]);
    y3[i] = 2*pow(power*TMath::Sin(0.5*Gamma34*x[i])/Gamma34,2.)*TMath::Exp(-gamma*x[i]);
    y4[i] = 0.5*(pow(TMath::Cos(0.5*Gamma34*x[i]),2.)+pow(detune*TMath::Sin(0.5*Gamma34*x[i])/Gamma34,2.))*TMath::Exp(-gamma*x[i]);
    y5[i] = 2*pow(power*TMath::Sin(0.5*Gamma56*x[i])/Gamma56,2.)*TMath::Exp(-gamma*x[i]);
    y6[i] = 0.5*(pow(TMath::Cos(0.5*Gamma56*x[i]),2.)+pow(detune*TMath::Sin(0.5*Gamma56*x[i])/Gamma34,2.))*TMath::Exp(-gamma*x[i]);
  }
  TGraph* gr1 = new TGraph(n,x_axis,y);
  gr1->SetTitle("P_{1} : #Delta#omega=0 [kHz] : #Gamma=250 [kHz]"); //SetTextColor(kRed);
  gr1->SetLineColor(2);
  gr1->SetLineWidth(2);
  gr1->SetFillStyle(0);
  TGraph* gr2 = new TGraph(n,x_axis,y2);
  gr2->SetTitle("P_{2} : #Delta#omega=0 [kHz] : #Gamma=250 [kHz]");
  gr2->SetLineColor(4);
  gr2->SetLineWidth(2);
  gr2->SetFillStyle(0);
  TGraph* gr3 = new TGraph(n,x_axis,y3);
  gr3->SetTitle("P_{1} : #Delta#omega=200 [kHz] : #Gamma=320 [kHz]");
  gr3->SetLineColor(3);
  gr3->SetLineWidth(2);
  gr3->SetFillStyle(0);
  TGraph* gr4 = new TGraph(n,x_axis,y4);
  gr4->SetTitle("P_{2} : #Delta#omega=200 [kHz] : #Gamma=320 [kHz]");
  gr4->SetLineColor(5);
  gr4->SetLineWidth(2);
  gr4->SetFillStyle(0);
  TGraph* gr5 = new TGraph(n,x_axis,y5);
  gr5->SetTitle("P_{1} : #Delta#omega=400 [kHz] : #Gamma=472 [kHz]");
  gr5->SetLineColor(7);
  gr5->SetLineWidth(2);
  gr5->SetFillStyle(0);
  TGraph* gr6 = new TGraph(n,x_axis,y6);
  gr6->SetTitle("P_{2} : #Delta#omega=400 [kHz] : #Gamma=472 [kHz]");
  gr6->SetLineColor(6);
  gr6->SetLineWidth(2);
  gr6->SetFillStyle(0);
 
  mg->Add(gr1);
  mg->Add(gr2);
  mg->Add(gr3);
  mg->Add(gr4);
  mg->Add(gr5);
  mg->Add(gr6);
  mg->Draw("AC");
  //mg->GetXaxis()->SetTitle("Time #times 10^{-6} [/#muSec]");
  mg->GetXaxis()->SetTitle("Time [/#muSec]");
  mg->GetYaxis()->SetTitle("Probability");
  c3->BuildLegend(0.33,0.64,0.88,0.86);
  c3->SaveAs("probability_detune=400.png");
}
