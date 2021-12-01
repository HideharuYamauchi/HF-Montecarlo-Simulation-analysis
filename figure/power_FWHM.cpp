#include <string>
#include <stdio.h>
void power_FWHM(){
  TCanvas* c3 = new TCanvas("c3","c3", 1200, 1200);
  c3->SetGrid();

  TMultiGraph* mg  = new TMultiGraph();
  double px1[41];
  double py1[41];
  double px2[41];
  double py2[41];
  double px3[41];
  double py3[41];
  double px4[41];
  double py4[41];
  double px5[41];
  double py5[41];
  double px6[41];
  double py6[41];
  double px[6][41];
  double py[6][41];
  for(int n=0; n<6;n++){
    /*  
	TF1 *fa1 = new TF1("fa1",
	"2*225*225*(1-cos(TMath::Sqrt(4*TMath::Pi()*TMath::Pi()*x*x+4*225*225)*0.001*9))/(4*TMath::Pi()*TMath::Pi()*x*x+4*225*225)",
	-400,400);
    */
    int p = 130+20*n;
    std::string power = std::to_string(p)+'*'+std::to_string(p);
    std::string start = "0.001*9";
    std::string end = "0.001*10";
    Double_t muon_life = 2.1969811*1.0e-6; // s
    Double_t gamma = (1/muon_life)*1.0e-3;
    //double px[5][41];
    //double py[5][41];
    std::string formula =
      power+"*2/(((2*TMath::Pi()*x)**2)+4*"+power+")*(TMath::Exp(-"+std::to_string(gamma)+'*'+start+")*(1-(TMath::Cos((TMath::Sqrt(((2*TMath::Pi()*x)**2)+4*"+power+"))*"+start+")-(TMath::Sqrt(((2*TMath::Pi()*x)**2)+4*"+power+"))/"+std::to_string(gamma)+"*TMath::Sin((TMath::Sqrt(((2*TMath::Pi()*x)**2)+4*"+power+"))*"+start+"))*"+std::to_string(gamma)+'*'+std::to_string(gamma)+"/((((2*TMath::Pi()*x)**2)+4*"+power+")+"+std::to_string(gamma)+'*'+std::to_string(gamma)+"))-TMath::Exp(-"+std::to_string(gamma)+'*'+end+")*(1-(TMath::Cos((TMath::Sqrt(((2*TMath::Pi()*x)**2)+4*"+power+"))*"+end+")-(TMath::Sqrt(((2*TMath::Pi()*x)**2)+4*"+power+"))/"+std::to_string(gamma)+"*TMath::Sin((TMath::Sqrt(((2*TMath::Pi()*x)**2)+4*"+power+"))*"+end+"))*"+std::to_string(gamma)+'*'+std::to_string(gamma)+"/((((2*TMath::Pi()*x)**2)+4*"+power+")+"+std::to_string(gamma)+'*'+std::to_string(gamma)+")))/(TMath::Exp(-"+std::to_string(gamma)+'*'+start+")-TMath::Exp(-"+std::to_string(gamma)+'*'+end+"))";
    TF1 *fa1 = new TF1("fa1",
		       formula.c_str(),
		       -400,400);
    for(int i=0; i<41; i++){
      px[n][i] = -400+20*i;
      py[n][i] = fa1->Eval(px[n][i]);
    }   
  }
  for(int i=0; i<41; i++){
    px1[i] = px[0][i];
    py1[i] = py[0][i];
  }
  for(int i=0; i<41; i++){
    px2[i] = px[1][i];
    py2[i] = py[1][i];
  }
  for(int i=0; i<41; i++){
    px3[i] = px[2][i];
    py3[i] = py[2][i];
  }
  for(int i=0; i<41; i++){
    px4[i] = px[3][i];
    py4[i] = py[3][i];
  }
  for(int i=0; i<41; i++){
    px5[i] = px[4][i];
    py5[i] = py[4][i];
  }
  for(int i=0; i<41; i++){
    px6[i] = px[5][i];
    py6[i] = py[5][i];
  }
      //fa1->Draw();
      /*
	Double_t HM_left = fa1->GetX(0.5*fa1->GetMaximum(), -100, 0.);
	std::cout << "HM_left:" << HM_left << std::endl;
	Double_t HM_right = fa1->GetX(0.5*fa1->GetMaximum(), 0., 100);
	std::cout << "HM_right:" << HM_right << std::endl;
	std::cout << "FWHM:" << HM_left-HM_right << std::endl;
      */
  TGraph* gr1 = new TGraph(41, px1, py1);
  gr1->SetTitle("b=130 [/kHz]");
  gr1->SetLineColor(2);
  gr1->SetLineWidth(2);
  gr1->SetFillStyle(0);
  //gr1->Draw("AC");
  TGraph* gr2 = new TGraph(41, px2, py2);
  gr2->SetTitle("b=150 [/kHz]");
  gr2->SetLineColor(1);
  gr2->SetLineWidth(2);
  gr2->SetFillStyle(0);
  //gr2->Draw("AC");
  TGraph* gr3 = new TGraph(41, px3, py3);
  gr3->SetTitle("b=170 [/kHz]");
  gr3->SetLineColor(3);
  gr3->SetLineWidth(2);
  gr3->SetFillStyle(0);
  //gr3->Draw("AC");
  TGraph* gr4 = new TGraph(41, px4, py4);
  gr4->SetTitle("b=190 [/kHz]");
  gr4->SetLineColor(4);
  gr4->SetLineWidth(2);
  gr4->SetFillStyle(0);
  //gr4->Draw("AC");
  TGraph* gr5 = new TGraph(41, px5, py5);
  gr5->SetTitle("b=210 [/kHz]");
  gr5->SetLineColor(5);
  gr5->SetLineWidth(2);
  gr5->SetFillStyle(0);
  //gr5->Draw("AC");
  TGraph* gr6 = new TGraph(41, px6, py6);
  gr6->SetTitle("b=230 [/kHz]");
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
  mg->GetXaxis()->SetTitle("Frequency Detuning [/kHz]");
  mg->GetYaxis()->SetTitle("Signal");
  c3->BuildLegend(0.65,0.70,0.80,0.85);
  c3->SaveAs("power_FWHM.png");
}
