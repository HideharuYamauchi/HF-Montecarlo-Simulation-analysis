void twomethod(){

  gStyle->SetOptFit(); // need this line to draw staic box
  
  std::ifstream ifs_con("./reso_conv.dat", std::ios::in);
  std::ifstream ifs_old("./reso_old.dat", std::ios::in);

  if(ifs_con.fail()||ifs_old.fail()) {
    std::cout << "Failed to open file." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  TCanvas* c3 = new TCanvas("c3","c3", 1200, 1200);
  c3->SetGrid();
  
  TMultiGraph* mg  = new TMultiGraph("","");
  TGraphErrors* con_curve = new TGraphErrors();
  TGraphErrors* old_curve = new TGraphErrors();
  
  const double muon_life = 2.1969811*1.0e-6;
  const Double_t gamma = (1/muon_life)*1.0e-3;
  Double_t detuning, Non, Noff, error;

  int n=0;
  for(int i=-40; i<41 ; i++){
    ifs_con >> detuning >> Non >> Noff >> error;
    con_curve->SetPoint(n, detuning, Non/Noff-1);
    con_curve->SetPointError(n, 0, error);
    
    ifs_old >> detuning >> Non >> Noff >> error;
    old_curve->SetPoint(n, detuning, Non/Noff-1);
    old_curve->SetPointError(n, 0, error);
    n++;
  }
  std::string fit_no_con = "[0]+[3]*2*[1]*[1]/(4*TMath::Pi()*TMath::Pi()*(x-[2])*(x-[2])+4*[1]*[1]+"+std::to_string(gamma)+'*'+std::to_string(gamma)+')';
  TF1* f1_conv = new TF1("f1_conv", fit_no_con.c_str(), -400, 400);
  
  f1_conv->FixParameter(0, 0); // offset
  f1_conv->SetParameter(1, 124.); // microwave power
  f1_conv->SetParameter(2, 0); // center
  f1_conv->SetParameter(3, 1); // scaling
  f1_conv->SetParNames("Offset", "b [/kHz]", "Center", "Scaling");
  f1_conv->SetLineColor(4);
  
  //con_curve->SetMarkerColor(kRed);
  con_curve->Fit("f1_conv","EM", "", -400, 400);
 
  mg->Add(con_curve);


  
  std::string start_str = std::to_string(6.5) + "e-3";
  std::string end_str = std::to_string(7.5) + "e-3";
  std::string fit_no_old = "[0]+([3]*2*[1]*[1]*(TMath::Exp(-" + std::to_string(gamma) + '*';
  std::string fit_no_old2 = ")*(1-(TMath::Cos((TMath::Sqrt(4*TMath::Pi()*TMath::Pi()*(x-[2])*(x-[2])+4*[1]*[1]))*";
  std::string fit_no_old3 = ")-(TMath::Sqrt(4*TMath::Pi()*TMath::Pi()*(x-[2])*(x-[2])+4*[1]*[1]))*TMath::Sin((TMath::Sqrt(4*TMath::Pi()*TMath::Pi()*(x-[2])*(x-[2])+4*[1]*[1]))*";
  std::string fit_no_old4 = ")/"+std::to_string(gamma)+")*"+std::to_string(gamma)+'*'+std::to_string(gamma)+"/(4*TMath::Pi()*TMath::Pi()*(x-[2])*(x-[2])+4*[1]*[1]+"+std::to_string(gamma)+'*'+std::to_string(gamma)+"))-TMath::Exp(-"+std::to_string(gamma)+'*';
  std::string fit_no_old5 = ")*(1-(TMath::Cos(TMath::Sqrt(4*TMath::Pi()*TMath::Pi()*(x-[2])*(x-[2])+4*[1]*[1])*";
  std::string fit_no_old6 = ")-TMath::Sqrt(4*TMath::Pi()*TMath::Pi()*(x-[2])*(x-[2])+4*[1]*[1])*TMath::Sin(TMath::Sqrt(4*TMath::Pi()*TMath::Pi()*(x-[2])*(x-[2])+4*[1]*[1])*";
  std::string fit_no_old7 = ")/"+std::to_string(gamma)+")*"+std::to_string(gamma)+'*'+std::to_string(gamma)+"/(4*TMath::Pi()*TMath::Pi()*(x-[2])*(x-[2])+4*[1]*[1]+"+std::to_string(gamma)+'*'+std::to_string(gamma)+")))/(4*TMath::Pi()*TMath::Pi()*(x-[2])*(x-[2])+4*[1]*[1])";
  std::string fit_no_old8 = ")/(TMath::Exp(-"+std::to_string(gamma)+'*'+start_str+")-TMath::Exp(-"+std::to_string(gamma)+'*'+end_str+"))";
  TF1* f1_old = new TF1("f1_old",(fit_no_old+start_str+fit_no_old2+start_str+fit_no_old3+start_str+fit_no_old4
				  +end_str+fit_no_old5+end_str+fit_no_old6+end_str+fit_no_old7+fit_no_old8).c_str(), -400, 400);
  
  f1_old->FixParameter(0, 0); // offset                                                                                                                                                                  
  f1_old->SetParameter(1, 124.); // microwave power                                                                                                                                                      
  f1_old->SetParameter(2, 0); // center                                                                                                                                                                  
  f1_old->SetParameter(3, 1); // scaling                                                                                                                                                                 
  f1_old->SetParNames("Offset", "b [/kHz]", "Center", "Scaling");
  f1_old->SetLineColor(2);

  //old_curve->SetMarkerColor(kBlue);
  old_curve->Fit("f1_old","EM", "", -400, 400);
  
  mg->Add(old_curve);  


  mg->Draw("AP");
  mg->GetXaxis()->SetTitle("Frequency Detuning [/kHz]");
  mg->GetYaxis()->SetTitle("Signal");
  

  c3->Update();
  
  TPaveStats* stats1 = (TPaveStats*) con_curve->GetListOfFunctions()->FindObject("stats");
  TPaveStats* stats2 = (TPaveStats*) old_curve->GetListOfFunctions()->FindObject("stats");
  stats1->SetTextColor(kBlue);
  stats2->SetTextColor(kRed);
  stats1->SetX1NDC(0.12); stats1->SetX2NDC(0.32); stats1->SetY1NDC(0.75);
  stats2->SetX1NDC(0.72); stats2->SetX2NDC(0.92); stats2->SetY1NDC(0.78);
  c3->Modified();
  
  c3->SaveAs("con_VS_old.png");
}
