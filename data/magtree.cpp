#include <stdio.h>
#include <vector>
#include <string>
#include "TTree.h"
#include "TFile.h"
#define N 3602

void magtree(const char* filename = "./BRECON_MOM_20200716_6.txt"){
  std::ifstream ifs;
  ifs.open(filename, std::ios::in); // read only
  if(!ifs.is_open()){
    std::cout << "Failed to open the file..." << std::endl;
    return;
  }

  TFile* file = TFile::Open("magnet.root","RECREATE");
  TTree* tree = new TTree("tree","tree of magnet moments");
  
  int i=0;
  int dummy;
  Double_t x, y, z;
  Double_t mx, my, mz;
  tree->Branch("x",&x,"x/D");
  tree->Branch("y",&y,"y/D");
  tree->Branch("z",&z,"z/D");
  tree->Branch("mx",&mx,"mx/D");
  tree->Branch("my",&my,"my/D");
  tree->Branch("mz",&mz,"mz/D");
  
  while(!ifs.eof()){
    if(i==3602) break;
    ifs >> dummy >> x >> y >> z >> mx >> my >> mz;
    i++;
    tree->Fill();
  }

  tree->Scan("*");
  //std::cout << "entries:" << tree->GetEntries() << std::endl;
  if(tree->Write()) std::cout  << "magtree.root is made." << std::endl;
  file->Close();
  ifs.close();
  return;
}
