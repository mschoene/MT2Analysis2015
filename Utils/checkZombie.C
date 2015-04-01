#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"

#include <iostream> 

using namespace std;

//Note (Boris 20/03/15): This macro should be extended with the extra checks that Bruno was talking about


int checkZombie(string input){
  cout << "============= testing: " << input << endl;
  TFile* f = TFile::Open( ("dcap://t3se01.psi.ch:22125/"+input).c_str() );
  if(f->IsZombie()){
    cout << "WARNING: file is a Zombie: " << input << endl;
    f->Close(); delete f; return 1;
  }

  TH1D* count = (TH1D*)f->Get("Count");
  if(count->IsZombie())  {
    cout << "WARNING: histo is zombie :-|" << endl;
    f->Close(); delete f; return 1;
  }
  //else
  //cout << "histo has #entries: " << count->GetBinContent(1) << endl;

  TTree *tree = (TTree*)f->Get("tree");
  if(tree->IsZombie()) {
    cout << "WARNING: tree is zombie :-|" << endl;
    f->Close(); delete f; return 1;
  }
  //else
  //cout << "tree has #entries: " << tree->GetEntries() << endl;


  f->Close();
  delete f;
  return 0;
}
