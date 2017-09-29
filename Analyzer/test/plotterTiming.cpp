#include "Analyzer.h"

int main( int argc, char* argv[] ) {

  //DrawTools::setStyle();

  std::string setup = "WLS";
  std::string energy = "50GeV";

  if( argc>1 ) {
    std::string setup_str(argv[1]);
    setup = setup_str;
    if( argc>2 ) {
      std::string energy_str(argv[2]);
      energy=energy_str;
    }
  } else {
    std::cout << "Usage:" << std::endl;
    std::cout << "./plotterTiming ([setup]) ([energy])" << std::endl;
    exit(12345);
  }

  std::cout<<"plotting for setup:"<<setup<<" energy:"<<energy<<std::endl;

  TString fileName = "files/totalOut_"+setup+"_"+energy+".root";

  TFile* inFile = TFile::Open(fileName);
  TTree* tree = (TTree*) inFile->Get("tree");
  Analyzer t(tree);
  t.Loop(setup, energy);

  return 0;

}
