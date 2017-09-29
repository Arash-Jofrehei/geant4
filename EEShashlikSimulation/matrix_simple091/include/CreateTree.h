#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"

class CreateTree
{
 private:
  
  TTree*  ftree;
  TString fname;
  
 public:
  
  
  CreateTree(TString name);
  TTree*             GetTree() const { return ftree; };
  TString            GetName() const { return fname; };
  int                Fill() { return this->GetTree()->Fill(); };
  bool               Write();
  void               Clear();
  static CreateTree* Instance() { return fInstance; };
  static CreateTree* fInstance;
  int Event;
  std::vector<float> Time_deposit;
  std::vector<float> Z_deposit;
  std::vector<float> Process_deposit;//1 is WLS, 2 Scintillation, 3 Cerenkov
  std::vector<float> Theta_deposit;
  std::vector<float> opPhoton_time;
  std::vector<int> opPhoton_process;
  float  Eabs;
  float  Eact;
  float  Eabs_3x3;
  float  Eact_3x3;

  float  Eabs_1x3;
  float  Eact_1x3;
  float Eact_CentralXtal;
  float Eabs_CentralXtal;
  int    NPhot_Act;
  int    NPhot_Fib;
  int    NPhot_Fib2;
  int    NPhot_Fib3;
  int    NPhot_Fib4;
  float  EfibrCore;
  float  EfibrClad;
  int    nLayers;

  int   Fibre_start_0;

  int   Fibre_0;
  int   Fibre_1;
  int   Fibre_2;
  int   Fibre_3;

  float  xPosition;
  float  yPosition;
  
  float  EOpt_0;
  float  EOpt_1;
  float  EOpt_2;
  float  EOpt_3;  
  std::vector<float> EAPD;
  std::vector<float> Time_deposit_APD;
  std::vector<float> E_deposit_APD;

  int nParticlesAPD;
};
