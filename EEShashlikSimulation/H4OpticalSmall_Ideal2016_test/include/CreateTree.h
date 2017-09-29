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

  float  Eabs_3x3;
  float  Eact_3x3;

  float  Eabs_1x3;
  float  Eact_1x3;

  int    nLayers;

  float Eact_CentralXtal;
  float Eabs_CentralXtal;

  float  xPosition;
  float  yPosition;
  

    
};
