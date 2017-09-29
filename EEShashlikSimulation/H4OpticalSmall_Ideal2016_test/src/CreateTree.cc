#include "CreateTree.h"



CreateTree* CreateTree::fInstance = NULL;



CreateTree::CreateTree(TString name)
{
  if( fInstance )
    {
      return;
    }

  this -> fInstance = this;
  this -> fname     = name;
  this -> ftree     = new TTree(name,name);

  this->GetTree()->Branch("Event",&this->Event,"Event/I");

  this->GetTree()->Branch("Eabs_3x3",&this->Eabs_3x3,"Eabs_3x3/F");
  this->GetTree()->Branch("Eact_3x3",&this->Eact_3x3,"Eact_3x3/F");

  this->GetTree()->Branch("Eabs_1x3",&this->Eabs_1x3,"Eabs_1x3/F");
  this->GetTree()->Branch("Eact_1x3",&this->Eact_1x3,"Eact_1x3/F");


  this->GetTree()->Branch("Eabs_CentralXtal",&this->Eabs_CentralXtal,"Eabs_CentralXtal");
  this->GetTree()->Branch("Eact_CentralXtal",&this->Eact_CentralXtal,"Eact_CentralXtal");

  this->GetTree()->Branch("nLayers",&this->nLayers,"nLayers/I");

  this->GetTree()->Branch("xPosition",&this->xPosition,"xPosition/F");
  this->GetTree()->Branch("yPosition",&this->yPosition,"yPosition/F");




}

bool CreateTree::Write()
{
  TString filename = this->GetName();
  filename+=".root";
  TFile* file = new TFile(filename,"RECREATE");
  this->GetTree()->Write();
  file->Write();
  file->Close();
  return true;
}


void CreateTree::Clear()
{
  Event = 0;

  Eabs_3x3=0;
  Eact_3x3=0;

  Eabs_1x3=0;
  Eact_1x3=0;

  nLayers=0;
  
  
  xPosition=0;
  yPosition=0;
  
  Eact_CentralXtal=0.;
  Eabs_CentralXtal=0.;

}
