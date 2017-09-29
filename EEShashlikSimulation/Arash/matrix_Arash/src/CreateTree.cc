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
  
  this->GetTree()->Branch("Time_deposit",&this->Time_deposit);
  this->GetTree()->Branch("Process_deposit",&this->Process_deposit);
  this->GetTree()->Branch("Z_deposit",&this->Z_deposit);
  this->GetTree()->Branch("Theta_deposit",&this->Theta_deposit);
  this->GetTree()->Branch("opPhoton_time",&this->opPhoton_time);    
  this->GetTree()->Branch("opPhoton_process",&this->opPhoton_process);    
  this->GetTree()->Branch("Eabs",&this->Eabs,"Eabs/F");
  this->GetTree()->Branch("Eact",&this->Eact,"Eact/F");
  this->GetTree()->Branch("NPhot_Act",&this->NPhot_Act,"NPhot_Act/I");
  this->GetTree()->Branch("NPhot_Fib",&this->NPhot_Fib,"NPhot_Fib/I");
  this->GetTree()->Branch("NPhot_Fib2",&this->NPhot_Fib2,"NPhot_Fib2/I");
  this->GetTree()->Branch("NPhot_Fib3",&this->NPhot_Fib3,"NPhot_Fib3/I");
  this->GetTree()->Branch("NPhot_Fib4",&this->NPhot_Fib4,"NPhot_Fib4/I");
  this->GetTree()->Branch("EfibrCore",&this->EfibrCore,"EfibrCore/F");
  this->GetTree()->Branch("EfibrClad",&this->EfibrClad,"EfibrClad/F");
  this->GetTree()->Branch("nLayers",&this->nLayers,"nLayers/I");

  this->GetTree()->Branch("Fibre_start_0",&this->Fibre_start_0,"Fibre_start_0/I");

  this->GetTree()->Branch("Fibre_0",&this->Fibre_0,"Fibre_0/I");
  this->GetTree()->Branch("Fibre_1",&this->Fibre_1,"Fibre_1/I");
  this->GetTree()->Branch("Fibre_2",&this->Fibre_2,"Fibre_2/I");
  this->GetTree()->Branch("Fibre_3",&this->Fibre_3,"Fibre_3/I");

  this->GetTree()->Branch("xPosition",&this->xPosition,"xPosition/F");
  this->GetTree()->Branch("yPosition",&this->yPosition,"yPosition/F");

  this->GetTree()->Branch("EOpt_0",&this->EOpt_0,"EOpt_0/F");
  this->GetTree()->Branch("EOpt_1",&this->EOpt_1,"EOpt_1/F");
  this->GetTree()->Branch("EOpt_2",&this->EOpt_2,"EOpt_2/F");
  this->GetTree()->Branch("EOpt_3",&this->EOpt_3,"EOpt_3/F");
  
  this->GetTree()->Branch("EAPD",&this->EAPD);
  this->GetTree()->Branch("Time_deposit_APD",&this->Time_deposit_APD);

  this->GetTree()->Branch("nParticlesAPD",&this->nParticlesAPD,"nParticlesAPD/I");

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
  Time_deposit.clear();
  Process_deposit.clear();
  Z_deposit.clear();
  Theta_deposit.clear();
  opPhoton_time.clear();
  opPhoton_process.clear();
  Eabs=0;
  Eact=0;
  Eabs_3x3=0;
  Eact_3x3=0;

  Eabs_1x3=0;
  Eact_1x3=0;
  NPhot_Act=0;
  NPhot_Fib=0;
  NPhot_Fib2=0;
  NPhot_Fib3=0;
  NPhot_Fib4=0;
  EfibrCore=0;
  EfibrClad=0;
  nLayers=0;
  
  Fibre_start_0=0;

  Fibre_0=0;
  Fibre_1=0;
  Fibre_2=0;
  Fibre_3=0;
  
  xPosition=0;
  yPosition=0;
  
  EOpt_0=0;
  EOpt_1=0;
  EOpt_2=0;
  EOpt_3=0;  
  
  Eact_CentralXtal=0.;
  Eabs_CentralXtal=0.;

  EAPD.clear();
  Time_deposit_APD.clear();

  nParticlesAPD=0;
}
