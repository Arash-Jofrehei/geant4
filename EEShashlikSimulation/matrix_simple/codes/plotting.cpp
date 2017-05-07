#include <iostream>
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TF1.h"
#include "string.h"
using namespace std;
float xPosition;
float yPosition;
int  Fibre_0;
int  Fibre_1;
int  Fibre_2;
int  Fibre_3;
int  NPhot_Fib;
/*int  NPhot_Fib2;
int  NPhot_Fib3;
int  NPhot_Fib4;*/
vector<float> *EAPD;
int nParticlesAPD;
double APD1;
int nbins_spike = 76;
int nbins_WLS = 76;

void plotting(){
  TCanvas *canvas = new TCanvas("plotting","plotting");
  TFile *file = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/geant4/EEShashlikSimulation/matrix_simple/cmake/ntuples/fiber19000.root");
  TTree *tree = (TTree*) file->Get("tree");
  tree->SetBranchAddress("xPosition", &xPosition);
  tree->SetBranchAddress("yPosition", &yPosition);
  tree->SetBranchAddress("Fibre_0", &Fibre_0);
  tree->SetBranchAddress("Fibre_1", &Fibre_1);
  tree->SetBranchAddress("Fibre_2", &Fibre_2);
  tree->SetBranchAddress("Fibre_3", &Fibre_3);
  tree->SetBranchAddress("NPhot_Fib", &NPhot_Fib);
  /*tree->SetBranchAddress("NPhot_Fib2", &NPhot_Fib2);
  tree->SetBranchAddress("NPhot_Fib3", &NPhot_Fib3);
  tree->SetBranchAddress("NPhot_Fib4", &NPhot_Fib4);*/
  tree->SetBranchAddress("EAPD", &EAPD);
  tree->SetBranchAddress("nParticlesAPD", &nParticlesAPD);
  TProfile2D *energy_spike = new TProfile2D("spike energy","# spike MIPs with geant4 simulation;X(mm);Y(mm)",nbins_spike,-9.5,9.5,nbins_spike,-9.5,9.5);
  TProfile2D *number_spike = new TProfile2D("spike count","# spiking particles with geant4 simulation;X(mm);Y(mm)",nbins_spike,-9.5,9.5,nbins_spike,-9.5,9.5);
  TProfile2D *count_WLS = new TProfile2D("WLS contribution","# WLS optical photons with geant4 simulation;X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  TProfile2D *count_fscint = new TProfile2D("fiber scintillation contribution","# fiber scintillation optical photons with geant4 simulation;X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  Long64_t nentries = tree->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++){
    Long64_t ientry = tree->LoadTree(jentry);
    nb = tree->GetEntry(jentry);   nbytes += nb;
    //APD1 = EAPD[0];
    energy_spike->Fill(xPosition+18.5,yPosition,EAPD->at(0)/0.0019);
    number_spike->Fill(xPosition+18.5,yPosition,nParticlesAPD);
    count_WLS->Fill(xPosition+18.5,yPosition,Fibre_0);
    count_fscint->Fill(xPosition+18.5,yPosition,NPhot_Fib);
  }
  //cout << "spike/WLS ratio center:    " << energy_spike->GetBinContent(20,20) / count_WLS->GetBinContent(10,11) << endl;
  //cout << "spike/WLS ratio near the fiber:    " << energy_spike->GetBinContent(5,5) / count_WLS->GetBinContent(3,3) << endl;
  energy_spike->Draw("colz");
  TCanvas *canvas2 = new TCanvas("plotting2","plotting2");
  number_spike->Draw("colz");
  TCanvas *canvas3 = new TCanvas("plotting3","plotting3");
  count_WLS->Draw("colz");
  TCanvas *canvas4 = new TCanvas("plotting4","plotting4");
  count_fscint->Draw("colz");
}