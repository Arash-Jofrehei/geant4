#include <iostream>
#include "TCanvas.h"
using namespace std;
float xPosition;
float yPosition;
int  Fibre_0;
int  Fibre_1;
int  Fibre_2;
int  Fibre_3;
int  NPhot_Fib;
int  NPhot_Fib2;
int  NPhot_Fib3;
int  NPhot_Fib4;
vector<float> *Process_deposit;
vector<float> *Time_deposit;

void time_arrival(){
  TCanvas *canvas = new TCanvas("time_arrival","time_arrival");
  TFile *file = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/geant4/EEShashlikSimulation/ntuples/center_fiber30000.root");
  TTree *tree = (TTree*) file->Get("tree");
  tree->SetBranchAddress("xPosition", &xPosition);
  tree->SetBranchAddress("yPosition", &yPosition);
  tree->SetBranchAddress("Fibre_0", &Fibre_0);
  tree->SetBranchAddress("Fibre_1", &Fibre_1);
  tree->SetBranchAddress("Fibre_2", &Fibre_2);
  tree->SetBranchAddress("Fibre_3", &Fibre_3);
  tree->SetBranchAddress("NPhot_Fib", &NPhot_Fib);
  tree->SetBranchAddress("NPhot_Fib2", &NPhot_Fib2);
  tree->SetBranchAddress("NPhot_Fib3", &NPhot_Fib3);
  tree->SetBranchAddress("NPhot_Fib4", &NPhot_Fib4);
  tree->SetBranchAddress("Process_deposit", &Process_deposit);
  tree->SetBranchAddress("Time_deposit", &Time_deposit);
  TH1F *WLS_center = new TH1F("WLS time arrival at center","WLS time arrival at center",1024,-0.1,204.7);
  TH1F *WLS_fiber = new TH1F("WLS time arrival near fiber","WLS time arrival near fiber",1024,-0.1,204.7);
  TH1F *fscint_center = new TH1F("fiber scintillation time arrival at center","fiber scintillation time arrival at center",1024,-0.1,204.7);
  TH1F *fscint_fiber = new TH1F("fiber scintillation time arrival near fiber","fiber scintillation time arrival near fiber",1024,-0.1,204.7);
  
  Long64_t nentries = tree->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++){
    Long64_t ientry = tree->LoadTree(jentry);
    nb = tree->GetEntry(jentry);   nbytes += nb;
    for (int i = 0;i < Process_deposit->size();i++){
      if (TMath::Abs(xPosition+18.5)<1 && TMath::Abs(yPosition)<1 && Process_deposit->at(i)==1) WLS_center->Fill(Time_deposit->at(i)-5.0);
      if (TMath::Abs(xPosition+25.5)<1 && TMath::Abs(yPosition+7)<1 && Process_deposit->at(i)==1) WLS_fiber->Fill(Time_deposit->at(i)-5.0);
      if (TMath::Abs(xPosition+18.5)<1 && TMath::Abs(yPosition)<1 && Process_deposit->at(i)==2) fscint_center->Fill(Time_deposit->at(i)-5.0);
      if (TMath::Abs(xPosition+25.5)<1 && TMath::Abs(yPosition+7)<1 && Process_deposit->at(i)==2) fscint_fiber->Fill(Time_deposit->at(i)-5.0);
    }
  }
  WLS_center->Scale(1.0/WLS_center->Integral());
  WLS_fiber->Scale(1.0/WLS_fiber->Integral());
  fscint_center->Scale(1.0/fscint_center->Integral());
  fscint_fiber->Scale(1.0/fscint_fiber->Integral());
  TFile *time_arrivals = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/time_arrivals.root","recreate");
  time_arrivals->cd();
  WLS_center->Write();
  WLS_fiber->Write();
  fscint_center->Write();
  fscint_fiber->Write();
  
  WLS_center->SetLineColor(4);
  WLS_fiber->SetLineColor(3);
  fscint_fiber->SetLineColor(2);
  fscint_fiber->Draw();
  WLS_fiber->Draw("same");
  fscint_center->Draw("same");
  WLS_center->Draw("same");
}