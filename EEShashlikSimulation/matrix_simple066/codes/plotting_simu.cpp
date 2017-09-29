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
float Eact_CentralXtal;
float Eabs_CentralXtal;
float xPosition;
float yPosition;
float equi_xPosition;
float equi_yPosition;
int  Fibre_0;
int  Fibre_1;
int  Fibre_2;
int  Fibre_3;
int  NPhot_Fib;
int  NPhot_Fib2;
int  NPhot_Fib3;
int  NPhot_Fib4;
vector<float> *Time_deposit;
vector<float> *Process_deposit;
vector<float> *Time_deposit_APD;
vector<float> *EAPD;
int nParticlesAPD;
double APD1;
int nbins_spike = 76;
int nbins_WLS = 76;

double FWHM(TH1F *h){
  int bin1 = h->FindFirstBinAbove(h->GetMaximum()/2);
  int bin2 = h->FindLastBinAbove(h->GetMaximum()/2);
  //int bin1 = h->GetMaximumBin();
  float fwhm = h->GetBinCenter(bin2) - h->GetBinCenter(bin1);
  return fwhm;
}

/*TH1F *conv;
TH1F *convolution(TH1F *f1,TH1F *f2,int size = 1024){
  double dummy;
  for (int sample = 0;sample < size;sample++){
    dummy = 0;
    for (int point = 0;point <= sample;point++){
      dummy += f1->GetBinContent(point+1)*f2->GetBinContent(1+sample-point);
    }
    conv->SetBinContent(sample+1,dummy);
  }
  return conv;
}*/

TH1F *conv;
TH1F *convolution(TH1F *f1,TH1F *f2,int size = 1024){
  conv->Reset();
  for (int sample = 0;sample < size;sample++){
    if (f1->GetBinContent(sample+1) == 0) continue;
    for (int point = sample;point < size;point++) conv->Fill(point/5.0,f1->GetBinContent(sample+1) * f2->GetBinContent(point+1-sample));
  }
  return conv;
}


void plotting_simu(){
  TCanvas *canvas = new TCanvas("plotting","plotting");
  TFile *convoluted_pulses = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/convoluted_pulses.root");
  apd_plus_electronics = (TH1F*) convoluted_pulses->Get("APD+electronics");
  conv = new TH1F("convolution","convolution",1024,-0.1,204.7);
  //TFile *file = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/geant4/EEShashlikSimulation/ntuples/full_yield.root");
  //TFile *file = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/geant4/EEShashlikSimulation/ntuples/full_face_119000.root");
  TFile *file = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/geant4/EEShashlikSimulation/matrix_simple/cmake/ntuples/runEEShashlik.root");
  TTree *tree = (TTree*) file->Get("tree");
  tree->SetBranchAddress("Eact_CentralXtal", &Eact_CentralXtal);
  tree->SetBranchAddress("Eabs_CentralXtal", &Eabs_CentralXtal);
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
  tree->SetBranchAddress("Time_deposit", &Time_deposit);
  tree->SetBranchAddress("Process_deposit", &Process_deposit);
  tree->SetBranchAddress("Time_deposit_APD", &Time_deposit_APD);
  tree->SetBranchAddress("EAPD", &EAPD);
  tree->SetBranchAddress("nParticlesAPD", &nParticlesAPD);
  TProfile2D *energy_spike = new TProfile2D("spike energy","# spike MIPs with geant4 simulation;X(mm);Y(mm)",nbins_spike,-9.5,9.5,nbins_spike,-9.5,9.5);
  TProfile2D *number_spike = new TProfile2D("spike count","# spiking particles with geant4 simulation;X(mm);Y(mm)",nbins_spike,-9.5,9.5,nbins_spike,-9.5,9.5);
  TProfile2D *count_WLS = new TProfile2D("# WLS optical photons with geant4 simulation (4APDs)","# WLS optical photons with geant4 simulation (4APDs);X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  TProfile2D *count_fscint = new TProfile2D("# fiber scintillation optical photons with geant4 simulation (4APDs)","# fiber scintillation optical photons with geant4 simulation (4APDs);X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  TProfile2D *act_deposit = new TProfile2D("Active Material Deposited Energy - Central Crystal","Active Material Deposited Energy - Central Crystal;X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  
  TProfile2D *WLS_waveform_amp[5];
  WLS_waveform_amp[0] = new TProfile2D("4APDs WLS waveform amplitude","4APDs WLS waveform amplitude;X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  WLS_waveform_amp[1] = new TProfile2D("APD1 WLS waveform amplitude","APD1 WLS waveform amplitude;X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  WLS_waveform_amp[2] = new TProfile2D("APD2 WLS waveform amplitude","APD2 WLS waveform amplitude;X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  WLS_waveform_amp[3] = new TProfile2D("APD3 WLS waveform amplitude","APD3 WLS waveform amplitude;X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  WLS_waveform_amp[4] = new TProfile2D("APD4 WLS waveform amplitude","APD4 WLS waveform amplitude;X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  
  TProfile2D *fscint_waveform_amp[5];
  fscint_waveform_amp[0] = new TProfile2D("4APDs fiber scintillation waveform amplitude","4APDs fiber scintillation waveform amplitude;X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  fscint_waveform_amp[1] = new TProfile2D("APD1 fiber scintillation waveform amplitude","APD1 fiber scintillation waveform amplitude;X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  fscint_waveform_amp[2] = new TProfile2D("APD2 fiber scintillation waveform amplitude","APD2 fiber scintillation waveform amplitude;X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  fscint_waveform_amp[3] = new TProfile2D("APD3 fiber scintillation waveform amplitude","APD3 fiber scintillation waveform amplitude;X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  fscint_waveform_amp[4] = new TProfile2D("APD4 fiber scintillation waveform amplitude","APD4 fiber scintillation waveform amplitude;X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  
  TProfile2D *overall_waveform_amp[5];
  overall_waveform_amp[0] = new TProfile2D("4APDs simulated waveform amplitude","4APDs simulated waveform amplitude;X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  overall_waveform_amp[1] = new TProfile2D("APD1 simulated waveform amplitude","APD1 simulated waveform amplitude;X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  overall_waveform_amp[2] = new TProfile2D("APD2 simulated waveform amplitude","APD2 simulated waveform amplitude;X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  overall_waveform_amp[3] = new TProfile2D("APD3 simulated waveform amplitude","APD3 simulated waveform amplitude;X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  overall_waveform_amp[4] = new TProfile2D("APD4 simulated waveform amplitude","APD4 simulated waveform amplitude;X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  
  TProfile2D *WLS_waveform_FWHM[5];
  WLS_waveform_FWHM[0] = new TProfile2D("4APDs WLS waveform FWHM","4APDs WLS waveform FWHM;X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  WLS_waveform_FWHM[1] = new TProfile2D("APD1 WLS waveform FWHM","APD1 WLS waveform FWHM;X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  WLS_waveform_FWHM[2] = new TProfile2D("APD2 WLS waveform FWHM","APD2 WLS waveform FWHM;X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  WLS_waveform_FWHM[3] = new TProfile2D("APD3 WLS waveform FWHM","APD3 WLS waveform FWHM;X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  WLS_waveform_FWHM[4] = new TProfile2D("APD4 WLS waveform FWHM","APD4 WLS waveform FWHM;X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  
  TProfile2D *fscint_waveform_FWHM[5];
  fscint_waveform_FWHM[0] = new TProfile2D("4APDs fiber scintillation waveform FWHM","4APDs fiber scintillation waveform FWHM;X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  fscint_waveform_FWHM[1] = new TProfile2D("APD1 fiber scintillation waveform FWHM","APD1 fiber scintillation waveform FWHM;X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  fscint_waveform_FWHM[2] = new TProfile2D("APD2 fiber scintillation waveform FWHM","APD2 fiber scintillation waveform FWHM;X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  fscint_waveform_FWHM[3] = new TProfile2D("APD3 fiber scintillation waveform FWHM","APD3 fiber scintillation waveform FWHM;X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  fscint_waveform_FWHM[4] = new TProfile2D("APD4 fiber scintillation waveform FWHM","APD4 fiber scintillation waveform FWHM;X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  
  TProfile2D *overall_waveform_FWHM[5];
  overall_waveform_FWHM[0] = new TProfile2D("4APDs simulated waveform FWHM","4APDs simulated waveform FWHM;X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  overall_waveform_FWHM[1] = new TProfile2D("APD1 simulated waveform FWHM","APD1 simulated waveform FWHM;X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  overall_waveform_FWHM[2] = new TProfile2D("APD2 simulated waveform FWHM","APD2 simulated waveform FWHM;X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  overall_waveform_FWHM[3] = new TProfile2D("APD3 simulated waveform FWHM","APD3 simulated waveform FWHM;X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  overall_waveform_FWHM[4] = new TProfile2D("APD4 simulated waveform FWHM","APD4 simulated waveform FWHM;X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  
  
  TH1F *time_arrival_WLS = new TH1F("WLS arrival time","WLS arrival time",1024,-0.1,204.7);
  TH1F *time_arrival_fscint = new TH1F("fiber scintillation arrival time","fiber scintillation arrival time",1024,-0.1,204.7);
  TH1F *time_arrival_overall = new TH1F("overall arrival time","overall arrival time",1024,-0.1,204.7);
  
  TH1F *WLS_waveform = new TH1F("WLS simulated waveform","WLS simulated waveform",1024,-0.1,204.7);
  TH1F *fscint_waveform = new TH1F("fiber scintillation simulated waveform","fiber scintillation simulated waveform",1024,-0.1,204.7);
  TH1F *overall_waveform = new TH1F("overall simulated waveform","overall simulated waveform",1024,-0.1,204.7);
  
  TH1F *WLS_waveform_4APDs = new TH1F("WLS simulated waveform 4APDs","WLS simulated waveform 4APDs",1024,-0.1,204.7);
  TH1F *fscint_waveform_4APDs = new TH1F("fiber scintillation simulated waveform 4APDs","fiber scintillation simulated waveform 4APDs",1024,-0.1,204.7);
  TH1F *overall_waveform_4APDs = new TH1F("overall simulated waveform 4APDs","overall simulated waveform 4APDs",1024,-0.1,204.7);
  
  float x_rotation[8] = {1,-1,1,-1,1,-1,1,-1};
  float y_rotation[8] = {1,-1,-1,-1,1,-1,-1,-1};
  
  //--------------------------------------------------------------------------------------------------------------------
  
  Long64_t nentries = tree->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++){
    Long64_t ientry = tree->LoadTree(jentry);
    nb = tree->GetEntry(jentry);   nbytes += nb;
    //if (jentry%1000 != 0) continue;
    if (jentry%1000 == 0) cout << jentry << " / " << nentries << endl;
    if(nParticlesAPD>0) energy_spike->Fill(xPosition+18.5,yPosition,EAPD->at(0)/0.0011);
    if(nParticlesAPD>0) number_spike->Fill(xPosition+18.5,yPosition,nParticlesAPD);
    //count_WLS->Fill(xPosition+18.5,yPosition,Fibre_0);
    //count_fscint->Fill(xPosition+18.5,yPosition,NPhot_Fib);
    count_WLS->Fill(xPosition+18.5,yPosition,Fibre_0+Fibre_1+Fibre_2+Fibre_3);
    count_fscint->Fill(xPosition+18.5,yPosition,NPhot_Fib+NPhot_Fib2+NPhot_Fib3+NPhot_Fib4);
    act_deposit->Fill(xPosition+18.5,yPosition,Eact_CentralXtal);
    if (NPhot_Fib2==0&&NPhot_Fib3==0&&NPhot_Fib4==0&&Fibre_1==0&&Fibre_2==0&&Fibre_3==0) continue;
    for (int APD = 1;APD < 5;APD++){
      time_arrival_WLS->Reset();
      time_arrival_fscint->Reset();
      time_arrival_overall->Reset();
      WLS_waveform->Reset();
      fscint_waveform->Reset();
      overall_waveform->Reset();
      if (APD == 1){
        WLS_waveform_4APDs->Reset();
        fscint_waveform_4APDs->Reset();
        overall_waveform->Reset();
      }
      for (int i = 0;i < Time_deposit->size();i++){
        if (Process_deposit->at(i) == (3*(APD-1)+1)) time_arrival_WLS->Fill(Time_deposit->at(i));
        if (Process_deposit->at(i) == (3*(APD-1)+2)) time_arrival_fscint->Fill(Time_deposit->at(i));
      }
      //WLS_waveform->Add(convolution(time_arrival_WLS,apd_plus_electronics));
      //fscint_waveform->Add(convolution(time_arrival_fscint,apd_plus_electronics));
      overall_waveform->Add(WLS_waveform);
      overall_waveform->Add(fscint_waveform);
      WLS_waveform_4APDs->Add(WLS_waveform);
      fscint_waveform_4APDs->Add(fscint_waveform);
      overall_waveform_4APDs->Add(overall_waveform);
      equi_xPosition = xPosition + 18.5;
      equi_yPosition = yPosition;
      for (int i = 0;i < 8;i++){
        if (i == 4){
          equi_xPosition = yPosition;
          equi_yPosition = xPosition + 18.5;
        }
        equi_xPosition *= x_rotation[i];
        equi_yPosition *= y_rotation[i];
        if (APD == 1) act_deposit->Fill(equi_xPosition,equi_yPosition,Eact_CentralXtal);
        WLS_waveform_amp[APD]->Fill(equi_xPosition,equi_yPosition,WLS_waveform->GetBinContent(WLS_waveform->GetMaximumBin()));
        fscint_waveform_amp[APD]->Fill(equi_xPosition,equi_yPosition,fscint_waveform->GetBinContent(fscint_waveform->GetMaximumBin()));
        overall_waveform_amp[APD]->Fill(equi_xPosition,equi_yPosition,overall_waveform->GetBinContent(overall_waveform->GetMaximumBin()));
        WLS_waveform_FWHM[APD]->Fill(equi_xPosition,equi_yPosition,FWHM(WLS_waveform));
        fscint_waveform_FWHM[APD]->Fill(equi_xPosition,equi_yPosition,FWHM(fscint_waveform));
        overall_waveform_FWHM[APD]->Fill(equi_xPosition,equi_yPosition,FWHM(overall_waveform));
        if (APD == 4){
          WLS_waveform_amp[0]->Fill(equi_xPosition,equi_yPosition,WLS_waveform_4APDs->GetBinContent(WLS_waveform_4APDs->GetMaximumBin()));
          fscint_waveform_amp[0]->Fill(equi_xPosition,equi_yPosition,fscint_waveform_4APDs->GetBinContent(fscint_waveform_4APDs->GetMaximumBin()));
          overall_waveform_amp[0]->Fill(equi_xPosition,equi_yPosition,overall_waveform_4APDs->GetBinContent(overall_waveform_4APDs->GetMaximumBin()));
          WLS_waveform_FWHM[0]->Fill(equi_xPosition,equi_yPosition,FWHM(WLS_waveform_4APDs));
          fscint_waveform_FWHM[0]->Fill(equi_xPosition,equi_yPosition,FWHM(fscint_waveform_4APDs));
          overall_waveform_FWHM[0]->Fill(equi_xPosition,equi_yPosition,FWHM(overall_waveform_4APDs));
        }
      }
    }
  }
  /*TFile *simulated_waveforms_position_profiles = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/full_yield_simulated_waveforms_position_profiles.root","recreate");
  //TFile *simulated_waveforms_position_profiles = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/simulated_waveforms_position_profiles.root","recreate");
  simulated_waveforms_position_profiles->cd();
  for (int i = 0;i < 5;i++){
    WLS_waveform_amp[i]->Write();
    fscint_waveform_amp[i]->Write();
    overall_waveform_amp[i]->Write();
    WLS_waveform_FWHM[i]->Write();
    fscint_waveform_FWHM[i]->Write();
    overall_waveform_FWHM[i]->Write();
  }*/
  /*
  TFile *simulated_WLS_fscint_deposit = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/simulated_WLS_fscint_deposit.root","recreate");
  simulated_WLS_fscint_deposit->cd();
  act_deposit->Write();
  count_WLS->Write();
  count_fscint->Write();
  */
  //TFile *simulated_waveforms_position_profiles = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/simulated_waveforms_position_profiles.root");
  /*TFile *simulated_waveforms_position_profiles = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/full_yield_simulated_waveforms_position_profiles.root");
  
  //-----------------------------------------------------------------------------------------
  energy_spike->Draw("colz");
  simulated_waveforms_position_profiles->Get("4APDs WLS waveform amplitude")->Draw("colz");
  TCanvas *canvas2 = new TCanvas("plotting2","plotting2");
  number_spike->Draw("colz");
  simulated_waveforms_position_profiles->Get("4APDs fiber scintillation waveform amplitude")->Draw("colz");
  TCanvas *canvas3 = new TCanvas("plotting3","plotting3");
  count_WLS->Draw("colz");
  simulated_waveforms_position_profiles->Get("4APDs simulated waveform amplitude")->Draw("colz");
  TCanvas *canvas4 = new TCanvas("plotting4","plotting4");
  count_fscint->Draw("colz");
  simulated_waveforms_position_profiles->Get("4APDs WLS waveform FWHM")->Draw("colz");
  TCanvas *canvas5 = new TCanvas("plotting5","plotting5");
  act_deposit->Draw("colz");
  simulated_waveforms_position_profiles->Get("4APDs fiber scintillation waveform FWHM")->Draw("colz");
  TCanvas *canvas6 = new TCanvas("plotting6","plotting6");
  simulated_waveforms_position_profiles->Get("4APDs simulated waveform FWHM")->Draw("colz");*/
  count_WLS->Draw("colz");
}