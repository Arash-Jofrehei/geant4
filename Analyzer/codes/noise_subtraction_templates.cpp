#include <iostream>
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
using namespace std;
float WF_val1[1024];
float WF_val2[1024];
float WF_val3[1024];
float WF_val6[1024];
float WF_val7[1024];
float WF_val8[1024];
float WF_val9[1024];
float WF_time10[1024];
float WF_val10[1024];
float WF_val11[1024];
float WF_val14[1024];
float WF_val15[1024];
float WF_val16[1024];
float WF_val17[1024];
float Time[18];
Float_t         x[2];
Float_t         y[2];
Int_t           nFibresOnX[2];
Int_t           nFibresOnY[2];
TBranch        *b_nFibresOnX;
TBranch        *b_nFibresOnY;
float channel_weights[18] = {0,1,1,1,0,0,0.4,0.4,0.4,0.4,1,1,0,0,1,1,1,0};//{0.521416,0.791515,0.524998,0.504333,0.506726,0.517119,0.50302,1,1.00179,0.512666,0.622622,0.899577};//{0.932688,1.56220,0.850238,0.492765,0.482391,0.606879,0.508326,1,1.37240,0.963999,1.44540,0.862974};
TH1F* apd_plus_electronics;
TH1F* convoluted_pulse_fscint;
TH1F* convoluted_pulse_WLS;
double fit_function(double *v,double *par)
{
  int n = int(v[0]*5);
  int shift = int(par[2]);
  return par[0]*(par[1]*convoluted_pulse_fscint->GetBinContent(n+1+shift)+(1-par[1])*convoluted_pulse_WLS->GetBinContent(n+1+shift))+par[3];
  //return par[0]*(par[1]*apd_plus_electronics->GetBinContent(n+1+shift)+(1-par[1])*convoluted_pulse_WLS->GetBinContent(n+1+shift))+par[3];
}
void noise_subtraction_templates(){
  TCanvas *canvas = new TCanvas("noise_subtraction_templates","noise_subtraction_templates");
  //TFile *final = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/merged_crystal11.root");
  TFile *final = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/analysis_4683.root");
  TTree *ftree = (TTree*) final->Get("h4");
  TFile *convoluted_pulses = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/convoluted_pulses.root");
  convoluted_pulse_fscint = (TH1F*) convoluted_pulses->Get("normalized fiber scintillation convoluted with APD+electronics");
  convoluted_pulse_WLS = (TH1F*) convoluted_pulses->Get("normalized WLS+CeF3 convoluted with APD+electronics");
  apd_plus_electronics = (TH1F*) convoluted_pulses->Get("normalized APD+electronics");
  TH1F *waveform[13];
  waveform[0] = new TH1F("waveform1","waveform crystal 2",1024,-0.1,204.7);
  waveform[1] = new TH1F("waveform2","waveform crystal 3",1024,-0.1,204.7);
  waveform[2] = new TH1F("waveform3","waveform crystal 4",1024,-0.1,204.7);
  waveform[3] = new TH1F("waveform6","waveform crystal4APD 1",1024,-0.1,204.7);
  waveform[4] = new TH1F("waveform7","waveform crystal4APD 2",1024,-0.1,204.7);
  waveform[5] = new TH1F("waveform8","waveform crystal4APD 3",1024,-0.1,204.7);
  waveform[6] = new TH1F("waveform9","waveform crystal4APD 4",1024,-0.1,204.7);
  waveform[7] = new TH1F("waveform10","waveform crystal 11",1024,-0.1,204.7);
  waveform[8] = new TH1F("waveform11","waveform crystal 12",1024,-0.1,204.7);
  waveform[9] = new TH1F("waveform14","waveform crystal 15",1024,-0.1,204.7);
  waveform[10] = new TH1F("waveform15","waveform crystal 16",1024,-0.1,204.7);
  waveform[11] = new TH1F("waveform16","waveform crystal 17",256,-0.1,204.7);
  waveform[12] = new TH1F("waveform17","waveform crystal 18",256,-0.1,204.7);
  TH1 *difference_hist = new TH1F("difference","difference",256,-0.1,204.7);
  double integrated_charge = 0;
  double calibrated_integrated_charge = 0;
  TH1D *integrated_charge_crystal11 = new TH1D("integrated charge crystal 11","crystal 11 integrated charge (single crystal)",175,0,600000);
  TH1D *calibrated_integrated_charge_crystal11 = new TH1D("calibrated integrated charge crystal 11","crystal 11 integrated charge (3*3 matrix)",200,0,750000);
  TF1 *crys = new TF1("crys","crystalball",0,500000);
  TF1 *func = new TF1("fit",fit_function,0,200,3);
  func->SetParNames("integrated charge","fiber scintillation ratio","shift","added constant");
  func->SetParLimits(1,0,1);
  func->SetParLimits(2,40,100);
  ftree->SetBranchAddress("WF_val1",WF_val1);
  ftree->SetBranchAddress("WF_val2",WF_val2);
  ftree->SetBranchAddress("WF_val3",WF_val3);
  ftree->SetBranchAddress("WF_val6",WF_val6);
  ftree->SetBranchAddress("WF_val7",WF_val7);
  ftree->SetBranchAddress("WF_val8",WF_val8);
  ftree->SetBranchAddress("WF_val9",WF_val9);
  ftree->SetBranchAddress("WF_time10",WF_time10);
  ftree->SetBranchAddress("WF_val10",WF_val10);
  ftree->SetBranchAddress("WF_val11",WF_val11);
  ftree->SetBranchAddress("WF_val14",WF_val14);
  ftree->SetBranchAddress("WF_val15",WF_val15);
  ftree->SetBranchAddress("WF_val16",WF_val16);
  ftree->SetBranchAddress("WF_val17",WF_val17);
  ftree->SetBranchAddress("Time",Time);
  ftree->SetBranchAddress("x", x);
  ftree->SetBranchAddress("y", y);
  ftree->SetBranchAddress("nFibresOnX", nFibresOnX, &b_nFibresOnX);
  ftree->SetBranchAddress("nFibresOnY", nFibresOnY, &b_nFibresOnY);
  Long64_t nentries = ftree->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<30;jentry++){
    Long64_t ientry = ftree->LoadTree(jentry);
    nb = ftree->GetEntry(jentry);   nbytes += nb;
    if (TMath::Abs(x[1]-224.0)<2 && TMath::Abs(y[1]-294.5)<2){
      for (int i = 0;i < 12;i++) waveform[i]->Reset();
      difference_hist->Reset();
      for (int sample = 0;sample < 1024;sample++){
        waveform[0]->Fill(30+WF_time10[sample]-Time[1],channel_weights[1]*WF_val1[sample]);
        waveform[1]->Fill(30+WF_time10[sample]-Time[2],channel_weights[2]*WF_val2[sample]);
        waveform[2]->Fill(30+WF_time10[sample]-Time[10],channel_weights[3]*WF_val3[sample]);
        waveform[3]->Fill(30+WF_time10[sample]-Time[6],channel_weights[6]*WF_val6[sample]);
        waveform[4]->Fill(30+WF_time10[sample]-Time[7],channel_weights[7]*WF_val7[sample]);
        waveform[5]->Fill(30+WF_time10[sample]-Time[8],channel_weights[8]*WF_val8[sample]);
        waveform[6]->Fill(30+WF_time10[sample]-Time[9],channel_weights[9]*WF_val9[sample]);
        waveform[7]->Fill(30+WF_time10[sample]-Time[10],channel_weights[10]*WF_val10[sample]);
        waveform[8]->Fill(30+WF_time10[sample]-Time[11],channel_weights[11]*WF_val11[sample]);
        waveform[9]->Fill(30+WF_time10[sample]-Time[14],channel_weights[14]*WF_val14[sample]);
        waveform[10]->Fill(30+WF_time10[sample]-Time[15],channel_weights[15]*WF_val15[sample]);
        waveform[11]->Fill(30+WF_time10[sample]-Time[10],channel_weights[16]*WF_val16[sample]);
        waveform[12]->Fill(30+WF_time10[sample]-Time[10],WF_val17[sample]);
      }
      calibrated_integrated_charge = 0;
      func->SetParameters(300000,0.4,70,0,0);
      for (int j = 0;j < 10;j++){
        //waveform[12]->Fit("fit","","",30,80);
      }
      difference_hist->Add(waveform[11],1);
      //difference_hist->Add(func,-1);
      difference_hist->Add(waveform[12],-1);
      //waveform[11]->Add(difference_hist,-1);
      waveform[11]->Draw();
      //waveform[11]->Add(waveform[12],-1);
      difference_hist->Draw("same");
      /*for (int i = 0;i < 12;i++){
        func->SetParameters(300000,0.4,70,0,0);
        for (int j = 0;j < 10;j++){
          waveform[i]->Fit("fit","","",30,80);
        }
        if (i == 7) integrated_charge = func->GetParameter(0);
        calibrated_integrated_charge += func->GetParameter(0);
      }*/
      integrated_charge_crystal11->Fill(integrated_charge);
      calibrated_integrated_charge_crystal11->Fill(calibrated_integrated_charge);
    }
  }
  crys->SetParameters(1000,integrated_charge_crystal11->GetMean(),(integrated_charge_crystal11->GetStdDev())/4,1.1,2.1);
  //integrated_charge_crystal11->Draw();
  //integrated_charge_crystal11->Fit("crys");
  crys->SetParameters(800,calibrated_integrated_charge_crystal11->GetMean(),(calibrated_integrated_charge_crystal11->GetStdDev())/4,1.1,2.1);
  //calibrated_integrated_charge_crystal11->Draw();
  //calibrated_integrated_charge_crystal11->Fit("crys");
}