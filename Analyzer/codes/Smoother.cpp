#include <iostream>
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
using namespace std;
float WF_val0[1024];
float WF_val1[1024];
float WF_val2[1024];
float WF_val3[1024];
float WF_val4[1024];
float WF_val5[1024];
float WF_val6[1024];
float WF_val7[1024];
float WF_val8[1024];
float WF_val9[1024];
float WF_time10[1024];
float WF_val10[1024];
float WF_val11[1024];
float WF_val12[1024];
float WF_val13[1024];
float WF_val14[1024];
float WF_val15[1024];
float WF_val16[1024];
float WF_val17[1024];
float Time[18];
Float_t         x[2];
Float_t         y[2];
float median_integrated_charge[18];
float median_xbar;
float median_ybar;
float median_x[2];
float median_y[2];
int median_nFibersX[2];
int median_nFibersY[2];
Int_t           nFibresOnX[2];
Int_t           nFibresOnY[2];
TBranch        *b_nFibresOnX;
TBranch        *b_nFibresOnY;
                                                    bool we_are_impinging_xtal11 = true;

float channel_weights_crystal4APD[18] = {2.41098,1.62535,1.49631,0,0,0.470364,0.954885,0.931846,1.15933,0.983237,1.3819,0,0,0.452602,1.86884,1.41348,0,0};
float channel_weights_crystal11[18] = {0,1,1,1,0,0,1,0.35,0.35,0.35,0.35,1,0,0,1,1,1,0};
float crystal_x[18] = {190.0,207.0,224.0,241.0,258.0,190.0,207.0,207.0,207.0,207.0,224.0,241.0,258.0,190.0,207.0,224.0,241.0,258.0};
float crystal_y[18] = {311.5,311.5,311.5,311.5,311.5,294.5,294.5,294.5,294.5,294.5,294.5,294.5,294.5,277.5,277.5,277.5,277.5,277.5};

float integrated_charge(TH1F* histogram){
  float baseline = 0;
  for (int i = 15; i<100;i++) baseline += histogram->GetBinContent(i+1)/85.0;
  int start = 0;
  for (int i = 100; i<500;i++){
    if (histogram->GetBinContent(i+1) - baseline > 15){
      start = i;
      break;
    }
  }
  int end = 0;
  for (int i = start; i<1024;i++){
    if (histogram->GetBinContent(i+1) - baseline < 15){
      end = i;
      break;
    }
  }
  //start = 35*5;
  end = start + 500;
  float charge = 0;
  for (int i = start; i<end;i++) charge += histogram->GetBinContent(i+1) - baseline;
  return charge;
}

TH1F *smoothed_histogram = new TH1F("smoothed histogram","smoothed histogram",1024,-0.1,204.7);

TH1F *smoothed(TH1F* histogram,int radius = 10,int iteration = 1,bool saveMAX = true){
  float interval[2*radius+1];
  float buffer = -1000;
  for (int epoch = 0; epoch < iteration ;epoch++){
    for (int i = 0; i < 1024; i++){
      if (i < radius || i > 1024-radius || (saveMAX == true && TMath::Abs(i+1 - histogram->GetMaximumBin()) < radius)) smoothed_histogram->SetBinContent(i+1,histogram->GetBinContent(i+1));
      else{
        for (int j = 0; j < 2*radius+1 ;j++){
          interval[j] = histogram->GetBinContent(i+1+j-radius);
        }
        for (int j = 0;j < 2*radius+1 ;j++){
          for (int k = j;k < 2*radius+1 ;k++){
            if (interval[j] > interval[k]){
              buffer = interval[j];
              interval[j] = interval[k];
              interval[k] = buffer;
            }
          }
        }
        smoothed_histogram->SetBinContent(i+1,interval[radius]);
      }
    }
  }
  return smoothed_histogram;
}

void Smoother(){
  TCanvas *canvas = new TCanvas("Smoother","Smoother");
  TFile *final;
  if (we_are_impinging_xtal11){
    //final = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/merged_crystal11.root","update");
    final = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/analysis_4683.root","update");
  }
  else{
    //final = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/merged_crystal4apd.root","update");
    final = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/analysis_4684.root","update");
  }
  TTree *ftree = (TTree*) final->Get("h4");
  ftree->SetBranchAddress("WF_val0",WF_val0);
  ftree->SetBranchAddress("WF_val1",WF_val1);
  ftree->SetBranchAddress("WF_val2",WF_val2);
  ftree->SetBranchAddress("WF_val3",WF_val3);
  ftree->SetBranchAddress("WF_val4",WF_val4);
  ftree->SetBranchAddress("WF_val5",WF_val5);
  ftree->SetBranchAddress("WF_val6",WF_val6);
  ftree->SetBranchAddress("WF_val7",WF_val7);
  ftree->SetBranchAddress("WF_val8",WF_val8);
  ftree->SetBranchAddress("WF_val9",WF_val9);
  ftree->SetBranchAddress("WF_time10",WF_time10);
  ftree->SetBranchAddress("WF_val10",WF_val10);
  ftree->SetBranchAddress("WF_val11",WF_val11);
  ftree->SetBranchAddress("WF_val12",WF_val12);
  ftree->SetBranchAddress("WF_val13",WF_val13);
  ftree->SetBranchAddress("WF_val14",WF_val14);
  ftree->SetBranchAddress("WF_val15",WF_val15);
  ftree->SetBranchAddress("WF_val16",WF_val16);
  ftree->SetBranchAddress("WF_val17",WF_val17);
  ftree->SetBranchAddress("Time",Time);
  ftree->SetBranchAddress("x", x);
  ftree->SetBranchAddress("y", y);
  ftree->SetBranchAddress("nFibresOnX", nFibresOnX, &b_nFibresOnX);
  ftree->SetBranchAddress("nFibresOnY", nFibresOnY, &b_nFibresOnY);
  TTree *median = new TTree("median","median");
  median->Branch("median_integrated_charge",&median_integrated_charge,"median_integrated_charge[18]/F");
  median->Branch("median_xbar",&median_xbar,"median_xbar/F");
  median->Branch("median_ybar",&median_ybar,"median_ybar/F");
  median->Branch("median_x",&median_x,"median_x[2]/F");
  median->Branch("median_y",&median_y,"median_y[2]/F");
  median->Branch("median_nFibersX",median_nFibersX,"median_nFibersX[2]/I");
  median->Branch("median_nFibersY",median_nFibersY,"median_nFibersY[2]/I");
  TH1F *waveform[18];
  waveform[0] = new TH1F("waveform0","waveform crystal 1",1024,-0.1,204.7);
  waveform[1] = new TH1F("waveform1","waveform crystal 2",1024,-0.1,204.7);
  waveform[2] = new TH1F("waveform2","waveform crystal 3",1024,-0.1,204.7);
  waveform[3] = new TH1F("waveform3","waveform crystal 4",1024,-0.1,204.7);
  waveform[4] = new TH1F("waveform4","waveform crystal 5",1024,-0.1,204.7);
  waveform[5] = new TH1F("waveform5","waveform crystal 6",1024,-0.1,204.7);
  waveform[6] = new TH1F("waveform6","waveform crystal4APD 1",1024,-0.1,204.7);
  waveform[7] = new TH1F("waveform7","waveform crystal4APD 2",1024,-0.1,204.7);
  waveform[8] = new TH1F("waveform8","waveform crystal4APD 3",1024,-0.1,204.7);
  waveform[9] = new TH1F("waveform9","waveform crystal4APD 4",1024,-0.1,204.7);
  waveform[10] = new TH1F("waveform10","waveform crystal 11",1024,-0.1,204.7);
  waveform[11] = new TH1F("waveform11","waveform crystal 12",1024,-0.1,204.7);
  waveform[12] = new TH1F("waveform12","waveform crystal 13",1024,-0.1,204.7);
  waveform[13] = new TH1F("waveform13","waveform crystal 14",1024,-0.1,204.7);
  waveform[14] = new TH1F("waveform14","waveform crystal 15",1024,-0.1,204.7);
  waveform[15] = new TH1F("waveform15","waveform crystal 16",1024,-0.1,204.7);
  waveform[16] = new TH1F("waveform16","waveform crystal 17",1024,-0.1,204.7);
  waveform[17] = new TH1F("waveform17","waveform crystal 18",1024,-0.1,204.7);
  for (int i = 0; i < 18; i++) waveform[i]->SetLineWidth(0);
  
  int counter = 0;
  Long64_t nentries = ftree->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++){
    if (jentry%500 == 0) cout << jentry << endl;
    Long64_t ientry = ftree->LoadTree(jentry);
    nb = ftree->GetEntry(jentry);   nbytes += nb;
    median_x[0] = x[0];
    median_x[1] = x[1];
    median_y[0] = y[0];
    median_y[1] = y[1];
    median_nFibersX[0] = nFibresOnX[0];
    median_nFibersX[1] = nFibresOnX[1];
    median_nFibersY[0] = nFibresOnY[0];
    median_nFibersY[1] = nFibresOnY[1];
    for (int i=0;i<18;i++) median_integrated_charge[i] = 0;
    median_xbar = 0;
    median_ybar = 0;
    //if (TMath::Abs(x[1]-207)>5||TMath::Abs(y[1]-294.5)>5||nFibresOnX[1]>2||nFibresOnY[1]>2){
    if (nFibresOnX[1]>2||nFibresOnY[1]>2){
      median->Fill();
      continue;
    }
    counter += 1;
    for (int i=0;i<18;i++) waveform[i]->Reset();
    for (int sample = 0;sample < 1024;sample++){
      waveform[0]->Fill(WF_time10[sample],WF_val0[sample]);
      waveform[1]->Fill(WF_time10[sample],WF_val1[sample]);
      waveform[2]->Fill(WF_time10[sample],WF_val2[sample]);
      waveform[3]->Fill(WF_time10[sample],WF_val3[sample]);
      waveform[4]->Fill(WF_time10[sample],WF_val4[sample]);
      waveform[5]->Fill(WF_time10[sample],WF_val5[sample]);
      waveform[6]->Fill(WF_time10[sample],WF_val6[sample]);
      waveform[7]->Fill(WF_time10[sample],WF_val7[sample]);
      waveform[8]->Fill(WF_time10[sample],WF_val8[sample]);
      waveform[9]->Fill(WF_time10[sample],WF_val9[sample]);
      waveform[10]->Fill(WF_time10[sample],WF_val10[sample]);
      waveform[11]->Fill(WF_time10[sample],WF_val11[sample]);
      waveform[12]->Fill(WF_time10[sample],WF_val12[sample]);
      waveform[13]->Fill(WF_time10[sample],WF_val13[sample]);
      waveform[14]->Fill(WF_time10[sample],WF_val14[sample]);
      waveform[15]->Fill(WF_time10[sample],WF_val15[sample]);
      waveform[16]->Fill(WF_time10[sample],WF_val16[sample]);
      waveform[17]->Fill(WF_time10[sample],WF_val17[sample]);
    }
    for (int i=0;i<18;i++) median_integrated_charge[i] = integrated_charge(smoothed(smoothed(waveform[i],15,2,false),150,1));
    /*float x_ratio = channel_weights_crystal4APD[10]*median_integrated_charge[10]/(channel_weights_crystal4APD[5]*median_integrated_charge[5]);
    median_xbar = 207.0 + 17 *TMath::Log10(x_ratio) / (6-TMath::Log10(2 + x_ratio + 1/x_ratio));
    float y_ratio = channel_weights_crystal4APD[1]*median_integrated_charge[1]/(channel_weights_crystal4APD[14]*median_integrated_charge[14]);
    median_ybar = 294.5 + 17 *TMath::Log10(y_ratio) / (6-TMath::Log10(2 + y_ratio + 1/y_ratio));
    if (x_ratio < 0 || y_ratio < 0){
      median_xbar = -100;
      median_ybar = -100;
    }*/
    float energy_sum;
    float median_position_weight[18];
    float weight_sum;
    if (we_are_impinging_xtal11){
      energy_sum = 0;
      for (int i = 0;i < 18;i++) energy_sum += channel_weights_crystal11[i]*median_integrated_charge[i];
      for (int i = 0;i < 18;i++){
        if (channel_weights_crystal11[i]*median_integrated_charge[i]/energy_sum < 0.001){
          median_position_weight[i] = 0;
          continue;
        }
        median_position_weight[i] = 3 + TMath::Log10(channel_weights_crystal11[i]*median_integrated_charge[i]/energy_sum);
      }
      weight_sum = 0;
      for (int i = 0;i < 18;i++) weight_sum += median_position_weight[i];
      median_xbar = 0;
      median_ybar = 0;
      for (int i = 0;i < 18;i++) median_xbar += median_position_weight[i] * crystal_x[i] / weight_sum;
      for (int i = 0;i < 18;i++) median_ybar += median_position_weight[i] * crystal_y[i] / weight_sum;
    }
    else{
      energy_sum = 0;
      for (int i = 0;i < 18;i++) energy_sum += channel_weights_crystal4APD[i]*median_integrated_charge[i];
      for (int i = 0;i < 18;i++){
        if (channel_weights_crystal4APD[i]*median_integrated_charge[i]/energy_sum < 0.001){
          median_position_weight[i] = 0;
          continue;
        }
        median_position_weight[i] = 3 + TMath::Log10(channel_weights_crystal4APD[i]*median_integrated_charge[i]/energy_sum);
      }
      weight_sum = 0;
      for (int i = 0;i < 18;i++) weight_sum += median_position_weight[i];
      median_xbar = 0;
      median_ybar = 0;
      for (int i = 0;i < 18;i++) median_xbar += median_position_weight[i] * crystal_x[i] / weight_sum;
      for (int i = 0;i < 18;i++) median_ybar += median_position_weight[i] * crystal_y[i] / weight_sum;
    }
    /*if (counter == 2){
      cout << x[1]-207 << "    " << y[1]-294.5 << endl;
      //smoothed(smoothed(waveform[6],15,2,false),150,1)->Draw();
      cout << integrated_charge(smoothed(smoothed(waveform[6],15,2,false),150,1)) << endl;
      break;
    }*/
    median->Fill();
  }
  median->Write();
}