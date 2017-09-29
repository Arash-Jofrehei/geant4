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
float WF_val0[1024];
float WF_val1[1024];
float WF_val2[1024];
float WF_val5[1024];
float WF_val6[1024];
float WF_val7[1024];
float WF_val8[1024];
float WF_val9[1024];
float WF_time10[1024];
float WF_val10[1024];
float WF_val13[1024];
float WF_val14[1024];
float WF_val15[1024];
float WF_val17[1024];
float Time[18];
Float_t         x[2];
Float_t         y[2];
Int_t           nFibresOnX[2];
Int_t           nFibresOnY[2];
TBranch        *b_nFibresOnX;
TBranch        *b_nFibresOnY;
int nbins = 7;//20;
float x_start = 205;//197;
float x_end = 209;//217;
float y_start = 292.5;//284.5;
float y_end = 296.5;//304.5;
TH1F* apd_plus_electronics;
TH1F* convoluted_pulse_fscint;
TH1F* convoluted_pulse_WLS;
double fit_function(double *v,double *par)
{
  int n = int(v[0]*5);
  int shift = int(par[3]);
  return par[0]*apd_plus_electronics->GetBinContent(n+1+shift)+par[1]*convoluted_pulse_fscint->GetBinContent(n+1+shift)+par[2]*convoluted_pulse_WLS->GetBinContent(n+1+shift)+par[4];
  //return par[0]*(par[1]*convoluted_pulse_fscint->GetBinContent(n+1+shift)+(1-par[1])*convoluted_pulse_WLS->GetBinContent(n+1+shift))+par[3];
  //return par[0]*(par[1]*apd_plus_electronics->GetBinContent(n+1+shift)+(1-par[1])*convoluted_pulse_WLS->GetBinContent(n+1+shift))+par[3];
}
void fscint_contribution(){
  TCanvas *canvas = new TCanvas("fscint_contribution","fscint_contribution");
  //TFile *final = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/merged_crystal4apd.root");
  //TFile *final = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/analysis_4686.root");
  TFile *final = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/analysis_4684.root");
  //TFile *final = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/merged_crystal11.root");
  //TFile *apd_profile_waveform_bins = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/apd_profile_waveform_bins_shift.root");
  TFile *channel_profile_waveform_7bins = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/crystal4apd_channel_profile_waveform_7bins.root");
  //TFile *apd_profile_waveform_bins = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/apd_profile_waveform_bins.root");
  //TFile *apd_profile_waveform_bins = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/crystal11_profile_waveform_bins.root");
  TTree *ftree = (TTree*) final->Get("h4");
  TFile *convoluted_pulses = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/convoluted_pulses.root");
  convoluted_pulse_fscint = (TH1F*) convoluted_pulses->Get("normalized fiber scintillation convoluted with APD+electronics");
  convoluted_pulse_WLS = (TH1F*) convoluted_pulses->Get("normalized WLS+CeF3 convoluted with APD+electronics");
  apd_plus_electronics = (TH1F*) convoluted_pulses->Get("normalized APD+electronics");
  TF1 *func = new TF1("fit",fit_function,-100,200,4);
  func->SetParNames("spike on APD","fiber scintillation","WLS","shift","vertical offset");
  func->SetParLimits(0,0,100000000);
  func->SetParLimits(1,0,100000000);
  func->SetParLimits(2,0,100000000);
  func->SetParLimits(3,20,80);
  double spike_ratio = 0;
  double fscint_ratio = 0.4;
  double WLS_ratio = 0;
  double integrated_charge = 1;
  TH2F *ratio_hist2D_apd[6];
  ratio_hist2D_apd[0] = new TH2F("ratio_hist2D_apd1","nuclear counter effect contribution for APD1",nbins,x_start,x_end,nbins,y_start,y_end);
  ratio_hist2D_apd[1] = new TH2F("ratio_hist2D_apd2","nuclear counter effect contribution for APD2",nbins,x_start,x_end,nbins,y_start,y_end);
  ratio_hist2D_apd[2] = new TH2F("ratio_hist2D_apd3","nuclear counter effect contribution for APD3",nbins,x_start,x_end,nbins,y_start,y_end);
  ratio_hist2D_apd[3] = new TH2F("ratio_hist2D_apd4","nuclear counter effect contribution for APD4",nbins,x_start,x_end,nbins,y_start,y_end);
  ratio_hist2D_apd[4] = new TH2F("calibrated template amplitude (single crystal)","calibrated template amplitude (single crystal)",nbins,x_start,x_end,nbins,y_start,y_end);
  ratio_hist2D_apd[5] = new TH2F("calibrated template amplitude (3*3 matrix)","calibrated template amplitude (3*3 matrix)",nbins,x_start,x_end,nbins,y_start,y_end);
  TProfile *horizontal_strip = new TProfile("horizontal strip","aggregate WLS of 4APDs",40,x_start,x_end);
  TProfile *vertical_strip = new TProfile("vertical strip","aggregate WLS of 4APDs",40,y_start,y_end);
  bool save = false;
  TProfile *channel[18][nbins][nbins];
  string name;
  for (int a = 0;a < 18;a++){
    if (a == 3 || a == 4 || a == 11 || a == 12 || a == 16 || a == 17) continue;
    for(int i = 0;i < nbins;i++){
      for(int j = 0;j < nbins;j++){
        name = "channel";
        //name = "APD";
        name.append(to_string(a));
        name.append("_");
        name.append(to_string(i));
        name.append("_");
        name.append(to_string(j));
        const char *APD_name = name.c_str();
        if (save) channel[a][i][j] = new TProfile(APD_name,APD_name,1024,-0.1,204.7);
        //else channel[a][i][j] = (TProfile*) apd_profile_waveform_bins->Get(APD_name);
        else channel[a][i][j] = (TProfile*) channel_profile_waveform_7bins->Get(APD_name);
      }
    }
  }
  TH1F *APD_WF[4];
  APD_WF[0] = new TH1F("APD1","APD1",1024,-0.1,204.7);
  APD_WF[1] = new TH1F("APD2","APD2",1024,-0.1,204.7);
  APD_WF[2] = new TH1F("APD3","APD3",1024,-0.1,204.7);
  APD_WF[3] = new TH1F("APD4","APD4",1024,-0.1,204.7);
  ftree->SetBranchAddress("WF_val0",WF_val0);
  ftree->SetBranchAddress("WF_val1",WF_val1);
  ftree->SetBranchAddress("WF_val2",WF_val2);
  ftree->SetBranchAddress("WF_val5",WF_val5);
  ftree->SetBranchAddress("WF_val6",WF_val6);
  ftree->SetBranchAddress("WF_val7",WF_val7);
  ftree->SetBranchAddress("WF_val8",WF_val8);
  ftree->SetBranchAddress("WF_val9",WF_val9);
  ftree->SetBranchAddress("WF_time10",WF_time10);
  ftree->SetBranchAddress("WF_val10",WF_val10);
  ftree->SetBranchAddress("WF_val13",WF_val13);
  ftree->SetBranchAddress("WF_val14",WF_val14);
  ftree->SetBranchAddress("WF_val15",WF_val15);
  ftree->SetBranchAddress("WF_val17",WF_val17);
  ftree->SetBranchAddress("Time",Time);
  ftree->SetBranchAddress("x", x);
  ftree->SetBranchAddress("y", y);
  ftree->SetBranchAddress("nFibresOnX", nFibresOnX, &b_nFibresOnX);
  ftree->SetBranchAddress("nFibresOnY", nFibresOnY, &b_nFibresOnY);
  
  Long64_t nentries = ftree->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  if (save){
    for (Long64_t jentry=0; jentry<nentries;jentry++){
      if (jentry%500 == 0) cout << jentry << endl;
      Long64_t ientry = ftree->LoadTree(jentry);
      nb = ftree->GetEntry(jentry);   nbytes += nb;
      if (x[1]<x_start||x[1]>=x_end||y[1]<y_start||y[1]>=y_end) continue;
      int x_index = int((nbins/(x_end-x_start))*(x[1]-x_start));
      int y_index = int((nbins/(y_end-y_start))*(y[1]-y_start));
      for (int sample = 0;sample < 1024;sample++){
        channel[0][x_index][y_index]->Fill(30+WF_time10[sample]-Time[0],WF_val0[sample]);
        channel[1][x_index][y_index]->Fill(30+WF_time10[sample]-Time[1],WF_val1[sample]);
        channel[2][x_index][y_index]->Fill(30+WF_time10[sample]-Time[2],WF_val2[sample]);
        channel[5][x_index][y_index]->Fill(30+WF_time10[sample]-Time[5],WF_val5[sample]);
        channel[6][x_index][y_index]->Fill(30+WF_time10[sample]-Time[6],WF_val6[sample]);
        channel[7][x_index][y_index]->Fill(30+WF_time10[sample]-Time[7],WF_val7[sample]);
        channel[8][x_index][y_index]->Fill(30+WF_time10[sample]-Time[8],WF_val8[sample]);
        channel[9][x_index][y_index]->Fill(30+WF_time10[sample]-Time[9],WF_val9[sample]);
        channel[10][x_index][y_index]->Fill(30+WF_time10[sample]-Time[10],WF_val10[sample]);
        channel[13][x_index][y_index]->Fill(30+WF_time10[sample]-Time[13],WF_val13[sample]);
        channel[14][x_index][y_index]->Fill(30+WF_time10[sample]-Time[14],WF_val14[sample]);
        channel[15][x_index][y_index]->Fill(30+WF_time10[sample]-Time[15],WF_val15[sample]);
      }
    }
  }
  /*Long64_t nentries = ftree->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  int counter = 0;
  float WLS_amp = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++){
    if (jentry%500 == 0) cout << jentry << endl;
    Long64_t ientry = ftree->LoadTree(jentry);
    nb = ftree->GetEntry(jentry);   nbytes += nb;
    if (x[1]<x_start||x[1]>=x_end||y[1]<y_start||y[1]>=y_end) continue;
    for (int sample = 0;sample < 1024;sample++){
      APD_WF[0]->SetBinContent(1+int(5*(30+WF_time10[sample]-Time[6])),0.954843*WF_val6[sample]);
      APD_WF[1]->SetBinContent(1+int(5*(30+WF_time10[sample]-Time[7])),0.929462*WF_val7[sample]);
      APD_WF[2]->SetBinContent(1+int(5*(30+WF_time10[sample]-Time[8])),1.1492*WF_val8[sample]);
      APD_WF[3]->SetBinContent(1+int(5*(30+WF_time10[sample]-Time[9])),0.993394*WF_val9[sample]);
    }
    WLS_amp = 0;
    for (int a = 0;a < 4;a++){
      func->SetParameters(50000,1,1,30,0,0);
      for (int k = 0;k < 20;k++) APD_WF[a]->Fit("fit","Q","",25,60);
      WLS_amp += func->GetParameter(2);
      APD_WF[a]->Reset();
    }
    //counter+=1;
    //horizontal_strip->Fill(x[1],WLS_amp);
    //if (counter>3) break;
    //vertical_strip->Fill(y[1],WLS_amp);
  }*/
  //TFile *apd_profile_waveform_bins = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/apd_profile_waveform_bins.root","recreate");
  //TFile *apd_profile_waveform_bins = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/apd_profile_waveform_bins_shift.root","recreate");
  //TFile *xtal11_profile_waveform_bins = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/crystal11_profile_waveform_bins.root","recreate");
  //TFile *channel_profile_waveform_7bins = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/crystal4apd_channel_profile_waveform_7bins.root","recreate");
  func->FixParameter(1,0);
  for (int i = 0;i < nbins;i++){
    for (int j = 0;j < nbins;j++){
      cout << endl << i+1 << "    " << j+1 << endl;
      //xtal11_profile_waveform_bins->cd();
      //apd[0][i][j]->Write();
      for (int a = 0;a < 18;a++){
        if (a == 3 || a == 4 || a == 11 || a == 12 || a == 16 || a == 17) continue;
        //channel_profile_waveform_7bins->cd();
        //channel[a][i][j]->Write();
        func->SetParameters(50000,0,50000,30,0,0);
        for (int k = 0;k < 100;k++){
          //channel[a][i][j]->Fit("fit","Q","",25,60);
        }
        integrated_charge = func->GetParameter(0)+func->GetParameter(1)+func->GetParameter(2);
        spike_ratio = func->GetParameter(0)/integrated_charge;
        fscint_ratio = func->GetParameter(1)/integrated_charge;
        WLS_ratio = func->GetParameter(2)/integrated_charge;
        cout << "channel" << a << ":    " << func->GetParameter(2) << endl;
        //if (a==0) ratio_hist2D_apd[0]->SetBinContent(i+1,j+1,func->GetParameter(2));
        //else ratio_hist2D_apd[0]->SetBinContent(i+1,j+1,func->GetParameter(2)+ratio_hist2D_apd[0]->GetBinContent(i+1,j+1));
      }
    }
  }
  for (int i = 0;i < nbins;i++){
    for (int j = 0;j < nbins;j++){
      //cout << endl << i+1 << "    " << j+1 << "      " << ratio_hist2D_apd[2]->GetBinContent(i+1,j+1);
    }
  }
  //cout << endl;/*
  func->SetParameters(150000,0,150000,45,0,0);
  for (int i = 0;i < 100; i++){
    channel[6][3][3]->SetLineWidth(0);
    //channel[0][0][6]->Smooth(1);
    channel[6][3][3]->Fit("fit","","",25,95);
  }
  integrated_charge = func->GetParameter(0)+func->GetParameter(1)+func->GetParameter(2);
  spike_ratio = func->GetParameter(0)/integrated_charge;
  fscint_ratio = func->GetParameter(1)/integrated_charge;
  WLS_ratio = func->GetParameter(2)/integrated_charge;
  cout << "spike:      " << spike_ratio*100 << endl;
  cout << "fscint:      " << fscint_ratio*100 << endl;
  cout << "WLS:      " << WLS_ratio*100 << endl;
  TH1F *spike_hist = new TH1F("spike hist","nuclear counter effect contribution",1024,-0.1,204.7);
  TH1F *fscint_hist = new TH1F("fscint hist","fiber scintillation contribution",1024,-0.1,204.7);
  TH1F *WLS_hist = new TH1F("WLS hist","WLS contribution",1024,-0.1,204.7);
  for (int sample = 0;sample<1024;sample++){
    spike_hist->SetBinContent(sample+1,(apd_plus_electronics->GetBinContent(sample+1+int(func->GetParameter(3))))*func->GetParameter(0));
    fscint_hist->SetBinContent(sample+1,(convoluted_pulse_fscint->GetBinContent(sample+1+int(func->GetParameter(3))))*func->GetParameter(1));
    WLS_hist->SetBinContent(sample+1,(convoluted_pulse_WLS->GetBinContent(sample+1+int(func->GetParameter(3))))*func->GetParameter(2));
  }
  spike_hist->SetLineColor(2);
  fscint_hist->SetLineColor(3);
  WLS_hist->SetLineColor(4);
  spike_hist->Draw("same");
  fscint_hist->Draw("same");
  WLS_hist->Draw("same");
  /*
  ratio_hist2D_apd[3]->Draw("colz");
  cout << (ratio_hist2D_apd[1]->GetBinContent(5,3) - ratio_hist2D_apd[1]->GetBinContent(1,3))/(ratio_hist2D_apd[0]->GetBinContent(3,5) - ratio_hist2D_apd[0]->GetBinContent(3,1)) << endl;
  cout << (ratio_hist2D_apd[1]->GetBinContent(5,3) - ratio_hist2D_apd[1]->GetBinContent(1,3))/(ratio_hist2D_apd[1]->GetBinContent(5,3) - ratio_hist2D_apd[1]->GetBinContent(1,3)) << endl;
  cout << (ratio_hist2D_apd[1]->GetBinContent(5,3) - ratio_hist2D_apd[1]->GetBinContent(1,3))/(ratio_hist2D_apd[2]->GetBinContent(1,3) - ratio_hist2D_apd[2]->GetBinContent(5,3)) << endl;
  cout << (ratio_hist2D_apd[1]->GetBinContent(5,3) - ratio_hist2D_apd[1]->GetBinContent(1,3))/(ratio_hist2D_apd[3]->GetBinContent(3,1) - ratio_hist2D_apd[3]->GetBinContent(3,5)) << endl;
  int single_sum = 0;
  for (int i = 0;i < nbins;i++){
    for (int j = 0;j < nbins;j++){
      single_sum = 0;
      for (int k = 0;k < 4;k++){
        if (k==3) single_sum += ((ratio_hist2D_apd[2]->GetBinContent(ratio_hist2D_apd[2]->GetMaximumBin()) - ratio_hist2D_apd[2]->GetBinContent(ratio_hist2D_apd[2]->GetMinimumBin()))/(ratio_hist2D_apd[k]->GetBinContent(ratio_hist2D_apd[k]->GetMaximumBin()) - ratio_hist2D_apd[k]->GetBinContent(ratio_hist2D_apd[2]->GetMinimumBin())))*ratio_hist2D_apd[k]->GetBinContent(i+1,j+1);
        else single_sum += ((ratio_hist2D_apd[2]->GetBinContent(ratio_hist2D_apd[2]->GetMaximumBin()) - ratio_hist2D_apd[2]->GetBinContent(ratio_hist2D_apd[2]->GetMinimumBin()))/(ratio_hist2D_apd[k]->GetBinContent(ratio_hist2D_apd[k]->GetMaximumBin()) - ratio_hist2D_apd[k]->GetBinContent(ratio_hist2D_apd[k]->GetMinimumBin())))*ratio_hist2D_apd[k]->GetBinContent(i+1,j+1);
      }
      ratio_hist2D_apd[4]->SetBinContent(i+1,j+1,single_sum);
    }
  }*/
  //horizontal_strip->Draw();
  //vertical_strip->Draw();
  //ratio_hist2D_apd[0]->Draw("colz");
  //apd_plus_electronics->Draw();
  //convoluted_pulse_fscint->Draw("same");
  //convoluted_pulse_WLS->Draw("same");
}