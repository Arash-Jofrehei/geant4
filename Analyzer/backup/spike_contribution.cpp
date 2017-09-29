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
int nbins = 20;
float x_start = 214;//197;
float x_end = 234;//217;
float y_start = 284.5;
float y_end = 304.5;
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
void spike_contribution(){
  TCanvas *canvas = new TCanvas("fscint_contribution","fscint_contribution");
  //TFile *final = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/merged_crystal4apd.root");
  //TFile *final = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/analysis_4686.root");
  //TFile *final = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/analysis_4684.root");
  TFile *final = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/merged_crystal11.root");
  //TFile *apd_profile_waveform_bins = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/apd_profile_waveform_bins_shift.root");
  //TFile *apd_profile_waveform_bins = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/apd_profile_waveform_bins.root");
  TFile *apd_profile_waveform_bins = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/crystal11_profile_waveform_bins.root");
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
  ratio_hist2D_apd[0] = new TH2F("ratio_hist2D_apd1","WLS contribution for crystal 11",nbins,x_start,x_end,nbins,y_start,y_end);
  ratio_hist2D_apd[1] = new TH2F("ratio_hist2D_apd2","nuclear counter effect contribution for APD2",nbins,x_start,x_end,nbins,y_start,y_end);
  ratio_hist2D_apd[2] = new TH2F("ratio_hist2D_apd3","nuclear counter effect contribution for APD3",nbins,x_start,x_end,nbins,y_start,y_end);
  ratio_hist2D_apd[3] = new TH2F("ratio_hist2D_apd4","nuclear counter effect contribution for APD4",nbins,x_start,x_end,nbins,y_start,y_end);
  ratio_hist2D_apd[4] = new TH2F("calibrated template amplitude (single crystal)","calibrated template amplitude (single crystal)",nbins,x_start,x_end,nbins,y_start,y_end);
  ratio_hist2D_apd[5] = new TH2F("calibrated template amplitude (3*3 matrix)","calibrated template amplitude (3*3 matrix)",nbins,x_start,x_end,nbins,y_start,y_end);
  bool save = false;
  TProfile *apd[4][nbins][nbins];
  string name;
  for (int a = 0;a < 4;a++){
    if (a != 0) continue;
    for(int i = 0;i < nbins;i++){
      for(int j = 0;j < nbins;j++){
        name = "crystal11";
        //name = "APD";
        //name.append(to_string(a+1));
        name.append("_");
        name.append(to_string(i));
        name.append("_");
        name.append(to_string(j));
        const char *APD_name = name.c_str();
        if (save) apd[a][i][j] = new TProfile(APD_name,APD_name,1024,-0.1,204.7);
        else apd[a][i][j] = (TProfile*) apd_profile_waveform_bins->Get(APD_name);
      }
    }
  }
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
  ftree->SetBranchAddress("nFibresOnY", nFibresOnY, &b_nFibresOnY);/*
  Long64_t nentries = ftree->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  if (save){
    for (Long64_t jentry=0; jentry<nentries;jentry++){
      Long64_t ientry = ftree->LoadTree(jentry);
      nb = ftree->GetEntry(jentry);   nbytes += nb;
      if (x[1]<x_start||x[1]>=x_end||y[1]<y_start||y[1]>=y_end) continue;
      int x_index = int((nbins/(x_end-x_start))*(x[1]-x_start));
      int y_index = int((nbins/(y_end-y_start))*(y[1]-y_start));
      for (int sample = 0;sample < 1009;sample++){
        apd[0][x_index][y_index]->Fill(WF_time6[sample]-28.6,WF_val6[sample]);
        apd[1][x_index][y_index]->Fill(WF_time7[sample]-28.6,WF_val7[sample]);
        apd[2][x_index][y_index]->Fill(WF_time8[sample]-28.6,WF_val8[sample]);
        apd[3][x_index][y_index]->Fill(WF_time9[sample]-28.6,WF_val9[sample]);
      }
      for (int sample = 866;sample < 1024;sample++){
        for (int a = 0;a < 4;a++){
          apd[a][x_index][y_index]->Fill(sample/5.0,0);
        }
      }
    }
  }*/
  /*
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
        apd[0][x_index][y_index]->Fill(30+WF_time10[sample]-Time[10],WF_val10[sample]);
        //apd[0][x_index][y_index]->Fill(30+WF_time10[sample]-Time[1],WF_val1[sample]);
        //apd[1][x_index][y_index]->Fill(30+WF_time10[sample]-Time[5],WF_val5[sample]);
        //apd[2][x_index][y_index]->Fill(30+WF_time10[sample]-Time[10],WF_val10[sample]);
        //apd[3][x_index][y_index]->Fill(30+WF_time10[sample]-Time[14],WF_val14[sample]);
        //apd[0][x_index][y_index]->Fill(30+WF_time10[sample]-Time[6],WF_val6[sample]);
        //apd[1][x_index][y_index]->Fill(30+WF_time10[sample]-Time[7],WF_val7[sample]);
        //apd[2][x_index][y_index]->Fill(30+WF_time10[sample]-Time[8],WF_val8[sample]);
        //apd[3][x_index][y_index]->Fill(30+WF_time10[sample]-Time[9],WF_val9[sample]);
      }
    }
  }*/
  //TFile *apd_profile_waveform_bins = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/apd_profile_waveform_bins.root","recreate");
  //TFile *apd_profile_waveform_bins = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/apd_profile_waveform_bins_shift.root","recreate");
  //TFile *xtal11_profile_waveform_bins = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/crystal11_profile_waveform_bins.root","recreate");
  for (int i = 0;i < nbins;i++){
    for (int j = 0;j < nbins;j++){
      cout << endl << i+1 << "    " << j+1 << endl;
      //xtal11_profile_waveform_bins->cd();
      //apd[0][i][j]->Write();
      for (int a = 0;a < 4;a++){
        //apd_profile_waveform_bins->cd();
        //apd[a][i][j]->Write();
        func->SetParameters(50000,1,0.0,30,0,0);
        for (int k = 0;k < 100;k++){
          if (a==0) apd[a][i][j]->Fit("fit","Q","",25,60);
        }
        integrated_charge = func->GetParameter(0)+func->GetParameter(1)+func->GetParameter(2);
        spike_ratio = func->GetParameter(0)/integrated_charge;
        fscint_ratio = func->GetParameter(1)/integrated_charge;
        WLS_ratio = func->GetParameter(2)/integrated_charge;
        ratio_hist2D_apd[a]->SetBinContent(i+1,j+1,func->GetParameter(2));
      }
    }
  }
  for (int i = 0;i < nbins;i++){
    for (int j = 0;j < nbins;j++){
      //cout << endl << i+1 << "    " << j+1 << "      " << ratio_hist2D_apd[2]->GetBinContent(i+1,j+1);
    }
  }
  //cout << endl;
  func->SetParameters(50000,1,0.0,30,0,0);
  //func->FixParameter(1,0);
  //func->SetParLimits(3,20,50);
  for (int i = 0;i < 100; i++){
    apd[0][10][10]->SetLineWidth(0);
    //apd[0][3][6]->Smooth(1);
    apd[0][10][10]->Fit("fit","","",25,60);
  }
  integrated_charge = func->GetParameter(0)+func->GetParameter(1)+func->GetParameter(2);
  spike_ratio = func->GetParameter(0)/integrated_charge;
  fscint_ratio = func->GetParameter(1)/integrated_charge;
  WLS_ratio = func->GetParameter(2)/integrated_charge;
  cout << "spike:      " << spike_ratio*100 << endl;
  cout << "fscint:      " << fscint_ratio*100 << endl;
  cout << "WLS:      " << WLS_ratio*100 << endl;
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
  ratio_hist2D_apd[0]->Draw("colz");
  //apd_plus_electronics->Draw();
  //convoluted_pulse_fscint->Draw("same");
  //convoluted_pulse_WLS->Draw("same");
}