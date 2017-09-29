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
float WF_time6[1024];
float WF_val6[1024];
float WF_time7[1024];
float WF_val7[1024];
float WF_time8[1024];
float WF_val8[1024];
float WF_time9[1024];
float WF_val9[1024];
Float_t         x[2];
Float_t         y[2];
Int_t           nFibresOnX[2];
Int_t           nFibresOnY[2];
TBranch        *b_nFibresOnX;
TBranch        *b_nFibresOnY;
int nbins = 20;
float x_start = 197;
float x_end = 217;
float y_start = 284.5;
float y_end = 304.5;
TH1F* convoluted_pulse_fscint;
TH1F* convoluted_pulse_WLS;
double fit_function(double *v,double *par)
{
  int n = int(v[0]*5);
  int WLS_shift = int(par[3]);
  int fscint_shift = int(par[4]);
  return par[0]*convoluted_pulse_fscint->GetBinContent(n+1+fscint_shift) + par[1]*convoluted_pulse_WLS->GetBinContent(n+1+WLS_shift) + par[2];
}
void fscint_contribution(){
  TCanvas *canvas = new TCanvas("fscint_contribution","fscint_contribution");
  TFile *final = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/merged_crystal4apd.root");
  //TFile *final = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/analysis_4684.root");
  TTree *ftree = (TTree*) final->Get("h4");
  TFile *convoluted_pulses = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/convoluted_pulses.root");
  convoluted_pulse_fscint = (TH1F*) convoluted_pulses->Get("convoluted pulse fiber scintillation");
  convoluted_pulse_WLS = (TH1F*) convoluted_pulses->Get("convoluted pulse WLS");
  TF1 *func = new TF1("fit",fit_function,0,132,2);
  func->SetParNames("fscint constant","WLS constant","adding constant","WLS shift","fscint shift");
  func->SetParLimits(0,0,1000);
  func->SetParLimits(1,0,1000);
  double WLS_ratio = 0.5; //    WLS/(WLS+fscint)
  double c_fscint = 1;
  double c_WLS = 1;
  TH2F *apd1_ratio_hist2D = new TH2F("apd1_ratio_hist2D","apd1_ratio_hist2D",nbins,x_start,x_end,nbins,y_start,y_end);
  /*TH2F ***apd1;
  string name;
  for(int i = 0;i < nbins;i++){
    for(int j = 0;j < nbins;j++){
      name = "APD1_";
      name.append(to_string(i));
      name.append("_");
      name.append(to_string(j));
      const char *APD_name = name.c_str();
      apd1[i][j] = new TH2F(APD_name,APD_name,1024,-0.1,204.7,124000,-44000,80000);
    }
  }*/
  TH2F *apd1 = new TH2F("apd1","apd1",1024,-0.1,204.7,124000,-44000,80000);
  TH2F *apd2 = new TH2F("apd2","apd2",1024,-0.1,204.7,124000,-44000,80000);
  TH2F *apd3 = new TH2F("apd3","apd3",1024,-0.1,204.7,124000,-44000,80000);
  TH2F *apd4 = new TH2F("apd4","apd4",1024,-0.1,204.7,124000,-44000,80000);
  ftree->SetBranchAddress("WF_time6",WF_time6);
  ftree->SetBranchAddress("WF_val6",WF_val6);
  ftree->SetBranchAddress("WF_time7",WF_time7);
  ftree->SetBranchAddress("WF_val7",WF_val7);
  ftree->SetBranchAddress("WF_time8",WF_time8);
  ftree->SetBranchAddress("WF_val8",WF_val8);
  ftree->SetBranchAddress("WF_time9",WF_time9);
  ftree->SetBranchAddress("WF_val9",WF_val9);
  ftree->SetBranchAddress("x", x);
  ftree->SetBranchAddress("y", y);
  ftree->SetBranchAddress("nFibresOnX", nFibresOnX, &b_nFibresOnX);
  ftree->SetBranchAddress("nFibresOnY", nFibresOnY, &b_nFibresOnY);
  Long64_t nentries = ftree->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  bool empty = true;
  for (int i = 0;i < nbins;i++){
    for (int j = 0;j < nbins;j++){
      apd1->Reset();3
      //apd2->Reset();
      //apd3->Reset();
      //apd4->Reset();
      empty = true;
      cout << endl << i << "    " << j << endl;
      for (Long64_t jentry=0; jentry<nentries;jentry++){
        Long64_t ientry = ftree->LoadTree(jentry);
        nb = ftree->GetEntry(jentry);   nbytes += nb;
        if (x[1]<x_start||x[1]>=x_end||y[1]<y_start||y[1]>=y_end) continue;
        int x_index = int((nbins/(x_end-x_start))*(x[1]-x_start));
        int y_index = int((nbins/(y_end-y_start))*(y[1]-y_start));
        if (x_index != i || y_index != j) continue;
        empty = false;
        for (int sample = 0;sample < 1009;sample++){
          apd1->Fill(WF_time6[sample]-28.6,20*WF_val6[sample]);
          //apd2->Fill(WF_time7[sample]-28.6,20*WF_val7[sample]);
          //apd3->Fill(WF_time8[sample]-28.6,20*WF_val8[sample]);
          //apd4->Fill(WF_time9[sample]-28.6,20*WF_val9[sample]);
        }
        for (int sample = 866;sample < 1024;sample++){
          apd1->Fill(sample/5.0,0);//-1000);
          //apd2->Fill(sample/5.0,0);//-1000);
          //apd3->Fill(sample/5.0,0);//-1000);
          //apd4->Fill(sample/5.0,0);//-1000);
        }
      }
      if (empty) continue;
      TProfile *apd1_prof = apd1->ProfileX("apd1_prof");
      //TProfile *apd2_prof = apd2->ProfileX("apd2_prof");
      //TProfile *apd3_prof = apd3->ProfileX("apd3_prof");
      //TProfile *apd4_prof = apd4->ProfileX("apd4_prof");
      TH1D *apd1_1D = apd1_prof->ProjectionX();
      //TH1D *apd2_1D = apd2_prof->ProjectionX();
      //TH1D *apd3_1D = apd3_prof->ProjectionX();
      //TH1D *apd4_1D = apd4_prof->ProjectionX();
      func->SetParameters(1,1,0,0,0);
      apd1_1D->Fit("fit","","",0,50);
      c_fscint = func->GetParameter(0);
      c_WLS = func->GetParameter(1);
      WLS_ratio = c_WLS/(c_WLS+c_fscint);
      cout << endl << WLS_ratio << endl;
      apd1_ratio_hist2D->SetBinContent(i+1,j+1,100*WLS_ratio);
    }
  }
  apd1_ratio_hist2D->Draw("colz");
  //TProfile *apd1_prof[nbins][nbins];
  //TH1D *apd1_1D[nbins][nbins];
  /*for(int i = 0;i < nbins;i++){
    for(int j = 0;j < nbins;j++){
      apd1_prof[i][j] = apd1[i][j]->ProfileX();
      apd1_1D[i][j] = apd1_prof[i][j]->ProjectionX();
    }
  }*/
  /*double WLS_ratio = 0.5; //    WLS/(WLS+fscint)
  double c_fscint = 1;
  double c_WLS = 1;
  func->SetParameters(1,1,0,0,0);
  //apd1_1D->Fit("fit","","",0,50);
  c_fscint = func->GetParameter(0);
  c_WLS = func->GetParameter(1);
  double adding_constant = func->GetParameter(2);
  double WLS_shift = func->GetParameter(3);
  double fscint_shift = func->GetParameter(4);
  WLS_ratio = c_WLS/(c_WLS+c_fscint);
  cout << endl << WLS_ratio << endl;
  cout << c_fscint << "      " << c_WLS << "      " << WLS_ratio << "      " << adding_constant << "      " << WLS_shift << "      " << fscint_shift << endl;*/
}