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
int nbins = 20;
float x_start = 197;
float x_end = 217;
float y_start = 284.5;
float y_end = 304.5;
TH1F* apd_plus_electronics;
TH1F* convoluted_pulse_fscint;
TH1F* convoluted_pulse_WLS;
TH1D* deconv;

TH1D *Deconvolution(TH1F *p,TH1F *h,int size = 1024){
  double dummy;
  for (int sample = 0;sample < size;sample++){
    dummy = p->GetBinContent(sample+1);
    for (int point = 0;point < sample;point++){
      dummy -= deconv->GetBinContent(point+1)*h->GetBinContent(sample-point+1);
    }
    dummy /= h->GetBinContent(1);
    cout << sample << "    " << dummy << endl;
    deconv->SetBinContent(sample+1,dummy);
  }
  return deconv;
}

void deconvolution(){
  TCanvas *canvas = new TCanvas("deconvolution","deconvolution");
  deconv = new TH1D("d","d",1024,-0.1,204.7);
  TH1F *test = new TH1F("test","test",1024,-0.1,204.7);
  TFile *apd_profile_waveform_bins = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/apd_profile_waveform_bins_shift.root");
  TFile *convoluted_pulses = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/convoluted_pulses.root");
  apd_plus_electronics = (TH1F*) convoluted_pulses->Get("normalized APD+electronics");
  int offset = 0;
  for(int i = 0;i < 1024;i++){
    if ((apd_plus_electronics->GetBinContent(i+1)/apd_plus_electronics->GetBinContent(apd_plus_electronics->GetMaximumBin())) < 0.02){
      offset += 1;
    }
    if ((apd_plus_electronics->GetBinContent(i+1)/apd_plus_electronics->GetBinContent(apd_plus_electronics->GetMaximumBin())) > 0.02) break;
  }
  for(int i = 0;i < 1024-offset;i++){
    apd_plus_electronics->SetBinContent(i+1,apd_plus_electronics->GetBinContent(i+1+offset));
  }
  for(int i = 1024-offset;i < 1024;i++){
    apd_plus_electronics->SetBinContent(i+1,0);
  }
  TH2F *spike_hist2D_apd[4];
  spike_hist2D_apd[0] = new TH2F("spike_hist2D_apd1","nuclear counter effect contribution for APD1;X (mm);Y (mm)",nbins,x_start,x_end,nbins,y_start,y_end);
  spike_hist2D_apd[1] = new TH2F("spike_hist2D_apd2","nuclear counter effect contribution for APD2;X (mm);Y (mm)",nbins,x_start,x_end,nbins,y_start,y_end);
  spike_hist2D_apd[2] = new TH2F("spike_hist2D_apd3","nuclear counter effect contribution for APD3;X (mm);Y (mm)",nbins,x_start,x_end,nbins,y_start,y_end);
  spike_hist2D_apd[3] = new TH2F("spike_hist2D_apd4","nuclear counter effect contribution for APD4;X (mm);Y (mm)",nbins,x_start,x_end,nbins,y_start,y_end);
  TH2F *fscint_hist2D_apd[4];
  fscint_hist2D_apd[0] = new TH2F("fscint_hist2D_apd1","fiber scintillation contribution for APD1;X (mm);Y (mm)",nbins,x_start,x_end,nbins,y_start,y_end);
  fscint_hist2D_apd[1] = new TH2F("fscint_hist2D_apd2","fiber scintillation contribution for APD2;X (mm);Y (mm)",nbins,x_start,x_end,nbins,y_start,y_end);
  fscint_hist2D_apd[2] = new TH2F("fscint_hist2D_apd3","fiber scintillation contribution for APD3;X (mm);Y (mm)",nbins,x_start,x_end,nbins,y_start,y_end);
  fscint_hist2D_apd[3] = new TH2F("fscint_hist2D_apd4","fiber scintillation contribution for APD4;X (mm);Y (mm)",nbins,x_start,x_end,nbins,y_start,y_end);
  TH2F *WLS_hist2D_apd[4];
  WLS_hist2D_apd[0] = new TH2F("WLS_hist2D_apd1","WLS contribution for APD1;X (mm);Y (mm)",nbins,x_start,x_end,nbins,y_start,y_end);
  WLS_hist2D_apd[1] = new TH2F("WLS_hist2D_apd2","WLS contribution for APD2;X (mm);Y (mm)",nbins,x_start,x_end,nbins,y_start,y_end);
  WLS_hist2D_apd[2] = new TH2F("WLS_hist2D_apd3","WLS contribution for APD3;X (mm);Y (mm)",nbins,x_start,x_end,nbins,y_start,y_end);
  WLS_hist2D_apd[3] = new TH2F("WLS_hist2D_apd4","WLS contribution for APD4;X (mm);Y (mm)",nbins,x_start,x_end,nbins,y_start,y_end);
  TProfile *apd[4][nbins][nbins];
  string name;
  for (int a = 0;a < 4;a++){           //reading the waveform profiles for each bin
    for(int i = 0;i < nbins;i++){
      for(int j = 0;j < nbins;j++){
        name = "APD";
        name.append(to_string(a+1));
        name.append("_");
        name.append(to_string(i));
        name.append("_");
        name.append(to_string(j));
        const char *APD_name = name.c_str();
        apd[a][i][j] = (TProfile*) apd_profile_waveform_bins->Get(APD_name);
      }
    }
  }
  offset = 0;
  for(int i = 0;i < 1024;i++){
    if ((apd[0][1][1]->GetBinContent(i+1)/apd[0][1][1]->GetBinContent(apd[0][1][1]->GetMaximumBin())) < 0.0025){
      offset += 1;
    }
    if ((apd[0][1][1]->GetBinContent(i+1)/apd[0][1][1]->GetBinContent(apd[0][1][1]->GetMaximumBin())) > 0.0025) break;
  }
  for(int i = 0;i < 1024-offset;i++){
    test->SetBinContent(i+1,apd[0][1][1]->GetBinContent(i+1+offset));
  }
  for(int i = 1024-offset;i < 1024;i++){
    test->SetBinContent(i+1,0);
  }
  test->Smooth(1000);
  apd_plus_electronics->Smooth(1000);
  //TH1D *dec_hist = Deconvolution(apd_plus_electronics,apd_plus_electronics);
  TH1D *dec_hist = Deconvolution(test,apd_plus_electronics);
  dec_hist->Draw();
  TCanvas *canvas2 = new TCanvas("deconvolution2","deconvolution2");
  apd_plus_electronics->Draw();
  TCanvas *canvas3 = new TCanvas("deconvolution3","deconvolution3");
  test->Draw();
  cout << endl << endl << endl << endl << test->GetBinContent(1) << "        " << test->GetBinContent(2) << endl << apd_plus_electronics->GetBinContent(1) << "        " << apd_plus_electronics->GetBinContent(2) << endl;
}