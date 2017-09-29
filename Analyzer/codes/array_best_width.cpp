#include <iostream>
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "string.h"
using namespace std;
float max_WLS_spike_ratio = 2;
float max_fscint_spike_ratio = 2;
float WLS_ratio_jump = 0.002;
float fscint_ratio_jump = 0.002;
TH1F* convoluted_pulse_spike;
TH1F* convoluted_pulse_fscint;
TH1F* convoluted_pulse_WLS;
float Stat[6];

float *stat(TH1F *h){
  Stat[1] = h->FindFirstBinAbove(h->GetMaximum()/3);
  Stat[2] = h->GetMaximumBin();
  Stat[3] = h->FindLastBinAbove(h->GetMaximum()/6);
  Stat[0] = h->GetBinContent(Stat[2]);
  Stat[4] = 100.0*h->GetBinContent((Stat[2]+Stat[3])/2)/Stat[0];
  Stat[5] = h->FindFirstBinAbove(h->GetMaximum()/2);
  return Stat;
}

float *stat(TProfile *h){
  Stat[1] = h->FindFirstBinAbove(h->GetMaximum()/3);
  Stat[2] = h->GetMaximumBin();
  Stat[3] = h->FindLastBinAbove(h->GetMaximum()/6);
  Stat[0] = h->GetBinContent(Stat[2]);
  Stat[4] = 100.0*h->GetBinContent((Stat[2]+Stat[3])/2)/Stat[0];
  Stat[5] = h->FindFirstBinAbove(h->GetMaximum()/2);
  return Stat;
}

void array_best_width(){
  TCanvas *canvas = new TCanvas("array best width","array best width");
  TFile *convoluted_pulses = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/convoluted_pulses.root");
  convoluted_pulse_spike = (TH1F*) convoluted_pulses->Get("normalized nuclear counter effect convoluted with APD+electronics");
  convoluted_pulse_fscint = (TH1F*) convoluted_pulses->Get("normalized fiber scintillation convoluted with APD+electronics");
  convoluted_pulse_WLS = (TH1F*) convoluted_pulses->Get("normalized WLS+CeF3 convoluted with APD+electronics");
  apd_plus_electronics = (TH1F*) convoluted_pulses->Get("normalized APD+electronics");
  TH1F *sum = new TH1F("sum","sum",1024,-0.1,204.7);
  short int best_array[350][300][750][2];
  for (int jump_WLS = 0;jump_WLS < int(max_WLS_spike_ratio/WLS_ratio_jump);jump_WLS++){
    if (jump_WLS%100 == 0) cout << jump_WLS/10 << "    percent" << endl;
    for (int jump_fscint = 0;jump_fscint < int(max_fscint_spike_ratio/fscint_ratio_jump);jump_fscint++){
      sum->Reset();
      sum->Add(convoluted_pulse_spike);
      sum->Add(convoluted_pulse_fscint,jump_fscint*fscint_ratio_jump);
      sum->Add(convoluted_pulse_WLS,jump_WLS*WLS_ratio_jump);
      int width1 = int(stat(sum)[2]-stat(sum)[1]);
      int width2 = int(stat(sum)[2]-stat(sum)[5]);
      int width3 = int(stat(sum)[3]-stat(sum)[2]);
      if (best_array[width1][width2][width3][0] != -1) cout << "again!" << endl; 
      best_array[width1][width2][width3][0] = jump_fscint*fscint_ratio_jump;
      best_array[width1][width2][width3][1] = jump_WLS*WLS_ratio_jump;
    }
  }
}  