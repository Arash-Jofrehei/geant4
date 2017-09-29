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
int nbins = 19;
float x_start = -9.5;
float x_end = 9.5;
float y_start = -9.5;
float y_end = 9.5;
float WF_val6[1024];
float WF_time6[1024];
float Time[18];
Float_t         x[2];
Float_t         y[2];
Int_t           nFibresOnX[2];
Int_t           nFibresOnY[2];
TBranch        *b_nFibresOnX;
TBranch        *b_nFibresOnY;
TH1F* apd_plus_electronics;
TH1F* convoluted_pulse_fscint;
TH1F* convoluted_pulse_WLS;

TH1F *smoothed_histogram = new TH1F("smoothed histogram","smoothed histogram",1024,-0.1,204.7);
TH1F *smoothed(TH1F* histogram,int radius = 10,int iteration = 1,bool saveMAX = true){
  if (iteration == 0){
    for (int i = 0; i < 1024; i++) smoothed_histogram->SetBinContent(i+1,histogram->GetBinContent(i+1));
    return smoothed_histogram;
  }
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

TH1F* deconv;

TH1F *Deconvolution(TH1F *p,TH1F *h,int horizon = 150,int size = 1024){
  float main[size];
  float templ[size];
  for (int i = 0;i < size;i++){
    main[i] = p->GetBinContent(i+1);
    templ[i] = h->GetBinContent(i+1);
  }
  float dummy;
  for (int sample = 0;sample < size-4;sample++){
    if (sample > 500) horizon = 50;
    dummy = 0;
    dummy = main[sample]/templ[0];
    //if (main[sample]>0 && main[sample+1]>=0 && main[sample+2]>=0 && main[sample+3]>=0 && main[sample+4]>=0) dummy = main[sample]/templ[0];
    for (int i = sample;i < size-6;i++){
      //if (main[i]>=0 && main[i+1]>=0 && main[i+2]>=0 && main[i+3]>=0 && main[i+4]>=0 && templ[i-sample]>=0 && main[i]/templ[i-sample] < dummy && i < (sample+250)) dummy = main[i]/(1.0*templ[i-sample]);
      //****if (TMath::Abs(main[i]/templ[i-sample]) < TMath::Abs(dummy) && i < (sample+horizon)) dummy = main[i]/(1.0*templ[i-sample]);
      if (main[i]<=0 && i < (sample+horizon)) dummy = 0;
      if (TMath::Abs(main[i]/templ[i-sample]) < TMath::Abs(dummy) && TMath::Abs(main[i+6]/templ[i+6-sample]) < TMath::Abs(dummy) && i < (sample+horizon) && main[i]>=0) dummy = main[i]/(1.0*templ[i-sample]);
      //if (sample == 500) cout << i-sample+1 << "    " << main[i] << "    " << templ[i-sample] << "    " << dummy << endl;
    }
    for (int i = sample;i < size-6;i++){
      main[i] -= dummy*templ[i-sample];
      if (main[i] < 0 && main[i-1] > 0 && main[i+1] > 0) main[i] = 0.5*(main[i-1]+main[i+1]);
    }
    deconv->SetBinContent(sample+1,dummy);
  }
  return deconv;
}

/*
TH1F* deconv;

TH1F *Deconvolution(TH1F *p,TH1F *h,int size = 1024){
  float dummy;
  for (int sample = 0;sample < size;sample++){
    dummy = p->GetBinContent(sample+1);
    for (int point = 0;point < sample;point++){
      dummy -= deconv->GetBinContent(point+1)*h->GetBinContent(sample-point+1);
    }
    dummy /= h->GetBinContent(1);
    cout << sample << "    " << dummy << endl;
    if (dummy < 0) dummy = 0;
    deconv->SetBinContent(sample+1,dummy);
  }
  return deconv;
}
*/
TH1F *conv;
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
}

void deconvolution(){
  TCanvas *canvas = new TCanvas("deconvolution","deconvolution");
  TFile *time_arrivals = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/time_arrivals.root");
  TH1F *fscint_hist_center = (TH1F*) time_arrivals->Get("fiber scintillation time arrival at center");
  TH1F *WLS_hist_center = (TH1F*) time_arrivals->Get("WLS time arrival at center");
  TH1F *fscint_hist_fiber = (TH1F*) time_arrivals->Get("fiber scintillation time arrival near fiber");
  TH1F *WLS_hist_fiber = (TH1F*) time_arrivals->Get("WLS time arrival near fiber");
  TFile *final = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/merged_crystal4apd.root");
  TTree *ftree = (TTree*) final->Get("h4");
  ftree->SetBranchAddress("WF_val6",WF_val6);
  ftree->SetBranchAddress("WF_time6",WF_time6);
  ftree->SetBranchAddress("Time",Time);
  ftree->SetBranchAddress("x", x);
  ftree->SetBranchAddress("y", y);
  ftree->SetBranchAddress("nFibresOnX", nFibresOnX, &b_nFibresOnX);
  ftree->SetBranchAddress("nFibresOnY", nFibresOnY, &b_nFibresOnY);
  deconv = new TH1F("d","d",1024,-0.1,204.7);
  conv = new TH1F("c","c",1024,-0.1,204.7);
  TH1F *test = new TH1F("APD1 fiber WF(black) - conv. of deconv. (red) valid for positive values","APD1 fiber WF(black) - conv. of deconv. (red) valid for positive values;time(ns)",1024,-0.1,204.7);
  TH1F *test2 = new TH1F("APD1 center of crystal(black) - conv. of deconv. (red) valid for positive values","APD1 center of crystal(black) - conv. of deconv. (red) valid for positive values;time(ns)",1024,-0.1,204.7);
  TFile *apd_profile_waveform_bins = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/apd_profile_waveform_19bins_shift.root");
  TFile *convoluted_pulses = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/convoluted_pulses.root");
  apd_plus_electronics = (TH1F*) convoluted_pulses->Get("APD+electronics");
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
  for(int i = 0;i < 1024;i++){
    test->SetBinContent(i+1,apd[0][2][2]->GetBinContent(i+1));
  }
  TH1F *smoothed_test = smoothed(smoothed(test,10,1,false),100,1,true);
  for (int i = 0;i < 1024;i++){
    test->SetBinContent(i+1,smoothed_test->GetBinContent(i+1));
  }
  for(int i = 0;i < 1024;i++){
    test2->SetBinContent(i+1,apd[0][9][9]->GetBinContent(i+1));
  }
  TH1F *smoothed_test2 = smoothed(smoothed(test2,10,1,false),100,1,true);
  for (int i = 0;i < 1024;i++){
    test2->SetBinContent(i+1,smoothed_test2->GetBinContent(i+1));
  }
  TH1F *smoothed_apd_plus_electronics = smoothed(apd_plus_electronics,10,1,true);
  smoothed_apd_plus_electronics->Draw();
  for (int i = 0;i < 1024;i++){
    apd_plus_electronics->SetBinContent(i+1,smoothed_apd_plus_electronics->GetBinContent(i+1));
  }
  int offset = 0;
  for(int i = 0;i < 1024;i++){
    if ((apd_plus_electronics->GetBinContent(i+1)/apd_plus_electronics->GetBinContent(apd_plus_electronics->GetMaximumBin())) < 0.001){
      offset += 1;
    }
    if ((apd_plus_electronics->GetBinContent(i+1)/apd_plus_electronics->GetBinContent(apd_plus_electronics->GetMaximumBin())) > 0.001) break;
  }
  for(int i = 0;i < 1024-offset;i++){
    /*if (apd_plus_electronics->GetBinContent(i+1+offset) > 0) */apd_plus_electronics->SetBinContent(i+1,apd_plus_electronics->GetBinContent(i+1+offset));
    //else apd_plus_electronics->SetBinContent(i+1,0);
  }
  for(int i = 1024-offset;i < 1024;i++){
    apd_plus_electronics->SetBinContent(i+1,0);
  }
  /*offset = 0;
  for(int i = 0;i < 1024;i++){
    if ((apd[0][2][2]->GetBinContent(i+1)/apd[0][2][2]->GetBinContent(apd[0][2][2]->GetMaximumBin())) < 0.001){
      offset += 1;
    }
    if ((apd[0][2][2]->GetBinContent(i+1)/apd[0][2][2]->GetBinContent(apd[0][2][2]->GetMaximumBin())) > 0.001) break;
  }
  for(int i = 0;i < 1024-offset;i++){
    test->SetBinContent(i+1,apd[0][2][2]->GetBinContent(i+1+offset));
  }
  for(int i = 1024-offset;i < 1024;i++){
    test->SetBinContent(i+1,0);
  }*/
  //test->Smooth(1000);
  //apd_plus_electronics->Smooth(1000);
  //TH1D *dec_hist = Deconvolution(apd_plus_electronics,apd_plus_electronics);
  
  TH1F *waveform = new TH1F("waveform","waveform",1024,-0.1,204.7);
  TH1F *smoothed_waveform = new TH1F("smoothed waveform","smoothed waveform",1024,-0.1,204.7);
  TH1F *total_deconv = new TH1F("overall deconvolution of signals from APD1","overall deconvolution of signals from APD1 near center;time (ns)",1024,-0.1,204.7);
  TH1F *average_waveform = new TH1F("average waveform","average waveform",1024,-0.1,204.7);
  Long64_t nentries = ftree->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  int counter = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++){
    if (jentry%1000 == 0) cout << jentry << " / " << nentries << endl;
    Long64_t ientry = ftree->LoadTree(jentry);
    nb = ftree->GetEntry(jentry);   nbytes += nb;
    //if (TMath::Abs(x[1]-207.0+7)>1||TMath::Abs(y[1]-294.5+7)>1) continue;
    if (TMath::Abs(x[1]-207.0)>1||TMath::Abs(y[1]-294.5)>1) continue;
    waveform->Reset();
    for (int sample = 0;sample < 1024;sample++){
      waveform->Fill(30+WF_time6[sample]-Time[6],WF_val6[sample]);
    }
    //if (waveform->GetBinContent(waveform->GetMaximumBin()) < 1500 || waveform->GetBinContent(waveform->GetMaximumBin()) > 2500) continue;
    counter++;
    //if (counter != 3) continue;
    smoothed_waveform->Reset();
    //smoothed_waveform->Add(smoothed(smoothed(waveform,10,1,false),300,1,true));
    smoothed_waveform->Add(smoothed(smoothed(waveform,10,1,false),150,0,true));
    /*offset = -50;
    for(int i = 0;i < 1024;i++){
      if ((smoothed_waveform->GetBinContent(i+1)/smoothed_waveform->GetBinContent(smoothed_waveform->GetMaximumBin())) < 0.005){
        offset += 1;
      }
      if ((smoothed_waveform->GetBinContent(i+1)/smoothed_waveform->GetBinContent(smoothed_waveform->GetMaximumBin())) > 0.005) break;
    }
    for(int i = 0;i < 1024-offset;i++){
      smoothed_waveform->SetBinContent(i+1,smoothed_waveform->GetBinContent(i+1+offset));
    }
    for(int i = 1024-offset;i < 1024;i++){
      smoothed_waveform->SetBinContent(i+1,0);
    }*/
    average_waveform->Add(smoothed_waveform);
    total_deconv->Add(Deconvolution(smoothed_waveform,apd_plus_electronics,150));
  }
  //counter = 1;
  //total_deconv = Deconvolution(average_waveform,apd_plus_electronics,150);
  average_waveform->Scale(1.0/counter);
  total_deconv->Scale(1.0/counter);
  
  
  TH1F *dec_hist = new TH1F("dec_hist","deconvolved APD1 WF : fiber;time(ns)",1024,-0.1,204.7);
  dec_hist->Add(Deconvolution(test,apd_plus_electronics));
  TH1F *dec_hist2 = new TH1F("dec_hist2","deconvolved APD1 WF : center;time(ns)",1024,-0.1,204.7);
  dec_hist2->Add(Deconvolution(test2,apd_plus_electronics));
  TH1F *integral_dec_hist2 = new TH1F("integral of deconvolved APD1 WF : center;time(ns)","",1024,-0.1,204.7);
  for (int i = 0;i < 1024;i++){
    integral_dec_hist2->SetBinContent(i+1,dec_hist2->Integral(1,i+1));
  }
  integral_dec_hist2->Smooth(0);
  for (int i = 1024;i > 1;i--){
    //integral_dec_hist2->SetBinContent(i,integral_dec_hist2->GetBinContent(i)-integral_dec_hist2->GetBinContent(i-1));
  }
  conv->SetLineColor(2);
  convolution(dec_hist2,apd_plus_electronics);
  //test->Add(conv,-1);
  test2->Draw();
  //conv->Smooth(10);
  conv->Draw("same");
  total_deconv->Draw();
  TCanvas *canvas2 = new TCanvas("deconvolution2","deconvolution2");
  apd_plus_electronics->Draw();
  integral_dec_hist2->Draw();
  TH1F *conv_center = convolution(total_deconv,apd_plus_electronics);
  conv_center->SetTitle("APD1 center of crystal(black) - overall conv. of deconv. per event(red) valid for positive values;time(ns)");
  conv_center->Draw();
  //average_waveform->Scale(conv_center->GetBinContent(conv_center->GetMaximumBin())/test2->GetBinContent(test2->GetMaximumBin()));
  cout << "number of events: " << counter << endl;
  average_waveform->Draw("same");
  //smoothed(smoothed(waveform,1,1,false),1,1,true)->Draw("same");
  //smoothed(smoothed(waveform,10,1,false),150,1,true)->Draw("same");
  TCanvas *canvas3 = new TCanvas("deconvolution3","deconvolution3");
  dec_hist2->Draw();
  fscint_hist_center->Draw();
  WLS_hist_center->Draw("same");
  smoothed_waveform->Add(smoothed_waveform,-1);
  smoothed_waveform->Add(average_waveform,1);
  smoothed_waveform->Add(conv_center,-1);
  smoothed_waveform->SetTitle("APD1 waveform at the center minus conv. of deconv. of it;time (ns)");
  smoothed_waveform->Draw();
}