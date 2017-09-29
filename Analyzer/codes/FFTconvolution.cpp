#include <iostream>
#include "TCanvas.h"
#include "TComplex.h"
using namespace std;
float WF_time2[1024];
float WF_val2[1024];
float Time[9];

TH1D *normalized(TH1D *histogram,int nsamplesNormalization = 450,int offset = 382){
  string strName = "normalized ";
  strName.append(string(histogram->GetName()));
  const char *name = strName.c_str();
  string strTitle = "normalized ";
  strTitle.append(string(histogram->GetTitle()));
  const char *title = strTitle.c_str();
  TH1D *normalized_histogram = new TH1D(name,title,1024,-0.1,204.7);
  nsamplesNormalization = histogram->GetMaximumBin() - 382;
  double integral = 0;
  for (int i = offset;i < nsamplesNormalization+offset;i++){
    integral += histogram->GetBinContent(i+1);
  }
  for (int i = 0;i < 1024;i++){
    normalized_histogram->SetBinContent(i+1,histogram->GetBinContent(i+1)/integral);
  }
  return normalized_histogram;
}

double tail(double *v,double *par){
  return par[0]*TMath::Exp(par[1]*v[0]);
}


double *convolution(double *f1,double *f2,int size = 1024){
  double conv[size];
  double dummy;
  for (int sample = 0;sample < size;sample++){
    dummy = 0;
    for (int point = 0;point <= sample;point++){
      dummy += f1[point]*f2[sample-point];
    }
    conv[sample] = dummy;
  }
  return conv;
}

void FFTconvolution(){
  TCanvas *canvas = new TCanvas("FFTconvolution","FFTconvolution");
  TFile *elecs = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/apd_plus_electronics_4444.root");
  TFile *WLS_file = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/4arash/plotterTiming_WLS_100GeV_qe.root");
  TFile *fscint_file = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/plotterTiming_SingleFibre_50GeV.root");
  TFile *spike_file = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/plotterTiming_Ideal2016_APD_100GeV.root");
  
  TFile *time_arrivals = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/time_arrivals.root");
  
  TTree *etree = (TTree*) elecs->Get("h4");
  TH1F *spike_hist0 = (TH1F*) spike_file->Get("timeArrival");
  spike_hist0->Scale(1.0/spike_hist0->Integral());
  TH1F *WLS_hist0 = (TH1F*) WLS_file->Get("timeArrival");
  TH1F *fscint_hist0 = (TH1F*) fscint_file->Get("timeArrival_scint");
  
  TH1F *fscint_hist = (TH1F*) time_arrivals->Get("fiber scintillation time arrival at center");
  TH1F *WLS_hist = (TH1F*) time_arrivals->Get("WLS time arrival at center");
  
  TF1 *tail = new TF1("tail","[0]*pow(e,[1]*x)",0,1000);
  tail->SetParameters(100000,-1);
  tail->SetParLimits(1,-1000,0);
  int nsamples = 1024;
  double spike[nsamples];
  TH1F *spike_hist = new TH1F("spike","spike",1024,-0.1,204.7);
  for (int i = 0;i < 400;i++){
    spike_hist->Fill(spike_hist0->GetBinLowEdge(i+1)-5.0,spike_hist0->GetBinContent(i+1));
  }
  double fscint[nsamples];
  /*TH1F *fscint_tail_hist = new TH1F("fscint tail hist","fscint tail hist",2048,-0.1,409.5);
  fscint_hist0->Fit("tail","","",13,48);
  fscint_tail_hist->Add(tail,1);
  TH1F *fscint_hist = new TH1F("fiber scintillation","fiber scintillation",1024,-0.1,204.7);
  for (int i = 0;i < 1024;i++){
    //if (i<70) fscint_hist->SetBinContent(i+1,fscint_hist0->GetBinContent(i+1+32));
    //else fscint_hist->SetBinContent(i+1,fscint_tail_hist->GetBinContent(i+1+32));
    fscint_hist->SetBinContent(fscint_hist0->FindBin(fscint_hist0->GetBinLowEdge(i+1)),fscint_hist0->GetBinContent(i+1));
  }*/
  double WLS[nsamples];
  /*TH1F *WLS_tail_hist = new TH1F("WLS tail hist","WLS tail hist",2048,-0.1,409.5);
  WLS_hist0->Fit("tail","","",70,190);
  WLS_tail_hist->Add(tail,1);
  TH1F *WLS_hist = new TH1F("WLS","WLS",1024,-0.1,204.7);
  for (int i = 0;i < 1024;i++){
    //if (i<750) WLS_hist->SetBinContent(i+1,WLS_hist0->GetBinContent(i+1+71));
    //else WLS_hist->SetBinContent(i+1,WLS_tail_hist->GetBinContent(i+1+71));
    WLS_hist->SetBinContent(WLS_hist0->FindBin(WLS_hist0->GetBinLowEdge(i+1)),WLS_hist0->GetBinContent(i+1));
  }*/
  for (int i = 0;i < 1024;i++){
    spike[i] = spike_hist->GetBinContent(i+1);
    fscint[i] = fscint_hist->GetBinContent(i+1);
    WLS[i] = WLS_hist->GetBinContent(i+1);
  }
  TProfile *apd_plus_electronics = new TProfile("APD+electronics response","APD+electronics response",1024,-0.1,204.7);
  TH1F *dummy = new TH1F("dummy","dummy",1024,-0.1,204.7);
  etree->SetBranchAddress("WF_time2",WF_time2);
  etree->SetBranchAddress("WF_val2",WF_val2);
  etree->SetBranchAddress("Time",Time);
  Long64_t nentries = etree->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = etree->LoadTree(jentry);
    nb = etree->GetEntry(jentry);   nbytes += nb;
    dummy->Reset();
    for (int sample = 0;sample < 1024;sample++){
      dummy->Fill(WF_time2[sample]-Time[2]+90,WF_val2[sample]);
    }
    for (int sample = 0;sample < 1024;sample++){
      apd_plus_electronics->Fill(WF_time2[sample]-Time[2]+90,WF_val2[sample]/dummy->GetBinContent(dummy->GetMaximumBin()));
    }
  }
  double apd_plus_electronics_array[nsamples];
  for (int i = 0;i < 1024;i++){
    apd_plus_electronics_array[i] = apd_plus_electronics->GetBinContent(i+1);
  }
  TVirtualFFT* fftr2c1 = TVirtualFFT::FFT(1, &nsamples, "R2C EX K");
  TVirtualFFT* fftr2c2 = TVirtualFFT::FFT(1, &nsamples, "R2C EX K");
  TVirtualFFT* fftr2c3 = TVirtualFFT::FFT(1, &nsamples, "R2C EX K");
  TVirtualFFT* fftc2r1 = TVirtualFFT::FFT(1, &nsamples, "C2R EX K");
  TVirtualFFT* fftc2r2 = TVirtualFFT::FFT(1, &nsamples, "C2R EX K");
  fftr2c1->SetPoints(fscint);
  fftr2c2->SetPoints(WLS);
  fftr2c3->SetPoints(apd_plus_electronics_array);
  fftr2c1->Transform();
  fftr2c2->Transform();
  fftr2c3->Transform();
  double r1,i1,r2,i2,r3,i3,re1,im1,re2,im2;
  for (int i = 0;i < nsamples;i++){
    fftr2c1->GetPointComplex(i,r1,i1);
    fftr2c2->GetPointComplex(i,r2,i2);
    fftr2c3->GetPointComplex(i,r3,i3);
    //if (i == 0) cout << r1 << "  " << r1 << "  " << r1 << "  " << r1 << "  " <<
    re1 = r1*r3 - i1*i3;
    im1 = r1*i3 + r3*i1;
    re2 = r2*r3 - i2*i3;
    im2 = r2*i3 + r3*i2;
    TComplex t1(re1,im1);
    TComplex t2(re2,im2);
    fftc2r1->SetPointComplex(i,t1);
    fftc2r2->SetPointComplex(i,t2);
  }
  fftc2r1->Transform();
  fftc2r2->Transform();
  double *conv_spike = convolution(spike,apd_plus_electronics_array);
  double convoluted_spike[1024];
  for (int i = 0;i < 1024;i++){
    convoluted_spike[i] = conv_spike[i];
  }
  double *conv_fscint = convolution(fscint,apd_plus_electronics_array);
  double convoluted_fscint[1024];
  for (int i = 0;i < 1024;i++){
    convoluted_fscint[i] = conv_fscint[i];
  }
  double *conv_WLS = convolution(WLS,apd_plus_electronics_array);
  double convoluted_WLS[1024];
  for (int i = 0;i < 1024;i++){
    convoluted_WLS[i] = conv_WLS[i];
  }
  TH1D *convoluted_pulse_spike = new TH1D("nuclear counter effect convoluted with APD+electronics","nuclear counter effect convoluted with APD+electronics",1024,-0.1,204.7);
  TH1D *convoluted_pulse_WLS = new TH1D("WLS+CeF3 convoluted with APD+electronics","WLS+CeF3 convoluted with APD+electronics",1024,-0.1,204.7);
  TH1D *convoluted_pulse_fscint = new TH1D("fiber scintillation convoluted with APD+electronics","fiber scintillation convoluted with APD+electronics",1024,-0.1,204.7);
  double dummy1;
  double dummy2;
  double dummy3;
  for (int i = 0;i < 1024;i++){
    //convoluted_pulse_fscint->SetBinContent(i+1,fftc2r1->GetPointReal(i));
    //convoluted_pulse_WLS->SetBinContent(i+1,fftc2r2->GetPointReal(i));
    convoluted_pulse_spike->SetBinContent(i+1,convoluted_spike[i]);
    convoluted_pulse_fscint->SetBinContent(i+1,convoluted_fscint[i]);
    convoluted_pulse_WLS->SetBinContent(i+1,convoluted_WLS[i]);
    /*dummy1 = 0;
    dummy2 = 0;
    dummy3 = 0;
    for (int point = 0;point <= i;point++){
      dummy1 += spike[point]*apd_plus_electronics_array[i-point];
      dummy2 += fscint[point]*apd_plus_electronics_array[i-point];
      dummy3 += WLS[point]*apd_plus_electronics_array[i-point];
    }
    convoluted_pulse_spike->SetBinContent(i+1,dummy1);
    convoluted_pulse_fscint->SetBinContent(i+1,dummy2);
    convoluted_pulse_WLS->SetBinContent(i+1,dummy3);*/
  }
  TH1D *apd_plus_electronics_1D = apd_plus_electronics->ProjectionX();
  apd_plus_electronics_1D->SetName("APD+electronics");
  apd_plus_electronics_1D->SetTitle("APD+electronics");
  /*TH1D *normalized_apd_plus_electronics = normalized(apd_plus_electronics_1D,382);
  TH1D *normalized_convoluted_pulse_spike = normalized(convoluted_pulse_spike,413);
  TH1D *normalized_convoluted_pulse_fscint = normalized(convoluted_pulse_fscint,424);
  TH1D *normalized_convoluted_pulse_WLS = normalized(convoluted_pulse_WLS,470);*/
  TH1D *normalized_apd_plus_electronics = apd_plus_electronics_1D;
  TH1D *normalized_convoluted_pulse_spike = convoluted_pulse_spike;
  TH1D *normalized_convoluted_pulse_fscint = convoluted_pulse_fscint;
  TH1D *normalized_convoluted_pulse_WLS = convoluted_pulse_WLS;
  //apd_plus_electronics->Draw();
  //fscint_hist->Draw();
  // fscint: bin# - 32 , WLS: bin# - 71
  //WLS_hist->Draw();
  //convoluted_pulse_fscint->Draw();
  //convoluted_pulse_WLS->Draw();
  for (int i = 0;i < 1024;i++){
    //cout << i << "    " << apd_plus_electronics->GetBinContent(i+1) << "    " << convoluted_pulse_spike->GetBinContent(i+1) << "    " << convoluted_pulse_fscint->GetBinContent(i+1) << "    " << convoluted_pulse_WLS->GetBinContent(i+1) << endl;
  }
  normalized_convoluted_pulse_spike->SetLineColor(2);
  normalized_convoluted_pulse_fscint->SetLineColor(3);
  normalized_convoluted_pulse_WLS->SetLineColor(4);
  
  TFile *convoluted_pulses = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/convoluted_pulses.root","recreate");
  convoluted_pulses->cd();
  normalized_apd_plus_electronics->Write();
  normalized_convoluted_pulse_spike->Write();
  normalized_convoluted_pulse_fscint->Write();
  normalized_convoluted_pulse_WLS->Write();
  
  normalized_apd_plus_electronics->Draw();
  normalized_convoluted_pulse_spike->Draw("same");
  normalized_convoluted_pulse_fscint->Draw("same");
  normalized_convoluted_pulse_WLS->Draw("same");
  //WLS_hist->Draw();
  TCanvas *canvas2 = new TCanvas("FFTconvolution2","FFTconvolution2");
  fscint_hist->Draw();
}