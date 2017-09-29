#include <iostream>
#include "string.h"
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TApplication.h"
#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooGenericPdf.h"
#include "RooFFTConvPdf.h"
#include "RooPlot.h"
#include "RooTFnBinding.h"
#include "RooTFnPdfBinding.h"
#include  "RooNumConvPdf.h"
#include  "RooBinning.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TVectorD.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/geant4/Analyzer/interface/Waveform.h"
#include "/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/geant4/Analyzer/interface/WaveformNew.h"
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include "TMath.h"
#include "/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/geant4/Analyzer/interface/FFTConvolution.h"

using namespace std;
float WF_time2[1024];
float WF_val2[1024];

TH1D *normalized(TH1F *histogram,int nsamplesNormalization = 450,int offset = 100){
  string strName = "normalized ";
  strName.append(string(histogram->GetName()));
  const char *name = strName.c_str();
  string strTitle = "normalized ";
  strTitle.append(string(histogram->GetTitle()));
  const char *title = strTitle.c_str();
  TH1D *normalized_histogram = new TH1D(name,title,1024,-0.1,204.7);
  double integral = 0;
  for (int i = 0;i < nsamplesNormalization;i++){
    integral += histogram->GetBinContent(i+1);
  }
  for (int i = 0;i < offset;i++){
    normalized_histogram->SetBinContent(i+1,0);
  }
  for (int i = offset;i < 1024;i++){
    normalized_histogram->SetBinContent(i+1,histogram->GetBinContent(i+1-offset)/integral);
  }
  return normalized_histogram;
}

int main(int argc, char **argv)
{
TApplication theApp("tapp", &argc, argv);
TCanvas *canvas = new TCanvas("convolve","convolve");

TFile *elecs = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/apd_plus_electronics.root");
TFile *WLS = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/4arash/plotterTiming_WLS_100GeV_qe.root");
TFile *fscint = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/plotterTiming_SingleFibre_50GeV.root");
TTree *etree = (TTree*) elecs->Get("h4");
TH1F *WLS_hist0 = (TH1F*) WLS->Get("timeArrival");
TH1F *fscint_hist0 = (TH1F*) fscint->Get("timeArrival_scint");
TH1F *WLS_hist = new TH1F("WLS_hist","WLS",1024,-0.1,204.7);
TH1F *fscint_hist = new TH1F("fscint_hist","fiber scintillation",1024,-0.1,204.7);
for (int i = 0;i < 951;i++){
  WLS_hist->SetBinContent(i,(WLS_hist0->GetBinContent(i+73)));
}
for (int i = 951;i < 1024;i++){
  WLS_hist->SetBinContent(i,0);
}
for (int i = 0;i < 217;i++){
  fscint_hist->SetBinContent(i,(fscint_hist0->GetBinContent(i+33)));
}
for (int i = 217;i < 517;i++){
  fscint_hist->SetBinContent(i,0);
}
for (int i = 517;i < 1024;i++){
  fscint_hist->SetBinContent(i,0);
}

TProfile *apd_plus_electronics = new TProfile("APD+electronics response","APD+electronics response",1024,-0.1,204.7);
etree->SetBranchAddress("WF_time2",WF_time2);
etree->SetBranchAddress("WF_val2",WF_val2);
/*
Long64_t nentries = etree->GetEntriesFast();
Long64_t nbytes = 0, nb = 0;
for (Long64_t jentry=0; jentry<nentries;jentry++) {
  Long64_t ientry = etree->LoadTree(jentry);
  nb = etree->GetEntry(jentry);   nbytes += nb;
  for (int sample = 0;sample < 1009;sample++){
    apd_plus_electronics->Fill(WF_time2[sample]-72.2,WF_val2[sample]);
  }
  for (int sample = 648;sample < 1024;sample++){
    apd_plus_electronics->Fill(sample/5.0,0);
  }
}*/

Long64_t nentries = etree->GetEntriesFast();
Long64_t nbytes = 0, nb = 0;
for (Long64_t jentry=0; jentry<nentries;jentry++) {
  Long64_t ientry = etree->LoadTree(jentry);
  nb = etree->GetEntry(jentry);   nbytes += nb;
  for (int sample = 0;sample < 1009;sample++){
    apd_plus_electronics->Fill(WF_time2[sample]-72.2,WF_val2[sample]);
  }
  for (int sample = 648;sample < 1024;sample++){
    apd_plus_electronics->Fill(sample/5.0,0);
  }
}

TH1D *apd_plus_electronics_1D = apd_plus_electronics->ProjectionX();
TGraph* apd_plus_electronics_graph = new TGraph(apd_plus_electronics_1D);
TGraph* WLS_graph = new TGraph(WLS_hist);
TGraph* fscint_graph = new TGraph(fscint_hist);
WaveformNew *apd_plus_electronics_WFN = new WaveformNew(apd_plus_electronics_graph->GetN(),apd_plus_electronics_graph->GetX(),apd_plus_electronics_graph->GetY());
WaveformNew *WLS_WFN = new WaveformNew(WLS_graph->GetN(),WLS_graph->GetX(),WLS_graph->GetY());
WaveformNew *fscint_WFN = new WaveformNew(fscint_graph->GetN(),fscint_graph->GetX(),fscint_graph->GetY());

//----------------------------------------------------------------------------------------------------------------------------

FFTConvolution* con = new FFTConvolution();
WaveformNew *convoluted_WLS = con->fftConvolute(WLS_WFN,apd_plus_electronics_WFN);
WaveformNew *convoluted_fscint = con->fftConvolute(fscint_WFN,apd_plus_electronics_WFN);
TH1F *convoluted_pulse_WLS = convoluted_WLS->get_histo("WLS+CeF3 convoluted with APD+electronics");
TH1F *convoluted_pulse_fscint = convoluted_fscint->get_histo("fiber scintillation convoluted with APD+electronics");
TH1D *normalized_convoluted_pulse_WLS = normalized(convoluted_pulse_WLS);
TH1D *normalized_convoluted_pulse_fscint = normalized(convoluted_pulse_fscint);
TH1D *normalized_apd_plus_electronics = normalized(apd_plus_electronics_WFN->get_histo("APD+electronics"));
//TFile *convoluted_pulses = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/convoluted_pulses.root","recreate");
//convoluted_pulses->cd();
//normalized_apd_plus_electronics->Write();
//normalized_convoluted_pulse_fscint->Write();
//normalized_convoluted_pulse_WLS->Write();
normalized_convoluted_pulse_fscint->Draw();
normalized_convoluted_pulse_WLS->Draw("same");
normalized_apd_plus_electronics->Draw("same");
//fscint_hist->Draw();
//WLS_hist->Draw("same");

//----------------------------------------------------------------------------------------------------------------------------
canvas->Update();
theApp.Run();

return 0;
}