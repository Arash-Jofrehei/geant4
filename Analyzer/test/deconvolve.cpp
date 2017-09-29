#include <iostream>
#include <string>
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
float WF_time6[1024];
float WF_val6[1024];
float WF_time7[1024];
float WF_val7[1024];
float WF_time8[1024];
float WF_val8[1024];
float WF_time9[1024];
float WF_val9[1024];
float time2[1024];
float val2[1024];
float time6[1024];
float val6[1024];
float time7[1024];
float val7[1024];
float time8[1024];
float val8[1024];
float time9[1024];
float val9[1024];
Float_t         x[2];
Float_t         y[2];
Int_t           nFibresOnX[2];
Int_t           nFibresOnY[2];

TBranch        *b_x;
TBranch        *b_y;
TBranch        *b_nFibresOnX;
TBranch        *b_nFibresOnY;

int main(int argc, char **argv)
{
TApplication theApp("tapp", &argc, argv);
TCanvas *canvas = new TCanvas("c","c");

TFile *final = new TFile("~/cern-summer-2016/H4Analysis/ntuples/analysis_4684.root");
TFile *elecs = new TFile("~/cern-summer-2016/H4Analysis/ntuples/apd_plus_electronics.root");
TFile *WLS = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/4arash/plotterTiming_WLS_100GeV_qe.root");
TTree *ftree = (TTree*) final->Get("h4");
TTree *etree = (TTree*) elecs->Get("reco");
TH1F *WLS_hist = (TH1F*) WLS->Get("timeArrival");
TH2F *e_time_val = new TH2F("e_time_val","e_time_val",1024,0,204.8,4000,-500,3500);
TH2F *apd1_time_val = new TH2F("apd1_time_val","apd1_time_val",1024,0,204.8,6200,-2200,4000);
TH2F *apd2_time_val = new TH2F("apd2_time_val","apd2_time_val",1024,0,204.8,6200,-2200,4000);
TH2F *apd3_time_val = new TH2F("apd3_time_val","apd3_time_val",1024,0,204.8,6200,-2200,4000);
TH2F *apd4_time_val = new TH2F("apd4_time_val","apd4_time_val",1024,0,204.8,6200,-2200,4000);

etree->SetBranchAddress("WF_time2",WF_time2);//,&b_WF_time2);
etree->SetBranchAddress("WF_val2",WF_val2);//,&b_WF_val2);
ftree->SetBranchAddress("WF_time6",WF_time6);//,&b_WF_time6);
ftree->SetBranchAddress("WF_val6",WF_val6);//,&b_WF_val6);
ftree->SetBranchAddress("WF_time7",WF_time7);//,&b_WF_time7);
ftree->SetBranchAddress("WF_val7",WF_val7);//,&b_WF_val7);
ftree->SetBranchAddress("WF_time8",WF_time8);//,&b_WF_time8);
ftree->SetBranchAddress("WF_val8",WF_val8);//,&b_WF_val8);
ftree->SetBranchAddress("WF_time9",WF_time9);//,&b_WF_time9);
ftree->SetBranchAddress("WF_val9",WF_val9);//,&b_WF_val9);
ftree->SetBranchAddress("x", x);//, &b_x);
ftree->SetBranchAddress("y", y);//, &b_y);
ftree->SetBranchAddress("nFibresOnX", nFibresOnX, &b_nFibresOnX);
ftree->SetBranchAddress("nFibresOnY", nFibresOnY, &b_nFibresOnY);

  //------------------------------------------------------------------------------------------------------------
  
  Long64_t nentries = etree->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = etree->LoadTree(jentry);
    nb = etree->GetEntry(jentry);   nbytes += nb;
    for (int sample = 0;sample < 1024;sample++){
      e_time_val->Fill(WF_time2[sample],WF_val2[sample]);
    }
  }
  
  nentries = ftree->GetEntriesFast();
  nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = ftree->LoadTree(jentry);
    nb = ftree->GetEntry(jentry);   nbytes += nb;
    if (TMath::Abs(x[0]-207)<2&&TMath::Abs(y[0]-294.5)<2)
    {
      for (int sample = 0;sample < 1024;sample++){
        apd1_time_val->Fill(WF_time6[sample],WF_val6[sample]);
        apd2_time_val->Fill(WF_time7[sample],WF_val7[sample]);
        apd3_time_val->Fill(WF_time8[sample],WF_val8[sample]);
        apd4_time_val->Fill(WF_time9[sample],WF_val9[sample]);
      }
    }
  }

TProfile *e_time_val_prof = e_time_val->ProfileX("e_time_val_prof");
TProfile *apd1_time_val_prof = apd1_time_val->ProfileX("apd1_time_val_prof");
TProfile *apd2_time_val_prof = apd2_time_val->ProfileX("apd2_time_val_prof");
TProfile *apd3_time_val_prof = apd3_time_val->ProfileX("apd3_time_val_prof");
TProfile *apd4_time_val_prof = apd4_time_val->ProfileX("apd4_time_val_prof");
TH1D *e_time_val_1D = e_time_val_prof->ProjectionX();
TH1D *apd1_time_val_1D = apd1_time_val_prof->ProjectionX();
TH1D *apd2_time_val_1D = apd2_time_val_prof->ProjectionX();
TH1D *apd3_time_val_1D = apd3_time_val_prof->ProjectionX();
TH1D *apd4_time_val_1D = apd4_time_val_prof->ProjectionX();
TGraph* e_time_val_graph = new TGraph(e_time_val_1D);
TGraph* WLS_graph = new TGraph(WLS_hist);
TGraph* apd1_time_val_graph = new TGraph(apd1_time_val_1D);
TGraph* apd2_time_val_graph = new TGraph(apd2_time_val_1D);
TGraph* apd3_time_val_graph = new TGraph(apd3_time_val_1D);
TGraph* apd4_time_val_graph = new TGraph(apd4_time_val_1D);
WaveformNew *e_time_val_WFN = new WaveformNew(e_time_val_graph->GetN(),e_time_val_graph->GetX(),e_time_val_graph->GetY());
WaveformNew *WLS_WFN = new WaveformNew(WLS_graph->GetN(),WLS_graph->GetX(),WLS_graph->GetY());
WaveformNew *apd1_time_val_WFN = new WaveformNew(apd1_time_val_graph->GetN(),apd1_time_val_graph->GetX(),apd1_time_val_graph->GetY());
WaveformNew *apd2_time_val_WFN = new WaveformNew(apd2_time_val_graph->GetN(),apd2_time_val_graph->GetX(),apd2_time_val_graph->GetY());
WaveformNew *apd3_time_val_WFN = new WaveformNew(apd3_time_val_graph->GetN(),apd3_time_val_graph->GetX(),apd3_time_val_graph->GetY());
WaveformNew *apd4_time_val_WFN = new WaveformNew(apd4_time_val_graph->GetN(),apd4_time_val_graph->GetX(),apd4_time_val_graph->GetY());

//----------------------------------------------------------------------------------------------------------------------------

FFTConvolution* dec = new FFTConvolution();
WaveformNew *deconvolved_apd1 = dec->fftConvolute(apd1_time_val_WFN,e_time_val_WFN,false);
WaveformNew *convolved_apd1 = dec->fftConvolute(deconvolved_apd1,e_time_val_WFN);
FFTConvolution* con = new FFTConvolution();
WaveformNew *convoluted = con->fftConvolute(WLS_WFN,e_time_val_WFN);
TH1F* hist_deconvolved_apd1 = deconvolved_apd1->get_histo("deconvolved_apd1");
TH1F* hist_convolved_apd1 = convolved_apd1->get_histo("convolved_apd1");
apd1_time_val_1D->Draw();
e_time_val_1D->Draw("SAME");
hist_convolved_apd1->Draw("SAME");
//hist_deconvolved_apd1->Draw();
//convoluted->Draw();
canvas->Update();
theApp.Run();

return 0;
}