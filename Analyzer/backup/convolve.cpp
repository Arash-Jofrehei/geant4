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

TBranch        *b_nFibresOnX;
TBranch        *b_nFibresOnY;

TH1F* convoluted_pulse_fscint;
TH1F* convoluted_pulse_WLS;
Double_t fit_function(Double_t *v,Double_t *par)
{
  int n = int(v[0]*5);
  int WLS_shift = int(par[3]);
  int fscint_shift = int(par[4]);
  return par[0]*convoluted_pulse_fscint->GetBinContent(n+1+fscint_shift) + par[1]*convoluted_pulse_WLS->GetBinContent(n+1+WLS_shift) + par[2];
}

int main(int argc, char **argv)
{
TApplication theApp("tapp", &argc, argv);
TCanvas *canvas = new TCanvas("convolve","convolve");

TFile *final = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/merged_crystal4apd.root");
TFile *elecs = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/apd_plus_electronics.root");
TFile *WLS = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/4arash/plotterTiming_WLS_100GeV_qe.root");
TFile *fscint = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/plotterTiming_SingleFibre_50GeV.root");
TTree *ftree = (TTree*) final->Get("h4");
TTree *etree = (TTree*) elecs->Get("h4");
TH1F *WLS_hist0 = (TH1F*) WLS->Get("timeArrival");
TH1F *fscint_hist0 = (TH1F*) fscint->Get("timeArrival_scint");
TH1F *WLS_hist = new TH1F("WLS_hist","WLS",1024,-0.1,204.7);
TH1F *fscint_hist = new TH1F("fscint_hist","fiber scintillation",1024,-0.1,204.7);
for (int i = 0;i < 951;i++){
  WLS_hist->SetBinContent(i,(WLS_hist0->GetBinContent(i+73))/10000000.0);
}
for (int i = 951;i < 1024;i++){
  WLS_hist->SetBinContent(i,0);//0.00011);
}
for (int i = 0;i < 217;i++){
  fscint_hist->SetBinContent(i,(fscint_hist0->GetBinContent(i+33))/50000000.0);
}
for (int i = 217;i < 517;i++){
  fscint_hist->SetBinContent(i,0);//0.000025);
}
for (int i = 517;i < 1024;i++){
  fscint_hist->SetBinContent(i,0);//0.000007);
}
TH2F *e_time_val = new TH2F("e_time_val","e_time_val",1024,-0.1,204.7,19500,-16.25,178.75);
TH2F *apd1_time_val = new TH2F("apd1_time_val","apd1_time_val",1024,-0.1,204.7,124000,-44000,80000);
TH2F *apd2_time_val = new TH2F("apd2_time_val","apd2_time_val",1024,-0.1,204.7,124000,-44000,80000);
TH2F *apd3_time_val = new TH2F("apd3_time_val","apd3_time_val",1024,-0.1,204.7,124000,-44000,80000);
TH2F *apd4_time_val = new TH2F("apd4_time_val","apd4_time_val",1024,-0.1,204.7,124000,-44000,80000);
/*TF1 *f = new  TF1("f", "pow(exp(1.)*[1]*(x)/[0],[0])*exp(-[1]*(x))",0.,200.);
f->SetParameters(3.0, 0.8);
RooRealVar t("t","t", 0, 200);
t.setBins(1024,"cache") ;
t.setBinning(RooBinning(1024,0,200));
RooAbsPdf* rfa1 = RooFit::bindPdf(f,t);
TH1* histoResponse = rfa1->createHistogram("t",1024);
TGraph* waveGraphResp = new TGraph(histoResponse);
WaveformNew* waveRespNew = new WaveformNew(histoResponse->GetNbinsX(),waveGraphResp->GetX(),waveGraphResp->GetY());*/

etree->SetBranchAddress("WF_time2",WF_time2);
etree->SetBranchAddress("WF_val2",WF_val2);
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

  //------------------------------------------------------------------------------------------------------------
  
  Long64_t nentries = etree->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = etree->LoadTree(jentry);
    nb = etree->GetEntry(jentry);   nbytes += nb;
    for (int sample = 0;sample < 1009;sample++){
      e_time_val->Fill(WF_time2[sample]-72.2,3.25*WF_val2[sample]/100.0);
    }
    for (int sample = 648;sample < 1024;sample++){
      e_time_val->Fill(sample/5.0,0);//-1.2);
    }
  }

  nentries = ftree->GetEntriesFast();
  nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries/10;jentry++) {
    Long64_t ientry = ftree->LoadTree(jentry);
    nb = ftree->GetEntry(jentry);   nbytes += nb;
    if (TMath::Abs(x[0]-207.0)<0.5&&TMath::Abs(y[0]-294.5)<0.5)
    {
      for (int sample = 0;sample < 1009;sample++){
        apd1_time_val->Fill(WF_time6[sample]-28.6,20*WF_val6[sample]);
        apd2_time_val->Fill(WF_time7[sample]-28.6,20*WF_val7[sample]);
        apd3_time_val->Fill(WF_time8[sample]-28.6,20*WF_val8[sample]);
        apd4_time_val->Fill(WF_time9[sample]-28.6,20*WF_val9[sample]);
      }
      for (int sample = 866;sample < 1024;sample++){
        apd1_time_val->Fill(sample/5.0,0);//-1000);
        apd2_time_val->Fill(sample/5.0,0);//-1000);
        apd3_time_val->Fill(sample/5.0,0);//-1000);
        apd4_time_val->Fill(sample/5.0,0);//-1000);
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
TGraph* fscint_graph = new TGraph(fscint_hist);
TGraph* apd1_time_val_graph = new TGraph(apd1_time_val_1D);
TGraph* apd2_time_val_graph = new TGraph(apd2_time_val_1D);
TGraph* apd3_time_val_graph = new TGraph(apd3_time_val_1D);
TGraph* apd4_time_val_graph = new TGraph(apd4_time_val_1D);
WaveformNew *e_time_val_WFN = new WaveformNew(e_time_val_graph->GetN(),e_time_val_graph->GetX(),e_time_val_graph->GetY());
WaveformNew *WLS_WFN = new WaveformNew(WLS_graph->GetN(),WLS_graph->GetX(),WLS_graph->GetY());
WaveformNew *fscint_WFN = new WaveformNew(fscint_graph->GetN(),fscint_graph->GetX(),fscint_graph->GetY());
WaveformNew *apd1_time_val_WFN = new WaveformNew(apd1_time_val_graph->GetN(),apd1_time_val_graph->GetX(),apd1_time_val_graph->GetY());
WaveformNew *apd2_time_val_WFN = new WaveformNew(apd2_time_val_graph->GetN(),apd2_time_val_graph->GetX(),apd2_time_val_graph->GetY());
WaveformNew *apd3_time_val_WFN = new WaveformNew(apd3_time_val_graph->GetN(),apd3_time_val_graph->GetX(),apd3_time_val_graph->GetY());
WaveformNew *apd4_time_val_WFN = new WaveformNew(apd4_time_val_graph->GetN(),apd4_time_val_graph->GetX(),apd4_time_val_graph->GetY());

//----------------------------------------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------------------------------------------

//FFTConvolution* dec = new FFTConvolution();
//WaveformNew *deconvolved_apd1 = dec->fftConvolute(apd1_time_val_WFN,e_time_val_WFN,false);
//WaveformNew *convolved_apd1 = dec->fftConvolute(deconvolved_apd1,e_time_val_WFN);
FFTConvolution* con = new FFTConvolution();
WaveformNew *convoluted_WLS = con->fftConvolute(WLS_WFN,e_time_val_WFN);
WaveformNew *convoluted_fscint = con->fftConvolute(fscint_WFN,e_time_val_WFN);
//WaveformNew *convoluted = con->fftConvolute(WLS_WFN,waveRespNew);
//TH1F* hist_deconvolved_apd1 = deconvolved_apd1->get_histo("deconvolved_apd1");
//TH1F* hist_convolved_apd1 = convolved_apd1->get_histo("convolved_apd1");
convoluted_pulse_WLS = convoluted_WLS->get_histo("convoluted pulse WLS");
convoluted_pulse_fscint = convoluted_fscint->get_histo("convoluted pulse fiber scintillation");
convoluted_pulse_WLS->SetLineWidth(0);
convoluted_pulse_fscint->SetLineWidth(0);
TFile *convoluted_pulses = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/convoluted_pulses.root","recreate");
convoluted_pulses->cd();
convoluted_pulse_WLS->Write();
convoluted_pulse_fscint->Write();
//apd1_time_val_1D->SetLineWidth(0);
apd1_time_val_1D->Draw();
convoluted_pulse_WLS->Draw("same");
convoluted_pulse_fscint->Draw("SAME");
//histoResponse->Draw();
//e_time_val_1D->Draw();
//hist_convolved_apd1->Draw();
//hist_deconvolved_apd1->Draw();
//convoluted->Draw();
WLS_hist->Draw();
//fscint_hist->Draw("same");

for (int i = 0; i < 1024 ; i++){
  //cout << (apd1_time_val_graph->GetX())[i] << "    " << (apd1_time_val_graph->GetY())[i] << "    " << (e_time_val_graph->GetX())[i] << "    " << (e_time_val_graph->GetY())[i] << endl;
  //cout << i << "      " << convoluted_pulse_WLS->GetBinContent(i) << "        " << convoluted_pulse_fscint->GetBinContent(i) << "        " << apd1_time_val_1D->GetBinContent(i) << endl;
}

//----------------------------------------------------------------------------------------------------------------------------

double WLS_ratio = 0.5; //    WLS/(WLS+fscint)
double c_fscint = 1;
double c_WLS = 1;
TF1 *func = new TF1("fit",fit_function,0,132,2);
func->SetParameters(1,1,0,0,0);
func->SetParNames("fscint constant","WLS constant","adding constant","WLS shift","fscint shift");
func->SetParLimits(0,0,1000);
func->SetParLimits(1,0,1000);
apd1_time_val_1D->Fit("fit","","",0,50);
c_fscint = func->GetParameter(0);
c_WLS = func->GetParameter(1);
double adding_constant = func->GetParameter(2);
double WLS_shift = func->GetParameter(3);
double fscint_shift = func->GetParameter(4);
WLS_ratio = c_WLS/(c_WLS+c_fscint);
cout << c_fscint << "      " << c_WLS << "      " << WLS_ratio << "      " << adding_constant << "      " << WLS_shift << "      " << fscint_shift << endl;

//----------------------------------------------------------------------------------------------------------------------------
canvas->Update();
theApp.Run();

return 0;
}