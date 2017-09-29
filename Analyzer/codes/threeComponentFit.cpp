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
TH1F* convoluted_pulse_spike;
TH1F* convoluted_pulse_fscint;
TH1F* convoluted_pulse_WLS;
double fit_function(double *v,double *par)
{
  int n = int(v[0]*5);
  int shift = int(par[3]);
  //int shift1 = int(par[3]);
  //int shift2 = int(par[4]);
  //int shift3 = int(par[5]);
  //the shift is the same for all three components. Fits don't work fine with different shifts and also fiber scintillation and WLS are the APD response convoluted with optical simulation WFs for which it was easy to find the rising point so I've already shifted the rising points to the same point.
  return par[0]*convoluted_pulse_spike->GetBinContent(n+1+shift)+par[1]*convoluted_pulse_fscint->GetBinContent(n+1+shift)+par[2]*convoluted_pulse_WLS->GetBinContent(n+1+shift)+par[4]; //I haven't used par[4] (vertical shift) here. It unstables the fits. For marginal crystals it might be useful since I noticed that due to the noise, the baseline is incorrect for many events.
  
  //return par[0]*apd_plus_electronics->GetBinContent(n+1+shift1)+par[1]*convoluted_pulse_fscint->GetBinContent(n+1+shift2)+par[2]*convoluted_pulse_WLS->GetBinContent(n+1+shift3)+par[6];
}

double FWHM(TProfile *h){
  int bin1 = h->FindFirstBinAbove(h->GetMaximum()/2);
  int bin2 = h->FindLastBinAbove(h->GetMaximum()/2);
  //int bin1 = h->GetMaximumBin();
  double fwhm = h->GetBinCenter(bin2) - h->GetBinCenter(bin1);
  return fwhm;
}
void threeComponentFit(){
  TCanvas *canvas = new TCanvas("threeComponentFit","threeComponentFit");
  TFile *apd_profile_waveform_bins = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/apd_profile_waveform_bins_shift.root");
  TFile *convoluted_pulses = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/convoluted_pulses.root");
  convoluted_pulse_spike = (TH1F*) convoluted_pulses->Get("normalized nuclear counter effect convoluted with APD+electronics");
  convoluted_pulse_fscint = (TH1F*) convoluted_pulses->Get("normalized fiber scintillation convoluted with APD+electronics");
  convoluted_pulse_WLS = (TH1F*) convoluted_pulses->Get("normalized WLS+CeF3 convoluted with APD+electronics");
  apd_plus_electronics = (TH1F*) convoluted_pulses->Get("normalized APD+electronics");
  TF1 *func = new TF1("fit",fit_function,-100,200,4);
  func->SetLineColor(6);
  func->SetParNames("spike on APD","fiber scintillation","WLS","spike shift","WLS & fiber scintillation shift","WLS shift","vertical offset");
  func->SetParLimits(0,0,100000000);
  func->SetParLimits(1,0,100000000);
  func->SetParLimits(2,0,100000000);
  func->SetParLimits(3,250,400);
  func->SetParLimits(4,250,420);
  func->SetParLimits(5,220,320);
  double spike_ratio = 0;
  double fscint_ratio = 0;
  double WLS_ratio = 0;
  double integrated_charge = 1;
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
  for (int i = 0;i < nbins;i++){
    for (int j = 0;j < nbins;j++){
      //cout << endl << i+1 << "    " << j+1 << endl;
      //for (int a = 0;a < 4;a++){
      for (int a = 0;a < 1;a++){
        apd[a][i][j]->SetLineWidth(0);
        //func->SetParameters(50000,1,1,30,0,0);
        func->SetParameters(20000,15000,20000,330);
        //func->SetParameters(30000,10000,30000,310,350,320);
        //fitting process depends a lot on the initial shift and more importantly the limits on the shift. Current ones work fine
        for (int k = 0;k < 30;k++){              //works with a lot less than 100 times. Sometimes fails to fit a few bins for less than ~20
          //apd[a][i][j]->Fit("fit","Q","",25,80);
        }
        for (int k = 0;k < 400;k++){              //works with a lot less than 100 times. Sometimes fails to fit a few bins for less than ~20
          //apd[a][i][j]->Fit("fit","Q","",10,60);
          //apd[a][i][j]->Fit("fit","Q","",25,80); ####
          //if (i == 9 && j == 9) apd[a][i][j]->Fit("fit","","",10,80);
          //else apd[a][i][j]->Fit("fit","Q","",10,80);
        }
        integrated_charge = func->GetParameter(0)+func->GetParameter(1)+func->GetParameter(2);
        spike_ratio = func->GetParameter(0)/integrated_charge;
        fscint_ratio = func->GetParameter(1)/integrated_charge;
        WLS_ratio = func->GetParameter(2)/integrated_charge;
        spike_hist2D_apd[a]->SetBinContent(i+1,j+1,func->GetParameter(0));
        fscint_hist2D_apd[a]->SetBinContent(i+1,j+1,func->GetParameter(1));
        WLS_hist2D_apd[a]->SetBinContent(i+1,j+1,func->GetParameter(2));
      }
    }
  }
  //func->SetParameters(30000,10000,30000,300);
  //func->SetParameters(30000,10000,30000,300,370,300);
  func->SetParameters(20000,15000,20000,330);
  //func->SetParameters(20000,15000,20000,310,350,320);
  for (int i = 0;i < 40; i++){
    //apd[0][9][9]->Fit("fit","Q","",10,100);
  }
  for (int i = 0;i < 400; i++){ // to check the fit for a specific bin
    //apd[0][9][9]->SetLineWidth(0);
    //apd[0][1][1]->Fit("fit","","",25,60);
    apd[0][1][1]->Fit("fit","Q","",25,80);
  }
  integrated_charge = func->GetParameter(0)+func->GetParameter(1)+func->GetParameter(2);
  spike_ratio = func->GetParameter(0)/integrated_charge;
  fscint_ratio = func->GetParameter(1)/integrated_charge;
  WLS_ratio = func->GetParameter(2)/integrated_charge;
  cout << "spike:      " << spike_ratio*100 << endl;
  cout << "fscint:      " << fscint_ratio*100 << endl;
  cout << "WLS:      " << WLS_ratio*100 << endl;
  
  //name = "APD1 center of crystal: ";
  name = "";
  name.append(to_string(spike_ratio*100));
  name.append("% spike(red) + ");
  name.append(to_string(fscint_ratio*100));
  name.append("% fiber scintillation(green) + ");
  name.append(to_string(WLS_ratio*100));
  name.append("% WLS(blue);25ns - 80ns fit with one time shift - impinging the APD");
  const char *hist_title = name.c_str();
  apd[0][1][1]->SetTitle(hist_title);
  apd[0][1][1]->Draw();
  
  TH1F *spike_hist = new TH1F("spike hist","nuclear counter effect contribution",1024,-0.1,204.7);
  TH1F *fscint_hist = new TH1F("fscint hist","fiber scintillation contribution",1024,-0.1,204.7);
  TH1F *WLS_hist = new TH1F("WLS hist","WLS contribution",1024,-0.1,204.7);
  for (int sample = 0;sample<1024;sample++){
    spike_hist->SetBinContent(sample+1,(convoluted_pulse_spike->GetBinContent(sample+1+int(func->GetParameter(3))))*func->GetParameter(0));
    fscint_hist->SetBinContent(sample+1,(convoluted_pulse_fscint->GetBinContent(sample+1+int(func->GetParameter(3))))*func->GetParameter(1));
    WLS_hist->SetBinContent(sample+1,(convoluted_pulse_WLS->GetBinContent(sample+1+int(func->GetParameter(3))))*func->GetParameter(2));
  }
  spike_hist->SetLineColor(2);
  fscint_hist->SetLineColor(3);
  WLS_hist->SetLineColor(4);
  spike_hist->Draw("same");
  fscint_hist->Draw("same");
  WLS_hist->Draw("same");
  
  //spike_hist2D_apd[0]->SetTitle("FWHM for APD1 pulse (ns)");
  //spike_hist2D_apd[0]->SetTitle("maximum amplitude - APD4");
  //spike_hist2D_apd[0]->Draw("colz");
  
  /*
  TCanvas *canvas2 = new TCanvas("threeComponentFit2","threeComponentFit2");
  fscint_hist2D_apd[0]->Draw("colz");
  func->SetParameters(20000,15000,20000,330);
  for (int i = 0;i < 400; i++){ // to check the fit for a specific bin
    //apd[0][16][16]->Fit("fit","Q","",10,100);
  }
  integrated_charge = func->GetParameter(0)+func->GetParameter(1)+func->GetParameter(2);
  spike_ratio = func->GetParameter(0)/integrated_charge;
  fscint_ratio = func->GetParameter(1)/integrated_charge;
  WLS_ratio = func->GetParameter(2)/integrated_charge;
  cout << endl << "spike:      " << spike_ratio*100 << endl;
  cout << "fscint:      " << fscint_ratio*100 << endl;
  cout << "WLS:      " << WLS_ratio*100 << endl;
  TH1F *spike_hist2 = new TH1F("spike hist2","nuclear counter effect contribution",1024,-0.1,204.7);
  TH1F *fscint_hist2 = new TH1F("fscint hist2","fiber scintillation contribution",1024,-0.1,204.7);
  TH1F *WLS_hist2 = new TH1F("WLS hist2","WLS contribution",1024,-0.1,204.7);
  for (int sample = 0;sample<1024;sample++){
    spike_hist2->SetBinContent(sample+1,(apd_plus_electronics->GetBinContent(sample+1+int(func->GetParameter(3))))*func->GetParameter(0));
    fscint_hist2->SetBinContent(sample+1,(convoluted_pulse_fscint->GetBinContent(sample+1+int(func->GetParameter(4))))*func->GetParameter(1));
    WLS_hist2->SetBinContent(sample+1,(convoluted_pulse_WLS->GetBinContent(sample+1+int(func->GetParameter(5))))*func->GetParameter(2));
  }
  spike_hist2->SetLineColor(2);
  fscint_hist2->SetLineColor(3);
  WLS_hist2->SetLineColor(4);
  spike_hist2->Draw("same");
  fscint_hist2->Draw("same");
  WLS_hist2->Draw("same");
  
  
  TCanvas *canvas3 = new TCanvas("threeComponentFit3","threeComponentFit3");
  WLS_hist2D_apd[0]->Draw("colz");
  func->SetParameters(20000,15000,20000,330);
  for (int i = 0;i < 400; i++){ // to check the fit for a specific bin
    //apd[0][1][1]->Fit("fit","Q","",10,100);
  }
  integrated_charge = func->GetParameter(0)+func->GetParameter(1)+func->GetParameter(2);
  spike_ratio = func->GetParameter(0)/integrated_charge;
  fscint_ratio = func->GetParameter(1)/integrated_charge;
  WLS_ratio = func->GetParameter(2)/integrated_charge;
  cout << endl << "spike:      " << spike_ratio*100 << endl;
  cout << "fscint:      " << fscint_ratio*100 << endl;
  cout << "WLS:      " << WLS_ratio*100 << endl;
  TH1F *spike_hist3 = new TH1F("spike hist3","nuclear counter effect contribution",1024,-0.1,204.7);
  TH1F *fscint_hist3 = new TH1F("fscint hist3","fiber scintillation contribution",1024,-0.1,204.7);
  TH1F *WLS_hist3 = new TH1F("WLS hist3","WLS contribution",1024,-0.1,204.7);
  for (int sample = 0;sample<1024;sample++){
    spike_hist3->SetBinContent(sample+1,(apd_plus_electronics->GetBinContent(sample+1+int(func->GetParameter(3))))*func->GetParameter(0));
    fscint_hist3->SetBinContent(sample+1,(convoluted_pulse_fscint->GetBinContent(sample+1+int(func->GetParameter(4))))*func->GetParameter(1));
    WLS_hist3->SetBinContent(sample+1,(convoluted_pulse_WLS->GetBinContent(sample+1+int(func->GetParameter(5))))*func->GetParameter(2));
  }
  spike_hist3->SetLineColor(2);
  fscint_hist3->SetLineColor(3);
  WLS_hist3->SetLineColor(4);
  spike_hist3->Draw("same");
  fscint_hist3->Draw("same");
  WLS_hist3->Draw("same");*/
  //apd_plus_electronics->Draw();
  //convoluted_pulse_fscint->Draw("same");
  //convoluted_pulse_WLS->Draw("same");
}