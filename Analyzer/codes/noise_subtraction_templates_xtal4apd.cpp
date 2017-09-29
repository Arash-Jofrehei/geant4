#include <iostream>
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
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
float WLS[18];
float spike[18];
float positionX[2];
float positionY[2];
float nFibersX[2];
float nFibersY[2];
Int_t           nFibresOnX[2];
Int_t           nFibresOnY[2];
TBranch        *b_nFibresOnX;
TBranch        *b_nFibresOnY;
float channel_weights[18] = {1,1,1,0,0,1,1,1,1,1,1,0,0,1,1,1,0,0};//{1.28673,1.37646,0.862157,0,0,0.38369,0.9088*0.951252,0.9088*0.932984,0.9088*1.16544,0.9088*0.981471,1.14679,0,0,0.234947,1.55411,0.795243,0,0};//{0.3*5.80376,0.7*3.52714,0.3*4.26522,0,0,0.7*1,0.951252,0.932984,1.16544,0.981471,0.7*2.79305,0,0,0.3*1,0.7*3.80542,0.3*3.21488,0,0};//{0,1.72464,0,0,0,0.470221,0.951252,0.932984,1.16544,0.981471,1.35214,0,0,0,1.80526,0,0,0};//{2.29894,1.72464,1.60079,0,0,0.470221,0.951252,0.932984,1.16544,0.981471,1.35214,0,0,0.451124,1.80526,1.3819,0,0};//{1,1,1,0,0,0.3,0.4,0.4,0.4,0.4,1,0,0,0.4,1,1,0,0};//{4.22939,3.03897,1.9653,0,0,1,0.810848,0.552276,0.541423,1,2.11745,0,0,1,3.73868,3.58402,0,0};//{0.521416,0.791515,0.524998,0.504333,0.506726,0.517119,0.50302,1,1.00179,0.512666,0.622622,0.899577};//{0.932688,1.56220,0.850238,0.492765,0.482391,0.606879,0.508326,1,1.37240,0.963999,1.44540,0.862974};
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
void noise_subtraction_templates_xtal4apd(){
  TCanvas *canvas = new TCanvas("noise_subtraction_templates","noise_subtraction_templates");
  TFile *convoluted_pulses = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/convoluted_pulses.root");
  convoluted_pulse_fscint = (TH1F*) convoluted_pulses->Get("normalized fiber scintillation convoluted with APD+electronics");
  convoluted_pulse_WLS = (TH1F*) convoluted_pulses->Get("normalized WLS+CeF3 convoluted with APD+electronics");
  apd_plus_electronics = (TH1F*) convoluted_pulses->Get("normalized APD+electronics");
  TFile *mean_noise_4APDs = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/mean_noise_4APDs");
  TH1 *mean_noise_4APD[4];
  //mean_noise_4APD[0] = new TH1F("mean noise APD1","mean noise APD1",1024,-0.1,204.7);
  //mean_noise_4APD[1] = new TH1F("mean noise APD2","mean noise APD2",1024,-0.1,204.7);
  //mean_noise_4APD[2] = new TH1F("mean noise APD3","mean noise APD3",1024,-0.1,204.7);
  //mean_noise_4APD[3] = new TH1F("mean noise APD4","mean noise APD4",1024,-0.1,204.7);
  mean_noise_4APD[0] = (TH1*) mean_noise_4APDs->Get("mean noise APD1");
  mean_noise_4APD[1] = (TH1*) mean_noise_4APDs->Get("mean noise APD2");
  mean_noise_4APD[2] = (TH1*) mean_noise_4APDs->Get("mean noise APD3");
  mean_noise_4APD[3] = (TH1*) mean_noise_4APDs->Get("mean noise APD4");
  //TFile *final = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/analysis_4787.root");
  TFile *final = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/analysis_4684.root","update");
  TTree *ftree = (TTree*) final->Get("h4");
  TTree *templates = new TTree("templates","templates");
  templates->Branch("WLS",&WLS,"WLS[18]/F");
  templates->Branch("spike",&spike,"spike[18]/F");
  templates->Branch("positionX",positionX,"positionX[2]/F");
  templates->Branch("positionY",positionY,"positionY[2]/F");
  templates->Branch("nFibersX",nFibersX,"nFibersX[2]/F");
  templates->Branch("nFibersY",nFibersY,"nFibersY[2]/F");
  TH1F *waveform[13];
  waveform[0] = new TH1F("waveform0","waveform crystal 1",1024,-0.1,204.7);
  waveform[1] = new TH1F("waveform1","waveform crystal 2",1024,-0.1,204.7);
  waveform[2] = new TH1F("waveform2","waveform crystal 3",1024,-0.1,204.7);
  waveform[3] = new TH1F("waveform5","waveform crystal 6",1024,-0.1,204.7);
  waveform[4] = new TH1F("waveform6","waveform crystal4APD 1",1024,-0.1,204.7);
  waveform[5] = new TH1F("waveform7","waveform crystal4APD 2",1024,-0.1,204.7);
  waveform[6] = new TH1F("waveform8","waveform crystal4APD 3",1024,-0.1,204.7);
  waveform[7] = new TH1F("waveform9","waveform crystal4APD 4",1024,-0.1,204.7);
  waveform[8] = new TH1F("waveform10","waveform crystal 11",1024,-0.1,204.7);
  waveform[9] = new TH1F("waveform13","waveform crystal 14",1024,-0.1,204.7);
  waveform[10] = new TH1F("waveform14","waveform crystal 15",1024,-0.1,204.7);
  waveform[11] = new TH1F("waveform15","waveform crystal 16",1024,-0.1,204.7);
  waveform[12] = new TH1F("waveform17","waveform crystal 18",1024,-0.1,204.7);
  TH1F *difference_hist = new TH1F("difference","difference",1024,-0.1,204.7);
  TH1F *difference_hist2 = new TH1F("difference2","difference2",1024,-0.1,204.7);
  TH1F *difference_hist3 = new TH1F("difference3","difference3",1024,-0.1,204.7);
  float E = 1;
  TH1F *integrated_charge_APD[4];
  //integrated_charge_APD[0] = new TH1F("integrated charge APD1","integrated charge APD1",100*E,25000*E,250000*E);
  //integrated_charge_APD[1] = new TH1F("integrated charge APD2","integrated charge APD2",100*E,25000*E,250000*E);
  //integrated_charge_APD[2] = new TH1F("integrated charge APD3","integrated charge APD3",100*E,25000*E,250000*E);
  //integrated_charge_APD[3] = new TH1F("integrated charge APD4","integrated charge APD4",100*E,25000*E,250000*E);
  integrated_charge_APD[0] = new TH1F("integrated charge APD1","integrated charge APD1",100*E,250*E,20000*E);
  integrated_charge_APD[1] = new TH1F("integrated charge APD2","integrated charge APD2",100*E,250*E,20000*E);
  integrated_charge_APD[2] = new TH1F("integrated charge APD3","integrated charge APD3",200*E,250*E,40000*E);
  integrated_charge_APD[3] = new TH1F("integrated charge APD4","integrated charge APD4",100*E,250*E,20000*E);
  double integrated_charge = 0;
  double calibrated_integrated_charge = 0;
  TH1D *integrated_charge_crystal11 = new TH1D("integrated charge crystal 4APD","crystal 4APD integrated charge (single crystal)",180*E,100000*E,1800000*E);
  TH1D *calibrated_integrated_charge_crystal11 = new TH1D("calibrated integrated charge crystal 4APD","crystal 4APD integrated charge (3*3 matrix)",180*E,100000,1800000*E);
  TF1 *crys = new TF1("crys","crystalball",0,5000000);
  TF1 *g = new TF1("g","gaus",0,5000000);
  TF1 *func = new TF1("fit",fit_function,-100,200,4);
  func->SetParNames("spike on APD","fiber scintillation","WLS","shift","vertical offset");
  func->SetParLimits(0,0,100000000);
  func->SetParLimits(1,0,100000000);
  func->SetParLimits(2,0,100000000);
  func->SetParLimits(3,20,180);
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
  ftree->SetBranchAddress("nFibresOnY", nFibresOnY, &b_nFibresOnY);
  float charge[12];
  int counter = 0;
  Long64_t nentries = ftree->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++){
    if (jentry%500 == 0) cout << jentry << endl;
    Long64_t ientry = ftree->LoadTree(jentry);
    nb = ftree->GetEntry(jentry);   nbytes += nb;
    for (int i = 0;i < 2;i++){
      positionX[i] = x[i];
      positionY[i] = y[i];
      nFibersX[i] = nFibresOnX[i];
      nFibersY[i] = nFibresOnY[i];
    }
    for (int i = 0;i < 18;i++) WLS[i] = 0;
    for (int i = 0;i < 18;i++) spike[i] = 0;
    if (TMath::Abs(x[1]-207)<2&&TMath::Abs(y[1]-294.5)<2&&nFibresOnX[1]<3&&nFibresOnY[1]<3){
      //cout << counter << endl;
      //cout << endl << endl << x[1]-207 << "    " << y[1]-294.5 << endl;
      for (int i = 0;i < 13;i++) waveform[i]->Reset();
      //difference_hist->Reset();
      for (int sample = 0;sample < 1024;sample++){
        //waveform[0]->Fill(30+WF_time10[sample]-Time[1],channel_weights[0]*WF_val0[sample]-WF_val17[sample]);
        waveform[0]->Fill(WF_time10[sample]-25,channel_weights[0]*WF_val0[sample]-WF_val17[sample]);
        waveform[1]->Fill(WF_time10[sample]-25,channel_weights[1]*WF_val1[sample]-WF_val17[sample]);
        waveform[2]->Fill(WF_time10[sample]-25,channel_weights[2]*WF_val2[sample]-WF_val17[sample]);
        //waveform[1]->Fill(30+WF_time10[sample]-Time[1],channel_weights[1]*WF_val1[sample]-WF_val17[sample]);
        //waveform[2]->Fill(30+WF_time10[sample]-Time[2],channel_weights[2]*WF_val2[sample]-WF_val17[sample]);
        waveform[3]->Fill(30+WF_time10[sample]-Time[5],channel_weights[5]*WF_val5[sample]);
        waveform[4]->Fill(30+WF_time10[sample]-Time[6],channel_weights[6]*WF_val6[sample]);
        waveform[5]->Fill(30+WF_time10[sample]-Time[7],channel_weights[7]*WF_val7[sample]);
        waveform[6]->Fill(30+WF_time10[sample]-Time[8],channel_weights[8]*WF_val8[sample]);
        waveform[7]->Fill(30+WF_time10[sample]-Time[9],channel_weights[9]*WF_val9[sample]);
        waveform[8]->Fill(30+WF_time10[sample]-Time[10],channel_weights[10]*WF_val10[sample]-WF_val17[sample]);
        waveform[9]->Fill(WF_time10[sample]-25,channel_weights[13]*WF_val13[sample]);
        //waveform[9]->Fill(30+WF_time10[sample]-Time[13],channel_weights[13]*WF_val13[sample]-WF_val17[sample]);
        waveform[10]->Fill(30+WF_time10[sample]-Time[14],channel_weights[14]*WF_val14[sample]-WF_val17[sample]);
        //waveform[11]->Fill(30+WF_time10[sample]-Time[14],channel_weights[15]*WF_val15[sample]-WF_val17[sample]);
        waveform[11]->Fill(WF_time10[sample]-25,channel_weights[15]*WF_val15[sample]-WF_val17[sample]);
        waveform[12]->Fill(30+WF_time10[sample]-Time[7],WF_val17[sample]);
      }
      //for (int i = 0;i < 12;i++) difference_hist->Add(waveform[i],1);
      integrated_charge = 0;
      calibrated_integrated_charge = 0;
      /*func->SetParameters(300000,0.4,70,0,0);
      for (int j = 0;j < 10;j++){
        //waveform[12]->Fit("fit","","",30,80);
      }
      difference_hist->Add(waveform[11],1);
      //difference_hist->Add(func,-1);
      difference_hist->Add(waveform[12],-1);
      waveform[11]->Draw();
      difference_hist->Draw("same");
      waveform[12]->Draw("same");*/
      func->FixParameter(1,0);
      for (int i = 0;i < 12;i++){
        waveform[i]->SetLineWidth(0);
        //waveform[i]->Smooth(100);
        if (i > 3 && i < 8){
          waveform[i]->Smooth(1);
          func->SetParameters(channel_weights[7]*500000,0,channel_weights[7]*10000,45,0,0);
          for (int j = 0;j < 30;j++){
            waveform[i]->Fit("fit","Q","",25,90);
          }
          integrated_charge += func->GetParameter(0) + func->GetParameter(2);
          calibrated_integrated_charge += func->GetParameter(0) + func->GetParameter(2);
          WLS[i+2] = func->GetParameter(2);
          spike[i+2] = func->GetParameter(0);
          //cout << i << "            " << func->GetParameter(0) + func->GetParameter(2) - 314000 << endl;
        }
        if (i == 1){
          waveform[i]->Smooth(1);
          func->SetParameters(channel_weights[1]*15000,0,channel_weights[1]*10000,45,0,0);
          //func->SetParameters(150000,0,150000,45,0,0);
          for (int j = 0;j < 30;j++){
            waveform[i]->Fit("fit","Q","",25,90);
          }
          calibrated_integrated_charge += func->GetParameter(0) + func->GetParameter(2);
          integrated_charge_APD[0]->Fill(func->GetParameter(0) + func->GetParameter(2));
          WLS[1] = func->GetParameter(2);
          spike[1] = func->GetParameter(0);
          //cout << i << "    " << func->GetParameter(0) + func->GetParameter(2) - 14500 << endl;
        }
        if (i == 10){
          waveform[i]->Smooth(1);
          func->SetParameters(channel_weights[14]*15000,0,channel_weights[14]*10000,45,0,0);
          //func->SetParameters(150000,0,150000,45,0,0);
          for (int j = 0;j < 30;j++){
            waveform[i]->Fit("fit","Q","",25,65);
          }
          /*difference_hist->Reset();
          difference_hist2->Reset();
          difference_hist->Add(waveform[i],1);
          difference_hist->Add(func,-1);
          for (int i = 0;i < 1024;i++){
            difference_hist2->SetBinContent(i+1,difference_hist->GetBinContent(i+42-int(func->GetParameter(3))));
          }
          //difference_hist2->Add(mean_noise_4APD[3],-1);
          //func->SetParameters(150000,0,150000,25,0,0);
          //func->SetParameters(100000,1,1,30,0,0);
          for (int j = 0;j < 50;j++){
            //difference_hist2->Fit("fit","Q","",25,125);
          }
          //difference_hist2->Add(func,-1);
          //difference_hist2->Draw();
          //waveform[i]->Add(func,-1);
          waveform[i]->Draw();
          //mean_noise_4APD[i-4]->Draw();
          //difference_hist->Add(waveform[i],1.0/3958);
          //difference_hist->Add(func,-1.0/3958);
          //for (int i = 0;i < 1024;i++){
            //difference_hist2->SetBinContent(i+1,difference_hist->GetBinContent(i+42-int(func->GetParameter(3))));
          //}
          //mean_noise_4APD[i-4]->Add(difference_hist2,1);
          //integrated_charge_APD[i-4]->Fill(func->GetParameter(2));*/
          //integrated_charge += func->GetParameter(0) + func->GetParameter(2);//100*func->GetMaximum();//
          calibrated_integrated_charge += func->GetParameter(0) + func->GetParameter(2);//100*func->GetMaximum();//
          WLS[14] = func->GetParameter(2);
          spike[14] = func->GetParameter(0);
          integrated_charge_APD[3]->Fill(func->GetParameter(0) + func->GetParameter(2));
          //cout << i << "    " << func->GetParameter(0) + func->GetParameter(2) - 14500 << endl;
        }
        /*else{
          func->SetParameters(40000,1,1,30,0,0);
          for (int j = 0;j < 30;j++){
            waveform[i]->Fit("fit","Q","",25,60);
          }
          calibrated_integrated_charge += func->GetParameter(2);
        }*/
        if (i == 8){
          waveform[i]->Smooth(1);
          func->SetParameters(channel_weights[10]*17000,0,channel_weights[10]*12000,45,0,0);
          for (int j = 0;j < 30;j++){
            waveform[i]->Fit("fit","Q","",25,90);
          }
          calibrated_integrated_charge += func->GetParameter(0) + func->GetParameter(2);
          integrated_charge_APD[2]->Fill(func->GetParameter(0) + func->GetParameter(2));
          WLS[10] = func->GetParameter(2);
          spike[10] = func->GetParameter(0);
          //cout << i << "    " << func->GetParameter(0) + func->GetParameter(2) - 14500 << endl;
        }
        if (i == 3){
          waveform[i]->Smooth(1);
          func->SetParameters(channel_weights[5]*4000,0,channel_weights[5]*1000,45,0,0);
          for (int j = 0;j < 30;j++){
            waveform[i]->Fit("fit","Q","",25,90);
          }
          calibrated_integrated_charge += func->GetParameter(0) + func->GetParameter(2);
          integrated_charge_APD[1]->Fill(func->GetParameter(0) + func->GetParameter(2));
          WLS[5] = func->GetParameter(2);
          spike[5] = func->GetParameter(0);
          //cout << i << "    " << func->GetParameter(0) + func->GetParameter(2) - 14500 << endl;
        }
        if (i == 2){
          waveform[i]->Smooth(1);
          func->SetParameters(channel_weights[2]*4000,0,channel_weights[2]*1000,45,0,0);
          for (int j = 0;j < 30;j++){
            waveform[i]->Fit("fit","Q","",25,90);
          }
          calibrated_integrated_charge += func->GetParameter(0) + func->GetParameter(2);
          integrated_charge_APD[1]->Fill(func->GetParameter(0) + func->GetParameter(2));
          WLS[2] = func->GetParameter(2);
          spike[2] = func->GetParameter(0);
          //cout << i << "    " << func->GetParameter(0) + func->GetParameter(2) << endl;
        }
        if (i == 9){
          waveform[i]->Smooth(1);
          func->SetParameters(channel_weights[13]*4000,0,channel_weights[13]*1000,50,0,0);
          for (int j = 0;j < 30;j++){
            waveform[i]->Fit("fit","Q","",25,70);
          }
          calibrated_integrated_charge += func->GetParameter(0) + func->GetParameter(2);
          integrated_charge_APD[2]->Fill(func->GetParameter(0) + func->GetParameter(2));
          WLS[13] = func->GetParameter(2);
          spike[13] = func->GetParameter(0);
          //cout << i << "    " << func->GetParameter(0) + func->GetParameter(2) << endl;
        }
        if (i == 11){
          waveform[i]->Smooth(1);
          func->SetParameters(channel_weights[15]*5500,0,channel_weights[15]*1000,45,0,0);
          for (int j = 0;j < 30;j++){
            waveform[i]->Fit("fit","Q","",25,80);
          }
          calibrated_integrated_charge += func->GetParameter(0) + func->GetParameter(2);
          integrated_charge_APD[3]->Fill(func->GetParameter(0) + func->GetParameter(2));
          WLS[15] = func->GetParameter(2);
          spike[15] = func->GetParameter(0);
          //cout << i << "    " << func->GetParameter(0) + func->GetParameter(2) << endl;
        }
        if (i == 0){
          waveform[i]->Smooth(100);
          func->SetParameters(channel_weights[15]*4000,0,channel_weights[15]*1000,55,0,0);
          for (int j = 0;j < 30;j++){
            waveform[i]->Fit("fit","Q","",25,80);
          }
          calibrated_integrated_charge += func->GetParameter(0) + func->GetParameter(2);
          integrated_charge_APD[0]->Fill(func->GetParameter(0) + func->GetParameter(2));
          WLS[0] = func->GetParameter(2);
          spike[0] = func->GetParameter(0);
          //cout << i << "    " << func->GetParameter(0) + func->GetParameter(2) << endl;
        }
        //cout << i << "    " << func->GetParameter(0) + func->GetParameter(2) << endl;
      }
      /*
      difference_hist->Reset();
      difference_hist->Add(waveform[0],1);
      difference_hist->Add(waveform[1],1);
      difference_hist->Add(waveform[2],1);
      difference_hist->Add(waveform[3],1);
      difference_hist->Add(waveform[8],1);
      difference_hist->Add(waveform[9],1);
      difference_hist->Add(waveform[10],1);
      difference_hist->Add(waveform[11],1);
      difference_hist->SetLineWidth(0);
      //difference_hist->Smooth(1000);
      func->SetParameters(300000,0.4,80,0,0);
      for (int j = 0;j < 10;j++){
        difference_hist->Fit("fit","Q","",20,90);
      }*/
      //integrated_charge = func->GetParameter(0);
      //calibrated_integrated_charge += 1.1*func->GetParameter(0);//100*func->GetMaximum();//
      //cout << integrated_charge << "      " << calibrated_integrated_charge << "      " << 100*(calibrated_integrated_charge-integrated_charge)/integrated_charge << endl;
      integrated_charge_crystal11->Fill(integrated_charge);
      calibrated_integrated_charge_crystal11->Fill(calibrated_integrated_charge);
      //difference_hist->Reset();
      //difference_hist2->Reset();
      difference_hist->SetLineWidth(0);
      if (2+2==5){
        difference_hist->Add(waveform[7],1);
        //difference_hist->Add(func,-1);
        for (int i = 0;i < 1024;i++){
          difference_hist2->SetBinContent(i+1,difference_hist->GetBinContent(i+42-int(func->GetParameter(3))));
        }
        difference_hist3->Add(difference_hist2,1);
        difference_hist->SetLineColor(4);
        //difference_hist->Smooth(1);
        for (int i = 0;i < 20;i++) difference_hist3->Fit("fit","","",20,130);
        difference_hist3->Add(func,-1);
        difference_hist3->Draw();
        //func->Draw("same");
      }
      if (2+2==5){
        difference_hist->Add(waveform[7],1.0/3958);
        difference_hist->Add(func,-1.0/3958);
        for (int i = 0;i < 1024;i++){
          difference_hist2->SetBinContent(i+1,difference_hist->GetBinContent(i+42-int(func->GetParameter(3))));
        }
        difference_hist3->Add(difference_hist2,1);
      }
      /*//if (TMath::Abs(calibrated_integrated_charge-577571) > 10000){
        cout << x[1]-207.0 << "    " << y[1]-294.5 << endl;
        cout << integrated_charge-485715 << endl;
        cout << calibrated_integrated_charge-577571 << endl;
        for (int i = 0;i < 12;i++){
          cout << i << "        " << charge[i] << endl;
        }
        waveform[11]->SetLineWidth(0);
        //waveform[11]->Smooth(100);
        waveform[8]->SetLineWidth(0);
        //waveform[8]->Smooth(100);
        difference_hist->Reset();
        difference_hist->Add(waveform[11],1);
        difference_hist->Add(waveform[8],1);
        difference_hist->Add(waveform[0],1);
        difference_hist->Add(waveform[1],1);
        difference_hist->Add(waveform[2],1);
        difference_hist->Add(waveform[9],1);
        difference_hist->Add(waveform[10],1);
        difference_hist->SetLineWidth(0);
        //difference_hist->Smooth(100);
        difference_hist->Draw();
        waveform[0]->SetLineWidth(0);
        waveform[4]->SetLineWidth(0);
        waveform[1]->SetLineWidth(0);
        waveform[2]->SetLineWidth(0);
        waveform[9]->SetLineWidth(0);
        waveform[10]->SetLineWidth(0);
        waveform[0]->Draw("same");
        waveform[1]->Draw("same");
        waveform[2]->Draw("same");
        waveform[8]->Draw("same");
        waveform[9]->Draw("same");
        waveform[10]->Draw("same");
        waveform[11]->Draw("same");
        waveform[4]->Draw("same");*/
        /*if (difference_hist->GetBinContent(50)<-40)*/ counter+=1; //3958 events in the 4*4 window
        //cout << endl << jentry << endl;
        //cout << integrated_charge - 1141620 << "    " << calibrated_integrated_charge - 1262960 << "    " << calibrated_integrated_charge-integrated_charge << endl;
        //waveform[0]->Draw();
        //if (counter > 40) break;
        //cout << counter << endl;
      //}
    }
    templates->Fill();
  }
  //TFile *mean_noise_4APDs = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/mean_noise_4APDs","recreate");
  //mean_noise_4APDs->cd();
  //for (int i = 0;i < 4;i++) mean_noise_4APD[i]->Write(); 
  //cout << counter << endl;
  /*double APD_mean[4];
  double APD_mean_error[4];
  double APD_sigma[4];
  double APD_sigma_error[4];
  double APD_resolution[4];
  double APD_resolution_error[4];
  g->SetParLimits(0,0,1000000);
  for (int i = 0;i < 4;i++){
    //crys->SetParameters(500,integrated_charge_APD[i]->GetMean(),(integrated_charge_APD[i]->GetStdDev())/4,1.1,2.1);
    g->SetParameters(50,integrated_charge_APD[i]->GetMean(),(integrated_charge_APD[i]->GetStdDev())/2);
    integrated_charge_APD[i]->Draw("same");
    integrated_charge_APD[i]->Fit("g");
    APD_mean[i] = g->GetParameter(1);
    APD_mean_error[i] = crys->GetParError(1);
    APD_sigma[i] = crys->GetParameter(2);
    APD_sigma_error[i] = crys->GetParError(2);
    APD_resolution[i] = 100.0 * APD_sigma[i] / APD_mean[i];
    APD_resolution_error[i] = APD_resolution[i] * TMath::Sqrt(TMath::Power(APD_mean_error[i]/APD_mean[i],2.0)+TMath::Power(APD_sigma_error[i]/APD_sigma[i],2.0));
    cout << "APD" << i+1 << " resolution:    " << APD_resolution[i] << " +- " << APD_resolution_error[i] << endl;
  }
  double APD_weight[4] = {1,1,1,1};
  for (int i = 0;i < 4;i++){
    //APD_weight[i] = (APD_mean[0]+APD_mean[1]+APD_mean[2]+APD_mean[3])/(4*APD_mean[i]);
    APD_weight[i] = (APD_mean[2])/(APD_mean[i]);
    cout << "APD" << i+1 << " weight:    " << APD_weight[i] << endl;
  }*/
  
  templates->Write();
  
  crys->SetParameters(1000,integrated_charge_crystal11->GetMean(),(integrated_charge_crystal11->GetStdDev())/4,1.1,2.1);
  integrated_charge_crystal11->Draw();
  integrated_charge_crystal11->Fit("crys");
  double weighted_mean = crys->GetParameter(1);
  double weighted_mean_error = crys->GetParError(1);
  double weighted_sigma = crys->GetParameter(2);
  double weighted_sigma_error = crys->GetParError(2);
  double weighted_resolution = 100.0 * weighted_sigma / weighted_mean;
  double weighted_resolution_error = weighted_resolution * TMath::Sqrt(TMath::Power(weighted_mean_error/weighted_mean,2.0)+TMath::Power(weighted_sigma_error/weighted_sigma,2.0));
  cout << "weighted resolution:  " << weighted_resolution << " +- " << weighted_resolution_error << endl;
  crys->SetParameters(800,calibrated_integrated_charge_crystal11->GetMean(),(calibrated_integrated_charge_crystal11->GetStdDev())/4,1.1,2.1);
  calibrated_integrated_charge_crystal11->Draw("same");
  calibrated_integrated_charge_crystal11->Fit("crys");
  double calibrated_mean = crys->GetParameter(1);
  double calibrated_mean_error = crys->GetParError(1);
  double calibrated_sigma = crys->GetParameter(2);
  double calibrated_sigma_error = crys->GetParError(2);
  double calibrated_resolution = 100.0 * calibrated_sigma / calibrated_mean;
  double calibrated_resolution_error = calibrated_resolution * TMath::Sqrt(TMath::Power(calibrated_mean_error/calibrated_mean,2.0)+TMath::Power(calibrated_sigma_error/calibrated_sigma,2.0));
  cout << "calibrated resolution:  " << calibrated_resolution << " +- " << calibrated_resolution_error << endl;
  crys->SetParameters(1000,integrated_charge_crystal11->GetMean(),(integrated_charge_crystal11->GetStdDev())/4,1.1,2.1);
}