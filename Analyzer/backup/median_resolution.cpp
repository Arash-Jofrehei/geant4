#include <iostream>
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
using namespace std;
float median_integrated_charge[18];
float median_x[2];
float median_y[2];
int median_nFibersX[2];
int median_nFibersY[2];
float channel_weights[18] = {2.41098,1.62535,1.49631,0,0,0.470364,0.954885,0.931846,1.15933,0.983237,1.3819,0,0,0.452602,1.86884,1.41348,0,0};//{1,1,1,0,0,1,1,1,1,1,1,0,0,1,1,1,0,0};
float E = 1;

void median_resolution(){
  TCanvas *canvas = new TCanvas("median_resolution","median_resolution");
  TFile *final = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/analysis_4684.root");
  TTree *ftree = (TTree*) final->Get("median");
  ftree->SetBranchAddress("median_integrated_charge",median_integrated_charge);
  ftree->SetBranchAddress("median_x",median_x);
  ftree->SetBranchAddress("median_y",median_y);
  ftree->SetBranchAddress("median_nFibersX",median_nFibersX);
  ftree->SetBranchAddress("median_nFibersY",median_nFibersY);
  TF1 *crys = new TF1("crys","crystalball",0,5000000);
  TF1 *g = new TF1("g","gaus",0,5000000);
  TH1F *integrated_charge_APD[4];
  integrated_charge_APD[0] = new TH1F("integrated charge APD1","integrated charge APD1",100*E,1000*E,500000*E);
  integrated_charge_APD[1] = new TH1F("integrated charge APD2","integrated charge APD2",100*E,1000*E,500000*E);
  integrated_charge_APD[2] = new TH1F("integrated charge APD3","integrated charge APD3",100*E,1000*E,500000*E);
  integrated_charge_APD[3] = new TH1F("integrated charge APD4","integrated charge APD4",100*E,1000*E,500000*E);
  TH1F *charge = new TH1F("charge","charge",200,50000,2000000);
  TH1F *calibrated_charge = new TH1F("calibrated_charge","calibrated_charge",200,50000,2000000);
  float integrated_charge;
  float calibrated_integrated_charge;
  Long64_t nentries = ftree->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++){
    Long64_t ientry = ftree->LoadTree(jentry);
    nb = ftree->GetEntry(jentry);   nbytes += nb;
    if (TMath::Abs(median_x[1]-207)>2||TMath::Abs(median_y[1]-294.5)>2||median_nFibersX[1]>2||median_nFibersY[1]>2) continue;
    integrated_charge = 0;
    for (int i = 6;i < 10;i++) integrated_charge += channel_weights[i] * median_integrated_charge[i];
    charge->Fill(integrated_charge);
    calibrated_integrated_charge = 0;
    for (int i = 0;i < 18;i++) calibrated_integrated_charge += channel_weights[i] * median_integrated_charge[i];
    calibrated_charge->Fill(calibrated_integrated_charge);
    //for (int i = 0;i < 4;i++) integrated_charge_APD[i]->Fill(median_integrated_charge[i+6]);
    /*integrated_charge_APD[0]->Fill(median_integrated_charge[0]);
    integrated_charge_APD[1]->Fill(median_integrated_charge[2]);
    integrated_charge_APD[2]->Fill(median_integrated_charge[13]);
    integrated_charge_APD[3]->Fill(median_integrated_charge[15]);*/
  }
  
  //    --------------------------------------------------------------------------------------------------------------
  /*
  double APD_mean[4];
  double APD_mean_error[4];
  double APD_sigma[4];
  double APD_sigma_error[4];
  double APD_resolution[4];
  double APD_resolution_error[4];
  for (int i = 0;i < 4;i++){
    crys->SetParameters(500,integrated_charge_APD[i]->GetMean(),(integrated_charge_APD[i]->GetStdDev())/4,1.1,2.1);
    g->SetParameters(50,integrated_charge_APD[i]->GetMean(),(integrated_charge_APD[i]->GetStdDev())/2);
    integrated_charge_APD[i]->Draw("same");
    //integrated_charge_APD[i]->Fit("g");
    integrated_charge_APD[i]->Fit("crys");
    //APD_mean[i] = g->GetParameter(1);
    APD_mean[i] = crys->GetParameter(1);
    //APD_mean_error[i] = g->GetParError(1);
    APD_mean_error[i] = crys->GetParError(1);
    //APD_sigma[i] = g->GetParameter(2);
    APD_sigma[i] = crys->GetParameter(2);
    //APD_sigma_error[i] = g->GetParError(2);
    APD_sigma_error[i] = crys->GetParError(2);
    APD_resolution[i] = 100.0 * APD_sigma[i] / APD_mean[i];
    APD_resolution_error[i] = APD_resolution[i] * TMath::Sqrt(TMath::Power(APD_mean_error[i]/APD_mean[i],2.0)+TMath::Power(APD_sigma_error[i]/APD_sigma[i],2.0));
    cout << "APD" << i+1 << " resolution:    " << APD_resolution[i] << " +- " << APD_resolution_error[i] << endl;
  }
  double APD_weight[4] = {1,1,1,1};
  for (int i = 0;i < 4;i++){
    APD_weight[i] = (APD_mean[0]+APD_mean[1]+APD_mean[2]+APD_mean[3])/(4*APD_mean[i]);
    //APD_weight[i] = (APD_mean[2])/(APD_mean[i]);
    cout << "APD" << i+1 << " weight:    " << APD_weight[i] << endl;
  }
  */
  //    --------------------------------------------------------------------------------------------------------------
  
  crys->SetParameters(500,charge->GetMean(),(charge->GetStdDev())/4,1.1,2.1);
  charge->Draw();
  charge->Fit("crys");
  double weighted_mean = crys->GetParameter(1);
  double weighted_mean_error = crys->GetParError(1);
  double weighted_sigma = crys->GetParameter(2);
  double weighted_sigma_error = crys->GetParError(2);
  double weighted_resolution = 100.0 * weighted_sigma / weighted_mean;
  double weighted_resolution_error = weighted_resolution * TMath::Sqrt(TMath::Power(weighted_mean_error/weighted_mean,2.0)+TMath::Power(weighted_sigma_error/weighted_sigma,2.0));
  cout << "weighted resolution:  " << weighted_resolution << " +- " << weighted_resolution_error << endl;
  crys->SetParameters(800,calibrated_charge->GetMean(),(calibrated_charge->GetStdDev())/4,1.1,2.1);
  calibrated_charge->Draw("same");
  calibrated_charge->Fit("crys");
  double calibrated_mean = crys->GetParameter(1);
  double calibrated_mean_error = crys->GetParError(1);
  double calibrated_sigma = crys->GetParameter(2);
  double calibrated_sigma_error = crys->GetParError(2);
  double calibrated_resolution = 100.0 * calibrated_sigma / calibrated_mean;
  double calibrated_resolution_error = calibrated_resolution * TMath::Sqrt(TMath::Power(calibrated_mean_error/calibrated_mean,2.0)+TMath::Power(calibrated_sigma_error/calibrated_sigma,2.0));
  cout << "calibrated resolution:  " << calibrated_resolution << " +- " << calibrated_resolution_error << endl;
}