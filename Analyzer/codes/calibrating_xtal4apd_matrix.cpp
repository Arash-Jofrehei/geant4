#include <iostream>
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "string.h"
#include "TMatrixDSym.h"
using namespace std;
float median_integrated_charge[18];
float median_x[2];
float median_y[2];
int median_nFibersX[2];
int median_nFibersY[2];
float charge_sig[18];
float x[2];
float y[2];
Int_t           nFibresOnX[2];
Int_t           nFibresOnY[2];
TBranch        *b_nFibresOnX;
TBranch        *b_nFibresOnY;
float WLS[18];
float spike[18];
float positionX[2];
float positionY[2];
float nFibersX[2];
float nFibersY[2];
float corner_percent = 0.5;
float side_percent = 1.78;
float channel_weights[18] = {2.41098,1.62535,1.49631,0,0,0.470364,0.954885,0.931846,1.15933,0.983237,1.3819,0,0,0.452602,1.86884,1.41348,0,0};//{1.33715,1.37695,0.862342,0,0,0.382559,0.9088*0.951252,0.9088*0.932984,0.9088*1.16544,0.9088*0.981471,1.14502,0,0,0.236246,1.55053,0.769922,0,0};//{3.08814,2.16522,2.06918,0,0,0.603557,0.84*0.951252,0.84*0.932984,0.84*1.16544,0.84*0.981471,1.80394,0,0,0.563872,2.44467,1.90858,0,0};//{3.08814,2.16522,2.06918,0,0,0.603557,0.84*0.951252,0.84*0.932984,0.84*1.16544,0.84*0.981471,1.80394,0,0,0.563872,2.44467,1.90858,0,0};//{3.16812,2.1125,2.05083,0,0,0.601527,0.84*0.951252,0.84*0.932984,0.84*1.16544,0.84*0.981471,1.69442,0,0,0.564436,2.31684,1.83204,0,0};//{0.35*5.80376,0.7*3.52714,0.35*4.26522,0,0,0.7*1,0.951252,0.932984,1.16544,0.981471,0.7*2.79305,0,0,0.35*1,0.7*3.80542,0.35*3.21488,0,0};//{21.2314,-23.811,10.5023,0,0,-0.178681,0.951252,0.932984,1.16544,0.981471,-1.01619,0,0,-17.4042,16.457,-37.735};
void calibrating_xtal4apd_matrix(){
  TCanvas *canvas = new TCanvas("calibrating xtal4apd matrix","calibrating xtal4apd matrix");
  TFile *file = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/analysis_4684.root");
  TTree *tree = (TTree*) file->Get("templates");
  TTree *main = (TTree*) file->Get("h4");
  TTree *ftree = (TTree*) file->Get("median");
  main->SetBranchAddress("x", x);
  main->SetBranchAddress("y", y);
  main->SetBranchAddress("nFibresOnX", nFibresOnX, &b_nFibresOnX);
  main->SetBranchAddress("nFibresOnY", nFibresOnY, &b_nFibresOnY);
  main->SetBranchAddress("charge_sig",charge_sig);
  tree->SetBranchAddress("WLS",WLS);
  tree->SetBranchAddress("spike",spike);
  tree->SetBranchAddress("positionX",positionX);
  tree->SetBranchAddress("positionY",positionY);
  tree->SetBranchAddress("nFibersX",nFibersX);
  tree->SetBranchAddress("nFibersY",nFibersY);
  ftree->SetBranchAddress("median_integrated_charge",median_integrated_charge);
  ftree->SetBranchAddress("median_x",median_x);
  ftree->SetBranchAddress("median_y",median_y);
  ftree->SetBranchAddress("median_nFibersX",median_nFibersX);
  ftree->SetBranchAddress("median_nFibersY",median_nFibersY);
  int x_index = -100;
  int y_index = -100;
  float E = 1;
  float template_amp;
  float calibrated_template_amp;
  TH1F *template_amplitude_crystal4apd = new TH1F("template amplitude crystal 4APD","crystal 4APD template amplitude (single crystal)",1800,100*E,1800000*E);
  TH1F *calibrated_template_amplitude_crystal4apd = new TH1F("calibrated template amplitude crystal 4APD","crystal 4APD template amplitude (3*3 matrix)",1800,100,1800000*E);
  TH1F *calibration_hists[100][100];
  string name;
  for (int i = 0; i < 100;i++){
    for (int j = 0; j < 100;j++){
      name = "calib_bin";
      name.append("_");
      name.append(to_string(i));
      name.append("_");
      name.append(to_string(j));
      const char *bin_name = name.c_str();
      calibration_hists[i][j] = new TH1F(bin_name,bin_name,100,100*E,2000000*E);
    }
  }
  TH2F *threeGroupsMin = new TH2F("calibration constant probe","resolution for 100GeV as a function of ratios of energy deposit (percent)",100,0,2,100,0,4);
  TH1D *bin_amp[9][9];
  for (int i = 0;i < 9;i++){
    for (int j = 0;j < 9;j++){
      name = "bin";
      name.append("_");
      name.append(to_string(i));
      name.append("_");
      name.append(to_string(j));
      const char *bin_name = name.c_str();
      if(j == 4) bin_amp[i][j] = new TH1D(bin_name,bin_name,130*E,100*E,1600000);
      else bin_amp[i][j] = new TH1D(bin_name,bin_name,100,1000*E,100000);
    }
  }/*
  Long64_t nentries = tree->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++){
    Long64_t ientry = tree->LoadTree(jentry);
    nb = tree->GetEntry(jentry);   nbytes += nb;
    if (TMath::Abs(positionX[1]-207)<4&&TMath::Abs(positionY[1]-294.5)<4&&nFibersX[1]<3&&nFibersY[1]<3){
      template_amp = 0;
      calibrated_template_amp = 0;
      for (int i = 0;i < 18;i++){
        calibrated_template_amp += channel_weights[i]*(WLS[i] + spike[i]);
        if (i >= 6 && i <= 9) template_amp += channel_weights[i]*(WLS[i] + spike[i]);
      }
      template_amplitude_crystal4apd->Fill(template_amp);
      calibrated_template_amplitude_crystal4apd->Fill(calibrated_template_amp);
      for (int i = 0;i < 100;i++){
        for (int j = 0;j < 100;j++){
          calibrated_template_amp = 0;
          for (int k = 0;k < 18;k++){
            if (k == 0 || k == 2 || k == 13 || k == 15) calibrated_template_amp += (0.02*i/corner_percent)*channel_weights[k]*(WLS[k] + spike[k]);
            if (k == 1 || k == 5 || k == 10 || k == 14) calibrated_template_amp += (0.04*j/side_percent)*channel_weights[k]*(WLS[k] + spike[k]);
            if (k == 6 || k == 7 || k == 8 || k == 9) calibrated_template_amp += (100-4*0.02*i-4*0.04*j)*0.01*channel_weights[k]*(WLS[k] + spike[k])/0.9088;
          }
          calibration_hists[i][j]->Fill(calibrated_template_amp);
        }
      }
      //x_index = int(positionX[1]-205.5);
      //y_index = int(positionY[1]-293);
      x_index = int((positionX[1]-204)/2);
      y_index = int((positionY[1]-291.5)/2);
      if (x_index < 0 || x_index > 2 || y_index < 0 || y_index > 2) continue;
      bin_amp[3*y_index + x_index][0]->Fill(WLS[0] + spike[0]);
      bin_amp[3*y_index + x_index][1]->Fill(WLS[1] + spike[1]);
      bin_amp[3*y_index + x_index][2]->Fill(WLS[2] + spike[2]);
      bin_amp[3*y_index + x_index][3]->Fill(WLS[5] + spike[5]);
      bin_amp[3*y_index + x_index][4]->Fill(0.951252*(WLS[6] + spike[6])+0.932984*(WLS[7] + spike[7])+1.16544*(WLS[8] + spike[8])+0.981471*(WLS[9] + spike[9]));
      bin_amp[3*y_index + x_index][5]->Fill(WLS[10] + spike[10]);
      bin_amp[3*y_index + x_index][6]->Fill(WLS[13] + spike[13]);
      bin_amp[3*y_index + x_index][7]->Fill(WLS[14] + spike[14]);
      bin_amp[3*y_index + x_index][8]->Fill(WLS[15] + spike[15]);
    }
  }*/
  //Long64_t nentries = main->GetEntriesFast();
  Long64_t nentries = ftree->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++){
    //Long64_t ientry = main->LoadTree(jentry);
    Long64_t ientry = ftree->LoadTree(jentry);
    //nb = main->GetEntry(jentry);   nbytes += nb;
    nb = ftree->GetEntry(jentry);   nbytes += nb;
    //if (TMath::Abs(x[1]-207)<4&&TMath::Abs(y[1]-294.5)<4&&nFibresOnX[1]<3&&nFibresOnY[1]<3){
    if (TMath::Abs(median_x[1]-207)<1&&TMath::Abs(median_y[1]-294.5)<1&&median_nFibersX[1]<3&&median_nFibersY[1]<3){
      template_amp = 0;
      calibrated_template_amp = 0;
      for (int i = 0;i < 18;i++){
        calibrated_template_amp += channel_weights[i]*(charge_sig[i]);
        if (i >= 6 && i <= 9) template_amp += channel_weights[i]*(charge_sig[i]);
      }
      template_amplitude_crystal4apd->Fill(template_amp);
      calibrated_template_amplitude_crystal4apd->Fill(calibrated_template_amp);
      for (int i = 0;i < 100;i++){
        for (int j = 0;j < 100;j++){
          calibrated_template_amp = 0;
          for (int k = 0;k < 18;k++){
            /*if (k == 0 || k == 2 || k == 13 || k == 15) calibrated_template_amp += (0.02*i/corner_percent)*channel_weights[k]*(charge_sig[k]);
            if (k == 1 || k == 5 || k == 10 || k == 14) calibrated_template_amp += (0.04*j/side_percent)*channel_weights[k]*(charge_sig[k]);
            if (k == 6 || k == 7 || k == 8 || k == 9) calibrated_template_amp += (100-4*0.02*i-4*0.04*j)*0.01*channel_weights[k]*(charge_sig[k])/0.9088;*/
            if (k == 0 || k == 2 || k == 13 || k == 15) calibrated_template_amp += (0.02*i/corner_percent)*channel_weights[k]*(median_integrated_charge[k]);
            if (k == 1 || k == 5 || k == 10 || k == 14) calibrated_template_amp += (0.04*j/side_percent)*channel_weights[k]*(median_integrated_charge[k]);
            if (k == 6 || k == 7 || k == 8 || k == 9) calibrated_template_amp += (100-4*0.02*i-4*0.04*j)*0.01*channel_weights[k]*(median_integrated_charge[k])/0.9088;
          }
          calibration_hists[i][j]->Fill(calibrated_template_amp);
        }
      }
      //x_index = int(positionX[1]-205.5);
      //y_index = int(positionY[1]-293);
      x_index = int((positionX[1]-204)/2);
      y_index = int((positionY[1]-291.5)/2);
      if (x_index < 0 || x_index > 2 || y_index < 0 || y_index > 2) continue;
      bin_amp[3*y_index + x_index][0]->Fill(charge_sig[0]);
      bin_amp[3*y_index + x_index][1]->Fill(charge_sig[1]);
      bin_amp[3*y_index + x_index][2]->Fill(charge_sig[2]);
      bin_amp[3*y_index + x_index][3]->Fill(charge_sig[5]);
      bin_amp[3*y_index + x_index][4]->Fill(0.951252*(charge_sig[6])+0.932984*(charge_sig[7])+1.16544*(charge_sig[8])+0.981471*(charge_sig[9]));
      bin_amp[3*y_index + x_index][5]->Fill(charge_sig[10]);
      bin_amp[3*y_index + x_index][6]->Fill(charge_sig[13]);
      bin_amp[3*y_index + x_index][7]->Fill(charge_sig[14]);
      bin_amp[3*y_index + x_index][8]->Fill(charge_sig[15]);
    }
  }
  TF1 *crys = new TF1("crys","crystalball",0,5000000);
  TF1 *g = new TF1("g","gaus",0,5000000);
  crys->SetParameters(600,template_amplitude_crystal4apd->GetMean(),(template_amplitude_crystal4apd->GetStdDev())/4,1.1,2.1);
  template_amplitude_crystal4apd->Draw();
  template_amplitude_crystal4apd->Fit("crys","Q");
  template_amplitude_crystal4apd->Fit("crys");
  double weighted_mean = crys->GetParameter(1);
  double weighted_mean_error = crys->GetParError(1);
  double weighted_sigma = crys->GetParameter(2);
  double weighted_sigma_error = crys->GetParError(2);
  double weighted_resolution = 100.0 * weighted_sigma / weighted_mean;
  double weighted_resolution_error = weighted_resolution * TMath::Sqrt(TMath::Power(weighted_mean_error/weighted_mean,2.0)+TMath::Power(weighted_sigma_error/weighted_sigma,2.0));
  cout << "weighted resolution:  " << weighted_resolution << " +- " << weighted_resolution_error << endl;
  crys->SetParameters(600,calibrated_template_amplitude_crystal4apd->GetMean(),(calibrated_template_amplitude_crystal4apd->GetStdDev())/4,1.1,2.1);
  calibrated_template_amplitude_crystal4apd->Draw("same");
  calibrated_template_amplitude_crystal4apd->Fit("crys","Q");
  calibrated_template_amplitude_crystal4apd->Fit("crys");
  double calibrated_mean = crys->GetParameter(1);
  double calibrated_mean_error = crys->GetParError(1);
  double calibrated_sigma = crys->GetParameter(2);
  double calibrated_sigma_error = crys->GetParError(2);
  double calibrated_resolution = 100.0 * calibrated_sigma / calibrated_mean;
  double calibrated_resolution_error = calibrated_resolution * TMath::Sqrt(TMath::Power(calibrated_mean_error/calibrated_mean,2.0)+TMath::Power(calibrated_sigma_error/calibrated_sigma,2.0));
  cout << "calibrated resolution:  " << calibrated_resolution << " +- " << calibrated_resolution_error << endl;
  for (int i = 0;i < 100;i++){
    for (int j = 0;j < 100;j++){
      crys->SetParameters(600,calibration_hists[i][j]->GetMean(),(calibration_hists[i][j]->GetStdDev())/4,1.1,2.1);
      for (int k = 0;k < 3;k++) calibration_hists[i][j]->Fit("crys","Q");
      threeGroupsMin->SetBinContent(i+1,j+1,100*crys->GetParameter(2)/crys->GetParameter(1));
    }
  }
  threeGroupsMin->Draw("colz");
  double center_bin_weights[9];
  double bin_means[81];
  for (int i = 0; i < 9;i++){
    for (int j = 0;j < 9;j++){
      crys->SetParameters(40,bin_amp[i][j]->GetMean(),(bin_amp[i][j]->GetStdDev())/4,1.1,2.1);
      if (j == 4){
        if (i == 8) crys->SetParameter(1,1260000);
        bin_amp[i][j]->Fit("crys","Q");
        bin_means[9*i+j] = crys->GetParameter(1);
      }
      else{
        bin_amp[i][j]->Fit("g","Q");
        bin_means[9*i+j] = g->GetParameter(1);
      }
      if (i==4){
        center_bin_weights[j] = bin_means[9*i+j];
      }
    }
  }
  /*TMatrixDSym binMeans(9,bin_means);
  binMeans.Print();
  binMeans.Invert();
  binMeans.Print();
  double dummy = 0;
  double calibration_constants[9];
  for (int i = 0; i < 9;i++) dummy += binMeans[4][i];
  for (int i = 0; i < 9;i++){
    calibration_constants[i] = 0;
    for (int j = 0;j < 9;j++){
      calibration_constants[i] += binMeans[i][j]/dummy;
    }
    cout << i << "    " << calibration_constants[i] << endl << endl;
  }
  for (int i = 0; i < 9;i++){
    if (i == 0 || i == 2 || i == 6 || i == 8) cout << center_bin_weights[i] << "        " << 0.01*corner_percent * center_bin_weights[4] / center_bin_weights[i] << endl;
    if (i == 1 || i == 3 || i == 5 || i == 7) cout << center_bin_weights[i] << "        " << 0.01*side_percent * center_bin_weights[4] / center_bin_weights[i] << endl;
    if (i == 4) cout << center_bin_weights[i] << "    "  << 1-0.04*corner_percent-0.04*side_percent << endl;
  }*/
}