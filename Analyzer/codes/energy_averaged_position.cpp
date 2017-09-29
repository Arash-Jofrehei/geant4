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
float median_charge[18];
float subtracted_median_charge[18];
float x[2];
float y[2];
int nFibresOnX[2];
int nFibresOnY[2];
double xbar;
double ybar;
float originalX[2];
float originalY[2];
double Xbar;
double Ybar;
double uniformed_calibrated_sum;
int x_index,y_index;
const int nbins = 44;
double c[nbins][nbins];
double calibrated_sum = 0;
//float weights[18] = {2.41098,1.62535,1.49631,0,0,0.470364,0.954885,0.931846,1.15933,0.983237,1.3819,0,0,0.452602,1.86884,1.41348,0,0};
float channel_weights_crystal11[18] = {0,1,1,1,0,0,0.327*0.954885,0.327*0.931846,0.327*1.15933,0.327*0.983237,1,0.957,0,0,1,1.06,1,0};

void energy_averaged_position(){
  TCanvas *canvas = new TCanvas("energy averaged position","energy averaged position");
  TH1F *uniformed_calibrated_response = new TH1F("uniformed calibrated response","uniformed calibrated response crystal 11",210,0.3,1.7);
  TH3F *Ybar_Xbar_uni = new TH3F("Ybar_Xbar_uni","Ybar_Xbar_uni",nbins,-11,11,nbins,-11,11,12000,0,12);
  TH3F *Y1_X1_uni = new TH3F("Y1_X1_uni","Y1_X1_uni",nbins,-11,11,nbins,-11,11,12000,0,12);
  TProfile *Y1_uni = new TProfile("Y1_uni","uniformed calibrated response;Y(mm)",33,-11,11);
  TFile *merged_file = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/merged_crystal11.root"/*,"UPDATE"*/);
  //TFile *merged_file = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/analysis_4684.root"/*,"UPDATE"*/);
  TTree *ftree = (TTree*) merged_file->Get("h4");
  TTree *uniformed = new TTree("uniformed","uniformed");
  uniformed->Branch("uniformed_calibrated_sum",&uniformed_calibrated_sum,"uniformed_calibrated_sum/D");
  uniformed->Branch("originalX",originalX,"originalX[2]/F");
  uniformed->Branch("originalY",originalY,"originalY[2]/F");
  uniformed->Branch("Xbar",&Xbar,"Xbar/D");
  uniformed->Branch("Ybar",&Ybar,"Ybar/D");
  ftree->SetBranchAddress("median_charge",median_charge);
  ftree->SetBranchAddress("subtracted_median_charge",subtracted_median_charge);
  ftree->SetBranchAddress("x",x);
  ftree->SetBranchAddress("y",y);
  ftree->SetBranchAddress("xbar",&xbar);
  ftree->SetBranchAddress("ybar",&ybar);
  ftree->SetBranchAddress("nFibresOnX",nFibresOnX);
  ftree->SetBranchAddress("nFibresOnY",nFibresOnY);
  
  TH1F *amp[nbins][nbins];
  string name;
  for(int i = 0;i < nbins;i++){
    for(int j = 0;j < nbins;j++){
      name = "amp";
      name.append("_");
      name.append(to_string(i));
      name.append("_");
      name.append(to_string(j));
      const char *amp_name = name.c_str();
      amp[i][j] = new TH1F(amp_name,amp_name,60,10000,1000000);
    }
  }
  
  Long64_t nentries = ftree->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = ftree->LoadTree(jentry);
    nb = ftree->GetEntry(jentry);   nbytes += nb;
    calibrated_sum = 0;
    for (int i = 0;i < 18;i++){
      if(i < 6 || i > 9) calibrated_sum += median_charge[i] * channel_weights_crystal11[i];
      else calibrated_sum += median_charge[i] * channel_weights_crystal11[i];
    }
    x_index = int((nbins/22.0)*(xbar + 11));
    y_index = int((nbins/22.0)*(ybar + 11));
    if (x_index<0||y_index<0||x_index>(nbins-1)||y_index>(nbins-1) || median_charge[10] < 0.1)
      continue;
    amp[x_index][y_index]->Fill(calibrated_sum);
  }
  
  for (int i = 0;i < nbins;i++){
    for (int j = 0;j < nbins;j++){
      c[i][j] = 1.0/(amp[i][j]->GetXaxis()->GetBinCenter(amp[i][j]->GetMaximumBin()));
      //cout << i << "    " << j << "    " << 1.0/c[i][j] << endl;
    }
  }
  
  nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++){
    Long64_t ientry = ftree->LoadTree(jentry);
    nb = ftree->GetEntry(jentry);   nbytes += nb;
    if (jentry%10000 == 0) cout << jentry << " / " << nentries << endl;
    originalX[0] = x[0];
    originalY[0] = y[0];
    originalX[1] = x[1];
    originalY[1] = y[1];
    Xbar = xbar;
    Ybar = ybar;
    calibrated_sum = 0;
    x_index = int((nbins/22.0)*(xbar + 11));
    y_index = int((nbins/22.0)*(ybar + 11));
    if (x_index<0||y_index<0||x_index>(nbins-1)||y_index>(nbins-1) || median_charge[10] < 0.1){
      uniformed_calibrated_sum = 0;
      //uniformed->Fill();
      continue;
    }
    for (int i = 0;i < 18;i++){
      if(i < 6 || i > 9) calibrated_sum += median_charge[i] * channel_weights_crystal11[i];
      else calibrated_sum += median_charge[i] * channel_weights_crystal11[i];
    }
    uniformed_calibrated_sum = calibrated_sum * c[x_index][y_index];
    if (TMath::Abs(x[1])<2 && TMath::Abs(y[1])<2) uniformed_calibrated_response->Fill(uniformed_calibrated_sum);
    if (x[1] > -500 && y[1] > -500){
      Ybar_Xbar_uni->Fill(xbar,ybar,uniformed_calibrated_sum);
      Y1_X1_uni->Fill(x[1],y[1],uniformed_calibrated_sum);
      if (TMath::Abs(x[1])<2) Y1_uni->Fill(y[1],uniformed_calibrated_sum);
    }
    //uniformed->Fill();
  }
  TProfile2D *Ybar_Xbar_uni_prof = Ybar_Xbar_uni->Project3DProfile("yx");
  Ybar_Xbar_uni_prof->SetTitle("uniformed calibrated response crystal 11;energy averaged X(mm);energy averaged Y(mm)");
  TProfile2D *Y1_X1_uni_prof = Y1_X1_uni->Project3DProfile("yx");
  Y1_X1_uni_prof->SetTitle("uniformed calibrated response crystal 11;X(mm);Y(mm)");
  Ybar_Xbar_uni_prof->Draw("colz");
  TCanvas *canvas2 = new TCanvas("energy averaged position2","energy averaged position2");
  Y1_X1_uni_prof->Draw("colz");
  TCanvas *canvas3 = new TCanvas("energy averaged position3","energy averaged position3");
  uniformed_calibrated_response->Draw();
  TCanvas *canvas4 = new TCanvas("energy averaged position4","energy averaged position4");
  Y1_uni->Draw();
  //uniformed->Write();
}