#include <iostream>
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
using namespace std;
int nbins = 19;
float x[2];
float y[2];
int nFibresOnX[2];
int nFibresOnY[2];
TBranch        *b_nFibresOnX;
TBranch        *b_nFibresOnY;
float charge_sig[18];
float maximum[18];

double FWHM(TProfile *h){
  int bin1 = h->FindFirstBinAbove(h->GetMaximum()/2);
  int bin2 = h->FindLastBinAbove(h->GetMaximum()/2);
  //int bin1 = h->GetMaximumBin();
  float fwhm = h->GetBinCenter(bin2) - h->GetBinCenter(bin1);
  return fwhm;
}

void plotting(){
  TCanvas *canvas = new TCanvas("plotting","plotting");
  //TFile *simulated_waveforms_position_profiles = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/simulated_waveforms_position_profiles.root");
  TFile *simulated_waveforms_position_profiles = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/ideal_surface_simulated_waveforms_position_profiles.root");
  //TFile *simulated_waveforms_position_profiles = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/full_yield_simulated_waveforms_position_profiles.root");
  TProfile2D *simulated_WLS_max_4APDs = (TProfile2D*) simulated_waveforms_position_profiles->Get("4APDs WLS waveform amplitude");
  TProfile2D *simulated_WLS_FWHM_4APDs = (TProfile2D*) simulated_waveforms_position_profiles->Get("4APDs WLS waveform FWHM");
  TFile *simulated_WLS_fscint_deposit = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/simulated_WLS_fscint_deposit_ideal_surface.root");
  TProfile2D *WLS_count = (TProfile2D*) simulated_WLS_fscint_deposit->Get("# WLS optical photons with geant4 simulation (4APDs)");
  TProfile2D *fscint_count = (TProfile2D*) simulated_WLS_fscint_deposit->Get("# fiber scintillation optical photons with geant4 simulation (4APDs)");
  TProfile2D *act_deposit = (TProfile2D*) simulated_WLS_fscint_deposit->Get("Active Material Deposited Energy - Central Crystal");
  TFile *data_waveforms_FWHM_file = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/data_waveforms_FWHM.root");
  TFile *apd_profile_waveform_bins = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/apd_profile_waveform_19bins_shift.root");
  TFile *final = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/merged_crystal4apd.root");
  TTree *ftree = (TTree*) final->Get("h4");
  ftree->SetBranchAddress("x",x);
  ftree->SetBranchAddress("y",y);
  ftree->SetBranchAddress("nFibresOnX",nFibresOnX);
  ftree->SetBranchAddress("nFibresOnY",nFibresOnY);
  ftree->SetBranchAddress("charge_sig",charge_sig);
  ftree->SetBranchAddress("maximum",maximum);
  
  TH2F *data_waveform_FWHM[5];
  data_waveform_FWHM[0] = new TH2F("4APDs data waveform FWHM","4APDs data waveform FWHM;X(mm);Y(mm)",nbins,-9.5,9.5,nbins,-9.5,9.5);
  data_waveform_FWHM[1] = new TH2F("APD1 data waveform FWHM","APD1 data waveform FWHM;X(mm);Y(mm)",nbins,-9.5,9.5,nbins,-9.5,9.5);
  data_waveform_FWHM[2] = new TH2F("APD2 data waveform FWHM","APD2 data waveform FWHM;X(mm);Y(mm)",nbins,-9.5,9.5,nbins,-9.5,9.5);
  data_waveform_FWHM[3] = new TH2F("APD3 data waveform FWHM","APD3 data waveform FWHM;X(mm);Y(mm)",nbins,-9.5,9.5,nbins,-9.5,9.5);
  data_waveform_FWHM[4] = new TH2F("APD4 data waveform FWHM","APD4 data waveform FWHM;X(mm);Y(mm)",nbins,-9.5,9.5,nbins,-9.5,9.5);
  
  TProfile *waveform_dummy = new TProfile("Something","Something",1024,-0.1,204.7);
  TProfile *apd[4][nbins][nbins];
  string name;
         //reading the waveform profiles for each bin
  for(int i = 0;i < nbins;i++){
    for(int j = 0;j < nbins;j++){
      for (int a = 0;a < 4;a++){    
        if (a == 0) waveform_dummy->Reset();
        name = "APD";
        name.append(to_string(a+1));
        name.append("_");
        name.append(to_string(i));
        name.append("_");
        name.append(to_string(j));
        const char *APD_name = name.c_str();
        apd[a][i][j] = (TProfile*) apd_profile_waveform_bins->Get(APD_name);
        waveform_dummy->Add(apd[a][i][j]);
        data_waveform_FWHM[a]->SetBinContent(i+1,j+1,FWHM(apd[a][i][j]));
        if (a == 3) data_waveform_FWHM[0]->SetBinContent(i+1,j+1,FWHM(waveform_dummy));
      }
    }
  }
  
  TProfile2D *charge = new TProfile2D("integrated charge from data","integrated charge from data;X(mm);Y(mm)",38,-9.5,9.5,38,-9.5,9.5);
  TProfile2D *amp_max = new TProfile2D("maximum amplitude from data","maximum amplitude from data;X(mm);Y(mm)",38,-9.5,9.5,38,-9.5,9.5);
  Long64_t nentries = ftree->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  float sum_max = 0;
  float sum_charge = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++){
    Long64_t ientry = ftree->LoadTree(jentry);
    nb = ftree->GetEntry(jentry);   nbytes += nb;
    if (jentry%10000 == 0) cout << jentry << "/" << nentries << endl;
    if (x[1]<0 || y[1]<0) continue;
    sum_charge = charge_sig[6]+charge_sig[7]+charge_sig[8]+charge_sig[9];
    sum_max = maximum[6]+maximum[7]+maximum[8]+maximum[9];
    if (sum_charge < 10000) continue;
    if (sum_max < 800) continue;
    charge->Fill(x[1]-207,y[1]-294.5,sum_charge);
    amp_max->Fill(x[1]-207,y[1]-294.5,sum_max);
  }
  
  //                                                plotting some amplitude and FWHM ratios
  
  TH2F *ratio_WLS_deposit = new TH2F("ratio of simulated WLS over deposited energy","ratio of simulated WLS over deposited energy;X(mm);Y(mm)",76,-9.5,9.5,76,-9.5,9.5);
  TH2F *ratio_amp_data_WLS = new TH2F("ratio of amplitudes - data over simulated WLS","ratio of amplitudes - data over simulated WLS;X(mm);Y(mm)",38,-9.5,9.5,38,-9.5,9.5);
  TH2F *ratio_charge_data_WLS = new TH2F("ratio of integrated charge - data over simulated WLS","ratio of integrated charge - data over simulated WLS;X(mm);Y(mm)",38,-9.5,9.5,38,-9.5,9.5);
  TH2F *ratio_FWHM_data_WLS = new TH2F("ratio of FWHM - data to simulated WLS","ratio of FWHM - data to simulated WLS;X(mm);Y(mm)",19,-9.5,9.5,19,-9.5,9.5);
  TH2F *ratio_amp_data_depos = new TH2F("ratio of data amplitude over simulated deposited energy","ratio of data amplitude over simulated deposited energy;X(mm);Y(mm)",38,-9.5,9.5,38,-9.5,9.5);
  TH2F *ratio_charge_data_depos = new TH2F("ratio of data integrated charge over simulated deposited energy","ratio of data integrated charge over simulated deposited energy;X(mm);Y(mm)",38,-9.5,9.5,38,-9.5,9.5);
  for (int i = 0;i < 38;i++){
    for (int j = 0;j < 38;j++){
      float average_max = 0;
      float average_charge = 0;
      float average_FWHM = 0;
      float average_depos = 0;
      float n = 0;
      if (simulated_WLS_max_4APDs->GetBinContent(2*i+1,2*j+1) > 0){
        n += 1;
        average_max += simulated_WLS_max_4APDs->GetBinContent(2*i+1,2*j+1);
      }
      if (simulated_WLS_max_4APDs->GetBinContent(2*i+2,2*j+1) > 0){
        n += 1;
        average_max += simulated_WLS_max_4APDs->GetBinContent(2*i+2,2*j+1);
      }
      if (simulated_WLS_max_4APDs->GetBinContent(2*i+1,2*j+2) > 0){
        n += 1;
        average_max += simulated_WLS_max_4APDs->GetBinContent(2*i+1,2*j+2);
      }
      if (simulated_WLS_max_4APDs->GetBinContent(2*i+2,2*j+2) > 0){
        n += 1;
        average_max += simulated_WLS_max_4APDs->GetBinContent(2*i+2,2*j+2);
      }
      
      if (WLS_count->GetBinContent(2*i+1,2*j+1) > 0){
        average_charge += WLS_count->GetBinContent(2*i+1,2*j+1);
      }
      if (WLS_count->GetBinContent(2*i+2,2*j+1) > 0){
        average_charge += WLS_count->GetBinContent(2*i+2,2*j+1);
      }
      if (WLS_count->GetBinContent(2*i+1,2*j+2) > 0){
        average_charge += WLS_count->GetBinContent(2*i+1,2*j+2);
      }
      if (WLS_count->GetBinContent(2*i+2,2*j+2) > 0){
        average_charge += WLS_count->GetBinContent(2*i+2,2*j+2);
      }
      /*if (i < 19 && j < 19){
        if (simulated_WLS_FWHM_4APDs->GetBinContent(4*i+1,4*j+1) > 0){
          average_FWHM += simulated_WLS_FWHM_4APDs->GetBinContent(4*i+1,4*j+1);
        }
        if (simulated_WLS_FWHM_4APDs->GetBinContent(4*i+2,4*j+1) > 0){
          average_FWHM += simulated_WLS_FWHM_4APDs->GetBinContent(4*i+2,4*j+1);
        }
        if (simulated_WLS_FWHM_4APDs->GetBinContent(4*i+1,4*j+2) > 0){
          average_FWHM += simulated_WLS_FWHM_4APDs->GetBinContent(4*i+1,4*j+2);
        }
        if (simulated_WLS_FWHM_4APDs->GetBinContent(4*i+2,4*j+2) > 0){
          average_FWHM += simulated_WLS_FWHM_4APDs->GetBinContent(4*i+2,4*j+2);
        }
      }*/
      
      if (act_deposit->GetBinContent(2*i+1,2*j+1) > 0){
        average_depos += act_deposit->GetBinContent(2*i+1,2*j+1);
      }
      if (act_deposit->GetBinContent(2*i+2,2*j+1) > 0){
        average_depos += act_deposit->GetBinContent(2*i+2,2*j+1);
      }
      if (act_deposit->GetBinContent(2*i+1,2*j+2) > 0){
        average_depos += act_deposit->GetBinContent(2*i+1,2*j+2);
      }
      if (act_deposit->GetBinContent(2*i+2,2*j+2) > 0){
        average_depos += act_deposit->GetBinContent(2*i+2,2*j+2);
      }
      average_max /= n;
      average_charge /= n;
      //average_FWHM /= n;
      average_depos /= n;
      if (average_max == 0 || average_charge == 0/* || average_FWHM == 0*/) continue;
      float dummy = amp_max->GetBinContent(i+1,j+1)/average_max;
      //if (dummy > 500) continue;
      ratio_amp_data_WLS->SetBinContent(i+1,j+1,dummy);
      ratio_charge_data_WLS->SetBinContent(i+1,j+1,charge->GetBinContent(i+1,j+1)/average_charge);
      //if (i < 19 && j < 19)ratio_FWHM_data_WLS->SetBinContent(i+1,j+1,data_waveform_FWHM[0]->GetBinContent(i+1,j+1)/average_FWHM);
      ratio_amp_data_depos->SetBinContent(i+1,j+1,amp_max->GetBinContent(i+1,j+1)/average_depos);
      ratio_charge_data_depos->SetBinContent(i+1,j+1,charge->GetBinContent(i+1,j+1)/average_depos);
    }
  }
  
  for (int i = 0;i < 76;i++){
    for (int j = 0;j < 76;j++){
      ratio_WLS_deposit->SetBinContent(i+1,j+1,simulated_WLS_max_4APDs->GetBinContent(i+1,j+1)/act_deposit->GetBinContent(i+1,j+1));
    }
  }
  
  simulated_WLS_FWHM_4APDs->Rebin2D(4,4);
  for (int i = 0;i < 19;i++){
    for (int j = 0;j < 19;j++){
      ratio_FWHM_data_WLS->SetBinContent(i+1,j+1,data_waveform_FWHM[0]->GetBinContent(i+1,j+1)/simulated_WLS_FWHM_4APDs->GetBinContent(i+1,j+1));
    }
  }
  
  
  ratio_WLS_deposit->Scale(1.0/ratio_WLS_deposit->GetBinContent(38,38));
  for (int i = 0;i < 76;i++) for (int j = 0;j < 76;j++) if (ratio_WLS_deposit->GetBinContent(i+1,j+1) > 3) ratio_WLS_deposit->SetBinContent(i+1,j+1,0);
  
  ratio_amp_data_WLS->Scale(1.0/ratio_amp_data_WLS->GetBinContent(19,19));
  for (int i = 0;i < 38;i++) for (int j = 0;j < 38;j++) if (ratio_amp_data_WLS->GetBinContent(i+1,j+1) > 3) ratio_amp_data_WLS->SetBinContent(i+1,j+1,0);
  
  ratio_charge_data_WLS->Scale(1.0/ratio_charge_data_WLS->GetBinContent(19,19));
  for (int i = 0;i < 38;i++) for (int j = 0;j < 38;j++) if (ratio_charge_data_WLS->GetBinContent(i+1,j+1) > 3) ratio_charge_data_WLS->SetBinContent(i+1,j+1,0);
  
  ratio_FWHM_data_WLS->Scale(1.0/ratio_FWHM_data_WLS->GetBinContent(10,10));
  for (int i = 0;i < 19;i++) for (int j = 0;j < 19;j++) if (ratio_FWHM_data_WLS->GetBinContent(i+1,j+1) > 1.2) ratio_FWHM_data_WLS->SetBinContent(i+1,j+1,0);
  
  ratio_amp_data_depos->Scale(1.0/ratio_amp_data_depos->GetBinContent(19,19));
  for (int i = 0;i < 38;i++) for (int j = 0;j < 38;j++) if (ratio_amp_data_depos->GetBinContent(i+1,j+1) > 3) ratio_amp_data_depos->SetBinContent(i+1,j+1,0);
  
  ratio_charge_data_depos->Scale(1.0/ratio_charge_data_depos->GetBinContent(19,19));
  for (int i = 0;i < 38;i++) for (int j = 0;j < 38;j++) if (ratio_charge_data_depos->GetBinContent(i+1,j+1) > 3) ratio_charge_data_depos->SetBinContent(i+1,j+1,0);
  
  /*TFile *data_waveforms_FWHM_file = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/data_waveforms_FWHM.root","recreate");
  data_waveforms_FWHM_file->cd();
  for (int i = 0;i < 5;i++) data_waveform_FWHM[i]->Write();*/
  
  ratio_WLS_deposit->Draw("colz");
  TCanvas *canvas2 = new TCanvas("plotting2","plotting2");
  ratio_amp_data_WLS->Draw("colz");
  TCanvas *canvas3 = new TCanvas("plotting3","plotting3");
  ratio_charge_data_WLS->Draw("colz");
  TCanvas *canvas4 = new TCanvas("plotting4","plotting4");
  ratio_FWHM_data_WLS->Draw("colz");
  TCanvas *canvas5 = new TCanvas("plotting5","plotting5");
  ratio_amp_data_depos->Draw("colz");
  data_waveform_FWHM[0]->Draw("colz");
  TCanvas *canvas6 = new TCanvas("plotting6","plotting6");
  ratio_charge_data_depos->Draw("colz");
  simulated_WLS_FWHM_4APDs->Draw("colz");
}