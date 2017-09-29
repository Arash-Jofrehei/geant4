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
int nbins_WLS = 76;

string ratio_plots(){
  TCanvas *canvas = new TCanvas("ratio","ratio");
  TFile *file1 = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/simulated_WLS_fscint_deposit_badlyPolished.root");
  TProfile2D *count_WLS_4APDs_1 = (TProfile2D*) file1->Get("# WLS optical photons with geant4 simulation (4APDs)");
  TFile *file2 = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/simulated_WLS_fscint_deposit_ideal_surface.root");
  TProfile2D *count_WLS_4APDs_2 = (TProfile2D*) file2->Get("# WLS optical photons with geant4 simulation (4APDs)");
  TH2F *ratio_count_WLS = new TH2F("ratio of WLS photons: after fixing surfaces / before fixing surfaces","ratio of WLS photons: after fixing surfaces / before fixing surfaces;X(mm);Y(mm)",nbins_WLS,-9.5,9.5,nbins_WLS,-9.5,9.5);
  for (int i = 0; i < nbins_WLS ;i++){
    for (int j = 0; j < nbins_WLS ;j++){
      if (count_WLS_4APDs_1->GetBinContent(i+1,j+1) != 0) ratio_count_WLS->SetBinContent(i+1,j+1,count_WLS_4APDs_2->GetBinContent(i+1,j+1)/count_WLS_4APDs_1->GetBinContent(i+1,j+1));
    }
  }
  ratio_count_WLS->Draw("colz");
  return "TLI";
}