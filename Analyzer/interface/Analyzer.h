//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Mar 20 17:23:21 2016 by ROOT version 5.34/18
// from TTree tree/tree
// found on file: files/totalOut_WLS_50GeV.root
//////////////////////////////////////////////////////////

#ifndef Analyzer_h
#define Analyzer_h

#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <sstream>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "Waveform.h"
//#include <Drawtools.h>

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TVectorD.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TProfile.h"
#include "TRandom3.h"

#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooFFTConvPdf.h"
#include "RooPlot.h"
#include "RooCBShape.h"
#include "RooCruijff.h"
#include "RooGenericPdf.h"
#include "RooTFnBinding.h"
#include "RooTFnPdfBinding.h"
#include  "RooNumConvPdf.h"


#include <TVectorD.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <string>
#include <map>
#include <algorithm>

//for fft
#include "/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/geant4/Analyzer/interface/Waveform.h"
#include "/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/geant4/Analyzer/interface/FFTConvolution.h"

// Fixed size dimensions of array or collections stored in the TTree if any.

class Analyzer {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   std::map<TString,TH1F*> histos_;
   std::map<TString,TH2F*> histos2D_;

   TString outDir_;
   std::string energy_;
   std::string setup_;

   TVectorD* meanValueEnergy;
   TVectorD* meanValueEnergyErr;

   TVectorD* resValueEnergy;
   TVectorD* resValueEnergyErr;

   TVectorD* meanValueEactLYScaled;
   TVectorD* meanValueEactLYScaledErr;

   TVectorD* resValueEactLYScaled;
   TVectorD* resValueEactLYScaledErr;


   TVectorD* resValueTime;
   TVectorD* resValueTimeErr;

   TVectorD* resValueTime_frac50;
   TVectorD* resValueTimeErr_frac50;

   TVectorD* resValueTime_frac50_shaped;
   TVectorD* resValueTimeErr_frac50_shaped;

   TVectorD* resValueTime_frac50_shaped_apd;
   TVectorD* resValueTimeErr_frac50_shaped_apd;
   

   // Declaration of leaf types
   Int_t           Event;
   std::vector<float>   *Time_deposit;
   std::vector<float>   *Eact_layer;
   std::vector<float>   *Process_deposit;
   std::vector<float>   *Z_deposit;
   std::vector<float>   *opPhoton_time;
   std::vector<int>     *opPhoton_process;
   Float_t         Eabs;
   Float_t         Eact;
   Float_t         EfibrCore;
   Float_t         EfibrClad;
   Int_t           nLayers;
   Int_t           Fibre_0;
   Int_t           Fibre_1;
   Int_t           Fibre_2;
   Int_t           Fibre_3;
   Float_t         xPosition;
   Float_t         yPosition;
   Float_t         EOpt_0;
   Float_t         EOpt_1;
   Float_t         EOpt_2;
   Float_t         EOpt_3;

   // List of branches
   TBranch        *b_Event;   //!
   TBranch        *b_Time_deposit;   //!
   TBranch        *b_Eact_layer;   //!
   TBranch        *b_Process_deposit;   //!
   TBranch        *b_Z_deposit;   //!
   TBranch        *b_opPhoton_time;   //!
   TBranch        *b_opPhoton_process;   //!
   TBranch        *b_Eabs;   //!
   TBranch        *b_Eact;   //!
   TBranch        *b_EfibrCore;   //!
   TBranch        *b_EfibrClad;   //!
   TBranch        *b_nLayers;   //!
   TBranch        *b_Fibre_0;   //!
   TBranch        *b_Fibre_1;   //!
   TBranch        *b_Fibre_2;   //!
   TBranch        *b_Fibre_3;   //!
   TBranch        *b_xPosition;   //!
   TBranch        *b_yPosition;   //!
   TBranch        *b_EOpt_0;   //!
   TBranch        *b_EOpt_1;   //!
   TBranch        *b_EOpt_2;   //!
   TBranch        *b_EOpt_3;   //!

   Analyzer(TTree *tree=0);
   virtual ~Analyzer();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(std::string setup, std::string energy);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   void addHisto(TString name, int nBins, float XLow, float XUp,TString XLabel);
   void addHisto2D(TString name, int nBinsX, float XLow, float XUp,TString XLabel,int nBinsY, float YLow, float YUp,TString YLabel);
   void createHistos();
   void drawHistos();
   void fitHisto(TH1F* histo, TVectorD* res, TVectorD* resErr,bool isEnergy=false,TVectorD* mean=0, TVectorD* meanErr=0);
   std::vector<float> createLY(int nLayers);
   template <typename T>  std::vector<size_t> sort_indexes(const std::vector<T> &v);
};

#endif

#ifdef Analyzer_cxx
Analyzer::Analyzer(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("files/totalOut_WLS_50GeV.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("files/totalOut_WLS_50GeV.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

Analyzer::~Analyzer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Analyzer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Analyzer::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Analyzer::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

  meanValueEnergy = new TVectorD(1);
  meanValueEnergyErr = new TVectorD(1);

  resValueEnergy = new TVectorD(1);
  resValueEnergyErr = new TVectorD(1);

  meanValueEactLYScaled = new TVectorD(1);
  meanValueEactLYScaledErr = new TVectorD(1);

  resValueEactLYScaled = new TVectorD(1);
  resValueEactLYScaledErr = new TVectorD(1);


  resValueTime = new TVectorD(1);
  resValueTimeErr = new TVectorD(1);

  resValueTime_frac50 = new TVectorD(1);
  resValueTimeErr_frac50 = new TVectorD(1);

  resValueTime_frac50_shaped = new TVectorD(1);
  resValueTimeErr_frac50_shaped = new TVectorD(1);

  resValueTime_frac50_shaped_apd = new TVectorD(1);
  resValueTimeErr_frac50_shaped_apd = new TVectorD(1);


   // Set object pointer
   Time_deposit = 0;
   Eact_layer = 0;
   Process_deposit = 0;
   Z_deposit = 0;
   opPhoton_time = 0;
   opPhoton_process = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Event", &Event, &b_Event);
   fChain->SetBranchAddress("Time_deposit", &Time_deposit, &b_Time_deposit);
   fChain->SetBranchAddress("Eact_layer", &Eact_layer, &b_Eact_layer);
   fChain->SetBranchAddress("Process_deposit", &Process_deposit, &b_Process_deposit);
   fChain->SetBranchAddress("Z_deposit", &Z_deposit, &b_Z_deposit);
   fChain->SetBranchAddress("opPhoton_time", &opPhoton_time, &b_opPhoton_time);
   fChain->SetBranchAddress("opPhoton_process", &opPhoton_process, &b_opPhoton_process);
   fChain->SetBranchAddress("Eabs", &Eabs, &b_Eabs);
   fChain->SetBranchAddress("Eact", &Eact, &b_Eact);
   fChain->SetBranchAddress("EfibrCore", &EfibrCore, &b_EfibrCore);
   fChain->SetBranchAddress("EfibrClad", &EfibrClad, &b_EfibrClad);
   fChain->SetBranchAddress("nLayers", &nLayers, &b_nLayers);
   fChain->SetBranchAddress("Fibre_0", &Fibre_0, &b_Fibre_0);
   fChain->SetBranchAddress("Fibre_1", &Fibre_1, &b_Fibre_1);
   fChain->SetBranchAddress("Fibre_2", &Fibre_2, &b_Fibre_2);
   fChain->SetBranchAddress("Fibre_3", &Fibre_3, &b_Fibre_3);
   fChain->SetBranchAddress("xPosition", &xPosition, &b_xPosition);
   fChain->SetBranchAddress("yPosition", &yPosition, &b_yPosition);
   fChain->SetBranchAddress("EOpt_0", &EOpt_0, &b_EOpt_0);
   fChain->SetBranchAddress("EOpt_1", &EOpt_1, &b_EOpt_1);
   fChain->SetBranchAddress("EOpt_2", &EOpt_2, &b_EOpt_2);
   fChain->SetBranchAddress("EOpt_3", &EOpt_3, &b_EOpt_3);
   Notify();
}

Bool_t Analyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Analyzer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Analyzer::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

template <typename T> 
std::vector<size_t> Analyzer::sort_indexes(const std::vector<T> &v) {

  // initialize original index locations
  std::vector<size_t> idx(v.size());
  for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;

}

#endif // #ifdef Analyzer_cxx
