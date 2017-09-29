#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <sstream>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

//#include "DrawTools.h"
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

#include <vector>
#include <string>
#include <map>
#include <algorithm>

float getRatioError( float num, float denom, float numErr, float denomErr ) {

  return sqrt( numErr*numErr/(denom*denom) + denomErr*denomErr*num*num/(denom*denom*denom*denom) );

}


int main( int argc, char* argv[] ) {

  //DrawTools::setStyle();

  std::string setup = "WLS";
  std::string energy = "50GeV";

  if( argc>1 ) {
    std::string setup_str(argv[1]);
    setup = setup_str;
  }else{
    std::cout << "Usage:" << std::endl;
    std::cout << "./compareTimingEnergy ([setup]) " << std::endl;
    exit(12345);
  }

  std::cout<<"plotting for setup:"<<setup<<std::endl;

  std::map<TString,float> energies;
  energies["20GeV"]=20;
  energies["50GeV"]=50;
  energies["100GeV"]=100;
  energies["200GeV"]=200;
  
  TH2F * axes = new TH2F("axes","axes",200,0.,250.,300,0,400);
 
  TH2F * axes_2 = new TH2F("axes_2","axes_2",200,0.,250.,300,0,4);
  
  TGraphErrors* resVsEnergy = new TGraphErrors(0);
  resVsEnergy->SetName("resVsEnergy");
  resVsEnergy->SetMarkerStyle(20);
  resVsEnergy->SetMarkerSize(1.6);
  resVsEnergy->SetMarkerColor(kBlue);

  TGraphErrors* resEnVsEnergy = new TGraphErrors(0);
  resEnVsEnergy->SetName("resEnVsEnergy");
  resEnVsEnergy->SetMarkerStyle(20);
  resEnVsEnergy->SetMarkerSize(1.6);
  resEnVsEnergy->SetMarkerColor(kBlue);


  axes->GetYaxis()->SetTitle("Timing Resolution [ps]");
  axes->GetXaxis()->SetTitle("Energy [GeV]");

  axes_2->GetYaxis()->SetTitle("Energy Resolution [%]");
  axes_2->GetXaxis()->SetTitle("Energy [GeV]");


  int i=0;
  for(std::map<TString,float>::const_iterator it = energies.begin();it!=energies.end(); it++){
    TFile* file = TFile::Open("outFiles/plotterTiming_"+setup+"_"+it->first+".root");
    TVectorD* resValue = (TVectorD*)file->Get("resValueTime_frac50");
    TVectorD* resValueErr = (TVectorD*) file->Get("resValueTimeErr_frac50");

    TVectorD* resValueEnergy = (TVectorD*)file->Get("resValueEnergy");
    TVectorD* resValueErrEnergy = (TVectorD*) file->Get("resValueEnergyErr");

    TVectorD* meanValueEnergy = (TVectorD*)file->Get("meanValueEnergy");
    TVectorD* meanValueErrEnergy = (TVectorD*) file->Get("meanValueEnergyErr");


    std::cout<<"processing file:"<<i<<std::endl;
    std::cout<< (*resValue)[0]<<" "<<(*resValueEnergy)[0]<<" "<<it->first<<std::endl;
    resVsEnergy->SetPoint(i,it->second, (*resValue)[0]);
    resVsEnergy->SetPointError(i,0, (*resValueErr)[0]);

     resEnVsEnergy->SetPoint(i,it->second, 100*(*resValueEnergy)[0]/(*meanValueEnergy)[0]);
     resEnVsEnergy->SetPointError(i,0, getRatioError((*resValueEnergy)[0],(*meanValueEnergy)[0],(*resValueErrEnergy)[0],(*meanValueErrEnergy)[0]));


    ++i;  
 }

  std::cout<<"daje"<<std::endl;

  TCanvas c1;
  TF1* f = new TF1("f","[0]/sqrt(x)",0,250);
  axes->Draw();
  resVsEnergy->Fit(f,"R",0,250);
  resVsEnergy->Draw("Psame");
  c1.SaveAs("outFiles/plots/resVsEnergy.pdf");
  c1.SaveAs("outFiles/plots/resVsEnergy.png");

  c1.Clear();
  TF1* f2 = new TF1("f2","[0]/sqrt(x)",0,250);
  axes_2->Draw();
  resEnVsEnergy->Fit(f2,"R",0,250);
  resEnVsEnergy->Draw("Psame");
  c1.SaveAs("outFiles/plots/resEnVsEnergy.pdf");
  c1.SaveAs("outFiles/plots/resEnVsEnergy.png");


  return 0;

}
