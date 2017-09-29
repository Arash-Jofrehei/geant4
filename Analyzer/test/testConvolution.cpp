#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <sstream>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

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


//#include <DrawTools.h>
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
#include "Waveform.h"
#include "WaveformNew.h"

#include <vector>
#include <string>
#include <map>
#include <algorithm>

#include "FFTConvolution.h"

int main( int argc, char* argv[] ) {

  //DrawTools::setStyle();

  std::string setup = "WLS";
  std::string energy = "50GeV";

  if( argc>1 ) {
    std::string setup_str(argv[1]);
    setup = setup_str;
  }else{
    std::cout << "Usage:" << std::endl;
    std::cout << "./testConvolution ([setup]) " << std::endl;
    exit(12345);
  }

  std::cout<<"plotting for setup:"<<setup<<std::endl;

  std::map<TString,float> energies;
  energies["100GeV"]=100;
  TH1F* waveform;

  Waveform* waveF;
  WaveformNew* waveFNew; //FIXME

  for(std::map<TString,float>::const_iterator it = energies.begin();it!=energies.end(); it++){
    TFile* file = TFile::Open("outFiles/plotterTiming_"+setup+"_"+it->first+".root");
    waveform = (TH1F*) file->Get("timeArrival");
    TGraph* waveGraph = new TGraph(waveform);
    waveF = new Waveform(waveGraph->GetN(),waveGraph->GetX(),waveGraph->GetY());
    waveFNew = new WaveformNew(waveGraph->GetN(),waveGraph->GetX(),waveGraph->GetY());

    //    waveform->Sumw2();
    //    waveform->Scale(1./waveform->Integral());
  }




  RooRealVar t("t","t", 0, 200);
  t.setBins(1024,"cache") ;
  t.setBinning(RooBinning(1024,0,200));
  RooDataHist wave("wave","dataset with t",t,RooFit::Import(*waveform) );

  RooHistPdf histpdf1("histpdf1","histpdf1",t,wave);


  TF1 *f = new  TF1("f", "pow(exp(1.)*[1]*(x)/[0],[0])*exp(-[1]*(x))",0.,200.);
  f->SetParameters(3.0, 0.8);



  RooAbsPdf* rfa1 = RooFit::bindPdf(f,t);
  //try with fft convolution in time is multiplication in frequency
  TH1* histoResponse = rfa1->createHistogram("t",1024);
  TGraph* waveGraphResp = new TGraph(histoResponse);
  Waveform* waveResp = new Waveform(histoResponse->GetNbinsX(),waveGraphResp->GetX(),waveGraphResp->GetY());
  WaveformNew* waveRespNew = new WaveformNew(histoResponse->GetNbinsX(),waveGraphResp->GetX(),waveGraphResp->GetY());

  std::cout<<"BEFORE"<<std::endl;
  //  TH1F* histConvBef = waveF->get_histo("waveF");//FIXME!check if it's ok before fft

  waveFNew->fft();
  TH1F* histConvBef = waveFNew->get_histo("waveF");//FIXME!check if it's ok before fft


  //  waveF->inv_fft();

  waveResp->fft();

  //  waveF->multiply_fft(waveResp);//FIXME implement correctly
  //  waveF->inv_fft(0,30);
  //  waveF->inv_fft();
  //  waveF->inv_fft_multiply();
  //  waveF->inv_fft(100,10);
  //  waveF->print();

  FFTConvolution* convProd = new FFTConvolution();
  WaveformNew* convWave = convProd->fftConvolute(waveFNew,waveRespNew);
  TH1F* histConvNew = convWave->get_histo("waveFNew");
  histConvNew->SetLineColor(kGreen+2);


  TFile* file2 = TFile::Open("circuitResponse.root");
  TGraph* graphCircuit = (TGraph*) file2->Get("WaveFormCircuit");
  WaveformNew* waveCircuit = new WaveformNew(graphCircuit->GetN(),graphCircuit->GetX(),graphCircuit->GetY());
  FFTConvolution* convProd2 = new FFTConvolution();
  WaveformNew* convWave2 = convProd2->fftConvolute(convWave,waveCircuit);
  //    WaveformNew* convWave2 = convProd2->fftConvolute(convWave,waveFNew);
  TH1F* histConvNew2 = convWave2->get_histo("waveFNew2");
  histConvNew2->SetLineColor(kRed);

  

  std::cout<<"------------------------AFTER-------------------------"<<std::endl;
  TH1F* histConv = waveF->get_histo("waveF");//FIXME!check if it's ok before fft
  histConv->SetLineColor(kViolet);
  //  RooRealVar tau("tau","tau",1.22451e+7,1.e+5,pow(10,30)) ;
  //  RooRealVar timeshift("timeshift","timeshift",0.03e-6,0,0.05e-6) ;

  //  RooGenericPdf crrc("crrc","(@0-0.05e-6)*@1*exp(-(@0-0.05e-6)*@1)",RooArgSet(t,tau)) ;


  //  TFile* outFile = TFile::Open("outFiles/testConvolution_"+setup+".root");

  std::string fileName(Form("outFiles/testConvolution_%s",setup.c_str()));
  fileName+=".root";
  TFile* outFile = TFile::Open(fileName.c_str(),"recreate");

  //  RooNumConvPdf convolutionShaping("convolutionShaping","exp (X) crrc",t,histpdf1,*rfa1) ;
  //  RooDataHist& cacheHist = histpdf1.dataHist() ;
  //  RooRealVar* histX = (RooRealVar*) cacheHist.get()->find("t") ;
  for(int i=0;i<1056;++i)  {
    //  histX->setBin(i);   
    //    std::cout<<"histX"<<histX->getBinning().binCenter(i)<<" ";
    std::cout<<"histX"<<t.getBinning().binNumber(0);
    t.setBin(i);   
    std::cout<<" -----histpdf1: "<<i<<" "<<histpdf1.getVal(wave.get())<<std::endl;
  }
  RooFFTConvPdf convolutionShaping("convolutionShaping","exp (X) crrc",t,histpdf1,*rfa1) ;
  RooRealVar center("center","center",50,0,400);
  RooRealVar width("width","width",40,0,400);

  //convolutionShaping.setConvolutionWindow(center,width,1);
  convolutionShaping.verboseEval(0);
  std::cout<<convolutionShaping.bufferStrategy()<<"-------------- strategy-----"<<std::endl;
  //  convolutionShaping.convIntConfig().setEpsAbs(0.001);//absolute precision
  //  convolutionShaping.convIntConfig().setEpsRel(0.001);//relative precision

  RooPlot* frame1 = t.frame(RooFit::Title("histogram pdf")) ;
  //  rfa1->plotOn(frame1);
//  //  std::cout<<"norm: "<<waveform->Integral()<<std::endl;
  wave.plotOn(frame1,RooFit::LineColor(kBlack),RooFit::DrawOption("L"),RooFit::LineWidth(2),RooFit::Precision(0.1));
  float binWidth = histConvNew->GetXaxis()->GetBinWidth(1);
  std::string ytitle = Form("Response [a.u.]",binWidth);
  frame1->SetYTitle(ytitle.c_str());
  //  convolutionShaping.plotOn(frame1,RooFit::LineColor(kRed)) ;

  TH1* histo_shaped =  convolutionShaping.createHistogram("t");
  histo_shaped->Sumw2();
  histo_shaped->SetName("histo_shaped");


  TH1* hdummy=wave.createHistogram("t",histConv->GetNbinsX());
  hdummy->SetLineColor(kRed);
  TCanvas c1("c_histpdf");

  //  frame1->Draw();
  TGaxis::SetMaxDigits(1); 
  waveform->GetXaxis()->SetNoExponent(kTRUE);  
  waveform->GetYaxis()->SetTitle(ytitle.c_str());
  waveform->Scale(1./waveform->GetMaximum());
  histConvNew->Scale(1./histConvNew->GetMaximum());
  histConvNew2->Scale(1./histConvNew2->GetMaximum());
  waveform->Draw("histL");
  std::cout<<"NBINS:"<<histConv->GetNbinsX()<<" "<<hdummy->GetNbinsX()<<std::endl;
  std::cout<<hdummy->Integral()<<" "<<histConv->Integral()<<std::endl;
  histConv->Scale(hdummy->Integral(0,200)/histConv->Integral(0,200));
  ///  histConvNew->Scale(hdummy->Integral(0,200)/histConvNew->Integral(0,200));
  //  hdummy->DrawNormalized("L");
  //  histo_shaped->DrawNormalized("Lsame");
  //    histConv->Draw("same");

  histConvNew->Draw("same");
  histConvNew2->Draw("same");

  c1.SaveAs("histpdf.png");
  c1.SaveAs("histpdf.pdf");

  c1.Write();
  outFile->Write();
  histo_shaped->Write();
  hdummy->Write("hdummy");
  histConv->Write("histConv");
  histConvNew->Write("histConvNew");
  histConvNew2->Write("histConvNew2");
  histConvBef->Write("histConvBef");
  outFile->Close();

  return 0;

}
