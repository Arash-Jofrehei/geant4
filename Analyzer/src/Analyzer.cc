#define Analyzer_cxx
#include "Analyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#define DOSHAPING 1
#define LIMIT_NENTRIES 1

void Analyzer::addHisto(TString name, int nBins, float XLow, float XUp, TString XLabel){

  TH1F* histo = new TH1F(name,name,nBins,XLow,XUp);
  histo->GetXaxis()->SetTitle(XLabel);
  histos_[name]=histo;
}

void Analyzer::addHisto2D(TString name, int nBinsX, float XLow, float XUp,TString XLabel,int nBinsY, float YLow, float YUp,TString YLabel){
  TH2F* histo = new TH2F(name,name,nBinsX,XLow,XUp,nBinsY,YLow,YUp);
  histo->GetXaxis()->SetTitle(XLabel);
  histo->GetYaxis()->SetTitle(YLabel);
  histos2D_[name]=histo;
}



void Analyzer::drawHistos(){


  TCanvas c1;
   for(std::map<TString,TH1F*>::const_iterator it = histos_.begin();it!=histos_.end(); it++){
     std::cout<<outDir_+it->first+".png"<<std::endl;
     it->second->Draw();

     c1.SaveAs(outDir_+it->first+".png");
     c1.SaveAs(outDir_+it->first+".pdf");
     it->second->Write();
     c1.Clear();
   }

   for(std::map<TString,TH2F*>::const_iterator it2D = histos2D_.begin();it2D!=histos2D_.end(); it2D++){
     it2D->second->Draw("colz");
     TProfile* prof =  it2D->second->ProfileX("prof");
     prof->Draw("same");
     //     it2D->second->Draw("profXsame");
     std::cout<<outDir_+it2D->first+".png"<<std::endl;
     c1.SaveAs(outDir_+it2D->first+".png");
     c1.SaveAs(outDir_+it2D->first+".pdf");
     it2D->second->Write();
     c1.Clear();
   }


}

void Analyzer::createHistos(){
  if(setup_=="SingleFibre" || setup_=="Ideal2016")  addHisto("timeArrival_avg",500,5,10,"time [ns]");
  else addHisto("timeArrival_avg",500,10,20,"time [ns]");
  addHisto("ZTimeArrival_avg",420,0,420,"Z [mm]");
  addHisto("processTimeArrival_avg",3,0.5,3.5,"process");
 
  addHisto("nPhotons",1000,0,50000,"WLS Photons");

  if(setup_=="SingleFibre" || setup_=="Ideal2016")   addHisto("timeArrival_avg_wls",500,5,10,"time [ns]");
  else   addHisto("timeArrival_avg_wls",500,10,20,"time [ns]");
  addHisto("ZTimeArrival_avg_wls",420,0,420,"Z [mm]");
  
  if(setup_=="SingleFibre" || setup_=="Ideal2016")   addHisto("timeArrival_avg_scint",500,5,10,"time [ns]");
  else   addHisto("timeArrival_avg_scint",500,10,20,"time [ns]");
  addHisto("ZTimeArrival_avg_scint",420,0,420,"Z [mm]");
  
  if(setup_=="SingleFibre" || setup_=="Ideal2016")   addHisto("timeArrival_avg_cher",500,5,10,"time [ns]");
  else   addHisto("timeArrival_avg_cher",500,10,20,"time [ns]");
  addHisto("ZTimeArrival_avg_cher",420,0,420,"Z [mm]");

  if(setup_=="SingleFibre" || setup_=="Ideal2016")   addHisto("timeArrival",500,5,10,"time [ns]");
  else   addHisto("timeArrival",1024,0,200,"time [ns]");
  addHisto("ZTimeArrival",420,0,420,"Z [mm]");
  addHisto("processTimeArrival",3,0.5,3.5,"process");

  if(setup_=="SingleFibre" || setup_=="Ideal2016")   addHisto("timeArrival_scint",500,5,10,"time [ns]");
  else   addHisto("timeArrival_scint",500,10,20,"time [ns]");
  addHisto("ZTimeArrival_scint",420,0,420,"Z [mm]");

  if(setup_=="SingleFibre" || setup_=="Ideal2016")   addHisto("timeArrival_cher",500,5,10,"time [ns]");
  else   addHisto("timeArrival_cher",500,10,20,"time [ns]");
  addHisto("ZTimeArrival_cher",420,0,420,"Z [mm]");

  if(setup_=="SingleFibre" || setup_=="Ideal2016")   addHisto("timeArrival_wls",500,5,10,"time [ns]");
  else   addHisto("timeArrival_wls",500,10,20,"time [ns]");
  addHisto("ZTimeArrival_wls",420,0,420,"Z [mm]");

 addHisto("waveform_total",200,0,100,"time [ns]");
 addHisto("waveform_scint",200,0,100,"time [ns]");
 addHisto("waveform_cher",200,0,100,"time [ns]");
 addHisto("waveform_wls",200,0,100,"time [ns]");

 addHisto("time_at_max",100,0,50,"time [ns]");

 addHisto("Z_nPhotTiming",420,0,420,"Z [mm]");
 addHisto("Z_nPhotTimingNorm",420,0,420,"Z [mm]");
 addHisto("Z_firstPhotons",420,0,420,"Z [mm]");
 addHisto("Z_firstPhotonsNorm",420,0,420,"Z [mm]");
 addHisto("Process_nPhotTiming",3,0.5,3.5,"process"); 
 addHisto("Time_nPhotTiming",100,0,50,"time [ns]"); 

 addHisto("time_frac50",300,10,20,"time[ns]");
 addHisto("time_frac50_shaped_apd",500,16,19,"time [ns]");
 addHisto("time_frac50_shaped",500,16,19,"time [ns]");

 addHisto2D("time_frac50_vs_Z",420,0,420,"Z [mm]",500,10,100,"time [ns]");
 addHisto2D("time_frac50_vs_Z_fast",420,0,420,"Z [mm]",500,10,30,"time [ns]");

 addHisto("EactLYScaled",100,37000,47000,"E [GeV]");

}


void Analyzer::Loop(std::string setup, std::string energy)
{
  //  gROOT->ProcessLine("gPrintViaErrorHandler = kTRUE");
  //  gROOT->ProcessLine("gErrorIgnoreLevel = 6000");

//   In a ROOT session, you can do:
//      Root > .L Analyzer.C
//      Root > Analyzer t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   std::vector<float> LY =createLY(12);
   TRandom3 gausGen;
   for(int i=0;i<LY.size();++i){
     LY[i]*= gausGen.Gaus(1,0.04);
     std::cout<<gausGen.Gaus(1,0.04);
   }

         std::sort(LY.begin(),LY.end());
            std::reverse(LY.begin(),LY.end());

   //      std::random_shuffle(LY.begin(),LY.end());
   //      std::random_shuffle(LY.begin(),LY.end());

   Long64_t nentries = fChain->GetEntries();
   if(LIMIT_NENTRIES)   nentries = 200;
   energy_=energy;
   setup_=setup;

   outDir_="outFiles/plots/"+energy+"/"+setup+"/";
   if(LIMIT_NENTRIES)outDir_+="test/";
   system(Form("mkdir -p %s", outDir_.Data()));
   TString fileName = "outFiles/plotterTiming_"+setup+"_"+energy;
   if(LIMIT_NENTRIES)fileName+="_test_";
   fileName+=".root";



  createHistos();

  int nPhotTiming=200;
  LoadTree(0);
  fChain->GetEntry(0);
  //  nPhotTiming=(int)(Time_deposit->size()*0.01);
  std::cout<<"averaging time on "<<nPhotTiming<<" photons"<<std::endl;


   Long64_t nbytes = 0, nb = 0;
   TH1F* wave = new TH1F("wave","wave",1024,0,200);


   //photodetector response

   ///MAPD response from ECAL fit of APD pulse shape
   RooRealVar t("t","t", 0, 200);
   TF1 *f = new  TF1("f", "pow(exp(1.)*[1]*(x)/[0],[0])*exp(-[1]*(x))",0.,200.);
   f->SetParameters(3.0, 0.8);
	
   RooAbsPdf* rfa1 = RooFit::bindPdf(f,t);
   TH1* MAPDResponse = rfa1->createHistogram("t",1024);
   TGraph* MAPDGraphResp = new TGraph(MAPDResponse);
   WaveformNew* waveRespNew = new WaveformNew(MAPDResponse->GetNbinsX(),MAPDGraphResp->GetX(),MAPDGraphResp->GetY());

   //Response of preamp board, from Ulf Roser simulations of circuit
   TFile* file2 = TFile::Open("circuitResponse.root");
   TGraph* graphCircuit = (TGraph*) file2->Get("WaveFormCircuit");
   WaveformNew* waveCircuit = new WaveformNew(graphCircuit->GetN(),graphCircuit->GetX(),graphCircuit->GetY());


  TFile* outFile = TFile::Open(fileName,"recreate");

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if(jentry%100 == 0)std::cout<<"Entry:"<<jentry<<"/"<<nentries<<std::endl;
      float averageTiming = 0;
      float averageZ = 0;
      float averageTiming_wls = 0;
      float averageZ_wls = 0;
      float averageTiming_scint = 0;
      float averageZ_scint = 0;
      float averageTiming_cher = 0;
      float averageZ_cher = 0;

      int nPhotTiming_wls=0;
      int nPhotTiming_scint=0;
      int nPhotTiming_cher=0;

      float EactTotal=0;

      for(int i=0;i<nLayers;++i){
	//	EactTotal += Eact_layer->at(i)*LY[i];   //FIXMEE
	//	std::cout<<Eact_layer->at(i)<<" "<<LY[i]<<std::endl;
      }
      //      std::cout<<EactTotal<<std::endl;

      histos_["EactLYScaled"]->Fill(EactTotal);


      histos_["nPhotons"]->Fill(Fibre_0);

      for(int j=0;j<Time_deposit->size();++j){

	float time = Time_deposit->at(j);
	int process = Process_deposit->at(j); 
	int Z = Z_deposit->at(j); 
	

	histos2D_["time_frac50_vs_Z"]->Fill(Z,time);
	histos2D_["time_frac50_vs_Z_fast"]->Fill(Z,time);

	histos_["timeArrival"]->Fill(time);
	histos_["ZTimeArrival"]->Fill(Z);
	histos_["processTimeArrival"]->Fill(process);
	if(jentry<10)	histos_["waveform_total"]->Fill(time);
	wave->Fill(time);
	if(process==1){
	histos_["timeArrival_wls"]->Fill(time);
	histos_["ZTimeArrival_wls"]->Fill(Z);
	if(jentry<10)	histos_["waveform_wls"]->Fill(time);
	}else if(process==2){
	histos_["timeArrival_scint"]->Fill(time);
	histos_["ZTimeArrival_scint"]->Fill(Z);
	if(jentry<10)	histos_["waveform_scint"]->Fill(time);
	}else  if(process==3){
	histos_["timeArrival_cher"]->Fill(time);
	histos_["ZTimeArrival_cher"]->Fill(Z);
	if(jentry<10)	histos_["waveform_cher"]->Fill(time);
	}
      }






      //      wave->Sumw2();
      
   
      TGraph* waveGraph = new TGraph(wave);
      WaveformNew* waveF = new WaveformNew(wave->GetNbinsX(),waveGraph->GetX(),waveGraph->GetY());
      WaveformNew::max_amplitude_informations wave_max_bare = waveF->max_amplitude(4,100,3);

      float time_frac50 = waveF->time_at_frac(0,100,0.5,wave_max_bare,3);

      histos_["time_frac50"]->Fill(time_frac50);

      

    
      if(DOSHAPING){
	//convolution with MAPD
	FFTConvolution* convProd = new FFTConvolution();
	WaveformNew* convWave = convProd->fftConvolute(waveF,waveRespNew);

	//convolution with preamp
	 FFTConvolution* convProd2 = new FFTConvolution();
	 WaveformNew* convWave2 = convProd2->fftConvolute(convWave,waveCircuit);



	 if(jentry==0){
	   TCanvas c1("c_histpdf");
	   TH1F* hist = waveF->get_histo("waveF");
	   TH1F* histConvNew = convWave->get_histo("waveFNew");
	   histConvNew->SetLineColor(kGreen+2);
	   TH1F* histConvNew2 = convWave2->get_histo("waveFNew2");
	   hist->Scale(1./hist->GetMaximum());
	   histConvNew->Scale(1./histConvNew->GetMaximum());
	   histConvNew2->Scale(1./histConvNew2->GetMaximum());

	   histConvNew2->SetLineColor(kRed);
	   hist->Draw("histL");
	   histConvNew->Draw("same");
	   histConvNew2->Draw("same");
	  c1.SaveAs("histpdf_analyzer.png");
	  c1.SaveAs("histpdf_analyzer.pdf");
	}
	

	 //	WaveformNew::max_amplitude_informations wave_max_bare_shaped = convWave2->max_amplitude(4,100,3);
	 WaveformNew::max_amplitude_informations wave_max_bare_shaped_apd = convWave->max_amplitude(4,100,3);
	 histos_["time_frac50_shaped_apd"]->Fill(convWave->time_at_frac(0,100,0.5,wave_max_bare_shaped_apd,3));


	 //	WaveformNew::max_amplitude_informations wave_max_bare_shaped = convWave2->max_amplitude(4,100,3);
	 WaveformNew::max_amplitude_informations wave_max_bare_shaped = convWave2->max_amplitude(4,100,3);
	 histos_["time_frac50_shaped"]->Fill(convWave2->time_at_frac(0,100,0.5,wave_max_bare_shaped,3));



	 //	std::cout<<" frac50_shaped:"<<convWave->time_at_frac(0,100,0.5,wave_max_bare_shaped,3)<<std::endl;

      }

      
      if(jentry==200){
	TCanvas dummy;
	TF1 f("f","pol2",0,50);
	wave->Fit("f","R","",wave->GetBinCenter(wave->GetMaximumBin())-2*wave->GetBinWidth(0),wave->GetBinCenter(wave->GetMaximumBin())+2*wave->GetBinWidth(0));
	wave->Draw();
	f.SetLineColor(kRed);
	f.Draw("same");
	dummy.SaveAs("dummy.png");
      }
	  wave->Reset();	


      if(setup_=="SingleFibre" || setup_=="Ideal2016")         nPhotTiming=(int)(Time_deposit->size()*0.01);
      else nPhotTiming=(int)(Time_deposit->size()*0.01);
      //      std::cout<<nPhotTiming<<std::endl;

      //      for (auto i: Analyzer::sort_indexes( *Time_deposit)) {
	//	std::cout << Time_deposit->at(i) << std::endl;
      //      } 


      int nPhot=0;      
      for (auto i: Analyzer::sort_indexes( *Time_deposit)) {
	nPhot++;
	if(nPhot>nPhotTiming)break;
	float time = Time_deposit->at(i);
	int process = Process_deposit->at(i); 
	int Z = Z_deposit->at(i); 
	//	std::cout<<nPhot<<" "<< time << std::endl;
	if(nPhot<Time_deposit->size()){
	  averageTiming+=time;
	  histos_["processTimeArrival_avg"]->Fill(process);
	  averageZ+=Z;
	  if(nPhot<4)	  histos_["Z_firstPhotons"]->Fill(Z);
	  histos_["Z_nPhotTiming"]->Fill(Z);
	  histos_["Process_nPhotTiming"]->Fill(process);
	  histos_["Time_nPhotTiming"]->Fill(time);
	  if(process==1){
	    nPhotTiming_wls++;
	    averageTiming_wls+=time;
	    averageZ_wls+=Z;
	  }else if(process==2){
	    nPhotTiming_scint++;
	    averageTiming_scint+=time;
	    averageZ_scint+=Z;
	  }else if(process==3){
	    nPhotTiming_cher++;
	    averageTiming_cher+=time;
	    averageZ_cher+=Z;
	  }
	}
      }
      histos_["timeArrival_avg"]->Fill(averageTiming/nPhotTiming);
      histos_["ZTimeArrival_avg"]->Fill(averageZ/nPhotTiming);
      
      histos_["timeArrival_avg_wls"]->Fill(averageTiming_wls/nPhotTiming_wls);
      histos_["ZTimeArrival_avg_wls"]->Fill(averageZ_wls/nPhotTiming_wls);
      
      histos_["timeArrival_avg_scint"]->Fill(averageTiming_scint/nPhotTiming_scint);
      histos_["ZTimeArrival_avg_scint"]->Fill(averageZ_scint/nPhotTiming_scint);
      
      histos_["timeArrival_avg_cher"]->Fill(averageTiming_cher/nPhotTiming_cher);
      histos_["ZTimeArrival_avg_cher"]->Fill(averageZ_cher/nPhotTiming_cher);
   }


   histos_["Z_nPhotTimingNorm"] = (TH1F*)histos_["Z_nPhotTiming"]->Clone();
   histos_["ZTimeArrival"]->Print();
   TH1F* dummy=(TH1F*)histos_["ZTimeArrival"]->Clone();
   histos_["Z_nPhotTimingNorm"] -> Divide(dummy);

   histos_["Z_firstPhotonsNorm"] = (TH1F*)histos_["Z_firstPhotons"]->Clone();
   histos_["Z_firstPhotonsNorm"] -> Divide(dummy);



   fitHisto(histos_["timeArrival_avg"],resValueTime,resValueTimeErr);
   fitHisto(histos_["time_frac50"],resValueTime_frac50,resValueTimeErr_frac50);
   fitHisto(histos_["nPhotons"],resValueEnergy,resValueEnergyErr,true,meanValueEnergy,meanValueEnergyErr);
   fitHisto(histos_["EactLYScaled"],resValueEactLYScaled,resValueEactLYScaledErr,true,meanValueEactLYScaled,meanValueEactLYScaledErr);
   if(DOSHAPING){
     fitHisto(histos_["time_frac50_shaped_apd"],resValueTime_frac50_shaped_apd,resValueTimeErr_frac50_shaped_apd);
     fitHisto(histos_["time_frac50_shaped"],resValueTime_frac50_shaped,resValueTimeErr_frac50_shaped);

   }
   drawHistos();

   TCanvas wave_canvas;
   histos_["waveform_total"]->Draw();
   histos_["waveform_scint"]->SetLineColor(kRed);
   histos_["waveform_scint"]->Draw("same");
   histos_["waveform_cher"]->SetLineColor(kBlue);
   histos_["waveform_wls"]->SetLineColor(kViolet);
   histos_["waveform_cher"]->Draw("same");
   histos_["waveform_wls"]->Draw("same");
   wave_canvas.SaveAs(outDir_+"waveform.png");
   wave_canvas.SaveAs(outDir_+"waveform.pdf");

   meanValueEnergy->Write("meanValueEnergy");
   meanValueEnergyErr->Write("meanValueEnergyErr");

   resValueEnergy->Write("resValueEnergy");
   resValueEnergyErr->Write("resValueEnergyErr");

   meanValueEactLYScaled->Write("meanValueEactLYScaled");
   meanValueEactLYScaledErr->Write("meanValueEactLYScaledErr");

   resValueEactLYScaled->Write("resValueEactLYScaled");
   resValueEactLYScaledErr->Write("resValueEactLYScaledErr");


   resValueTime->Write("resValueTime");
   resValueTimeErr->Write("resValueTimeErr");

   resValueTime_frac50->Write("resValueTime_frac50");
   resValueTimeErr_frac50->Write("resValueTimeErr_frac50");

   resValueTime_frac50_shaped_apd->Write("resValueTime_frac50_shaped_apd");
   resValueTimeErr_frac50_shaped_apd->Write("resValueTimeErr_frac50_shaped_apd");
   
   resValueTime_frac50_shaped->Write("resValueTime_frac50_shaped");
   resValueTimeErr_frac50_shaped->Write("resValueTimeErr_frac50_shaped");


   outFile->Write();
   outFile->Close();

}


void Analyzer::fitHisto(TH1F* histo, TVectorD* res, TVectorD* resErr,bool isEnergy,TVectorD* mean, TVectorD* meanErr){
//-----------------fit with cruijff ------------------------
  double peakpos = histo->GetBinCenter(histo->GetMaximumBin());
  double sigma = histo->GetRMS();

  double fitmin;
  double fitmax;
  
  
  fitmin = peakpos-4*sigma;
  fitmax = peakpos+4*sigma;
  
  RooRealVar x("x","deltaT", fitmin, fitmax);
  RooDataHist data("data","dataset with x",x,RooFit::Import(*histo) );

  RooRealVar meanr("meanr","Mean",peakpos+sigma,peakpos-3*sigma, peakpos+3*sigma);
  RooRealVar widthL("widthL","#sigmaL",sigma , 0, 5*sigma);
  RooRealVar widthR("widthR","#sigmaR",sigma , 0, 5*sigma);
  RooRealVar alphaL("alphaL","#alpha",5.08615e-02 , 0., 1.);
  RooRealVar alphaR("alphaR","#alpha",5.08615e-02, 0., 1.);
  int ndf;

  TGaxis::SetMaxDigits(2);
  RooPlot* frame;

  RooCruijff fit_fct("fit_fct","fit_fct",x,meanr,widthL,widthR,alphaL,alphaR); ndf = 5;
  fit_fct.fitTo(data);
  
  frame = x.frame("Title");
  if(!isEnergy)  frame->SetXTitle("time [ns]");
  else if(setup_=="Ideal2016")frame->SetXTitle("Energy [MeV]");
  else frame->SetXTitle("nPhotons");

  std::string ytitle = Form("Events");
  frame->SetYTitle(ytitle.c_str());
  
  data.plotOn(frame,RooFit::Precision(1));  //this will show histogram data points on canvas 
  fit_fct.plotOn(frame);//this will show fit overlay on canvas  

  double rms,rmsErr;
  rms = (widthL.getVal()+widthR.getVal())/2;
  rmsErr = 0.5*sqrt(widthL.getError()*widthL.getError()+widthR.getError()*widthR.getError());


  TCanvas* cans = new TCanvas();
  cans->cd();
  frame->Draw();
  TLegend* lego = new TLegend(0.57, 0.75, 0.89, 0.9); 
  lego->SetTextAlign(32);
  lego->SetTextSize(0.036);
  if(!isEnergy){
    lego->AddEntry(  (TObject*)0 ,Form("#mu = %.0f #pm %.0f ps", meanr.getVal()*1.e3, meanr.getError()*1.e3), "");
    lego->AddEntry(  (TObject*)0 ,Form("#sigma = %.0f #pm %.0f ps", rms*1.e3, rmsErr*1.e3), "");

    res[0]=rms*1.e3;
    resErr[0]=rmsErr*1.e3;

  }else{
    lego->AddEntry(  (TObject*)0 ,Form("#mu = %.0f #pm %.0f ps", meanr.getVal(), meanr.getError()), "");
    lego->AddEntry(  (TObject*)0 ,Form("#sigma = %.0f #pm %.0f ps", rms, rmsErr), "");

    mean[0]=meanr.getVal();
    meanErr[0]=meanr.getError();
    
    res[0]=rms;
    resErr[0]=rmsErr;


  }



  lego->SetFillColor(0);
  lego->Draw("same");
  
  //TPaveText* pave = DrawTools::getLabelTop_expOnXaxis(energy_+" Electron Beam");
  //pave->Draw("same");

  cans->SaveAs(outDir_+histo->GetName()+"_cruijff_fit.png");
  cans->SaveAs(outDir_+histo->GetName()+"_cruijff_fit.pdf");

  //  TGaxis::SetMaxDigits(5);//default value

}



std::vector<float> Analyzer::createLY(int nLayers){
  ///FIXME for the moment just quick initialization since in a hurry make it configurable

  std::vector<float> LYvalues;

//  LYvalues.push_back(1.052);
//  LYvalues.push_back(1.096);
//  LYvalues.push_back(0.973);
//  LYvalues.push_back(1.082);
//  LYvalues.push_back(1.044);
//  LYvalues.push_back(0.977);
//  LYvalues.push_back(1.024);
//  LYvalues.push_back(0.992);
//  LYvalues.push_back(1.037);
//  LYvalues.push_back(1.007);



  LYvalues.push_back(1.029);
  LYvalues.push_back(0.996);
  LYvalues.push_back(0.989);
  LYvalues.push_back(1.040);
  LYvalues.push_back(1.011);
  LYvalues.push_back(1.000);
  LYvalues.push_back(1.016);
  LYvalues.push_back(1.028);
  LYvalues.push_back(0.979);
  LYvalues.push_back(1.021);
  LYvalues.push_back(1.001);
  LYvalues.push_back(1.018);

  return LYvalues;

}
