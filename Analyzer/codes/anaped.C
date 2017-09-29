#include <TVirtualFFT.h>
#include <math.h>
#include <dirent.h>
#include <TFile.h>
#include <TGraph.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TChain.h>
#include <TCanvas.h>
#include "Event.h"
#define NGRP 4
#define NCH 36
//#define NGSPS 5.0 // 5.0 Gs/s
#define ADCSCALE 0.25e-3    // 250 uV/adc count
#define NSAMPLE 1024
#define NFILEMAX 2000    // Max number of files to read
#define UNDER_SAMPLING 1
#define IS0 0
#include "ped_25NGSPS_ch20.dat"
double NGSPS;
double shape[36][1024],dshape[36][1024],dnorm[36];
int cur_ch;

double fshape(double *x, double *par)
{
  double t0=par[0];
  double qmax=par[1];
  double t=x[0]-t0;
  int i=(int)(t*NGSPS);
  double y=0.;
  if(i>0 && i<NSAMPLE-50)
  {
    double t1=((double)i)/NGSPS;
    double t2=((double)i+1)/NGSPS;
    double y1=shape[cur_ch][i];
    double y2=shape[cur_ch][i+1];
    y=y1+(y2-y1)/(t2-t1)*(t-t1);
  } 
  return y*qmax;
}

void anaped (int numrun=5297, int BW4=100, double ngsps=5.0, int nfile_max=NFILEMAX) 
{
  //NGSPS=5.0; // 5.0 Gsps by default. Should be overwritten at first event since the actual value is in the data
  //NGSPS=2.5; // 5.0 Gsps by default. Should be overwritten at first event since the actual value is in the data
  int ch_ref=20;
  NGSPS=ngsps;
  int n_LPF = 4;
  n_LPF = n_LPF/2;
  double s = NGSPS*1.e9;
  double f = 100.e6;
  double a = tan(M_PI*f/s);
  double aa = a*a;
  double r;
  double *A = (double *)malloc(n_LPF*sizeof(double));
  double *d1 = (double *)malloc(n_LPF*sizeof(double));
  double *d2 = (double *)malloc(n_LPF*sizeof(double));
  double *w0 = (double *)calloc(n_LPF, sizeof(double));
  double *w1 = (double *)calloc(n_LPF, sizeof(double));
  double *w2 = (double *)calloc(n_LPF, sizeof(double));
  double x;
  for(int i=0; i<36; i++)
  {
    char fname[80];
    int ch=i; 
    if(i!=18 && i!=24 && i!=25 && i!=20)ch=24;
    sprintf(fname,"shape_%d.dat",ch);
    FILE *fd=fopen(fname,"r");
    double qmax=0.;
    for(int is=0; is<1024; is++)
    {
      fscanf(fd,"%lf",&shape[i][is]);
      if(shape[i][is]>qmax)qmax=shape[i][is];
    }
    for(int is=0; is<1024; is++)shape[i][is]/=qmax;
    for(int is=0; is<1023; is++)dshape[i][is]=(shape[i][is+1]-shape[i][is])/(1./NGSPS);
    dshape[i][1023]=0.;
    dnorm[i]=0.;
    for(int is=0; is<1024; is++){if(dshape[i][is]>dnorm[i])dnorm[i]=dshape[i][is];};
    printf("ch %d : dnorm : %e\n",i,dnorm[i]);
  }

  for(int i=0; i<n_LPF; i++)
  {
    r = sin(M_PI*(2.0*i+1.0)/(4.0*n_LPF));
    s = aa + 2.0*a*r + 1.0;
    A[i] = aa/s;
    d1[i] = 2.0*(1-aa)/s;
    d2[i] = -(aa - 2.0*a*r + 1.0)/s;
  }

  double pi=2.*asin(1.);
  double dt=1./(NGSPS*1.e9);
// Butterworth filter order 4
  double o1=0.7654, o2=1.8478;
// Cut-off frequency :
  double oc=2.*pi*BW4*1.e6;
  double a1=1+o1*oc*dt/2.+(oc*dt/2.)*(oc*dt/2.);
  double b1=2.*((oc*dt/2.)*(oc*dt/2.)-1.);
  double c1=1-o1*oc*dt/2.+(oc*dt/2.)*(oc*dt/2.);
  double a2=1+o2*oc*dt/2.+(oc*dt/2.)*(oc*dt/2.);
  double b2=2.*((oc*dt/2.)*(oc*dt/2.)-1.);
  double c2=1-o2*oc*dt/2.+(oc*dt/2.)*(oc*dt/2.);
  double ad=a1*a2;
  double bd=a1*b2+b1*a2;
  double cd=a1*c2+b1*b2+c1*a2;
  double dd=b1*c2+c1*b2;
  double ed=c1*c2;
  double an=1.;
  double bn=4.;
  double cn=6.;
  double dn=4.;
  double en=1.;
  double g=pow(dt*oc,4)/16.;
  an=g*an/ad;
  bn=g*bn/ad;
  cn=g*cn/ad;
  dn=g*dn/ad;
  en=g*en/ad;
  bd=bd/ad;
  cd=cd/ad;
  dd=dd/ad;
  ed=ed/ad;
  ad=1.;
  double y4=0.,y3=0.,y2=0.,y1=0.;
  double x4=0.,x3=0.,x2=0.,x1=0.;
  printf("poly 1 : %e %e %e\n",a1,b1,c1);
  printf("poly 2 : %e %e %e\n",a2,b2,c2);
  printf("poly d : %e %e %e %e %e\n",ad,bd,cd,dd,ed);
  printf("poly n : %e %e %e %e %e\n",an,bn,cn,dn,en);
  printf("gain %e\n",g);

  double yref=1.;
// Step response :
  for(int i=0; i<200; i++)
  {
    double value=0.;
    if(i>10)value=1.;
    double y=x4*en+x3*dn+x2*cn*x1*bn+value*an-y4*ed-y3*dd-y2*cd-y1*bd;
    x4=x3; x3=x2; x2=x1; x1=value;
    y4=y3; y3=y2; y2=y1; y1=y;
    yref=y;
  }
  y4=0.;y3=0.;y2=0.;y1=0.;
  x4=0.;x3=0.;x2=0.;x1=0.;
  printf("Butterworth total gain : %f\n",yref);

  int nfile=0;
  TF1 *f1=new TF1("shape",fshape,0.,1024./NGSPS,2);

// Scan directory for run numrun and register files to analyze :
  char flist[NFILEMAX][500];
  char gname[500],dname[500];
  sprintf(dname,"%d.list",numrun);
  FILE *fl=fopen(dname,"r");
  if(fl==NULL)
  {
    sprintf(dname,"data/digi/%d/",numrun);
    struct dirent *dirlist;
    DIR *dird = opendir(dname);
    if( dird == NULL)
    {
      printf(" error while opening rawdata directory : %s\n",dname);
      exit(-1);
    }
    int iret=1;
    while(iret > 0)  
    {    
      dirlist = readdir( dird );  
      if (dirlist !=NULL)  
      {    
        char *cpos1=strstr(dirlist->d_name,".root");
        if(cpos1!=NULL)
        {
          //printf("%s\n", dirlist->d_name);
          sprintf(flist[nfile++],"%s%s",dname,dirlist->d_name);
        }
      }
      else
        iret = 0;
    }
    iret = closedir( dird );
    printf(" Closing directory scan : %d, %d valid files found\n", iret,nfile);
  }
  else
  {
    int eof=0;
    char ctmp[500];
    while(eof!=EOF)
    {
      eof=fscanf(fl,"%s",ctmp);
      if(eof==EOF)continue;
      //printf("%d %s\n",nfile,ctmp);
      if(strstr(ctmp,"eoscms")==NULL)
        sprintf(flist[nfile++],"root://eoscms.cern.ch//eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Jul2016/raw/DataTree/%d/%s",
              numrun,ctmp);
      else
        sprintf(flist[nfile++],"%s",ctmp);
    }
    fclose(fl);
  }
  if(nfile>nfile_max)nfile=nfile_max;

  printf("Start with file %s\n",flist[0]);
  TChain *inputTree_ = new TChain ("H4tree","H4tree") ;
  for(int ifile=0; ifile<nfile; ifile++)
  {
    printf("Add file %d : %s\n",ifile,flist[ifile]);
    inputTree_->Add(flist[ifile]) ;
  }

  treeStructData treeStruct_;

  inputTree_->SetBranchAddress("evtNumber"    ,&treeStruct_.evtNumber);
  inputTree_->SetBranchAddress("evtTimeDist"  ,&treeStruct_.evtTimeDist);
  inputTree_->SetBranchAddress("evtTimeStart" ,&treeStruct_.evtTimeStart);

  inputTree_->SetBranchAddress("nEvtTimes"    ,&treeStruct_.nEvtTimes);
  inputTree_->SetBranchAddress("evtTime"      ,&treeStruct_.evtTime);
  inputTree_->SetBranchAddress("evtTimeBoard" ,&treeStruct_.evtTimeBoard);

  inputTree_->SetBranchAddress("nAdcChannels" ,&treeStruct_.nAdcChannels);
  inputTree_->SetBranchAddress("adcBoard"     ,treeStruct_.adcBoard);
  inputTree_->SetBranchAddress("adcChannel"   ,treeStruct_.adcChannel);
  inputTree_->SetBranchAddress("adcData"      ,treeStruct_.adcData);

  inputTree_->SetBranchAddress("nTdcChannels" ,&treeStruct_.nTdcChannels);
  inputTree_->SetBranchAddress("tdcBoard"     ,treeStruct_.tdcBoard);
  inputTree_->SetBranchAddress("tdcChannel"   ,treeStruct_.tdcChannel);
  inputTree_->SetBranchAddress("tdcData"      ,treeStruct_.tdcData);

  inputTree_->SetBranchAddress("nDigiSamples" ,&treeStruct_.nDigiSamples);
  inputTree_->SetBranchAddress("digiGroup"    ,treeStruct_.digiGroup);
  inputTree_->SetBranchAddress("digiChannel"  ,treeStruct_.digiChannel);
  inputTree_->SetBranchAddress("digiFrequency",treeStruct_.digiFrequency);
  inputTree_->SetBranchAddress("digiSampleIndex",treeStruct_.digiSampleIndex);
  inputTree_->SetBranchAddress("digiSampleValue",treeStruct_.digiSampleValue);
  inputTree_->SetBranchAddress("digiBoard"    ,treeStruct_.digiBoard);

  inputTree_->SetBranchAddress("nScalerWords" ,&treeStruct_.nScalerWords);
  inputTree_->SetBranchAddress("scalerWord"   ,treeStruct_.scalerWord);
  inputTree_->SetBranchAddress("scalerBoard"  ,treeStruct_.scalerBoard);

  inputTree_->SetBranchAddress("nPatterns"    ,&treeStruct_.nPatterns);
  inputTree_->SetBranchAddress("pattern"      ,treeStruct_.pattern);
  inputTree_->SetBranchAddress("patternBoard" ,treeStruct_.patternBoard);
  inputTree_->SetBranchAddress("patternChannel",treeStruct_.patternChannel);

  inputTree_->SetBranchAddress("nTriggerWords",&treeStruct_.nTriggerWords);
  inputTree_->SetBranchAddress("triggerWords" ,treeStruct_.triggerWords);
  inputTree_->SetBranchAddress("triggerWordsBoard",treeStruct_.triggerWordsBoard);

  int nEntries = inputTree_->GetEntries();
  printf("Will process %d events\n",nEntries);

  TProfile *pped=new TProfile("ped_ave","ped_ave",NSAMPLE,0.,NSAMPLE);
  TProfile *pmod=new TProfile("mod_ave","mod_ave",NSAMPLE,0.,NGSPS);
  TProfile *pmod_LPF=new TProfile("mod_ave_LPF","mod_ave_LPF",NSAMPLE,0.,NGSPS);
  TProfile *pmod_BW4=new TProfile("mod_ave_BW4","mod_ave_BW4",NSAMPLE,0.,NGSPS);
  TH1D *hmod=new TH1D("NDS","NDS",10*(NSAMPLE/2-1),0.,NGSPS/2.*1000.);
  hmod->SetLineColor(kBlue);
  hmod->SetLineWidth(2);
  TH1D *hmod_LPF=new TH1D("NDS_LPF","NDS_LPF",10*(NSAMPLE/2-1),0.,NGSPS/2.*1000.);
  hmod_LPF->SetLineColor(kGreen);
  hmod_LPF->SetLineWidth(2);
  TH1D *hmod_BW4=new TH1D("NDS_BW4","NDS_BW4",10*(NSAMPLE/2-1),0.,NGSPS/2.*1000.);
  hmod_BW4->SetLineColor(kRed);
  hmod_BW4->SetLineWidth(2);

  int size=NSAMPLE;
  TVirtualFFT *fft_f = TVirtualFFT::FFT(1, &size, "C2CF M K");

  //TGraph *tgpulse[NCH];
  char hname[80];
  double cont[NCH][NSAMPLE],cor[NCH][NSAMPLE];
  double cont_LPF[NCH][NSAMPLE];
  double cont_BW4[NCH][NSAMPLE];
  sprintf(hname,"data/graph/ped_5297.root");
  if(numrun>=5740) sprintf(hname,"data/graph/ped_5744.root");
  TFile *tfcor = new TFile(hname);
  TH1D *hcor;
  for(int ich=0; ich<NCH; ich++)
  {
    sprintf(hname,"Digitizer_ch%d_raw_mean",ich);
    gDirectory->GetObject(hname,hcor);
    for(int is=0; is<NSAMPLE; is++)
    {
      if((ich%9)!=8)
        cor[ich][is]=hcor->GetBinContent(is+1);
      else
        cor[ich][is]=0.;
    }
  }
  tfcor->Close();

  double ped;
  for(int ich=0; ich<NCH; ich++)
  {
    ped=0.;
    //for(int is=0; is<10/NGSPS; ped+=cor[ich][is++]);
    //ped/=(10./NGSPS);
    for(int is=0; is<NSAMPLE-20; ped+=cor[ich][is++]);
    ped/=NSAMPLE-20;
    for(int is=0; is<NSAMPLE; cor[ich][is++]-=ped);
    //for(int is=0; is<NSAMPLE; cor[ich][is++]=0.);
  }
  if(NGSPS==2.5){for(int is=0; is<NSAMPLE; is++){cor[20][is]=ped_25_20[is];}}

  char fout[80];
  sprintf(fout,"data/graph/%d.%4.4d_%3.3dMHz.root",numrun,nfile,BW4);
  TFile *tfo=new TFile(fout,"recreate");

// First event to :
  int freq=-1;
  ped=0.;
  inputTree_->GetEntry(0);
  //printf("Number of samples : %d\n",treeStruct_.nDigiSamples);
  for(unsigned int is=0; is<treeStruct_.nDigiSamples; is++)
  {
    int igrp=treeStruct_.digiGroup[is];
    //if(igrp>0)continue;
    int ich=igrp*9+treeStruct_.digiChannel[is];
    int idx=treeStruct_.digiSampleIndex[is];
    if(freq==-1)
    {
      freq=treeStruct_.digiFrequency[is];
      if(freq==0)     NGSPS=5.0;
      else if(freq==1)NGSPS=2.5;
      else if(freq==2)NGSPS=1.;
      else printf("Strange Frequency : %d . Should not be used !\n",freq);
      printf("New Frequency : %d : %f\n",freq,NGSPS);
    }
    double val=treeStruct_.digiSampleValue[is];
    cont[ich][idx]=val-cor[ich][is];
  }
  ped=0.;
  for(int is=0; is<100; ped+=cont[26][is++]);
  ped/=100.;

  dt=1./NGSPS;
  double tmax=dt*NSAMPLE;
  double fmax=1./dt;
  double df=fmax/NSAMPLE;

// Now, we now know the sampling frequency. We can book the histo :
  //for(int ich=0; ich<NCH; ich++)
  //{
  //  sprintf(hname,"Digitizer_ch%d_event",ich);
  //  tgpulse[ich]=new TGraph();
  //  tgpulse[ich]->SetName(hname);
  //  tgpulse[ich]->SetTitle(hname);
  //  tgpulse[ich]->SetMarkerStyle(21);
  //  tgpulse[ich]->SetMarkerSize(0.5);
  //  tgpulse[ich]->SetMarkerColor(kRed);
  //  sprintf(hname,"Digitizer_ch%d_event_under_sampled",ich);
  //}

// Second pass : analyze events and apply corrections :
  double  noise_t[NCH], ave_t[NCH];
  double  noise_t_LPF[NCH], ave_t_LPF[NCH];
  double  noise_t_BW4[NCH], ave_t_BW4[NCH];
  for(int ich=0; ich<NCH; ich++)
  {
    noise_t[ich]=0.;
    noise_t_LPF[ich]=0.;
    noise_t_BW4[ich]=0.;
    ave_t[ich]=0.;
    ave_t_LPF[ich]=0.;
    ave_t_BW4[ich]=0.;
  }
  int nevt_good=0;
  for(int iEntry=0; iEntry<nEntries; iEntry++)
  //for(int iEntry=0; iEntry<1; iEntry++)
  {
    inputTree_->GetEntry(iEntry);
    //for(int ich=0; ich<NCH; ich++)
    //{
    //  tgpulse[ich]->Clear();
    //}

    for(unsigned int is=0; is<treeStruct_.nDigiSamples; is++)
    {
      int igrp=treeStruct_.digiGroup[is];
      int ich=igrp*9+treeStruct_.digiChannel[is];
      int idx=treeStruct_.digiSampleIndex[is];
      double val=treeStruct_.digiSampleValue[is];
      //if(trig_type==0)hraw_mean[ich]->Fill((double(idx)+0.5)/NGSPS,val);
      cont[ich][idx]=val-cor[ich][idx];
    }

    for(int i=0; i<NSAMPLE; i++)
    {
      pped->Fill(0.5+i,cont[ch_ref][i]);
    }

// remove pedestal (computed on first 10 ns)
    for(int ich=0; ich<NCH; ich++)
    {
// remove events with more than 5 mV spread in first 10 ns :
      ped=0.;
      //for(int is=0; is<10/NGSPS; is++) ped+=cont[ich][is];
      //ped/=(10./NGSPS);
      for(int is=0; is<NSAMPLE-20; is++) ped+=cont[ich][is];
      ped/=NSAMPLE-20;

      for(int is=0; is<NSAMPLE; is++)
      {   
        cont[ich][is]-=ped;
        cont[ich][is]*=ADCSCALE; // mV
        if(is>NSAMPLE-20)
        {
          cont[ich][is]=cont[ich][is-20];
        }
      }

//LPF :
      for(int i=0; i<n_LPF; i++)
      {
        w2[i] = 0.;
        w1[i] = 0.;
      }
      for(int is=0; is<NSAMPLE; is++)
      {   
        for(int i=0; i<n_LPF; i++)
        {
          w0[i] = d1[i]*w1[i] + d2[i]*w2[i] + cont[ich][is];
          cont_LPF[ich][is] = A[i]*(w0[i] + 2.0*w1[i] + w2[i]);
          w2[i] = w1[i];
          w1[i] = w0[i];
        }
      }
// BW4 :
      y4=0.,y3=0.,y2=0.,y1=0.;
      x4=0.,x3=0.,x2=0.,x1=0.;
      for(int is=0; is<NSAMPLE; is++)
      {
        double value=cont[ich][is];
        double y=(x4*en+x3*dn+x2*cn*x1*bn+value*an)/yref-y4*ed-y3*dd-y2*cd-y1*bd; 
        x4=x3; x3=x2; x2=x1; x1=value;
        y4=y3; y3=y2; y2=y1; y1=y;
        cont_BW4[ich][is]=y;
      }

      for(int is=0; is<NSAMPLE; is++)
      {
        //tgpulse[ich]->SetPoint(is,t,cont[ich][is]); 
        noise_t[ich]+=cont[ich][is]*cont[ich][is];
        noise_t_LPF[ich]+=cont_LPF[ich][is]*cont_LPF[ich][is];
        noise_t_BW4[ich]+=cont_BW4[ich][is]*cont_BW4[ich][is];
        ave_t[ich]+=cont[ich][is];
        ave_t_LPF[ich]+=cont_LPF[ich][is];
        ave_t_BW4[ich]+=cont_BW4[ich][is];
      }
    }

// For selected channels, do the FFT and take the average
    double rex[NSAMPLE],    imx[NSAMPLE],    rey[NSAMPLE],    imy[NSAMPLE],    mod[NSAMPLE];
    double rex_LPF[NSAMPLE],imx_LPF[NSAMPLE],rey_LPF[NSAMPLE],imy_LPF[NSAMPLE],mod_LPF[NSAMPLE];
    double rex_BW4[NSAMPLE],imx_BW4[NSAMPLE],rey_BW4[NSAMPLE],imy_BW4[NSAMPLE],mod_BW4[NSAMPLE];
    //int n=tgpulse[ch_ref]->GetN();
    for(int i=0; i<NSAMPLE; i++)
    {
      rex[i]=0.;
      rex_LPF[i]=0.;
      rex_BW4[i]=0.;
      imx[i]=0.;
      imx_LPF[i]=0.;
      imx_BW4[i]=0.;
      rex[i]=cont[ch_ref][i];
      rex_LPF[i]=cont_LPF[ch_ref][i];
      rex_BW4[i]=cont_BW4[ch_ref][i];
    }

    fft_f->SetPointsComplex(rex, imx);
    fft_f->Transform();
    fft_f->GetPointsComplex(rey, imy);
    fft_f->SetPointsComplex(rex_LPF, imx_LPF);
    fft_f->Transform();
    fft_f->GetPointsComplex(rey_LPF, imy_LPF);
    fft_f->SetPointsComplex(rex_BW4, imx_BW4);
    fft_f->Transform();
    fft_f->GetPointsComplex(rey_BW4, imy_BW4);
    for(int i=0; i<NSAMPLE; i++)
    {   
      rey[i]/=NSAMPLE;
      imy[i]/=NSAMPLE;
      rey_LPF[i]/=NSAMPLE;
      imy_LPF[i]/=NSAMPLE;
      rey_BW4[i]/=NSAMPLE;
      imy_BW4[i]/=NSAMPLE;
      mod[i]=sqrt(rey[i]*rey[i]+imy[i]*imy[i]);
      mod_LPF[i]=sqrt(rey_LPF[i]*rey_LPF[i]+imy_LPF[i]*imy_LPF[i]);
      mod_BW4[i]=sqrt(rey_BW4[i]*rey_BW4[i]+imy_BW4[i]*imy_BW4[i]);
      pmod->Fill(i*df,mod[i]);
      pmod_LPF->Fill(i*df,mod_LPF[i]);
      pmod_BW4->Fill(i*df,mod_BW4[i]);
    }

    nevt_good++;
  }

  for(int ich=0; ich<NCH; ich++)
  {
    noise_t[ich]/=NSAMPLE*nevt_good;
    noise_t_LPF[ich]/=NSAMPLE*nevt_good;
    noise_t_BW4[ich]/=NSAMPLE*nevt_good;
    ave_t[ich]/=NSAMPLE*nevt_good;
    ave_t_LPF[ich]/=NSAMPLE*nevt_good;
    ave_t_BW4[ich]/=NSAMPLE*nevt_good;
    noise_t[ich]=sqrt(noise_t[ich]-ave_t[ich]*ave_t[ich]);
    noise_t_LPF[ich]=sqrt(noise_t_LPF[ich]-ave_t_LPF[ich]*ave_t_LPF[ich]);
    noise_t_BW4[ich]=sqrt(noise_t_BW4[ich]-ave_t_BW4[ich]*ave_t_BW4[ich]);
  }
  double noise_f=0, noise_f_LPF=0., noise_f_BW4=0.;
  double loc_noise;
  loc_noise=pmod->GetBinContent(1);
  //noise_f+=2.*loc_noise*loc_noise;
  loc_noise=pmod_LPF->GetBinContent(1);
  //noise_f_LPF+=2.*loc_noise*loc_noise;
  loc_noise=pmod_BW4->GetBinContent(1);
  //noise_f_BW4+=2.*loc_noise*loc_noise;
  for(int i=1; i<=NSAMPLE/2; i++)
  {
    loc_noise=pmod->GetBinContent(i+1);
    noise_f+=loc_noise*loc_noise;
    for(int j=10*(i-1)+1; j<10*(i-1)+11; j++) hmod->SetBinContent(j,loc_noise/sqrt(df*1.e9)*1.e9); // F in Hz and V in nV
    loc_noise=pmod->GetBinContent(NSAMPLE-i+1);
    noise_f+=loc_noise*loc_noise;

    loc_noise=pmod_LPF->GetBinContent(i+1);
    noise_f_LPF+=loc_noise*loc_noise;
    for(int j=10*(i-1)+1; j<10*(i-1)+11; j++) hmod_LPF->SetBinContent(j,loc_noise/sqrt(df*1.e9)*1.e9);
    loc_noise=pmod_LPF->GetBinContent(NSAMPLE-i+1);
    noise_f_LPF+=loc_noise*loc_noise;

    loc_noise=pmod_BW4->GetBinContent(i+1);
    noise_f_BW4+=loc_noise*loc_noise;
    for(int j=10*(i-1)+1; j<10*(i-1)+11; j++) hmod_BW4->SetBinContent(j,loc_noise/sqrt(df*1.e9)*1.e9);
    loc_noise=pmod_BW4->GetBinContent(NSAMPLE-i+1);
    noise_f_BW4+=loc_noise*loc_noise;
  }
  pped->Write();
  pmod->Write();
  pmod_LPF->Write();
  pmod_BW4->Write();
  hmod->Write();
  hmod_LPF->Write();
  hmod_BW4->Write();
  tfo->Close();
  printf("End of analysis : %d good evt\n",nevt_good);
  printf("No filter noise from spectrum = %e, from samples =%e V\n",sqrt(noise_f),noise_t[ch_ref]);
  printf("LPF filter noise from spectrum = %e, from samples =%e V\n",sqrt(noise_f_LPF),noise_t_LPF[ch_ref]);
  printf("BW4 filter noise from spectrum = %e, from samples =%e V\n",sqrt(noise_f_BW4),noise_t_BW4[ch_ref]);
  for(int ich=0; ich<NCH;ich++)
  {
    printf("ich %d : noise: no_filter %e, LPF %e, BW4 %e\n",ich,noise_t[ich],noise_t_LPF[ich],noise_t_BW4[ich]);
  }
  return;
}
