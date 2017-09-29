#include "interface/FFTConvolution.h"

#include "TGraph.h"
#include "TF1.h"
#include "TMath.h"
#include <assert.h>

#include "TComplex.h"

using namespace std;


////////////////////////////////////////////////////////////////////////////////
 /// 'Extend' means is that the input p.d.f convolution observable range is widened to include the buffer range
 /// 'Flat' means that the buffer is filled with the p.d.f. value at the boundary of the observable range
 /// 'Mirror' means that the buffer is filled with a ,irror image of the p.d.f. around the convolution observable boundary 
 ///
 /// The default strategy is extend. If one of the input p.d.f.s is a RooAddPdf, it is configured so that the interpretation
 /// range of the fraction coefficients is kept at the nominal convolutions observable range (instead of interpreting coefficients
 /// in the widened range including the buffer)
void FFTConvolution::setBufferStrategy(BufStrat bs) 
 {
  _bufStrat = bs ;
 }


/// Change the size of the buffer on either side of the observable range to frac times the
/// size of the range of the convolution observable
 
 void FFTConvolution::setBufferFraction(Double_t frac)
{
   if (frac<0) {
     std::cout<< "FFTConvolution::setBufferFraction(" << frac << ") fraction should be greater than or equal to zero" << endl ;
     return ;
   }
   _bufFrac = frac ;
}

//convolution of two waveforms using fft  
WaveformNew* FFTConvolution::fftConvolute(WaveformNew*  wave1,WaveformNew*  wave2,bool con){


  int N=wave1->_samples.size();//n samples of original waveform
  //assert(N == wave2->_samples.size());

  int NBuf=static_cast<int>((N*bufferFraction())/2+0.5);//number of buffer samples needed for fft
  int N2=N+2*NBuf;//number of samples + buffer zones

  int zeroBin1,zeroBin2;

  Double_t* input1 = sampleWaveform(wave1,N,NBuf,N2,zeroBin1);
  Double_t* input2 = sampleWaveform(wave2,N,NBuf,N2,zeroBin2);

  //TVirtualFFT* fftr2c1 = TVirtualFFT::FFT(1, &N2, "R2CK");
  //TVirtualFFT* fftr2c2 = TVirtualFFT::FFT(1, &N2, "R2CK");
  //TVirtualFFT* fftc2r = TVirtualFFT::FFT(1, &N2, "C2RK");
  
  TVirtualFFT* fftr2c1 = TVirtualFFT::FFT(1, &N2, "R2C EX K");
  TVirtualFFT* fftr2c2 = TVirtualFFT::FFT(1, &N2, "R2C EX K");
  TVirtualFFT* fftc2r = TVirtualFFT::FFT(1, &N2, "C2R EX K");

  // Real->Complex FFT Transform on p.d.f. 1 sampling
  fftr2c1->SetPoints(input1);
  fftr2c1->Transform();
  
  // Real->Complex FFT Transform on p.d.f 2 sampling
  fftr2c2->SetPoints(input2);
  fftr2c2->Transform();
  
  // Loop over first half +1 of complex output results, multiply 
  // and set as input of reverse transform
  for (Int_t i=0 ; i<N2 ; i++) {
  //for (Int_t i=0 ; i<N2/2+1 ; i++) {
    Double_t re1,re2,im1,im2 ;
    fftr2c1->GetPointComplex(i,re1,im1) ;
    fftr2c2->GetPointComplex(i,re2,im2) ;
    //cout << re1 << "    " << im1 << "    " << re2 << "    " << im2 << endl;
    Double_t re=0,im=0;
    //if (i<5)
    {
      if (con){
        re = re1*re2 - im1*im2 ;
        im = re1*im2 + re2*im1 ;
      }
      else{
        re = (re1*re2 + im1*im2)/(re2*re2+im2*im2);
        im = (-re1*im2 + re2*im1)/(re2*re2+im2*im2);
      }
    }
    TComplex t(re,im) ;
    //cout << re << "        " << im << endl;
    fftc2r->SetPointComplex(i,t) ;
  }
  
  // Reverse Complex->Real FFT transform product
  fftc2r->Transform() ;

  Int_t totalShift = zeroBin1 + (N2-N)/2 ;

  WaveformNew* convolutedWave = new WaveformNew();
  
  for (Int_t i =0 ; i<N ; i++) {
    
    // Cyclically shift array back so that bin containing zero is back in zeroBin
    Int_t j = i + totalShift ;
    while (j<0) j+= N2 ;
    while (j>=N2) j-= N2 ;
    convolutedWave->addSample(fftc2r->GetPointReal(j),(wave1->_times[1]-wave1->_times[0]));//FIXME do dynamic
  }


  // Delete input arrays
   delete[] input1 ;
   delete[] input2 ;

   return   convolutedWave;

}

Double_t* FFTConvolution::sampleWaveform(WaveformNew*  wave, int N,int NBuf,int N2, int& zeroBin){

  // Allocate array of sampling size plus optional buffer zones
  Double_t* array = new Double_t[N2] ;

  zeroBin=0;
  
  // Find bin ID that contains zero value
  for(int i=0;i<N;i++){
    if(wave->_times[i]<0)zeroBin++; 
    else break;
  }
  

  //  int binShift = int(N2/(wave->_times[N-1]-wave->_times[0]));
  //  zeroBin+=binShift;
  zeroBin+=NBuf-1;
  
  while(zeroBin>=N2) zeroBin-= N2 ;
  while(zeroBin<0) zeroBin+= N2 ;

  // First scan hist into temp array 
  Double_t *tmp = new Double_t[N2] ;
  Int_t k(0) ;

	
  switch(_bufStrat) {
    
  case Extend:
    // Sample entire extended range (N2 samples)
    for (k=0 ; k<N2 ; k++) {
      if (k<N)      tmp[k]= wave->_samples[k];
      else tmp[k]=0;
    }  
    break;

  case Flat:    
    // Sample original range (N samples) and fill lower and upper buffer
    // bins with p.d.f. value at respective boundary
    for (k=0 ; k<NBuf ; k++) {
      tmp[k] = wave->_samples[0];
    }
    for(k=0 ; k<N ; k++) {
      tmp[k+NBuf] = wave->_samples[k];
    }
    for (k=0;k<NBuf;k++){
      tmp[N+NBuf+k] = wave->_samples[N-1];
    }

  case Mirror:
    // Sample original range (N samples) and fill lower and upper buffer
    // bins with mirror image of sampled range
    for (k=0 ; k<N ; k++) {
      tmp[k+NBuf] = wave->_samples[k];
    }  
    for (k=1 ; k<=NBuf ; k++) {
      tmp[NBuf-k] = wave->_samples[k];
      tmp[NBuf+N+k-1] = wave->_samples[N-k];
    }  
    break ;
    
  }

  // Scan function and store values in array
  for (Int_t i=0 ; i<N2 ; i++) {
    // Cyclically shift writing location by zero bin position    
    Int_t j = i - (zeroBin) ;
    if (j<0) j+= N2 ;
    if (j>=N2) j-= N2 ;
    array[i] = tmp[j] ;
  }  
  
  // Cleanup 
  delete[] tmp ;
  return array ;


}
