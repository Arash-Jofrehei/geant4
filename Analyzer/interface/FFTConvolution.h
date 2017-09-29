#ifndef FFTConvolution_h
#define FFTConvolution_h

#include "WaveformNew.h"


class FFTConvolution : public WaveformNew
{

 public:
  WaveformNew* fftConvolute(WaveformNew*  wave1, WaveformNew*  wave2, bool con = true);//convolution of two waveforms using fft  
  Double_t*  sampleWaveform(WaveformNew*  wave, int N,int NBuf,int N2,int& zeroBin);//performs needed operation for fft: buffer zones, time shifting....

  Double_t bufferFraction() const { 
    // Return value of buffer fraction applied in FFT calculation array beyond either
    // end of the observable domain to reduce cyclical effects
    return _bufFrac ; 
  }
  
  enum BufStrat { Extend=0, Mirror=1, Flat=2 } ;

  Double_t _bufFrac ; // Sampling buffer size as fraction of domain size 
  BufStrat _bufStrat ; // Strategy to fill the buffer


  BufStrat bufferStrategy() const {
    // Return the strategy currently used to fill the buffer: 
    // 'Extend' means is that the input p.d.f convolution observable range is widened to include the buffer range
    // 'Flat' means that the buffer is filled with the p.d.f. value at the boundary of the observable range
    // 'Mirror' means that the buffer is filled with a mirror image of the p.d.f. around the convolution observable boundary 
    return _bufStrat ;
  }
  void setBufferStrategy(BufStrat bs) ;
  void setBufferFraction(Double_t frac) ;
  

  FFTConvolution(BufStrat strategy=Extend, Double_t frac=0.1){//default values
    setBufferStrategy(strategy);
    setBufferFraction(frac);
  }

};

#endif
