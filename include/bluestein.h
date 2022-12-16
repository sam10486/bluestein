#ifndef _BLUESTEIN_H_
#define _BLUESTEIN_H_

#include <vector>
#include <iostream>
#include <NTL/ZZ.h>
#include <math.h>

#include "NTT.h"

using namespace NTL;

class bluestein{
public:
  std::vector<ZZ> h2_freq;  // stroing m-th twiddle  factor at frequency domain
  std::vector<ZZ> h2_time; 
  unsigned long n;  //original length
  unsigned long M;
  //samll prime
  ZZ prime;
  ZZ ROU;
  ZZ IROU;
  //Power of 2 prime 
  ZZ Prime_pof2;
  ZZ ROU_pof2; 
  //
  NTT FFT; 
  
  void init(ZZ P,ZZ root,unsigned long length);
  void BFFT(std::vector<ZZ> &A);
  void N_ROU(ZZ prime,unsigned long n,ZZ &root);
  int  PROU(ZZ prime,ZZ root,unsigned long n);
};

#endif