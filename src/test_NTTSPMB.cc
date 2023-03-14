// Copyright (C) IBM, All Rights Reserved

#include <cstring>
#include <string>
#include <iomanip>
#include <ctime>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <NTL/tools.h>
#include <NTL/ZZ.h>
#include <time.h>
#include <bluestein.h>
#include <BitOperate.h>

using namespace std;

NTL_CLIENT
#include "NTTSPMB.h"
#include "NTT.h"
void test_NTTSPMB()
{
  NTTSPMB test;
  NTT     NTT_test;
  long difference_length;
  long difference_16;
  unsigned long fft_point;
  int band_memory_size; 
  int radix;

  //----------bluestein--------------
  bluestein blue;
  ZZ tmp_prime;
  conv(tmp_prime, "97");
  ZZ ROU;
  unsigned int u_n = 16;
  blue.N_ROU(tmp_prime, u_n, ROU);
  cout << "ROU = " << ROU << std::endl;
  //---------------------------------

  radix = 16;
  ZZ fft_prime;
  ZZ fft_twiddle;
  ZZ fft_twiddle_16;
  ZZ fft_IW;
  ZZ fft_twiddle_65536;
  
  fft_point         = 65536;//16
  difference_length = 65536 / fft_point;
  difference_16     = fft_point / 16;
  band_memory_size  = fft_point / 32;
  conv(fft_prime,"18446744069414584321");
  conv(fft_twiddle_65536,"14603442835287214144");  //65536-th twiddle factor
  //-------test--------
  //conv(fft_prime,"197");
  //conv(fft_twiddle_65536,"8");  //65536-th twiddle factor
  //-------------------
  
  PowerMod(fft_twiddle,fft_twiddle_65536,difference_length,fft_prime);
  std::cout << "difference_length = " << difference_length << ", fft_twiddle = " << fft_twiddle << std::endl;
   
  InvMod(fft_IW,fft_twiddle,fft_prime);
 
  std::cout << "test NTTSPMB Init!!! \n";
  //---
  test.init(fft_point,fft_prime,fft_twiddle,radix);
  NTT_test.NTT_init(fft_point,fft_prime,fft_twiddle);
  //----

  std::vector<ZZ> A;
  std::vector<ZZ> A_1;
  std::vector<ZZ> A_NTT_B0R0;
  std::vector<ZZ> A_NTT_B0R1;
  std::vector<ZZ> A_NTT_B0R2;
  std::vector<ZZ> A_NTT_B0R3;
  std::vector<ZZ> A_NTT_B0R4;
  std::vector<ZZ> A_NTT_B0R5;
  std::vector<ZZ> A_NTT_B0R6;
  std::vector<ZZ> A_NTT_B0R7;
  std::vector<ZZ> A_NTT_B0R8;
  std::vector<ZZ> A_NTT_B0R9;
  std::vector<ZZ> A_NTT_B0R10;
  std::vector<ZZ> A_NTT_B0R11;
  std::vector<ZZ> A_NTT_B0R12;
  std::vector<ZZ> A_NTT_B0R13;
  std::vector<ZZ> A_NTT_B0R14;
  std::vector<ZZ> A_NTT_B0R15;
  std::vector<ZZ> A_NTT_B1R0;
  std::vector<ZZ> A_NTT_B1R1;
  std::vector<ZZ> A_NTT_B1R2;
  std::vector<ZZ> A_NTT_B1R3;
  std::vector<ZZ> A_NTT_B1R4;
  std::vector<ZZ> A_NTT_B1R5;
  std::vector<ZZ> A_NTT_B1R6;
  std::vector<ZZ> A_NTT_B1R7;
  std::vector<ZZ> A_NTT_B1R8;
  std::vector<ZZ> A_NTT_B1R9;
  std::vector<ZZ> A_NTT_B1R10;
  std::vector<ZZ> A_NTT_B1R11;
  std::vector<ZZ> A_NTT_B1R12;
  std::vector<ZZ> A_NTT_B1R13;
  std::vector<ZZ> A_NTT_B1R14;
  std::vector<ZZ> A_NTT_B1R15;
 
  //--------------------
  //modify 2020/08/12
  A.resize(fft_point);
  A_1.resize(fft_point);
  
  for(int i = 0;i < fft_point;i++){
		  A[i]   = i;
		  A_1[i] = i;
  }

  switch(fft_point){
    case 65536:
      test.NTT_radix16(A);
      break;
    case 32768:
      test.NTT_r16_r8(A,A_NTT_B0R0,A_NTT_B0R1,A_NTT_B0R2,A_NTT_B0R3,
			    A_NTT_B0R4,A_NTT_B0R5,A_NTT_B0R6,A_NTT_B0R7,
                  A_NTT_B0R8,A_NTT_B0R9,A_NTT_B0R10,A_NTT_B0R11,
                  A_NTT_B0R12,A_NTT_B0R13,A_NTT_B0R14,A_NTT_B0R15,
                  A_NTT_B1R0,A_NTT_B1R1,A_NTT_B1R2,A_NTT_B1R3,
				A_NTT_B1R4,A_NTT_B1R5,A_NTT_B1R6,A_NTT_B1R7,
                  A_NTT_B1R8,A_NTT_B1R9,A_NTT_B1R10,A_NTT_B1R11,
				A_NTT_B1R12,A_NTT_B1R13,A_NTT_B1R14,A_NTT_B1R15);
      break;
    case 16384:
      test.NTT_r16_r4(A,A_NTT_B0R0,A_NTT_B0R1,A_NTT_B0R2,A_NTT_B0R3,
			    A_NTT_B0R4,A_NTT_B0R5,A_NTT_B0R6,A_NTT_B0R7,
                  A_NTT_B0R8,A_NTT_B0R9,A_NTT_B0R10,A_NTT_B0R11,
                  A_NTT_B0R12,A_NTT_B0R13,A_NTT_B0R14,A_NTT_B0R15,
                  A_NTT_B1R0,A_NTT_B1R1,A_NTT_B1R2,A_NTT_B1R3,
				A_NTT_B1R4,A_NTT_B1R5,A_NTT_B1R6,A_NTT_B1R7,
                  A_NTT_B1R8,A_NTT_B1R9,A_NTT_B1R10,A_NTT_B1R11,
				A_NTT_B1R12,A_NTT_B1R13,A_NTT_B1R14,A_NTT_B1R15);
      break;
    case 8192:
      test.NTT_r16_r2(A,A_NTT_B0R0,A_NTT_B0R1,A_NTT_B0R2,A_NTT_B0R3,
			    A_NTT_B0R4,A_NTT_B0R5,A_NTT_B0R6,A_NTT_B0R7,
                  A_NTT_B0R8,A_NTT_B0R9,A_NTT_B0R10,A_NTT_B0R11,
                  A_NTT_B0R12,A_NTT_B0R13,A_NTT_B0R14,A_NTT_B0R15,
                  A_NTT_B1R0,A_NTT_B1R1,A_NTT_B1R2,A_NTT_B1R3,
				A_NTT_B1R4,A_NTT_B1R5,A_NTT_B1R6,A_NTT_B1R7,
                  A_NTT_B1R8,A_NTT_B1R9,A_NTT_B1R10,A_NTT_B1R11,
				A_NTT_B1R12,A_NTT_B1R13,A_NTT_B1R14,A_NTT_B1R15);
      break;
    case 256:
      test.NTT_radix4(A);
      break;
    case 128:
      test.NTT_r4_r2(A,A_NTT_B0R0,A_NTT_B0R1,A_NTT_B0R2,A_NTT_B0R3,
              A_NTT_B1R0,A_NTT_B1R1,A_NTT_B1R2,A_NTT_B1R3);
      break;
    case 16:
      test.NTT_radix2(A);
      break;
  }
  
    //multipler
  NTT_test.NTT_t(A_1);
  
  std::ofstream A_o("./A_output.txt");
  std::ofstream A_1_o("./A_1_output.txt");
  std::ofstream A_INTT_o("./A_INTT_output.txt");
  //std::ofstream Mult_data_o("./MULT_DATA_o.txt");
  

  int error = 0;
   
  for(int i = 0; i < fft_point;i++){
	 A_o << A[i];  
	 A_1_o << A_1[i];  
     A_o << "\n";
     A_1_o << "\n";
	 if(A[i] != A_1[i]) {
		 std::cout << "error index: " << i <<"\n";
		 error = error + 1;
	 }
  }
  std::cout << "error : " << error << "\n";
 
  
  std::cout << "------------------------------\n";
  std::cout << " INTT Start!!                 \n";

  //test.INTT_r16_r2(A,A_NTT_B0R0,A_NTT_B0R1,A_NTT_B0R2,A_NTT_B0R3,
	//			     A_NTT_B0R4,A_NTT_B0R5,A_NTT_B0R6,A_NTT_B0R7,
  //                   A_NTT_B0R8,A_NTT_B0R9,A_NTT_B0R10,A_NTT_B0R11,
  //                   A_NTT_B0R12,A_NTT_B0R13,A_NTT_B0R14,A_NTT_B0R15,
  //                   A_NTT_B1R0,A_NTT_B1R1,A_NTT_B1R2,A_NTT_B1R3,
	//			     A_NTT_B1R4,A_NTT_B1R5,A_NTT_B1R6,A_NTT_B1R7,
  //                   A_NTT_B1R8,A_NTT_B1R9,A_NTT_B1R10,A_NTT_B1R11,
	//			     A_NTT_B1R12,A_NTT_B1R13,A_NTT_B1R14,A_NTT_B1R15);
  
  //test.NTT_radix16(A);
  /*test.INTT_r4_r2(A,A_NTT_B0R0,A_NTT_B0R1,A_NTT_B0R2,A_NTT_B0R3,
                A_NTT_B1R0,A_NTT_B1R1,A_NTT_B1R2,A_NTT_B1R3);*/
  
  for(int i = 0; i < fft_point;i++){
	A_INTT_o << A[i];  
	A_INTT_o << "\n";  
  }
  
}
