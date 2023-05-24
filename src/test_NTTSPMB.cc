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
#include "DIT_NTTSPMB.h"
#include "DIF_NTTSPMB.h" 
#include "NWC_Algo.h"
#include "DIF_INWC.h"

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
  int radix_r1;
  int radix_r2;
  //----------bluestein--------------
  //bluestein blue;
  //ZZ tmp_prime;
  //conv(tmp_prime, "97");
  //ZZ ROU;
  //unsigned int u_n = 16;
  //blue.N_ROU(tmp_prime, u_n, ROU);
  //cout << "ROU = " << ROU << std::endl;
  //---------------------------------

  radix_r1 = 16;
  radix_r2 = 16;
  ZZ fft_prime;
  ZZ fft_twiddle;
  ZZ fft_twiddle_16;
  ZZ fft_IW;
  ZZ fft_twiddle_65536;
  
  fft_point         = pow(radix_r1, 3) * radix_r2;//16
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
  //-------------NWC PART---------------------
  std::ofstream INWC_golden_o("./NWC_PrintData/INWC_golden.txt"); 
  ZZ Phi, InvPhi, IW;
  SqrRootMod(Phi, fft_twiddle, fft_prime);
  InvMod(InvPhi, Phi, fft_prime);
  InvMod(IW, fft_twiddle, fft_prime);
  std::cout << "Phi = " << Phi << ", InvPhi = " << InvPhi << endl;
  NWC_Algo nwc_algo(radix_r1, radix_r2, fft_point, fft_prime);
  nwc_algo.setValue(radix_r1, radix_r2, fft_point, fft_prime, Phi, InvPhi, fft_twiddle, IW);
  //nwc_algo.showInfo();
  vector<ZZ > NWC_arr, NWC_golden, INWC_golden;
  NWC_arr.resize(fft_point);
  NWC_golden.resize(fft_point);
  INWC_golden.resize(fft_point);
  for (int i = 0; i < fft_point; i++){
    NWC_golden[i] = i;
    NWC_arr[i] = i;
    INWC_golden[i] = i;
  }
  nwc_algo.NWC(NWC_arr);
  nwc_algo.INWC(NWC_arr);
  nwc_algo.INWC(INWC_golden);
  int err = 0;
  for (int i = 0; i < fft_point; i++){
    if (NWC_golden[i] != NWC_arr[i]){
      err++;
    }
  }
  std::cout << "--------------NWC function test----------------" << endl;
  std::cout << "err = " << err << endl;
  std::cout << "--------------NWC function test fin------------" << endl;
  
  for (int i = 0; i < fft_point; i++){
    INWC_golden_o << INWC_golden[i] << endl;
  }
  

  //------------------------------------------


  std::cout << "test NTTSPMB Init!!! \n";
  //---
  test.init(fft_point,fft_prime,fft_twiddle,radix_r1, Phi);
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

  switch(fft_point){
    case 32768:
      test.INTT_r16_r8(A,A_NTT_B0R0,A_NTT_B0R1,A_NTT_B0R2,A_NTT_B0R3,
    	    A_NTT_B0R4,A_NTT_B0R5,A_NTT_B0R6,A_NTT_B0R7,
                  A_NTT_B0R8,A_NTT_B0R9,A_NTT_B0R10,A_NTT_B0R11,
                  A_NTT_B0R12,A_NTT_B0R13,A_NTT_B0R14,A_NTT_B0R15,
                  A_NTT_B1R0,A_NTT_B1R1,A_NTT_B1R2,A_NTT_B1R3,
    		A_NTT_B1R4,A_NTT_B1R5,A_NTT_B1R6,A_NTT_B1R7,
                  A_NTT_B1R8,A_NTT_B1R9,A_NTT_B1R10,A_NTT_B1R11,
    		A_NTT_B1R12,A_NTT_B1R13,A_NTT_B1R14,A_NTT_B1R15);
      break;
    case 16384:
      test.INTT_r16_r4(A,A_NTT_B0R0,A_NTT_B0R1,A_NTT_B0R2,A_NTT_B0R3,
    	    A_NTT_B0R4,A_NTT_B0R5,A_NTT_B0R6,A_NTT_B0R7,
                  A_NTT_B0R8,A_NTT_B0R9,A_NTT_B0R10,A_NTT_B0R11,
                  A_NTT_B0R12,A_NTT_B0R13,A_NTT_B0R14,A_NTT_B0R15,
                  A_NTT_B1R0,A_NTT_B1R1,A_NTT_B1R2,A_NTT_B1R3,
    		A_NTT_B1R4,A_NTT_B1R5,A_NTT_B1R6,A_NTT_B1R7,
                  A_NTT_B1R8,A_NTT_B1R9,A_NTT_B1R10,A_NTT_B1R11,
    		A_NTT_B1R12,A_NTT_B1R13,A_NTT_B1R14,A_NTT_B1R15);
      break;
    case 8192:
      test.INTT_r16_r2(A,A_NTT_B0R0,A_NTT_B0R1,A_NTT_B0R2,A_NTT_B0R3,
    	    A_NTT_B0R4,A_NTT_B0R5,A_NTT_B0R6,A_NTT_B0R7,
                  A_NTT_B0R8,A_NTT_B0R9,A_NTT_B0R10,A_NTT_B0R11,
                  A_NTT_B0R12,A_NTT_B0R13,A_NTT_B0R14,A_NTT_B0R15,
                  A_NTT_B1R0,A_NTT_B1R1,A_NTT_B1R2,A_NTT_B1R3,
    		A_NTT_B1R4,A_NTT_B1R5,A_NTT_B1R6,A_NTT_B1R7,
                  A_NTT_B1R8,A_NTT_B1R9,A_NTT_B1R10,A_NTT_B1R11,
    		A_NTT_B1R12,A_NTT_B1R13,A_NTT_B1R14,A_NTT_B1R15);
      break;
    case 128:
      test.INTT_r4_r2(A,A_NTT_B0R0,A_NTT_B0R1,A_NTT_B0R2,A_NTT_B0R3,
              A_NTT_B1R0,A_NTT_B1R1,A_NTT_B1R2,A_NTT_B1R3);
      break;
    default:
      cout << "NO INTT function exist!!" << std::endl;
      break;
  }

  for(int i = 0; i < fft_point;i++){
	  A_INTT_o << A[i];  
	  A_INTT_o << "\n";  
  }

  std::cout << "--------------------DIT FFT part----------------" << std::endl;
  int debug_for_DIT = 1;
  if(debug_for_DIT){
    error = 0;
    DIT_NTTSPMB DIT_spmb;
    BitOperate BitRev;
    vector<ZZ> B;
    vector<ZZ> B_NTT_B0R0;
    vector<ZZ> B_NTT_B0R1;
    vector<ZZ> B_NTT_B0R2;
    vector<ZZ> B_NTT_B0R3;
    vector<ZZ> B_NTT_B0R4;
    vector<ZZ> B_NTT_B0R5;
    vector<ZZ> B_NTT_B0R6;
    vector<ZZ> B_NTT_B0R7;
    vector<ZZ> B_NTT_B0R8;
    vector<ZZ> B_NTT_B0R9;
    vector<ZZ> B_NTT_B0R10;
    vector<ZZ> B_NTT_B0R11;
    vector<ZZ> B_NTT_B0R12;
    vector<ZZ> B_NTT_B0R13;
    vector<ZZ> B_NTT_B0R14;
    vector<ZZ> B_NTT_B0R15;
    vector<ZZ> B_NTT_B1R0;
    vector<ZZ> B_NTT_B1R1;
    vector<ZZ> B_NTT_B1R2;
    vector<ZZ> B_NTT_B1R3;
    vector<ZZ> B_NTT_B1R4;
    vector<ZZ> B_NTT_B1R5;
    vector<ZZ> B_NTT_B1R6;
    vector<ZZ> B_NTT_B1R7;
    vector<ZZ> B_NTT_B1R8;
    vector<ZZ> B_NTT_B1R9;
    vector<ZZ> B_NTT_B1R10;
    vector<ZZ> B_NTT_B1R11;
    vector<ZZ> B_NTT_B1R12;
    vector<ZZ> B_NTT_B1R13;
    vector<ZZ> B_NTT_B1R14;
    vector<ZZ> B_NTT_B1R15;
    B.resize(fft_point);
    std::cout << "fft_twiddle = " << fft_twiddle << endl;

    for(int i = 0;i < fft_point;i++){
      B[i]   = i;
    }
    DIT_spmb.init(fft_point,fft_prime,fft_twiddle,radix_r1, Phi);

    int no_display = 0;
    switch(fft_point){
      case 65536:
        DIT_spmb.DIT_NTT_radix16(B);
        break;
      case 256:
        DIT_spmb.DIT_NTT_radix4(B);
        //DIT_spmb.test_radix4(B);
        break;
      /*case 128:
        DIT_spmb.DIT_NTT_r4_r2(
          B,
          B_NTT_B0R0,B_NTT_B0R1,B_NTT_B0R2,B_NTT_B0R3,
          B_NTT_B1R0,B_NTT_B1R1,B_NTT_B1R2,B_NTT_B1R3);
        break;*/
      case 16:
        DIT_spmb.DIT_NTT_radix2(B);
        //DIT_spmb.test_radix2(B);
        break;
      default:
        no_display = 1;
        std::cout << "No this DIT selection" << endl;
        break;
    }
    
    std::ofstream B_o("./B_output.txt");
    for(int i = 0; i < fft_point;i++){
	    B_o << B[i];  
      B_o << "\n";
	    if(B[i] != A_1[i]) {
	  	  if(!no_display) std::cout << "error index: " << i <<"\n";
	  	  error = error + 1;
      }else if (B[i] == A_1[i])
      {
        //cout << "B[" << i << "] = " << B[i] << endl;
      }
      
    }
    if(!no_display) std::cout << "error : " << error << "\n";
  }
  
  std::cout << "------------------DIF FFT test-----------------" << endl;
  int debug_for_DIF = 1;
  if(debug_for_DIF){
    error = 0;
    DIF_NTTSPMB DIF_spmb;
    BitOperate BitRev;
    std::vector<ZZ> C;
    std::vector<ZZ> C_NTT_B0R0;
    std::vector<ZZ> C_NTT_B0R1;
    std::vector<ZZ> C_NTT_B0R2;
    std::vector<ZZ> C_NTT_B0R3;
    std::vector<ZZ> C_NTT_B0R4;
    std::vector<ZZ> C_NTT_B0R5;
    std::vector<ZZ> C_NTT_B0R6;
    std::vector<ZZ> C_NTT_B0R7;
    std::vector<ZZ> C_NTT_B0R8;
    std::vector<ZZ> C_NTT_B0R9;
    std::vector<ZZ> C_NTT_B0R10;
    std::vector<ZZ> C_NTT_B0R11;
    std::vector<ZZ> C_NTT_B0R12;
    std::vector<ZZ> C_NTT_B0R13;
    std::vector<ZZ> C_NTT_B0R14;
    std::vector<ZZ> C_NTT_B0R15;
    std::vector<ZZ> C_NTT_B1R0;
    std::vector<ZZ> C_NTT_B1R1;
    std::vector<ZZ> C_NTT_B1R2;
    std::vector<ZZ> C_NTT_B1R3;
    std::vector<ZZ> C_NTT_B1R4;
    std::vector<ZZ> C_NTT_B1R5;
    std::vector<ZZ> C_NTT_B1R6;
    std::vector<ZZ> C_NTT_B1R7;
    std::vector<ZZ> C_NTT_B1R8;
    std::vector<ZZ> C_NTT_B1R9;
    std::vector<ZZ> C_NTT_B1R10;
    std::vector<ZZ> C_NTT_B1R11;
    std::vector<ZZ> C_NTT_B1R12;
    std::vector<ZZ> C_NTT_B1R13;
    std::vector<ZZ> C_NTT_B1R14;
    std::vector<ZZ> C_NTT_B1R15;

    C.resize(fft_point);
    std::cout << "fft_twiddle = " << fft_twiddle << endl;

    for(int i = 0;i < fft_point;i++){
      C[i]   = i;
    }
    DIF_spmb.init(fft_point,fft_prime,fft_twiddle,radix_r1, Phi);

    switch(fft_point){
      case 65536:
        DIF_spmb.DIF_NTT_radix16(C);
        break;
      case 32768:
        DIF_spmb.DIF_NTT_r16_r8(
          C,
          C_NTT_B0R0, C_NTT_B0R1, C_NTT_B0R2, C_NTT_B0R3,
			    C_NTT_B0R4, C_NTT_B0R5, C_NTT_B0R6, C_NTT_B0R7,
          C_NTT_B0R8, C_NTT_B0R9, C_NTT_B0R10,C_NTT_B0R11,
          C_NTT_B0R12,C_NTT_B0R13,C_NTT_B0R14,C_NTT_B0R15,
          C_NTT_B1R0, C_NTT_B1R1, C_NTT_B1R2, C_NTT_B1R3,
				  C_NTT_B1R4, C_NTT_B1R5, C_NTT_B1R6, C_NTT_B1R7,
          C_NTT_B1R8, C_NTT_B1R9, C_NTT_B1R10,C_NTT_B1R11,
				  C_NTT_B1R12,C_NTT_B1R13,C_NTT_B1R14,C_NTT_B1R15
        );
        break;
      case 16384:
        DIF_spmb.DIF_NTT_r16_r4(
          C,
          C_NTT_B0R0,   C_NTT_B0R1, C_NTT_B0R2, C_NTT_B0R3,
    	    C_NTT_B0R4,   C_NTT_B0R5, C_NTT_B0R6, C_NTT_B0R7,
          C_NTT_B0R8,   C_NTT_B0R9, C_NTT_B0R10,C_NTT_B0R11,
          C_NTT_B0R12,  C_NTT_B0R13,C_NTT_B0R14,C_NTT_B0R15,
          C_NTT_B1R0,   C_NTT_B1R1, C_NTT_B1R2, C_NTT_B1R3,
    		  C_NTT_B1R4,   C_NTT_B1R5, C_NTT_B1R6, C_NTT_B1R7,
          C_NTT_B1R8,   C_NTT_B1R9, C_NTT_B1R10,C_NTT_B1R11,
    		  C_NTT_B1R12,  C_NTT_B1R13,C_NTT_B1R14,C_NTT_B1R15
        );
        break;
      case 8192:
        DIF_spmb.DIF_NTT_r16_r2(
          C,          C_NTT_B0R0, C_NTT_B0R1, C_NTT_B0R2,A_NTT_B0R3,
			    C_NTT_B0R4, C_NTT_B0R5, C_NTT_B0R6, C_NTT_B0R7,
          C_NTT_B0R8, C_NTT_B0R9, C_NTT_B0R10,C_NTT_B0R11,
          C_NTT_B0R12,C_NTT_B0R13,C_NTT_B0R14,C_NTT_B0R15,
          C_NTT_B1R0, C_NTT_B1R1, C_NTT_B1R2, C_NTT_B1R3,
				  C_NTT_B1R4, C_NTT_B1R5, C_NTT_B1R6, C_NTT_B1R7,
          C_NTT_B1R8, C_NTT_B1R9, C_NTT_B1R10,C_NTT_B1R11,
				  C_NTT_B1R12,C_NTT_B1R13,C_NTT_B1R14,C_NTT_B1R15);
        break;
      case 256:
        DIF_spmb.DIF_NTT_radix4(C);
        break;
      case 128:
        DIF_spmb.DIF_NTT_r4_r2(C,C_NTT_B0R0,C_NTT_B0R1,C_NTT_B0R2,C_NTT_B0R3,
              C_NTT_B1R0,C_NTT_B1R1,C_NTT_B1R2,C_NTT_B1R3);
        break;
      case 16:
        DIF_spmb.DIF_NTT_radix2(C);
        break;
      default:
        std::cout << "No this DIF selection" << endl;
        break;
    }
    
    std::ofstream C_o("./C_output.txt");
    for(int i = 0; i < fft_point;i++){
	    C_o << C[i];  
      C_o << "\n";
	    if(C[i] != A_1[i]) {
	  	  std::cout << "error index: " << i <<"\n";
	  	  error = error + 1;
      }else{
        //std::cout << "DIF_arr[" << i  << "] = " <<  C[i] << "\n";
      }
    }
    std::cout << "error : " << error << "\n";
  }
  

  cout << "------------------DIF NWC test-----------------" << endl;
  std::vector<ZZ> INWC_B0R0;
  std::vector<ZZ> INWC_B0R1;
  std::vector<ZZ> INWC_B0R2;
  std::vector<ZZ> INWC_B0R3;
  std::vector<ZZ> INWC_B0R4;
  std::vector<ZZ> INWC_B0R5;
  std::vector<ZZ> INWC_B0R6;
  std::vector<ZZ> INWC_B0R7;
  std::vector<ZZ> INWC_B0R8;
  std::vector<ZZ> INWC_B0R9;
  std::vector<ZZ> INWC_B0R10;
  std::vector<ZZ> INWC_B0R11;
  std::vector<ZZ> INWC_B0R12;
  std::vector<ZZ> INWC_B0R13;
  std::vector<ZZ> INWC_B0R14;
  std::vector<ZZ> INWC_B0R15;

  std::vector<ZZ> INWC_B1R0;
  std::vector<ZZ> INWC_B1R1;
  std::vector<ZZ> INWC_B1R2;
  std::vector<ZZ> INWC_B1R3;
  std::vector<ZZ> INWC_B1R4;
  std::vector<ZZ> INWC_B1R5;
  std::vector<ZZ> INWC_B1R6;
  std::vector<ZZ> INWC_B1R7;
  std::vector<ZZ> INWC_B1R8;
  std::vector<ZZ> INWC_B1R9;
  std::vector<ZZ> INWC_B1R10;
  std::vector<ZZ> INWC_B1R11;
  std::vector<ZZ> INWC_B1R12;
  std::vector<ZZ> INWC_B1R13;
  std::vector<ZZ> INWC_B1R14;
  std::vector<ZZ> INWC_B1R15;


  DIF_INWC DIF_inwc;
  std::vector<ZZ> INWC_arr;
  INWC_arr.resize(fft_point);
  error = 0;
  std::cout << "fft_IW = " << fft_IW << endl;
  for(int i = 0;i < fft_point;i++){
    INWC_arr[i]   = i;
  }
  DIF_inwc.init(fft_point,fft_prime,fft_twiddle,radix_r1, Phi);
  switch (fft_point){
    case 65536:
      DIF_inwc.DIF_INWC_radix16(INWC_arr);
      break;
    case 32768:
       DIF_inwc.DIF_INWC_r16_r8(INWC_arr, INWC_B0R0, INWC_B0R1, INWC_B0R2,INWC_B0R3,
			    INWC_B0R4, INWC_B0R5, INWC_B0R6, INWC_B0R7,
          INWC_B0R8, INWC_B0R9, INWC_B0R10,INWC_B0R11,
          INWC_B0R12,INWC_B0R13,INWC_B0R14,INWC_B0R15,
          INWC_B1R0, INWC_B1R1, INWC_B1R2, INWC_B1R3,
				  INWC_B1R4, INWC_B1R5, INWC_B1R6, INWC_B1R7,
          INWC_B1R8, INWC_B1R9, INWC_B1R10,INWC_B1R11,
				  INWC_B1R12,INWC_B1R13,INWC_B1R14,INWC_B1R15);
      break;
    case 16384:
      DIF_inwc.DIF_INWC_r16_r4(INWC_arr, INWC_B0R0, INWC_B0R1, INWC_B0R2,INWC_B0R3,
			    INWC_B0R4, INWC_B0R5, INWC_B0R6, INWC_B0R7,
          INWC_B0R8, INWC_B0R9, INWC_B0R10,INWC_B0R11,
          INWC_B0R12,INWC_B0R13,INWC_B0R14,INWC_B0R15,
          INWC_B1R0, INWC_B1R1, INWC_B1R2, INWC_B1R3,
				  INWC_B1R4, INWC_B1R5, INWC_B1R6, INWC_B1R7,
          INWC_B1R8, INWC_B1R9, INWC_B1R10,INWC_B1R11,
				  INWC_B1R12,INWC_B1R13,INWC_B1R14,INWC_B1R15);
      break;
    case 8192:
      DIF_inwc.DIF_INWC_r16_r2(INWC_arr, INWC_B0R0, INWC_B0R1, INWC_B0R2,INWC_B0R3,
			    INWC_B0R4, INWC_B0R5, INWC_B0R6, INWC_B0R7,
          INWC_B0R8, INWC_B0R9, INWC_B0R10,INWC_B0R11,
          INWC_B0R12,INWC_B0R13,INWC_B0R14,INWC_B0R15,
          INWC_B1R0, INWC_B1R1, INWC_B1R2, INWC_B1R3,
				  INWC_B1R4, INWC_B1R5, INWC_B1R6, INWC_B1R7,
          INWC_B1R8, INWC_B1R9, INWC_B1R10,INWC_B1R11,
				  INWC_B1R12,INWC_B1R13,INWC_B1R14,INWC_B1R15);
      break;
    case 256:
      DIF_inwc.DIF_INWC_radix4(INWC_arr);
      break;
    case 128:
      DIF_inwc.DIF_INWC_r4_r2(INWC_arr, INWC_B0R0, INWC_B0R1, INWC_B0R2, INWC_B0R3,
                                INWC_B1R0, INWC_B1R1, INWC_B1R2, INWC_B1R3);
      break;
    case 16:
      DIF_inwc.DIF_INWC_radix2(INWC_arr);
      break;
    default:
      break;
  }
  std::ofstream DIF_INWC_o("./NWC_PrintData/DIF_INWC_output.txt");
    for(int i = 0; i < fft_point;i++){
	    DIF_INWC_o << INWC_arr[i];  
      DIF_INWC_o << "\n";
	    if(INWC_golden[i] != INWC_arr[i]) {
	  	  std::cout << "error index: " << i <<"\n";
	  	  error = error + 1;
      }else{
        //std::cout << "INWC_arr[" << i  << "] = " <<  INWC_arr[i] << "\n";
      }
    }
    std::cout << "error : " << error << "\n";


  //std::cout << "------------------DIF NWC MergeFactor test-----------------" << endl;
  //std::vector<ZZ> INWC_MergeFactor_B0R0;
  //std::vector<ZZ> INWC_MergeFactor_B0R1;
  //std::vector<ZZ> INWC_MergeFactor_B0R2;
  //std::vector<ZZ> INWC_MergeFactor_B0R3;
  //std::vector<ZZ> INWC_MergeFactor_B0R4;
  //std::vector<ZZ> INWC_MergeFactor_B0R5;
  //std::vector<ZZ> INWC_MergeFactor_B0R6;
  //std::vector<ZZ> INWC_MergeFactor_B0R7;
  //std::vector<ZZ> INWC_MergeFactor_B0R8;
  //std::vector<ZZ> INWC_MergeFactor_B0R9;
  //std::vector<ZZ> INWC_MergeFactor_B0R10;
  //std::vector<ZZ> INWC_MergeFactor_B0R11;
  //std::vector<ZZ> INWC_MergeFactor_B0R12;
  //std::vector<ZZ> INWC_MergeFactor_B0R13;
  //std::vector<ZZ> INWC_MergeFactor_B0R14;
  //std::vector<ZZ> INWC_MergeFactor_B0R15;
//
  //std::vector<ZZ> INWC_MergeFactor_B1R0;
  //std::vector<ZZ> INWC_MergeFactor_B1R1;
  //std::vector<ZZ> INWC_MergeFactor_B1R2;
  //std::vector<ZZ> INWC_MergeFactor_B1R3;
  //std::vector<ZZ> INWC_MergeFactor_B1R4;
  //std::vector<ZZ> INWC_MergeFactor_B1R5;
  //std::vector<ZZ> INWC_MergeFactor_B1R6;
  //std::vector<ZZ> INWC_MergeFactor_B1R7;
  //std::vector<ZZ> INWC_MergeFactor_B1R8;
  //std::vector<ZZ> INWC_MergeFactor_B1R9;
  //std::vector<ZZ> INWC_MergeFactor_B1R10;
  //std::vector<ZZ> INWC_MergeFactor_B1R11;
  //std::vector<ZZ> INWC_MergeFactor_B1R12;
  //std::vector<ZZ> INWC_MergeFactor_B1R13;
  //std::vector<ZZ> INWC_MergeFactor_B1R14;
  //std::vector<ZZ> INWC_MergeFactor_B1R15;
//
//
  //DIF_INWC DIF_inwc_MergeFactor;
  //std::vector<ZZ> INWC_MergeFactor_arr;
  //INWC_MergeFactor_arr.resize(fft_point);
  //error = 0;
  //std::cout << "fft_IW = " << fft_IW << endl;
  //for(int i = 0;i < fft_point;i++){
  //  INWC_MergeFactor_arr[i]   = i;
  //}
  //DIF_inwc_MergeFactor.init(fft_point,fft_prime,fft_twiddle,radix_r1, Phi);
  //switch (fft_point){
  //  case 65536:
  //    
  //    break;
  //  case 32768:
  //     
  //    break;
  //  case 16384:
  //    
  //    break;
  //  case 8192:
  //    
  //    break;
  //  case 256:
  //    
  //    break;
  //  case 128:
  //    
  //    break;
  //  case 16:
  //    DIF_inwc_MergeFactor.DIF_INWC_MergeFactor_radix2(INWC_MergeFactor_arr);
  //    break;
  //  default:
  //    break;
  //}
  //std::ofstream DIF_INWC_MergeFactor_o("./NWC_PrintData/DIF_INWC_MergeFactor_output.txt");
  //  for(int i = 0; i < fft_point;i++){
	//    DIF_INWC_MergeFactor_o << INWC_MergeFactor_arr[i];  
  //    DIF_INWC_MergeFactor_o << "\n";
	//    if(INWC_golden[i] != INWC_MergeFactor_arr[i]) {
	//  	  std::cout << "error index: " << i <<"\n";
	//  	  error = error + 1;
  //    }else{
  //      //std::cout << "INWC_MergeFactor_arr[" << i  << "] = " <<  INWC_MergeFactor_arr[i] << "\n";
  //    }
  //  }
  //  std::cout << "error : " << error << "\n";

  std::cout << "------------------DIF NWC seperateInvN test-----------------" << endl;
  std::vector<ZZ> INWC_seperateInvN_B0R0;
  std::vector<ZZ> INWC_seperateInvN_B0R1;
  std::vector<ZZ> INWC_seperateInvN_B0R2;
  std::vector<ZZ> INWC_seperateInvN_B0R3;
  std::vector<ZZ> INWC_seperateInvN_B0R4;
  std::vector<ZZ> INWC_seperateInvN_B0R5;
  std::vector<ZZ> INWC_seperateInvN_B0R6;
  std::vector<ZZ> INWC_seperateInvN_B0R7;
  std::vector<ZZ> INWC_seperateInvN_B0R8;
  std::vector<ZZ> INWC_seperateInvN_B0R9;
  std::vector<ZZ> INWC_seperateInvN_B0R10;
  std::vector<ZZ> INWC_seperateInvN_B0R11;
  std::vector<ZZ> INWC_seperateInvN_B0R12;
  std::vector<ZZ> INWC_seperateInvN_B0R13;
  std::vector<ZZ> INWC_seperateInvN_B0R14;
  std::vector<ZZ> INWC_seperateInvN_B0R15;

  std::vector<ZZ> INWC_seperateInvN_B1R0;
  std::vector<ZZ> INWC_seperateInvN_B1R1;
  std::vector<ZZ> INWC_seperateInvN_B1R2;
  std::vector<ZZ> INWC_seperateInvN_B1R3;
  std::vector<ZZ> INWC_seperateInvN_B1R4;
  std::vector<ZZ> INWC_seperateInvN_B1R5;
  std::vector<ZZ> INWC_seperateInvN_B1R6;
  std::vector<ZZ> INWC_seperateInvN_B1R7;
  std::vector<ZZ> INWC_seperateInvN_B1R8;
  std::vector<ZZ> INWC_seperateInvN_B1R9;
  std::vector<ZZ> INWC_seperateInvN_B1R10;
  std::vector<ZZ> INWC_seperateInvN_B1R11;
  std::vector<ZZ> INWC_seperateInvN_B1R12;
  std::vector<ZZ> INWC_seperateInvN_B1R13;
  std::vector<ZZ> INWC_seperateInvN_B1R14;
  std::vector<ZZ> INWC_seperateInvN_B1R15;


  DIF_INWC DIF_inwc_seperateInvN;
  std::vector<ZZ> INWC_seperateInvN_arr;
  INWC_seperateInvN_arr.resize(fft_point);
  error = 0;
  cout << "fft_IW = " << fft_IW << endl;
  for(int i = 0;i < fft_point;i++){
    INWC_seperateInvN_arr[i]   = i;
  }
  DIF_inwc_seperateInvN.init(fft_point,fft_prime,fft_twiddle,radix_r1, Phi);
  switch (fft_point){
    case 65536:
      DIF_inwc_seperateInvN.DIF_INWC_seperateInvN_radix16(INWC_seperateInvN_arr);
      break;
    case 32768:
       DIF_inwc_seperateInvN.DIF_INWC_seperateInvN_r16_r8(INWC_seperateInvN_arr, 
          INWC_seperateInvN_B0R0, INWC_seperateInvN_B0R1, INWC_seperateInvN_B0R2, INWC_seperateInvN_B0R3,
			    INWC_seperateInvN_B0R4, INWC_seperateInvN_B0R5, INWC_seperateInvN_B0R6, INWC_seperateInvN_B0R7,
          INWC_seperateInvN_B0R8, INWC_seperateInvN_B0R9, INWC_seperateInvN_B0R10,INWC_seperateInvN_B0R11,
          INWC_seperateInvN_B0R12,INWC_seperateInvN_B0R13,INWC_seperateInvN_B0R14,INWC_seperateInvN_B0R15,
          INWC_seperateInvN_B1R0, INWC_seperateInvN_B1R1, INWC_seperateInvN_B1R2, INWC_seperateInvN_B1R3,
				  INWC_seperateInvN_B1R4, INWC_seperateInvN_B1R5, INWC_seperateInvN_B1R6, INWC_seperateInvN_B1R7,
          INWC_seperateInvN_B1R8, INWC_seperateInvN_B1R9, INWC_seperateInvN_B1R10,INWC_seperateInvN_B1R11,
				  INWC_seperateInvN_B1R12,INWC_seperateInvN_B1R13,INWC_seperateInvN_B1R14,INWC_seperateInvN_B1R15);
      break;
    case 16384:
      DIF_inwc_seperateInvN.DIF_INWC_seperateInvN_r16_r4(INWC_seperateInvN_arr, 
          INWC_seperateInvN_B0R0, INWC_seperateInvN_B0R1, INWC_seperateInvN_B0R2, INWC_seperateInvN_B0R3,
			    INWC_seperateInvN_B0R4, INWC_seperateInvN_B0R5, INWC_seperateInvN_B0R6, INWC_seperateInvN_B0R7,
          INWC_seperateInvN_B0R8, INWC_seperateInvN_B0R9, INWC_seperateInvN_B0R10,INWC_seperateInvN_B0R11,
          INWC_seperateInvN_B0R12,INWC_seperateInvN_B0R13,INWC_seperateInvN_B0R14,INWC_seperateInvN_B0R15,
          INWC_seperateInvN_B1R0, INWC_seperateInvN_B1R1, INWC_seperateInvN_B1R2, INWC_seperateInvN_B1R3,
				  INWC_seperateInvN_B1R4, INWC_seperateInvN_B1R5, INWC_seperateInvN_B1R6, INWC_seperateInvN_B1R7,
          INWC_seperateInvN_B1R8, INWC_seperateInvN_B1R9, INWC_seperateInvN_B1R10,INWC_seperateInvN_B1R11,
				  INWC_seperateInvN_B1R12,INWC_seperateInvN_B1R13,INWC_seperateInvN_B1R14,INWC_seperateInvN_B1R15);
      break;
    case 8192:
      DIF_inwc_seperateInvN.DIF_INWC_seperateInvN_r16_r2(INWC_seperateInvN_arr, 
          INWC_seperateInvN_B0R0, INWC_seperateInvN_B0R1, INWC_seperateInvN_B0R2, INWC_seperateInvN_B0R3,
			    INWC_seperateInvN_B0R4, INWC_seperateInvN_B0R5, INWC_seperateInvN_B0R6, INWC_seperateInvN_B0R7,
          INWC_seperateInvN_B0R8, INWC_seperateInvN_B0R9, INWC_seperateInvN_B0R10,INWC_seperateInvN_B0R11,
          INWC_seperateInvN_B0R12,INWC_seperateInvN_B0R13,INWC_seperateInvN_B0R14,INWC_seperateInvN_B0R15,
          INWC_seperateInvN_B1R0, INWC_seperateInvN_B1R1, INWC_seperateInvN_B1R2, INWC_seperateInvN_B1R3,
				  INWC_seperateInvN_B1R4, INWC_seperateInvN_B1R5, INWC_seperateInvN_B1R6, INWC_seperateInvN_B1R7,
          INWC_seperateInvN_B1R8, INWC_seperateInvN_B1R9, INWC_seperateInvN_B1R10,INWC_seperateInvN_B1R11,
				  INWC_seperateInvN_B1R12,INWC_seperateInvN_B1R13,INWC_seperateInvN_B1R14,INWC_seperateInvN_B1R15);
      break;
    case 256:
      DIF_inwc_seperateInvN.DIF_INWC_seperateInvN_radix4(INWC_seperateInvN_arr);
      break;
    case 128:
      cout << "No this one!" << endl;
      break;
    case 16:
      DIF_inwc_seperateInvN.DIF_INWC_seperateInvN_radix2(INWC_seperateInvN_arr);
      break;
    default:
      break;
  }
  std::ofstream DIF_INWC_seperateInvN_o("./NWC_PrintData/DIF_INWC_seperateInvN_output.txt");
    for(int i = 0; i < fft_point;i++){
	    DIF_INWC_seperateInvN_o << INWC_seperateInvN_arr[i];  
      DIF_INWC_seperateInvN_o << "\n";
	    if(INWC_golden[i] != INWC_seperateInvN_arr[i]) {
	  	  std::cout << "error index: " << i <<"\n";
	  	  error = error + 1;
      }else{
        //std::cout << "INWC_seperateInvN_arr[" << i  << "] = " <<  INWC_seperateInvN_arr[i] << "\n";
      }
    }
    std::cout << "error : " << error << "\n";

  //---------------test R16 BU---------------
  int test_N = 8;
  vector<ZZ > R16_arr;
  R16_arr.resize(test_N);
   vector<ZZ > R16_golden;
  R16_golden.resize(test_N);
  for(int i=0; i<test_N; i++){
    R16_arr[i] = i;
    R16_golden[i] = i;
  }

  //DIF_inwc_seperateInvN.INWC_seperateInvN_Radix8_BU(
  //  R16_arr[0], R16_arr[1], R16_arr[2], R16_arr[3], 
  //  R16_arr[4], R16_arr[5], R16_arr[6], R16_arr[7]
  //);

  //DIF_inwc_seperateInvN.INWC_seperateInvN_Radix16_BU(
  //  R16_arr[0], R16_arr[1], R16_arr[2], R16_arr[3], 
  //  R16_arr[4], R16_arr[5], R16_arr[6], R16_arr[7], 
  //  R16_arr[8], R16_arr[9], R16_arr[10], R16_arr[11],
  //  R16_arr[12], R16_arr[13], R16_arr[14], R16_arr[15]
  //);

  //for (int i = 0; i < test_N; i++){
  //  cout << "R16_arr[" << i << "] = " << R16_arr[i] << endl;
  //}
  //
  //NTT R16_test_golden;
  //R16_test_golden.NTT_init(test_N,(ZZ)97,(ZZ)8);
  //R16_test_golden.NTT_t(R16_golden);
  //error = 0;
  //for(int i = 0; i < test_N;i++){
	// if(R16_arr[i] != R16_golden[i]) {
	//	 std::cout << "error index: " << i <<"\n";
  //   cout << "R16_arr = " << R16_arr[i] << " != " << "R16_golden = " << R16_golden[i] << endl;
	//	 error = error + 1;
	// }
  //}
  //std::cout << "error : " << error << "\n";
}
