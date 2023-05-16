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

  radix_r1 = 2;
  radix_r2 = 2;
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
  ZZ Phi, InvPhi, IW;
  SqrRootMod(Phi, fft_twiddle, fft_prime);
  InvMod(InvPhi, Phi, fft_prime);
  InvMod(IW, fft_twiddle, fft_prime);
  cout << "Phi = " << Phi << ", InvPhi = " << InvPhi << endl;
  NWC_Algo nwc_algo(radix_r1, radix_r2, fft_point, fft_prime);
  nwc_algo.setValue(radix_r1, radix_r2, fft_point, fft_prime, Phi, InvPhi, fft_twiddle, IW);
  nwc_algo.showInfo();
  vector<ZZ > NWC_arr, NWC_golden;
  NWC_arr.resize(fft_point);
  NWC_golden.resize(fft_point);
  for (int i = 0; i < fft_point; i++){
    NWC_golden[i] = i;
    NWC_arr[i] = i;
  }
  nwc_algo.NWC(NWC_arr);
  nwc_algo.INWC(NWC_arr);
  int err = 0;
  for (int i = 0; i < fft_point; i++){
    if (NWC_golden[i] != NWC_arr[i]){
      err++;
    }
  }
  cout << "err = " << err << endl;
  
  //------------------------------------------

  //----------For INTT-----------
  std::ofstream INTT_output("./NWC_PrintData/INTT_output.txt");
  NTT     INTT_instance;
  vector<ZZ > INTT_arr;
  INTT_arr.resize(fft_point);
  for(int i=0; i<fft_point; i++){
    INTT_arr[i] = i;
  }
  INTT_instance.NTT_init(fft_point,fft_prime,fft_IW);
  INTT_instance.NTT_t(INTT_arr);
  for (int i = 0; i < fft_point; i++){
    INTT_output << INTT_arr[i] << endl;
  }
  
  //-----------------------------

  std::cout << "test NTTSPMB Init!!! \n";
  //---
  test.init(fft_point,fft_prime,fft_twiddle,radix_r1);
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

  cout << "--------------------DIT FFT part----------------" << std::endl;
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
    cout << "fft_twiddle = " << fft_twiddle << endl;

    for(int i = 0;i < fft_point;i++){
      B[i]   = i;
    }
    DIT_spmb.init(fft_point,fft_prime,fft_twiddle,radix_r1);

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
        cout << "No this DIT selection" << endl;
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
  
  cout << "------------------DIF FFT test-----------------" << endl;
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
    cout << "fft_twiddle = " << fft_twiddle << endl;

    for(int i = 0;i < fft_point;i++){
      C[i]   = i;
    }
    DIF_spmb.init(fft_point,fft_prime,fft_twiddle,radix_r1);

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
        cout << "No this DIF selection" << endl;
        break;
    }
    
    std::ofstream C_o("./C_output.txt");
    for(int i = 0; i < fft_point;i++){
	    C_o << C[i];  
      C_o << "\n";
	    if(C[i] != A_1[i]) {
	  	  std::cout << "error index: " << i <<"\n";
	  	  error = error + 1;
      }
    }
    std::cout << "error : " << error << "\n";
  }
  
}
