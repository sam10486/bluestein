// Copyright (C) IBM, All Rights Reserved
/*
  Transform verilog code to C code 
*/
#include <cstring>
#include <string>
#include <iomanip>
#include <ctime>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>

#include "FFTC.h"
/*
change FFT_point with radix-16  will change the rom data 
*/
void my_test(int argc, char *argv[])
{
  /******************************************************/
  //argv[0]:file name => ./test_file.cc
  //argv[1]: IsConfig
  //argv[2]: m-th cyclotomic polynomial 
  //argv[3]: cyclotomic polynomial coefficient bit size
  //argv[4]: FFT processor radix
  /*****************************************************/
  FFTC test_fftc;
  int IsRandom = 1;
  int IsConfig = 0;
  int IsRandom_cyclotomic_prime = 1;
  int Freq;
  std::cout << "----------------------------------------------------------\n";
  std::cout << "basic parameter:\n";
  if(argc > 1)std::cout  << "Is Configureable:  "   << argv[1] << "\n";
  if(argc > 2)std::cout  << "m-th: "                << argv[2] << "\n";
  if(argc > 3)std::cout  << "coefficient bit size:" << argv[3] << "\n";
  if(argc > 4)std::cout  << "FFTP radix is :"       << argv[4] << "\n";
  if(argc > 5)std::cout  << "Frequency : "          << argv[5] << "\n";
  if(argc > 6){std::cout << "File name : "          << argv[6] << "\n"; IsRandom = 0;}
  if(argc > 7){std::cout << "cyclotomic prime : "   << argv[7] << "\n"; IsRandom_cyclotomic_prime = 0;}
  std::cout << "------------------------------------------------------------\n";
  //
  long m ; //m = 1705
  unsigned long OP_width ; // = 22;
  unsigned long radix;
  std::string string_tmp = "./";
  std::string Data_in;
  IsConfig = atol(argv[1]);
  m        = atol(argv[2]);
  OP_width = atol(argv[3]);
  radix    = atol(argv[4]);
  Freq     = atol(argv[5]);
  if(argc > 6 ) {
	  Data_in = argv[6];
	  Data_in = string_tmp + Data_in;
  }
  long CP_tmp;
  if(argc > 7){
	CP_tmp = atol(argv[7]);  	  
  }
  //calculate FFT point
  unsigned long L1;
  double L2;
  double M; 
  L1 = 2*m - 2;
  L2 = log2(L1);
  M  = ceil(L2);
  M  = pow(2,M); // Power of 2 length
  std::cout << "need FFT_Point: " << M <<"\n";
  unsigned long FFT_Point ;
  FFT_Point = (unsigned long) M;
  
  double Stage;         // need stage for radix-r
  double Stage_ceil;
  double radix_bit;     // number of radix bit
  double FFT_Point_bit; // number of fft_point bits

  radix_bit     = log2((double)radix);
  FFT_Point_bit = log2(M);
  Stage = FFT_Point_bit / radix_bit;
  Stage_ceil = ceil(Stage);
  std::cout << "Stage : " << Stage <<" , Stage_ceil: " << Stage_ceil <<"\n";  
 
  //==================================================
  test_fftc.parameter_in(radix,FFT_Point,OP_width,m,CP_tmp,IsRandom_cyclotomic_prime,Freq);
  std::string        string_buf1 = "./BFFTP";
  std::string        string_in;
  std::stringstream  ss;
  ss << string_buf1;
  std::cout << "file fold : " << ss.str() << "\n"; 
  string_in = ss.str();
  if(IsConfig == 0){
     std::cout <<"-----------------------------------------------\n";
     std::cout <<"         Test file generate!!!!!!              \n";
     test_fftc.testingfile_gen(string_in,Data_in,IsRandom);
     std::cout <<"-----------------------------------------------\n";
     std::cout <<"         Verilog Module generate!!!!!!         \n";
     test_fftc.gen(string_in);
     std::cout <<"-----------------------------------------------\n";
  }else {
     //****************************************************************
     //Configureable generate function
     std::cout <<"-----------------------------------------------\n";
     std::cout <<"         Test file generate!!!!!!              \n";	 
     test_fftc.testingfile_Reconfigure_gen(string_in,Data_in,IsRandom);
	 std::cout <<"-----------------------------------------------\n";
     std::cout <<"         Verilog Module generate!!!!!!         \n";
	 test_fftc.FFTP_Reconfigure_gen(string_in);
	 test_fftc.testfftp_Reconfigure_r16_for_coverage("./Reconfigure");
     //*****************************************************************
  std::cout <<"-----------------------------------------------\n";
  }

}