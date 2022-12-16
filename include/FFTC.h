#ifndef _FFTC_H_
#define _FFTC_H_
/*
  Using special prime to design FFT processor 
  prime = 2^64 - 2^32 + 1
  Sum radix-16 out for special prime FFT 
*/
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <sstream>
#include <NTL/ZZ.h>
//change parameter with fft point
#include "CenCtrl.h"
//radix-16
#include "CSA.h"
#include "Pipe.h"
#include "Sum.h"
//fixed 
//except Mux1
#include "Mux.h"
#include "CLA.h"
#include "Mod.h"
//
#include "bluestein.h"
#include "PhiFuc.h"
#include "NTT.h"
#include "SPMB.h"
//configurable BFFT 
#include "configurable_BFFT.h"

using namespace NTL;

//FFT Complier generate
class FFTC{
public:
  unsigned long radix;     //radix
  unsigned long fft_point; //fft_point 
  unsigned long CP_width;  // cyclotomic polynomial prime width
  unsigned long PRE_width;  // pre-computing width
  ZZ cyclotomic_prime;
  ZZ pre_computing;
  unsigned long addr_width;
  unsigned long ROM_width;
  unsigned long IReROM_width;
  unsigned long data_cnt_width;
  // Frequency
  int Freq;
  
  //cyclotomic polynomial parameter
  long m;    //m-th cyclotomic polynomial
  long phi_m; //phim
  long m_2;  //2m  
  
  ZZ IN;
  //parameter set
  //r is radix , N is fft_point
  void parameter_in(unsigned long r,unsigned long N,unsigned long CP_w,long m_th,long CP_in,int IsRandom_CP,int frequency);
  void cyclotpoly_parameter_set(long m_th,long CP_i,int IsRandom_CP); //frist call
  //synthesis script (design vision)
  void syn_script_gen();
  void syn_script_r4();
  void syn_script_r8();
  void syn_script_r16();
  //generate all module of FFT
  void gen(std::string string_in);
  void testingfile_gen(std::string string_in,std::string Data_in,int IsRandom);
  void testingData_gen(std::string Data_in , int random_data_option);
  void testingData_gen_for_coverage();
  void GoldenData_o();
  void mkflags_gen();
  //FFT processor top module generate
  void FFTP_gen(std::string string_in);
  //radix-4 FFT top moduler verilog file
  void FFTP_r4(std::string string_in);
  void testfftp_r4(std::string string_in);
  void testfftp_r4_R2P(std::string string_in);
  void testfftp_syn_r4();
  //Mixed radix-4 radix-2 
  void FFTP_r4_r2(std::string string_in);
  //radix-8 top modular verilog file
  void FFTP_r8(std::string string_in);
  void testfftp_r8(std::string string_in);
  void testfftp_r8_R2P(std::string string_in);
  void testfftp_syn_r8();  
  //radix-16 FFT top moduler verilog file
  void FFTP_r16(std::string string_in);
  void FFTP_r16_Mixed_radix(std::string string_in);
  void testfftp_r16(std::string string_in);
  void testfftp_r16_R2P(std::string string_in);
  void testfftp_syn_r16();
  //lc_shell script
  void lc_script_gen();
  void lc_script_gen_r4();
  void lc_script_gen_r8();
  void lc_script_gen_r16();
  //----------------------------
  //ReConfigurable BFFT generate function
  void FFTP_Reconfigure_gen(std::string string_in);
  void syn_script_Reconfigure_r16();
  void testingfile_Reconfigure_gen(std::string string_in,std::string Data_in,int IsRandom);
  void GoldenData_Reconfigure_o();
  void mkflags_Reconfigure_gen();
  void testfftp_Reconfigure_r16(std::string string_in);
  void testfftp_Reconfigure_r16_for_coverage(std::string string_in);
  void testfftp_Reconfigure_syn_r16();
  void lc_script_gen_Reconfigure_r16();
  //----------------------------
  std::string ZZtohex(ZZ zz_tmp,unsigned long data_bit_length); // using in cyclotomic_prime

};
#endif