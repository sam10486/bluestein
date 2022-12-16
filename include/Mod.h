#ifndef _MOD_H_
#define _MOD_H_
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

#include "SPMB.h"

using namespace NTL;

class Mod{
public:
  //int tmp;
  long m_2;
  unsigned long r;
  unsigned long cyclotomic_prime; // cyclotomic_prime prime 
  unsigned long CP_width; //cyclotomic prime width
  //Barrett reduction parameter
  long rf_fri;  //right shift frist  time
  long rf_sec;  //right shift second time
  long pre_computing;
  double data_mult_pre_par;
  long data_fri_rs_width;
  long data_mult_pre_width;
  
  //Generate Mod arithmetic verilog file
  void gen(unsigned long radix,std::string string_in,unsigned long CP_w,unsigned long cp_i,long m_2_i);
  void gen_configurable(unsigned long radix,std::string string_in,unsigned long CP_w,unsigned long cp_i,long m_2_i);
  void BM_parameter();  //barrett  reduction parameter generate
  //original prime modulus arithmetic
  void Mul(std::string string_in);
  void MulMod(std::string string_in);
  void MulMod_configurable(std::string string_in);
  void CP_CLA(std::string string_in);
  void CP_CLA_clg(std::string string_in,unsigned long Num_digits);
  void BR(std::string string_in);
  void BR_configurable(std::string string_in);
  //for radix-4
  //2020/03/10 modify   
  void Mod192_r4(std::string string_in);
  void Mod96_r4(std::string string_in);
  void MulMod128_r4(std::string string_in);
  void R4_TMulMod(std::string string_in);
  //radix-8 
  void Mod192_r8(std::string string_in);
  void Mod96_r8(std::string string_in);
  void MulMod128_r8(std::string string_in);
  //these parameter of any Mod module are fixed 
  //for radix-16
  void Mod96(std::string string_in);
  void Mod96PD(std::string string_in);
  void Mod192(std::string string_in);
  void Mod192PD(std::string string_in);
  void ModMux(std::string string_in);
  void Mul64(std::string string_in);
  void MulMod128(std::string string_in);
  void MulMod128PD(std::string string_in);
};
#endif