#ifndef _MUX_H_
#define _MUX_H_
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

using namespace NTL;

class Mux{
public:
  unsigned long r; // radix
  unsigned long addr_width;
  //original prime , small prime for HeLib
  unsigned long CP_width;
  ZZ IN; // inverse N (mod p)
  
  //generate mux verilog file and init parameter
  void gen(unsigned long addr_w,unsigned long radix,unsigned long CP_w,std::string string_in ,ZZ inverse_N);
  //
  void parameter_init(unsigned long addr_w,unsigned long CP_w);
  // TWIMux
  void TWIMux_r4(std::string string_in);
  void TWIMux_r8(std::string string_in);
  void TWIMux_r16(std::string string_in);
  //for radix-4
  void Mux1_r4(std::string string_in);
  void Mux2_r4(std::string string_in);
  void Mux3_r4(std::string string_in);
  void Mux4_r4(std::string string_in);
  //bluestein's fft,new module     
  void Mux5_r4(std::string string_in);
  void Mux6_r4(std::string string_in);
  //==================================
  void MuxMA_r4(std::string string_in);
  void MuxROMA_r4(std::string string_in);
  //==================================
  //radix-8
  void Mux1_r8(std::string string_in);
  void Mux2_r8(std::string string_in);
  void Mux3_r8(std::string string_in);
  void Mux4_r8(std::string string_in);
  //bluestein's fft,new module     
  void Mux5_r8(std::string string_in);
  void Mux6_r8(std::string string_in);  
  //one parameter addr_width in mux1
  //for radix-16
  void Mux1(std::string string_in);
  void Mux2(std::string string_in);
  void Mux3(std::string string_in);
  void Mux4(std::string string_in);
  //bluestein's fft,new module
  void Mux5_r16(std::string string_in);
  void Mux6_r16(std::string string_in);
};
#endif