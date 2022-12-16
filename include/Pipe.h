#ifndef _PIPE_H_
#define _PIPE_H_
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

class Pipe{
public:
  unsigned long fft_point; //fft_point
  unsigned long r; //radix
  unsigned long addr_width;
  //generate radix-16 pipeline verilog file
  void gen(unsigned long fft_p,unsigned long radix,std::string string_in);
  //radix-4  pipeline verilog file
  //twiddle factor mux and pipe
  void TWIpipe_r4(std::string string_in);
  void R4_NPipeReg1(std::string string_in);
  void R4_NPipeReg2(std::string string_in);
  void R4_NPipeReg3(std::string string_in);
  void R4_PipeReg1(std::string string_in);
  void R4_PipeReg2(std::string string_in);
  void R4_PipeReg3(std::string string_in);
  void R4_PipeReg4(std::string string_in);
  void R4_PipeReg4_2(std::string string_in);
  void R4_PipeReg5_1(std::string string_in);
  void R4_PipeReg5_2(std::string string_in);
  void Radix4_Pipe(std::string string_in);
  void Radix4_R2_Pipe(std::string string_in);
  // radix-2^(2)
  void BU_R4_S0_R2P(std::string string_in);
  void BU_R4_R2P(std::string string_in);
  void Pipe_R4_R2P(std::string string_in);
  void Radix4_Pipe_R2P(std::string string_in);
  //radix-8  pipeline verilog file
  //twiddle factor mux and pipe
  void TWIpipe_r8(std::string string_in);
  void R8_NPipeReg1(std::string string_in);  
  void R8_NPipeReg2(std::string string_in);  
  void R8_NPipeReg3(std::string string_in);  
  void R8_PipeReg1(std::string string_in);  
  void R8_PipeReg2(std::string string_in);  
  void R8_PipeReg3(std::string string_in);  
  void R8_PipeReg4(std::string string_in);  
  void R8_PipeReg4_2(std::string string_in);  
  void R8_PipeReg5_1(std::string string_in);  
  void R8_PipeReg5_2(std::string string_in);  
  void Radix8_Pipe(std::string string_in);  
  // Radix-2^(4) 
  // Radix-2^(4) butterfly unit  
  void BU_R8_S0_R2P(std::string string_in);
  void BU_R8_R2P(std::string string_in);
  void Pipe_R8_R2P(std::string string_in);
  void Radix8_Pipe_R2P(std::string string_in);
  void TWIpipe_r8_R2P(std::string string_in);  
  
  //radix-16 pipeline verilog file
  //twiddle factor mux and pipe
  void TWIpipe_r16(std::string string_in);
  void R16_NPipeReg1(std::string string_in);
  void R16_NPipeReg2(std::string string_in);
  void R16_NPipeReg3(std::string string_in);
  void R16_PipeReg1(std::string string_in);
  void R16_PipeReg2(std::string string_in);
  void R16_PipeReg3(std::string string_in);
  void R16_PipeReg4(std::string string_in);
  void R16_PipeReg4_2(std::string string_in);
  void R16_PipeReg5_1(std::string string_in);
  void R16_PipeReg5_2(std::string string_in);
  void Radix16_Pipe(std::string string_in);
  void Radix16_Mixed_Radix_Pipe(std::string string_in);
  // Radix-2^(4) 
  // Radix-2^(4) butterfly unit
  void BU_R16_S0_R2P(std::string string_in);
  void BU_R16_R2P(std::string string_in);
  void Pipe_R2P(std::string string_in);
  void Radix16_Pipe_R2P(std::string string_in);
  void Radix16_Pipe_R2P_Mixed_Radix(std::string string_in);
  void TWIpipe_r16_R2P(std::string string_in);


};
#endif