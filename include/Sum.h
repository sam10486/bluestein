#ifndef _SUM_H_
#define _SUM_H_
/*
  Using special prime to design FFT processor 
  prime = 2^64 - 2^32 + 1
  Sum radix-16 or radix-4 out for special prime FFT 
*/
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <sstream>


class Sum{
public:
  unsigned long fft_point; //fft_point
  unsigned long r; //radix
  
  //generate radix-16 sum_out verilog file
  void gen(unsigned long fft_p,unsigned long radix,std::string string_in);
 //Sum_out for radix-4
  void Sum4_out0(std::string string_in);
  void Sum4_out1(std::string string_in);
  void Sum4_out2(std::string string_in);
  void Sum4_out3(std::string string_in);
  //Mixed-radix
  void Sum4_r2_out0(std::string string_in);
  void Sum4_r2_out1(std::string string_in);
  void Sum4_r2_out2(std::string string_in);
  void Sum4_r2_out3(std::string string_in);
 //Sum_out for radix-8
  void Sum8_out0(std::string string_in); 
  void Sum8_out1(std::string string_in); 
  void Sum8_out2(std::string string_in); 
  void Sum8_out3(std::string string_in); 
  void Sum8_out4(std::string string_in); 
  void Sum8_out5(std::string string_in); 
  void Sum8_out6(std::string string_in); 
  void Sum8_out7(std::string string_in); 
  
 //Sum_out for radix-16 
  void Sum16_out0(std::string string_in);
  void Sum16_out1(std::string string_in);
  void Sum16_out2(std::string string_in);
  void Sum16_out3(std::string string_in);
  void Sum16_out4(std::string string_in);
  void Sum16_out5(std::string string_in);
  void Sum16_out6(std::string string_in);
  void Sum16_out7(std::string string_in);
  void Sum16_out8(std::string string_in);
  void Sum16_out9(std::string string_in);
  void Sum16_out10(std::string string_in);
  void Sum16_out11(std::string string_in);
  void Sum16_out12(std::string string_in);
  void Sum16_out13(std::string string_in);
  void Sum16_out14(std::string string_in);
  void Sum16_out15(std::string string_in);
 //Sum_out for Mixed radix-16 
  void Sum16_out0_Mixed(std::string string_in);
  void Sum16_out1_Mixed(std::string string_in);
  void Sum16_out2_Mixed(std::string string_in);
  void Sum16_out3_Mixed(std::string string_in);
  void Sum16_out4_Mixed(std::string string_in);
  void Sum16_out5_Mixed(std::string string_in);
  void Sum16_out6_Mixed(std::string string_in);
  void Sum16_out7_Mixed(std::string string_in);
  void Sum16_out8_Mixed(std::string string_in);
  void Sum16_out9_Mixed(std::string string_in);
  void Sum16_out10_Mixed(std::string string_in);
  void Sum16_out11_Mixed(std::string string_in);
  void Sum16_out12_Mixed(std::string string_in);
  void Sum16_out13_Mixed(std::string string_in);
  void Sum16_out14_Mixed(std::string string_in);
  void Sum16_out15_Mixed(std::string string_in);  
};
#endif