#ifndef _CSA_H_
#define _CSA_H_
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

class CSA{
public:
  unsigned long fft_point; //fft_point
  unsigned long r; //radix
  
  void gen(unsigned long fft_p,unsigned long radix , std::string string_in);
 //Carry Save adder for radix-4
 //Sum_CSAout for radix-4
 //string_in must be "./R16_ or R4_ after FFTPoint_FFT"
 //for example string_in = "./R16_65536P_FFT"
  void Sum4_CSAout0(std::string string_in);
  void Sum4_CSAout1(std::string string_in);
  void Sum4_CSAout2(std::string string_in);
  void Sum4_CSAout3(std::string string_in);
 //Mixed radix-4 radix-2 
  void Sum4_r2_CSAout0(std::string string_in);
  void Sum4_r2_CSAout1(std::string string_in);
  void Sum4_r2_CSAout2(std::string string_in);
  void Sum4_r2_CSAout3(std::string string_in);
 //Carry Save adder for radix-4
 //Sum_CSAout for radix-4
  void Sum8_CSAout0(std::string string_in);
  void Sum8_CSAout1(std::string string_in);
  void Sum8_CSAout2(std::string string_in);
  void Sum8_CSAout3(std::string string_in);
  void Sum8_CSAout4(std::string string_in);
  void Sum8_CSAout5(std::string string_in);
  void Sum8_CSAout6(std::string string_in);
  void Sum8_CSAout7(std::string string_in);
  
 //Carry save adder for radix-16
 //Sum_CSAout for radix-16
  void Sum16_CSAout0(std::string string_in);
  void Sum16_CSAout1(std::string string_in);
  void Sum16_CSAout2(std::string string_in);
  void Sum16_CSAout3(std::string string_in);
  void Sum16_CSAout4(std::string string_in);
  void Sum16_CSAout5(std::string string_in);
  void Sum16_CSAout6(std::string string_in);
  void Sum16_CSAout7(std::string string_in);
  void Sum16_CSAout8(std::string string_in);
  void Sum16_CSAout9(std::string string_in);
  void Sum16_CSAout10(std::string string_in);
  void Sum16_CSAout11(std::string string_in);
  void Sum16_CSAout12(std::string string_in);
  void Sum16_CSAout13(std::string string_in);
  void Sum16_CSAout14(std::string string_in);
  void Sum16_CSAout15(std::string string_in);
  //-----------------
  //Carry save adder for radix-16
  //Sum_CSAout for radix-16 and radix-2
  void Sum16_R2_CSAout0(std::string string_in);
  void Sum16_R2_CSAout1(std::string string_in);
  void Sum16_R2_CSAout2(std::string string_in);
  void Sum16_R2_CSAout3(std::string string_in);
  void Sum16_R2_CSAout4(std::string string_in);
  void Sum16_R2_CSAout5(std::string string_in);
  void Sum16_R2_CSAout6(std::string string_in);
  void Sum16_R2_CSAout7(std::string string_in);
  void Sum16_R2_CSAout8(std::string string_in);
  void Sum16_R2_CSAout9(std::string string_in);
  void Sum16_R2_CSAout10(std::string string_in);
  void Sum16_R2_CSAout11(std::string string_in);
  void Sum16_R2_CSAout12(std::string string_in);
  void Sum16_R2_CSAout13(std::string string_in);
  void Sum16_R2_CSAout14(std::string string_in);
  void Sum16_R2_CSAout15(std::string string_in);  
  //-----------------
  //Carry save adder for radix-16
  //Sum_CSAout for radix-16 and radix-4
  void Sum16_R4_CSAout0(std::string string_in);
  void Sum16_R4_CSAout1(std::string string_in);
  void Sum16_R4_CSAout2(std::string string_in);
  void Sum16_R4_CSAout3(std::string string_in);
  void Sum16_R4_CSAout4(std::string string_in);
  void Sum16_R4_CSAout5(std::string string_in);
  void Sum16_R4_CSAout6(std::string string_in);
  void Sum16_R4_CSAout7(std::string string_in);
  void Sum16_R4_CSAout8(std::string string_in);
  void Sum16_R4_CSAout9(std::string string_in);
  void Sum16_R4_CSAout10(std::string string_in);
  void Sum16_R4_CSAout11(std::string string_in);
  void Sum16_R4_CSAout12(std::string string_in);
  void Sum16_R4_CSAout13(std::string string_in);
  void Sum16_R4_CSAout14(std::string string_in);
  void Sum16_R4_CSAout15(std::string string_in); 
  //-----------------
  //Carry save adder for radix-16
  //Sum_CSAout for radix-16 and radix-8
  void Sum16_R8_CSAout0(std::string string_in);
  void Sum16_R8_CSAout1(std::string string_in);
  void Sum16_R8_CSAout2(std::string string_in);
  void Sum16_R8_CSAout3(std::string string_in);
  void Sum16_R8_CSAout4(std::string string_in);
  void Sum16_R8_CSAout5(std::string string_in);
  void Sum16_R8_CSAout6(std::string string_in);
  void Sum16_R8_CSAout7(std::string string_in);
  void Sum16_R8_CSAout8(std::string string_in);
  void Sum16_R8_CSAout9(std::string string_in);
  void Sum16_R8_CSAout10(std::string string_in);
  void Sum16_R8_CSAout11(std::string string_in);
  void Sum16_R8_CSAout12(std::string string_in);
  void Sum16_R8_CSAout13(std::string string_in);
  void Sum16_R8_CSAout14(std::string string_in);
  void Sum16_R8_CSAout15(std::string string_in); 

};
#endif