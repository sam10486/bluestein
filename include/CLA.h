#ifndef _CLA_H_
#define _CLA_H_
/*
  Using special prime to design FFT processor 
  prime = 2^64 - 2^32 + 1
  Carry lookahead generator for special prime FFT 
*/
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <sstream>

class CLA{
public:
  int r; //radix
  
  //generate carry lookahead generator verilog file
  void gen(std::string string_in);
  //Carry lookahead generator
  void CLA1(std::string string_in);
  void CLA2(std::string string_in);
  void CLA3(std::string string_in);
  void CLA4(std::string string_in);
  void CLA6(std::string string_in);
  void CLA16(std::string string_in);
  void CLA16clg(std::string string_in);
  void CLA24(std::string string_in);
  void CLA24clg(std::string string_in);
  void CLA32(std::string string_in);
  void CLA32clg(std::string string_in);
  void CLA64(std::string string_in);
  void CLA64_co(std::string string_in);
  void CLA64clg(std::string string_in);
  //need to modify function 2020/3/4
  void CLA64clg_co(std::string string_in);
  //
  void CLA65(std::string string_in);
  void CLA65clg(std::string string_in);
  void CLA96(std::string string_in);
  void CLA96clg(std::string string_in);
  void CLA144(std::string string_in);
  void CLA144clg(std::string string_in);
  void CLA192(std::string string_in);
  void CLA192clg(std::string string_in);
  void CLA432(std::string string_in);
  void CLA432clg(std::string string_in);  
};
#endif
