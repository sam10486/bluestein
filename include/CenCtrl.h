#ifndef _CENCTRL_H_
#define _CNECTRL_H_
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

class CenCtrl{
public:
  unsigned long radix; //radix 
  unsigned long fft_point; //fft_point
  unsigned long addr_width;
  unsigned long data_cnt_width;
  long          m;
  long          phim;
  /************************************************/
  //cyclotomic polynomial prime or small prime , this prime is used in HeLib
  unsigned long CP_width;
  /************************************************/
  double  IsMixed;
  
  /******parameter for center controll unit*****/
  unsigned long No_stage; // number of fft stage
  //number of load data.//For verilog parameter EXTMA_V1
  unsigned long NO_LD; 
  //number of using butterfly unit per stage. for verilog parameter DCNT_V2
  unsigned long NO_U_BUPS;
  //number of delay for Muldata. for verilog parameter DCNT_V3
  unsigned long NO_D_Muldata;
  unsigned long ROM_width;
  unsigned long IReROM_width; //reorder rom width
  //FFT Stage width
  unsigned long stage_width;
  //butterfly unit counter width
  unsigned long BU_counter_width;
  //data counter bit position1 for MSB of j_value
  //radix-r MSB bit position
  unsigned long DCNT_BP1;
  //data counter bit position2 for LSB of i_value
  //LSB of relocation groups width
  unsigned long DCNT_BP2;
  //data counter bit position3 for MSB of i_value
  //MSB of butterfly unit counter
  unsigned long DCNT_BP3;
  //data counter bit position4 for LSB of SC_value
  //LSB of Stage counter
  unsigned long DCNT_BP4;
  //FFT data counter. Using it to determine where stage of the FFT
  //FFT_dc[last one]  is for verilog parameter DCNT_V1
  //other in order for veilog parameter DCNT_V4~V*; * = 4 + (number of stage -1);
  std::vector<unsigned long> FFT_dc;
  //DCNT_FS1     = (FFT_Point / radix)*(stage-1) + 2
  //DCNT_FS(j+1) = (FFT_Point / radix)*(stage-1) + 2 + j      , for j = 0 ~ (radix-3)
  //DCNT_FS(K+1)  //data counter value(K+3) for FFT_Smode_sel , for K = 0 ~ (radix-3)
  std::vector<unsigned long> DCNT_FS; // no use for radix-16
  
  
  //Memory-based
  void parameter_set(unsigned long r,long m_th,long phi_m,unsigned long fft_p,unsigned long CP_w);
  //generate file
  void gen(unsigned long r,std::string string_in);
  //for radix-4
  //new module is added, in order to do bluestein's fft
  void IReorderMA_pip(std::string string_in);
  void Rectrl(std::string string_in);

  //===============================================
  //radix-4 
  void R4_InpipeReg(std::string string_in);
  void R4_BU_outpipe(std::string string_in);
  void order_ROMReg_r4(std::string string_in);
  void CenCtrl_r4(std::string string_in);
  void Ctrl_PipeReg1_r4(std::string string_in);
  void R4_ROMPipeReg1(std::string string_in);
  void R4_AGU(std::string string_in);
  void R4_WAddr(std::string string_in);
  void R4_WD_buf(std::string string_in);
  void R4_DC(std::string string_in);
  
  //radix-4 ,final stage using radix-2
  void Ctrl_PipeReg1_r4_r2(std::string string_in);
  void R4_R2_AGU(std::string string_in);
  void R4_R2_DC(std::string string_in);
  
  //radix-8
  void R8_InpipeReg(std::string string_in);
  void R8_BU_outpipe(std::string string_in);
  void order_ROMReg_r8(std::string string_in);
  void CenCtrl_r8(std::string string_in);
  void Ctrl_PipeReg1_r8(std::string string_in);
  void R8_ROMPipeReg1(std::string string_in);
  void R8_AGU(std::string string_in);
  void R8_WAddr(std::string string_in);
  void R8_WD_buf(std::string string_in);
  void R8_DC(std::string string_in);
  
  void Ctrl_PipeReg1_r8_R2P(std::string string_in);
  void R8_WD_buf_R2P(std::string string_in);
  void R8_ROMPipeReg1_R2P(std::string string_in);  

  //for radix-16
  void order_ROMReg_r16(std::string string_in);
  void R16_InpipeReg(std::string string_in);
  void R16_BU_outpipe(std::string string_in);
  void CenCtrl_16(std::string string_in);
  void Ctrl_PipeReg1(std::string string_in);
  void R16_ROMPipeReg1(std::string string_in);
  void R16_AGU(std::string string_in);
  void R16_WAddr(std::string string_in);
  void R16_WD_buf(std::string string_in);
  void R16_DC(std::string string_in);
  
  //Mixed radix-16 , final stage using radix-2 
  void Ctrl_PipeReg1_r16_r2(std::string string_in);
  void R16_R2_AGU(std::string string_in);
  void R16_R2_DC(std::string string_in);
  //Mixed radix-16 , final stage using radix-4 
  void Ctrl_PipeReg1_r16_r4(std::string string_in);
  void R16_R4_AGU(std::string string_in);
  void R16_R4_DC(std::string string_in);
  //Mixed radix-16 , final stage using radix-8 
  void Ctrl_PipeReg1_r16_r8(std::string string_in);
  void R16_R8_AGU(std::string string_in);
  void R16_R8_DC(std::string string_in);
  
  void Ctrl_PipeReg1_R2P(std::string string_in);
  void Ctrl_PipeReg1_R2P_Mixed_Radix(std::string string_in);
  void R16_WD_buf_R2P(std::string string_in);
  void R16_ROMPipeReg1_R2P(std::string string_in);
  
};
#endif