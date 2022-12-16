#ifndef _CONFIGURABLE_BFFT_H_
#define _CONFIGURABLE_BFFT_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <sstream>
#include <NTL/ZZ.h>

using namespace NTL;

class configurable_BFFT{
public:
     unsigned long fft_point; //fft_point
     unsigned long CP_width;  // cyclotomic polynomial prime width
	 unsigned long PRE_width;
	 
   
     void gen(std::string string_in,unsigned long fft_p,unsigned long cp_width,unsigned long pre_width);
	 //TOP module 
	 void FFTP(std::string string_in);
	 //Cenctrl
	 void Config_Reg(std::string string_in); // configurable new module 
	 void Rectrl(std::string string_in);
	 void IReorderMA_pip(std::string string_in);
	 void R16_InpipeReg(std::string string_in);
	 void R16_BU_outpipe(std::string string_in);
     void CenCtrl_16(std::string string_in);
	 //no use
	 void CSM_4096(std::string string_in);
	 void CSM_8192(std::string string_in);
	 void CSM_16384(std::string string_in);
	 void CSM_32768(std::string string_in);
	 void CSM_65536(std::string string_in);
	 //
     void Ctrl_PipeReg1(std::string string_in);
     void R16_ROMPipeReg1(std::string string_in);
	 void order_ROMReg_r16(std::string string_in);
     void R16_AGU(std::string string_in);
     void R16_AGU_4096(std::string string_in);
     void R16_AGU_8192(std::string string_in);
     void R16_AGU_16384(std::string string_in);
     void R16_AGU_32768(std::string string_in);
     void R16_AGU_65536(std::string string_in);
     void R16_WD_buf(std::string string_in);
	 void R16_WAddr(std::string string_in);
     void R16_DC(std::string string_in);
	 
	 //Mux
	 //for radix-16
	 void MAMux_HSRAM(std::string string_in);
	 void MAMux_RESRAM(std::string string_in);
	 void TWIMux(std::string string_in);
     void Mux1(std::string string_in);
     void Mux2(std::string string_in);
     void Mux3(std::string string_in);
     void Mux4(std::string string_in);
     //bluestein's fft,new module
     void Mux5(std::string string_in);
     void Mux6(std::string string_in);	 
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
	 //Sum 
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
	 //CSA
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
	 // R2P Butterfly unit
	 void Ctrl_PipeReg1_R2P(std::string string_in); // Need to configre as Mixed radix
	 void R16_WD_buf_R2P(std::string string_in);
	 void R16_ROMPipeReg1_R2P(std::string string_in);
	 void BU_R16_S0_R2P(std::string string_in);
	 void BU_R16_R2P(std::string string_in);
	 void Pipe_R2P(std::string string_in);
	 void Radix16_Pipe_R2P(std::string string_in);
	 void TWIpipe_r16_R2P(std::string string_in);
     //CLA
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
     void CLA64clg_co(std::string string_in);
     void CLA65(std::string string_in);
     void CLA65clg(std::string string_in);
     void CLA96(std::string string_in);
     void CLA96clg(std::string string_in);
     void CLA192(std::string string_in);
     void CLA192clg(std::string string_in);   
	 

};
#endif