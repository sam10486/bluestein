#ifndef _SPMB_H_
#define _SPMB_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <NTL/ZZ.h>

using namespace NTL;

class SPMB{
public:
    unsigned long fft_point;
    unsigned long radix;
    unsigned long counter_iteration;
	unsigned long CP_width;
	double        IsMixed;
	long 		  m;
	long          phim;
    int           bc_width;
    int           offset;
    int           group;
    int           fft_point_bit_size;
    ZZ            cyclotomic_prime; //calculate test for BFFT
    
    //FFT SPMB 
    std::vector<int> ma;
    std::vector<int> bit_array_tmp;
    std::vector<int> bn;
    //Inverse FFT SPMB
    std::vector<int> ima;
    std::vector<int> ibit_array_tmp;
    std::vector<int> ibit_br; //bit reverse 
    std::vector<int> ibn;
	
    //parameter init
	//Calculate SPMB for FFT ,calculate bn and ma.
	void init(unsigned long fft_p,unsigned long r,int bc_width,unsigned long CP_w,ZZ cyclotomic_p,long m_th,long phi_m);
	void init_Reconfigure(unsigned long fft_p,unsigned long r,int bc_width,unsigned long CP_w,ZZ cyclotomic_p,long m_th,long phi_m);
	void time_o(std::vector<ZZ> time_data,std::string string_in);
	void H_freq_o(std::vector<ZZ> B_NTT);
	void re_order_factor(ZZ m_2_rou);
	
	//radix-4
    void init_r4(unsigned long fft_p , unsigned long r,int bc_width,unsigned long CP_w,
	ZZ cyclotomic_p,long m_th,long phi_m); 
    void init_r4_r2(unsigned long fft_p , unsigned long r,int bc_width,unsigned long CP_w,
	ZZ cyclotomic_p,long m_th,long phi_m); 
    void R4_R2_Pointwise_Mult_Index_AGU(int index,int &MA,int &BN);
    void time_o_r4(std::vector<ZZ> time_data,std::string string_in);
    void H_freq_o_r4(std::vector<ZZ> B_NTT);  //synthesis mux-8 //output the frequency tpye of B_NTT, and it is binary type
    void H_freq_o_r4_r2(std::vector<ZZ> B_NTT);  //synthesis mux-8 //output the frequency tpye of B_NTT, and it is binary type
    void H_freq_o_r4_Mux16(std::vector<ZZ> B_NTT); //synthesis mux-16 //output the frequency tpye of B_NTT, and it is binary type
    void H_freq_o_r4_r2_Mux16(std::vector<ZZ> B_NTT); //synthesis mux-16 //output the frequency tpye of B_NTT, and it is binary type
	void re_order_factor_r4(ZZ m_2_rou);  //reorder factor  output
    void r4_FFT_TW_ROM_D4();
	void r4_FFT_TW_ROM();   //
    void r4_IFFT_TW_ROM_D4();
	void r4_IFFT_TW_ROM();
	//radix-8
    void init_r8(unsigned long fft_p , unsigned long r,int bc_width,unsigned long CP_w,
	ZZ cyclotomic_p,long m_th,long phi_m);  
    void time_o_r8(std::vector<ZZ> time_data,std::string string_in);
    void H_freq_o_r8(std::vector<ZZ> B_NTT);  //synthesis mux-8 //output the frequency tpye of B_NTT, and it is binary type
    void re_order_factor_r8(ZZ m_2_rou);  //reorder factor  output
    void r8_FFT_TW_ROM();   //
    void r8_IFFT_TW_ROM();
	//void r8_FFT_TW_ROM_D4();
    //void r8_IFFT_TW_ROM_D4();
    //void H_freq_o_r8_Mux16(std::vector<ZZ> B_NTT); //synthesis mux-16 //output the frequency tpye of B_NTT, and it is binary type	
	//radix-16
    void init_r16(unsigned long fft_p , unsigned long r,int bc_width,unsigned long CP_w,
	ZZ cyclotomic_p,long m_th,long phi_m,int IsReconfig);
	void init_r16_Mixed_radix(unsigned long fft_p , unsigned long r,int bc_width,unsigned long CP_w,
	ZZ cyclotomic_p,long m_th,long phi_m,int IsReconfig);
	//Mixed radix point wise multplication index agu
	void R16_Mixed_Radix_Pointwise_Mult_Index_AGU(int index,int &MA,int &BN);
	void R16_R2_Index_AGU(int index,int &MA,int &BN);
	void R16_R4_Index_AGU(int index,int &MA,int &BN);
	void R16_R8_Index_AGU(int index,int &MA,int &BN);
    //time domain data using SPMB memory addressing	
    void time_o_r16(std::vector<ZZ> time_data,std::string string_in);
    void H_freq_o_r16(std::vector<ZZ> H_NTT);   //output the frequency tpye of B_NTT, and it is binary type
    void H_freq_o_r16_r2(std::vector<ZZ> H_NTT);   //output the frequency tpye of B_NTT, and it is binary type
    void H_freq_o_r16_r4(std::vector<ZZ> H_NTT);   //output the frequency tpye of B_NTT, and it is binary type
    void H_freq_o_r16_r8(std::vector<ZZ> H_NTT);   //output the frequency tpye of B_NTT, and it is binary type
    void re_order_factor_r16(ZZ m_2_rou);  //reorder factor  output
    void r16_FFT_TW_ROM();
    void r16_IFFT_TW_ROM();
    void r16_FFT_TW_ROM_Reconfig();
    void r16_IFFT_TW_ROM_Reconfig();	
	// 
	void radix16_Reconfigure_DATA(std::vector<ZZ> B_NTT,ZZ m_2_rou);
	void H_freq_Reconfigure_r16(std::vector<ZZ> H_NTT);
	void H_freq_Reconfigure_r16_r2(std::vector<ZZ> H_NTT);
	void H_freq_Reconfigure_r16_r4(std::vector<ZZ> H_NTT);
	void H_freq_Reconfigure_r16_r8(std::vector<ZZ> H_NTT);
	void re_order_factor_Reconfigure_r16(ZZ m_2_rou);
    std::string ZZtohex(ZZ zz_tmp);
	std::string ZZtohex_cp(ZZ zz_tmp); // using in cyclotomic_prime
	std::string BITtohex(std::vector<int> A_in,int length);
};

#endif