#ifndef _NTTSPMB_H_
#define _NTTSPMB_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <NTL/ZZ.h>

using namespace NTL;

class NTTSPMB{
public:
	unsigned long  N; // NTT ponit
	int radix;

	ZZ N_ZZ; 
	ZZ IN; // inverse of N (mod p)

	ZZ  p; // modulus of NTT
	       // must be a prime
	ZZ  W; // primitive nth root of unity in Zp
	       // W^(N) (mod p) = 1
	ZZ  IW;// inverse of W  (mod p)
		
	void init(unsigned long n, ZZ prime, ZZ root,int r); //init parameters  
	void Radix2_BU(ZZ &a, ZZ &b); //Radix-2 buterfly unit
	void Radix4_BU(ZZ &a, ZZ &b,ZZ &c,ZZ &d); //Radix-4 buterfly unit
	void Radix4_BU_INTT(ZZ &a, ZZ &b,ZZ &c,ZZ &d); //Radix-4 buterfly unit
	void Radix8_BU(ZZ &a_r0, ZZ &a_r1,ZZ &a_r2,ZZ &a_r3,
	               ZZ &a_r4, ZZ &a_r5,ZZ &a_r6,ZZ &a_r7); //Radix-8 buterfly unit
    void Radix8_BU_INTT(ZZ &a_r0, ZZ &a_r1,ZZ &a_r2,ZZ &a_r3,				   
	                    ZZ &a_r4, ZZ &a_r5,ZZ &a_r6,ZZ &a_r7); //Radix-8 buterfly unit INTT			   
	void Radix16_BU(ZZ &a_r0, ZZ &a_r1,ZZ &a_r2,ZZ &a_r3,
					ZZ &a_r4, ZZ &a_r5,ZZ &a_r6,ZZ &a_r7,
					ZZ &a_r8, ZZ &a_r9,ZZ &a_r10,ZZ &a_r11,
					ZZ &a_r12, ZZ &a_r13,ZZ &a_r14,ZZ &a_r15);
    void Radix16_BU_INTT(ZZ &a_r0, ZZ &a_r1,ZZ &a_r2,ZZ &a_r3,
					     ZZ &a_r4, ZZ &a_r5,ZZ &a_r6,ZZ &a_r7,
					     ZZ &a_r8, ZZ &a_r9,ZZ &a_r10,ZZ &a_r11,				
					     ZZ &a_r12, ZZ &a_r13,ZZ &a_r14,ZZ &a_r15);	
    void NTT_radix2(std::vector<ZZ> &A); //radix-2 NTT
	void NTT_radix4(std::vector<ZZ> &A); //radix-4 NTT
	void NTT_radix16(std::vector<ZZ> &A); //radix-4 NTT
	void NTT_r4_r2(std::vector<ZZ> &A,std::vector<ZZ> &B0R0,std::vector<ZZ> &B0R1,std::vector<ZZ> &B0R2,std::vector<ZZ> &B0R3
	,std::vector<ZZ> &B1R0,std::vector<ZZ> &B1R1,std::vector<ZZ> &B1R2,std::vector<ZZ> &B1R3); //radix-4,radix-2 NTT
	void INTT_r4_r2(std::vector<ZZ> &A,std::vector<ZZ> &B0R0,std::vector<ZZ> &B0R1,std::vector<ZZ> &B0R2,std::vector<ZZ> &B0R3,
	std::vector<ZZ> &B1R0,std::vector<ZZ> &B1R1,std::vector<ZZ> &B1R2,std::vector<ZZ> &B1R3); //radix-4,radix-2 INTT
	//*************************************************************
	//Mixed radix-16 ,and the final stage is radix-2 FFT 
	//for 512-point , 8192 point ..... , but 131072 this parameters need to change twiddle factor.
	//*************************************************************
	//radix-16 and radix-2 Mixed radix NTT
	void NTT_r16_r2(std::vector<ZZ> &A,std::vector<ZZ> &B0R0,std::vector<ZZ> &B0R1,std::vector<ZZ> &B0R2,std::vector<ZZ> &B0R3
	                ,std::vector<ZZ> &B0R4,std::vector<ZZ> &B0R5,std::vector<ZZ> &B0R6,std::vector<ZZ> &B0R7,std::vector<ZZ> &B0R8
					,std::vector<ZZ> &B0R9,std::vector<ZZ> &B0R10,std::vector<ZZ> &B0R11,std::vector<ZZ> &B0R12,std::vector<ZZ> &B0R13
					,std::vector<ZZ> &B0R14,std::vector<ZZ> &B0R15,std::vector<ZZ> &B1R0,std::vector<ZZ> &B1R1,std::vector<ZZ> &B1R2
					,std::vector<ZZ> &B1R3,std::vector<ZZ> &B1R4,std::vector<ZZ> &B1R5,std::vector<ZZ> &B1R6,std::vector<ZZ> &B1R7
					,std::vector<ZZ> &B1R8,std::vector<ZZ> &B1R9,std::vector<ZZ> &B1R10,std::vector<ZZ> &B1R11,std::vector<ZZ> &B1R12
					,std::vector<ZZ> &B1R13,std::vector<ZZ> &B1R14,std::vector<ZZ> &B1R15); 
	//radix-16 and radix-2 Mixed radix NTT
	void INTT_r16_r2(std::vector<ZZ> &A,std::vector<ZZ> &B0R0,std::vector<ZZ> &B0R1,std::vector<ZZ> &B0R2,std::vector<ZZ> &B0R3
	                ,std::vector<ZZ> &B0R4,std::vector<ZZ> &B0R5,std::vector<ZZ> &B0R6,std::vector<ZZ> &B0R7,std::vector<ZZ> &B0R8
					,std::vector<ZZ> &B0R9,std::vector<ZZ> &B0R10,std::vector<ZZ> &B0R11,std::vector<ZZ> &B0R12,std::vector<ZZ> &B0R13
					,std::vector<ZZ> &B0R14,std::vector<ZZ> &B0R15,std::vector<ZZ> &B1R0,std::vector<ZZ> &B1R1,std::vector<ZZ> &B1R2
					,std::vector<ZZ> &B1R3,std::vector<ZZ> &B1R4,std::vector<ZZ> &B1R5,std::vector<ZZ> &B1R6,std::vector<ZZ> &B1R7
					,std::vector<ZZ> &B1R8,std::vector<ZZ> &B1R9,std::vector<ZZ> &B1R10,std::vector<ZZ> &B1R11,std::vector<ZZ> &B1R12
					,std::vector<ZZ> &B1R13,std::vector<ZZ> &B1R14,std::vector<ZZ> &B1R15); 
	//----------------------------------------------------------------------------------------------
	//radix-16 and radix-4 Mixed radix NTT
	void NTT_r16_r4(std::vector<ZZ> &A,std::vector<ZZ> &B0R0,std::vector<ZZ> &B0R1,std::vector<ZZ> &B0R2,std::vector<ZZ> &B0R3
	                ,std::vector<ZZ> &B0R4,std::vector<ZZ> &B0R5,std::vector<ZZ> &B0R6,std::vector<ZZ> &B0R7,std::vector<ZZ> &B0R8
					,std::vector<ZZ> &B0R9,std::vector<ZZ> &B0R10,std::vector<ZZ> &B0R11,std::vector<ZZ> &B0R12,std::vector<ZZ> &B0R13
					,std::vector<ZZ> &B0R14,std::vector<ZZ> &B0R15,std::vector<ZZ> &B1R0,std::vector<ZZ> &B1R1,std::vector<ZZ> &B1R2
					,std::vector<ZZ> &B1R3,std::vector<ZZ> &B1R4,std::vector<ZZ> &B1R5,std::vector<ZZ> &B1R6,std::vector<ZZ> &B1R7
					,std::vector<ZZ> &B1R8,std::vector<ZZ> &B1R9,std::vector<ZZ> &B1R10,std::vector<ZZ> &B1R11,std::vector<ZZ> &B1R12
					,std::vector<ZZ> &B1R13,std::vector<ZZ> &B1R14,std::vector<ZZ> &B1R15);
    void INTT_r16_r4(std::vector<ZZ> &A,std::vector<ZZ> &B0R0,std::vector<ZZ> &B0R1,std::vector<ZZ> &B0R2,std::vector<ZZ> &B0R3
	                ,std::vector<ZZ> &B0R4,std::vector<ZZ> &B0R5,std::vector<ZZ> &B0R6,std::vector<ZZ> &B0R7,std::vector<ZZ> &B0R8
					,std::vector<ZZ> &B0R9,std::vector<ZZ> &B0R10,std::vector<ZZ> &B0R11,std::vector<ZZ> &B0R12,std::vector<ZZ> &B0R13
					,std::vector<ZZ> &B0R14,std::vector<ZZ> &B0R15,std::vector<ZZ> &B1R0,std::vector<ZZ> &B1R1,std::vector<ZZ> &B1R2
					,std::vector<ZZ> &B1R3,std::vector<ZZ> &B1R4,std::vector<ZZ> &B1R5,std::vector<ZZ> &B1R6,std::vector<ZZ> &B1R7
					,std::vector<ZZ> &B1R8,std::vector<ZZ> &B1R9,std::vector<ZZ> &B1R10,std::vector<ZZ> &B1R11,std::vector<ZZ> &B1R12
					,std::vector<ZZ> &B1R13,std::vector<ZZ> &B1R14,std::vector<ZZ> &B1R15);	
	//----------------------------------------------------------------------------------------------
	void NTT_r16_r8(std::vector<ZZ> &A,std::vector<ZZ> &B0R0,std::vector<ZZ> &B0R1,std::vector<ZZ> &B0R2,std::vector<ZZ> &B0R3
	                ,std::vector<ZZ> &B0R4,std::vector<ZZ> &B0R5,std::vector<ZZ> &B0R6,std::vector<ZZ> &B0R7,std::vector<ZZ> &B0R8
					,std::vector<ZZ> &B0R9,std::vector<ZZ> &B0R10,std::vector<ZZ> &B0R11,std::vector<ZZ> &B0R12,std::vector<ZZ> &B0R13
					,std::vector<ZZ> &B0R14,std::vector<ZZ> &B0R15,std::vector<ZZ> &B1R0,std::vector<ZZ> &B1R1,std::vector<ZZ> &B1R2
					,std::vector<ZZ> &B1R3,std::vector<ZZ> &B1R4,std::vector<ZZ> &B1R5,std::vector<ZZ> &B1R6,std::vector<ZZ> &B1R7
					,std::vector<ZZ> &B1R8,std::vector<ZZ> &B1R9,std::vector<ZZ> &B1R10,std::vector<ZZ> &B1R11,std::vector<ZZ> &B1R12
					,std::vector<ZZ> &B1R13,std::vector<ZZ> &B1R14,std::vector<ZZ> &B1R15);
    void INTT_r16_r8(std::vector<ZZ> &A,std::vector<ZZ> &B0R0,std::vector<ZZ> &B0R1,std::vector<ZZ> &B0R2,std::vector<ZZ> &B0R3
	                ,std::vector<ZZ> &B0R4,std::vector<ZZ> &B0R5,std::vector<ZZ> &B0R6,std::vector<ZZ> &B0R7,std::vector<ZZ> &B0R8
					,std::vector<ZZ> &B0R9,std::vector<ZZ> &B0R10,std::vector<ZZ> &B0R11,std::vector<ZZ> &B0R12,std::vector<ZZ> &B0R13
					,std::vector<ZZ> &B0R14,std::vector<ZZ> &B0R15,std::vector<ZZ> &B1R0,std::vector<ZZ> &B1R1,std::vector<ZZ> &B1R2
					,std::vector<ZZ> &B1R3,std::vector<ZZ> &B1R4,std::vector<ZZ> &B1R5,std::vector<ZZ> &B1R6,std::vector<ZZ> &B1R7
					,std::vector<ZZ> &B1R8,std::vector<ZZ> &B1R9,std::vector<ZZ> &B1R10,std::vector<ZZ> &B1R11,std::vector<ZZ> &B1R12
					,std::vector<ZZ> &B1R13,std::vector<ZZ> &B1R14,std::vector<ZZ> &B1R15);
	//----------------------------------------------------------------------------------------------
	//radix-2
	void RR_R2(int BC,int shift_bit,int &result); //circulant right shift
	void AGU_R2(int BC,int &BN,int &MA); //radix-2 address generate unit
	void BR_R2(int BC,int &result); //bit reverse radix-2
	//radix-4
	void RR_R4(int BC,int shift_bit,int &result); //circulant right shift 
	void AGU_R4(int BC,int &BN,int &MA); //radix-4 address generate unit 
	void RR_R4_R2(int BC,int shift_bit,int &result); //circulant right shift for R4_R2 NTT	
	void BR_R4(int BC,int &result); //bit reverse radix-4
	void BR_R4_R2(int BC,int &result);
	void REORDERBC_R4_R2_OUT(int BC,int &result); // Reorder bc for output index
	void BC_IFFT_Reorder_R4_R2(int BC,int &result);
	void INTT_REORDERBC_R4_R2_OUT(int BC,int &result); // INTT Reorder bc for output index
	//radix-16
	void RR_R16(int BC,int shift_bit,int &result); //circulant right shift 
	void AGU_R16(int BC,int &BN,int &MA); //radix-16 address generate unit 
	void BR_R16(int BC,int &result); //bit reverse radix-16
	void RR_R16_R2(int BC,int shift_bit,int &result); //circulant right shift for R16_R2 NTT	
	void RR_R16_R4(int BC,int shift_bit,int &result); //circulant right shift for R16_R4 NTT	
	void RR_R16_R8(int BC,int shift_bit,int &result); //circulant right shift for R16_R8 NTT	
	// Mixed radix reindex for NTT output
	void NTT_REORDERINDEX_R16_R2_OUT(int BC,int &result); // Reorder INDEX for output index
	void NTT_REORDERINDEX_R16_R4_OUT(int BC,int &result); // Reorder INDEX for output index
	void NTT_REORDERINDEX_R16_R8_OUT(int BC,int &result); // Reorder INDEX for output index
	//Mixed radix-2
	void BC_IFFT_Reorder_R16_R2(int BC,int &result);
	void BC_IFFT_Reorder_R16_R2_OUT(int BC,int &result);
	//Mixed radix-4
	void BC_IFFT_Reorder_R16_R4(int BC,int &result);
	void BC_IFFT_Reorder_R16_R4_OUT(int BC,int &result);
	//Mixed radix-8
	void BC_IFFT_Reorder_R16_R8(int BC,int &result);
	void BC_IFFT_Reorder_R16_R8_OUT(int BC,int &result);
	
	//Gray code transform
	int  Gray(int index,int group);
};

#endif
