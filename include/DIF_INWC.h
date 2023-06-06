#ifndef _DIF_INWC_H_
#define _DIF_INWC_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <NTL/ZZ.h>
#include "NTTSPMB.h"
#include "BitOperate.h"

using namespace NTL;

class DIF_INWC : public NTTSPMB, public BitOperate{
public:
    void INWC_Radix2_BU(ZZ &a,ZZ &b);
    void INWC_Radix4_BU(ZZ &a,ZZ &b,ZZ &c,ZZ &d);
    void INWC_Radix8_BU(ZZ &a_r0,ZZ &a_r1,ZZ &a_r2,ZZ &a_r3,ZZ &a_r4,ZZ &a_r5,ZZ &a_r6,ZZ &a_r7);
    void INWC_Radix16_BU(   ZZ &a_r0,ZZ &a_r1,ZZ &a_r2,ZZ &a_r3,ZZ &a_r4,
                            ZZ &a_r5,ZZ &a_r6,ZZ &a_r7,
                            ZZ &a_r8,ZZ &a_r9,ZZ &a_r10,ZZ &a_r11,
                            ZZ &a_r12,ZZ &a_r13,ZZ &a_r14,ZZ &a_r15);
    void DIF_INWC_radix2(std::vector<ZZ> &A); //radix-2 NTT
    void DIF_INWC_radix4(std::vector<ZZ> &A);
    void DIF_INWC_radix16(std::vector<ZZ> &A);
    void DIF_INWC_r4_r2(std::vector<ZZ> &A,
	std::vector<ZZ> &B0R0,std::vector<ZZ> &B0R1,std::vector<ZZ> &B0R2,std::vector<ZZ> &B0R3,
	std::vector<ZZ> &B1R0,std::vector<ZZ> &B1R1,std::vector<ZZ> &B1R2,std::vector<ZZ> &B1R3);
    void DIF_INWC_r16_r2(std::vector<ZZ> &A,std::vector<ZZ> &B0R0,std::vector<ZZ> &B0R1,std::vector<ZZ> &B0R2,std::vector<ZZ> &B0R3,
	std::vector<ZZ> &B0R4,std::vector<ZZ> &B0R5,std::vector<ZZ> &B0R6,std::vector<ZZ> &B0R7,std::vector<ZZ> &B0R8,std::vector<ZZ> &B0R9,
	std::vector<ZZ> &B0R10,std::vector<ZZ> &B0R11,std::vector<ZZ> &B0R12,std::vector<ZZ> &B0R13,std::vector<ZZ> &B0R14,std::vector<ZZ> &B0R15,
	std::vector<ZZ> &B1R0,std::vector<ZZ> &B1R1,std::vector<ZZ> &B1R2,std::vector<ZZ> &B1R3,std::vector<ZZ> &B1R4,std::vector<ZZ> &B1R5,
	std::vector<ZZ> &B1R6,std::vector<ZZ> &B1R7,std::vector<ZZ> &B1R8,std::vector<ZZ> &B1R9,std::vector<ZZ> &B1R10,std::vector<ZZ> &B1R11,
	std::vector<ZZ> &B1R12,std::vector<ZZ> &B1R13,std::vector<ZZ> &B1R14,std::vector<ZZ> &B1R15);
    void DIF_INWC_r16_r4(std::vector<ZZ> &A,std::vector<ZZ> &B0R0,std::vector<ZZ> &B0R1,std::vector<ZZ> &B0R2,std::vector<ZZ> &B0R3,
	std::vector<ZZ> &B0R4,std::vector<ZZ> &B0R5,std::vector<ZZ> &B0R6,std::vector<ZZ> &B0R7,std::vector<ZZ> &B0R8,std::vector<ZZ> &B0R9,
	std::vector<ZZ> &B0R10,std::vector<ZZ> &B0R11,std::vector<ZZ> &B0R12,std::vector<ZZ> &B0R13,std::vector<ZZ> &B0R14,std::vector<ZZ> &B0R15,
	std::vector<ZZ> &B1R0,std::vector<ZZ> &B1R1,std::vector<ZZ> &B1R2,std::vector<ZZ> &B1R3,std::vector<ZZ> &B1R4,std::vector<ZZ> &B1R5,
	std::vector<ZZ> &B1R6,std::vector<ZZ> &B1R7,std::vector<ZZ> &B1R8,std::vector<ZZ> &B1R9,std::vector<ZZ> &B1R10,std::vector<ZZ> &B1R11,
	std::vector<ZZ> &B1R12,std::vector<ZZ> &B1R13,std::vector<ZZ> &B1R14,std::vector<ZZ> &B1R15);
    void DIF_INWC_r16_r8(std::vector<ZZ> &A,std::vector<ZZ> &B0R0,std::vector<ZZ> &B0R1,std::vector<ZZ> &B0R2,std::vector<ZZ> &B0R3,
	std::vector<ZZ> &B0R4,std::vector<ZZ> &B0R5,std::vector<ZZ> &B0R6,std::vector<ZZ> &B0R7,std::vector<ZZ> &B0R8,std::vector<ZZ> &B0R9,
	std::vector<ZZ> &B0R10,std::vector<ZZ> &B0R11,std::vector<ZZ> &B0R12,std::vector<ZZ> &B0R13,std::vector<ZZ> &B0R14,std::vector<ZZ> &B0R15,
	std::vector<ZZ> &B1R0,std::vector<ZZ> &B1R1,std::vector<ZZ> &B1R2,std::vector<ZZ> &B1R3,std::vector<ZZ> &B1R4,std::vector<ZZ> &B1R5,
	std::vector<ZZ> &B1R6,std::vector<ZZ> &B1R7,std::vector<ZZ> &B1R8,std::vector<ZZ> &B1R9,std::vector<ZZ> &B1R10,std::vector<ZZ> &B1R11,
	std::vector<ZZ> &B1R12,std::vector<ZZ> &B1R13,std::vector<ZZ> &B1R14,std::vector<ZZ> &B1R15);

	//----------INWC MergeFacor--------------
	void INWC_MergeFactor_Radix2_BU(ZZ &a,ZZ &b);
	void INWC_MergeFactor_Radix4_BU(ZZ &a,ZZ &b,ZZ &c,ZZ &d);
	void INWC_MergeFactor_Radix8_BU(ZZ &a_r0,ZZ &a_r1,ZZ &a_r2,ZZ &a_r3,ZZ &a_r4,ZZ &a_r5,ZZ &a_r6,ZZ &a_r7);
	void INWC_MergeFactor_Radix16_BU(ZZ &a_r0,ZZ &a_r1,ZZ &a_r2,ZZ &a_r3,ZZ &a_r4,ZZ &a_r5,ZZ &a_r6,ZZ &a_r7,
		ZZ &a_r8,ZZ &a_r9,ZZ &a_r10,ZZ &a_r11,ZZ &a_r12,ZZ &a_r13,ZZ &a_r14,ZZ &a_r15);

	void DIF_INWC_MergeFactor_radix2(std::vector<ZZ> &A); //radix-2 NTT
	void DIF_INWC_MergeFactor_radix4(std::vector<ZZ> &A);
	void DIF_INWC_MergeFactor_radix16(std::vector<ZZ> &A);
	void DIF_INWC_MergeFactor_r16_r2(std::vector<ZZ> &A,std::vector<ZZ> &B0R0,std::vector<ZZ> &B0R1,std::vector<ZZ> &B0R2,std::vector<ZZ> &B0R3,
	std::vector<ZZ> &B0R4,std::vector<ZZ> &B0R5,std::vector<ZZ> &B0R6,std::vector<ZZ> &B0R7,std::vector<ZZ> &B0R8,std::vector<ZZ> &B0R9,
	std::vector<ZZ> &B0R10,std::vector<ZZ> &B0R11,std::vector<ZZ> &B0R12,std::vector<ZZ> &B0R13,std::vector<ZZ> &B0R14,std::vector<ZZ> &B0R15,
	std::vector<ZZ> &B1R0,std::vector<ZZ> &B1R1,std::vector<ZZ> &B1R2,std::vector<ZZ> &B1R3,std::vector<ZZ> &B1R4,std::vector<ZZ> &B1R5,
	std::vector<ZZ> &B1R6,std::vector<ZZ> &B1R7,std::vector<ZZ> &B1R8,std::vector<ZZ> &B1R9,std::vector<ZZ> &B1R10,std::vector<ZZ> &B1R11,
	std::vector<ZZ> &B1R12,std::vector<ZZ> &B1R13,std::vector<ZZ> &B1R14,std::vector<ZZ> &B1R15);
	void DIF_INWC_MergeFactor_r16_r4(std::vector<ZZ> &A,std::vector<ZZ> &B0R0,std::vector<ZZ> &B0R1,std::vector<ZZ> &B0R2,std::vector<ZZ> &B0R3,
	std::vector<ZZ> &B0R4,std::vector<ZZ> &B0R5,std::vector<ZZ> &B0R6,std::vector<ZZ> &B0R7,std::vector<ZZ> &B0R8,std::vector<ZZ> &B0R9,
	std::vector<ZZ> &B0R10,std::vector<ZZ> &B0R11,std::vector<ZZ> &B0R12,std::vector<ZZ> &B0R13,std::vector<ZZ> &B0R14,std::vector<ZZ> &B0R15,
	std::vector<ZZ> &B1R0,std::vector<ZZ> &B1R1,std::vector<ZZ> &B1R2,std::vector<ZZ> &B1R3,std::vector<ZZ> &B1R4,std::vector<ZZ> &B1R5,
	std::vector<ZZ> &B1R6,std::vector<ZZ> &B1R7,std::vector<ZZ> &B1R8,std::vector<ZZ> &B1R9,std::vector<ZZ> &B1R10,std::vector<ZZ> &B1R11,
	std::vector<ZZ> &B1R12,std::vector<ZZ> &B1R13,std::vector<ZZ> &B1R14,std::vector<ZZ> &B1R15);
	void DIF_INWC_MergeFactor_r16_r8(std::vector<ZZ> &A,std::vector<ZZ> &B0R0,std::vector<ZZ> &B0R1,std::vector<ZZ> &B0R2,std::vector<ZZ> &B0R3,
	std::vector<ZZ> &B0R4,std::vector<ZZ> &B0R5,std::vector<ZZ> &B0R6,std::vector<ZZ> &B0R7,std::vector<ZZ> &B0R8,std::vector<ZZ> &B0R9,
	std::vector<ZZ> &B0R10,std::vector<ZZ> &B0R11,std::vector<ZZ> &B0R12,std::vector<ZZ> &B0R13,std::vector<ZZ> &B0R14,std::vector<ZZ> &B0R15,
	std::vector<ZZ> &B1R0,std::vector<ZZ> &B1R1,std::vector<ZZ> &B1R2,std::vector<ZZ> &B1R3,std::vector<ZZ> &B1R4,std::vector<ZZ> &B1R5,
	std::vector<ZZ> &B1R6,std::vector<ZZ> &B1R7,std::vector<ZZ> &B1R8,std::vector<ZZ> &B1R9,std::vector<ZZ> &B1R10,std::vector<ZZ> &B1R11,
	std::vector<ZZ> &B1R12,std::vector<ZZ> &B1R13,std::vector<ZZ> &B1R14,std::vector<ZZ> &B1R15);
	//---------------------------------------

	//---------------INWC seperate N^-1-----------------
	void INWC_seperateInvN_Radix2_BU(ZZ &a,ZZ &b, ZZ InvTwo);
	void INWC_seperateInvN_Radix4_BU(ZZ &a,ZZ &b,ZZ &c,ZZ &d, ZZ InvTwo);
	void INWC_seperateInvN_Radix8_BU(ZZ &a_r0,ZZ &a_r1,ZZ &a_r2,ZZ &a_r3,ZZ &a_r4,ZZ &a_r5,ZZ &a_r6,ZZ &a_r7, ZZ InvTwo);
	void INWC_seperateInvN_Radix16_BU(ZZ &a_r0,ZZ &a_r1,ZZ &a_r2,ZZ &a_r3,ZZ &a_r4,ZZ &a_r5,ZZ &a_r6,ZZ &a_r7,
	ZZ &a_r8,ZZ &a_r9,ZZ &a_r10,ZZ &a_r11,ZZ &a_r12,ZZ &a_r13,ZZ &a_r14,ZZ &a_r15, ZZ InvTwo);
	void DIF_INWC_seperateInvN_radix2(std::vector<ZZ> &A);
	void DIF_INWC_seperateInvN_radix4(std::vector<ZZ> &A);
	void DIF_INWC_seperateInvN_radix16(std::vector<ZZ> &A);

	void DIF_INWC_seperateInvN_r16_r2(std::vector<ZZ> &A,std::vector<ZZ> &B0R0,std::vector<ZZ> &B0R1,std::vector<ZZ> &B0R2,std::vector<ZZ> &B0R3,
	std::vector<ZZ> &B0R4,std::vector<ZZ> &B0R5,std::vector<ZZ> &B0R6,std::vector<ZZ> &B0R7,std::vector<ZZ> &B0R8,std::vector<ZZ> &B0R9,
	std::vector<ZZ> &B0R10,std::vector<ZZ> &B0R11,std::vector<ZZ> &B0R12,std::vector<ZZ> &B0R13,std::vector<ZZ> &B0R14,std::vector<ZZ> &B0R15,
	std::vector<ZZ> &B1R0,std::vector<ZZ> &B1R1,std::vector<ZZ> &B1R2,std::vector<ZZ> &B1R3,std::vector<ZZ> &B1R4,std::vector<ZZ> &B1R5,
	std::vector<ZZ> &B1R6,std::vector<ZZ> &B1R7,std::vector<ZZ> &B1R8,std::vector<ZZ> &B1R9,std::vector<ZZ> &B1R10,std::vector<ZZ> &B1R11,
	std::vector<ZZ> &B1R12,std::vector<ZZ> &B1R13,std::vector<ZZ> &B1R14,std::vector<ZZ> &B1R15);
	void DIF_INWC_seperateInvN_r16_r4(std::vector<ZZ> &A,std::vector<ZZ> &B0R0,std::vector<ZZ> &B0R1,std::vector<ZZ> &B0R2,std::vector<ZZ> &B0R3,
	std::vector<ZZ> &B0R4,std::vector<ZZ> &B0R5,std::vector<ZZ> &B0R6,std::vector<ZZ> &B0R7,std::vector<ZZ> &B0R8,std::vector<ZZ> &B0R9,
	std::vector<ZZ> &B0R10,std::vector<ZZ> &B0R11,std::vector<ZZ> &B0R12,std::vector<ZZ> &B0R13,std::vector<ZZ> &B0R14,std::vector<ZZ> &B0R15,
	std::vector<ZZ> &B1R0,std::vector<ZZ> &B1R1,std::vector<ZZ> &B1R2,std::vector<ZZ> &B1R3,std::vector<ZZ> &B1R4,std::vector<ZZ> &B1R5,
	std::vector<ZZ> &B1R6,std::vector<ZZ> &B1R7,std::vector<ZZ> &B1R8,std::vector<ZZ> &B1R9,std::vector<ZZ> &B1R10,std::vector<ZZ> &B1R11,
	std::vector<ZZ> &B1R12,std::vector<ZZ> &B1R13,std::vector<ZZ> &B1R14,std::vector<ZZ> &B1R15);
	void DIF_INWC_seperateInvN_r16_r8(std::vector<ZZ> &A,std::vector<ZZ> &B0R0,std::vector<ZZ> &B0R1,std::vector<ZZ> &B0R2,std::vector<ZZ> &B0R3,
	std::vector<ZZ> &B0R4,std::vector<ZZ> &B0R5,std::vector<ZZ> &B0R6,std::vector<ZZ> &B0R7,std::vector<ZZ> &B0R8,std::vector<ZZ> &B0R9,
	std::vector<ZZ> &B0R10,std::vector<ZZ> &B0R11,std::vector<ZZ> &B0R12,std::vector<ZZ> &B0R13,std::vector<ZZ> &B0R14,std::vector<ZZ> &B0R15,
	std::vector<ZZ> &B1R0,std::vector<ZZ> &B1R1,std::vector<ZZ> &B1R2,std::vector<ZZ> &B1R3,std::vector<ZZ> &B1R4,std::vector<ZZ> &B1R5,
	std::vector<ZZ> &B1R6,std::vector<ZZ> &B1R7,std::vector<ZZ> &B1R8,std::vector<ZZ> &B1R9,std::vector<ZZ> &B1R10,std::vector<ZZ> &B1R11,
	std::vector<ZZ> &B1R12,std::vector<ZZ> &B1R13,std::vector<ZZ> &B1R14,std::vector<ZZ> &B1R15);
	//--------------------------------------------------
};
#endif