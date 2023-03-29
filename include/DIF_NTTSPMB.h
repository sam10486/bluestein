#ifndef _DIF_NTTSPMB_H_
#define _DIF_NTTSPMB_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <NTL/ZZ.h>
#include "NTTSPMB.h"
#include "BitOperate.h"

using namespace NTL;

class DIF_NTTSPMB : public NTTSPMB, public BitOperate{
public:
    void DIF_NTT_radix2(std::vector<ZZ> &A); //radix-2 NTT
    void DIF_NTT_radix4(std::vector<ZZ> &A); //radix-4 NTT
    void DIF_NTT_radix16(std::vector<ZZ> &A);
    void DIF_NTT_r4_r2(std::vector<ZZ> &A,
        std::vector<ZZ> &B0R0,std::vector<ZZ> &B0R1,std::vector<ZZ> &B0R2,std::vector<ZZ> &B0R3,
        std::vector<ZZ> &B1R0,std::vector<ZZ> &B1R1,std::vector<ZZ> &B1R2,std::vector<ZZ> &B1R3);

    void DIF_NTT_r16_r2(
        vector<ZZ> &A,      vector<ZZ> &B0R0,   vector<ZZ> &B0R1,   vector<ZZ> &B0R2,   vector<ZZ> &B0R3,
	    vector<ZZ> &B0R4,   vector<ZZ> &B0R5,   vector<ZZ> &B0R6,   vector<ZZ> &B0R7,   vector<ZZ> &B0R8,std::vector<ZZ> &B0R9,
	    vector<ZZ> &B0R10,  vector<ZZ> &B0R11,  vector<ZZ> &B0R12,  vector<ZZ> &B0R13,  vector<ZZ> &B0R14,std::vector<ZZ> &B0R15,
	    vector<ZZ> &B1R0,   vector<ZZ> &B1R1,   vector<ZZ> &B1R2,   vector<ZZ> &B1R3,   vector<ZZ> &B1R4,std::vector<ZZ> &B1R5,
	    vector<ZZ> &B1R6,   vector<ZZ> &B1R7,   vector<ZZ> &B1R8,   vector<ZZ> &B1R9,   vector<ZZ> &B1R10,std::vector<ZZ> &B1R11,
	    vector<ZZ> &B1R12,  vector<ZZ> &B1R13,  vector<ZZ> &B1R14,  vector<ZZ> &B1R15);

    void DIF_NTT_r16_r4(
        vector<ZZ> &A,      vector<ZZ> &B0R0,   vector<ZZ> &B0R1,   vector<ZZ> &B0R2,   vector<ZZ> &B0R3,
	    vector<ZZ> &B0R4,   vector<ZZ> &B0R5,   vector<ZZ> &B0R6,   vector<ZZ> &B0R7,   vector<ZZ> &B0R8,   vector<ZZ> &B0R9,
	    vector<ZZ> &B0R10,  vector<ZZ> &B0R11,  vector<ZZ> &B0R12,  vector<ZZ> &B0R13,  vector<ZZ> &B0R14,  vector<ZZ> &B0R15,
	    vector<ZZ> &B1R0,   vector<ZZ> &B1R1,   vector<ZZ> &B1R2,   vector<ZZ> &B1R3,   vector<ZZ> &B1R4,   vector<ZZ> &B1R5,
	    vector<ZZ> &B1R6,   vector<ZZ> &B1R7,   vector<ZZ> &B1R8,   vector<ZZ> &B1R9,   vector<ZZ> &B1R10,  vector<ZZ> &B1R11,
	    vector<ZZ> &B1R12,  vector<ZZ> &B1R13,  vector<ZZ> &B1R14,  vector<ZZ> &B1R15);
    
    void DIF_NTT_r16_r8(
        vector<ZZ> &A,
        vector<ZZ> &B0R0,   vector<ZZ> &B0R1,   vector<ZZ> &B0R2,   vector<ZZ> &B0R3,
	    vector<ZZ> &B0R4,   vector<ZZ> &B0R5,   vector<ZZ> &B0R6,   vector<ZZ> &B0R7,   vector<ZZ> &B0R8,   vector<ZZ> &B0R9,
	    vector<ZZ> &B0R10,  vector<ZZ> &B0R11,  vector<ZZ> &B0R12,  vector<ZZ> &B0R13,  vector<ZZ> &B0R14,  vector<ZZ> &B0R15,
	    vector<ZZ> &B1R0,   vector<ZZ> &B1R1,   vector<ZZ> &B1R2,   vector<ZZ> &B1R3,   vector<ZZ> &B1R4,   vector<ZZ> &B1R5,
	    vector<ZZ> &B1R6,   vector<ZZ> &B1R7,   vector<ZZ> &B1R8,   vector<ZZ> &B1R9,   vector<ZZ> &B1R10,  vector<ZZ> &B1R11,
	    vector<ZZ> &B1R12,  vector<ZZ> &B1R13,  vector<ZZ> &B1R14,  vector<ZZ> &B1R15);
};
#endif