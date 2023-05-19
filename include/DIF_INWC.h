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
    void INWC_Radix16_BU(   ZZ &a_r0,ZZ &a_r1,ZZ &a_r2,ZZ &a_r3,ZZ &a_r4,
                            ZZ &a_r5,ZZ &a_r6,ZZ &a_r7,
                            ZZ &a_r8,ZZ &a_r9,ZZ &a_r10,ZZ &a_r11,
                            ZZ &a_r12,ZZ &a_r13,ZZ &a_r14,ZZ &a_r15);
    void DIF_INWC_radix2(std::vector<ZZ> &A); //radix-2 NTT
    void DIF_INWC_radix4(std::vector<ZZ> &A);
    void DIF_INWC_radix16(std::vector<ZZ> &A);
};
#endif