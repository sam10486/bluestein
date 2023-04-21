#ifndef _DIT_NTTSPMB_H_
#define _DIT_NTTSPMB_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <NTL/ZZ.h>
#include "NTTSPMB.h"
#include "BitOperate.h"

using namespace NTL;

class DIT_NTTSPMB : public NTTSPMB, public BitOperate{
public:
    void DIT_NTT_radix2(std::vector<ZZ> &A); //radix-2 NTT
    void DIT_NTT_radix4(std::vector<ZZ> &A); //radix-4 NTT
    void DIT_NTT_radix16(std::vector<ZZ> &A);

    void DIT_NTT_r4_r2(
        vector<ZZ> &A,
	    vector<ZZ> &B0R0,   vector<ZZ> &B0R1,   vector<ZZ> &B0R2,   vector<ZZ> &B0R3,
	    vector<ZZ> &B1R0,   vector<ZZ> &B1R1,   vector<ZZ> &B1R2,   vector<ZZ> &B1R3);
    
    void test_radix2(std::vector<ZZ> &A);
    void test_radix4(std::vector<ZZ> &A);
};
#endif