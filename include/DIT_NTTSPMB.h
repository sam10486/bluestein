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
};
#endif