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
    void DIF_INWC_radix2(std::vector<ZZ> &A); //radix-2 NTT
    void DIF_INWC_radix4(std::vector<ZZ> &A);
};
#endif