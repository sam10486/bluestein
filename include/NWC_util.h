#ifndef _NWC_util_H_
#define _NWC_util_H_

#include <iostream>
#include <NTL/ZZ.h>
#include <vector>

using namespace NTL;

class NWC_util{

private:
    ZZ Radix_r1;
    ZZ Radix_r2;
    ZZ N;           // polynomial degree
    ZZ Modular;     // polynomial modulus
    ZZ W;           // polynomial twiddle factor
    ZZ IW;          // polynomial inverse twiddle factor

public:
    NWC_util();
    ~NWC_util();
};


#endif
