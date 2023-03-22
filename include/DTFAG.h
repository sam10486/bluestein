#ifndef _DTFAG_H_
#define _DTFAG_H_

#include <iostream>
#include "NTT.h"
#include "NTTSPMB.h"
#include "BitOperate.h"

using namespace std;

class DTFAG : public BitOperate {
public:
    void DTFAG_DIT();
    void DTFAG_DIF();
    void DTFAG_test();
    void DTFAG_DIF_MixedRadix();
    void DTFAG_verify();
};

#endif