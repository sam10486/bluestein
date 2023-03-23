#ifndef _DTFAG_H_
#define _DTFAG_H_

#include <iostream>
#include "NTT.h"
#include "NTTSPMB.h"
#include "BitOperate.h"

using namespace std;

class DTFAG : public BitOperate {
public:
    void DTFAG_SPMB_DIT( int stage, 
                    vector<ZZ > &st0_Tw, vector<ZZ > &st1_Tw, vector<ZZ > &st2_Tw, 
                    int DTFAG_t, int DTFAG_i, int DTFAG_j);
    void DTFAG_DIF();
    void DTFAG_DIT();
    void DTFAG_DIF_MixedRadix();
    void DTFAG_verify();
};

#endif