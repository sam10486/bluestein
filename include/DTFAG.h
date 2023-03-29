#ifndef _DTFAG_H_
#define _DTFAG_H_

#include <iostream>
#include "NTT.h"
#include "NTTSPMB.h"
#include "BitOperate.h"

using namespace std;

class DTFAG : public BitOperate {
public:
    void DTFAG_SPMB_DIT( 
                    int stage, int fft_point, int radix_r1, int radix_r2, int debug,
                    vector<vector<ZZ > > &ROM0,  vector<vector<ZZ > > &ROM1,  vector<vector<ZZ > > &ROM2,
                    vector<ZZ > &st0_Tw, vector<ZZ > &st1_Tw, vector<ZZ > &st2_Tw, 
                    int DTFAG_t, int DTFAG_i, int DTFAG_j);
    void DTFAG_DIF();
    void DTFAG_DIT();
    void DTFAG_DIF_MixedRadix();
    void DTFAG_verify();
    void DTFAG_ROM_init(
    int radix_r1, int radix_r2, ZZ fft_twiddle, ZZ fft_prime, int debug,
    vector<vector<ZZ > > &ROM0,  vector<vector<ZZ > > &ROM1,  vector<vector<ZZ > > &ROM2);

    void DTFAG_SPMB_DIF_MR (
        int stage, int fft_point, int radix_r1, int radix_r2, int debug,
        vector<vector<ZZ > > &ROM0,  vector<vector<ZZ > > &ROM1,  vector<vector<ZZ > > &ROM2,
        vector<ZZ > &st0_Tw, vector<ZZ > &st1_Tw, vector<ZZ > &st2_Tw,
        int DTFAG_i, int DTFAG_t, int DTFAG_j);
};

#endif