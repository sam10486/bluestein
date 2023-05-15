#ifndef _NWC_util_H_
#define _NWC_util_H_

#include <iostream>
#include <NTL/ZZ.h>
#include <vector>

using namespace NTL;

class NWC_util{

private:
    long long Radix_r1;
    long long Radix_r2;
    long long N;           // polynomial degree
    ZZ Modular;     // polynomial modulus
    ZZ Phi;           // polynomial twiddle factor
    ZZ InvPhi;          // polynomial inverse twiddle factor

public:
    NWC_util(long long Radix_r1, long long Radix_r2, long long N, ZZ Modular);
    ~NWC_util();
    void setValue(long long Radix_r1, long long Radix_r2, long long N, ZZ Modular, ZZ Phi, ZZ InvPhi);
    void getValue(long long *Radix_r1, long long *Radix_r2, long long *N, ZZ *Modular, ZZ *Phi, ZZ *InvPhi);
    ZZ find_n_rou(ZZ base, long long N, ZZ Modular);
    ZZ find_phi(long long N, ZZ Modular);
    bool check_prou(ZZ n_rou, long long N, ZZ Modular);
    void showInfo();
};


#endif
