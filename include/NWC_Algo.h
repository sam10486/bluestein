#ifndef _NWC_Algo_H_
#define _NWC_Algo_H_

#include <iostream>
#include <NTL/ZZ.h>
#include <vector>
#include "NWC_util.h"

using namespace NTL;
using namespace std;

class NWC_Algo: public NWC_util{
private:
    /* data */
public:
    NWC_Algo(long long Radix_r1, long long Radix_r2, long long N, ZZ Modular);
    ~NWC_Algo();
    vector<ZZ > NWC(vector<ZZ > &arr);
    vector<ZZ > INWC(vector<ZZ > &arr);
};





#endif