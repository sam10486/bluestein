#ifndef _forward_NWC_H_
#define _forward_NWC_H_

#include <iostream>
#include <NTL/ZZ.h>
#include <vector>
#include "NWC_util.h"

using namespace NTL;
using namespace std;

class forward_NWC : public NWC_util{

public:
    vector<ZZ > NWC(vector<ZZ > arr, long long N, ZZ phi, ZZ Modular);

};


#endif