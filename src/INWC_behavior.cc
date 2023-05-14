#include <iostream>
#include <NTL/ZZ.h>
#include <vector>
#include "NWC_util.h"
#include "assert.h"
#include "NWC_util.h"

using namespace std;
using namespace NTL;

int INWC_behavior(){
    long long Radix_r1 = 2;
    long long Radix_r2 = 2;
    long long N = 16;
    ZZ Modular;
    conv(Modular, "97");

    NWC_util NWC(Radix_r1, Radix_r2, N, Modular);
    NWC.showInfo();
    return 0;
}