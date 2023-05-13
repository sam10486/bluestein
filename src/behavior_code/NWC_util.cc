#include <iostream>
#include <NTL/ZZ.h>
#include <vector>
#include "NWC_util.h"


using namespace NTL;

NWC_util::NWC_util(){
    Radix_r1 = 0;
    Radix_r2 = 0;
    N = 0;
    Modular = 0;
    W = 0;
}

NWC_util::~NWC_util(){
}