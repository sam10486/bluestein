#include <iostream>
#include <NTL/ZZ.h>
#include <vector>
#include "forward_NWC.h"
#include "BitOperate.h"
#include "assert.h"

using namespace NTL;
using namespace std;

vector<ZZ > forward_NWC::NWC(vector<ZZ > arr, long long N, ZZ phi, ZZ Modular){
    long long t = N;
    BitOperate Bitrev;
    for(long long m = 1; m < N; m = m << 1){
        t = t / 2;
        for(long long i = 0; i < m; i++){
            long long j1 = 2*i*t;
            long long j2 = j1 + t - 1;
            long long bit_num = ceil(log2(N));
            long long index = Bitrev.BitReserve((m+i), bit_num);
            ZZ S;
            PowerMod(S, phi, index, Modular);
            for(long long j = j1; j <= j2; j++){
                ZZ U, V;
                U = arr.at(j);
                V = MulMod(arr.at(j+t), S, Modular);
                arr.at(j) = AddMod(U, V, Modular);
                arr.at(j+t) = SubMod(U, V, Modular);
            }
        }
    }
    return arr;
}