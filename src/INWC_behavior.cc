#include <iostream>
#include <NTL/ZZ.h>
#include <vector>
#include "NWC_util.h"
#include "assert.h"
#include "NWC_util.h"
#include "NWC_Algo.h"

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

    cout << "---------------" << endl;
    vector<ZZ > arr;
    arr.resize(N);
    NWC_Algo nwc_algo(Radix_r1, Radix_r2, N, Modular);
    for (int i = 0; i < N; i++){
        arr[i] = i;
    }    
    nwc_algo.NWC(arr);
    for (int i = 0; i < N; i++){
        cout << "arr_NWC[" << i << "] = " << arr[i] << endl;
    }
    cout << "---------------" << endl;
    nwc_algo.INWC(arr);
    for (int i = 0; i < N; i++){
        cout << "arr_INWC[" << i << "] = " << arr[i] << endl;
    }
    
    return 0;
}