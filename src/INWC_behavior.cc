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
    long long N = 4096;
    ZZ Modular;
    ZZ Phi, InvPhi, W, IW;
    conv(Modular, "18446744069414584321");
    //conv(Phi, "")

    NWC_util NWC(Radix_r1, Radix_r2, N, Modular);
    NWC.showInfo();

    cout << "---------------" << endl;
    vector<ZZ > arr, golden_arr;
    arr.resize(N);
    golden_arr.resize(N);
    NWC_Algo nwc_algo(Radix_r1, Radix_r2, N, Modular);
    //nwc_algo.setValue(Radix_r1, Radix_r2, N, Modular, )
    for (int i = 0; i < N; i++){
        golden_arr[i] = i;
        arr[i] = i;
    }    
    nwc_algo.NWC(arr);
    //for (int i = 0; i < N; i++){
    //    cout << "arr_NWC[" << i << "] = " << arr[i] << endl;
    //}
    cout << "---------------" << endl;
    nwc_algo.INWC(arr);
    //for (int i = 0; i < N; i++){
    //    cout << "arr_INWC[" << i << "] = " << arr[i] << endl;
    //}

    int err = 0;
    for (int i = 0; i < N; i++){
        if (arr[i] != golden_arr[i]){
            err++;
            cout << "arr[" << i << "] = " << arr[i] << " != " << golden_arr[i] << endl;
        }
    }
    if (err == 0){
        cout << "ALL pass!" << endl;
    }else{
        cout << "Something error..." << endl;
    }
    
    

    
    
    return 0;
}