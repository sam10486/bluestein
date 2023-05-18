#include <iostream>
#include <NTL/ZZ.h>
#include <vector>
#include "NWC_Algo.h"
#include "BitOperate.h"
#include "assert.h"

using namespace NTL;
using namespace std;


NWC_Algo::NWC_Algo(long long Radix_r1, long long Radix_r2, long long N, ZZ Modular): NWC_util(Radix_r1, Radix_r2, N, Modular){
    //TODO
}

NWC_Algo::~NWC_Algo(){
    //TODO
}

vector<ZZ > NWC_Algo::NWC(vector<ZZ > &arr){
    long long Radix_r1, Radix_r2, N;
    ZZ Modular, Phi, InvPhi;
    ZZ W, IW;
    getValue(&Radix_r1, &Radix_r2, &N, &Modular , &Phi, &InvPhi, &W, &IW);
    BitOperate Bitrev;
    vector<ZZ> Arr_scramble;
    Arr_scramble.resize(N);
    for (int i = 0; i < N; i++){
        long long bit_num = ceil(log2(N));
        long long rev_index = Bitrev.BitReserve(i, bit_num);
        Arr_scramble[i] = arr[rev_index];
    }
    for (int s = 1; s <= log2(N); s++){
        //cout << "stage = " << s << endl;
        int m = pow(2, s);
        for (int j = 0; j <= (m/2)-1; j++){
            ZZ TF;
            int TF_deg = (2*j+1) * (N/m);
            TF = PowerMod(Phi, TF_deg, Modular);
            for (int k = 0; k <= (N/m)-1; k++){
                ZZ u, t;
                u = Arr_scramble[k*m+j];
                t = MulMod(TF, Arr_scramble[k*m+j+(m/2)], Modular);
                //cout << "u[" << k*m+j << "] = " << u << endl;
                //cout << "t[" << k*m+j+(m/2) << "] = " << t << ", TF = " << TF << endl;
                Arr_scramble[k*m+j] = AddMod(u, t, Modular);
                Arr_scramble[k*m+j+(m/2)] = SubMod(u, t, Modular);
                //cout << "Arr_scramble[" << k*m+j << "] = " << Arr_scramble[k*m+j] << endl;
                //cout << "Arr_scramble[" << k*m+j+(m/2) << "] = " << Arr_scramble[k*m+j+(m/2)] << endl;
            }
        }
    }
    for (int i = 0; i < N; i++){
        arr[i] = Arr_scramble[i];
    }
    
    return arr;
}

vector<ZZ > NWC_Algo::INWC(vector<ZZ> &arr){
    long long Radix_r1, Radix_r2, N;
    ZZ Modular, Phi, InvPhi;
    ZZ W, IW;
    getValue(&Radix_r1, &Radix_r2, &N, &Modular , &Phi, &InvPhi, &W, &IW);
    
    vector<ZZ > A_arr;
    A_arr.resize(N);
    for (int i = 0; i < N; i++){
        A_arr[i] = arr[i];
    }
    long long bit_num = ceil(log2(N));

    BitOperate Bitrev;
    ZZ InvTwo;
    InvMod(InvTwo, (ZZ)2, Modular);
    //cout << "InvTwo = " << InvTwo << endl;
    for (int s = log2(N); s >= 1; s--){
        //cout << "stage = " << s << endl;
        int m = pow(2, s);
        for (int j = 0; j <= (m/2)-1; j++){
            //ZZ TF;
            //int TF_deg = (2*j+1) * (N/m);
            //TF = PowerMod(InvPhi, TF_deg, Modular);
            int Phi_deg, W_deg;
            Phi_deg = (N/m);
            W_deg = (j) * (N/m);
            ZZ InvPhi_Order = PowerMod(InvPhi, Phi_deg, Modular);
            ZZ IW_Order = PowerMod(IW, W_deg, Modular);
            //cout << "Phi_deg = " << Phi_deg << ", W_deg = " << W_deg << endl;
            //cout << "TF = " << TF << ", InvPhi_Order = " << InvPhi_Order << ", IW = " << IW_Order << endl;
            //cout << "---------------" << endl;
            for (int k = 0; k <= (N/m)-1 ; k++){
                ZZ u, t;
                u = A_arr[k*m+j];
                t = A_arr[k*m+j+(m/2)];
                //cout << "u[" << k*m+j << "] = " << u << endl;
                //cout << "t[" << k*m+j+(m/2) << "] = " << t << endl;
                // upper
                A_arr[k*m+j] = AddMod(u, t, Modular);
                A_arr[k*m+j] = MulMod(A_arr[k*m+j], InvTwo, Modular);
                // down
                ZZ InvPhi_dot_IW = MulMod(InvPhi_Order, IW_Order, Modular);
                ZZ InvPhi_dot_IW_dot_Inv_two = MulMod(InvPhi_dot_IW, InvTwo, Modular);
                A_arr[k*m+j+(m/2)] = SubMod(u, t, Modular);
                A_arr[k*m+j+(m/2)] = MulMod(A_arr[k*m+j+(m/2)], InvPhi_dot_IW_dot_Inv_two, Modular);
                //cout << "A_arr[" << k*m+j << "] = " << A_arr[k*m+j] << endl;
                //cout << "A_arr[" << k*m+j+(m/2) << "] = " << A_arr[k*m+j+(m/2)] << endl;
            }
        }
    }
    for (int i = 0; i < N; i++){
        long long rev_idx = Bitrev.BitReserve(i, bit_num);
        arr[i] = A_arr[rev_idx];
    }
    return arr;
}

