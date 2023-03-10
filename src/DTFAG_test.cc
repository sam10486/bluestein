#include <iostream>
#include "NTT.h"
#include "NTTSPMB.h"

using namespace std;

vector<int> DecToBin_test(int data, int bit_width);
int VecToInt_test(vector<int > bit_array, int N);
int number_complement_test(int i, int radix_r);


void DTFAG_test () {
    //------- radix sel-----
    int radix_r = 4;
    ZZ P ;
    ZZ W ;
    conv(P, "18446744069414584321");
    conv(W, "14603442835287214144");
    //-----------------------

    int MA0 = 0;
    int MA1 = 0;
    int MA2 = 0;
    int arr_size = radix_r * radix_r;

    //int v1, v2;
    //int v0[radix_r] = {0};
    ZZ v1, v2;
    vector<ZZ > v0;

    v0.resize(radix_r);

    //int Tw0[radix_r] = {0};
    //int Tw1, Tw2;
    ZZ Tw1, Tw2;
    vector<ZZ > Tw0;

    Tw0.resize(radix_r);

    NTTSPMB NTTSPMB;
    std::ofstream DTFAG("./my_print_data/DTFAG.txt");

    int Tw2_display = 0;
    int Tw1_display = 0;
    int Tw0_display = false;
    int Tw_th = 1;

    //int ROM0[radix_r][radix_r]; // ROM0[k][n]
    //int ROM1[arr_size]; // ROM1[group][n]
    //int ROM2[arr_size]; // ROM2[group][n]
    vector<vector<ZZ > > ROM0;
    vector<ZZ > ROM1, ROM2;
    ROM0.resize(radix_r);
    for(int i=0; i<radix_r; i++){
        ROM0[i].resize(radix_r);
    }
    ROM1.resize(arr_size);
    ROM2.resize(arr_size);


    DTFAG << "****************initial ROM********************" << endl;
    for(int n=0; n<radix_r; n++){
        DTFAG << "ROM0[" << n <<  "] = [";
        for(int k=0; k<radix_r; k++){
            ROM0[n][k] = (radix_r * radix_r) * k * n;
            ZZ ROM0_dg; 
            ROM0_dg = (radix_r * radix_r) * k * n;
            //PowerMod(ROM0[k][n], W, ROM0_dg, P);
            DTFAG << ROM0[n][k] << ", ";
        }
        DTFAG << "]\n";
    }
    DTFAG << "----------------------------------" << endl;
    for(int group=0; group<radix_r; group++){
        DTFAG << "ROM1[" << group <<  "] = [";
        for(int n=0; n<radix_r; n++){
            ROM1[group*radix_r+n] = group * radix_r * n;
            //ROM1[group*radix_r+n] = group * n;
            ZZ ROM1_dg;
            ROM1_dg = group * radix_r * n;
            //PowerMod(ROM1[group*radix_r+n], W, ROM1_dg, P);
            DTFAG <<  ROM1[group*radix_r+n] << ", ";
        }
        DTFAG << "]\n";
    }
    DTFAG << "----------------------------------" << endl;
    for(int group=0; group<radix_r; group++){
        DTFAG << "ROM2[" << group <<  "] = [";
        for(int n=0; n<radix_r; n++){
            ROM2[group*radix_r+n] = group * n;
            //ROM2[group*radix_r+n] = group * radix_r * n;
            ZZ ROM2_dg;
            ROM2_dg = group * n;
            //PowerMod(ROM2[group*radix_r+n], W, ROM2_dg, P);
            DTFAG << ROM2[group*radix_r+n] << ", ";
        }
        DTFAG << "]\n";
    }
    DTFAG << "**************intital fin****************" << endl;

    for(int t=0; t<radix_r; t++){
        DTFAG << "     t = " << t << endl;
        for(int i=0; i<radix_r; i++){
            DTFAG << "     i = " << i << endl;
            for(int j=0; j<radix_r; j++){
                DTFAG << "----------------------------" << endl;
                if(Tw2_display || Tw1_display || Tw0_display){
                    if(j == Tw_th) DTFAG << "     j = " << j << endl;
                }else{
                    DTFAG << "     j = " << j << " (TF" << j << ")"<< endl;
                }
                
                MA0 = j;
                DTFAG << "     MA0 = " << MA0 << endl;
                if(t % 2 == 0){
                    MA1 = radix_r * NTTSPMB.Gray(i,radix_r) + j;
                    //MA1 = radix_r * i + j;
                    DTFAG << "t = " << t << ", G(i) = " << NTTSPMB.Gray(i,radix_r) << ", MA1 : " << MA1 << " = " << radix_r << " * " << 
                                    NTTSPMB.Gray(i,radix_r) << " + " << j << endl;
                }else{
                    int i_complement = number_complement_test(i, radix_r);
                    MA1 = radix_r * NTTSPMB.Gray(i_complement,radix_r) + j; 
                    //MA1 = radix_r * i_complement + j;
                    DTFAG << "t = " << t << ", G(~i) = " << NTTSPMB.Gray(i_complement,radix_r) << ", MA1 : " << MA1 << " = " << radix_r << " * " << 
                                    NTTSPMB.Gray(i_complement,radix_r) << " + " << j << endl;
                }
                //**********************
                MA2 = radix_r * NTTSPMB.Gray(t,radix_r) + j;
                DTFAG << "MA2 : " << MA2 << " = " << radix_r << " * " << 
                                    NTTSPMB.Gray(t,radix_r) << " + " << j << endl; 
                //*******************               

                for(int i=0; i<radix_r; i++){
                    v0[i] = ROM0[MA0][i];
                }
                v1 = ROM1[MA1];
                v2 = ROM2[MA2];

                //----------TW compute-----------
                for(int i=0; i<radix_r; i++){
                    if(Tw0_display){
                        if(j == Tw_th){
                            //DTFAG << "(len): v0[" << i << "] = " << v0[i] << ", ";
                            Tw0[i] = v0[i];
                            //DTFAG << "Tw0: " << Tw0[i] << " = " << v0[i] << endl;
                        }
                    }else {
                        DTFAG << "(len): v0[" << i << "] = " << v0[i] << ", ";
                        Tw0[i] = v0[i];
                        DTFAG << "Tw0: " << Tw0[i] << " = " << v0[i] << endl;
                    }
                    
                }
                
                for(int index=0; index<radix_r; index++){
                    if(Tw1_display){
                        if(j == Tw_th){
                            /*DTFAG << "v0[" << index << "] = " << v0[index] << ", v1 = " << v1 << ", ";
                            Tw1 = v0[index] + v1;
                            DTFAG << "Tw1: " << Tw1 << " = " << v0[index] << " + " << v1 << endl;*/

                            //DTFAG << "(len): " << "v0[" << index << "] = " << v0[index] << ", v1 = " << v1 << ", ";
                            MulMod(Tw1, v0[index], v1, P);
                            //DTFAG << "Tw1: " << Tw1 << " = " << v0[index] << " * " << v1 << endl;

                        }
                    }else {
                        DTFAG << "v0[" << index << "] = " << v0[index] << ", v1 = " << v1 << ", ";
                        Tw1 = v0[index] + v1;
                        DTFAG << "Tw1: " << Tw1 << " = " << v0[index] << " + " << v1 << endl;

                        //DTFAG << "(len): " << "v0[" << index << "] = " << v0[index] << ", v1 = " << v1 << ", ";
                        MulMod(Tw1, v0[index], v1, P);
                        //DTFAG << " Tw1: " << Tw1 << " = " << v0[index] << " * " << v1 << endl;
                    }
                    
                }
                for(int index=0; index<radix_r; index++){
                    if(Tw2_display){
                        if(j == Tw_th) {
                            DTFAG << "v0[" << index << "] = " << v0[index] << ", v1 = " << v1 << ", v2 = " << v2 << ", ";
                            Tw2 = v0[index] + v1 + v2;
                            DTFAG << "Tw2: " << Tw2 << " = " << v0[index] << " + " << v1 << " + " << v2 << endl;

                            //DTFAG << "v0[" << index << "] = " << v0[index] << ", v1 = " << v1 << ", v2 = " << v2 << ", ";
                            ZZ tmp;
                            MulMod(tmp, v1, v2, P);
                            MulMod(Tw2, v0[index], tmp, P);
                            //DTFAG << "Tw2: " << Tw2 << " = " << v0[index] << " * " << v1 << " * " << v2 
                            //        << " ( index = " << index << " ) " << endl;
                        }
                    }else{
                        DTFAG << "v0[" << index << "] = " << v0[index] << ", v1 = " << v1 << ", v2 = " << v2 << ", ";
                        Tw2 = v0[index] + v1 + v2;
                        DTFAG << "Tw2: " << Tw2 << " = " << v0[index] << " + " << v1 << " + " << v2 << endl;

                        //DTFAG << "v0[" << index << "] = " << v0[index] << ", v1 = " << v1 << ", v2 = " << v2 << ", ";
                        //ZZ tmp;
                        //MulMod(tmp, v1, v2, P);
                        //MulMod(Tw2, v0[index], tmp, P);
                        //DTFAG << "Tw2: " << Tw2 << " = " << v0[index] << " * " << v1 << " * " << v2 
                        //        << " ( index = " << index << " ) " << endl;
                    }
                    
                    
                }

                //--------------------------------
            }
        }
    }
}

int number_complement_test(int i, int radix_r){
    
    vector<int > i_complement_array;
    int i_complement;
    int bit_width = (int)ceil(log2(radix_r));
    i_complement_array.resize(bit_width);

    //cout << "i = " << i << endl;

    vector<int > i_array = DecToBin_test(i, bit_width);
    int tmp;
    for(int i=0; i<bit_width; i++){
        tmp = i_array[i];
        if(tmp){
            i_complement_array[i] = 0;
        }else{
            i_complement_array[i] = 1;
        }
    }
    i_complement = VecToInt_test(i_complement_array, radix_r);

    //cout << "i_complement = " << i_complement << endl;
    return i_complement;
}

vector<int> DecToBin_test(int data, int bit_width){
    vector<int> BinVec(bit_width);
    for(int j=0; j<bit_width; j++){
        BinVec.at(j) = (data >> j) & 1;
    }
    return BinVec;
}
int VecToInt_test(vector<int > bit_array, int N){
    int bit_width = (int)ceil(log2(N));
    int integer = 0;
    for(int j=0; j < bit_width; j++){
        integer += bit_array[j] << j;
    }
    return integer;
}