#include <iostream>
#include "NTT.h"
#include "NTTSPMB.h"

using namespace std;

vector<int> DecToBin(int data, int bit_width);
int VecToInt(vector<int > bit_array, int N);
int number_complement(int i, int radix_r1);


void DTFAG (ZZ P, ZZ W) {
    //------- radix sel-----
    int radix_r1 = 4;
    int radix_r2 = 4;
    //-----------------------

    int MA0 = 0;
    int MA1 = 0;
    int MA2 = 0;
    int arr_size = radix_r1 * radix_r1;

    ZZ v1, v2;
    vector<ZZ > v0;

    v0.resize(radix_r1);

    //int Tw0[radix_r1] = {0};
    //int Tw1, Tw2;
    ZZ Tw1, Tw2;
    vector<ZZ > Tw0;

    Tw0.resize(radix_r1);

    NTTSPMB NTTSPMB;
    std::ofstream DTFAG("./my_print_data/DTFAG.txt");
    std::ofstream DTFAG_pattern_Tw0("./my_print_data/DTFAG_pattern_Tw0.txt");
    std::ofstream DTFAG_pattern_Tw1("./my_print_data/DTFAG_pattern_Tw1.txt");
    std::ofstream DTFAG_pattern_Tw2("./my_print_data/DTFAG_pattern_Tw2.txt");

    int Tw2_display = 0;
    int Tw1_display = 0;
    int Tw0_display = 0;

    int Addr_display = 1;
    int Tw_th = 1;

    vector<vector<ZZ > > ROM0;
    vector<ZZ > ROM1, ROM2;
    ROM0.resize(radix_r1);
    for(int i=0; i<radix_r1; i++){
        ROM0[i].resize(radix_r1);
    }
    ROM1.resize(arr_size);
    ROM2.resize(arr_size);

    //-------------test pattern array-------------
    vector<vector<vector<vector<ZZ > > > > Tw0_ROM;
    vector<vector<vector<vector<ZZ > > > > Tw1_ROM;
    vector<vector<vector<vector<ZZ > > > > Tw2_ROM;

    Tw0_ROM.resize(radix_r1);
    for(int t=0; t<radix_r1; t++){
        Tw0_ROM[t].resize(radix_r1);
        for(int i=0; i<radix_r1; i++){
            Tw0_ROM[t][i].resize(radix_r1);
            for(int len_idx=0; len_idx<radix_r1; len_idx++){
                Tw0_ROM[t][i][len_idx].resize(radix_r1);
            }
        }
    }

    Tw1_ROM.resize(radix_r1);
    for(int t=0; t<radix_r1; t++){
        Tw1_ROM[t].resize(radix_r1);
        for(int i=0; i<radix_r1; i++){
            Tw1_ROM[t][i].resize(radix_r1);
            for(int len_idx=0; len_idx<radix_r1; len_idx++){
                Tw1_ROM[t][i][len_idx].resize(radix_r1);
            }
        }
    }

    Tw2_ROM.resize(radix_r1);
    for(int t=0; t<radix_r1; t++){
        Tw2_ROM[t].resize(radix_r1);
        for(int i=0; i<radix_r1; i++){
            Tw2_ROM[t][i].resize(radix_r1);
            for(int len_idx=0; len_idx<radix_r1; len_idx++){
                Tw2_ROM[t][i][len_idx].resize(radix_r1);
            }
        }
    }
    //--------------------------------------------


    DTFAG << "****************initial ROM********************" << endl;
    for(int n=0; n<radix_r1; n++){
        DTFAG << "n = " << n << ", ROM0[" << n <<  "] = [";
        for(int k=0; k<radix_r1; k++){
            ZZ ROM0_dg; 
            ROM0_dg = (radix_r1 * radix_r1) * n * k;
            PowerMod(ROM0[n][k], W, ROM0_dg, P);
            //------------------addr-------------------
            if(Addr_display) ROM0[n][k] = (radix_r1 * radix_r1) * k * n;
            //------------------addr-------------------
            DTFAG << ROM0[n][k] << ", ";
        }
        DTFAG << "]\n";
    }
    DTFAG << "----------------------------------" << endl;
    for(int group=0; group<radix_r1; group++){
        DTFAG << "ROM1[" << group <<  "] = [";
        for(int n=0; n<radix_r1; n++){
            ZZ ROM1_dg;
            ROM1_dg = group * radix_r1 * n;
            PowerMod(ROM1[group*radix_r1+n], W, ROM1_dg, P);
            //------------------addr-------------------
            if(Addr_display) ROM1[group*radix_r1+n] = group * radix_r1 * n;
            //------------------addr-------------------
            DTFAG <<  ROM1[group*radix_r1+n] << ", ";
        }
        DTFAG << "]\n";
    }
    DTFAG << "----------------------------------" << endl;
    for(int group=0; group<radix_r1; group++){
        DTFAG << "ROM2[" << group <<  "] = [";
        for(int n=0; n<radix_r1; n++){
            ZZ ROM2_dg;
            ROM2_dg = group * n;
            PowerMod(ROM2[group*radix_r1+n], W, ROM2_dg, P);
            //------------------addr-------------------
            if(Addr_display) ROM2[group*radix_r1+n] = group * n;
            //------------------addr-------------------
            DTFAG << ROM2[group*radix_r1+n] << ", ";
        }
        DTFAG << "]\n";
    }
    DTFAG << "**************intital fin****************" << endl;

    for(int t=0; t<radix_r1; t++){
        DTFAG << "     t = " << t << endl;
        for(int i=0; i<radix_r1; i++){
            DTFAG << "     i = " << i << endl;
            // length_idx indicate the length in NTTSPMB
            for(int len_idx=0; len_idx<radix_r1; len_idx++){
                // j indicate the sequence of TFs (ex: j=1 => TF1, j=2 => TF2 ... )
                DTFAG << "-----------------crossbar for len--------------------"<< endl;
                DTFAG << "  len_idx = " << len_idx << endl;
                for(int j=0; j<radix_r1; j++){
                    if(Tw2_display || Tw1_display || Tw0_display){
                        if(j == Tw_th) DTFAG << "     j = " << j << endl;
                    }else{
                        DTFAG << "     j = " << j << endl;
                    }
                    MA0 = j;
                    //*************test*****
                    MA1 = radix_r1 * NTTSPMB.Gray(t,radix_r1) + j;
                    if(Tw2_display || Tw1_display || Tw0_display){
                        if( j == Tw_th) {
                            DTFAG << "t = " << t << ", G(t) = " << NTTSPMB.Gray(t,radix_r1) 
                            << ", i = " << i << ", MA1 : " << MA1 << " = " << radix_r1 << " * "  
                            << NTTSPMB.Gray(t,radix_r1) << " + " << j << endl;
                        } 
                    }else{
                        //DTFAG << "t = " << t << ", G(t) = " << NTTSPMB.Gray(t,radix_r1) 
                        //<< ", i = " << i << ", MA1 : " << MA1 << " = " << radix_r1 << " * "  
                        //<< NTTSPMB.Gray(t,radix_r1) << " + " << j << endl;
                    }
                    //**********************

                    //******test********
                    if(t % 2 == 0) {
                        MA2 = radix_r1 * NTTSPMB.Gray(i,radix_r1) + j;
                        if(Tw2_display || Tw1_display || Tw0_display){
                            if(j == Tw_th){
                                //DTFAG << "MA2 : " << MA2 << " = " << radix_r1 << " * " << 
                                //            NTTSPMB.Gray(i,radix_r1) << " + " << j << endl; 
                            }
                        }else {
                            //DTFAG << "MA2 : " << MA2 << " = " << radix_r1 << " * " << 
                            //                NTTSPMB.Gray(i,radix_r1) << " + " << j << endl; 
                        }
                    } else {
                        int i_complement = number_complement(i, radix_r1);
                        MA2 = radix_r1 * NTTSPMB.Gray(i_complement,radix_r1) + j;
                        if(Tw2_display || Tw1_display || Tw0_display){
                            if(j == Tw_th){
                                //DTFAG << "MA2 : " << MA2 << " = " << radix_r1 << " * " << 
                                //            NTTSPMB.Gray(i_complement,radix_r1) << " + " << j << endl; 
                            }
                        }else {
                            //DTFAG << "MA2 : " << MA2 << " = " << radix_r1 << " * " << 
                            //                NTTSPMB.Gray(i_complement,radix_r1) << " + " << j << endl; 
                        }
                    }
                    //*******************               
                    for(int i=0; i<radix_r1; i++){
                        //int Tw0_len_idx = i % radix_r2 ;
                        //v0[i] = ROM0[MA0][Tw0_len_idx];
                        v0[i] = ROM0[MA0][i];
                    }
                    
                    v1 = ROM1[MA1];
                    v2 = ROM2[MA2];

                    //int Tw0_len_idx = len_idx % radix_r2 ;
                    //int Tw1_len_idx = len_idx % radix_r1 ;
                    //int Tw2_len_idx = len_idx % radix_r1 ;

                    if(Tw0_display){
                        if(j == Tw_th){
                            //DTFAG << "(TF" << j << "): " << "len_idx = " << len_idx << ", ";
                            Tw0[len_idx] = v0[len_idx];
                            //DTFAG << "Tw0: " << Tw0[len_idx] << " = " << v0[len_idx] << endl;
                        }
                    }else {
                        //DTFAG << "(TF" << j << "): " << "len_idx = " << len_idx << ", ";
                        Tw0[len_idx] = v0[len_idx];
                        //DTFAG << "Tw0: " << Tw0[len_idx] << " = " << v0[len_idx] << endl;
                    }

                    if(Tw1_display){
                        if(j == Tw_th){   
                            if(Addr_display){
                                //----------------addr---------------------
                                DTFAG << "(TF" << j << "): " << "len_idx = " << len_idx << ", ";
                                Tw1 = v0[len_idx] + v1;
                                DTFAG << "Tw1: " << Tw1 << " = " << v0[len_idx] << " + " << v1 << endl;     
                                //-----------------------------------------
                            }else{
                                //DTFAG << "(TF" << j << "): " << "len_idx = " << len_idx << ", ";
                                MulMod(Tw1, v0[len_idx], v1, P);
                                //DTFAG << "Tw1: " << Tw1 << " = " << v0[len_idx] << " * " << v1 << endl;     
                            }
                        }
                    }else {
                        if(Addr_display){
                            //----------------addr---------------------
                            //DTFAG << "(TF" << j << "): " << "len_idx = " << len_idx << ", ";
                            Tw1 = v0[len_idx] + v1;
                            //DTFAG << "Tw1: " << Tw1 << " = " << v0[len_idx] << " + " << v1 << endl;
                            //-----------------------------------------
                        }else{
                            //DTFAG << "(TF" << j << "): " << "len_idx = " << len_idx << ", ";
                            MulMod(Tw1, v0[len_idx], v1, P);
                            //DTFAG << "Tw1: " << Tw1 << " = " << v0[len_idx] << " * " << v1 << endl;     
                        }
                    }

                    if(Tw2_display){
                        if(j == Tw_th) {
                            if(Addr_display){
                                //----------------addr---------------------
                                //DTFAG << "(TF" << j << "): " << "len_idx = " << len_idx << ", ";
                                Tw2 = v0[len_idx] + v1 + v2;
                                //DTFAG << "Tw2: " << Tw2 << " = " << v0[len_idx] << " + " << v1 << " + " << v2 << endl;
                                //-----------------------------------------
                            }else{
                                //DTFAG << "(TF" << j << "): " << "len_idx = " << len_idx << ", ";
                                ZZ tmp;
                                MulMod(tmp, v1, v2, P);
                                MulMod(Tw2, v0[len_idx], tmp, P);
                                //DTFAG << "Tw2: " << Tw2 << " = " << v0[len_idx] << " * " << v1 << " * " << v2 << endl;
                            }
                        }
                    }else{
                        if(Addr_display){
                            //----------------addr---------------------
                            DTFAG << "(TF" << j << "): " << "len_idx = " << len_idx << ", ";
                            Tw2 = v0[len_idx] + v1 + v2;
                            DTFAG << "Tw2: " << Tw2 << " = " << v0[len_idx] << " + " << v1 << " + " << v2 << endl;
                            //-----------------------------------------
                        }else{
                            //DTFAG << "(TF" << j << "): "  << "len_idx = " << len_idx << ", ";
                            ZZ tmp;
                            MulMod(tmp, v1, v2, P);
                            MulMod(Tw2, v0[len_idx], tmp, P);
                            //DTFAG << "Tw2: " << Tw2 << " = " << v0[len_idx] << " * " << v1 << " * " << v2 << endl;
                        }
                        
                    }

                    Tw0_ROM[t][i][len_idx][j] = Tw0[len_idx];
                    Tw1_ROM[t][i][len_idx][j] = Tw1;
                    Tw2_ROM[t][i][len_idx][j] = Tw2;

                    //cout << "st0_ROM[" << t << "][" << i << "][" << len_idx << "][" << j << "] = " <<  st0_ROM[t][i][len_idx][j] << endl;
                }          
            }
        }
    }

    for(int t=0; t<radix_r1; t++){
        for(int i=0; i<radix_r1; i++){
            for(int len_idx=0; len_idx<radix_r1; len_idx++){
                for(int j=0; j<radix_r1; j++){
                    //cout << "st0_ROM[" << t << "][" << i << "][" << len_idx << "][" << j << "] = " <<  st0_ROM[t][i][len_idx][j] << endl;
                    DTFAG_pattern_Tw0 << Tw0_ROM[t][i][len_idx][j] << "\n";
                    DTFAG_pattern_Tw1 << Tw1_ROM[t][i][len_idx][j] << "\n";
                    DTFAG_pattern_Tw2 << Tw2_ROM[t][i][len_idx][j] << "\n";
                }
            }
        }
    }

}

int number_complement(int i, int radix_r1){
    
    vector<int > i_complement_array;
    int i_complement;
    int bit_width = (int)ceil(log2(radix_r1));
    i_complement_array.resize(bit_width);

    //cout << "i = " << i << endl;

    vector<int > i_array = DecToBin(i, bit_width);
    int tmp;
    for(int i=0; i<bit_width; i++){
        tmp = i_array[i];
        if(tmp){
            i_complement_array[i] = 0;
        }else{
            i_complement_array[i] = 1;
        }
    }
    i_complement = VecToInt(i_complement_array, radix_r1);

    //cout << "i_complement = " << i_complement << endl;
    return i_complement;
}

vector<int> DecToBin(int data, int bit_width){
    vector<int> BinVec(bit_width);
    for(int j=0; j<bit_width; j++){
        BinVec.at(j) = (data >> j) & 1;
    }
    return BinVec;
}
int VecToInt(vector<int > bit_array, int N){
    int bit_width = (int)ceil(log2(N));
    int integer = 0;
    for(int j=0; j < bit_width; j++){
        integer += bit_array[j] << j;
    }
    return integer;
}