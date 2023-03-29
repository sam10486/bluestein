#include <iostream>
#include "NTT.h"
#include "NTTSPMB.h"
#include "BitOperate.h"
#include "DTFAG.h"

using namespace std;

void DTFAG::DTFAG_DIT() {
    //------- radix sel-----
    int radix_r1 = 2;
    int radix_r2 = 2;
    unsigned long fft_point = 16;
    ZZ fft_prime ;
    ZZ fft_twiddle_65536 ;
    ZZ fft_twiddle ;
    long difference_length = 65536 / fft_point ;
    conv(fft_prime, "18446744069414584321"); // prime number
    conv(fft_twiddle_65536, "14603442835287214144"); // twiddle factor based setting by main.cc
    
    PowerMod(fft_twiddle,fft_twiddle_65536,difference_length,fft_prime);
    cout << "fft_twiddle = " << fft_twiddle << ", fft_prime = " << fft_prime << endl;
    //-----------------------
    BitOperate number_complement;

    int MA0 = 0;
    int MA1 = 0;
    int MA2 = 0;
    int arr_size = radix_r1 * radix_r1;

    ZZ v1, v2;
    vector<ZZ > v0;
    v0.resize(radix_r1);


    vector<ZZ > Tw0;
    Tw0.resize(radix_r1);
    ZZ Tw1, Tw2;
    
    
    
    NTTSPMB NTTSPMB;
    std::ofstream DTFAG_DIT("./my_print_data/DTFAG_DIT.txt");
    //std::ofstream DTFAG_pattern_Tw0("./my_print_data/DTFAG_pattern_Tw0.txt");
    //std::ofstream DTFAG_pattern_Tw1("./my_print_data/DTFAG_pattern_Tw1.txt");
    //std::ofstream DTFAG_pattern_Tw2("./my_print_data/DTFAG_pattern_Tw2.txt");

    int Tw2_display = 0;
    int Tw1_display = 0;
    int Tw0_display = 0;

    int debug = 1;
    int Tw_th = 1;

    vector<vector<ZZ > > ROM0;
    vector<ZZ > ROM1, ROM2;
    ROM0.resize(radix_r1);
    for(int i=0; i<radix_r1; i++){
        ROM0[i].resize(radix_r1);
    }
    ROM1.resize(arr_size);
    ROM2.resize(arr_size);

    DTFAG_DIT << "****************initial ROM********************" << endl;
    for(int k=0; k<radix_r1; k++){
        DTFAG_DIT << "ROM0[k=" << k <<  "] = [";
        for(int n=0; n<radix_r1; n++){
            ZZ ROM0_dg; 
            ROM0_dg = (radix_r1 * radix_r1) * n * k;
            PowerMod(ROM0[k][n], fft_twiddle, ROM0_dg, fft_prime);
            //------------------for debug-------------------
            if(debug) ROM0[k][n] = (radix_r1 * radix_r1) * k * n;
            //----------------------------------------------
            DTFAG_DIT << ROM0[k][n] << ", ";
        }
        DTFAG_DIT << "]\n";
    }
    DTFAG_DIT << "----------------------------------" << endl;
    for(int group=0; group<radix_r1; group++){
        DTFAG_DIT << "ROM1[g=" << group <<  "] = [";
        for(int n=0; n<radix_r1; n++){
            ZZ ROM1_dg;
            ROM1_dg = group * radix_r1 * n;
            PowerMod(ROM1[group*radix_r1+n], fft_twiddle, ROM1_dg, fft_prime);
            //------------------for debug-------------------
            if(debug) ROM1[group*radix_r1+n] = group * radix_r1 * n;
            //----------------------------------------------
            DTFAG_DIT <<  ROM1[group*radix_r1+n] << ", ";
        }
        DTFAG_DIT << "]\n";
    }
    DTFAG_DIT << "----------------------------------" << endl;
    for(int group=0; group<radix_r1; group++){
        DTFAG_DIT << "ROM2[g=" << group <<  "] = [";
        for(int n=0; n<radix_r1; n++){
            ZZ ROM2_dg;
            ROM2_dg = group * n;
            PowerMod(ROM2[group*radix_r1+n], fft_twiddle, ROM2_dg, fft_prime);
            //------------------for debug-------------------
            if(debug) ROM2[group*radix_r1+n] = group * n;
            //----------------------------------------------
            DTFAG_DIT << ROM2[group*radix_r1+n] << ", ";
        }
        DTFAG_DIT << "]\n";
    }
    DTFAG_DIT << "**************intital fin****************" << endl;

    for(int t=0; t<radix_r1; t++){
        DTFAG_DIT << "     t = " << t << endl;
        for(int i=0; i<radix_r1; i++){
            DTFAG_DIT << "-----------------crossbar for len--------------------"<< endl;
            DTFAG_DIT << "     i = " << i << endl;
            for(int j=0; j<radix_r1; j++){
                if(Tw2_display || Tw1_display || Tw0_display){
                    if(j == Tw_th) DTFAG_DIT << "     j = " << j << endl;
                }else{
                    DTFAG_DIT << "     j = " << j << endl;
                }
                MA0 = j;
                DTFAG_DIT << "MA0 = " << MA0 << endl;
                if(t % 2== 0) {
                    MA1 = radix_r1 * NTTSPMB.Gray(i,radix_r1) + j;
                    DTFAG_DIT << "MA1 = " << MA1 << endl;
                }else {
                    int i_complement = number_complement.number_complement(i, radix_r1);
                    MA1 = radix_r1 * NTTSPMB.Gray(i_complement,radix_r1) + j;
                    DTFAG_DIT << "MA1 = " << MA1 << endl;
                }
                MA2 = radix_r1 * NTTSPMB.Gray(t, radix_r1) + j;
                DTFAG_DIT << "MA2 = " << MA2 << endl;
                
                for(int i=0; i<radix_r1; i++){
                    v0[i] = ROM0[MA0][i];
                }
                v1 = ROM1[MA1];
                v2 = ROM2[MA2];                

                
                for(int idx=0; idx<radix_r1; idx++){
                    DTFAG_DIT << "v0[" << idx << "] = " << v0[idx];
                    Tw0[idx] = v0[idx];
                    DTFAG_DIT << ", Tw0[" << idx << "] = " << Tw0[idx] << endl;
                }
                

                for(int idx=0; idx<radix_r1; idx++){
                    //------------for debug--------------------
                    if(debug) DTFAG_DIT << "v0[" << idx << "] = " << v0[idx] << ", v1 = " << v1;
                    if(debug) Tw1 = v0[idx] + v1;
                    if(debug) DTFAG_DIT << ", Tw1[" << idx << "] = " << Tw1 << endl;
                    //------------------------------------------

                    //--------------for compute--------------------
                    if(!debug) DTFAG_DIT << "v0[" << idx << "] = " << v0[idx] << ", v1 = " << v1;
                    if(!debug) MulMod(Tw1, v0[idx], v1, fft_prime);
                    if(!debug) DTFAG_DIT << ", Tw1[" << idx << "] = " << Tw1 << endl;
                    //---------------------------------------------
                }

                for(int idx=0; idx<radix_r1; idx++){
                    //------------for debug--------------------
                    if(debug) DTFAG_DIT << "v0[" << idx << "] = " << v0[idx] << ", v1 = " << v1 << ", v2 = " << v2;
                    if(debug) Tw2 = v0[idx] + v1 + v2;
                    if(debug) DTFAG_DIT << ", Tw2[" << idx << "] = " << Tw2 << endl;
                    //------------------------------------------

                    //--------------for compute--------------------
                    ZZ tmp;
                    if(!debug) DTFAG_DIT << "v0[" << i << "] = " << v0[idx] << ", v1 = " << v1 << ", v2 = " << v2;
                    if(!debug) MulMod(tmp, v0[idx], v1, fft_prime);
                    if(!debug) MulMod(Tw2, tmp, v2, fft_prime);
                    if(!debug) DTFAG_DIT << ", Tw2[" << idx << "] = " << Tw2 << endl;
                    //---------------------------------------------
                }
            }          
        }
    }
}
