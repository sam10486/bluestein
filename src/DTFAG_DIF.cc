#include <iostream>
#include "NTT.h"
#include "NTTSPMB.h"
#include "BitOperate.h"
#include "DTFAG.h"

using namespace std;

void DTFAG::DTFAG_DIF () {
    //------- radix sel-----
    int radix_r = 2;
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
    int arr_size = radix_r * radix_r;

    //int v1, v2;
    //int v0[radix_r] = {0};
    vector<ZZ > v0;
    vector<ZZ > v1, v2;

    v0.resize(radix_r);
    v1.resize(radix_r);
    v2.resize(radix_r);

    //int Tw0[radix_r] = {0};
    //int Tw1, Tw2;
    vector<ZZ > Tw0;
    vector<ZZ > Tw1, Tw2;

    Tw0.resize(radix_r);
    Tw1.resize(radix_r);
    Tw2.resize(radix_r);

    NTTSPMB NTTSPMB;
    std::ofstream DTFAG_DIF("./my_print_data/DTFAG_DIF.txt");
    std::ofstream DTFAG_TestPattern_Tw0("./SPMB_tw/DTFAG_TestPattern_Tw0.txt");
    std::ofstream DTFAG_TestPattern_Tw1("./SPMB_tw/DTFAG_TestPattern_Tw1.txt");
    std::ofstream DTFAG_TestPattern_Tw2("./SPMB_tw/DTFAG_TestPattern_Tw2.txt");

    int Tw2_display = 0;
    int Tw1_display = 0;
    int Tw0_display = false;
    int Tw_th = 1;
    int debug = 0;

    vector<vector<ZZ > > ROM0;
    vector<vector<ZZ > > ROM1, ROM2;
    ROM0.resize(radix_r);
    for(int i=0; i<radix_r; i++){
        ROM0[i].resize(radix_r);
    }
    ROM1.resize(arr_size);
    for(int i=0; i<radix_r; i++){
        ROM1[i].resize(radix_r);
    }
    ROM2.resize(arr_size);
    for(int i=0; i<radix_r; i++){
        ROM2[i].resize(radix_r);
    }


    DTFAG_DIF << "****************initial ROM********************" << endl;
    for(int n=0; n<radix_r; n++){
        DTFAG_DIF << "ROM0[" << n <<  "] = [";
        for(int k=0; k<radix_r; k++){
            ROM0[n][k] = (radix_r * radix_r) * k * n;
            ZZ ROM0_dg; 
            ROM0_dg = (radix_r * radix_r) * k * n;
            if(!debug) PowerMod(ROM0[n][k], fft_twiddle, ROM0_dg, fft_prime);
            DTFAG_DIF << ROM0[n][k] << ", ";
        }
        DTFAG_DIF << "]\n";
    }
    DTFAG_DIF << "----------------------------------" << endl;
    for(int n=0; n<radix_r; n++){
        DTFAG_DIF << "ROM1[" << n <<  "] = [";
        for(int k=0; k<radix_r; k++){
            ROM1[n][k] = (radix_r) * k * n;
            ZZ ROM1_dg; 
            ROM1_dg = (radix_r) * k * n;
            if(!debug) PowerMod(ROM1[n][k], fft_twiddle, ROM1_dg, fft_prime);
            DTFAG_DIF << ROM1[n][k] << ", ";
        }
        DTFAG_DIF << "]\n";
    }
    DTFAG_DIF << "----------------------------------" << endl;
    for(int n=0; n<radix_r; n++){
        DTFAG_DIF << "ROM2[" << n <<  "] = [";
        for(int k=0; k<radix_r; k++){
            ROM2[n][k] = k * n;
            ZZ ROM2_dg; 
            ROM2_dg = k * n;
            if(!debug) PowerMod(ROM2[n][k], fft_twiddle, ROM2_dg, fft_prime);
            DTFAG_DIF << ROM2[n][k] << ", ";
        }
        DTFAG_DIF << "]\n";
    }
    DTFAG_DIF << "**************intital fin****************" << endl;

    for(int i=0; i<radix_r; i++){
        DTFAG_DIF << "  i = " << i << endl;
        for(int t=0; t<radix_r; t++){
            DTFAG_DIF << "  t = " << t << endl;
            for(int j=0; j<radix_r; j++){
                DTFAG_DIF << "----------------------------" << endl;
                if(Tw2_display || Tw1_display || Tw0_display){
                    if(j == Tw_th) DTFAG_DIF << "       j = " << j << endl;
                }else{
                    DTFAG_DIF << "      j = " << j << " (TF" << j << ")"<< endl;
                }
                
                MA0 = j;
                DTFAG_DIF << "      MA0 = " << MA0 << endl;
                MA1 = NTTSPMB.Gray(i,radix_r);
                DTFAG_DIF << "      MA1 = " << MA1 << endl;
                if(i % 2 == 1){
                    int t_complement = number_complement.number_complement(t, radix_r);
                    MA2 = NTTSPMB.Gray(t_complement,radix_r);
                    DTFAG_DIF << "      MA2 : " << MA2 << " = " << "G(" << t_complement << ")"<< endl; 
                }else{
                    MA2 = NTTSPMB.Gray(t,radix_r);
                    DTFAG_DIF << "      MA2 : " << MA2 << " = " << "G(" << t << ")"<< endl; 
                }          

                for(int i=0; i<radix_r; i++){
                    v0[i] = ROM0[MA0][i];
                }
                for(int i=0; i<radix_r; i++){
                    v1[i] = ROM1[MA1][i];
                }
                for(int i=0; i<radix_r; i++){
                    v2[i] = ROM2[MA2][i];
                }

                //----------TW compute-----------
                for(int i=0; i<radix_r; i++){
                    if(Tw0_display){
                        if(j == Tw_th){
                            //-----------for debug-------------------
                            if(debug) DTFAG_DIF << "v0[" << i << "] = " << v0[i] << ", v1[" << i << "] = " << v1[i] << ", v2[" << i << "] = " << v2[i];
                            if(debug) Tw0[i] = v0[i] + v1[i] + v2[i];
                            if(debug) DTFAG_DIF << ", Tw0[" << i << "] = " << Tw0[i] << endl;
                            //----------------------------------------

                            //--------real compute------------
                            ZZ tmp;
                            if(!debug) DTFAG_DIF << "v0[" << i << "] = " << v0[i] << ", v1[" << i << "] = " << v1[i] << ", v2[" << i << "] = " << v2[i];
                            if(!debug) MulMod(tmp, v1[i], v2[i], fft_prime);
                            if(!debug) MulMod(Tw0[i], v0[i], tmp, fft_prime);
                            if(!debug) DTFAG_DIF << ", Tw0[" << i << "] = " << Tw0[i] << endl;
                            DTFAG_TestPattern_Tw0 << Tw0[i] << endl;
                            //---------------------------------
                        }
                    }else {
                        //-----------for debug-------------------
                        if(debug) DTFAG_DIF << "v0[" << i << "] = " << v0[i] << ", v1[" << i << "] = " << v1[i] << ", v2[" << i << "] = " << v2[i];
                        if(debug) Tw0[i] = v0[i] + v1[i] + v2[i];
                        if(debug) DTFAG_DIF << ", Tw0[" << i << "] = " << Tw0[i] << endl;
                        //----------------------------------------

                        //--------real compute------------
                        ZZ tmp;
                        if(!debug) DTFAG_DIF << "v0[" << i << "] = " << v0[i] << ", v1[" << i << "] = " << v1[i] << ", v2[" << i << "] = " << v2[i];
                        if(!debug) MulMod(tmp, v1[i], v2[i], fft_prime);
                        if(!debug) MulMod(Tw0[i], v0[i], tmp, fft_prime);
                        if(!debug) DTFAG_DIF << ", Tw0[" << i << "] = " << Tw0[i] << endl;
                        DTFAG_TestPattern_Tw0 << Tw0[i] << endl;
                        //---------------------------------
                    }
                    
                }
                  
                for(int i=0; i<radix_r; i++){
                   if(Tw1_display){
                       if(j == Tw_th){
                            //-----------for debug-------------------
                            if(debug) DTFAG_DIF << "v0[" << i << "] = " << v0[i] << ", v1[" << i << "] = " << v1[i];
                            if(debug) Tw1[i] = v0[i] + v1[i];
                            if(debug) DTFAG_DIF << ", Tw1[" << i << "] = " << Tw1[i] << endl;
                            //----------------------------------------

                            //--------real compute------------
                            if(!debug) DTFAG_DIF << "v0[" << i << "] = " << v0[i] << ", v1[" << i << "] = " << v1[i];
                            if(!debug) MulMod(Tw1[i], v1[i], v0[i], fft_prime);
                            if(!debug) DTFAG_DIF << ", Tw1[" << i << "] = " << Tw1[i] << endl;
                            DTFAG_TestPattern_Tw1 << Tw1[i] << endl;
                            //---------------------------------
                       }
                   }else {
                        //-----------for debug-------------------
                        if(debug) DTFAG_DIF << "v0[" << i << "] = " << v0[i] << ", v1[" << i << "] = " << v1[i];
                        if(debug) Tw1[i] = v0[i] + v1[i];
                        if(debug) DTFAG_DIF << ", Tw1[" << i << "] = " << Tw1[i] << endl;
                        //----------------------------------------

                        //--------real compute------------
                        if(!debug) DTFAG_DIF << "v0[" << i << "] = " << v0[i] << ", v1[" << i << "] = " << v1[i];
                        if(!debug) MulMod(Tw1[i], v1[i], v0[i], fft_prime);
                        if(!debug) DTFAG_DIF << ", Tw1[" << i << "] = " << Tw1[i] << endl;
                        DTFAG_TestPattern_Tw1 << Tw1[i] << endl;
                        //---------------------------------
                   }
                   
                }  

                for(int i=0; i<radix_r; i++){
                   if(Tw2_display){
                       if(j == Tw_th){
                            //--------real compute------------
                            DTFAG_DIF << "v0[" << i << "] = " << v0[i];
                            Tw2[i] = v0[i];
                            DTFAG_DIF << ", Tw2[" << i << "] = " << Tw2[i] << endl;
                            DTFAG_TestPattern_Tw2 << Tw2[i] << endl;
                            //---------------------------------
                       }
                   }else {
                        //--------real compute------------
                        DTFAG_DIF << "v0[" << i << "] = " << v0[i];
                        Tw2[i] = v0[i];
                        DTFAG_DIF << ", Tw2[" << i << "] = " << Tw2[i] << endl;
                        DTFAG_TestPattern_Tw2 << Tw2[i] << endl;
                        //---------------------------------
                   }
                   
                } 
                
                //--------------------------------
            }
        }
    }
}
