#include <iostream>
#include "NTT.h"
#include "NTTSPMB.h"
#include "SPMB.h"
#include "BitOperate.h"
#include "DTFAG.h"

using namespace std;

void DTFAG::DTFAG_DIF_MixedRadix () {
    //------- radix sel-----
    int radix_r1 = 2;
    int radix_r2 = 2;
    unsigned long fft_point = pow(radix_r1, 3) * radix_r2;
    ZZ fft_prime ;
    ZZ fft_twiddle_65536 ;
    ZZ fft_twiddle ;
    long difference_length = 65536 / fft_point ;
    conv(fft_prime, "18446744069414584321"); // prime number
    conv(fft_twiddle_65536, "14603442835287214144"); // twiddle factor based setting by main.cc
    
    PowerMod(fft_twiddle,fft_twiddle_65536,difference_length,fft_prime);
    cout << "fft_twiddle = " << fft_twiddle << ", fft_prime = " << fft_prime << endl;
    //-----------------------
    
    SPMB spmb;
    BitOperate number_complement;

    int MA0 = 0;
    int MA1 = 0;
    int MA2 = 0;
    int MA_MixedRadix = 0;
    int arr_size = radix_r1 * radix_r1;

    //int v1, v2;
    //int v0[radix_r1] = {0};
    vector<ZZ > v0;
    vector<ZZ > v1, v2;
    vector<ZZ > v_MixedRadix;

    v0.resize(radix_r1);
    v1.resize(radix_r1);
    v2.resize(radix_r1);
    v_MixedRadix.resize(radix_r1);

    //int Tw0[radix_r1] = {0};
    //int Tw1, Tw2;
    vector<ZZ > Tw0;
    vector<ZZ > Tw1, Tw2;

    Tw0.resize(radix_r1);
    Tw1.resize(radix_r1);
    Tw2.resize(radix_r1);

    NTTSPMB NTTSPMB;
    std::ofstream DTFAG_DIF_MixedRadix("./my_print_data/DTFAG_DIF_MixedRadix.txt");
    std::ofstream DTFAG_TestPattern_Tw0("./SPMB_tw/MixR_DTFAG_DIF_TestPattern_Tw0.txt");
    std::ofstream DTFAG_TestPattern_Tw1("./SPMB_tw/MixR_DTFAG_DIF_TestPattern_Tw1.txt");
    std::ofstream DTFAG_TestPattern_Tw2("./SPMB_tw/MixR_DTFAG_DIF_TestPattern_Tw2.txt");
    //---------------ZZ to hex-----------------
    std::ofstream hex_DTFAG_MixedRadix("./SPMB_tw/hex_DTFAG_DIF_MixedRadix.txt");
    //-----------------------------------------

    int Tw2_display = 0;
    int Tw1_display = 0;
    int Tw0_display = false;
    int Tw_th = 1;
    int debug = 1;

    vector<vector<ZZ > > ROM0;
    vector<vector<ZZ > > ROM1, ROM2;
    ROM0.resize(radix_r1);
    for(int i=0; i<radix_r1; i++){
        ROM0[i].resize(radix_r1);
    }
    ROM1.resize(arr_size);
    for(int i=0; i<radix_r1; i++){
        ROM1[i].resize(radix_r1);
    }
    ROM2.resize(arr_size);
    for(int i=0; i<radix_r1; i++){
        ROM2[i].resize(radix_r1);
    }


    DTFAG_DIF_MixedRadix << "****************initial ROM********************" << endl;
    for(int k=0; k<radix_r1; k++){
        DTFAG_DIF_MixedRadix << "ROM0[" << "k=" << k <<  "] = [";
        for(int n=0; n<radix_r1; n++){
            ROM0[k][n] = (radix_r1 * radix_r2) * k * n;
            ZZ ROM0_dg; 
            ROM0_dg = (radix_r1 * radix_r2) * k * n;
            if(!debug) PowerMod(ROM0[k][n], fft_twiddle, ROM0_dg, fft_prime);
            DTFAG_DIF_MixedRadix << ROM0[k][n] << ", ";
            hex_DTFAG_MixedRadix << "ROM0_array[" << k << "][" << n << "] = "<< "64'h" << spmb.ZZtohex(ROM0[k][n]) << std::endl;
        }
        DTFAG_DIF_MixedRadix << "]\n";
    }
    DTFAG_DIF_MixedRadix << "----------------------------------" << endl;
    for(int k=0; k<radix_r2; k++){
        DTFAG_DIF_MixedRadix << "ROM1[" << "k=" << k <<  "] = [";
        for(int n=0; n<radix_r1; n++){
            ROM1[k][n] = (radix_r1) * k * n;
            ZZ ROM1_dg; 
            ROM1_dg = (radix_r1) * k * n;
            if(!debug) PowerMod(ROM1[k][n], fft_twiddle, ROM1_dg, fft_prime);
            DTFAG_DIF_MixedRadix << ROM1[k][n] << ", ";
        }
        DTFAG_DIF_MixedRadix << "]\n";
    }
    DTFAG_DIF_MixedRadix << "----------------------------------" << endl;
    for(int k=0; k<radix_r1; k++){
        DTFAG_DIF_MixedRadix << "ROM2[" << "k=" << k <<  "] = [";
        for(int n=0; n<radix_r1; n++){
            ROM2[k][n] = k * n;
            ZZ ROM2_dg; 
            ROM2_dg = k * n;
            if(!debug) PowerMod(ROM2[k][n], fft_twiddle, ROM2_dg, fft_prime);
            DTFAG_DIF_MixedRadix << ROM2[k][n] << ", ";
        }
        DTFAG_DIF_MixedRadix << "]\n";
    }
    DTFAG_DIF_MixedRadix << "**************intital fin****************" << endl;

    for(int i=0; i<radix_r2; i++){
        DTFAG_DIF_MixedRadix << "  i = " << i << endl;
        for(int t=0; t<radix_r1; t++){
            DTFAG_DIF_MixedRadix << "  t = " << t << endl;
            for(int j=0; j<radix_r1; j++){
                DTFAG_DIF_MixedRadix << "----------------------------" << endl;
                if(Tw2_display || Tw1_display || Tw0_display){
                    if(j == Tw_th) DTFAG_DIF_MixedRadix << "       j = " << j << endl;
                }else{
                    DTFAG_DIF_MixedRadix << "      j = " << j << " (TF" << j << ")"<< endl;
                }
                
                MA0 = j;
                DTFAG_DIF_MixedRadix << "      MA0 = " << MA0 << endl;
                MA1 = NTTSPMB.Gray(i,radix_r1);
                DTFAG_DIF_MixedRadix << "      MA1 = " << MA1 << endl;
                if(i % 2 == 1){
                    int t_complement = number_complement.number_complement(t, radix_r1);
                    MA2 = NTTSPMB.Gray(t_complement,radix_r1);
                    DTFAG_DIF_MixedRadix << "      MA2 : " << MA2 << " = " << "G(" << t_complement << ")"<< endl; 
                }else{
                    MA2 = NTTSPMB.Gray(t,radix_r1);
                    DTFAG_DIF_MixedRadix << "      MA2 : " << MA2 << " = " << "G(" << t << ")"<< endl; 
                }        
                if(radix_r1 == radix_r2) {
                    MA_MixedRadix = j;
                    DTFAG_DIF_MixedRadix << "      MA_MixedRadix = " << MA_MixedRadix << endl;
                }else {
                    MA_MixedRadix =  (j * (radix_r1/radix_r2)) % radix_r1;
                    DTFAG_DIF_MixedRadix << "      MA_MixedRadix = " << MA_MixedRadix << endl;
                }
                
                for(int i=0; i<radix_r1; i++){
                    v0[i] = ROM0[MA0][i];
                }
                for(int i=0; i<radix_r1; i++){
                    v1[i] = ROM1[MA1][i];
                }
                for(int i=0; i<radix_r1; i++){
                    v2[i] = ROM2[MA2][i];
                }
                for(int i=0; i<radix_r1; i++){
                    v_MixedRadix[i] = ROM0[MA_MixedRadix][i];
                }
                
                //----------TW compute-----------
                for(int i=0; i<radix_r1; i++){
                    if(Tw0_display){
                        if(j == Tw_th){
                            //-----------for debug-------------------
                            if(debug) DTFAG_DIF_MixedRadix << "v0[" << i << "] = " << v0[i] << ", v1[" << i << "] = " << v1[i] << ", v2[" << i << "] = " << v2[i];
                            if(debug) Tw0[i] = v0[i] + v1[i] + v2[i];
                            if(debug) DTFAG_DIF_MixedRadix << ", Tw0[" << i << "] = " << Tw0[i] << endl;
                            //----------------------------------------

                            //--------real compute------------
                            ZZ tmp;
                            if(!debug) DTFAG_DIF_MixedRadix << "v0[" << i << "] = " << v0[i] << ", v1[" << i << "] = " << v1[i] << ", v2[" << i << "] = " << v2[i];
                            if(!debug) MulMod(tmp, v1[i], v2[i], fft_prime);
                            if(!debug) MulMod(Tw0[i], v0[i], tmp, fft_prime);
                            if(!debug) DTFAG_DIF_MixedRadix << ", Tw0[" << i << "] = " << Tw0[i] << endl;
                            DTFAG_TestPattern_Tw0 << Tw0[i] << endl;
                            //---------------------------------
                        }
                    }else {
                        //-----------for debug-------------------
                        if(debug) DTFAG_DIF_MixedRadix << "v0[" << i << "] = " << v0[i] << ", v1[" << i << "] = " << v1[i] << ", v2[" << i << "] = " << v2[i];
                        if(debug) Tw0[i] = v0[i] + v1[i] + v2[i];
                        if(debug) DTFAG_DIF_MixedRadix << ", Tw0[" << i << "] = " << Tw0[i] << endl;
                        //----------------------------------------

                        //--------real compute------------
                        ZZ tmp;
                        if(!debug) DTFAG_DIF_MixedRadix << "v0[" << i << "] = " << v0[i] << ", v1[" << i << "] = " << v1[i] << ", v2[" << i << "] = " << v2[i];
                        if(!debug) MulMod(tmp, v1[i], v2[i], fft_prime);
                        if(!debug) MulMod(Tw0[i], v0[i], tmp, fft_prime);
                        if(!debug) DTFAG_DIF_MixedRadix << ", Tw0[" << i << "] = " << Tw0[i] << endl;
                        DTFAG_TestPattern_Tw0 << Tw0[i] << endl;
                        //---------------------------------
                    }
                    
                }
                  
                for(int i=0; i<radix_r1; i++){
                   if(Tw1_display){
                       if(j == Tw_th){
                            //-----------for debug-------------------
                            if(debug) DTFAG_DIF_MixedRadix << "v0[" << i << "] = " << v0[i] << ", v1[" << i << "] = " << v1[i];
                            if(debug) Tw1[i] = v0[i] + v1[i];
                            if(debug) DTFAG_DIF_MixedRadix << ", Tw1[" << i << "] = " << Tw1[i] << endl;
                            //----------------------------------------

                            //--------real compute------------
                            if(!debug) DTFAG_DIF_MixedRadix << "v0[" << i << "] = " << v0[i] << ", v1[" << i << "] = " << v1[i];
                            if(!debug) MulMod(Tw1[i], v1[i], v0[i], fft_prime);
                            if(!debug) DTFAG_DIF_MixedRadix << ", Tw1[" << i << "] = " << Tw1[i] << endl;
                            DTFAG_TestPattern_Tw1 << Tw1[i] << endl;
                            //---------------------------------
                       }
                   }else {
                        //-----------for debug-------------------
                        if(debug) DTFAG_DIF_MixedRadix << "v0[" << i << "] = " << v0[i] << ", v1[" << i << "] = " << v1[i];
                        if(debug) Tw1[i] = v0[i] + v1[i];
                        if(debug) DTFAG_DIF_MixedRadix << ", Tw1[" << i << "] = " << Tw1[i] << endl;
                        //----------------------------------------

                        //--------real compute------------
                        if(!debug) DTFAG_DIF_MixedRadix << "v0[" << i << "] = " << v0[i] << ", v1[" << i << "] = " << v1[i];
                        if(!debug) MulMod(Tw1[i], v1[i], v0[i], fft_prime);
                        if(!debug) DTFAG_DIF_MixedRadix << ", Tw1[" << i << "] = " << Tw1[i] << endl;
                        DTFAG_TestPattern_Tw1 << Tw1[i] << endl;
                        //---------------------------------
                   }
                   
                }  

                for(int i=0; i<radix_r1; i++){
                   if(Tw2_display){
                       if(j == Tw_th){
                            //--------real compute------------
                            if(radix_r1 == radix_r2){
                                DTFAG_DIF_MixedRadix << "v0[" << i << "] = " << v0[i];
                                Tw2[i] = v0[i];
                                DTFAG_DIF_MixedRadix << ", Tw2[" << i << "] = " << Tw2[i] << endl;
                                DTFAG_TestPattern_Tw2 << Tw2[i] << endl;
                            }else{
                                DTFAG_DIF_MixedRadix << "v_MixedRadix[" << i << "] = " << v_MixedRadix[i];
                                Tw2[i] = v_MixedRadix[i];
                                DTFAG_DIF_MixedRadix << ", Tw2[" << i << "] = " << Tw2[i] << endl;
                                DTFAG_TestPattern_Tw2 << Tw2[i] << endl;
                            }
                            //---------------------------------
                       }
                   }else {
                        //--------real compute------------
                        if(radix_r1 == radix_r2){
                            DTFAG_DIF_MixedRadix << "v0[" << i << "] = " << v0[i];
                            Tw2[i] = v0[i];
                            DTFAG_DIF_MixedRadix << ", Tw2[" << i << "] = " << Tw2[i] << endl;
                            DTFAG_TestPattern_Tw2 << Tw2[i] << endl;
                        }else{
                            DTFAG_DIF_MixedRadix << "v_MixedRadix[" << i << "] = " << v_MixedRadix[i];
                            Tw2[i] = v_MixedRadix[i];
                            DTFAG_DIF_MixedRadix << ", Tw2[" << i << "] = " << Tw2[i] << endl;
                            DTFAG_TestPattern_Tw2 << Tw2[i] << endl;
                        }
                        //---------------------------------
                   }
                   
                } 
                
                //--------------------------------
            }
        }
    }
}
