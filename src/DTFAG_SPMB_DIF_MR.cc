#include <iostream>
#include "NTT.h"
#include "NTTSPMB.h"
#include "BitOperate.h"
#include "DTFAG.h"

using namespace std;

void DTFAG::DTFAG_SPMB_DIF_MR (
        int stage, int fft_point, int radix_r1, int radix_r2, int debug,
        vector<vector<ZZ > > &ROM0,  vector<vector<ZZ > > &ROM1,  vector<vector<ZZ > > &ROM2,
        vector<ZZ > &st0_Tw, vector<ZZ > &st1_Tw, vector<ZZ > &st2_Tw,
        int DTFAG_i, int DTFAG_t, int DTFAG_j
) {
    //------- radix sel-----
    ZZ fft_prime ;
    ZZ fft_twiddle_65536 ;
    ZZ fft_twiddle ;
    long difference_length = 65536 / fft_point ;
    conv(fft_prime, "18446744069414584321"); // prime number
    conv(fft_twiddle_65536, "14603442835287214144"); // twiddle factor based setting by main.cc
    
    PowerMod(fft_twiddle,fft_twiddle_65536,difference_length,fft_prime);
    //-----------------------

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
    std::ofstream DTFAG_SPMB_DIF_MR("./my_print_data/DTFAG_SPMB_DIF_MR.txt");


    int Tw2_display = 0;
    int Tw1_display = 0;
    int Tw0_display = false;
    int Tw_th = 1;
    
    DTFAG_SPMB_DIF_MR << "     DTFAG_i = " << DTFAG_i << endl;
    DTFAG_SPMB_DIF_MR << "-----------------crossbar for len--------------------"<< endl;
    DTFAG_SPMB_DIF_MR << "     DTFAG_t = " << DTFAG_t << endl;

    
    DTFAG_SPMB_DIF_MR << "----------------------------" << endl;
    if(Tw2_display || Tw1_display || Tw0_display){
        if(DTFAG_j == Tw_th) DTFAG_SPMB_DIF_MR << "       DTFAG_j = " << DTFAG_j << endl;
    }else{
        DTFAG_SPMB_DIF_MR << "      DTFAG_j = " << DTFAG_j << " (TF" << DTFAG_j << ")"<< endl;
    }
    
    MA0 = DTFAG_j;
    DTFAG_SPMB_DIF_MR << "      MA0 = " << MA0 << endl;
    MA1 = NTTSPMB.Gray(DTFAG_i,radix_r1);
    DTFAG_SPMB_DIF_MR << "      MA1 = " << MA1 << endl;
    if(DTFAG_i % 2 == 1){
        int t_complement = number_complement.number_complement(DTFAG_t, radix_r1);
        MA2 = NTTSPMB.Gray(t_complement,radix_r1);
        DTFAG_SPMB_DIF_MR << "      MA2 : " << MA2 << " = " << "G(" << t_complement << ")"<< endl; 
    }else{
        MA2 = NTTSPMB.Gray(DTFAG_t,radix_r1);
        DTFAG_SPMB_DIF_MR << "      MA2 : " << MA2 << " = " << "G(" << DTFAG_t << ")"<< endl; 
    }        
    if(radix_r1 == radix_r2) {
        MA_MixedRadix = DTFAG_j;
        DTFAG_SPMB_DIF_MR << "      MA_MixedRadix = " << MA_MixedRadix << endl;
    }else {
        MA_MixedRadix =  (DTFAG_j * (radix_r1/radix_r2)) % radix_r1;
        DTFAG_SPMB_DIF_MR << "      MA_MixedRadix = " << MA_MixedRadix << endl;
    }
    
    for(int i=0; i<radix_r1; i++){
        v0[i] = ROM0[MA0][i];
    }
    //if(stage==1) cout << "MA1 = " << MA1 << endl;
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
        //-----------for debug-------------------
        if(debug) DTFAG_SPMB_DIF_MR << "v0[" << i << "] = " << v0[i] << ", v1[" << i << "] = " << v1[i] << ", v2[" << i << "] = " << v2[i];
        if(debug) Tw0[i] = v0[i] + v1[i] + v2[i];
        if(debug) DTFAG_SPMB_DIF_MR << ", Tw0[" << i << "] = " << Tw0[i] << endl;
        //----------------------------------------
        //--------real compute------------
        ZZ tmp;
        if(!debug) DTFAG_SPMB_DIF_MR << "v0[" << i << "] = " << v0[i] << ", v1[" << i << "] = " << v1[i] << ", v2[" << i << "] = " << v2[i];
        if(!debug) MulMod(tmp, v1[i], v2[i], fft_prime);
        if(!debug) MulMod(Tw0[i], v0[i], tmp, fft_prime);
        if(!debug) DTFAG_SPMB_DIF_MR << ", Tw0[" << i << "] = " << Tw0[i] << endl;
        //---------------------------------
        if(stage == 0){
            st0_Tw[i] = Tw0[i];
        }
    }
      
    for(int i=0; i<radix_r1; i++){
        //-----------for debug-------------------
        if(debug) DTFAG_SPMB_DIF_MR << "v0[" << i << "] = " << v0[i] << ", v1[" << i << "] = " << v1[i];
        if(debug) Tw1[i] = v0[i] + v1[i];
        if(debug) DTFAG_SPMB_DIF_MR << ", Tw1[" << i << "] = " << Tw1[i] << endl;
        //----------------------------------------
        //--------real compute------------
        if(!debug) DTFAG_SPMB_DIF_MR << "v0[" << i << "] = " << v0[i] << ", v1[" << i << "] = " << v1[i];
        if(!debug) MulMod(Tw1[i], v1[i], v0[i], fft_prime);
        if(!debug) DTFAG_SPMB_DIF_MR << ", Tw1[" << i << "] = " << Tw1[i] << endl;
        //---------------------------------
        if(stage == 1){
            st1_Tw[i] = Tw1[i];
        }
    }
    for(int i=0; i<radix_r1; i++){
        //--------real compute------------
        if(radix_r1 == radix_r2){
            DTFAG_SPMB_DIF_MR << "v0[" << i << "] = " << v0[i];
            Tw2[i] = v0[i];
            DTFAG_SPMB_DIF_MR << ", Tw2[" << i << "] = " << Tw2[i] << endl;
        }else{
            DTFAG_SPMB_DIF_MR << "v_MixedRadix[" << i << "] = " << v_MixedRadix[i];
            Tw2[i] = v_MixedRadix[i];
            DTFAG_SPMB_DIF_MR << ", Tw2[" << i << "] = " << Tw2[i] << endl;
        }
        //---------------------------------
        if(stage == 2){
            st2_Tw[i] = Tw2[i];
        }
    }           
}
