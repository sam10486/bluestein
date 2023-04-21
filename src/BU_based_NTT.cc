#include <iostream>
#include "NTT.h"
#include "NTTSPMB.h"
#include "SPMB.h"
#include "BitOperate.h"
#include "DTFAG.h"
#include "assert.h"
#include <NTL/tools.h>
#include <NTL/ZZ.h>
#include <math.h>


using namespace std;


void DIT(vector<ZZ > &DFT_data, vector<ZZ > &data_in, long long n, ZZ prou, ZZ modular);

void BU_based_NTT(){
    long long fft_point = 64;
    long long difference_length = 65536 / fft_point;
    long long difference_16     = fft_point / 16;
    long long band_memory_size  = fft_point / 32;
    ZZ fft_prime;
    ZZ fft_twiddle_65536;
    ZZ fft_twiddle;
    conv(fft_prime,"18446744069414584321");
    conv(fft_twiddle_65536,"14603442835287214144");  //65536-th twiddle factor

    conv(fft_prime,"197");
    conv(fft_twiddle_65536,"8");  //65536-th twiddle factor

    PowerMod(fft_twiddle,fft_twiddle_65536,difference_length,fft_prime);
    cout << "fft_twiddle = " << fft_twiddle << endl;
    //--------NTT golden----------
    NTT     NTT_test;
    std::vector<ZZ> A_1;
    A_1.resize(fft_point);
    for(int i = 0;i < fft_point;i++){
		  A_1[i] = i;
    }
    NTT_test.NTT_init(fft_point,fft_prime,fft_twiddle);
    NTT_test.NTT_t(A_1);
    
    //----------------------------
    

    vector<ZZ > data_in;
    vector<ZZ > data_out;
    data_in.resize(fft_point);
    data_out.resize(fft_point);
    BitOperate BR;
    

    for (int i = 0; i < fft_point; i++)
    {
        int idx_BR = BR.BitReserve(i, log2(fft_point));
        //data_in.at(i) = idx_BR;
        data_in.at(i) = i;
    }
    DIT(data_out, data_in, fft_point, fft_twiddle, fft_prime);

    for (int i = 0; i < fft_point; i++)
    {
        //cout << "data_out[" << i << "] = " << data_out[i] << endl;
    }

    int err=0;
    for (int i = 0; i < fft_point; i++)
    {
        if (data_out[i] != A_1[i])
        {
            err++;
            cout << "data_out[" << i << "] = " << data_out[i] << endl; 
        }else{
            //cout << "data_out "<< i << " pass" << endl;
        }
        
    }
    cout << "err = " << err << endl;
    

}



void DIT(vector<ZZ > &NWC_ans, vector<ZZ > &data_in, long long n, ZZ prou, ZZ modular)
{
    std::ofstream BU_based_NTT("./BU_based_NTT.txt");
    vector<ZZ > tmp_arr;
    tmp_arr.resize(n);
    BitOperate BR;
    ZZ up_correct_tw, down_correct_tw;
    ZZ correct_tw;
    int debug = 0;
    long long bit_width = log2(n);
    int radix = 4;
    int stage = log(n)/log(radix);
    for(long long i=0; i<n; i++){
        tmp_arr[i] = data_in[i];
    }
    
    for (int s = log2(n); s > 0 ; s--)
    {
        BU_based_NTT << "stage = " << s << endl;
        long long m = pow(2, s);
        for (int k = 0; k < (n/m); k++)
        {
           long long idx_BR = BR.BitReserve(k, log2(n)-s);
           long long twiddle_dg = idx_BR * (m/2);
           ZZ twiddle; 
           if(!debug) PowerMod(twiddle, prou, twiddle_dg, modular);
           //----------debug-------------
           if(debug) twiddle = twiddle_dg;
           //----------------------------
           BU_based_NTT << "-----------------------------------------------------"<< endl;
           BU_based_NTT << "    k = " << k << endl;
           for (int j = 0; j < m/2; j++)
           {
                BU_based_NTT << "-----------------"<< endl;
                BU_based_NTT << "   j = " << j << endl;
                ZZ u = tmp_arr[k*m+j];
                ZZ t ;

                if(!debug) t = MulMod(twiddle, tmp_arr[k*m+j+(m/2)], modular);


                if (ceil(s/2) != stage)
                {
                    if (s % 2 == 0 && j == 1){
                        int br = BR.BitReserve(k, log2((n/m)));
                        //--------debug---------
                        if(debug) up_correct_tw = br;
                        if(debug) down_correct_tw = br;
                        //----------------------
                        if(!debug) PowerMod(up_correct_tw, prou, br, modular);
                        if(!debug) PowerMod(down_correct_tw, prou, br, modular);
                        if(!debug) MulMod(u, u, up_correct_tw, modular);
                        if(!debug) MulMod(t, t, down_correct_tw, modular);
                    }
                    if (s % 2 ==  1 ){    
                        int k_shift = k >> 1;
                        int br = BR.BitReserve(k_shift, log2( n/(2*m) ) );
                        //--------debug---------
                        if(debug) correct_tw = br;
                        //----------------------
                        if(!debug) PowerMod(correct_tw, prou, br, modular);
                         if(!debug) InvMod(correct_tw, correct_tw, modular);
                        if(!debug) MulMod(t, t, correct_tw, modular);
                    }           
                }

                BU_based_NTT << "Before butterfly unit operation!" << endl;
                ZZ tw_up;
                ZZ tw_down;
                switch (n)
                {
                 case 16:
                    //--------------debug---------------
                    if(s == 2 && j==1){
                        if(debug) tw_up = up_correct_tw;
                        if(debug) tw_down = twiddle + down_correct_tw;
                    }else if (s == 1){
                        if(debug) tw_up = 0;
                        if(debug) tw_down = twiddle - correct_tw;
                    }else{
                        if(debug) tw_up = 0;
                        if(debug) tw_down = twiddle;
                    }
                    //----------------------------------
                    if(!debug) tw_up = up_correct_tw;
                    if(s == 2){
                        if(!debug) MulMod(tw_down,twiddle, down_correct_tw, modular);
                    }else if (s == 1){
                        if(!debug) MulMod(tw_down,twiddle, correct_tw, modular);
                    }
                    BU_based_NTT << "    idx_up[" << k*m+j << "] = " << u << ", tw = " << tw_up << endl;
                    BU_based_NTT << "    idx_down[" << k*m+j+(m/2) << "] = " << t << ", tw = " << tw_down << endl;
                    break;
                 case 32:
                    //--------------debug---------------
                    if(s == 2 && j==1){
                        if(debug) tw_up = up_correct_tw;
                        if(debug) tw_down = twiddle + down_correct_tw;
                    }else if (s == 1){
                        if(debug) tw_up = 0;
                        if(debug) tw_down = twiddle - correct_tw;
                    }else{
                        if(debug) tw_up = 0;
                        if(debug) tw_down = twiddle;
                    }
                    //----------------------------------
                    if(!debug) tw_up = up_correct_tw;
                    if(s == 2){
                        if(!debug) MulMod(tw_down,twiddle, down_correct_tw, modular);
                    }else if (s == 1){
                        if(!debug) MulMod(tw_down,twiddle, correct_tw, modular);
                    }
                    BU_based_NTT << "    idx_up[" << k*m+j << "] = " << u << ", tw = " << tw_up << endl;
                    BU_based_NTT << "    idx_down[" << k*m+j+(m/2) << "] = " << t << ", tw = " << tw_down << endl;
                    break;
                 default:
                    break;
               }
                //--------radix 2 BU---------------
                tmp_arr[k*m+j] = AddMod(u, t, modular);
                tmp_arr[k*m+j+(m/2)] = SubMod(u, t, modular);
                //---------------------------------

                BU_based_NTT << "After butterfly unit operation!" << endl;
                switch (n)
                {
                 case 16:
                    BU_based_NTT << "    idx_up[" << k*m+j << "] = " <<  tmp_arr[k*m+j] << endl;
                    BU_based_NTT << "    idx_down[" << k*m+j+(m/2) << "] = " << tmp_arr[k*m+j+(m/2)] << endl;
                    break;
                 case 32:
                     BU_based_NTT << "    idx_up[" << k*m+j << "] = " <<  tmp_arr[k*m+j] << endl;
                     BU_based_NTT << "    idx_down[" << k*m+j+(m/2) << "] = " << tmp_arr[k*m+j+(m/2)] << endl;
                     break;
                 default:
                    break;
               }
           }
        }
    }
    for (int i = 0; i < n; i++)
    {   
        int idx_BR = BR.BitReserve(i, log2(n));
        NWC_ans[i] = tmp_arr[idx_BR];
    }
    

}
