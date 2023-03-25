#include <iostream>
#include "NTT.h"
#include "NTTSPMB.h"
#include "BitOperate.h"
#include "DTFAG.h"

using namespace std;

void DTFAG::DTFAG_ROM_init(
    int radix_r1, int radix_r2, ZZ fft_twiddle, ZZ fft_prime, int debug,
    vector<vector<ZZ > > &ROM0,  vector<ZZ > &ROM1,  vector<ZZ > &ROM2){

    std::ofstream DTFAG_ROM_init("./my_print_data/DTFAG_ROM_init.txt");

    DTFAG_ROM_init << "****************initial ROM********************" << endl;
    for(int k=0; k<radix_r1; k++){
        DTFAG_ROM_init << "ROM0[k=" << k <<  "] = [";
        for(int n=0; n<radix_r1; n++){
            ZZ ROM0_dg; 
            ROM0_dg = (radix_r1 * radix_r1) * n * k;
            PowerMod(ROM0[k][n], fft_twiddle, ROM0_dg, fft_prime);
            //------------------for debug-------------------
            if(debug) ROM0[k][n] = (radix_r1 * radix_r1) * k * n;
            //----------------------------------------------
            DTFAG_ROM_init << ROM0[k][n] << ", ";
        }
        DTFAG_ROM_init << "]\n";
    }

    DTFAG_ROM_init << "----------------------------------" << endl;
    for(int group=0; group<radix_r1; group++){
        DTFAG_ROM_init << "ROM1[g=" << group <<  "] = [";
        for(int n=0; n<radix_r1; n++){
            ZZ ROM1_dg;
            ROM1_dg = group * radix_r1 * n;
            PowerMod(ROM1[group*radix_r1+n], fft_twiddle, ROM1_dg, fft_prime);
            //------------------for debug-------------------
            if(debug) ROM1[group*radix_r1+n] = group * radix_r1 * n;
            //----------------------------------------------
            DTFAG_ROM_init <<  ROM1[group*radix_r1+n] << ", ";
        }
        DTFAG_ROM_init << "]\n";
    }

    DTFAG_ROM_init << "----------------------------------" << endl;
    for(int group=0; group<radix_r1; group++){
        DTFAG_ROM_init << "ROM2[g=" << group <<  "] = [";
        for(int n=0; n<radix_r1; n++){
            ZZ ROM2_dg;
            ROM2_dg = group * n;
            PowerMod(ROM2[group*radix_r1+n], fft_twiddle, ROM2_dg, fft_prime);
            //------------------for debug-------------------
            if(debug) ROM2[group*radix_r1+n] = group * n;
            //----------------------------------------------
            DTFAG_ROM_init << ROM2[group*radix_r1+n] << ", ";
        }
        DTFAG_ROM_init << "]\n";
    }
    DTFAG_ROM_init << "**************intital fin****************" << endl;
}

