#include <iostream>
#include "NTT.h"
#include "NTTSPMB.h"
#include "BitOperate.h"
#include "DTFAG.h"

using namespace std;

void DTFAG::DTFAG_INWC_MergeFactor_ROM_init(
    int radix_r1, int radix_r2, ZZ fft_twiddle, ZZ fft_prime, int debug, int stage, ZZ InvPhi,
    vector<vector<ZZ > > &ROM0,  vector<vector<ZZ > > &ROM1,  vector<vector<ZZ > > &ROM2){


    BitOperate DecToBin;
    std::ofstream DTFAG_ROM_init("./NWC_PrintData/DTFAG_INWC_ROM_init.txt");

    ZZ InvPhi_deg = PowerMod((ZZ)radix_r1, stage, fft_prime);
	ZZ InvPhi_Order;
    ZZ InvTwo;
    InvMod(InvTwo, (ZZ)2, fft_prime);

    DTFAG_ROM_init << "****************initial ROM********************" << endl;
    for(int k=0; k<radix_r1; k++){
        DTFAG_ROM_init << "ROM0[" << "k=" << k <<  "] = [";
        //cout << "stage" << stage <<  endl;
        //cout << "ROM0[" << "k=" << k <<  "] = [";
        for(int n=0; n<radix_r1; n++){
            ZZ ROM0_dg; 
            switch (stage){
            case 0:
                ROM0[k][n] = (radix_r1 * radix_r2) * k * n * 2;
                ROM0_dg = (radix_r1 * radix_r2) * k * n * 2;
                InvPhi_Order = PowerMod(InvPhi, InvPhi_deg*n, fft_prime);
                if(!debug) PowerMod(ROM0[k][n], InvPhi, ROM0_dg, fft_prime);
                //--------------for INWC-------------------
                if(!debug) MulMod(ROM0[k][n], ROM0[k][n], InvPhi_Order, fft_prime); else ROM0[k][n] = ROM0[k][n] + InvPhi_deg*n;
                if(!debug) MulMod(ROM0[k][n], ROM0[k][n], InvTwo, fft_prime);
                //-----------------------------------------
                DTFAG_ROM_init << ROM0[k][n] << ", ";
                //cout << ROM0[k][n] << ", ";
                break;
            case 1:
                ROM0[k][n] = (radix_r1 * radix_r2) * k * n * 2;
                ROM0_dg = (radix_r1 * radix_r2) * k * n * 2;
                InvPhi_Order = PowerMod(InvPhi, InvPhi_deg*n, fft_prime);
                if(!debug) PowerMod(ROM0[k][n], InvPhi, ROM0_dg, fft_prime);
                //--------------for INWC-------------------
                if(!debug) MulMod(ROM0[k][n], ROM0[k][n], InvPhi_Order, fft_prime); else ROM0[k][n] = ROM0[k][n] + InvPhi_deg*n;
                if(!debug) MulMod(ROM0[k][n], ROM0[k][n], InvTwo, fft_prime);
                //-----------------------------------------
                DTFAG_ROM_init << ROM0[k][n] << ", ";
                //cout << ROM0[k][n] << ", ";
                break;
            case 2:
                ROM0[k][n] = (radix_r1 * radix_r2) * k * n * 2;
                ROM0_dg = (radix_r1 * radix_r2) * k * n * 2;
                InvPhi_Order = PowerMod(InvPhi, InvPhi_deg*n, fft_prime);
                if(!debug) PowerMod(ROM0[k][n], InvPhi, ROM0_dg, fft_prime);
                //--------------for INWC-------------------
                if(!debug) MulMod(ROM0[k][n], ROM0[k][n], InvPhi_Order, fft_prime); else ROM0[k][n] = ROM0[k][n] + InvPhi_deg*n;
                if(!debug) MulMod(ROM0[k][n], ROM0[k][n], InvTwo, fft_prime);
                //-----------------------------------------
                DTFAG_ROM_init << ROM0[k][n] << ", ";
                //cout << ROM0[k][n] << ", ";
                break;
            case 3:
                ROM0[k][n] = (radix_r1 * radix_r2) * k * n * 2;
                ROM0_dg = (radix_r1 * radix_r2) * k * n * 2;
                InvPhi_Order = PowerMod(InvPhi, InvPhi_deg*n, fft_prime);
                if(!debug) PowerMod(ROM0[k][n], InvPhi, ROM0_dg, fft_prime);
                //--------------for INWC-------------------
                if(!debug) MulMod(ROM0[k][n], ROM0[k][n], InvPhi_Order, fft_prime); else ROM0[k][n] = ROM0[k][n] + InvPhi_deg*n;
                if(!debug) MulMod(ROM0[k][n], ROM0[k][n], InvTwo, fft_prime);
                //-----------------------------------------
                DTFAG_ROM_init << ROM0[k][n] << ", ";
                //cout << ROM0[k][n] << ", ";
                break;
            default:
                break;
            }
        }
        DTFAG_ROM_init << "]\n";
        //cout << "]\n";
    }

    DTFAG_ROM_init << "----------------------------------" << endl;
    for(int k=0; k<radix_r2; k++){
        DTFAG_ROM_init << "ROM1[" << "k=" << k <<  "] = [";
        for(int n=0; n<radix_r1; n++){
            ROM1[k][n] = (radix_r1) * k * n * 2;
            ZZ ROM1_dg; 
            ROM1_dg = (radix_r1) * k * n * 2;
            if(!debug) PowerMod(ROM1[k][n], InvPhi, ROM1_dg, fft_prime);
            DTFAG_ROM_init << ROM1[k][n] << ", ";
        }
        DTFAG_ROM_init << "]\n";
    }

    DTFAG_ROM_init << "----------------------------------" << endl;
    for(int k=0; k<radix_r1; k++){
        DTFAG_ROM_init << "ROM2[" << "k=" << k <<  "] = [";
        for(int n=0; n<radix_r1; n++){
            ROM2[k][n] = k * n * 2;
            ZZ ROM2_dg; 
            ROM2_dg = k * n * 2;
            if(!debug) PowerMod(ROM2[k][n], InvPhi, ROM2_dg, fft_prime);
            DTFAG_ROM_init << ROM2[k][n] << ", ";
        }
        DTFAG_ROM_init << "]\n";
    }
    DTFAG_ROM_init << "**************intital fin****************" << endl;


}

