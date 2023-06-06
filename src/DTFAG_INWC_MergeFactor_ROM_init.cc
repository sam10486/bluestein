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
    std::ofstream DTFAG_ROM_init("./NWC_PrintData/ROM_Data/DTFAG_INWC_ROM_init.txt");
    std::ofstream DTFAG_ROM0_0_init("./NWC_PrintData/ROM_Data/DTFAG_INWC_ROM0_0_init.txt");
    std::ofstream DTFAG_ROM0_1_init("./NWC_PrintData/ROM_Data/DTFAG_INWC_ROM0_1_init.txt");
    std::ofstream DTFAG_ROM0_2_init("./NWC_PrintData/ROM_Data/DTFAG_INWC_ROM0_2_init.txt");
    std::ofstream DTFAG_ROM0_3_init("./NWC_PrintData/ROM_Data/DTFAG_INWC_ROM0_3_init.txt");

    std::ofstream DTFAG_ROM1_init("./NWC_PrintData/ROM_Data/DTFAG_INWC_ROM1_init.txt");
    std::ofstream DTFAG_ROM2_init("./NWC_PrintData/ROM_Data/DTFAG_INWC_ROM2_init.txt");


    ZZ InvPhi_deg_st0 = PowerMod((ZZ)radix_r1, 0, fft_prime);
    ZZ InvPhi_deg_st1 = PowerMod((ZZ)radix_r1, 1, fft_prime);
    ZZ InvPhi_deg_st2 = PowerMod((ZZ)radix_r1, 2, fft_prime);
    ZZ InvPhi_deg_st3 = PowerMod((ZZ)radix_r1, 3, fft_prime);
	ZZ InvPhi_Order_st0;
    ZZ InvPhi_Order_st1;
    ZZ InvPhi_Order_st2;
    ZZ InvN;
    InvMod(InvN, (ZZ)radix_r1, fft_prime);
    vector<vector<ZZ > > ROM0_st0;
    vector<vector<ZZ > > ROM0_st1;
    vector<vector<ZZ > > ROM0_st2;
    vector<vector<ZZ > > ROM0_st3;
    ROM0_st0.resize(radix_r1);
    for(int i=0; i<radix_r1; i++){
        ROM0_st0[i].resize(radix_r1);
    }
    ROM0_st1.resize(radix_r1);
	for(int i=0; i<radix_r1; i++){
        ROM0_st1[i].resize(radix_r1);
    }
    ROM0_st2.resize(radix_r1);
	for(int i=0; i<radix_r1; i++){
        ROM0_st2[i].resize(radix_r1);
    }
    ROM0_st3.resize(radix_r1);
	for(int i=0; i<radix_r1; i++){
        ROM0_st3[i].resize(radix_r1);
    }

    DTFAG_ROM_init << "****************initial ROM********************" << endl;
    for(int k=0; k<radix_r1; k++){
        DTFAG_ROM_init << "ROM0[" << "k=" << k <<  "] = [";
        for(int n=0; n<radix_r1; n++){
            ZZ ROM0_dg; 
            ROM0[k][n] = (radix_r1 * radix_r2) * k * n * 2;
            ROM0_dg = (radix_r1 * radix_r2) * k * n * 2;
            InvPhi_Order_st0 = PowerMod(InvPhi, InvPhi_deg_st0*n, fft_prime);
            InvPhi_Order_st1 = PowerMod(InvPhi, InvPhi_deg_st1*n, fft_prime);
            InvPhi_Order_st2 = PowerMod(InvPhi, InvPhi_deg_st2*n, fft_prime);
            if(!debug) PowerMod(ROM0[k][n], InvPhi, ROM0_dg, fft_prime);
            //--------------for INWC-------------------
            if(!debug) MulMod(ROM0_st0[k][n], ROM0[k][n], InvPhi_Order_st0, fft_prime); else ROM0[k][n] = ROM0[k][n] + InvPhi_deg_st0*n;
            if(!debug) MulMod(ROM0_st1[k][n], ROM0[k][n], InvPhi_Order_st1, fft_prime); else ROM0[k][n] = ROM0[k][n] + InvPhi_deg_st1*n;
            if(!debug) MulMod(ROM0_st2[k][n], ROM0[k][n], InvPhi_Order_st2, fft_prime); else ROM0[k][n] = ROM0[k][n] + InvPhi_deg_st2*n;
            if(!debug) MulMod(ROM0_st0[k][n], ROM0_st0[k][n], InvN, fft_prime);
            if(!debug) MulMod(ROM0_st1[k][n], ROM0_st1[k][n], InvN, fft_prime);
            if(!debug) MulMod(ROM0_st2[k][n], ROM0_st2[k][n], InvN, fft_prime);
            //-----------------------------------------
            DTFAG_ROM_init << ROM0[k][n] << ", ";
        }
        DTFAG_ROM_init << "]\n";
    }
    
    for(int k=0; k<radix_r1; k++){
        for(int n=0; n<radix_r1; n++){
            switch(stage){
                case 0:
                    ROM0[k][n] = ROM0_st0[k][n];
                    break;
                case 1:
                    ROM0[k][n] = ROM0_st1[k][n];
                    break;
                case 2:
                    ROM0[k][n] = ROM0_st2[k][n];
                    break;
                case 3:
                    ROM0_st3[k][n] = PowerMod(InvPhi, n, fft_prime);
                    ROM0_st3[k][n] = PowerMod(ROM0_st3[k][n], InvPhi_deg_st3, fft_prime);
                    ROM0_st3[k][n] = MulMod(ROM0_st3[k][n], InvN, fft_prime);
                    ROM0[k][n] = ROM0_st3[k][n];
                    //cout << " ROM0[" << k << "][" << n << "] = " << ROM0[k][n] << endl;
                    break;
                default:
                    break;
            }
            DTFAG_ROM0_0_init << "ROM0_st0_arr[" << k << "][" << n << "] <= 64'd" << ROM0_st0[k][n] << ";\n";
            DTFAG_ROM0_1_init << "ROM0_st1_arr[" << k << "][" << n << "] <= 64'd" << ROM0_st1[k][n] << ";\n";
            DTFAG_ROM0_2_init << "ROM0_st2_arr[" << k << "][" << n << "] <= 64'd" << ROM0_st2[k][n] << ";\n";
            DTFAG_ROM0_3_init << "ROM0_st3_arr[" << k << "][" << n << "] <= 64'd" << ROM0_st3[k][n] << ";\n";
        }
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
            DTFAG_ROM1_init << "ROM1_arr[" << k << "][" << n << "] <= 64'd" << ROM1[k][n] << ";\n";
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
            DTFAG_ROM2_init << "ROM2_arr[" << k << "][" << n << "] <= 64'd" << ROM2[k][n] << ";\n";
        }
        DTFAG_ROM_init << "]\n";
    }
    DTFAG_ROM_init << "**************intital fin****************" << endl;
}

