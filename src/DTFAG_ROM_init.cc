#include <iostream>
#include "NTT.h"
#include "NTTSPMB.h"
#include "BitOperate.h"
#include "DTFAG.h"

using namespace std;

void DTFAG::DTFAG_ROM_init(
    int radix_r1, int radix_r2, ZZ fft_twiddle, ZZ fft_prime, int debug,
    vector<vector<ZZ > > &ROM0,  vector<vector<ZZ > > &ROM1,  vector<vector<ZZ > > &ROM2){


    BitOperate DecToBin;
    std::ofstream DTFAG_ROM_init("./my_print_data/DTFAG_ROM_init.txt");
    //---------------ZZ to bin-----------------
    // ROM0
    std::ofstream Bin_DTFAG_DIF_ROM0_B0("./ROM_Data/ROM0/Bin_DTFAG_DIF_ROM0_B0.txt");
    std::ofstream Bin_DTFAG_DIF_ROM0_B1("./ROM_Data/ROM0/Bin_DTFAG_DIF_ROM0_B1.txt");
    std::ofstream Bin_DTFAG_DIF_ROM0_B2("./ROM_Data/ROM0/Bin_DTFAG_DIF_ROM0_B2.txt");
    std::ofstream Bin_DTFAG_DIF_ROM0_B3("./ROM_Data/ROM0/Bin_DTFAG_DIF_ROM0_B3.txt");
    std::ofstream Bin_DTFAG_DIF_ROM0_B4("./ROM_Data/ROM0/Bin_DTFAG_DIF_ROM0_B4.txt");
    std::ofstream Bin_DTFAG_DIF_ROM0_B5("./ROM_Data/ROM0/Bin_DTFAG_DIF_ROM0_B5.txt");
    std::ofstream Bin_DTFAG_DIF_ROM0_B6("./ROM_Data/ROM0/Bin_DTFAG_DIF_ROM0_B6.txt");
    std::ofstream Bin_DTFAG_DIF_ROM0_B7("./ROM_Data/ROM0/Bin_DTFAG_DIF_ROM0_B7.txt");
    // ROM1
    std::ofstream Bin_DTFAG_DIF_ROM1_B0("./ROM_Data/ROM1/Bin_DTFAG_DIF_ROM1_B0.txt");
    std::ofstream Bin_DTFAG_DIF_ROM1_B1("./ROM_Data/ROM1/Bin_DTFAG_DIF_ROM1_B1.txt");
    std::ofstream Bin_DTFAG_DIF_ROM1_B2("./ROM_Data/ROM1/Bin_DTFAG_DIF_ROM1_B2.txt");
    std::ofstream Bin_DTFAG_DIF_ROM1_B3("./ROM_Data/ROM1/Bin_DTFAG_DIF_ROM1_B3.txt");
    std::ofstream Bin_DTFAG_DIF_ROM1_B4("./ROM_Data/ROM1/Bin_DTFAG_DIF_ROM1_B4.txt");
    std::ofstream Bin_DTFAG_DIF_ROM1_B5("./ROM_Data/ROM1/Bin_DTFAG_DIF_ROM1_B5.txt");
    std::ofstream Bin_DTFAG_DIF_ROM1_B6("./ROM_Data/ROM1/Bin_DTFAG_DIF_ROM1_B6.txt");
    std::ofstream Bin_DTFAG_DIF_ROM1_B7("./ROM_Data/ROM1/Bin_DTFAG_DIF_ROM1_B7.txt");
    // ROM2
    std::ofstream Bin_DTFAG_DIF_ROM2_B0("./ROM_Data/ROM2/Bin_DTFAG_DIF_ROM2_B0.txt");
    std::ofstream Bin_DTFAG_DIF_ROM2_B1("./ROM_Data/ROM2/Bin_DTFAG_DIF_ROM2_B1.txt");
    std::ofstream Bin_DTFAG_DIF_ROM2_B2("./ROM_Data/ROM2/Bin_DTFAG_DIF_ROM2_B2.txt");
    std::ofstream Bin_DTFAG_DIF_ROM2_B3("./ROM_Data/ROM2/Bin_DTFAG_DIF_ROM2_B3.txt");
    std::ofstream Bin_DTFAG_DIF_ROM2_B4("./ROM_Data/ROM2/Bin_DTFAG_DIF_ROM2_B4.txt");
    std::ofstream Bin_DTFAG_DIF_ROM2_B5("./ROM_Data/ROM2/Bin_DTFAG_DIF_ROM2_B5.txt");
    std::ofstream Bin_DTFAG_DIF_ROM2_B6("./ROM_Data/ROM2/Bin_DTFAG_DIF_ROM2_B6.txt");
    std::ofstream Bin_DTFAG_DIF_ROM2_B7("./ROM_Data/ROM2/Bin_DTFAG_DIF_ROM2_B7.txt");
    //-----------------------------------------

    DTFAG_ROM_init << "****************initial ROM********************" << endl;
    for(int k=0; k<radix_r1; k++){
        DTFAG_ROM_init << "ROM0[" << "k=" << k <<  "] = [";
        for(int n=0; n<radix_r1; n++){
            ROM0[k][n] = (radix_r1 * radix_r2) * k * n;
            ZZ ROM0_dg; 
            ROM0_dg = (radix_r1 * radix_r2) * k * n;
            if(!debug) PowerMod(ROM0[k][n], fft_twiddle, ROM0_dg, fft_prime);
            DTFAG_ROM_init << ROM0[k][n] << ", ";
        }
        DTFAG_ROM_init << "]\n";
    }

    DTFAG_ROM_init << "----------------------------------" << endl;
    for(int k=0; k<radix_r2; k++){
        DTFAG_ROM_init << "ROM1[" << "k=" << k <<  "] = [";
        for(int n=0; n<radix_r1; n++){
            ROM1[k][n] = (radix_r1) * k * n;
            ZZ ROM1_dg; 
            ROM1_dg = (radix_r1) * k * n;
            if(!debug) PowerMod(ROM1[k][n], fft_twiddle, ROM1_dg, fft_prime);
            DTFAG_ROM_init << ROM1[k][n] << ", ";
        }
        DTFAG_ROM_init << "]\n";
    }

    DTFAG_ROM_init << "----------------------------------" << endl;
    for(int k=0; k<radix_r1; k++){
        DTFAG_ROM_init << "ROM2[" << "k=" << k <<  "] = [";
        for(int n=0; n<radix_r1; n++){
            ROM2[k][n] = k * n;
            ZZ ROM2_dg; 
            ROM2_dg = k * n;
            if(!debug) PowerMod(ROM2[k][n], fft_twiddle, ROM2_dg, fft_prime);
            DTFAG_ROM_init << ROM2[k][n] << ", ";
        }
        DTFAG_ROM_init << "]\n";
    }
    DTFAG_ROM_init << "**************intital fin****************" << endl;


    // -----------------------ZZ To Bin-------------------------
    // ROM0
    for(int k=0; k<radix_r1; k++){
        vector<ZZ > ROM0_D1;
        vector<ZZ > ROM0_D2;
        vector<ZZ > ROM0_D3;
        vector<ZZ > ROM0_D4;
        vector<ZZ > ROM0_D5;
        vector<ZZ > ROM0_D6;
        vector<ZZ > ROM0_D7;
        vector<ZZ > ROM0_D8;
        vector<ZZ > ROM0_D9;
        vector<ZZ > ROM0_D10;
        vector<ZZ > ROM0_D11;
        vector<ZZ > ROM0_D12;
        vector<ZZ > ROM0_D13;
        vector<ZZ > ROM0_D14;
        vector<ZZ > ROM0_D15;
        int bit_width = 64;
        ROM0_D1.resize(bit_width);
        ROM0_D2.resize(bit_width);
        ROM0_D3.resize(bit_width);
        ROM0_D4.resize(bit_width);
        ROM0_D5.resize(bit_width);
        ROM0_D6.resize(bit_width);
        ROM0_D7.resize(bit_width);
        ROM0_D8.resize(bit_width);
        ROM0_D9.resize(bit_width);
        ROM0_D10.resize(bit_width);
        ROM0_D11.resize(bit_width);
        ROM0_D12.resize(bit_width);
        ROM0_D13.resize(bit_width);
        ROM0_D14.resize(bit_width);
        ROM0_D15.resize(bit_width);
        ROM0_D1 = DecToBin.ZZ_DecToBin(ROM0[k][1], bit_width);
        ROM0_D2 = DecToBin.ZZ_DecToBin(ROM0[k][2], bit_width);
        ROM0_D3 = DecToBin.ZZ_DecToBin(ROM0[k][3], bit_width);
        ROM0_D4 = DecToBin.ZZ_DecToBin(ROM0[k][4], bit_width);
        ROM0_D5 = DecToBin.ZZ_DecToBin(ROM0[k][5], bit_width);
        ROM0_D6 = DecToBin.ZZ_DecToBin(ROM0[k][6], bit_width);
        ROM0_D7 = DecToBin.ZZ_DecToBin(ROM0[k][7], bit_width);
        ROM0_D8 = DecToBin.ZZ_DecToBin(ROM0[k][8], bit_width);
        ROM0_D9 = DecToBin.ZZ_DecToBin(ROM0[k][9], bit_width);
        ROM0_D10 = DecToBin.ZZ_DecToBin(ROM0[k][10], bit_width);
        ROM0_D11 = DecToBin.ZZ_DecToBin(ROM0[k][11], bit_width);
        ROM0_D12 = DecToBin.ZZ_DecToBin(ROM0[k][12], bit_width);
        ROM0_D13 = DecToBin.ZZ_DecToBin(ROM0[k][13], bit_width);
        ROM0_D14 = DecToBin.ZZ_DecToBin(ROM0[k][14], bit_width);
        ROM0_D15 = DecToBin.ZZ_DecToBin(ROM0[k][15], bit_width);
        for (int i = 0; i < bit_width; i++){
            // 64 bits
            Bin_DTFAG_DIF_ROM0_B0 << ROM0_D1[bit_width-1-i];
            // 128 bits
            Bin_DTFAG_DIF_ROM0_B1 << ROM0_D2[bit_width-1-i];
            Bin_DTFAG_DIF_ROM0_B2 << ROM0_D4[bit_width-1-i];
            Bin_DTFAG_DIF_ROM0_B3 << ROM0_D6[bit_width-1-i];
            Bin_DTFAG_DIF_ROM0_B4 << ROM0_D8[bit_width-1-i];
            Bin_DTFAG_DIF_ROM0_B5 << ROM0_D10[bit_width-1-i];
            Bin_DTFAG_DIF_ROM0_B6 << ROM0_D12[bit_width-1-i];
            Bin_DTFAG_DIF_ROM0_B7 << ROM0_D14[bit_width-1-i];
        }
        Bin_DTFAG_DIF_ROM0_B1 << "_";
        Bin_DTFAG_DIF_ROM0_B2 << "_";
        Bin_DTFAG_DIF_ROM0_B3 << "_";
        Bin_DTFAG_DIF_ROM0_B4 << "_";
        Bin_DTFAG_DIF_ROM0_B5 << "_";
        Bin_DTFAG_DIF_ROM0_B6 << "_";
        Bin_DTFAG_DIF_ROM0_B7 << "_";
        for (int i = 0; i < bit_width; i++){
            // 128 bits
            Bin_DTFAG_DIF_ROM0_B1 << ROM0_D3[bit_width-1-i];
            Bin_DTFAG_DIF_ROM0_B2 << ROM0_D5[bit_width-1-i];
            Bin_DTFAG_DIF_ROM0_B3 << ROM0_D7[bit_width-1-i];
            Bin_DTFAG_DIF_ROM0_B4 << ROM0_D9[bit_width-1-i];
            Bin_DTFAG_DIF_ROM0_B5 << ROM0_D11[bit_width-1-i];
            Bin_DTFAG_DIF_ROM0_B6 << ROM0_D13[bit_width-1-i];
            Bin_DTFAG_DIF_ROM0_B7 << ROM0_D15[bit_width-1-i];
        }
        Bin_DTFAG_DIF_ROM0_B0 << endl;
        Bin_DTFAG_DIF_ROM0_B1 << endl;
        Bin_DTFAG_DIF_ROM0_B2 << endl;
        Bin_DTFAG_DIF_ROM0_B3 << endl;
        Bin_DTFAG_DIF_ROM0_B4 << endl;
        Bin_DTFAG_DIF_ROM0_B5 << endl;
        Bin_DTFAG_DIF_ROM0_B6 << endl;
        Bin_DTFAG_DIF_ROM0_B7 << endl;
    }
    
    // ROM1
    for(int k=0; k<radix_r1; k++){
        vector<ZZ > ROM1_D1;
        vector<ZZ > ROM1_D2;
        vector<ZZ > ROM1_D3;
        vector<ZZ > ROM1_D4;
        vector<ZZ > ROM1_D5;
        vector<ZZ > ROM1_D6;
        vector<ZZ > ROM1_D7;
        vector<ZZ > ROM1_D8;
        vector<ZZ > ROM1_D9;
        vector<ZZ > ROM1_D10;
        vector<ZZ > ROM1_D11;
        vector<ZZ > ROM1_D12;
        vector<ZZ > ROM1_D13;
        vector<ZZ > ROM1_D14;
        vector<ZZ > ROM1_D15;
        int bit_width = 64;
        ROM1_D1.resize(bit_width);
        ROM1_D2.resize(bit_width);
        ROM1_D3.resize(bit_width);
        ROM1_D4.resize(bit_width);
        ROM1_D5.resize(bit_width);
        ROM1_D6.resize(bit_width);
        ROM1_D7.resize(bit_width);
        ROM1_D8.resize(bit_width);
        ROM1_D9.resize(bit_width);
        ROM1_D10.resize(bit_width);
        ROM1_D11.resize(bit_width);
        ROM1_D12.resize(bit_width);
        ROM1_D13.resize(bit_width);
        ROM1_D14.resize(bit_width);
        ROM1_D15.resize(bit_width);
        ROM1_D1 = DecToBin.ZZ_DecToBin(ROM1[k][1], bit_width);
        ROM1_D2 = DecToBin.ZZ_DecToBin(ROM1[k][2], bit_width);
        ROM1_D3 = DecToBin.ZZ_DecToBin(ROM1[k][3], bit_width);
        ROM1_D4 = DecToBin.ZZ_DecToBin(ROM1[k][4], bit_width);
        ROM1_D5 = DecToBin.ZZ_DecToBin(ROM1[k][5], bit_width);
        ROM1_D6 = DecToBin.ZZ_DecToBin(ROM1[k][6], bit_width);
        ROM1_D7 = DecToBin.ZZ_DecToBin(ROM1[k][7], bit_width);
        ROM1_D8 = DecToBin.ZZ_DecToBin(ROM1[k][8], bit_width);
        ROM1_D9 = DecToBin.ZZ_DecToBin(ROM1[k][9], bit_width);
        ROM1_D10 = DecToBin.ZZ_DecToBin(ROM1[k][10], bit_width);
        ROM1_D11 = DecToBin.ZZ_DecToBin(ROM1[k][11], bit_width);
        ROM1_D12 = DecToBin.ZZ_DecToBin(ROM1[k][12], bit_width);
        ROM1_D13 = DecToBin.ZZ_DecToBin(ROM1[k][13], bit_width);
        ROM1_D14 = DecToBin.ZZ_DecToBin(ROM1[k][14], bit_width);
        ROM1_D15 = DecToBin.ZZ_DecToBin(ROM1[k][15], bit_width);
        for (int i = 0; i < bit_width; i++){
            // 64 bits
            Bin_DTFAG_DIF_ROM1_B0 << ROM1_D1[bit_width-1-i];
            // 128 bits
            Bin_DTFAG_DIF_ROM1_B1 << ROM1_D2[bit_width-1-i];
            Bin_DTFAG_DIF_ROM1_B2 << ROM1_D4[bit_width-1-i];
            Bin_DTFAG_DIF_ROM1_B3 << ROM1_D6[bit_width-1-i];
            Bin_DTFAG_DIF_ROM1_B4 << ROM1_D8[bit_width-1-i];
            Bin_DTFAG_DIF_ROM1_B5 << ROM1_D10[bit_width-1-i];
            Bin_DTFAG_DIF_ROM1_B6 << ROM1_D12[bit_width-1-i];
            Bin_DTFAG_DIF_ROM1_B7 << ROM1_D14[bit_width-1-i];
        }
        Bin_DTFAG_DIF_ROM1_B1 << "_";
        Bin_DTFAG_DIF_ROM1_B2 << "_";
        Bin_DTFAG_DIF_ROM1_B3 << "_";
        Bin_DTFAG_DIF_ROM1_B4 << "_";
        Bin_DTFAG_DIF_ROM1_B5 << "_";
        Bin_DTFAG_DIF_ROM1_B6 << "_";
        Bin_DTFAG_DIF_ROM1_B7 << "_";
        for (int i = 0; i < bit_width; i++){
            // 128 bits
            Bin_DTFAG_DIF_ROM1_B1 << ROM1_D3[bit_width-1-i];
            Bin_DTFAG_DIF_ROM1_B2 << ROM1_D5[bit_width-1-i];
            Bin_DTFAG_DIF_ROM1_B3 << ROM1_D7[bit_width-1-i];
            Bin_DTFAG_DIF_ROM1_B4 << ROM1_D9[bit_width-1-i];
            Bin_DTFAG_DIF_ROM1_B5 << ROM1_D11[bit_width-1-i];
            Bin_DTFAG_DIF_ROM1_B6 << ROM1_D13[bit_width-1-i];
            Bin_DTFAG_DIF_ROM1_B7 << ROM1_D15[bit_width-1-i];
        }
        Bin_DTFAG_DIF_ROM1_B0 << endl;
        Bin_DTFAG_DIF_ROM1_B1 << endl;
        Bin_DTFAG_DIF_ROM1_B2 << endl;
        Bin_DTFAG_DIF_ROM1_B3 << endl;
        Bin_DTFAG_DIF_ROM1_B4 << endl;
        Bin_DTFAG_DIF_ROM1_B5 << endl;
        Bin_DTFAG_DIF_ROM1_B6 << endl;
        Bin_DTFAG_DIF_ROM1_B7 << endl;
    }
    // ROM2
    for(int k=0; k<radix_r1; k++){
        vector<ZZ > ROM2_D1;
        vector<ZZ > ROM2_D2;
        vector<ZZ > ROM2_D3;
        vector<ZZ > ROM2_D4;
        vector<ZZ > ROM2_D5;
        vector<ZZ > ROM2_D6;
        vector<ZZ > ROM2_D7;
        vector<ZZ > ROM2_D8;
        vector<ZZ > ROM2_D9;
        vector<ZZ > ROM2_D10;
        vector<ZZ > ROM2_D11;
        vector<ZZ > ROM2_D12;
        vector<ZZ > ROM2_D13;
        vector<ZZ > ROM2_D14;
        vector<ZZ > ROM2_D15;
        int bit_width = 64;
        ROM2_D1.resize(bit_width);
        ROM2_D2.resize(bit_width);
        ROM2_D3.resize(bit_width);
        ROM2_D4.resize(bit_width);
        ROM2_D5.resize(bit_width);
        ROM2_D6.resize(bit_width);
        ROM2_D7.resize(bit_width);
        ROM2_D8.resize(bit_width);
        ROM2_D9.resize(bit_width);
        ROM2_D10.resize(bit_width);
        ROM2_D11.resize(bit_width);
        ROM2_D12.resize(bit_width);
        ROM2_D13.resize(bit_width);
        ROM2_D14.resize(bit_width);
        ROM2_D15.resize(bit_width);
        ROM2_D1 = DecToBin.ZZ_DecToBin(ROM2[k][1], bit_width);
        ROM2_D2 = DecToBin.ZZ_DecToBin(ROM2[k][2], bit_width);
        ROM2_D3 = DecToBin.ZZ_DecToBin(ROM2[k][3], bit_width);
        ROM2_D4 = DecToBin.ZZ_DecToBin(ROM2[k][4], bit_width);
        ROM2_D5 = DecToBin.ZZ_DecToBin(ROM2[k][5], bit_width);
        ROM2_D6 = DecToBin.ZZ_DecToBin(ROM2[k][6], bit_width);
        ROM2_D7 = DecToBin.ZZ_DecToBin(ROM2[k][7], bit_width);
        ROM2_D8 = DecToBin.ZZ_DecToBin(ROM2[k][8], bit_width);
        ROM2_D9 = DecToBin.ZZ_DecToBin(ROM2[k][9], bit_width);
        ROM2_D10 = DecToBin.ZZ_DecToBin(ROM2[k][10], bit_width);
        ROM2_D11 = DecToBin.ZZ_DecToBin(ROM2[k][11], bit_width);
        ROM2_D12 = DecToBin.ZZ_DecToBin(ROM2[k][12], bit_width);
        ROM2_D13 = DecToBin.ZZ_DecToBin(ROM2[k][13], bit_width);
        ROM2_D14 = DecToBin.ZZ_DecToBin(ROM2[k][14], bit_width);
        ROM2_D15 = DecToBin.ZZ_DecToBin(ROM2[k][15], bit_width);
        for (int i = 0; i < bit_width; i++){
            // 64 bits
            Bin_DTFAG_DIF_ROM2_B0 << ROM2_D1[bit_width-1-i];
            // 128 bits
            Bin_DTFAG_DIF_ROM2_B1 << ROM2_D2[bit_width-1-i];
            Bin_DTFAG_DIF_ROM2_B2 << ROM2_D4[bit_width-1-i];
            Bin_DTFAG_DIF_ROM2_B3 << ROM2_D6[bit_width-1-i];
            Bin_DTFAG_DIF_ROM2_B4 << ROM2_D8[bit_width-1-i];
            Bin_DTFAG_DIF_ROM2_B5 << ROM2_D10[bit_width-1-i];
            Bin_DTFAG_DIF_ROM2_B6 << ROM2_D12[bit_width-1-i];
            Bin_DTFAG_DIF_ROM2_B7 << ROM2_D14[bit_width-1-i];
        }
        Bin_DTFAG_DIF_ROM2_B1 << "_";
        Bin_DTFAG_DIF_ROM2_B2 << "_";
        Bin_DTFAG_DIF_ROM2_B3 << "_";
        Bin_DTFAG_DIF_ROM2_B4 << "_";
        Bin_DTFAG_DIF_ROM2_B5 << "_";
        Bin_DTFAG_DIF_ROM2_B6 << "_";
        Bin_DTFAG_DIF_ROM2_B7 << "_";
        for (int i = 0; i < bit_width; i++){
            // 128 bits
            Bin_DTFAG_DIF_ROM2_B1 << ROM2_D3[bit_width-1-i];
            Bin_DTFAG_DIF_ROM2_B2 << ROM2_D5[bit_width-1-i];
            Bin_DTFAG_DIF_ROM2_B3 << ROM2_D7[bit_width-1-i];
            Bin_DTFAG_DIF_ROM2_B4 << ROM2_D9[bit_width-1-i];
            Bin_DTFAG_DIF_ROM2_B5 << ROM2_D11[bit_width-1-i];
            Bin_DTFAG_DIF_ROM2_B6 << ROM2_D13[bit_width-1-i];
            Bin_DTFAG_DIF_ROM2_B7 << ROM2_D15[bit_width-1-i];
        }
        Bin_DTFAG_DIF_ROM2_B0 << endl;
        Bin_DTFAG_DIF_ROM2_B1 << endl;
        Bin_DTFAG_DIF_ROM2_B2 << endl;
        Bin_DTFAG_DIF_ROM2_B3 << endl;
        Bin_DTFAG_DIF_ROM2_B4 << endl;
        Bin_DTFAG_DIF_ROM2_B5 << endl;
        Bin_DTFAG_DIF_ROM2_B6 << endl;
        Bin_DTFAG_DIF_ROM2_B7 << endl;
    }
    //------------------------------------------------------------
}

