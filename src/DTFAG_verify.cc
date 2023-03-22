#include <iostream>
#include "NTT.h"
#include "NTTSPMB.h"
#include "DTFAG.h"

using namespace std;

void DTFAG::DTFAG_verify(){
    //-------setting-------------
    int radix_r1 = 16;
    int radix_r2 = 16;
    int Tw_display = 2;
    int FFT_version = 1; // 0 is DIT, 1 is DIF
    //---------------------------

    ifstream DTFAG_pattern_Tw0;
    ifstream DTFAG_pattern_Tw1;
    ifstream DTFAG_pattern_Tw2;

    //--------TestPattern fft point ---------------
    DTFAG_pattern_Tw0.open("./SPMB_tw/MixR_DTFAG_DIF_TestPattern_Tw0.txt");
    DTFAG_pattern_Tw1.open("./SPMB_tw/MixR_DTFAG_DIF_TestPattern_Tw1.txt");
    DTFAG_pattern_Tw2.open("./SPMB_tw/MixR_DTFAG_DIF_TestPattern_Tw2.txt");
    //-----------------------------------------

    ifstream DTFAG_golden_st0;
	ifstream DTFAG_golden_st1;
	ifstream DTFAG_golden_st2;

    //--------Golden fft point---------
    switch(radix_r1){
        case 16: 
            cout << "radix_r1 = " << radix_r1 << endl;
            cout << "16" << endl;
            switch(radix_r2){
                case 16:
                    // radix_r1=16, radix_r2=16, n=65536
                    DTFAG_golden_st0.open("./SPMB_tw/DTFAG_golden_st0_65536.txt");
                    DTFAG_golden_st1.open("./SPMB_tw/DTFAG_golden_st1_65536.txt");
                    DTFAG_golden_st2.open("./SPMB_tw/DTFAG_golden_st2_65536.txt");
                    break;
                case 8:
                    // radix_r1=16, radix_r2=8, n=32768
                    DTFAG_golden_st0.open("./SPMB_tw/DTFAG_golden_st0_32768.txt");
                    DTFAG_golden_st1.open("./SPMB_tw/DTFAG_golden_st1_32768.txt");
                    DTFAG_golden_st2.open("./SPMB_tw/DTFAG_golden_st2_32768.txt");
                    break;
                case 4:
                    // radix_r1=16, radix_r2=4, n=16384
                    DTFAG_golden_st0.open("./SPMB_tw/DTFAG_golden_st0_16384.txt");
                    DTFAG_golden_st1.open("./SPMB_tw/DTFAG_golden_st1_16384.txt");
                    DTFAG_golden_st2.open("./SPMB_tw/DTFAG_golden_st2_16384.txt");
                    break;
                case 2:
                    // radix_r1=16, radix_r2=2, n=8192
                    DTFAG_golden_st0.open("./SPMB_tw/DTFAG_golden_st0_8192.txt");
                    DTFAG_golden_st1.open("./SPMB_tw/DTFAG_golden_st1_8192.txt");
                    DTFAG_golden_st2.open("./SPMB_tw/DTFAG_golden_st2_8192.txt");
                    break;
                default:
                    cout << "DTFAG_verify radix_r2 error!" << endl;
                    cout << "radix_r2 = " << radix_r2 << endl;
                    break;
            }
            break;
        case 4: 
            cout << "radix_r1 = " << radix_r1 << endl;
            cout << "4" << endl;
            switch(radix_r2){
                case 4:
                    // radix_r1=4, radix_r2=4, n=256
                    DTFAG_golden_st0.open("./SPMB_tw/DTFAG_golden_st0_256.txt");
                    DTFAG_golden_st1.open("./SPMB_tw/DTFAG_golden_st1_256.txt");
                    DTFAG_golden_st2.open("./SPMB_tw/DTFAG_golden_st2_256.txt");
                    break;
                case 2:
                    cout << "12312313213" << endl;
                    // radix_r1=4, radix_r2=2, n=128
                    DTFAG_golden_st0.open("./SPMB_tw/DTFAG_golden_st0_128.txt");
                    DTFAG_golden_st1.open("./SPMB_tw/DTFAG_golden_st1_128.txt");
                    DTFAG_golden_st2.open("./SPMB_tw/DTFAG_golden_st2_128.txt");
                    break;
                default:
                    cout << "DTFAG_verify radix_r2 error!" << endl;
                    cout << "radix_r2 = " << radix_r2 << endl;
                    break;
            }
            break;
        case 2: 
            cout << "2" << endl;
            switch(radix_r2){
                case 2:
                    // radix_r1=2, radix_r2=2, n=16
                    DTFAG_golden_st0.open("./SPMB_tw/DTFAG_golden_st0_16.txt");
                    DTFAG_golden_st1.open("./SPMB_tw/DTFAG_golden_st1_16.txt");
                    DTFAG_golden_st2.open("./SPMB_tw/DTFAG_golden_st2_16.txt");
                    break;
                default:
                    cout << "DTFAG_verify radix_r2 error!" << endl;
                    cout << "radix_r2 = " << radix_r2 << endl;
                    break;
            }
            break;
        default: 
            cout << "DTFAG_verify radix_r1 error!" << endl;
            cout << "radix_r1 = " << radix_r1 << endl;
            break;
    }
    
    //----------------------------------------------
    
    //-------------test pattern array initial-------------------
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

    Tw1_ROM.resize(radix_r2);
    for(int t=0; t<radix_r2; t++){
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

    //------------test pattern initial fin---------------------


    //--------------golden data array initial-------------------------
    vector<vector<vector<vector<ZZ > > > > st0_golden;
    vector<vector<vector<vector<ZZ > > > > st1_golden;
    vector<vector<vector<vector<ZZ > > > > st2_golden;
    st0_golden.resize(radix_r1);
    for(int t=0; t<radix_r1; t++){
        st0_golden[t].resize(radix_r1);
        for(int i=0; i<radix_r1; i++){
            st0_golden[t][i].resize(radix_r1);
            for(int len_idx=0; len_idx<radix_r1; len_idx++){
                st0_golden[t][i][len_idx].resize(radix_r1);
            }
        }
    }

    st1_golden.resize(radix_r2);
    for(int t=0; t<radix_r2; t++){
        st1_golden[t].resize(radix_r1);
        for(int i=0; i<radix_r1; i++){
            st1_golden[t][i].resize(radix_r1);
            for(int len_idx=0; len_idx<radix_r1; len_idx++){
                st1_golden[t][i][len_idx].resize(radix_r1);
            }
        }
    }

    st2_golden.resize(radix_r1);
    for(int t=0; t<radix_r1; t++){
        st2_golden[t].resize(radix_r1);
        for(int i=0; i<radix_r1; i++){
            st2_golden[t][i].resize(radix_r1);
            for(int len_idx=0; len_idx<radix_r1; len_idx++){
                st2_golden[t][i][len_idx].resize(radix_r1);
            }
        }
    }
    //---------------golden data array fin----------------------------


    //-------------stored data into ROM-------------
    if(!DTFAG_pattern_Tw0.is_open() || !DTFAG_pattern_Tw1.is_open() || !DTFAG_pattern_Tw2.is_open()){
        cout << "failed to open file.\n" << endl;
    }else {

        for(int i=0; i<radix_r2; i++){
            for(int t=0; t<radix_r1; t++){
                for(int j=0; j<radix_r1; j++){
                    for(int idx=0; idx<radix_r1; idx++){
                        DTFAG_pattern_Tw0 >> Tw0_ROM[i][t][j][idx] ;
                        DTFAG_pattern_Tw1 >> Tw1_ROM[i][t][j][idx] ;
                        DTFAG_pattern_Tw2 >> Tw2_ROM[i][t][j][idx] ;
                    }
                }
            }
        }
    }
    //-------------stored data fin--------------------


    //-------------stored golden -----------------------
    if(!DTFAG_golden_st0.is_open() || !DTFAG_golden_st1.is_open() || !DTFAG_golden_st2.is_open()){
        cout << "failed to open file.\n" << endl;
    }else {
        for(int i=0; i<radix_r2; i++){
            for(int t=0; t<radix_r1; t++){
                for(int j=0; j<radix_r1; j++){
                    for(int idx=0; idx<radix_r1; idx++){
                        DTFAG_golden_st0 >> st0_golden[i][t][j][idx] ;
                        DTFAG_golden_st1 >> st1_golden[i][t][j][idx] ;
                        DTFAG_golden_st2 >> st2_golden[i][t][j][idx] ;
                    }
                }
            }
        }
    }
    //-------------stored golden fin--------------------

    
    //------------------verify DTFAG correct or not---------------
    int Tw0_error = 0;
    int Tw1_error = 0;
    int Tw2_error = 0;
    if(FFT_version) { // DIF version compare
        for(int i=0; i<radix_r2; i++){
            cout << "   i = " << i << endl; 
            for(int t=0; t<radix_r1; t++){
                cout << "   t = " << t << endl;
                for(int j=0; j<radix_r1; j++){
                    cout << "   j = " << j << endl;
                    for(int idx=0; idx<radix_r1; idx++){
                        if(Tw0_ROM[i][t][j][idx] != st0_golden[i][t][j][idx]) {
                            Tw0_error++;
                            if(Tw_display == 0){
                                cout << "i = " << i << ", t = " << t << ", j = " << j;
                                cout << ", Tw0_ROM[" << i << "][" << t << "][" << j << "][" << idx << "] = " << Tw0_ROM[i][t][j][idx]
                                << ", st0_golden[" << i << "][" << t << "][" << j << "][" << idx << "] = " << st0_golden[i][t][j][idx] 
                                << ", error = " << Tw0_error << endl;
                            }
                        }else {
                            if(Tw_display == 0) cout << Tw0_ROM[i][t][j][idx] << " = " << st0_golden[i][t][j][idx]  << endl;
                        }
                        if(Tw1_ROM[i][t][j][idx] != st1_golden[i][t][j][idx]) {
                            Tw1_error++;
                            if(Tw_display == 1){
                                cout << "i = " << i << ", t = " << t << ", j = " << j;
                                cout << ", Tw1_ROM[" << i << "][" << t << "][" << j << "][" << idx << "] = " << Tw1_ROM[i][t][j][idx]
                                << ", st1_golden[" << i << "][" << t << "][" << j << "][" << idx << "] = " << st1_golden[i][t][j][idx] 
                                << ", error = " << Tw1_error << endl;
                            } 
                        }else {
                            if(Tw_display == 1) cout << Tw1_ROM[i][t][j][idx] << " = " << st1_golden[i][t][j][idx]  << endl;
                        }
                        if(Tw2_ROM[i][t][j][idx] != st2_golden[i][t][j][idx]) {
                            Tw2_error++;
                            if(Tw_display == 2){
                                cout << "i = " << i << ", t = " << t << ", j = " << j;
                                cout << ", Tw2_ROM[" << i << "][" << t << "][" << j << "][" << idx << "] = " << Tw2_ROM[i][t][j][idx]
                                << ", st2_golden[" << i << "][" << t << "][" << j << "][" << idx << "] = " << st2_golden[i][t][j][idx] 
                                << ", error = " << Tw0_error << endl;
                            }
                        }else {
                            if(Tw_display == 2) cout << Tw2_ROM[i][t][j][idx] << " = " << st2_golden[i][t][j][idx]  << endl;
                        }
                    }    
                }
            }
        }
    }else { // DIT version compare
        
    }
    

    if(Tw0_error != 0){
        cout << "DTFAG_Tw0 is error!, and error = " << Tw0_error << endl;
    }else{
        cout << "DTFAG_Tw0 correct !! "<< endl;
    }
    if(Tw1_error != 0){
        cout << "DTFAG_Tw1 is error!, and error = " << Tw1_error << endl;
    }else{
        cout << "DTFAG_Tw1 correct !! "<< endl;
    }
    if(Tw2_error != 0){
        cout << "DTFAG_Tw2 is error!, and error = " << Tw2_error << endl;
    }else{
        cout << "DTFAG_Tw2 correct !! "<< endl;
    }
    //------------------verify fin--------------------------------
    
}