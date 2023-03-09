#include <iostream>
#include "NTT.h"
#include "NTTSPMB.h"

using namespace std;

void DTFAG_verify(){
    //-------setting-------------
    int radix_r = 16;
    int Tw_display = 2;
    //---------------------------

    ifstream DTFAG_pattern_Tw0;
    ifstream DTFAG_pattern_Tw1;
    ifstream DTFAG_pattern_Tw2;

    //--------TestPattern fft point ---------------
    DTFAG_pattern_Tw0.open("./SPMB_tw/DTFAG_TestPattern_Tw0.txt");
    DTFAG_pattern_Tw1.open("./SPMB_tw/DTFAG_TestPattern_Tw1.txt");
    DTFAG_pattern_Tw2.open("./SPMB_tw/DTFAG_TestPattern_Tw2.txt");
    //-----------------------------------------

    ifstream DTFAG_golden_st0;
	ifstream DTFAG_golden_st1;
	ifstream DTFAG_golden_st2;

    //--------Golden fft point---------
    if(radix_r==16) DTFAG_golden_st0.open("./SPMB_tw/DTFAG_golden_st0_65536.txt");
    if(radix_r==16) DTFAG_golden_st1.open("./SPMB_tw/DTFAG_golden_st1_65536.txt");
    if(radix_r==16) DTFAG_golden_st2.open("./SPMB_tw/DTFAG_golden_st2_65536.txt");
    //----------------------------------------------

    //--------Golden fft point-----------
    if(radix_r==4) DTFAG_golden_st0.open("./SPMB_tw/DTFAG_golden_st0_256.txt");
    if(radix_r==4) DTFAG_golden_st1.open("./SPMB_tw/DTFAG_golden_st1_256.txt");
    if(radix_r==4) DTFAG_golden_st2.open("./SPMB_tw/DTFAG_golden_st2_256.txt");
    //----------------------------------------------

    //--------Golden fft point-----------
    if(radix_r==2) DTFAG_golden_st0.open("./SPMB_tw/DTFAG_golden_st0_16.txt");
    if(radix_r==2) DTFAG_golden_st1.open("./SPMB_tw/DTFAG_golden_st1_16.txt");
    if(radix_r==2) DTFAG_golden_st2.open("./SPMB_tw/DTFAG_golden_st2_16.txt");
    //----------------------------------------------
    

    //-------------test pattern array initial-------------------
    vector<vector<vector<vector<ZZ > > > > Tw0_ROM;
    vector<vector<vector<vector<ZZ > > > > Tw1_ROM;
    vector<vector<vector<vector<ZZ > > > > Tw2_ROM;
    Tw0_ROM.resize(radix_r);
    for(int t=0; t<radix_r; t++){
        Tw0_ROM[t].resize(radix_r);
        for(int i=0; i<radix_r; i++){
            Tw0_ROM[t][i].resize(radix_r);
            for(int len_idx=0; len_idx<radix_r; len_idx++){
                Tw0_ROM[t][i][len_idx].resize(radix_r);
            }
        }
    }

    Tw1_ROM.resize(radix_r);
    for(int t=0; t<radix_r; t++){
        Tw1_ROM[t].resize(radix_r);
        for(int i=0; i<radix_r; i++){
            Tw1_ROM[t][i].resize(radix_r);
            for(int len_idx=0; len_idx<radix_r; len_idx++){
                Tw1_ROM[t][i][len_idx].resize(radix_r);
            }
        }
    }

    Tw2_ROM.resize(radix_r);
    for(int t=0; t<radix_r; t++){
        Tw2_ROM[t].resize(radix_r);
        for(int i=0; i<radix_r; i++){
            Tw2_ROM[t][i].resize(radix_r);
            for(int len_idx=0; len_idx<radix_r; len_idx++){
                Tw2_ROM[t][i][len_idx].resize(radix_r);
            }
        }
    }

    //------------test pattern initial fin---------------------


    //--------------golden data array initial-------------------------
    vector<vector<vector<vector<ZZ > > > > st0_golden;
    vector<vector<vector<vector<ZZ > > > > st1_golden;
    vector<vector<vector<vector<ZZ > > > > st2_golden;
    st0_golden.resize(radix_r);
    for(int t=0; t<radix_r; t++){
        st0_golden[t].resize(radix_r);
        for(int i=0; i<radix_r; i++){
            st0_golden[t][i].resize(radix_r);
            for(int len_idx=0; len_idx<radix_r; len_idx++){
                st0_golden[t][i][len_idx].resize(radix_r);
            }
        }
    }

    st1_golden.resize(radix_r);
    for(int t=0; t<radix_r; t++){
        st1_golden[t].resize(radix_r);
        for(int i=0; i<radix_r; i++){
            st1_golden[t][i].resize(radix_r);
            for(int len_idx=0; len_idx<radix_r; len_idx++){
                st1_golden[t][i][len_idx].resize(radix_r);
            }
        }
    }

    st2_golden.resize(radix_r);
    for(int t=0; t<radix_r; t++){
        st2_golden[t].resize(radix_r);
        for(int i=0; i<radix_r; i++){
            st2_golden[t][i].resize(radix_r);
            for(int len_idx=0; len_idx<radix_r; len_idx++){
                st2_golden[t][i][len_idx].resize(radix_r);
            }
        }
    }
    //---------------golden data array fin----------------------------


    //-------------stored data into ROM-------------
    if(!DTFAG_pattern_Tw0.is_open() || !DTFAG_pattern_Tw1.is_open() || !DTFAG_pattern_Tw2.is_open()){
        cout << "failed to open file.\n" << endl;
    }else {

        for(int i=0; i<radix_r; i++){
            for(int t=0; t<radix_r; t++){
                for(int j=0; j<radix_r; j++){
                    for(int idx=0; idx<radix_r; idx++){
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
        for(int i=0; i<radix_r; i++){
            for(int t=0; t<radix_r; t++){
                for(int j=0; j<radix_r; j++){
                    for(int idx=0; idx<radix_r; idx++){
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
    for(int i=0; i<radix_r; i++){
        cout << "   i = " << i << endl; 
        for(int t=0; t<radix_r; t++){
            cout << "   t = " << t << endl;
            for(int j=0; j<radix_r; j++){
                cout << "   j = " << j << endl;
                for(int idx=0; idx<radix_r; idx++){
                    if(Tw0_ROM[i][t][j][idx] != st0_golden[i][t][j][idx]) {
                        Tw0_error++;
                        if(Tw_display == 0){
                            cout << "i = " << i << ", t = " << t << ", j = " << j;
                            cout << ", Tw0_ROM[" << i << "][" << t << "][" << j << "][" << idx << "] = " << Tw0_ROM[i][t][j][idx]
                            << ", st2_golden[" << i << "][" << t << "][" << j << "][" << idx << "] = " << st0_golden[i][t][j][idx] 
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
                            << ", st0_golden[" << i << "][" << t << "][" << j << "][" << idx << "] = " << st2_golden[i][t][j][idx] 
                            << ", error = " << Tw0_error << endl;
                        }
                    }else {
                        if(Tw_display == 2) cout << Tw2_ROM[i][t][j][idx] << " = " << st2_golden[i][t][j][idx]  << endl;
                    }

                }
            }
        }
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