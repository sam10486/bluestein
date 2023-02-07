#include <iostream>
#include "NTT.h"
#include "NTTSPMB.h"

using namespace std;

void DTFAG_verify(){

    ifstream DTFAG_pattern_Tw0;
    ifstream DTFAG_pattern_Tw1;
    ifstream DTFAG_pattern_Tw2;

    DTFAG_pattern_Tw0.open("./my_print_data/DTFAG_pattern_Tw0.txt");
    DTFAG_pattern_Tw1.open("./my_print_data/DTFAG_pattern_Tw1.txt");
    DTFAG_pattern_Tw2.open("./my_print_data/DTFAG_pattern_Tw2.txt");

    ifstream DTFAG_golden_st0;
	ifstream DTFAG_golden_st1;
	ifstream DTFAG_golden_st2;

    DTFAG_golden_st0.open("./my_print_data/DTFAG_golden_st0.txt");
    DTFAG_golden_st1.open("./my_print_data/DTFAG_golden_st1.txt");
    DTFAG_golden_st2.open("./my_print_data/DTFAG_golden_st2.txt");

    int radix_r = 16;
    int Tw_display = 2;

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

        for(int t=0; t<radix_r; t++){
            for(int i=0; i<radix_r; i++){
                for(int len_idx=0; len_idx<radix_r; len_idx++){
                    for(int j=0; j<radix_r; j++){
                        DTFAG_pattern_Tw0 >> Tw0_ROM[t][i][len_idx][j] ;
                        DTFAG_pattern_Tw1 >> Tw1_ROM[t][i][len_idx][j] ;
                        DTFAG_pattern_Tw2 >> Tw2_ROM[t][i][len_idx][j] ;
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
        for(int t=0; t<radix_r; t++){
            for(int i=0; i<radix_r; i++){
                for(int len_idx=0; len_idx<radix_r; len_idx++){
                    for(int j=0; j<radix_r; j++){
                        DTFAG_golden_st0 >> st0_golden[t][i][len_idx][j] ;
                        DTFAG_golden_st1 >> st1_golden[t][i][len_idx][j] ;
                        DTFAG_golden_st2 >> st2_golden[t][i][len_idx][j] ;
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
    for(int t=0; t<radix_r; t++){
        cout << "   t = " << t << endl; 
        for(int i=0; i<radix_r; i++){
            cout << "   i = " << i << endl;
            for(int len_idx=0; len_idx<radix_r; len_idx++){
                cout << "   len_idx = " << len_idx << endl;
                for(int j=0; j<radix_r; j++){
                    if(Tw0_ROM[t][i][len_idx][j] != st2_golden[t][i][len_idx][j]) {
                        Tw0_error++;
                        if(Tw_display == 0){
                            cout << "t = " << t << ", i = " << i << ", len_idx = " << len_idx;
                            cout << ", Tw0_ROM[" << i << "][" << i << "][" << len_idx << "][" << j << "] = " << Tw0_ROM[t][i][len_idx][j]
                            << ", st2_golden[" << i << "][" << i << "][" << len_idx << "][" << j << "] = " << st2_golden[t][i][len_idx][j] 
                            << ", error = " << Tw0_error << endl;
                        }
                    }else {
                        if(Tw_display == 0) cout << "TF" << j << ", " << Tw0_ROM[t][i][len_idx][j] << " = " << st2_golden[t][i][len_idx][j]  << endl;
                    }
                    if(Tw1_ROM[t][i][len_idx][j] != st1_golden[t][i][len_idx][j]) {
                        Tw1_error++;
                        if(Tw_display == 1){
                            cout << "t = " << t << ", i = " << i << ", len_idx = " << len_idx;
                            cout << ", Tw1_ROM[" << i << "][" << i << "][" << len_idx << "][" << j << "] = " << Tw1_ROM[t][i][len_idx][j]
                            << ", st1_golden[" << i << "][" << i << "][" << len_idx << "][" << j << "] = " << st1_golden[t][i][len_idx][j] 
                            << ", error = " << Tw1_error << endl;
                        } 
                        
                    }else {
                        if(Tw_display == 1) cout << "TF" << j << ", " << Tw1_ROM[t][i][len_idx][j] << " = " << st1_golden[t][i][len_idx][j]  << endl;
                    }
                    if(Tw2_ROM[t][i][len_idx][j] != st0_golden[t][i][len_idx][j]) {
                        Tw2_error++;
                        if(Tw_display == 2){
                            cout << "t = " << t << ", i = " << i << ", len_idx = " << len_idx;
                            cout << ", Tw2_ROM[" << i << "][" << i << "][" << len_idx << "][" << j << "] = " << Tw2_ROM[t][i][len_idx][j]
                            << ", st0_golden[" << i << "][" << i << "][" << len_idx << "][" << j << "] = " << st0_golden[t][i][len_idx][j] 
                            << ", error = " << Tw0_error << endl;
                        }
                    }else {
                        if(Tw_display == 2) cout << "TF" << j << ", " << Tw2_ROM[t][i][len_idx][j] << " = " << st0_golden[t][i][len_idx][j]  << endl;
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