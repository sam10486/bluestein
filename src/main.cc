#include <iostream>
#include "NTT.h"
#include "NTTSPMB.h"
#include "DTFAG.h"
#include "DIT_NTTSPMB.h"
#include <vector>
void test_NTTSPMB();
void my_test(int argc, char *argv[]);


using namespace std;

int main(int argc, char *argv[]){
    cout << "+---------------------------------------------------------+" << endl;
    cout << "| The following examples should be executed while reading |" << endl;
    cout << "| comments in associated files in src.                    |" << endl;
    cout << "+---------------------------------------------------------+" << endl;
    cout << "| Examples                   | Source Files               |" << endl;
    cout << "+----------------------------+----------------------------+" << endl;
    cout << "| 1. test_NTTSPMB            | test-NTTSPMB.cc            |" << endl;
    cout << "| 2. my_test                 | my_test.cpp                |" << endl;
    cout << "| 3. DTFAG_SPMB_DIT          | DTFAG_SPMB_DIT.cpp         |" << endl;
    cout << "| 4. DTFAG_verify            | DTFAG_verify.cpp           |" << endl;
    cout << "| 5. DTFAG_DIT               | DTFAG_DIT.cpp              |" << endl;
    cout << "| 6. DTFAG_DIF_MixedRadix    | DTFAG_DIF_MixedRadix.cpp   |" << endl;
    cout << "+----------------------------+----------------------------+" << endl;
    
    //---------class delcare------------
    DTFAG DTFAG_SPMB_DIT;
    DTFAG DTFAG_DIT;
    DTFAG DTFAG_DIF_MixedRadix;
    DTFAG DTFAG_verify;
    //---------class delcare fin---------
    

    //-----------test declare------------
    int stage = 0;
    int radix_r1 = 2;
    vector<ZZ > st0_Tw, st1_Tw, st2_Tw;
    st0_Tw.resize(radix_r1);
    st1_Tw.resize(radix_r1);
    st2_Tw.resize(radix_r1);
    int DTFAG_t = 0;
    int DTFAG_i = 0;
    int DTFAG_j = 0;
    //------------------------------------
    
    
    int input_parameter = 6;
    int selection = 0;
    bool valid = true;
    //-------for DTFAG--------------
    unsigned long fft_point = 65536;
    ZZ fft_prime ;
    ZZ fft_twiddle_65536 ;
    ZZ fft_twiddle ;
    long difference_length = 65536 / fft_point ;
    conv(fft_prime, "18446744069414584321"); // prime number
    conv(fft_twiddle_65536, "14603442835287214144"); // twiddle factor based setting by main.cc
    
    PowerMod(fft_twiddle,fft_twiddle_65536,difference_length,fft_prime);
    
    //------------------------------
    do
    {
        cout << endl << "> Run example (1 ~ " << input_parameter << ") or exit (0): ";
        if (!(cin >> (selection))) {
            valid = false;
        }
        else if (selection < 0 || selection > input_parameter) {
            valid = false;
        }
        else {
            valid = true;
        }
        if (!valid) {
            cout << "  [Beep~~] valid option: type 0 ~ " << input_parameter << endl;
            cin.clear();
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
        }
    } while (!valid);
    switch (selection) {
        case 1:
            test_NTTSPMB();
            break;
        case 2:
            my_test(argc, argv);
            break;
        case 3:
            DTFAG_SPMB_DIT.DTFAG_SPMB_DIT(stage, st0_Tw, st1_Tw, st2_Tw, DTFAG_t, DTFAG_i, DTFAG_j);
            break;
        case 4:
            DTFAG_verify.DTFAG_verify();
            break;
        case 5:
            DTFAG_DIT.DTFAG_DIT();
            break;
        case 6:
            DTFAG_DIF_MixedRadix.DTFAG_DIF_MixedRadix();
            break;    
        case 0:
            return 0;
    }
    return 0;
}