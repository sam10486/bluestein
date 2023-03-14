#include <iostream>
#include "NTT.h"
#include "NTTSPMB.h"
#include "DTFAG.h"

void test_NTTSPMB();
void my_test(int argc, char *argv[]);
void DTFAG_verify();
void DTFAG_test();
void DTFAG_DIF();
void DTFAG_MixedRadix();

using namespace std;

int main(int argc, char *argv[]){
    //while (true)
    //{
        cout << "+---------------------------------------------------------+" << endl;
        cout << "| The following examples should be executed while reading |" << endl;
        cout << "| comments in associated files in src.                    |" << endl;
        cout << "+---------------------------------------------------------+" << endl;
        cout << "| Examples                   | Source Files               |" << endl;
        cout << "+----------------------------+----------------------------+" << endl;
        cout << "| 1. test_NTTSPMB            | test-NTTSPMB.cc            |" << endl;
        cout << "| 2. my_test                 | my_test.cpp                |" << endl;
        cout << "| 3. DTFAG                   | DTFAG.cpp                  |" << endl;
        cout << "| 4. DTFAG_verify            | DTFAG_verify.cpp           |" << endl;
        cout << "| 5. DTFAG_test              | DTFAG_test.cpp             |" << endl;
        cout << "| 6. DTFAG_DIF               | DTFAG_DIF.cpp              |" << endl;
        cout << "| 7. DTFAG_MixedRadix        | DTFAG_MixedRadix.cpp       |" << endl;
        cout << "+----------------------------+----------------------------+" << endl;

        int input_parameter = 7;
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
                DTFAG(fft_prime, fft_twiddle);
                break;
            case 4:
                DTFAG_verify();
                break;
            case 5:
                DTFAG_test();
                break;
            case 6:
                DTFAG_DIF();
                break;
            case 7:
                DTFAG_MixedRadix();
                break;    
            case 0:
                return 0;
        }
    //}
    return 0;
}