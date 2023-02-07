#include <iostream>
#include "NTT.h"
#include "NTTSPMB.h"
#include "DTFAG.h"

void test_NTTSPMB();
void my_test(int argc, char *argv[]);

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
        cout << "+----------------------------+----------------------------+" << endl;

        int input_parameter = 3;
        int selection = 0;
        bool valid = true;
        //-------for DTFAG--------------
        ZZ P ;
        ZZ W ;
        conv(P, "18446744069414584321"); // prime number
        conv(W, "14603442835287214144"); // twiddle factor based setting by main.cc
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
            DTFAG(P, W);
            break;
        case 0:
            return 0;
        }
    //}
    return 0;
}