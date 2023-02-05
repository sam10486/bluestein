#include <iostream>
#include "NTT.h"
#include "NTTSPMB.h"

void test_NTTSPMB();
void my_test(int argc, char *argv[]);
void DTFAG();


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
            DTFAG();
            break;
        case 0:
            return 0;
        }
    //}
    return 0;
}