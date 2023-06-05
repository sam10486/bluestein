#include <iostream>
#include <NTL/ZZ.h>
#include <vector>
#include "NWC_Algo.h"
#include "BitOperate.h"
#include "assert.h"
#include <sstream>

using namespace NTL;
using namespace std;


NWC_Algo::NWC_Algo(long long Radix_r1, long long Radix_r2, long long N, ZZ Modular): NWC_util(Radix_r1, Radix_r2, N, Modular){
    //TODO
}

NWC_Algo::~NWC_Algo(){
    //TODO
}

vector<ZZ > NWC_Algo::NWC(vector<ZZ > &arr){
    long long Radix_r1, Radix_r2, N;
    ZZ Modular, Phi, InvPhi;
    ZZ W, IW;
    getValue(&Radix_r1, &Radix_r2, &N, &Modular , &Phi, &InvPhi, &W, &IW);
    BitOperate Bitrev;
    vector<ZZ> Arr_scramble;
    Arr_scramble.resize(N);
    for (int i = 0; i < N; i++){
        long long bit_num = ceil(log2(N));
        long long rev_index = Bitrev.BitReserve(i, bit_num);
        Arr_scramble[i] = arr[rev_index];
    }
    for (int s = 1; s <= log2(N); s++){
        //cout << "-------------stage = " << s << "---------------" << endl;
        int m = pow(2, s);
        for (int j = 0; j <= (m/2)-1; j++){
            ZZ TF;
            int TF_deg = (2*j+1) * (N/m);
            TF = PowerMod(Phi, TF_deg, Modular);
            for (int k = 0; k <= (N/m)-1; k++){
                ZZ u, t;
                u = Arr_scramble[k*m+j];
                t = MulMod(TF, Arr_scramble[k*m+j+(m/2)], Modular);
                //cout << "u[" << k*m+j << "] = " << u << endl;
                //cout << "t[" << k*m+j+(m/2) << "] = " << t << ", TF = " << TF << endl;
                Arr_scramble[k*m+j] = AddMod(u, t, Modular);
                Arr_scramble[k*m+j+(m/2)] = SubMod(u, t, Modular);
                //cout << "Arr_scramble[" << k*m+j << "] = " << Arr_scramble[k*m+j] << endl;
                //cout << "Arr_scramble[" << k*m+j+(m/2) << "] = " << Arr_scramble[k*m+j+(m/2)] << endl;
            }
        }
    }
    for (int i = 0; i < N; i++){
        arr[i] = Arr_scramble[i];
    }
    
    return arr;
}

vector<ZZ > NWC_Algo::INWC(vector<ZZ> &arr){
    long long Radix_r1, Radix_r2, N;
    ZZ Modular, Phi, InvPhi;
    ZZ W, IW;
    getValue(&Radix_r1, &Radix_r2, &N, &Modular , &Phi, &InvPhi, &W, &IW);
    
    vector<ZZ > A_arr;
    A_arr.resize(N);
    for (int i = 0; i < N; i++){
        A_arr[i] = arr[i];
    }
    long long bit_num = ceil(log2(N));

    BitOperate Bitrev;
    ZZ InvTwo;
    InvMod(InvTwo, (ZZ)2, Modular);
    //cout << "InvTwo = " << InvTwo << endl;
    for (int s = log2(N); s >= 1; s--){
        //cout << "stage = " << s << endl;
        int m = pow(2, s);
        for (int j = 0; j <= (m/2)-1; j++){
            //ZZ TF;
            //int TF_deg = (2*j+1) * (N/m);
            //TF = PowerMod(InvPhi, TF_deg, Modular);
            int Phi_deg, W_deg;
            Phi_deg = (N/m);
            W_deg = (j) * (N/m);
            ZZ InvPhi_Order = PowerMod(InvPhi, Phi_deg, Modular);
            ZZ IW_Order = PowerMod(IW, W_deg, Modular);
            //cout << "Phi_deg = " << Phi_deg << ", W_deg = " << W_deg << endl;
            //cout << "TF = " << TF << ", InvPhi_Order = " << InvPhi_Order << ", IW = " << IW_Order << endl;
            //cout << "---------------" << endl;
            for (int k = 0; k <= (N/m)-1 ; k++){
                ZZ u, t;
                u = A_arr[k*m+j];
                t = A_arr[k*m+j+(m/2)];
                //cout << "u[" << k*m+j << "] = " << u << endl;
                //cout << "t[" << k*m+j+(m/2) << "] = " << t << endl;
                // upper
                A_arr[k*m+j] = AddMod(u, t, Modular);
                A_arr[k*m+j] = MulMod(A_arr[k*m+j], InvTwo, Modular);
                // down
                ZZ InvPhi_dot_IW = MulMod(InvPhi_Order, IW_Order, Modular);
                ZZ InvPhi_dot_IW_dot_Inv_two = MulMod(InvPhi_dot_IW, InvTwo, Modular);
                A_arr[k*m+j+(m/2)] = SubMod(u, t, Modular);
                A_arr[k*m+j+(m/2)] = MulMod(A_arr[k*m+j+(m/2)], InvPhi_dot_IW_dot_Inv_two, Modular);
                //cout << "A_arr[" << k*m+j << "] = " << A_arr[k*m+j] << endl;
                //cout << "A_arr[" << k*m+j+(m/2) << "] = " << A_arr[k*m+j+(m/2)] << endl;
            }
        }
    }
    for (int i = 0; i < N; i++){
        long long rev_idx = Bitrev.BitReserve(i, bit_num);
        arr[i] = A_arr[rev_idx];
    }
    return arr;
}

void NWC_Algo::time_o_r16(std::vector<ZZ> time_data,std::string string_in){
    long long Radix_r1, Radix_r2, N;
    ZZ Modular, Phi, InvPhi;
    ZZ W, IW;
    getValue(&Radix_r1, &Radix_r2, &N, &Modular , &Phi, &InvPhi, &W, &IW);


    std::vector<int> ma;
    std::vector<int> bn;
    std::vector<int> bit_array_tmp;
    unsigned long fft_point;
    unsigned long radix;

    fft_point = N;
    radix = Radix_r1;

    int offset             =  fft_point/radix;
    int group              =  fft_point/(radix*radix);
    int fft_point_bit_size =  log2(fft_point); 
    unsigned long counter_iteration = fft_point/(radix*2);
    ma.resize(fft_point);
    bn.resize(fft_point);
    bit_array_tmp.resize(fft_point_bit_size);

    int ma_tmp;
    int bit_tmp;
    int bn_tmp;
    int tmp = 0;
    int addr_tmp = 0;
    
    for(int i=0; i < group; i++){        
        for(int ss=0 ; ss < radix; ss++){
            bn_tmp = 0;
            ma_tmp = 0;
            tmp = ss * group + i;    
            //bit calculate
            for(int j=0; j < fft_point_bit_size;j++){
                bit_tmp = tmp % 2;
                tmp = tmp >> 1;
                bit_array_tmp[j] = bit_tmp;
            } 
            for(int rs = 0; rs < fft_point_bit_size; rs++){
                if((bit_array_tmp[rs] == 1) && (rs != 0)) ma_tmp = ma_tmp + exp2((rs-1)); 
                bn_tmp = bn_tmp ^ bit_array_tmp[rs];
            }
            for(int kk=0; kk<radix; kk++){
                addr_tmp = ss*group + i + kk * (offset);
                if(addr_tmp >= fft_point)std::cout << addr_tmp <<"\n";
                ma[addr_tmp] = ma_tmp;
                bn[addr_tmp] = bn_tmp;
            }
        }    
    }


    ofstream b0radix0(string_in   +  "b0radix0.txt");
    ofstream b0radix1(string_in   +  "b0radix1.txt");
    ofstream b0radix2(string_in   +  "b0radix2.txt");
    ofstream b0radix3(string_in   +  "b0radix3.txt");
    ofstream b0radix4(string_in   +  "b0radix4.txt");
    ofstream b0radix5(string_in   +  "b0radix5.txt");
    ofstream b0radix6(string_in   +  "b0radix6.txt");
    ofstream b0radix7(string_in   +  "b0radix7.txt");
    ofstream b0radix8(string_in   +  "b0radix8.txt");
    ofstream b0radix9(string_in   +  "b0radix9.txt");
    ofstream b0radix10(string_in  +  "b0radix10.txt");
    ofstream b0radix11(string_in  +  "b0radix11.txt");
    ofstream b0radix12(string_in  +  "b0radix12.txt");
    ofstream b0radix13(string_in  +  "b0radix13.txt");
    ofstream b0radix14(string_in  +  "b0radix14.txt");
    ofstream b0radix15(string_in  +  "b0radix15.txt");
    ofstream b1radix0(string_in   +  "b1radix0.txt");
    ofstream b1radix1(string_in   +  "b1radix1.txt");
    ofstream b1radix2(string_in   +  "b1radix2.txt");
    ofstream b1radix3(string_in   +  "b1radix3.txt");
    ofstream b1radix4(string_in   +  "b1radix4.txt");
    ofstream b1radix5(string_in   +  "b1radix5.txt");
    ofstream b1radix6(string_in   +  "b1radix6.txt");
    ofstream b1radix7(string_in   +  "b1radix7.txt");
    ofstream b1radix8(string_in   +  "b1radix8.txt");
    ofstream b1radix9(string_in   +  "b1radix9.txt");
    ofstream b1radix10(string_in  +  "b1radix10.txt");
    ofstream b1radix11(string_in  +  "b1radix11.txt");
    ofstream b1radix12(string_in  +  "b1radix12.txt");
    ofstream b1radix13(string_in  +  "b1radix13.txt");
    ofstream b1radix14(string_in  +  "b1radix14.txt");
    ofstream b1radix15(string_in  +  "b1radix15.txt");
   
    std::string string_buf; 
    
    ZZ time_data_tmp; 

    int counter =0;
    int loop_index = 0;
    int loop_address = 0;
 
    //bank0
    while(counter < counter_iteration){ 
       if((ma[loop_index] == counter )&& (bn[loop_index] == 0) ){
          for(int j=0;j<radix;j++){
            loop_address   =  loop_index + offset*j;
            time_data_tmp  =  time_data[loop_address];
            string_buf     =  ZZtohex_cp(time_data_tmp);
            
            if(j==0 ) { b0radix0  << string_buf; b0radix0  << "\n"; }     
            if(j==1 ) { b0radix1  << string_buf; b0radix1  << "\n"; }
            if(j==2 ) { b0radix2  << string_buf; b0radix2  << "\n"; } 
            if(j==3 ) { b0radix3  << string_buf; b0radix3  << "\n"; }
            if(j==4 ) { b0radix4  << string_buf; b0radix4  << "\n"; }
            if(j==5 ) { b0radix5  << string_buf; b0radix5  << "\n"; }
            if(j==6 ) { b0radix6  << string_buf; b0radix6  << "\n"; }
            if(j==7 ) { b0radix7  << string_buf; b0radix7  << "\n"; }
            if(j==8 ) { b0radix8  << string_buf; b0radix8  << "\n"; }
            if(j==9 ) { b0radix9  << string_buf; b0radix9  << "\n"; }
            if(j==10) { b0radix10 << string_buf; b0radix10 << "\n"; }
            if(j==11) { b0radix11 << string_buf; b0radix11 << "\n"; }
            if(j==12) { b0radix12 << string_buf; b0radix12 << "\n"; }
            if(j==13) { b0radix13 << string_buf; b0radix13 << "\n"; }
            if(j==14) { b0radix14 << string_buf; b0radix14 << "\n"; }
            if(j==15) { b0radix15 << string_buf; b0radix15 << "\n"; }

          }           
           loop_index = 0;
           counter = counter + 1;
       }
       else loop_index = loop_index + 1;
    }
  
    counter      = 0;
    loop_index   = 0;
    loop_address = 0;
    //bank 1
    while(counter < counter_iteration){ 
        if((ma[loop_index] == counter )&& (bn[loop_index] == 1) ){
            for(int j=0;j<radix;j++){
            loop_address   = loop_index + offset*j;
            time_data_tmp  = time_data[loop_address];
            string_buf     = ZZtohex_cp(time_data_tmp);

            if(j==0 ) { b1radix0  << string_buf; b1radix0  << "\n"; }
            if(j==1 ) { b1radix1  << string_buf; b1radix1  << "\n"; }
            if(j==2 ) { b1radix2  << string_buf; b1radix2  << "\n"; }
            if(j==3 ) { b1radix3  << string_buf; b1radix3  << "\n"; }
            if(j==4 ) { b1radix4  << string_buf; b1radix4  << "\n"; }
            if(j==5 ) { b1radix5  << string_buf; b1radix5  << "\n"; }
            if(j==6 ) { b1radix6  << string_buf; b1radix6  << "\n"; }
            if(j==7 ) { b1radix7  << string_buf; b1radix7  << "\n"; }
            if(j==8 ) { b1radix8  << string_buf; b1radix8  << "\n"; }
            if(j==9 ) { b1radix9  << string_buf; b1radix9  << "\n"; }
            if(j==10) { b1radix10 << string_buf; b1radix10 << "\n"; }
            if(j==11) { b1radix11 << string_buf; b1radix11 << "\n"; }
            if(j==12) { b1radix12 << string_buf; b1radix12 << "\n"; }
            if(j==13) { b1radix13 << string_buf; b1radix13 << "\n"; }
            if(j==14) { b1radix14 << string_buf; b1radix14 << "\n"; }
            if(j==15) { b1radix15 << string_buf; b1radix15 << "\n"; }
            
            }           
             loop_index = 0;
             counter    = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
    
    //time_of.close();
    b0radix0.close();
	b0radix1.close();
	b0radix2.close();
	b0radix3.close();
	b0radix4.close();
	b0radix5.close();
	b0radix6.close();
	b0radix7.close();
	b0radix8.close();
	b0radix9.close();
	b0radix10.close();
	b0radix11.close();
	b0radix12.close();
	b0radix13.close();
	b0radix14.close();
	b0radix15.close();
	b1radix0.close();
	b1radix1.close();
	b1radix2.close();
	b1radix3.close();
	b1radix4.close();
	b1radix5.close();
	b1radix6.close();
	b1radix7.close();
	b1radix8.close();
	b1radix9.close();
	b1radix10.close();
	b1radix11.close();
	b1radix12.close();
	b1radix13.close();
	b1radix14.close();
	b1radix15.close();
	
}

std::string NWC_Algo::ZZtohex_cp(ZZ zz_tmp){
    std::string string_tmp;
    std::vector<char> tmp_hex; 
    std::stringstream ss;
    
    long tmp;
    int length;
	length = 6;
    tmp_hex.resize(length);
    
    for(int i =0; i < length; i++){
        tmp  =  to_long(zz_tmp % 16);
        
        if(tmp == 0) tmp_hex[i] = '0'; 
        if(tmp == 1) tmp_hex[i] = '1'; 
        if(tmp == 2) tmp_hex[i] = '2'; 
        if(tmp == 3) tmp_hex[i] = '3'; 
        if(tmp == 4) tmp_hex[i] = '4'; 
        if(tmp == 5) tmp_hex[i] = '5'; 
        if(tmp == 6) tmp_hex[i] = '6'; 
        if(tmp == 7) tmp_hex[i] = '7'; 
        if(tmp == 8) tmp_hex[i] = '8'; 
        if(tmp == 9) tmp_hex[i] = '9'; 
        if(tmp == 10) tmp_hex[i] = 'a'; 
        if(tmp == 11) tmp_hex[i] = 'b'; 
        if(tmp == 12) tmp_hex[i] = 'c'; 
        if(tmp == 13) tmp_hex[i] = 'd'; 
        if(tmp == 14) tmp_hex[i] = 'e'; 
        if(tmp == 15) tmp_hex[i] = 'f'; 
        
        zz_tmp = zz_tmp >> 4;
    }
    for(int i = 1;i <= length; i++){
        ss << tmp_hex[length-i];
    }
    
    string_tmp = ss.str();
    
    return string_tmp;
}

