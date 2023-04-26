#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <random>
#include <bitset>
#include "BitOperate.h"
#include "math.h"
#include <NTT.h>

using namespace std;



long long BitOperate::BitReserve(long long DataToReverse, long long BitLength){
    long long result = 0;
    for(long long i = 0; i < BitLength; i++){
        if((DataToReverse >> i) & 1){
            result |= 1 << (BitLength - 1 -i);
        }
    }
    return result;
}

vector<long long> BitOperate::DecToBin(long long data, long long bit_width){
    vector<long long> BinVec(bit_width);
    for(long long int j=0; j<bit_width; j++){
        BinVec.at(j) = (data >> j) & 1;
    }
    return BinVec;
}

void BitOperate::IntToVec(long long integer, long long N, vector<long long> &bit_array){
    long long bit_width = (long long)ceil(log2(N));
    bit_array.resize(bit_width);
    for(long long j=0; j < bit_width; j++){
        bit_array[j] = (integer >> j) & 1;
        //cout << bit_array[j] << endl;
    }
    
}

long long BitOperate::VecToInt(vector<long long> bit_array, long long N){
    long long bit_width = (long long)ceil(log2(N));
    long long integer = 0;
    for(long long j=0; j < bit_width; j++){
        integer += bit_array[j] << j;
        //cout << "bit_array[" << j << "] = " << bit_array[j] << endl;
    }
    return integer;
}

long long BitOperate::RR(long long BC, long long shift_bit, long long N, long long r){
    long long RR_out = 0;
    long long bit_width = (long long)ceil(log2(N/r));
    vector<long long> bit_array(bit_width);
    BitOperate DecToBin;
    bit_array = DecToBin.DecToBin(BC, bit_width);

    rotate(bit_array.begin(), bit_array.begin()+shift_bit, bit_array.end());
    for(long long j=0; j < bit_width; j++){
        RR_out += bit_array[j] << j;
        //cout << RR_out << endl ;
    }
    return RR_out;
}

long long BitOperate::Gray_code(long long index, long long group){
    long long group_bit_width = (long long)ceil(log2(group));
    vector<long long> index_bit_array(group_bit_width);
    BitOperate DecToBin, VecToInt;
    index_bit_array = DecToBin.DecToBin(index, group_bit_width);
    
    for(long long i=0; i<group_bit_width; i++){
        if(i == (group_bit_width-1))
            index_bit_array.at(i) = index_bit_array.at(i);
        else
            index_bit_array.at(i) = index_bit_array.at(i) ^ index_bit_array.at(i+1);      
    } 
    long long Gray_num = VecToInt.VecToInt(index_bit_array, group);
    return Gray_num;
}

long long BitOperate::unary_xor(long long data_in, long long bit_width){
    BitOperate DecToBin;
    vector<long long> bit_array = DecToBin.DecToBin(data_in, bit_width);
    long long xor_out = 0;
    for(long long i=0; i<bit_width; i++){
        xor_out += bit_array[i];
    }
    xor_out %= 2;
    return xor_out;
}

long long BitOperate::left_rotate(long long input, long long shift_bit, long long N){
    long long RR_out = 0;
    long long bit_width = (long long)ceil(log2(N));
    //cout << "input = " << input << ", " << "shift_bit = " << shift_bit << ", " << "N = " << N << endl;
    vector<long long> bit_array(bit_width);
    BitOperate DecToBin;
    bit_array = DecToBin.DecToBin(input, bit_width);

    rotate(bit_array.begin(), bit_array.begin()+shift_bit, bit_array.end());
    //cout << "----------" << endl;
    for(long long j=0; j < bit_width; j++){
        RR_out += bit_array[j] << j;
        //cout << RR_out << endl ;
    }
    return RR_out;
}


//------------------for memory initial------------------------
long long BitOperate::DecToBin_mem_init(long long *BinVec, long long data, long long bit_width){
    //cout << endl;
    //cout << data << " ";
    for(long long int j=0; j<bit_width; j++){
        BinVec[j] = (data >> j) & 1;
        //cout << BinVec[j];
    }
    return 0;
}

long long BitOperate::unary_xor_mem_init(long long data_in, long long bit_width_m, long long s){
    BitOperate DecToBin;
    vector<long long> bit_array = DecToBin.DecToBin(data_in, bit_width_m);
    long long xor_out = 0;
    for(long long i=0; i<bit_width_m-s; i++){
        xor_out += bit_array[i];
    }
    xor_out %= 2;
    return xor_out;
}

long long BitOperate::VecToInt_mem_init(long long data_in, long long bit_width_m, long long s){
    BitOperate DecToBin;
    vector<long long> bit_array = DecToBin.DecToBin(data_in, bit_width_m);
    long long integer = 0;
    for(long long j=1; j < bit_width_m-s; j++){
        integer += bit_array[j] << j-1;
        //cout << "bit_array[" << j << "] = " << bit_array[j] << endl;
    }
    return integer;
}

long long BitOperate::number_complement(long long i, long long radix_r1){

    BitOperate DecToBin, VecToInt;
    
    vector<long long > i_complement_array;
    long long i_complement;
    long long bit_width = (long long)ceil(log2(radix_r1));
    i_complement_array.resize(bit_width);

    //cout << "i = " << i << endl;

    vector<long long > i_array = DecToBin.DecToBin(i, bit_width);
    long long tmp;
    for(long long i=0; i<bit_width; i++){
        tmp = i_array[i];
        if(tmp){
            i_complement_array[i] = 0;
        }else{
            i_complement_array[i] = 1;
        }
    }
    i_complement = VecToInt.VecToInt(i_complement_array, radix_r1);

    //cout << "i_complement = " << i_complement << endl;
    return i_complement;
}

vector<ZZ > BitOperate::ZZ_DecToBin(ZZ data, long long bit_width){
    vector<ZZ > BinVec(bit_width);
    for(long long int j=0; j<bit_width; j++){
        BinVec.at(j) = (data >> j) & 1;
    }
    return BinVec;
}