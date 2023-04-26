#ifndef _BitOperate_
#define _BitOperate_

#include <iostream>
#include <vector>
#include <NTT.h>
using namespace std; 

class BitOperate{
public:
    long long BitReserve(long long DataToReverse, long long BitLength);
    vector<long long> DecToBin(long long data, long long bit_width);
    void IntToVec(long long integer, long long N, vector<long long> &bit_array);
    long long VecToInt(vector<long long> bit_array, long long N);
    long long RR(long long BC, long long shift_bit, long long N, long long r);
    long long Gray_code(long long index, long long group);
    long long unary_xor(long long data_in, long long bit_width);
    long long DecToBin_mem_init(long long *BinVec, long long data, long long bit_width);
    long long unary_xor_mem_init(long long data_in, long long bit_width_m, long long s);
    long long VecToInt_mem_init(long long data_in, long long bit_width_m, long long s);
    long long left_rotate(long long input, long long shift_bit, long long N);
    long long number_complement(long long i, long long radix_r1);
    vector<ZZ > ZZ_DecToBin(ZZ data, long long bit_width);
};

#endif