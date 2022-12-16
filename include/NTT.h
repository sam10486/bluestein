#ifndef _NTT_H_
#define _NTT_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <NTL/ZZ.h>

using namespace NTL;

class NTT{
public:
	unsigned long  N; // NTT ponit
	ZZ N_ZZ; 
	ZZ IN; // inverse of N (mod p)

	ZZ  p; // modulus of NTT
	       // must be a prime
	ZZ  W; // primitive nth root of unity in Zp
	       // W^(N) (mod p) = 1
	ZZ  IW;// inverse of W  (mod p)
	
	void NTT_init(unsigned long n, ZZ prime, ZZ root); //init parameters  
	void NTT_BU(ZZ &a, ZZ &b); //buterfly unit
	void NTT_Radix4_BU(ZZ &a,ZZ &b,ZZ &c,ZZ &d);
	
    void NTT_t(std::vector<ZZ> &A); //transform time domain to frequency domain
    void NTT_t_r4(std::vector<ZZ> &A); //transform time domain to frequency domain
	void NTT_t_r4_r2(std::vector<ZZ> &A); //mix radix 
	void NTT_pointwise_add(std::vector<ZZ> &A, std::vector<ZZ> b, std::vector<ZZ> c);
	void NTT_pointwise_sub(std::vector<ZZ> &A, std::vector<ZZ> b, std::vector<ZZ> c);
	void NTT_pointwise_mult(std::vector<ZZ> &A, std::vector<ZZ> b, std::vector<ZZ> c);
	void NTT_slow_t(std::vector<ZZ> &a, std::vector<ZZ> &A);
	void INTT_Radix4_BU(ZZ &a,ZZ &b,ZZ &c,ZZ &d);
	void INTT_t(std::vector<ZZ> &A);// Inverse NTT, frequecy domain to time domain
	void INTT_t_r4_r2(std::vector<ZZ> &A); //mix radix 
	void INTT_slow_t(std::vector<ZZ> &A, std::vector<ZZ> &a);
};

#endif
