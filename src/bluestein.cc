#include <iostream>
#include <NTL/ZZ.h>
#include <time.h>
#include <math.h>
#include <fstream>
#include "bluestein.h"

using namespace NTL;

//root: 2n-th roots of unity
void bluestein::init(ZZ P,ZZ root,unsigned long length){
	std::vector<ZZ> h;
	prime = P;
	ROU   = root;
	n     = length;
	InvMod(IROU, ROU, prime);  //calculate inverse of ROU
	//calculate h
	h.resize(n); //resize the vector as length n
	
    std::ofstream h2_o;
    h2_o.open("h2_pre_reorder.txt");
	for(int i = 0; i < n;i++){    
		unsigned long exp;
		ZZ tmp;
		exp  = pow(i,2); // i ^ 2 
		exp  = exp % (2*n);
		PowerMod(tmp,IROU,exp,prime);
		h[i] = tmp;
        h2_o << h[i];
        h2_o << "\n";
	}
    h2_o.close();
	//for(int i=0;i <n ;i++){
	//   std::cout<<" h["<<i<<"]:"<<h[i];
	//}

	unsigned long L1;
	double L2;
	L1 = 2*n - 2;
	L2 = log2(L1);
	M  = ceil(L2);
	M  = pow(2,M); // Power of 2 length

	h2_freq.resize(M);
    h2_time.resize(M);
    
	for(int j=0; j < M;j++){
		if(j < n){
			h2_freq[j] = h[j];
            h2_time[j] = h[j];
		}
		else if( j > M-n){  //M-n+2 paper or M-n+1
			h2_freq[j] = h[M-j];
            h2_time[j] = h[M-j];
		}
	}
    
	//Init power of 2 parameter
    conv(Prime_pof2,"18446744069414584321");
    ZZ ROU_tmp_65536;
	conv(ROU_tmp_65536,"14603442835287214144");
    unsigned long  order;
    //order of 65536-th ROU
    order = 65536 / M ; //M - th root of unity generator parameter 
    ROU_pof2 = PowerMod(ROU_tmp_65536,order,Prime_pof2); //M-th root of unity
    
	
	FFT.NTT_init(M,Prime_pof2,ROU_pof2);
	FFT.NTT_t(h2_freq);
	
}

void bluestein::BFFT(std::vector<ZZ> &A){
     std::vector<ZZ> y;
	 y.resize(M);
	 for(int i = 0; i < n;i++){
		unsigned long exp;
		ZZ tmp;
		exp  = pow(i,2); // i ^ 2
		exp  = exp % (2*n);
		PowerMod(tmp,ROU,exp,prime);
		MulMod(tmp,A[i],tmp,prime);
		y[i] = tmp;		 
	 }

	 FFT.NTT_t(y);
	 
	 std::vector<ZZ> z;
	 z.resize(M);
	 FFT.NTT_pointwise_mult(z,y,h2_freq);
	 FFT.INTT_t(z);
	 for(int i = 0; i < n ; i++){
		 z[i] = z[i] % prime;
	 }
	 for(int i = 0; i < n;i++){
		 unsigned long exp;
		 ZZ tmp;
		 exp  = pow(i,2); // i ^ 2
		 PowerMod(tmp,ROU,exp,prime);
		 MulMod(tmp,z[i],tmp,prime);
		 A[i] = tmp;
	 }
	 
}

void bluestein::N_ROU(ZZ prime,unsigned long n,ZZ &root){
	ZZ exp;
	ZZ tmp;
	ZZ i;
	int k;
	exp = (prime-1) / n;
    for(i=2;i < prime; i=i+1){
		PowerMod(tmp,i,exp,prime);
		k = PROU(prime,tmp,n);
		if(k==0){
			root = tmp;
			i = prime;
		}
	}
}

int bluestein::PROU(ZZ prime,ZZ root,unsigned long n){
	ZZ tmp;
	int k = 0;
	for(unsigned long i = 0; i < n ;i++){
		PowerMod(tmp,root,i,prime);
		if(tmp == 1 && i!=0){
			k = 1; // it is not primitive n-th roots of unity
			i = n; // break the loop
		}
	}
	return k;
}