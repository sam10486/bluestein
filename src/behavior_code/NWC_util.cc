#include <iostream>
#include <NTL/ZZ.h>
#include <vector>
#include "NWC_util.h"
#include "assert.h"

using namespace NTL;
using namespace std;

NWC_util::NWC_util(long long Radix_r1, long long Radix_r2, long long N, ZZ Modular){
    this->Radix_r1 = 0;
    this->Radix_r2 = 0;
    this->N = 0;
    this->Modular = 0;
    this->W = 0;
    ZZ W, IW;
    W = find_phi(N, Modular);
	InvMod(IW, W, Modular);
    setValue(Radix_r1, Radix_r2, N, Modular, W, IW);
}

NWC_util::~NWC_util(){
	// TO DO
}

void NWC_util::setValue(long long Radix_r1, long long Radix_r2, long long N, ZZ Modular, ZZ W, ZZ IW){
    this->Radix_r1 = Radix_r1;
    this->Radix_r2 = Radix_r2;
    this->N = N;
    this->Modular = Modular;
	this->W = W;
	this->IW = IW;

}
void NWC_util::getValue(long long *Radix_r1, long long *Radix_r2, long long *N, ZZ *Modular, ZZ *W, ZZ *IW){
	*Radix_r1 = this->Radix_r1;
	*Radix_r2 = this->Radix_r2;
	*N = this->N;
	*Modular = this->Modular;
	*W = this->W;
	*IW = this->IW;
}

void NWC_util::showInfo(){
	long long Radix_r1, Radix_r2, N;
	ZZ Modular, W, IW;
	getValue(&Radix_r1, &Radix_r2, &N, &Modular, &W, &IW);
	cout << "+-----------------------------" << endl;
	cout << "| The Radix_r1 = " << Radix_r1 << endl;
	cout << "| The Radix_r2 = " << Radix_r2 << endl;
	cout << "| The degree N = " << N << endl;
	cout << "| The Modular  = " << Modular << endl;
	cout << "| The twiddle factor = " << W << endl;
	cout << "| The Inverse twiddle factor = " << IW << endl;
	cout << "+-----------------------------" << endl;
}

// a^(p-1) = 1 (mod p)  ---> base^(Modular-1) = 1 (mod Modular)
ZZ NWC_util::find_n_rou(ZZ base, long long N, ZZ Modular) {
	assert(( Modular % N ) == 1);
	ZZ i;
	ZZ n_rou;
	i = (Modular-1)/N ;   // base^(Modular - 1) = base^( n * i ) = (base^i)^n = 1 (mod Modular)
	PowerMod(n_rou, base, i, Modular);
	return n_rou;
}

ZZ NWC_util::find_phi(long long N, ZZ Modular){   
	bool is_prou = false;
	ZZ i = (ZZ)2 ;
	ZZ n_rou;
	ZZ prou;
    long long degree = 2 * N;
	while(is_prou == false) {
		n_rou = find_n_rou(i, degree, Modular);
		is_prou = check_prou(n_rou, degree, Modular);
		i = i + 1;		
	}
	prou = n_rou;
	return prou;
}

//check if n_rou^1, n_rou^2,...,n_rou^(N-1) is not equal 1;
bool NWC_util::check_prou(ZZ n_rou, long long N, ZZ Modular){ 
	bool is_prou = true;
	ZZ tmp;
	for(int i = 1; i < N; i++){
		PowerMod(tmp, n_rou, i, Modular);
		if(tmp == 1){
			is_prou = false;
			break;
		}
	}
	return is_prou;
}

