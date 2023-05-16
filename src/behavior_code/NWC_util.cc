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
    this->Phi = 0;
	this->W = 0;
	this->IW = 0;
    Phi = find_phi(N, Modular);
	InvMod(InvPhi, Phi, Modular);
	PowerMod(W, Phi, 2, Modular);
	InvMod(IW, W, Modular);
    setValue(Radix_r1, Radix_r2, N, Modular, Phi, InvPhi, W, IW);
}

NWC_util::~NWC_util(){
	// TO DO
}

void NWC_util::setValue(long long Radix_r1, long long Radix_r2, long long N, 
						ZZ Modular, ZZ Phi, ZZ InvPhi, ZZ W, ZZ IW){
    this->Radix_r1 = Radix_r1;
    this->Radix_r2 = Radix_r2;
    this->N = N;
    this->Modular = Modular;
	this->Phi = Phi;
	this->InvPhi = InvPhi;
	this->W = W;
	this->IW = IW;

}
void NWC_util::getValue(long long *Radix_r1, long long *Radix_r2, long long *N, ZZ *Modular, 
						ZZ *Phi, ZZ *InvPhi, ZZ *W, ZZ *IW){
	*Radix_r1 = this->Radix_r1;
	*Radix_r2 = this->Radix_r2;
	*N = this->N;
	*Modular = this->Modular;
	*Phi = this->Phi;
	*InvPhi = this->InvPhi;
	*W = this->W;
	*IW = this->IW;
}

void NWC_util::showInfo(){
	long long Radix_r1, Radix_r2, N;
	ZZ Modular, Phi, InvPhi;
	ZZ W, IW;
	ZZ W1, W2, W3, W4, W5, W6, W7;
	ZZ IW1, IW2, IW3, IW4, IW5, IW6, IW7;
	getValue(&Radix_r1, &Radix_r2, &N, &Modular, &Phi, &InvPhi, &W, &IW);
	PowerMod(W1, W, 1, Modular);
	PowerMod(W2, W, 2, Modular);
	PowerMod(W3, W, 3, Modular);
	PowerMod(W4, W, 4, Modular);
	PowerMod(W5, W, 5, Modular);
	PowerMod(W6, W, 6, Modular);
	PowerMod(W7, W, 7, Modular);
	
	PowerMod(IW1, IW, 1, Modular);
	PowerMod(IW2, IW, 2, Modular);
	PowerMod(IW3, IW, 3, Modular);
	PowerMod(IW4, IW, 4, Modular);
	PowerMod(IW5, IW, 5, Modular);
	PowerMod(IW6, IW, 6, Modular);
	PowerMod(IW7, IW, 7, Modular);
	cout << "W1 = " << W1 << endl;
	cout << "W2 = " << W2 << endl;
	cout << "W3 = " << W3 << endl;
	cout << "W4 = " << W4 << endl;
	cout << "W5 = " << W5 << endl;
	cout << "W6 = " << W6 << endl;
	cout << "W7 = " << W7 << endl;

	cout << "IW1 = " << IW1 << endl;
	cout << "IW2 = " << IW2 << endl;
	cout << "IW3 = " << IW3 << endl;
	cout << "IW4 = " << IW4 << endl;
	cout << "IW5 = " << IW5 << endl;
	cout << "IW6 = " << IW6 << endl;
	cout << "IW7 = " << IW7 << endl;

	ZZ Phi1, Phi2, Phi3, Phi4, Phi5, Phi6, Phi7;
	ZZ InvPhi1, InvPhi2, InvPhi3, InvPhi4, InvPhi5, InvPhi6, InvPhi7, InvPhi8;
	PowerMod(Phi1, Phi, 1, Modular);
	PowerMod(Phi2, Phi, 2, Modular);
	PowerMod(Phi3, Phi, 3, Modular);
	PowerMod(Phi4, Phi, 4, Modular);
	PowerMod(Phi5, Phi, 5, Modular);
	PowerMod(Phi6, Phi, 6, Modular);
	PowerMod(Phi7, Phi, 7, Modular);

	PowerMod(InvPhi1, InvPhi, 1, Modular);
	PowerMod(InvPhi2, InvPhi, 2, Modular);
	PowerMod(InvPhi3, InvPhi, 3, Modular);
	PowerMod(InvPhi4, InvPhi, 4, Modular);
	PowerMod(InvPhi5, InvPhi, 5, Modular);
	PowerMod(InvPhi6, InvPhi, 6, Modular);
	PowerMod(InvPhi7, InvPhi, 7, Modular);
	PowerMod(InvPhi8, InvPhi, 8, Modular);
	cout << "------------------------" << endl;
	cout << "Phi1 = " << Phi1 << endl;
	cout << "Phi2 = " << Phi2 << endl;
	cout << "Phi3 = " << Phi3 << endl;
	cout << "Phi4 = " << Phi4 << endl;
	cout << "Phi5 = " << Phi5 << endl;
	cout << "Phi6 = " << Phi6 << endl;
	cout << "Phi7 = " << Phi7 << endl;

	cout << "InvPhi1 = " << InvPhi1 << endl;
	cout << "InvPhi2 = " << InvPhi2 << endl;
	cout << "InvPhi3 = " << InvPhi3 << endl;
	cout << "InvPhi4 = " << InvPhi4 << endl;
	cout << "InvPhi5 = " << InvPhi5 << endl;
	cout << "InvPhi6 = " << InvPhi6 << endl;
	cout << "InvPhi7 = " << InvPhi7 << endl;
	cout << "InvPhi8 = " << InvPhi8 << endl;
	
	cout << "+-----------------------------" << endl;
	cout << "| The Radix_r1 = " << Radix_r1 << endl;
	cout << "| The Radix_r2 = " << Radix_r2 << endl;
	cout << "| The degree N = " << N << endl;
	cout << "| The Modular  = " << Modular << endl;
	cout << "| The Phi factor = " << Phi << endl;
	cout << "| The Inverse Phi = " << InvPhi << endl;
	cout << "| The W = " << W << endl;
	cout << "| The IW = " << IW << endl;
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

