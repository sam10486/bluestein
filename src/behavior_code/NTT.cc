#include <iostream>
#include <fstream>
#include <NTL/ZZ.h>
#include <time.h>

#include "NTT.h"

using namespace NTL;

void NTT::NTT_init(unsigned long n, ZZ prime, ZZ root){
	N = n;
	ZZ N_ZZ;
	N_ZZ = N;
	p = prime;
	W = root;
	InvMod(IW, W, p); //calculate inverse of w
	InvMod(IN, N_ZZ, p); //calculate Inverse of N
}

void NTT::NTT_BU(ZZ &a,ZZ &b){
	ZZ tmp_a;
	ZZ tmp_b;
	AddMod(tmp_a, a, b, p);
	if (b < 0)b = b + p;
	SubMod(tmp_b, a, b, p);
	a = tmp_a;
	b = tmp_b;
}

void NTT::NTT_Radix4_BU(ZZ &a,ZZ &b,ZZ &c,ZZ &d){
	//a: x[n] b:[n+N/4] c:[n+N/2] d: 
	ZZ twiddle_1_4; // W^(N/4)
	ZZ twiddle_3_4; // W^(3N/4)
	unsigned long len_1_4; // N/4
	unsigned long len_3_4; // N/4
	len_1_4 = N / 4;
	PowerMod(twiddle_1_4, W, len_1_4, p);
	
	ZZ tmp_a;
	ZZ tmp_b;
    ZZ tmp_c;
    ZZ tmp_d;
	
	//stage 0
	AddMod(tmp_a,a,c,p);
	SubMod(tmp_c,a,c,p);
	AddMod(tmp_b,b,d,p);
	SubMod(tmp_d,b,d,p);
	MulMod(tmp_d,tmp_d,twiddle_1_4,p);
	
	
	ZZ tmp_a_1;
	ZZ tmp_b_1;
	ZZ tmp_c_1;
	ZZ tmp_d_1;
	
	//stage 1
	AddMod(tmp_a_1,tmp_a,tmp_b,p);
	SubMod(tmp_b_1,tmp_a,tmp_b,p);
	AddMod(tmp_c_1,tmp_c,tmp_d,p);
	SubMod(tmp_d_1,tmp_c,tmp_d,p);
	
	//data relocation  
	//bit-reverse
	a = tmp_a_1;
	c = tmp_b_1;
	b = tmp_c_1;
	d = tmp_d_1;
	
}

void NTT::INTT_Radix4_BU(ZZ &a,ZZ &b,ZZ &c,ZZ &d){
	//a: x[n] b:[n+N/4] c:[n+N/2] d: 
	ZZ twiddle_1_4; // W^(N/4)
	ZZ twiddle_3_4; // W^(3N/4)
	unsigned long len_1_4; // N/4
	unsigned long len_3_4; // N/4
	len_1_4 = N / 4;
	PowerMod(twiddle_1_4, IW, len_1_4, p);
	
	ZZ tmp_a;
	ZZ tmp_b;
    ZZ tmp_c;
    ZZ tmp_d;
	
	//stage 0
	AddMod(tmp_a,a,c,p);
	SubMod(tmp_c,a,c,p);
	AddMod(tmp_b,b,d,p);
	SubMod(tmp_d,b,d,p);
	MulMod(tmp_d,tmp_d,twiddle_1_4,p);
	
	
	ZZ tmp_a_1;
	ZZ tmp_b_1;
	ZZ tmp_c_1;
	ZZ tmp_d_1;
	
	//stage 1
	AddMod(tmp_a_1,tmp_a,tmp_b,p);
	SubMod(tmp_b_1,tmp_a,tmp_b,p);
	AddMod(tmp_c_1,tmp_c,tmp_d,p);
	SubMod(tmp_d_1,tmp_c,tmp_d,p);
	
	//data relocation  
	//bit-reverse
	a = tmp_a_1;
	c = tmp_b_1;
	b = tmp_c_1;
	d = tmp_d_1;
	
}
void NTT::NTT_t(std::vector<ZZ> &A) {
	unsigned long  S;//stage
	unsigned int   bias;
	unsigned int   nc; //number of class
	unsigned int   ne; //number of element for each class
	S = (unsigned long)ceil(log2(N));
	std::vector<std::vector<ZZ> > A_t;
	ZZ factor;   //base factor
	ZZ factor_t; //acctually mul factor
	bias = N;
	
	for (int i = 0; i < S; i++) { //DIT 
		nc = N / bias;
		ne = bias;
		A_t.resize(nc); //resize
		for (int k = 0; k < nc; k++) {
			A_t[k].resize(ne);
		}                              //resize
		for (unsigned int j = 0; j < nc; j++) {
			for (unsigned int jj = 0; jj < (ne); jj++) {
				A_t[j][jj] = A[j * ne + jj];
			} //given value
		}
		bias = bias >> 1; //divide by 2
		if (i == 0) factor = W;
		else SqrMod(factor, factor, p); 

		for (unsigned int s = 0; s < nc; s++) {
			for (unsigned int ss = 0; ss < (ne/2); ss++) {
				NTT_BU(A_t[s][ss], A_t[s][ss + bias]);
				PowerMod(factor_t, factor, ss, p);
				MulMod(A_t[s][ss + bias], A_t[s][ss + bias], factor_t, p);
			}
		}
		for (unsigned int x_index = 0; x_index < nc; x_index++) {
			for (unsigned int y_index = 0; y_index < ne; y_index++) {
				A[x_index*ne + y_index] = A_t[x_index][y_index];  //write back
			}
		}
	}
	//data relocation
	ZZ tmp;
	int exchange_position = 0;
	int bit=0;
	int bit_weight = 0;
	int p_tmp = 0;
	for (int i = 0; i < N; i++) {
		p_tmp = i;
		exchange_position = 0;
		for (int ls = 0; ls < S; ls++) {
			bit = p_tmp % 2;
			bit_weight = 1 << ((S-1) - ls);
			if (bit == 1) exchange_position = exchange_position + bit_weight;
			else  exchange_position = exchange_position;
			p_tmp = p_tmp >> 1;
		}
		// exchange_posistion > i , then data exchange 
		if (exchange_position > i) {
			tmp = A[i];
			A[i] = A[exchange_position];
			A[exchange_position] = tmp;
		}
	}

}

void NTT::NTT_t_r4(std::vector<ZZ> &A) {
	unsigned long  S;//stage
	unsigned int   bias;
	unsigned int   nc; //number of class
	unsigned int   ne; //number of element for each class
	S = (unsigned long)ceil(log2(N));
	S = S / 2 ; // stage for radix-4
	std::vector<std::vector<ZZ> > A_t;
	ZZ factor;   //base factor
	ZZ factor_t; //acctually mul factor
	ZZ factor_2t; //acctually mul factor
	ZZ factor_3t; //acctually mul factor
	bias = N;

    std::ofstream DataRecord("./NTT_R4_NOSPMB.txt");
	DataRecord << "Stage: "<< S<<"\n";
	for (int i = 0; i < S; i++) { //DIT 
		nc = N / bias;
		ne = bias;
		A_t.resize(nc); //resize
		for (int k = 0; k < nc; k++) {
			A_t[k].resize(ne);
		}                              //resize
		for (unsigned int j = 0; j < nc; j++) {
			for (unsigned int jj = 0; jj < (ne); jj++) {
				A_t[j][jj] = A[j * ne + jj];
			} //given value
		}
		bias = bias >> 2; //divide by 4
		if (i == 0) factor = W;
		else {
		    SqrMod(factor, factor, p); 
		    SqrMod(factor, factor, p); 
		}
		DataRecord <<"----------------------------------\n";
		DataRecord <<"now stage: " << i <<" \n";
		DataRecord <<"twiddle factor : " << factor <<" \n";
		DataRecord <<"*********\n";
		for (unsigned int s = 0; s < nc; s++) {
			for (unsigned int ss = 0; ss < (ne/4); ss++) {
				DataRecord <<"Before butterfly unit operation! \n";
				DataRecord <<"A_t["<< s <<"]["<< ss <<"]: " << A_t[s][ss] <<" \n";
				DataRecord <<"A_t["<< s <<"]["<< ss + 1 * bias <<"]: " << A_t[s][ss + 1 * bias]<<" \n";
				DataRecord <<"A_t["<< s <<"]["<< ss + 2 * bias <<"]: " << A_t[s][ss + 2 * bias] <<" \n";
				DataRecord <<"A_t["<< s <<"]["<< ss + 3 * bias <<"]: " << A_t[s][ss + 3 * bias] <<" \n";
				NTT_Radix4_BU(A_t[s][ss],A_t[s][ss + bias],A_t[s][ss + 2 * bias],A_t[s][ss + 3 * bias]);
				PowerMod(factor_t,  factor, 1 * ss, p);
				PowerMod(factor_2t, factor, 2 * ss, p);
				PowerMod(factor_3t, factor, 3 * ss, p);
				DataRecord <<"length:  " << ss <<" \n";
				DataRecord <<"twiddle factor_1 : " << factor_t <<" \n";
				DataRecord <<"twiddle factor_2 : " << factor_2t <<" \n";
				DataRecord <<"twiddle factor_3 : " << factor_3t <<" \n";
				MulMod(A_t[s][ss + 1 * bias], A_t[s][ss + 1 * bias], factor_t, p);
				MulMod(A_t[s][ss + 2 * bias], A_t[s][ss + 2 * bias], factor_2t, p);
				MulMod(A_t[s][ss + 3 * bias], A_t[s][ss + 3 * bias], factor_3t, p);
				DataRecord <<"A_t["<< s <<"]["<< ss <<"]: " << A_t[s][ss] <<" \n";
				DataRecord <<"A_t["<< s <<"]["<< ss + 1 * bias <<"]: " << A_t[s][ss + 1 * bias]<<" \n";
				DataRecord <<"A_t["<< s <<"]["<< ss + 2 * bias <<"]: " << A_t[s][ss + 2 * bias] <<" \n";
				DataRecord <<"A_t["<< s <<"]["<< ss + 3 * bias <<"]: " << A_t[s][ss + 3 * bias] <<" \n";
				DataRecord <<"-----------------------------------------------------------------------\n";
			}
		}
		for (unsigned int x_index = 0; x_index < nc; x_index++) {
			for (unsigned int y_index = 0; y_index < ne; y_index++) {
				A[x_index*ne + y_index] = A_t[x_index][y_index];  //write back
			}
		}
	}
	
	
	//data relocation
	int N_bit_size;  
	double d_tmp;
	d_tmp = (double)ceil(log2(N));
	N_bit_size = (int)d_tmp;
	
	std::cout << "N_bit_size : " << N_bit_size <<"\n";
	std::vector<int> bit_rep;
	std::vector<int> bit_rep_br;
	bit_rep.resize(N_bit_size);
	bit_rep_br.resize(N_bit_size);
	
	ZZ tmp;
	int exchange_position = 0;
	int bit=0;
	int bit_weight=0;
	int weight_tmp;
	int p_tmp = 0;
	
	for (int i = 0; i < N; i++) {
		p_tmp = i;
		exchange_position = 0;
		for(int bit_index = 0;bit_index < N_bit_size; bit_index++){
		   bit = p_tmp % 2;
		   bit_rep[bit_index] = bit;
		   p_tmp = p_tmp >> 1;
		}
        
		for(int bit_index=0; bit_index < (N_bit_size/2);bit_index++){
        	bit_rep_br[(N_bit_size-1) - (2 * bit_index + 0)] = bit_rep[2*bit_index + 1];
        	bit_rep_br[(N_bit_size-1) - (2 * bit_index + 1)] = bit_rep[2*bit_index + 0];		 
		}
	    
		for(int bit_index = 0; bit_index < N_bit_size;bit_index++){
			weight_tmp = 1 << bit_index;
			if(bit_rep_br[bit_index] == 1) bit_weight = weight_tmp;
			else bit_weight = 0;		
			exchange_position = exchange_position + bit_weight;
		}

		//std::cout <<"**********************************\n";
		//exchange_position > i , then data exchange 
		if (exchange_position > i) {
			tmp = A[i];
			A[i] = A[exchange_position];
			A[exchange_position] = tmp;
		}
		
	}
    
}
//Mixed radix-r4 and radix-r2
//final stage for r2
void NTT::NTT_t_r4_r2(std::vector<ZZ> &A) {
	unsigned long  S;//stage
	double         S_tmp;
	unsigned int   bias;
	unsigned int   nc; //number of class
	unsigned int   ne; //number of element for each class
	std::ofstream DataRecord("./NTT_R4_R2_NOSPMB.txt");
	//radix-4 
	S = (unsigned long)(log2(N));
	S_tmp = S/2;
	S_tmp = floor(S_tmp);
	S = (unsigned long) S_tmp; // stage for radix-4
	std::vector<std::vector<ZZ> > A_t;
	ZZ factor;   //base factor
	ZZ factor_t; //acctually mul factor
	ZZ factor_2t; //acctually mul factor
	ZZ factor_3t; //acctually mul factor
	bias = N;

	//radix-4 DIT
	for (int i = 0; i < S; i++) { //DIT 
		nc = N / bias;
		ne = bias;
		A_t.resize(nc); //resize
		for (int k = 0; k < nc; k++) {
			A_t[k].resize(ne);
		}                              //resize
		for (unsigned int j = 0; j < nc; j++) {
			for (unsigned int jj = 0; jj < (ne); jj++) {
				A_t[j][jj] = A[j * ne + jj];
			} //given value
		}
		bias = bias >> 2; //divide by 4
		if (i == 0) factor = W;
		else {
		    SqrMod(factor, factor, p); 
		    SqrMod(factor, factor, p); 
		}
		DataRecord << "--------------------------------------------\n";
		DataRecord << "stage : " << i <<"\n";
		DataRecord << "factor : " << factor <<"\n";
		DataRecord << "********************************\n";
		for (unsigned int s = 0; s < nc; s++) {
			for (unsigned int ss = 0; ss < (ne/4); ss++) {
				DataRecord <<"Before butterfly unit operation! \n";
				DataRecord <<"A_t["<< s <<"]["<< ss <<"]: " << A_t[s][ss] <<" \n";
				DataRecord <<"A_t["<< s <<"]["<< ss + 1 * bias <<"]: " << A_t[s][ss + 1 * bias]<<" \n";
				DataRecord <<"A_t["<< s <<"]["<< ss + 2 * bias <<"]: " << A_t[s][ss + 2 * bias] <<" \n";
				DataRecord <<"A_t["<< s <<"]["<< ss + 3 * bias <<"]: " << A_t[s][ss + 3 * bias] <<" \n";
				NTT_Radix4_BU(A_t[s][ss],A_t[s][ss + bias],A_t[s][ss + 2 * bias],A_t[s][ss + 3 * bias]);
				PowerMod(factor_t,  factor, 1 * ss, p);
				PowerMod(factor_2t, factor, 2 * ss, p);
				PowerMod(factor_3t, factor, 3 * ss, p);
                DataRecord <<"length:  " << ss <<" \n";
                DataRecord <<"twiddle factor_1 : " << factor_t <<" \n";
                DataRecord <<"twiddle factor_2 : " << factor_2t <<" \n";
                DataRecord <<"twiddle factor_3 : " << factor_3t <<" \n";			
				MulMod(A_t[s][ss + 1 * bias], A_t[s][ss + 1 * bias], factor_t, p);
				MulMod(A_t[s][ss + 2 * bias], A_t[s][ss + 2 * bias], factor_2t, p);
				MulMod(A_t[s][ss + 3 * bias], A_t[s][ss + 3 * bias], factor_3t, p);
                DataRecord <<"A_t["<< s <<"]["<< ss <<"]: " << A_t[s][ss] <<" \n";
                DataRecord <<"A_t["<< s <<"]["<< ss + 1 * bias <<"]: " << A_t[s][ss + 1 * bias]<<" \n";				
                DataRecord <<"A_t["<< s <<"]["<< ss + 2 * bias <<"]: " << A_t[s][ss + 2 * bias] <<" \n";				
                DataRecord <<"A_t["<< s <<"]["<< ss + 3 * bias <<"]: " << A_t[s][ss + 3 * bias] <<" \n";				
                DataRecord <<"-----------------------------------------------------------------------\n";				
			}
		}
		for (unsigned int x_index = 0; x_index < nc; x_index++) {
			for (unsigned int y_index = 0; y_index < ne; y_index++) {
				A[x_index*ne + y_index] = A_t[x_index][y_index];  //write back
			}
		}
	}
    for(int i = 0;i < N/2;i++){
		DataRecord <<"A["<< 2*i    <<"]:"<< A[2*i]    <<" \n";
		DataRecord <<"A["<< 2*i +1 <<"]:"<< A[2*i + 1]<<" \n";					
		NTT_BU(A[2*i],A[2*i+1]);
		DataRecord <<"A["<< 2*i    <<"]:"<< A[2*i]    <<" \n";
		DataRecord <<"A["<< 2*i +1 <<"]:"<< A[2*i + 1]<<" \n";
        DataRecord <<"-----------------------------------------------------------------------\n";						
	}
	
	//data relocation
	int N_bit_size;  
	double d_tmp;
	d_tmp = (double)ceil(log2(N));
	N_bit_size = (int)d_tmp;
	
	std::cout << "N_bit_size : " << N_bit_size <<"\n";
	std::vector<int> bit_rep;
	std::vector<int> bit_rep_br;
	bit_rep.resize(N_bit_size);
	bit_rep_br.resize(N_bit_size);
	
	ZZ tmp;
	int exchange_position = 0;
	int bit=0;
	int bit_weight=0;
	int weight_tmp;
	int p_tmp = 0;
	
	std::vector<ZZ> A_tmp;
	A_tmp.resize(N);
	
	for (int i = 0; i < N; i++) {
		p_tmp = i;
		exchange_position = 0;
		for(int bit_index = 0;bit_index < N_bit_size; bit_index++){
		   bit = p_tmp % 2;
		   bit_rep[bit_index] = bit;
		   p_tmp = p_tmp >> 1;
		}
        
		//radix-4 and final_stage radix-2 DATA relocation as 
		//output order lsb = new MSB
		//other bit reverse
		//example [a b c d e] =====> [ e c d a b]
		for(int bit_index=0; bit_index < ((N_bit_size-1)/2);bit_index++){
        	bit_rep_br[(N_bit_size-1) - (2 * bit_index + 0)-1] = bit_rep[2*bit_index + 2];
        	bit_rep_br[(N_bit_size-1) - (2 * bit_index + 1)-1] = bit_rep[2*bit_index + 1];		 
		}
	    bit_rep_br[N_bit_size-1] = bit_rep[0];
		for(int bit_index = 0; bit_index < N_bit_size;bit_index++){
			weight_tmp = 1 << bit_index;
			if(bit_rep_br[bit_index] == 1) bit_weight = weight_tmp;
			else bit_weight = 0;		
			exchange_position = exchange_position + bit_weight;
		}
		A_tmp[exchange_position] = A[i];
		//std::cout <<"**********************************\n";
		//exchange_position > i , then data exchange 
		
		//if (exchange_position > i) {
		//	tmp = A[i];
		//	A[i] = A[exchange_position];
		//	A[exchange_position] = tmp;
		//}	
	}	
	for(int i=0;i<N;i++){
		A[i] = A_tmp[i];
	}
    
}

void NTT::INTT_t(std::vector<ZZ> &A) {
	unsigned long  S;//stage
	unsigned int   bias;
	unsigned int   nc; //number of class
	unsigned int   ne; //number of element for each class
	S = (unsigned long)ceil(log2(N));
	std::vector<std::vector<ZZ> > A_t;
	std::ofstream DataRecord("./INTT_output.txt");
	ZZ factor;   //base factor
	ZZ factor_t; //acctually mul factor
	bias = N;

	for (int i = 0; i < S; i++) {
		nc = N / bias;
		ne = bias;
		A_t.resize(nc); //resize
		DataRecord << "*********************************************\n";
		DataRecord << "stage: " << i << "\n";
		for (int k = 0; k < nc; k++) {
			A_t[k].resize(ne);
		}                              //resize
		for (unsigned int j = 0; j < nc; j++) {
			for (unsigned int jj = 0; jj < (ne); jj++) {
				A_t[j][jj] = A[j * ne + jj];
			} //given value
		}

		bias = bias >> 1; //divide by 2
		if (i == 0) factor = IW;
		else SqrMod(factor, factor, p);

		for (unsigned int s = 0; s < nc; s++) {
			for (unsigned int ss = 0; ss < (ne / 2); ss++) {
				DataRecord << "Before computing!!! \n";
                DataRecord << "A_t["<< s << "]["<< ss <<"]: "<< A_t[s][ss] << "\n";
                DataRecord << "A_t["<< s << "]["<< ss + bias <<"]: "<< A_t[s][ss + bias] << "\n";
				NTT_BU(A_t[s][ss], A_t[s][ss + bias]);
				DataRecord << "After Butterfly unit computing!!! \n";
                DataRecord << "A_t["<< s << "]["<< ss <<"]: "<< A_t[s][ss] << "\n";
                DataRecord << "A_t["<< s << "]["<< ss + bias <<"]: "<< A_t[s][ss + bias] << "\n";				
				PowerMod(factor_t, factor, ss, p);
				DataRecord << "ss: "<< ss << "\n";
				DataRecord << "factor: "<< factor << "\n";
				DataRecord << "factor_t: "<< factor_t << "\n";
				MulMod(A_t[s][ss + bias], A_t[s][ss + bias], factor_t, p);
				DataRecord << "After computing!!! \n";
				DataRecord << "A_t["<< s << "]["<< ss <<"]: "<< A_t[s][ss] << "\n";
				DataRecord << "A_t["<< s << "]["<< ss + bias <<"]: "<< A_t[s][ss + bias] << "\n";
				DataRecord << "--------------------------------\n";
			}
		}
		for (unsigned int x_index = 0; x_index < nc; x_index++) {
			for (unsigned int y_index = 0; y_index < ne; y_index++) {
				A[x_index*ne + y_index] = A_t[x_index][y_index];  //write back
			}
		}
	}
	//data relocation
	ZZ tmp;
	int exchange_position = 0;
	int bit = 0;
	int bit_weight = 0;
	int p_tmp = 0;
	for (int i = 0; i < N; i++) {
		p_tmp = i;
		exchange_position = 0;
		for (int ls = 0; ls < S; ls++) {
			bit = p_tmp % 2;
			bit_weight = 1 << ((S - 1) - ls);
			if (bit == 1) exchange_position = exchange_position + bit_weight;
			else  exchange_position = exchange_position;
			p_tmp = p_tmp >> 1;
		}
		// exchange_posistion > i , then data exchange 
		if (exchange_position > i) {
			tmp = A[i];
			A[i] = A[exchange_position];
			A[exchange_position] = tmp;
		}
	}
	for (int i = 0; i < N; i++) {
		MulMod(A[i], A[i], IN, p);
	}
}

void NTT::INTT_t_r4_r2(std::vector<ZZ> &A) {
	unsigned long  S;//stage
	double         S_tmp;
	unsigned int   bias;
	unsigned int   nc; //number of class
	unsigned int   ne; //number of element for each class
	std::ofstream  DataRecord("./INTT_R4_R2_NOSPMB.txt");
	//radix-4 
	S = (unsigned long)(log2(N));
	S_tmp = S/2;
	S_tmp = floor(S_tmp);
	S = (unsigned long) S_tmp; // stage for radix-4
	std::vector<std::vector<ZZ> > A_t;
	ZZ factor;   //base factor
	ZZ factor_t; //acctually mul factor
	ZZ factor_2t; //acctually mul factor
	ZZ factor_3t; //acctually mul factor
	bias = N;

	//radix-4 DIT
	for (int i = 0; i < S; i++) { //DIT 
		nc = N / bias;
		ne = bias;
		A_t.resize(nc); //resize
		for (int k = 0; k < nc; k++) {
			A_t[k].resize(ne);
		}                              //resize
		for (unsigned int j = 0; j < nc; j++) {
			for (unsigned int jj = 0; jj < (ne); jj++) {
				A_t[j][jj] = A[j * ne + jj];
			} //given value
		}
		bias = bias >> 2; //divide by 4
		if (i == 0) factor = IW;
		else {
		    SqrMod(factor, factor, p); 
		    SqrMod(factor, factor, p); 
		}
		DataRecord << "--------------------------------------------\n";
		DataRecord << "stage : " << i <<"\n";
		DataRecord << "factor : " << factor <<"\n";
		DataRecord << "********************************\n";
		for (unsigned int s = 0; s < nc; s++) {
			for (unsigned int ss = 0; ss < (ne/4); ss++) {
				DataRecord <<"Before butterfly unit operation! \n";
				DataRecord <<"A_t["<< s <<"]["<< ss <<"]: " << A_t[s][ss] <<" \n";
				DataRecord <<"A_t["<< s <<"]["<< ss + 1 * bias <<"]: " << A_t[s][ss + 1 * bias]<<" \n";
				DataRecord <<"A_t["<< s <<"]["<< ss + 2 * bias <<"]: " << A_t[s][ss + 2 * bias] <<" \n";
				DataRecord <<"A_t["<< s <<"]["<< ss + 3 * bias <<"]: " << A_t[s][ss + 3 * bias] <<" \n";
				INTT_Radix4_BU(A_t[s][ss],A_t[s][ss + bias],A_t[s][ss + 2 * bias],A_t[s][ss + 3 * bias]);
				DataRecord <<"After radix-4 butterfly unit  \n";
			    DataRecord <<"A_t["<< s <<"]["<< ss <<"]: " << A_t[s][ss] <<" \n";
				DataRecord <<"A_t["<< s <<"]["<< ss + 1 * bias <<"]: " << A_t[s][ss + 1 * bias]<<" \n";
				DataRecord <<"A_t["<< s <<"]["<< ss + 2 * bias <<"]: " << A_t[s][ss + 2 * bias] <<" \n";
				DataRecord <<"A_t["<< s <<"]["<< ss + 3 * bias <<"]: " << A_t[s][ss + 3 * bias] <<" \n";
				PowerMod(factor_t,  factor, 1 * ss, p);
				PowerMod(factor_2t, factor, 2 * ss, p);
				PowerMod(factor_3t, factor, 3 * ss, p);
				DataRecord <<"factor:  " << factor <<" \n";
                DataRecord <<"length:  " << ss <<" \n";
                DataRecord <<"twiddle factor_1 : " << factor_t <<" \n";
                DataRecord <<"twiddle factor_2 : " << factor_2t <<" \n";
                DataRecord <<"twiddle factor_3 : " << factor_3t <<" \n";			
				MulMod(A_t[s][ss + 1 * bias], A_t[s][ss + 1 * bias], factor_t, p);
				MulMod(A_t[s][ss + 2 * bias], A_t[s][ss + 2 * bias], factor_2t, p);
				MulMod(A_t[s][ss + 3 * bias], A_t[s][ss + 3 * bias], factor_3t, p);
                DataRecord <<"After Mult twiddle factor  \n";
				DataRecord <<"A_t["<< s <<"]["<< ss <<"]: " << A_t[s][ss] <<" \n";
                DataRecord <<"A_t["<< s <<"]["<< ss + 1 * bias <<"]: " << A_t[s][ss + 1 * bias]<<" \n";				
                DataRecord <<"A_t["<< s <<"]["<< ss + 2 * bias <<"]: " << A_t[s][ss + 2 * bias] <<" \n";				
                DataRecord <<"A_t["<< s <<"]["<< ss + 3 * bias <<"]: " << A_t[s][ss + 3 * bias] <<" \n";				
                DataRecord <<"-----------------------------------------------------------------------\n";				
			}
		}
		for (unsigned int x_index = 0; x_index < nc; x_index++) {
			for (unsigned int y_index = 0; y_index < ne; y_index++) {
				A[x_index*ne + y_index] = A_t[x_index][y_index];  //write back
			}
		}
	}
    for(int i = 0;i < N/2;i++){
		DataRecord <<"A["<< 2*i    <<"]:"<< A[2*i]    <<" \n";
		DataRecord <<"A["<< 2*i +1 <<"]:"<< A[2*i + 1]<<" \n";					
		NTT_BU(A[2*i],A[2*i+1]);
		DataRecord <<"A["<< 2*i    <<"]:"<< A[2*i]    <<" \n";
		DataRecord <<"A["<< 2*i +1 <<"]:"<< A[2*i + 1]<<" \n";
        DataRecord <<"-----------------------------------------------------------------------\n";						
	}
	//divide N
	for(int i = 0; i < N;i++){
		MulMod(A[i],A[i],IN,p);
	}
	
	
	//data relocation
	int N_bit_size;  
	double d_tmp;
	d_tmp = (double)ceil(log2(N));
	N_bit_size = (int)d_tmp;
	
	std::cout << "N_bit_size : " << N_bit_size <<"\n";
	std::vector<int> bit_rep;
	std::vector<int> bit_rep_br;
	bit_rep.resize(N_bit_size);
	bit_rep_br.resize(N_bit_size);
	
	ZZ tmp;
	int exchange_position = 0;
	int bit=0;
	int bit_weight=0;
	int weight_tmp;
	int p_tmp = 0;
	
	std::vector<ZZ> A_tmp;
	A_tmp.resize(N);
	
	for (int i = 0; i < N; i++) {
		p_tmp = i;
		exchange_position = 0;
		for(int bit_index = 0;bit_index < N_bit_size; bit_index++){
		   bit = p_tmp % 2;
		   bit_rep[bit_index] = bit;
		   p_tmp = p_tmp >> 1;
		}
        
		//radix-4 and final_stage radix-2 DATA relocation as 
		//output order lsb = new MSB
		//other bit reverse
		//example [a b c d e] =====> [ e c d a b]
		for(int bit_index=0; bit_index < ((N_bit_size-1)/2);bit_index++){
        	bit_rep_br[(N_bit_size-1) - (2 * bit_index + 0)-1] = bit_rep[2*bit_index + 2];
        	bit_rep_br[(N_bit_size-1) - (2 * bit_index + 1)-1] = bit_rep[2*bit_index + 1];		 
		}
	    bit_rep_br[N_bit_size-1] = bit_rep[0];
		for(int bit_index = 0; bit_index < N_bit_size;bit_index++){
			weight_tmp = 1 << bit_index;
			if(bit_rep_br[bit_index] == 1) bit_weight = weight_tmp;
			else bit_weight = 0;		
			exchange_position = exchange_position + bit_weight;
		}
		A_tmp[exchange_position] = A[i];
		//std::cout <<"**********************************\n";
		//exchange_position > i , then data exchange 
		
		//if (exchange_position > i) {
		//	tmp = A[i];
		//	A[i] = A[exchange_position];
		//	A[exchange_position] = tmp;
		//}	
	}	
	for(int i=0;i<N;i++){
		A[i] = A_tmp[i];
	}   
}

void NTT::NTT_pointwise_add(std::vector<ZZ> &A, std::vector<ZZ> b, std::vector<ZZ> c) {
	if(A.size() != N)A.resize(N);
	for (unsigned long i = 0; i < N; i++) {
		AddMod(A[i], b[i], c[i], p);
	}
}

void NTT::NTT_pointwise_sub(std::vector<ZZ> &A, std::vector<ZZ> b, std::vector<ZZ> c) {
	if(A.size() != N)A.resize(N);
	for (unsigned long i = 0; i < N; i++) {
		SubMod(A[i], b[i], c[i], p);
	}
}

void NTT::NTT_pointwise_mult(std::vector<ZZ> &A, std::vector<ZZ> b, std::vector<ZZ> c) {
	if(A.size() != N) A.resize(N); //resize A arrary to N point
	for (long i = 0; i < N; i++) {
		MulMod(A[i], b[i], c[i], p); 
	}
}

void NTT::NTT_slow_t(std::vector<ZZ> &a, std::vector<ZZ> &A) {
	ZZ tw;
	ZZ factor;
	ZZ tmp;
	A.resize(N); //resize
	for (long i = 0; i < N; i++) {
		PowerMod(tw, W, i, p);
		A[i] = 0;
		tmp  = 0 ;
		for (long j = 0; j < N; j++) {
			PowerMod(factor, tw, j, p);
			MulMod(tmp, a[j], factor, p);
			AddMod(A[i], A[i], tmp, p);
		}
	}
}

void NTT::INTT_slow_t(std::vector<ZZ> &A, std::vector<ZZ> &a) {
	ZZ tw;
	ZZ factor;
	ZZ tmp;
	a.resize(N); //resize
	for (long i = 0; i < N; i++) {
		PowerMod(tw, IW, i, p);
		A[i] = 0;
		if ((i % 1024) == 0)std::cout << "run" << std::endl;
		for (long j = 0; j < N; j++) {
			PowerMod(factor, tw, j, p);
			MulMod(tmp, A[j], factor, p);
			AddMod(a[i], a[i], tmp, p);
		}
		MulMod(a[i], a[i], IN, p);
	}
}