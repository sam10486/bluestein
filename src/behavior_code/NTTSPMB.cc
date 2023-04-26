 #include <iostream>
#include <fstream>
#include <NTL/ZZ.h>
#include <time.h>

#include "NTTSPMB.h"
#include "SPMB.h"

#include <vector>
#include <algorithm>
#include "BitOperate.h"

using namespace NTL;

void NTTSPMB::init(unsigned long n, ZZ prime, ZZ root,int r){
	N = n;
	radix = r;
	ZZ N_ZZ;
	N_ZZ = N;
	p = prime;
	W = root;
	InvMod(IW, W, p); //calculate inverse of w
	InvMod(IN, N_ZZ, p); //calculate Inverse of N
}

void NTTSPMB::Radix2_BU(ZZ &a,ZZ &b){
	ZZ tmp_a;
	ZZ tmp_b;
	AddMod(tmp_a, a, b, p);
	if (b < 0)b = b + p;
	SubMod(tmp_b, a, b, p);
	a = tmp_a;
	b = tmp_b;
}

//radix 2^(2) DIF
void NTTSPMB::Radix4_BU(ZZ &a,ZZ &b,ZZ &c,ZZ &d){
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
	
	//std::cout << "tmp_a = " << tmp_a << ", a = " << a << ", c = " << c << std::endl;
	//std::cout << "tmp_c = " << tmp_c << ", a = " << a << ", c = " << c << std::endl;
	//std::cout << "tmp_b = " << tmp_b << ", b = " << b << ", d = " << d << std::endl;
	//std::cout << "tmp_a = " << tmp_d << ", b = " << b << ", d = " << d << std::endl;

	
	
	ZZ tmp_a_1;
	ZZ tmp_b_1;
	ZZ tmp_c_1;
	ZZ tmp_d_1;
	
	//stage 1
	AddMod(tmp_a_1,tmp_a,tmp_b,p);
	SubMod(tmp_b_1,tmp_a,tmp_b,p);
	AddMod(tmp_c_1,tmp_c,tmp_d,p);
	SubMod(tmp_d_1,tmp_c,tmp_d,p);

	//std::cout << "tmp_a_1 = " << tmp_a_1 << ", mp_a = " << tmp_a << ", tmp_b = " << tmp_b << std::endl;
	//std::cout << "tmp_b_1 = " << tmp_b_1 << ", mp_a = " << tmp_a << ", tmp_b = " << tmp_b << std::endl;
	//std::cout << "tmp_c_1 = " << tmp_c_1 << ", mp_c = " << tmp_c << ", tmp_d = " << tmp_d << std::endl;
	//std::cout << "tmp_d_1 = " << tmp_d_1 << ", mp_c = " << tmp_c << ", tmp_d = " << tmp_d << std::endl;
	
	//data relocation  
	//bit-reverse
	a = tmp_a_1;
	c = tmp_b_1;
	b = tmp_c_1;
	d = tmp_d_1;
	
	//a = tmp_a_1;
	//b = tmp_b_1;
	//c = tmp_c_1;
	//d = tmp_d_1;
	
}
//radix 2^(2) DIF
void NTTSPMB::Radix4_BU_INTT(ZZ &a,ZZ &b,ZZ &c,ZZ &d){
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

void NTTSPMB::Radix8_BU(ZZ &a_r0,ZZ &a_r1,ZZ &a_r2,ZZ &a_r3,ZZ &a_r4,ZZ &a_r5,ZZ &a_r6,ZZ &a_r7){
	ZZ twiddle_1_8; // W^(N/8)
	ZZ twiddle_r8_nk;
	ZZ data_tmp;
	unsigned long len_1_8; // N/16
	
	len_1_8 = N/8;
	PowerMod(twiddle_1_8,W,len_1_8,p);
	std::vector<ZZ> A_in;
	std::vector<ZZ> A_out;
		
	A_in.resize(8);
	A_out.resize(8);
	//data input 
	A_in[0]  = a_r0;
	A_in[1]  = a_r1;
	A_in[2]  = a_r2;
	A_in[3]  = a_r3;
	A_in[4]  = a_r4;
	A_in[5]  = a_r5;
	A_in[6]  = a_r6;
	A_in[7]  = a_r7;
	
	long  nk_exp;
	nk_exp   = 0;
	data_tmp = 0;
	twiddle_r8_nk = 0;
	
	for(int i=0;i<8;i++){
		A_out[i] = 0;
	}
	
	for(long i=0;i<8;i++){
		for(long j=0;j<8;j++){
		    nk_exp = i * j;
			nk_exp = nk_exp % 8;
			PowerMod(twiddle_r8_nk,twiddle_1_8,nk_exp, p);
			MulMod(data_tmp,A_in[j],twiddle_r8_nk,p);
			AddMod(A_out[i],A_out[i],data_tmp,p);
		}
	}
	
	//data output
	a_r0 = A_out[0];
	a_r1 = A_out[1];
	a_r2 = A_out[2];
	a_r3 = A_out[3];
	a_r4 = A_out[4];
	a_r5 = A_out[5];
	a_r6 = A_out[6];
	a_r7 = A_out[7];
}

void NTTSPMB::Radix8_BU_INTT(ZZ &a_r0,ZZ &a_r1,ZZ &a_r2,ZZ &a_r3,ZZ &a_r4,ZZ &a_r5,ZZ &a_r6,ZZ &a_r7){
	ZZ twiddle_1_8; // W^(N/8)
	ZZ twiddle_r8_nk;
	ZZ data_tmp;
	unsigned long len_1_8; // N/16
	
	len_1_8 = N/8;
	PowerMod(twiddle_1_8,IW,len_1_8,p);
	std::vector<ZZ> A_in;
	std::vector<ZZ> A_out;
		
	A_in.resize(8);
	A_out.resize(8);
	//data input 
	A_in[0]  = a_r0;
	A_in[1]  = a_r1;
	A_in[2]  = a_r2;
	A_in[3]  = a_r3;
	A_in[4]  = a_r4;
	A_in[5]  = a_r5;
	A_in[6]  = a_r6;
	A_in[7]  = a_r7;
	
	long  nk_exp;
	nk_exp = 0;
	
	for(int i=0;i<8;i++){
		A_out[i] = 0;
	}
	
	for(long i=0;i<8;i++){
		for(long j=0;j<8;j++){
		    nk_exp = i * j;
			nk_exp = nk_exp % 8;
			PowerMod(twiddle_r8_nk,twiddle_1_8,nk_exp, p);
			MulMod(data_tmp,A_in[j],twiddle_r8_nk,p);
			AddMod(A_out[i],A_out[i],data_tmp,p);
		}
	}
	
	//data output
	a_r0 = A_out[0];
	a_r1 = A_out[1];
	a_r2 = A_out[2];
	a_r3 = A_out[3];
	a_r4 = A_out[4];
	a_r5 = A_out[5];
	a_r6 = A_out[6];
	a_r7 = A_out[7];
}
void NTTSPMB::Radix16_BU(ZZ &a_r0,ZZ &a_r1,ZZ &a_r2,ZZ &a_r3,ZZ &a_r4,ZZ &a_r5,ZZ &a_r6,ZZ &a_r7,
ZZ &a_r8,ZZ &a_r9,ZZ &a_r10,ZZ &a_r11,ZZ &a_r12,ZZ &a_r13,ZZ &a_r14,ZZ &a_r15){
	ZZ twiddle_1_16; // W^(N/8)
	ZZ twiddle_r16_nk;
	ZZ data_tmp;
	unsigned long len_1_16; // N/16
	
	len_1_16 = N/16;
	PowerMod(twiddle_1_16,W,len_1_16,p);
	std::vector<ZZ> A_in;
	std::vector<ZZ> A_out;
		
	A_in.resize(16);
	A_out.resize(16);
	//data input 
	A_in[0]  = a_r0;
	A_in[1]  = a_r1;
	A_in[2]  = a_r2;
	A_in[3]  = a_r3;
	A_in[4]  = a_r4;
	A_in[5]  = a_r5;
	A_in[6]  = a_r6;
	A_in[7]  = a_r7;
	A_in[8]  = a_r8;
	A_in[9]  = a_r9;
	A_in[10] = a_r10;
	A_in[11] = a_r11;
	A_in[12] = a_r12;
	A_in[13] = a_r13;
	A_in[14] = a_r14;
	A_in[15] = a_r15;
	
	long  nk_exp;
	nk_exp = 0;
	
	for(int i=0;i<16;i++){
		A_out[i] = 0;
	}
	
	for(long i=0;i<16;i++){
		for(long j=0;j<16;j++){
		    nk_exp = i * j;
			nk_exp = nk_exp % 16;
			PowerMod(twiddle_r16_nk,twiddle_1_16,nk_exp, p);
			MulMod(data_tmp,A_in[j],twiddle_r16_nk,p);
			AddMod(A_out[i],A_out[i],data_tmp,p);
		}
	}
	
	//data output
	a_r0 = A_out[0];
	a_r1 = A_out[1];
	a_r2 = A_out[2];
	a_r3 = A_out[3];
	a_r4 = A_out[4];
	a_r5 = A_out[5];
	a_r6 = A_out[6];
	a_r7 = A_out[7];
	a_r8 = A_out[8];
	a_r9 = A_out[9];
	a_r10 = A_out[10];
	a_r11 = A_out[11];
	a_r12 = A_out[12];
	a_r13 = A_out[13];
	a_r14 = A_out[14];
	a_r15 = A_out[15];
}
void NTTSPMB::Radix16_BU_INTT(ZZ &a_r0,ZZ &a_r1,ZZ &a_r2,ZZ &a_r3,ZZ &a_r4,ZZ &a_r5,ZZ &a_r6,ZZ &a_r7,
ZZ &a_r8,ZZ &a_r9,ZZ &a_r10,ZZ &a_r11,ZZ &a_r12,ZZ &a_r13,ZZ &a_r14,ZZ &a_r15){
	ZZ twiddle_1_16; // W^(N/8)
	ZZ twiddle_r16_nk;
	ZZ data_tmp;
	unsigned long len_1_16; // N/16
	
	len_1_16 = N/16;
	PowerMod(twiddle_1_16,W,len_1_16,p);
	std::vector<ZZ> A_in;
	std::vector<ZZ> A_out;
		
	A_in.resize(16);
	A_out.resize(16);
	//data input 
	A_in[0]  = a_r0;
	A_in[1]  = a_r1;
	A_in[2]  = a_r2;
	A_in[3]  = a_r3;
	A_in[4]  = a_r4;
	A_in[5]  = a_r5;
	A_in[6]  = a_r6;
	A_in[7]  = a_r7;
	A_in[8]  = a_r8;
	A_in[9]  = a_r9;
	A_in[10] = a_r10;
	A_in[11] = a_r11;
	A_in[12] = a_r12;
	A_in[13] = a_r13;
	A_in[14] = a_r14;
	A_in[15] = a_r15;
	
	long  nk_exp;
	nk_exp = 0;
	
	for(int i=0;i<16;i++){
		A_out[i] = 0;
	}
	
	for(long i=0;i<16;i++){
		for(long j=0;j<16;j++){
		    nk_exp = i * j;
			nk_exp = nk_exp % 16;
			PowerMod(twiddle_r16_nk,twiddle_1_16,nk_exp, p);
			MulMod(data_tmp,A_in[j],twiddle_r16_nk,p);
			AddMod(A_out[i],A_out[i],data_tmp,p);
		}
	}
	
	//data output
	a_r0 = A_out[0];
	a_r1 = A_out[15];
	a_r2 = A_out[14];
	a_r3 = A_out[13];
	a_r4 = A_out[12];
	a_r5 = A_out[11];
	a_r6 = A_out[10];
	a_r7 = A_out[9];
	a_r8 = A_out[8];
	a_r9 = A_out[7];
	a_r10 = A_out[6];
	a_r11 = A_out[5];
	a_r12 = A_out[4];
	a_r13 = A_out[3];
	a_r14 = A_out[2];
	a_r15 = A_out[1];
}
void NTTSPMB::NTT_radix2(std::vector<ZZ> &A){
	int            Stage;//stage	
	int            word_size;
	int            offset;
	int            BC;
	int            BC_tmp;
    int            group;
    int            ma_tmp;
	int            bn_tmp;
    int            bit_tmp;
	int            bn0_bc_tmp;
	int            bn1_bc_tmp;
	int            bn0_ma_reg;
	int            bn1_ma_reg;
	int            gray_i;
    int            BC_WIDTH;	
    std::vector<int> bit_array_tmp;

	std::ofstream DATARECORD("./NTT_R2_SPMB.txt");
	std::ofstream spmb_radix2("./SPMB_tw/spmb_radix2.txt");
	std::ofstream DTFAG_golden_st0("./SPMB_tw/DTFAG_golden_st0_16.txt");
	std::ofstream DTFAG_golden_st1("./SPMB_tw/DTFAG_golden_st1_16.txt");
	std::ofstream DTFAG_golden_st2("./SPMB_tw/DTFAG_golden_st2_16.txt");
	//-------------Bit Reverse--------------
	BitOperate RR, IntToVec, VecToInt;
	long long bit_width = log2(N);
	long long bit_width_s = log(radix)/log(2);
	vector<long long> bit_array;
	//----------------------------------------

	
    Stage = (unsigned long)ceil(log2(N));
	BC_WIDTH  = (int)ceil(log2(N/2));
	offset    =  (int)N /radix;
	word_size =  (int)N / (2 * radix);
	group     =  (int)N / (radix * radix);
	bit_array_tmp.resize(BC_WIDTH);
	DATARECORD <<"Stage: "<< Stage<<"\n";

	//std::cout << "-------------------------------------\n";
	//SRAM
	ZZ            data_tmp;
	std::vector<ZZ>  A_B0R0;
	std::vector<ZZ>  A_B0R1;
	std::vector<ZZ>  A_B1R0;
	std::vector<ZZ>  A_B1R1;
	
	A_B0R0.resize(word_size);
	A_B0R1.resize(word_size);
	A_B1R0.resize(word_size);
	A_B1R1.resize(word_size);
	int length;
	ZZ  factor;   //base factor
	ZZ  factor_t; //acctually mul factor
    
	//init load data
    for(int i = 0; i < group; i++){
		for(int j = 0 ; j < radix ; j++){
			bn_tmp = 0;
			ma_tmp = 0;
			BC = j * group + i;
            //bit calculate
            for(int j=0; j < BC_WIDTH;j++){
                bit_tmp = BC % 2;
                BC = BC >> 1;
                bit_array_tmp[j] = bit_tmp;
            } 
			BC = j * group + i;
            for(int rs = 0; rs < BC_WIDTH; rs++){
                if((bit_array_tmp[rs] == 1) && (rs != 0)) ma_tmp = ma_tmp + exp2((rs-1)); 
                bn_tmp = bn_tmp ^ bit_array_tmp[rs];
            }
			//std::cout <<"**************\n";
			//std::cout << "  BC: " <<BC<<"\n";
			if(bn_tmp == 0){
				A_B0R0[ma_tmp] = A[BC];
				A_B0R1[ma_tmp] = A[BC + offset];
				//std::cout <<"A_B0R0["<<ma_tmp<<"]"<<A[BC]<<"\n";
				//std::cout <<"A_B0R1["<<ma_tmp<<"]"<<A[BC+ offset]<<"\n";
			}else {
				A_B1R0[ma_tmp] = A[BC];
				A_B1R1[ma_tmp] = A[BC + offset];
				//std::cout <<"A_B1R0["<<ma_tmp<<"]"<<A[BC]<<"\n";
				//std::cout <<"A_B1R1["<<ma_tmp<<"]"<<A[BC+ offset]<<"\n";
			}
		}
	}
	
	ma_tmp = 0;
	bn_tmp = 0;
	BC     = 0;
	int tw_degree = 1; // siang
	int difference = 2;
	std::cout << "init load over! \n";
	//need modify to mult by twiddle factor
	for(int s = 0; s < Stage; s++){
		if(s == 0)factor = W;
		else {
			SqrMod(factor,factor,p);
			tw_degree = tw_degree * 2;
			difference = difference * 2;
		}
		DATARECORD <<"Stage: "<< s<<"\n";
		DATARECORD << "factor "<< factor <<"\n";
		DATARECORD << "****************************\n";

		spmb_radix2 <<"Now Stage: "<< s <<"\n";
		spmb_radix2 <<"twiddle factor : "<< factor <<"\n";
		for(int i = 0 ;i < group;i++){
			spmb_radix2 << "----------------i =" << i << " ----------------" << std::endl;
			DATARECORD << "----------------i =" << i << " ----------------" << std::endl;
			for(int j = 0;j < radix;j++){
				DATARECORD << "---j =" << j << " ---" << std::endl;
				DATARECORD <<"twiddle factor : "<< factor <<"\n";
				DATARECORD <<"p : "<< p <<"\n";
	    		DATARECORD << "********\n";
				gray_i  = Gray(i,group);
				BC_tmp  = j * group + gray_i;
				DATARECORD << "BC_tmp = " << BC_tmp << ", BC = " << BC << ", s = " << s << endl;
				RR_R2(BC_tmp,s,BC);
				//std::cout << "BC: " << BC <<"\n";
				length = BC_tmp >> s;
				DATARECORD << "length: " <<  length <<"\n";
				PowerMod(factor_t,factor,length,p);
				DATARECORD << "factor_t: "<<factor_t<<"\n";

				AGU_R2(BC,bn_tmp,ma_tmp);

				//-----------compute data idx-------------
				spmb_radix2 /*<< "BC = " << BC */<< "BC_tmp = " << BC_tmp << std::endl;
				
				IntToVec.IntToVec(BC_tmp, N, bit_array);
				rotate(bit_array.begin(), bit_array.begin()+ bit_width_s*s , bit_array.end());
				long long Data = VecToInt.VecToInt(bit_array, N);
				spmb_radix2 << "Data_index = ";
                spmb_radix2 << "( " ;					
				for(int k = 0; k < radix ; k++ ){
					spmb_radix2 << Data + k*(1<<(bit_width-bit_width_s-bit_width_s*s)) <<" ";	
				}
				spmb_radix2 << ") " ;
				spmb_radix2 << ", (w^" << 0 << ", w^" << tw_degree * length << ")" <<std::endl;
				//-----------------------------------------
				
				if(bn_tmp == 0){
					DATARECORD << "bn_tmp = " << bn_tmp << ", ma_tmp = " << ma_tmp << std::endl;
					bn0_bc_tmp = BC_tmp;
					DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp]<<"\n";
					DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp]<< ", factor_t = " << factor_t << ", w^" << tw_degree*length << "\n";
					Radix2_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp]);
					MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],factor_t,p);
					DATARECORD << "---after BU compute---" << std::endl;
				    DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp]<<"\n";
					DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp]<<"\n";
					bn0_ma_reg = ma_tmp;

					//----------------DTFAG golden pattern------------------
					switch(s){
						case 0:
							DTFAG_golden_st0 << 1 		  << " \n ";
							DTFAG_golden_st0 << factor_t  << " \n ";
							break;
						case 1:
							DTFAG_golden_st1 << 1 		  << " \n ";
							DTFAG_golden_st1 << factor_t  << " \n ";
							break;
						case 2:
							DTFAG_golden_st2 << 1 		  << " \n ";
							DTFAG_golden_st2 << factor_t  << " \n ";
							break;
						default:
							break;
					}
					//------------------------------------------------------
				}
			    else {
					DATARECORD << "bn_tmp = " << bn_tmp << ", ma_tmp = " << ma_tmp << std::endl;
					bn1_bc_tmp = BC_tmp;
					DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<<"\n";
					DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<< ", factor_t = " << factor_t << ", w^" << tw_degree*length << "\n";
					Radix2_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp]);
					MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],factor_t,p);
					DATARECORD << "---after BU compute---" << std::endl;
				    DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<<"\n";
					DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<<"\n";
					bn1_ma_reg = ma_tmp;
					//----------------DTFAG golden pattern------------------
					switch(s){
						case 0:
							DTFAG_golden_st0 << 1 		  << " \n ";
							DTFAG_golden_st0 << factor_t  << " \n ";
							break;
						case 1:
							DTFAG_golden_st1 << 1 		  << " \n ";
							DTFAG_golden_st1 << factor_t  << " \n ";
							break;
						case 2:
							DTFAG_golden_st2 << 1 		  << " \n ";
							DTFAG_golden_st2 << factor_t  << " \n ";
							break;
						default:
							break;
					}
					//------------------------------------------------------
				}
				DATARECORD <<"--------------------------------------------------------------------\n";
			}
		//data relocation
		if(s < Stage-1){
			if(bn1_bc_tmp > bn0_bc_tmp){
				DATARECORD << "bn0_bc_tmp = " << bn0_bc_tmp << ", bn1_bc_tmp = " << bn1_bc_tmp << std::endl;
				data_tmp = A_B0R1[bn0_ma_reg];
				A_B0R1[bn0_ma_reg] = A_B1R0[bn1_ma_reg];
				A_B1R0[bn1_ma_reg] = data_tmp;
			}else {
				DATARECORD << "bn0_bc_tmp = " << bn0_bc_tmp << ", bn1_bc_tmp = " << bn1_bc_tmp << std::endl;
				data_tmp = A_B1R1[bn1_ma_reg];
				A_B1R1[bn1_ma_reg] = A_B0R0[bn0_ma_reg];
				A_B0R0[bn0_ma_reg] = data_tmp;
			}
		}
		//std::cout <<"A_B0R0["<<bn0_ma_reg<<"]: "<<A_B0R0[bn0_ma_reg]<<"\n";
		//std::cout <<"A_B0R1["<<bn0_ma_reg<<"]: "<<A_B0R1[bn0_ma_reg]<<"\n";
		//std::cout <<"A_B1R0["<<bn1_ma_reg<<"]: "<<A_B1R0[bn1_ma_reg]<<"\n";
		//std::cout <<"A_B1R1["<<bn1_ma_reg<<"]: "<<A_B1R1[bn1_ma_reg]<<"\n";
        //std::cout <<"------------------------------------------\n";		 
		}
	}
	
	int index0;
    int index1;
	//data output
	//bit reverse
	for(int i = 0;i < group;i++){
		for(int j = 0;j < radix;j++){
			gray_i  = Gray(i,group);
			BC_tmp  = j * group + gray_i;
			RR_R2(BC_tmp,Stage-1,BC);
			AGU_R2(BC,bn_tmp,ma_tmp);
			if(bn_tmp == 0){
			   //std::cout <<"BC_tmp*2: "<< BC_tmp*2<<"\n";
			   //std::cout <<"Bit reverse: "<< index_r0<<"\n";
			   BR_R2(2 * BC_tmp,index0);
			   BR_R2(2 * BC_tmp + 1,index1);
               A[index0]     = A_B0R0[ma_tmp];
			   A[index1] = A_B0R1[ma_tmp];
			}
			else {
			   //std::cout <<"BC_tmp*2: "<< BC_tmp*2<<"\n";
			   //std::cout <<"Bit reverse: "<< index_r0<<"\n";
			   BR_R2(2 * BC_tmp,index0);
			   BR_R2(2 * BC_tmp + 1,index1);
               A[index0]     = A_B1R0[ma_tmp];
               A[index1] = A_B1R1[ma_tmp];
			}			
		}		
	}
}
void NTTSPMB::NTT_radix4(std::vector<ZZ> &A){
	int            Stage;//stage	
	int            word_size;
	int            offset;
	int            BC;
	int            BC_tmp;
    int            group;
    int            ma_tmp;
	int            bn_tmp;
    int            bit_tmp;
	int            bn0_bc_tmp;//frist in bc data
	int            bn1_bc_tmp;//frist in bc data
	int            bn0_ma_reg1;
	int            bn0_ma_reg2;
	int            bn1_ma_reg1;
	int            bn1_ma_reg2;
	int            gray_i;
    int            BC_WIDTH;	
    std::vector<int> bit_array_tmp;
	
	std::ofstream DATARECORD("./NTT_R4_SPMB.txt");
	std::ofstream spmb_radix4("./SPMB_tw/spmb_radix4.txt");
	std::ofstream DTFAG_golden_st0("./SPMB_tw/DTFAG_golden_st0_256.txt");
	std::ofstream DTFAG_golden_st1("./SPMB_tw/DTFAG_golden_st1_256.txt");
	std::ofstream DTFAG_golden_st2("./SPMB_tw/DTFAG_golden_st2_256.txt");
	//-------------Bit Reverse--------------
	BitOperate RR, IntToVec, VecToInt, BR;
	long long bit_width = log2(N);
	long long bit_width_s = log(radix)/log(2);
	vector<long long> bit_array;
	//----------------------------------------
	
    Stage = (unsigned long)ceil(log2(N));
	Stage = (unsigned long)Stage / 2;
	BC_WIDTH  = (int)ceil(log2(N/4));
	offset    =  (int)N /radix;
	word_size =  (int)N / (2 * radix);
	group     =  (int)N / (radix * radix);
	bit_array_tmp.resize(BC_WIDTH);
	DATARECORD <<"Stage: "<< Stage<<"\n";
	//std::cout <<"BC_WIDTH: "<< BC_WIDTH<<"\n";
	//std::cout <<"word_size: "<< word_size<<"\n";
	//std::cout << "-------------------------------------\n";
	//SRAM
	ZZ               data_tmp;
	ZZ               data_tmp_1;
	ZZ               data_tmp_2;
	std::vector<ZZ>  A_B0R0;
	std::vector<ZZ>  A_B0R1;
	std::vector<ZZ>  A_B0R2;
	std::vector<ZZ>  A_B0R3;
	std::vector<ZZ>  A_B1R0;
	std::vector<ZZ>  A_B1R1;
	std::vector<ZZ>  A_B1R2;
	std::vector<ZZ>  A_B1R3;
	
	A_B0R0.resize(word_size);
	A_B0R1.resize(word_size);
	A_B0R2.resize(word_size);
	A_B0R3.resize(word_size);
	A_B1R0.resize(word_size);
	A_B1R1.resize(word_size);
	A_B1R2.resize(word_size);
	A_B1R3.resize(word_size);
	int length;
	ZZ  factor;   //base factor
	ZZ  factor_t; //acctually mul factor
	ZZ  factor_2t;
	ZZ  factor_3t;
    
	//init load data
    for(int i = 0; i < group; i++){
		for(int j = 0 ; j < radix ; j++){
			bn_tmp = 0;
			ma_tmp = 0;
			gray_i  = Gray(i,group);
			BC = j * group + gray_i;
            //bit calculate
            for(int j=0; j < BC_WIDTH;j++){
                bit_tmp = BC % 2;
                BC = BC >> 1;
                bit_array_tmp[j] = bit_tmp;
            } 
			gray_i  = Gray(i,group);
			BC = j * group + gray_i;
            for(int rs = 0; rs < BC_WIDTH; rs++){
                if((bit_array_tmp[rs] == 1) && (rs != 0)) ma_tmp = ma_tmp + exp2((rs-1)); 
                bn_tmp = bn_tmp ^ bit_array_tmp[rs];
            }
			//std::cout <<"**************\n";
			//std::cout << "  BC: " <<BC<<"\n";
			if(bn_tmp == 0){
				A_B0R0[ma_tmp] = A[BC];
				A_B0R1[ma_tmp] = A[BC + offset];
				A_B0R2[ma_tmp] = A[BC + 2 * offset];
				A_B0R3[ma_tmp] = A[BC + 3 * offset];
				//std::cout <<"A_B0R0["<<ma_tmp<<"]"<<A[BC]<<"\n";
				//std::cout <<"A_B0R1["<<ma_tmp<<"]"<<A[BC+ offset]<<"\n";
				//std::cout <<"A_B0R2["<<ma_tmp<<"]"<<A[BC+ 2*offset]<<"\n";
				//std::cout <<"A_B0R3["<<ma_tmp<<"]"<<A[BC+ 3*offset]<<"\n";
			}else {
				A_B1R0[ma_tmp] = A[BC];
				A_B1R1[ma_tmp] = A[BC + offset];
				A_B1R2[ma_tmp] = A[BC + 2 * offset];
				A_B1R3[ma_tmp] = A[BC + 3 * offset];
				//std::cout <<"A_B1R0["<<ma_tmp<<"]"<<A[BC]<<"\n";
				//std::cout <<"A_B1R1["<<ma_tmp<<"]"<<A[BC+ offset]<<"\n";
				//std::cout <<"A_B1R2["<<ma_tmp<<"]"<<A[BC+ 2*offset]<<"\n";
				//std::cout <<"A_B1R3["<<ma_tmp<<"]"<<A[BC+ 3*offset]<<"\n";
			}
		}
	}
	std::cout <<"-----------------------------------------\n";
	ma_tmp = 0;
	bn_tmp = 0;
	BC     = 0;
	int tw_degree = 1; // siang
	std::cout << "init load over! \n";
	//need modify to mult by twiddle factor
	for(int s = 0; s < Stage; s++){
		if(s == 0){
			factor = W;
			std::cout << "W = " << W << std::endl;
		}
		else {
			SqrMod(factor,factor,p);
			SqrMod(factor,factor,p);
			tw_degree = tw_degree * 4;
		}
		DATARECORD <<"---------------------------------\n";
		DATARECORD <<"Now Stage: "<< s <<"\n";
		

		spmb_radix4 <<"Now Stage: "<< s <<"\n";
		spmb_radix4 <<"twiddle factor : "<< factor <<"\n";
	    spmb_radix4 << "********\n";
		for(int i = 0 ;i < group;i++){
			DATARECORD <<"twiddle factor : "<< factor <<"\n";
			DATARECORD <<"p : "<< p <<"\n";
	    	DATARECORD << "********\n";
			bn0_bc_tmp  = 0;
			bn1_bc_tmp  = 0;
			spmb_radix4 <<"--------------i = " << i << "----------------\n";
			for(int j = 0;j < radix;j++){
				gray_i  = Gray(i,group);
				BC_tmp  = j * group + gray_i;
				DATARECORD << "i: " << i <<"\n";
				DATARECORD << "gray_i: " << gray_i <<"\n";
				DATARECORD << "BC_tmp: " << BC_tmp <<"\n";
				RR_R4(BC_tmp,s,BC);
				DATARECORD << "BC: " << BC <<"\n";
				length = BC_tmp >> (2*s);
				DATARECORD << "length: " <<  length <<"\n";

				//-----------compute data idx-------------
				spmb_radix4 << "BC_tmp = " << BC_tmp << std::endl;
				
				IntToVec.IntToVec(BC_tmp, N, bit_array);
				rotate(bit_array.begin(), bit_array.begin()+ bit_width_s*s , bit_array.end());
				long long Data = VecToInt.VecToInt(bit_array, N);
				spmb_radix4 << "Data_index = ";
                spmb_radix4 << "( " ;		
				int spmb_radix4_arr[radix];		
				for(long long k = 0; k < radix ; k++ ){
					long long BR_tmp = BR.BitReserve(k, log2(radix));
					spmb_radix4_arr[BR_tmp] = Data + k*(1<<(bit_width-bit_width_s-bit_width_s*s));	
				}
				for(int k=0; k<radix; k++) {
					spmb_radix4 << spmb_radix4_arr[k] << ", ";
				}
				spmb_radix4 << ") " ;
				spmb_radix4 << ", (w^" << 0 << ", w^" << tw_degree * length 
				<< ", w^"  << tw_degree * length * 2 << ", w^" << tw_degree * length * 3 << ")" << std::endl;
				//-----------------------------------------
				
				PowerMod(factor_t,factor,length,p);
				DATARECORD << "factor_t: "<<factor_t<<"\n";
				AGU_R4(BC,bn_tmp,ma_tmp);
				//spmb_radix4 << "bn_tmp = " << bn_tmp << ", ma_tmp = " << ma_tmp << std::endl;
				DATARECORD << "bn_tmp: "<<bn_tmp<<"\n";
				int debug = 0;
				if(bn_tmp == 0){
					if(j < 2)bn0_bc_tmp = BC_tmp;
					PowerMod(factor_2t,factor_t,2,p);
					PowerMod(factor_3t,factor_t,3,p);
					if(debug) DATARECORD <<"Before butterfly unit operation! \n";
					if(debug) DATARECORD <<" A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp] << std::endl;
				    if(debug) DATARECORD <<" A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp] << std::endl;
					if(debug) DATARECORD <<" A_B0R2["<<ma_tmp<<"]: "<<A_B0R2[ma_tmp] << std::endl;
					if(debug) DATARECORD <<" A_B0R3["<<ma_tmp<<"]: "<<A_B0R3[ma_tmp] << std::endl;
					if(!debug) Radix4_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp],A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
					DATARECORD <<" -------------------" << std::endl;
					if(debug) DATARECORD <<" A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp] << ", factor_0 : w^" << tw_degree*0<<"\n";
					if(debug) DATARECORD <<" A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp] << ", facotr_1 : w^" << tw_degree*length  <<"\n";
					if(debug) DATARECORD <<" A_B0R2["<<ma_tmp<<"]: "<<A_B0R2[ma_tmp] << ", facotr_2 : w^" << tw_degree*length*2 <<"\n";
					if(debug) DATARECORD <<" A_B0R3["<<ma_tmp<<"]: "<<A_B0R3[ma_tmp] << ", facotr_3 : w^" << tw_degree*length*3 <<"\n";
					if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],factor_t,p);
					if(!debug) MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],factor_2t,p);
					if(!debug) MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],factor_3t,p);
					DATARECORD <<"facotr_1 : " << factor_t  << ", length = " << length * 1 << " \n";
					DATARECORD <<"facotr_2 : " << factor_2t << ", length = " << length * 2 << " \n";
					DATARECORD <<"facotr_3 : " << factor_3t << ", length = " << length * 3 << " \n";
					
					//----------------DTFAG golden pattern------------------
					switch(s){
						case 0:
							DTFAG_golden_st0 << 1 		  << " \n ";
							DTFAG_golden_st0 << factor_t  << " \n ";
							DTFAG_golden_st0 << factor_2t << " \n ";
							DTFAG_golden_st0 << factor_3t << " \n ";
							break;
						case 1:
							DTFAG_golden_st1 << 1 		  << " \n ";
							DTFAG_golden_st1 << factor_t  << " \n ";
							DTFAG_golden_st1 << factor_2t << " \n ";
							DTFAG_golden_st1 << factor_3t << " \n ";
							break;
						case 2:
							DTFAG_golden_st2 << 1 		  << " \n ";
							DTFAG_golden_st2 << factor_t  << " \n ";
							DTFAG_golden_st2 << factor_2t << " \n ";
							DTFAG_golden_st2 << factor_3t << " \n ";
							break;
						default:
							break;
					}
					//------------------------------------------------------
				    //DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp]<<"\n";
				    //DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp]<<"\n";
					//DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<<A_B0R2[ma_tmp]<<"\n";
					//DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<<A_B0R3[ma_tmp]<<"\n";
					DATARECORD <<"--------------------------------------------------------------------\n";
					if(j <  2)bn0_ma_reg1 = ma_tmp;
					if(j >= 2)bn0_ma_reg2 = ma_tmp;
				}
			    else {
					if(j < 2)bn1_bc_tmp = BC_tmp;
					PowerMod(factor_2t,factor_t,2,p);
					PowerMod(factor_3t,factor_t,3,p);
					if(debug) DATARECORD <<"Before butterfly unit operation! \n";
					if(debug) DATARECORD <<" A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp] << std::endl;
					if(debug) DATARECORD <<" A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp] << std::endl;
					if(debug) DATARECORD <<" A_B1R2["<<ma_tmp<<"]: "<<A_B1R2[ma_tmp] << std::endl;
					if(debug) DATARECORD <<" A_B1R3["<<ma_tmp<<"]: "<<A_B1R3[ma_tmp] << std::endl;
					if(!debug) Radix4_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp],A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
					DATARECORD <<" -------------------" << std::endl;
					if(debug) DATARECORD <<" A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<< ", factor_0 : w^" << tw_degree*0<<"\n";
					if(debug) DATARECORD <<" A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<< ", facotr_1 : w^" << tw_degree*length  <<"\n";
					if(debug) DATARECORD <<" A_B1R2["<<ma_tmp<<"]: "<<A_B1R2[ma_tmp]<< ", facotr_2 : w^" << tw_degree*length*2 <<"\n";
					if(debug) DATARECORD <<" A_B1R3["<<ma_tmp<<"]: "<<A_B1R3[ma_tmp]<< ", facotr_3 : w^" << tw_degree*length*3 <<"\n";
					if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],factor_t,p);
					if(!debug) MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],factor_2t,p);
					if(!debug) MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],factor_3t,p);
					DATARECORD <<"facotr_1 : " << factor_t  << ", length = " << length * 1 << " \n";
					DATARECORD <<"facotr_2 : " << factor_2t << ", length = " << length * 2 << " \n";
					DATARECORD <<"facotr_3 : " << factor_3t << ", length = " << length * 3 << " \n";
					//----------------DTFAG golden pattern------------------
					switch(s){
						case 0:
							DTFAG_golden_st0 << 1 		  << " \n ";
							DTFAG_golden_st0 << factor_t  << " \n ";
							DTFAG_golden_st0 << factor_2t << " \n ";
							DTFAG_golden_st0 << factor_3t << " \n ";
							break;
						case 1:
							DTFAG_golden_st1 << 1 		  << " \n ";
							DTFAG_golden_st1 << factor_t  << " \n ";
							DTFAG_golden_st1 << factor_2t << " \n ";
							DTFAG_golden_st1 << factor_3t << " \n ";
							break;
						case 2:
							DTFAG_golden_st2 << 1 		  << " \n ";
							DTFAG_golden_st2 << factor_t  << " \n ";
							DTFAG_golden_st2 << factor_2t << " \n ";
							DTFAG_golden_st2 << factor_3t << " \n ";
							break;
						default:
							break;
					}
					//------------------------------------------------------	
				    //DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp] <<"\n";
				    //DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp] <<"\n";
					//DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<<A_B1R2[ma_tmp] <<"\n";
					//DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<<A_B1R3[ma_tmp] <<"\n";
					DATARECORD <<"--------------------------------------------------------------------\n";
					if(j < 2) bn1_ma_reg1 = ma_tmp;
					if(j >= 2)bn1_ma_reg2 = ma_tmp;
				}
			}
		//data relocation
		 if(s < Stage-1){
		  if(bn1_bc_tmp > bn0_bc_tmp){
             if(bn0_ma_reg1 > bn0_ma_reg2){
				 ma_tmp = bn0_ma_reg1;
				 bn0_ma_reg1 = bn0_ma_reg2;
				 bn0_ma_reg2 = ma_tmp;
			 }
			 if(bn1_ma_reg1 > bn1_ma_reg2){
				 ma_tmp = bn1_ma_reg1;
				 bn1_ma_reg1 = bn1_ma_reg2;
				 bn1_ma_reg2 = ma_tmp;
			 } 
			 data_tmp   = A_B0R1[bn0_ma_reg1];
			 data_tmp_1 = A_B0R2[bn0_ma_reg1];
			 data_tmp_2 = A_B0R3[bn0_ma_reg1];
			 A_B0R1[bn0_ma_reg1] = A_B1R0[bn1_ma_reg1];
			 A_B0R2[bn0_ma_reg1] = A_B1R0[bn1_ma_reg2];
			 A_B0R3[bn0_ma_reg1] = A_B0R0[bn0_ma_reg2];
			 A_B1R0[bn1_ma_reg1] = data_tmp;
			 A_B1R0[bn1_ma_reg2] = data_tmp_1;
			 A_B0R0[bn0_ma_reg2] = data_tmp_2;
			 data_tmp   = A_B1R2[bn1_ma_reg1];
			 data_tmp_1 = A_B1R3[bn1_ma_reg1];
			 A_B1R2[bn1_ma_reg1] = A_B1R1[bn1_ma_reg2];
			 A_B1R3[bn1_ma_reg1] = A_B0R1[bn0_ma_reg2];
			 A_B1R1[bn1_ma_reg2] = data_tmp;
			 A_B0R1[bn0_ma_reg2] = data_tmp_1;
			 data_tmp = A_B1R3[bn1_ma_reg2];
			 A_B1R3[bn1_ma_reg2] = A_B0R2[bn0_ma_reg2];
			 A_B0R2[bn0_ma_reg2] = data_tmp;
		  }else {
             if(bn0_ma_reg1 > bn0_ma_reg2){
				 ma_tmp = bn0_ma_reg1;
				 bn0_ma_reg1 = bn0_ma_reg2;
				 bn0_ma_reg2 = ma_tmp;
			 }
			 if(bn1_ma_reg1 > bn1_ma_reg2){
				 ma_tmp = bn1_ma_reg1;
				 bn1_ma_reg1 = bn1_ma_reg2;
				 bn1_ma_reg2 = ma_tmp;
			 } 	 	 
		 	 data_tmp   = A_B1R1[bn1_ma_reg1];
		 	 data_tmp_1 = A_B1R2[bn1_ma_reg1];
		 	 data_tmp_2 = A_B1R3[bn1_ma_reg1];
			 A_B1R1[bn1_ma_reg1] = A_B0R0[bn0_ma_reg1];
			 A_B1R2[bn1_ma_reg1] = A_B0R0[bn0_ma_reg2];
			 A_B1R3[bn1_ma_reg1] = A_B1R0[bn1_ma_reg2];
			 A_B0R0[bn0_ma_reg1] = data_tmp;
			 A_B0R0[bn0_ma_reg2] = data_tmp_1;
			 A_B1R0[bn1_ma_reg2] = data_tmp_2;
			 data_tmp   = A_B0R2[bn0_ma_reg1];
			 data_tmp_1 = A_B0R3[bn0_ma_reg1];
			 A_B0R2[bn0_ma_reg1]  = A_B0R1[bn0_ma_reg2];
			 A_B0R3[bn0_ma_reg1]  = A_B1R1[bn1_ma_reg2];
			 A_B0R1[bn0_ma_reg2]  = data_tmp;
		 	 A_B1R1[bn1_ma_reg2]  = data_tmp_1;
			 data_tmp = A_B0R3[bn0_ma_reg2];
             A_B0R3[bn0_ma_reg2] = A_B1R2[bn1_ma_reg2];
			 A_B1R2[bn1_ma_reg2] = data_tmp;
			 
		  }
		 }	 
		}
	}
	
	int index0;
    int index1;
    int index2;
    int index3;
	//data output
	//bit reverse
	for(int i = 0;i < group;i++){
		for(int j = 0;j < radix;j++){
			gray_i  = Gray(i,group);
			BC_tmp  = j * group + gray_i;
			RR_R4(BC_tmp,Stage-1,BC);
			AGU_R4(BC,bn_tmp,ma_tmp);
			if(bn_tmp == 0){
			   //std::cout <<"BC_tmp*2: "<< BC_tmp*2<<"\n";
			   //std::cout <<"Bit reverse: "<< index_r0<<"\n";
			   BR_R4(4 * BC_tmp,index0);
			   BR_R4(4 * BC_tmp + 1,index1);
			   BR_R4(4 * BC_tmp + 2,index2);
			   BR_R4(4 * BC_tmp + 3,index3);
               A[index0] = A_B0R0[ma_tmp];
			   A[index1] = A_B0R1[ma_tmp];
			   A[index2] = A_B0R2[ma_tmp];
			   A[index3] = A_B0R3[ma_tmp];
			}
			else {
			   //std::cout <<"BC_tmp*2: "<< BC_tmp*2<<"\n";
			   //std::cout <<"Bit reverse: "<< index_r0<<"\n";
			   BR_R4(4 * BC_tmp,index0);
			   BR_R4(4 * BC_tmp + 1,index1);
			   BR_R4(4 * BC_tmp + 2,index2);
			   BR_R4(4 * BC_tmp + 3,index3);
               A[index0]     = A_B1R0[ma_tmp];
               A[index1]     = A_B1R1[ma_tmp];
               A[index2]     = A_B1R2[ma_tmp];
               A[index3]     = A_B1R3[ma_tmp];
			}			
		}		
	}
}
void NTTSPMB::NTT_radix16(std::vector<ZZ> &A){
	int            Stage;//stage	
	int            word_size;
	int            offset;
	int            BC;
	int            BC_tmp;
    int            group;
    int            ma_tmp;
	int            bn_tmp;
    int            bit_tmp;
	int            bn0_bc_tmp;//frist in bc data
	int            bn1_bc_tmp;//frist in bc data
	int            bn0_ma_reg1;
	int            bn0_ma_reg2;
	int            bn0_ma_reg3;
	int            bn0_ma_reg4;
	int            bn0_ma_reg5;
	int            bn0_ma_reg6;
	int            bn0_ma_reg7;
	int            bn0_ma_reg8;
	int            bn1_ma_reg1;
	int            bn1_ma_reg2;
	int            bn1_ma_reg3;
	int            bn1_ma_reg4;
	int            bn1_ma_reg5;
	int            bn1_ma_reg6;
	int            bn1_ma_reg7;
	int            bn1_ma_reg8;
	int            gray_i;
    int            BC_WIDTH;	
    std::vector<int> bit_array_tmp;
	
	std::ofstream DATARECORD("./NTT_R16_SPMB.txt");
	std::ofstream siang_record("./my_print_data/single_radix16_SPMB.txt");
	std::ofstream TF_record("./my_print_data/single_radix16_SPMB_TF_dg.txt");
	std::ofstream DTFAG_golden_st0("./SPMB_tw/DTFAG_golden_st0_65536.txt");
	std::ofstream DTFAG_golden_st1("./SPMB_tw/DTFAG_golden_st1_65536.txt");
	std::ofstream DTFAG_golden_st2("./SPMB_tw/DTFAG_golden_st2_65536.txt");

    Stage = (unsigned long)ceil(log2(N));
	Stage = (unsigned long)Stage / 4;
	BC_WIDTH  = (int)ceil(log2(N/16));
	offset    =  (int)N /radix;
	word_size =  (int)N / (2 * radix);
	group     =  (int)N / (radix * radix);
	bit_array_tmp.resize(BC_WIDTH);
	DATARECORD <<"Stage: "<< Stage<<"\n";
	siang_record <<"Stage: "<< Stage<<"\n";
	//std::cout <<"BC_WIDTH: "<< BC_WIDTH<<"\n";
	//std::cout <<"word_size: "<< word_size<<"\n";
	//std::cout << "-------------------------------------\n";
	//SRAM

	ZZ               data_tmp_1;
	ZZ               data_tmp_2;
	ZZ               data_tmp_3;
	ZZ               data_tmp_4;
	ZZ               data_tmp_5;
	ZZ               data_tmp_6;
	ZZ               data_tmp_7;
	ZZ               data_tmp_8;
	ZZ               data_tmp_9;
	ZZ               data_tmp_10;
	ZZ               data_tmp_11;
	ZZ               data_tmp_12;
	ZZ               data_tmp_13;
	ZZ               data_tmp_14;
	ZZ               data_tmp_15;
	std::vector<ZZ>  A_B0R0;
	std::vector<ZZ>  A_B0R1;
	std::vector<ZZ>  A_B0R2;
	std::vector<ZZ>  A_B0R3;
	std::vector<ZZ>  A_B0R4;
	std::vector<ZZ>  A_B0R5;
	std::vector<ZZ>  A_B0R6;
	std::vector<ZZ>  A_B0R7;
	std::vector<ZZ>  A_B0R8;
	std::vector<ZZ>  A_B0R9;
	std::vector<ZZ>  A_B0R10;
	std::vector<ZZ>  A_B0R11;
	std::vector<ZZ>  A_B0R12;
	std::vector<ZZ>  A_B0R13;
	std::vector<ZZ>  A_B0R14;
	std::vector<ZZ>  A_B0R15;
	std::vector<ZZ>  A_B1R0;
	std::vector<ZZ>  A_B1R1;
	std::vector<ZZ>  A_B1R2;
	std::vector<ZZ>  A_B1R3;
	std::vector<ZZ>  A_B1R4;
	std::vector<ZZ>  A_B1R5;
	std::vector<ZZ>  A_B1R6;
	std::vector<ZZ>  A_B1R7;
	std::vector<ZZ>  A_B1R8;
	std::vector<ZZ>  A_B1R9;
	std::vector<ZZ>  A_B1R10;
	std::vector<ZZ>  A_B1R11;
	std::vector<ZZ>  A_B1R12;
	std::vector<ZZ>  A_B1R13;
	std::vector<ZZ>  A_B1R14;
	std::vector<ZZ>  A_B1R15;
	
	A_B0R0.resize(word_size);
	A_B0R1.resize(word_size);
	A_B0R2.resize(word_size);
	A_B0R3.resize(word_size);
	A_B0R4.resize(word_size);
	A_B0R5.resize(word_size);
	A_B0R6.resize(word_size);
	A_B0R7.resize(word_size);
	A_B0R8.resize(word_size);
	A_B0R9.resize(word_size);
	A_B0R10.resize(word_size);
	A_B0R11.resize(word_size);
	A_B0R12.resize(word_size);
	A_B0R13.resize(word_size);
	A_B0R14.resize(word_size);
	A_B0R15.resize(word_size);
	A_B1R0.resize(word_size);
	A_B1R1.resize(word_size);
	A_B1R2.resize(word_size);
	A_B1R3.resize(word_size);
	A_B1R4.resize(word_size);
	A_B1R5.resize(word_size);
	A_B1R6.resize(word_size);
	A_B1R7.resize(word_size);
	A_B1R8.resize(word_size);
	A_B1R9.resize(word_size);
	A_B1R10.resize(word_size);
	A_B1R11.resize(word_size);
	A_B1R12.resize(word_size);
	A_B1R13.resize(word_size);
	A_B1R14.resize(word_size);
	A_B1R15.resize(word_size);
	int length;
	ZZ  factor;   //base factor
	ZZ  factor_t; //acctually mul factor
	ZZ  factor_2t;
	ZZ  factor_3t;
	ZZ  factor_4t;
	ZZ  factor_5t;
	ZZ  factor_6t;
	ZZ  factor_7t;
	ZZ  factor_8t;
	ZZ  factor_9t;
	ZZ  factor_10t;
	ZZ  factor_11t;
	ZZ  factor_12t;
	ZZ  factor_13t;
	ZZ  factor_14t;
	ZZ  factor_15t;
    
	//init load data
    for(int i = 0; i < group; i++){
		for(int j = 0 ; j < radix ; j++){
			bn_tmp = 0;
			ma_tmp = 0;
			gray_i  = Gray(i,group);
			BC = j * group + gray_i;
            //bit calculate
            for(int j=0; j < BC_WIDTH;j++){
                bit_tmp = BC % 2;
                BC = BC >> 1;
                bit_array_tmp[j] = bit_tmp;
            } 
			gray_i  = Gray(i,group);
			BC = j * group + gray_i;
            for(int rs = 0; rs < BC_WIDTH; rs++){
                if((bit_array_tmp[rs] == 1) && (rs != 0)) ma_tmp = ma_tmp + exp2((rs-1)); 
                bn_tmp = bn_tmp ^ bit_array_tmp[rs];
            }
			//std::cout <<"**************\n";
			//std::cout << "  BC: " <<BC<<"\n";
			if(bn_tmp == 0){
				A_B0R0[ma_tmp]  = A[BC];
				A_B0R1[ma_tmp]  = A[BC +      offset];
				A_B0R2[ma_tmp]  = A[BC + 2  * offset];
				A_B0R3[ma_tmp]  = A[BC + 3  * offset];
				A_B0R4[ma_tmp]  = A[BC + 4  * offset];
				A_B0R5[ma_tmp]  = A[BC + 5  * offset];
				A_B0R6[ma_tmp]  = A[BC + 6  * offset];
				A_B0R7[ma_tmp]  = A[BC + 7  * offset];
				A_B0R8[ma_tmp]  = A[BC + 8  * offset];
				A_B0R9[ma_tmp]  = A[BC + 9  * offset];
				A_B0R10[ma_tmp] = A[BC + 10 * offset];
				A_B0R11[ma_tmp] = A[BC + 11 * offset];
				A_B0R12[ma_tmp] = A[BC + 12 * offset];
				A_B0R13[ma_tmp] = A[BC + 13 * offset];
				A_B0R14[ma_tmp] = A[BC + 14 * offset];
				A_B0R15[ma_tmp] = A[BC + 15 * offset];
			}else {
				A_B1R0[ma_tmp]  = A[BC];
				A_B1R1[ma_tmp]  = A[BC +     offset];
				A_B1R2[ma_tmp]  = A[BC + 2 * offset];
				A_B1R3[ma_tmp]  = A[BC + 3 * offset];
				A_B1R4[ma_tmp]  = A[BC + 4 * offset];
				A_B1R5[ma_tmp]  = A[BC + 5 * offset];
				A_B1R6[ma_tmp]  = A[BC + 6 * offset];
				A_B1R7[ma_tmp]  = A[BC + 7 * offset];
				A_B1R8[ma_tmp]  = A[BC + 8 * offset];
				A_B1R9[ma_tmp]  = A[BC + 9 * offset];
				A_B1R10[ma_tmp] = A[BC + 10 * offset];
				A_B1R11[ma_tmp] = A[BC + 11 * offset];
				A_B1R12[ma_tmp] = A[BC + 12 * offset];
				A_B1R13[ma_tmp] = A[BC + 13 * offset];
				A_B1R14[ma_tmp] = A[BC + 14 * offset];
				A_B1R15[ma_tmp] = A[BC + 15 * offset];
			}
		}
	}
	//std::cout <<"-----------------------------------------\n";
	ma_tmp = 0;
	bn_tmp = 0;
	BC     = 0;
	int tw_degree = 1; // siang
	std::cout << "init load over! \n";
	//need modify to mult by twiddle factor
	for(int s = 0; s < Stage; s++){
		if(s == 0) {
			factor = W;
		}
		else {
			SqrMod(factor,factor,p);
			SqrMod(factor,factor,p);
			SqrMod(factor,factor,p);
			SqrMod(factor,factor,p);
			tw_degree = tw_degree * 16;
		}
		DATARECORD <<"---------------------------------\n";
		DATARECORD <<"Now Stage: "<< s <<"\n";
		DATARECORD <<"twiddle factor : "<< factor <<"\n";
		siang_record <<"---------------------------------\n"; // siang
		siang_record <<"Now Stage: "<< s <<"\n";				// siang
		siang_record <<"twiddle factor : "<< factor <<"\n";	// siang
		TF_record <<"---------------------------------\n";
		TF_record <<"Now Stage: "<< s <<"\n";				// siang
		std::cout << "twiddle factor : "<< factor <<"\n";
	    DATARECORD << "********\n";
		siang_record << "********\n";	// siang
		for(int i = 0 ;i < group;i++){
			bn0_bc_tmp  = 0;
			bn1_bc_tmp  = 0;
			TF_record << "----------------------------------" <<"\n";	// TF_record
			TF_record <<"i: "<< i << "\n";	
			for(int j = 0;j < radix;j++){
				siang_record <<"twiddle factor : "<< factor <<"\n";	// siang
				gray_i  = Gray(i,group);
				BC_tmp  = j * group + gray_i;
				DATARECORD << "i: " << i <<"\n";
				siang_record << "i: " << i << ", j: " << j << "\n";
				DATARECORD << "gray_i: " << gray_i <<"\n";
				DATARECORD << "BC_tmp: " << BC_tmp <<"\n";
				siang_record << "BC_tmp: " << BC_tmp <<"\n";	// siang
				RR_R16(BC_tmp,s,BC);
				DATARECORD << "BC: " << BC <<"\n";
				siang_record << "BC: " << BC <<"\n";	// siang
				length = BC_tmp >> (4*s);
				DATARECORD << "length: " <<  length <<"\n";
				siang_record << "length: " <<  length <<"\n";	// siang
				TF_record << "BC = " << BC << ", BC_tmp = " << BC_tmp 
				<< ", length: " << length << ", tw_dg: " <<  tw_degree * length << " = " << tw_degree << " * " << length <<"\n";	// TF_record
				
				PowerMod(factor_t,factor,length,p);
				DATARECORD << "factor_t: "<<factor_t<<"\n";
				siang_record << "factor_t: "<<factor_t<<"\n";	// siang
				AGU_R16(BC,bn_tmp,ma_tmp);
				DATARECORD << "bn_tmp: "<<bn_tmp<<"\n";
				if(bn_tmp == 0){
					if(j < 2)bn0_bc_tmp = BC_tmp;
					PowerMod(factor_2t,factor_t,2,p);
					PowerMod(factor_3t,factor_t,3,p);
					PowerMod(factor_4t,factor_t,4,p);
					PowerMod(factor_5t,factor_t,5,p);
					PowerMod(factor_6t,factor_t,6,p);
					PowerMod(factor_7t,factor_t,7,p);
					PowerMod(factor_8t,factor_t,8,p);
					PowerMod(factor_9t,factor_t,9,p);
					PowerMod(factor_10t,factor_t,10,p);
					PowerMod(factor_11t,factor_t,11,p);
					PowerMod(factor_12t,factor_t,12,p);
					PowerMod(factor_13t,factor_t,13,p);
					PowerMod(factor_14t,factor_t,14,p);
					PowerMod(factor_15t,factor_t,15,p);
					DATARECORD <<" Before butterfly unit operation! \n";
					DATARECORD <<" A_B0R0["<<ma_tmp<<"]:  "<<A_B0R0[ma_tmp]<<"\n";
				    DATARECORD <<" A_B0R1["<<ma_tmp<<"]:  "<<A_B0R1[ma_tmp]<<"\n";
					DATARECORD <<" A_B0R2["<<ma_tmp<<"]:  "<<A_B0R2[ma_tmp]<<"\n";
					DATARECORD <<" A_B0R3["<<ma_tmp<<"]:  "<<A_B0R3[ma_tmp]<<"\n";
					DATARECORD <<" A_B0R4["<<ma_tmp<<"]:  "<<A_B0R4[ma_tmp]<<"\n";
					DATARECORD <<" A_B0R5["<<ma_tmp<<"]:  "<<A_B0R5[ma_tmp]<<"\n";
					DATARECORD <<" A_B0R6["<<ma_tmp<<"]:  "<<A_B0R6[ma_tmp]<<"\n";
					DATARECORD <<" A_B0R7["<<ma_tmp<<"]:  "<<A_B0R7[ma_tmp]<<"\n";
					DATARECORD <<" A_B0R8["<<ma_tmp<<"]:  "<<A_B0R8[ma_tmp]<<"\n";
					DATARECORD <<" A_B0R9["<<ma_tmp<<"]:  "<<A_B0R9[ma_tmp]<<"\n";
					DATARECORD <<" A_B0R10["<<ma_tmp<<"]: "<<A_B0R10[ma_tmp]<<"\n";
					DATARECORD <<" A_B0R11["<<ma_tmp<<"]: "<<A_B0R11[ma_tmp]<<"\n";
					DATARECORD <<" A_B0R12["<<ma_tmp<<"]: "<<A_B0R12[ma_tmp]<<"\n";
					DATARECORD <<" A_B0R13["<<ma_tmp<<"]: "<<A_B0R13[ma_tmp]<<"\n";
					DATARECORD <<" A_B0R14["<<ma_tmp<<"]: "<<A_B0R14[ma_tmp]<<"\n";
					DATARECORD <<" A_B0R15["<<ma_tmp<<"]: "<<A_B0R15[ma_tmp]<<"\n";
					Radix16_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp], A_B0R2[ma_tmp], A_B0R3[ma_tmp],A_B0R4[ma_tmp],
							   A_B0R5[ma_tmp],A_B0R6[ma_tmp], A_B0R7[ma_tmp], A_B0R8[ma_tmp],A_B0R9[ma_tmp],
							   A_B0R10[ma_tmp],A_B0R11[ma_tmp],A_B0R12[ma_tmp],A_B0R13[ma_tmp],A_B0R14[ma_tmp],
							   A_B0R15[ma_tmp]);
					DATARECORD <<" after butterfly unit operation! \n";
					DATARECORD <<" A_B0R0["<<ma_tmp<<"]:  "<<A_B0R0[ma_tmp]<<"\n";
				    DATARECORD <<" A_B0R1["<<ma_tmp<<"]:  "<<A_B0R1[ma_tmp]<<"\n";
					DATARECORD <<" A_B0R2["<<ma_tmp<<"]:  "<<A_B0R2[ma_tmp]<<"\n";
					DATARECORD <<" A_B0R3["<<ma_tmp<<"]:  "<<A_B0R3[ma_tmp]<<"\n";
					DATARECORD <<" A_B0R4["<<ma_tmp<<"]:  "<<A_B0R4[ma_tmp]<<"\n";
					DATARECORD <<" A_B0R5["<<ma_tmp<<"]:  "<<A_B0R5[ma_tmp]<<"\n";
					DATARECORD <<" A_B0R6["<<ma_tmp<<"]:  "<<A_B0R6[ma_tmp]<<"\n";
					DATARECORD <<" A_B0R7["<<ma_tmp<<"]:  "<<A_B0R7[ma_tmp]<<"\n";
					DATARECORD <<" A_B0R8["<<ma_tmp<<"]:  "<<A_B0R8[ma_tmp]<<"\n";
					DATARECORD <<" A_B0R9["<<ma_tmp<<"]:  "<<A_B0R9[ma_tmp]<<"\n";
					DATARECORD <<" A_B0R10["<<ma_tmp<<"]: "<<A_B0R10[ma_tmp]<<"\n";
					DATARECORD <<" A_B0R11["<<ma_tmp<<"]: "<<A_B0R11[ma_tmp]<<"\n";
					DATARECORD <<" A_B0R12["<<ma_tmp<<"]: "<<A_B0R12[ma_tmp]<<"\n";
					DATARECORD <<" A_B0R13["<<ma_tmp<<"]: "<<A_B0R13[ma_tmp]<<"\n";
					DATARECORD <<" A_B0R14["<<ma_tmp<<"]: "<<A_B0R14[ma_tmp]<<"\n";
					DATARECORD <<" A_B0R15["<<ma_tmp<<"]: "<<A_B0R15[ma_tmp]<<"\n";
                    MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],factor_t,p);
					MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],factor_2t,p);
					MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],factor_3t,p);
					MulMod(A_B0R4[ma_tmp],A_B0R4[ma_tmp],factor_4t,p);
					MulMod(A_B0R5[ma_tmp],A_B0R5[ma_tmp],factor_5t,p);
					MulMod(A_B0R6[ma_tmp],A_B0R6[ma_tmp],factor_6t,p);
					MulMod(A_B0R7[ma_tmp],A_B0R7[ma_tmp],factor_7t,p);
					MulMod(A_B0R8[ma_tmp],A_B0R8[ma_tmp],factor_8t,p);
					MulMod(A_B0R9[ma_tmp],A_B0R9[ma_tmp],factor_9t,p);
					MulMod(A_B0R10[ma_tmp],A_B0R10[ma_tmp],factor_10t,p);
					MulMod(A_B0R11[ma_tmp],A_B0R11[ma_tmp],factor_11t,p);
					MulMod(A_B0R12[ma_tmp],A_B0R12[ma_tmp],factor_12t,p);
					MulMod(A_B0R13[ma_tmp],A_B0R13[ma_tmp],factor_13t,p);
					MulMod(A_B0R14[ma_tmp],A_B0R14[ma_tmp],factor_14t,p);
					MulMod(A_B0R15[ma_tmp],A_B0R15[ma_tmp],factor_15t,p);
					DATARECORD <<"facotr_1  : " << factor_t  << " \n";
					DATARECORD <<"facotr_2  : " << factor_2t << " \n";
					DATARECORD <<"facotr_3  : " << factor_3t << " \n";
					DATARECORD <<"facotr_4  : " << factor_4t << " \n";
					DATARECORD <<"facotr_5  : " << factor_5t << " \n";
					DATARECORD <<"facotr_6  : " << factor_6t << " \n";
					DATARECORD <<"facotr_7  : " << factor_7t << " \n";
					DATARECORD <<"facotr_8  : " << factor_8t << " \n";
					DATARECORD <<"facotr_9  : " << factor_9t << " \n";
					DATARECORD <<"facotr_10 : " << factor_10t << " \n";
					DATARECORD <<"facotr_11 : " << factor_11t << " \n";
					DATARECORD <<"facotr_12 : " << factor_12t << " \n";
					DATARECORD <<"facotr_13 : " << factor_13t << " \n";
					DATARECORD <<"facotr_14 : " << factor_14t << " \n";
					DATARECORD <<"facotr_15 : " << factor_15t << " \n";

					//----------------DTFAG golden pattern------------------
					switch(s){
						case 0:
							DTFAG_golden_st0 << 1 		  << " \n ";
							DTFAG_golden_st0 << factor_t  << " \n ";
							DTFAG_golden_st0 << factor_2t << " \n ";
							DTFAG_golden_st0 << factor_3t << " \n ";
							DTFAG_golden_st0 << factor_4t << " \n ";
							DTFAG_golden_st0 << factor_5t << " \n ";
							DTFAG_golden_st0 << factor_6t << " \n ";
							DTFAG_golden_st0 << factor_7t << " \n ";
							DTFAG_golden_st0 << factor_8t << " \n ";
							DTFAG_golden_st0 << factor_9t << " \n ";
							DTFAG_golden_st0 << factor_10t << " \n ";
							DTFAG_golden_st0 << factor_11t << " \n ";
							DTFAG_golden_st0 << factor_12t << " \n ";
							DTFAG_golden_st0 << factor_13t << " \n ";
							DTFAG_golden_st0 << factor_14t << " \n ";
							DTFAG_golden_st0 << factor_15t << " \n ";
							break;
						case 1:
							DTFAG_golden_st1 << 1 		  << " \n ";
							DTFAG_golden_st1 << factor_t  << " \n ";
							DTFAG_golden_st1 << factor_2t << " \n ";
							DTFAG_golden_st1 << factor_3t << " \n ";
							DTFAG_golden_st1 << factor_4t << " \n ";
							DTFAG_golden_st1 << factor_5t << " \n ";
							DTFAG_golden_st1 << factor_6t << " \n ";
							DTFAG_golden_st1 << factor_7t << " \n ";
							DTFAG_golden_st1 << factor_8t << " \n ";
							DTFAG_golden_st1 << factor_9t << " \n ";
							DTFAG_golden_st1 << factor_10t << " \n ";
							DTFAG_golden_st1 << factor_11t << " \n ";
							DTFAG_golden_st1 << factor_12t << " \n ";
							DTFAG_golden_st1 << factor_13t << " \n ";
							DTFAG_golden_st1 << factor_14t << " \n ";
							DTFAG_golden_st1 << factor_15t << " \n ";
							break;
						case 2:
							DTFAG_golden_st2 << 1 		  << " \n ";
							DTFAG_golden_st2 << factor_t  << " \n ";
							DTFAG_golden_st2 << factor_2t << " \n ";
							DTFAG_golden_st2 << factor_3t << " \n ";
							DTFAG_golden_st2 << factor_4t << " \n ";
							DTFAG_golden_st2 << factor_5t << " \n ";
							DTFAG_golden_st2 << factor_6t << " \n ";
							DTFAG_golden_st2 << factor_7t << " \n ";
							DTFAG_golden_st2 << factor_8t << " \n ";
							DTFAG_golden_st2 << factor_9t << " \n ";
							DTFAG_golden_st2 << factor_10t << " \n ";
							DTFAG_golden_st2 << factor_11t << " \n ";
							DTFAG_golden_st2 << factor_12t << " \n ";
							DTFAG_golden_st2 << factor_13t << " \n ";
							DTFAG_golden_st2 << factor_14t << " \n ";
							DTFAG_golden_st2 << factor_15t << " \n ";
							break;
						default:
							break;
					}
					//------------------------------------------------------

					//-------------------siang record----------------------
					//siang_record << "factor_t: "<< factor_t <<"\n";	// siang
					siang_record << "p: "<< p <<"\n";	// siang
					siang_record <<"facotr_1  : " << factor_t  << " \n";
					siang_record <<"facotr_2_facotr_3 : " << factor_2t << "_" << factor_3t << " \n";
					siang_record <<"facotr_4_facotr_5 : " << factor_4t << "_" << factor_5t << " \n";
					siang_record <<"facotr_6_facotr_7 : " << factor_6t << "_" << factor_7t << " \n";
					siang_record <<"facotr_8_facotr_9 : " << factor_8t << "_" << factor_9t << " \n";
					siang_record <<"facotr_10_facotr_11 : " << factor_10t << "_" << factor_11t << " \n";
					siang_record <<"facotr_12_facotr_13 : " << factor_12t << "_" << factor_13t << " \n";
					siang_record <<"facotr_14_facotr_15 : " << factor_14t << "_" << factor_15t << " \n";
					//-------------------------------------------------------

				    /*DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp]<<"\n";
				    DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp]<<"\n";
					DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<<A_B0R2[ma_tmp]<<"\n";
					DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<<A_B0R3[ma_tmp]<<"\n";
					DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<<A_B0R4[ma_tmp]<<"\n";
					DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<<A_B0R5[ma_tmp]<<"\n";
					DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<<A_B0R6[ma_tmp]<<"\n";
					DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<<A_B0R7[ma_tmp]<<"\n";
					DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<<A_B0R8[ma_tmp]<<"\n";
					DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<<A_B0R9[ma_tmp]<<"\n";
					DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<<A_B0R10[ma_tmp]<<"\n";
					DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<<A_B0R11[ma_tmp]<<"\n";
					DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<<A_B0R12[ma_tmp]<<"\n";
					DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<<A_B0R13[ma_tmp]<<"\n";
					DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<<A_B0R14[ma_tmp]<<"\n";
					DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<<A_B0R15[ma_tmp]<<"\n";*/
					DATARECORD <<"--------------------------------------------------------------------\n";
					siang_record <<"--------------------------------------------------------------------\n";
					if(j <  2)bn0_ma_reg1 = ma_tmp;
					if((j >= 2)  && (j < 4))bn0_ma_reg2 = ma_tmp;
					if((j >= 4)  && (j < 6))bn0_ma_reg3 = ma_tmp;
					if((j >= 6)  && (j < 8))bn0_ma_reg4 = ma_tmp;
					if((j >= 8)  && (j < 10))bn0_ma_reg5 = ma_tmp;
					if((j >= 10) && (j < 12))bn0_ma_reg6 = ma_tmp;
					if((j >= 12) && (j < 14))bn0_ma_reg7 = ma_tmp;
					if((j >= 14) && (j < 16))bn0_ma_reg8 = ma_tmp;
				}
			    else {
					if(j < 2)bn1_bc_tmp = BC_tmp;
					PowerMod(factor_2t,factor_t,2,p);
					PowerMod(factor_3t,factor_t,3,p);
					PowerMod(factor_4t,factor_t,4,p);
					PowerMod(factor_5t,factor_t,5,p);
					PowerMod(factor_6t,factor_t,6,p);
					PowerMod(factor_7t,factor_t,7,p);
					PowerMod(factor_8t,factor_t,8,p);
					PowerMod(factor_9t,factor_t,9,p);
					PowerMod(factor_10t,factor_t,10,p);
					PowerMod(factor_11t,factor_t,11,p);
					PowerMod(factor_12t,factor_t,12,p);
					PowerMod(factor_13t,factor_t,13,p);
					PowerMod(factor_14t,factor_t,14,p);
					PowerMod(factor_15t,factor_t,15,p);
					DATARECORD <<" Before butterfly unit operation! \n";
					/*DATARECORD <<" A_B1R0["<<ma_tmp<<"]:  "<<A_B1R0[ma_tmp]<<"\n";
				    DATARECORD <<" A_B1R1["<<ma_tmp<<"]:  "<<A_B1R1[ma_tmp]<<"\n";
					DATARECORD <<" A_B1R2["<<ma_tmp<<"]:  "<<A_B1R2[ma_tmp]<<"\n";
					DATARECORD <<" A_B1R3["<<ma_tmp<<"]:  "<<A_B1R3[ma_tmp]<<"\n";
					DATARECORD <<" A_B1R4["<<ma_tmp<<"]:  "<<A_B1R4[ma_tmp]<<"\n";
					DATARECORD <<" A_B1R5["<<ma_tmp<<"]:  "<<A_B1R5[ma_tmp]<<"\n";
					DATARECORD <<" A_B1R6["<<ma_tmp<<"]:  "<<A_B1R6[ma_tmp]<<"\n";
					DATARECORD <<" A_B1R7["<<ma_tmp<<"]:  "<<A_B1R7[ma_tmp]<<"\n";
					DATARECORD <<" A_B1R8["<<ma_tmp<<"]:  "<<A_B1R8[ma_tmp]<<"\n";
					DATARECORD <<" A_B1R9["<<ma_tmp<<"]:  "<<A_B1R9[ma_tmp]<<"\n";
					DATARECORD <<" A_B1R10["<<ma_tmp<<"]: "<<A_B1R10[ma_tmp]<<"\n";
					DATARECORD <<" A_B1R11["<<ma_tmp<<"]: "<<A_B1R11[ma_tmp]<<"\n";
					DATARECORD <<" A_B1R12["<<ma_tmp<<"]: "<<A_B1R12[ma_tmp]<<"\n";
					DATARECORD <<" A_B1R13["<<ma_tmp<<"]: "<<A_B1R13[ma_tmp]<<"\n";
					DATARECORD <<" A_B1R14["<<ma_tmp<<"]: "<<A_B1R14[ma_tmp]<<"\n";
					DATARECORD <<" A_B1R15["<<ma_tmp<<"]: "<<A_B1R15[ma_tmp]<<"\n";*/
					Radix16_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp], A_B1R2[ma_tmp], A_B1R3[ma_tmp],A_B1R4[ma_tmp],
							   A_B1R5[ma_tmp],A_B1R6[ma_tmp], A_B1R7[ma_tmp], A_B1R8[ma_tmp],A_B1R9[ma_tmp],
							   A_B1R10[ma_tmp],A_B1R11[ma_tmp],A_B1R12[ma_tmp],A_B1R13[ma_tmp],A_B1R14[ma_tmp],
							   A_B1R15[ma_tmp]);
					MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],factor_t,p);
					MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],factor_2t,p);
					MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],factor_3t,p);
					MulMod(A_B1R4[ma_tmp],A_B1R4[ma_tmp],factor_4t,p);
					MulMod(A_B1R5[ma_tmp],A_B1R5[ma_tmp],factor_5t,p);
					MulMod(A_B1R6[ma_tmp],A_B1R6[ma_tmp],factor_6t,p);
					MulMod(A_B1R7[ma_tmp],A_B1R7[ma_tmp],factor_7t,p);
					MulMod(A_B1R8[ma_tmp],A_B1R8[ma_tmp],factor_8t,p);
					MulMod(A_B1R9[ma_tmp],A_B1R9[ma_tmp],factor_9t,p);
					MulMod(A_B1R10[ma_tmp],A_B1R10[ma_tmp],factor_10t,p);
					MulMod(A_B1R11[ma_tmp],A_B1R11[ma_tmp],factor_11t,p);
					MulMod(A_B1R12[ma_tmp],A_B1R12[ma_tmp],factor_12t,p);
					MulMod(A_B1R13[ma_tmp],A_B1R13[ma_tmp],factor_13t,p);
					MulMod(A_B1R14[ma_tmp],A_B1R14[ma_tmp],factor_14t,p);
					MulMod(A_B1R15[ma_tmp],A_B1R15[ma_tmp],factor_15t,p);
					DATARECORD <<"facotr_1  : " << factor_t  << " \n";
					DATARECORD <<"facotr_2  : " << factor_2t << " \n";
					DATARECORD <<"facotr_3  : " << factor_3t << " \n";
					DATARECORD <<"facotr_4  : " << factor_4t << " \n";
					DATARECORD <<"facotr_5  : " << factor_5t << " \n";
					DATARECORD <<"facotr_6  : " << factor_6t << " \n";
					DATARECORD <<"facotr_7  : " << factor_7t << " \n";
					DATARECORD <<"facotr_8  : " << factor_8t << " \n";
					DATARECORD <<"facotr_9  : " << factor_9t << " \n";
					DATARECORD <<"facotr_10 : " << factor_10t << " \n";
					DATARECORD <<"facotr_11 : " << factor_11t << " \n";
					DATARECORD <<"facotr_12 : " << factor_12t << " \n";
					DATARECORD <<"facotr_13 : " << factor_13t << " \n";
					DATARECORD <<"facotr_14 : " << factor_14t << " \n";
					DATARECORD <<"facotr_15 : " << factor_15t << " \n";

					//----------------DTFAG golden pattern------------------
					switch(s){
						case 0:
							DTFAG_golden_st0 << 1 		  << " \n ";
							DTFAG_golden_st0 << factor_t  << " \n ";
							DTFAG_golden_st0 << factor_2t << " \n ";
							DTFAG_golden_st0 << factor_3t << " \n ";
							DTFAG_golden_st0 << factor_4t << " \n ";
							DTFAG_golden_st0 << factor_5t << " \n ";
							DTFAG_golden_st0 << factor_6t << " \n ";
							DTFAG_golden_st0 << factor_7t << " \n ";
							DTFAG_golden_st0 << factor_8t << " \n ";
							DTFAG_golden_st0 << factor_9t << " \n ";
							DTFAG_golden_st0 << factor_10t << " \n ";
							DTFAG_golden_st0 << factor_11t << " \n ";
							DTFAG_golden_st0 << factor_12t << " \n ";
							DTFAG_golden_st0 << factor_13t << " \n ";
							DTFAG_golden_st0 << factor_14t << " \n ";
							DTFAG_golden_st0 << factor_15t << " \n ";
							break;
						case 1:
							DTFAG_golden_st1 << 1 		  << " \n ";
							DTFAG_golden_st1 << factor_t  << " \n ";
							DTFAG_golden_st1 << factor_2t << " \n ";
							DTFAG_golden_st1 << factor_3t << " \n ";
							DTFAG_golden_st1 << factor_4t << " \n ";
							DTFAG_golden_st1 << factor_5t << " \n ";
							DTFAG_golden_st1 << factor_6t << " \n ";
							DTFAG_golden_st1 << factor_7t << " \n ";
							DTFAG_golden_st1 << factor_8t << " \n ";
							DTFAG_golden_st1 << factor_9t << " \n ";
							DTFAG_golden_st1 << factor_10t << " \n ";
							DTFAG_golden_st1 << factor_11t << " \n ";
							DTFAG_golden_st1 << factor_12t << " \n ";
							DTFAG_golden_st1 << factor_13t << " \n ";
							DTFAG_golden_st1 << factor_14t << " \n ";
							DTFAG_golden_st1 << factor_15t << " \n ";
							break;
						case 2:
							DTFAG_golden_st2 << 1 		  << " \n ";
							DTFAG_golden_st2 << factor_t  << " \n ";
							DTFAG_golden_st2 << factor_2t << " \n ";
							DTFAG_golden_st2 << factor_3t << " \n ";
							DTFAG_golden_st2 << factor_4t << " \n ";
							DTFAG_golden_st2 << factor_5t << " \n ";
							DTFAG_golden_st2 << factor_6t << " \n ";
							DTFAG_golden_st2 << factor_7t << " \n ";
							DTFAG_golden_st2 << factor_8t << " \n ";
							DTFAG_golden_st2 << factor_9t << " \n ";
							DTFAG_golden_st2 << factor_10t << " \n ";
							DTFAG_golden_st2 << factor_11t << " \n ";
							DTFAG_golden_st2 << factor_12t << " \n ";
							DTFAG_golden_st2 << factor_13t << " \n ";
							DTFAG_golden_st2 << factor_14t << " \n ";
							DTFAG_golden_st2 << factor_15t << " \n ";
							break;
						default:
							break;
					}
					//------------------------------------------------------

					//-------------------siang record----------------------.
					//siang_record << "factor_t: "<< factor_t <<"\n";	// siang
					siang_record << "p: "<< p <<"\n";	// siang
					siang_record <<"facotr_1  : " << factor_t  << " \n";
					siang_record <<"facotr_2_facotr_3 : " << factor_2t << "_" << factor_3t << " \n";
					siang_record <<"facotr_4_facotr_5 : " << factor_4t << "_" << factor_5t << " \n";
					siang_record <<"facotr_6_facotr_7 : " << factor_6t << "_" << factor_7t << " \n";
					siang_record <<"facotr_8_facotr_9 : " << factor_8t << "_" << factor_9t << " \n";
					siang_record <<"facotr_10_facotr_11 : " << factor_10t << "_" << factor_11t << " \n";
					siang_record <<"facotr_12_facotr_13 : " << factor_12t << "_" << factor_13t << " \n";
					siang_record <<"facotr_14_facotr_15 : " << factor_14t << "_" << factor_15t << " \n";
					//-------------------------------------------------------



				    /*DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<<"\n";
				    DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<<"\n";
					DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<<A_B1R2[ma_tmp]<<"\n";
					DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<<A_B1R3[ma_tmp]<<"\n";
					DATARECORD <<"A_B1R4["<<ma_tmp<<"]: "<<A_B1R4[ma_tmp]<<"\n";
					DATARECORD <<"A_B1R5["<<ma_tmp<<"]: "<<A_B1R5[ma_tmp]<<"\n";
					DATARECORD <<"A_B1R6["<<ma_tmp<<"]: "<<A_B1R6[ma_tmp]<<"\n";
					DATARECORD <<"A_B1R7["<<ma_tmp<<"]: "<<A_B1R7[ma_tmp]<<"\n";
					DATARECORD <<"A_B1R8["<<ma_tmp<<"]: "<<A_B1R8[ma_tmp]<<"\n";
					DATARECORD <<"A_B1R9["<<ma_tmp<<"]: "<<A_B1R9[ma_tmp]<<"\n";
					DATARECORD <<"A_B1R10["<<ma_tmp<<"]: "<<A_B1R10[ma_tmp]<<"\n";
					DATARECORD <<"A_B1R11["<<ma_tmp<<"]: "<<A_B1R11[ma_tmp]<<"\n";
					DATARECORD <<"A_B1R12["<<ma_tmp<<"]: "<<A_B1R12[ma_tmp]<<"\n";
					DATARECORD <<"A_B1R13["<<ma_tmp<<"]: "<<A_B1R13[ma_tmp]<<"\n";
					DATARECORD <<"A_B1R14["<<ma_tmp<<"]: "<<A_B1R14[ma_tmp]<<"\n";
					DATARECORD <<"A_B1R15["<<ma_tmp<<"]: "<<A_B1R15[ma_tmp]<<"\n";*/
					DATARECORD <<"--------------------------------------------------------------------\n";
					siang_record <<"--------------------------------------------------------------------\n";
                    if(j <  2)bn1_ma_reg1 = ma_tmp;					
                    if((j >= 2)  && (j < 4))bn1_ma_reg2 = ma_tmp;
                    if((j >= 4)  && (j < 6))bn1_ma_reg3 = ma_tmp;
                    if((j >= 6)  && (j < 8))bn1_ma_reg4 = ma_tmp;
                    if((j >= 8)  && (j < 10))bn1_ma_reg5 = ma_tmp;
                    if((j >= 10) && (j < 12))bn1_ma_reg6 = ma_tmp;
                    if((j >= 12) && (j < 14))bn1_ma_reg7 = ma_tmp;
                    if((j >= 14) && (j < 16))bn1_ma_reg8 = ma_tmp;					
				}
			}
		//data relocation
		 if(s < Stage-1){
		  if(bn1_bc_tmp > bn0_bc_tmp){
			 data_tmp_1  = A_B0R1[bn0_ma_reg1];
			 data_tmp_2  = A_B0R2[bn0_ma_reg1];
			 data_tmp_3  = A_B0R3[bn0_ma_reg1];
			 data_tmp_4  = A_B0R4[bn0_ma_reg1];
			 data_tmp_5  = A_B0R5[bn0_ma_reg1];
			 data_tmp_6  = A_B0R6[bn0_ma_reg1];
			 data_tmp_7  = A_B0R7[bn0_ma_reg1];
			 data_tmp_8  = A_B0R8[bn0_ma_reg1];
			 data_tmp_9  = A_B0R9[bn0_ma_reg1];
			 data_tmp_10 = A_B0R10[bn0_ma_reg1];
			 data_tmp_11 = A_B0R11[bn0_ma_reg1];
			 data_tmp_12 = A_B0R12[bn0_ma_reg1];
			 data_tmp_13 = A_B0R13[bn0_ma_reg1];
			 data_tmp_14 = A_B0R14[bn0_ma_reg1];
			 data_tmp_15 = A_B0R15[bn0_ma_reg1];
			 A_B0R1[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg1];
			 A_B0R2[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg2];
			 A_B0R3[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg2];
			 A_B0R4[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg3];
			 A_B0R5[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg3];
			 A_B0R6[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg1] = A_B0R0[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg1] = A_B1R0[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg1] = A_B0R0[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg1] = A_B1R0[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg1] = A_B1R0[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg1] = A_B0R0[bn0_ma_reg8];
			 A_B1R0[bn1_ma_reg1]  = data_tmp_1; 
			 A_B1R0[bn1_ma_reg2]  = data_tmp_2; 
			 A_B0R0[bn0_ma_reg2]  = data_tmp_3; 
			 A_B1R0[bn1_ma_reg3]  = data_tmp_4; 
			 A_B0R0[bn0_ma_reg3]  = data_tmp_5; 
			 A_B0R0[bn0_ma_reg4]  = data_tmp_6; 
			 A_B1R0[bn1_ma_reg4]  = data_tmp_7; 
			 A_B1R0[bn1_ma_reg5]  = data_tmp_8; 
			 A_B0R0[bn0_ma_reg5]  = data_tmp_9; 
			 A_B0R0[bn0_ma_reg6]  = data_tmp_10;
			 A_B1R0[bn1_ma_reg6]  = data_tmp_11;
			 A_B0R0[bn0_ma_reg7]  = data_tmp_12;
			 A_B1R0[bn1_ma_reg7]  = data_tmp_13;
			 A_B1R0[bn1_ma_reg8]  = data_tmp_14;
			 A_B0R0[bn0_ma_reg8]  = data_tmp_15;
			 /*********************************************/
			 data_tmp_1  = A_B1R2[bn1_ma_reg1];
			 data_tmp_2  = A_B1R3[bn1_ma_reg1];
			 data_tmp_3  = A_B1R4[bn1_ma_reg1];
			 data_tmp_4  = A_B1R5[bn1_ma_reg1];
			 data_tmp_5  = A_B1R6[bn1_ma_reg1];
			 data_tmp_6  = A_B1R7[bn1_ma_reg1];
			 data_tmp_7  = A_B1R8[bn1_ma_reg1];
			 data_tmp_8  = A_B1R9[bn1_ma_reg1];
			 data_tmp_9  = A_B1R10[bn1_ma_reg1];
			 data_tmp_10 = A_B1R11[bn1_ma_reg1];
			 data_tmp_11 = A_B1R12[bn1_ma_reg1];
			 data_tmp_12 = A_B1R13[bn1_ma_reg1];
			 data_tmp_13 = A_B1R14[bn1_ma_reg1];
			 data_tmp_14 = A_B1R15[bn1_ma_reg1];
             A_B1R2[bn1_ma_reg1]  =	 A_B1R1[bn1_ma_reg2];
             A_B1R3[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg2];
             A_B1R4[bn1_ma_reg1]  =  A_B1R1[bn1_ma_reg3];
             A_B1R5[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg3];
             A_B1R6[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg4];
             A_B1R7[bn1_ma_reg1]  =  A_B1R1[bn1_ma_reg4];
             A_B1R8[bn1_ma_reg1]  =  A_B1R1[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg1] =  A_B0R1[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg1] =  A_B1R1[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg1] =  A_B0R1[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg1] =  A_B1R1[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg1] =  A_B1R1[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg1] =  A_B0R1[bn0_ma_reg8];
			 A_B1R1[bn1_ma_reg2]  =  data_tmp_1; 
			 A_B0R1[bn0_ma_reg2]  =  data_tmp_2; 
			 A_B1R1[bn1_ma_reg3]  =  data_tmp_3; 
			 A_B0R1[bn0_ma_reg3]  =  data_tmp_4; 
			 A_B0R1[bn0_ma_reg4]  =  data_tmp_5; 
			 A_B1R1[bn1_ma_reg4]  =  data_tmp_6; 
			 A_B1R1[bn1_ma_reg5]  =  data_tmp_7; 
			 A_B0R1[bn0_ma_reg5]  =  data_tmp_8; 
			 A_B0R1[bn0_ma_reg6]  =  data_tmp_9; 
			 A_B1R1[bn1_ma_reg6]  =  data_tmp_10;
			 A_B0R1[bn0_ma_reg7]  =  data_tmp_11;
			 A_B1R1[bn1_ma_reg7]  =  data_tmp_12;
			 A_B1R1[bn1_ma_reg8]  =  data_tmp_13;
			 A_B0R1[bn0_ma_reg8]  =  data_tmp_14;
			/************************************************************/ 
			 data_tmp_1  =   A_B1R3[bn1_ma_reg2];
			 data_tmp_2  =   A_B1R4[bn1_ma_reg2];
			 data_tmp_3  =   A_B1R5[bn1_ma_reg2];
			 data_tmp_4  =   A_B1R6[bn1_ma_reg2];
			 data_tmp_5  =   A_B1R7[bn1_ma_reg2];
			 data_tmp_6  =   A_B1R8[bn1_ma_reg2];
			 data_tmp_7  =   A_B1R9[bn1_ma_reg2];
			 data_tmp_8  =   A_B1R10[bn1_ma_reg2];
			 data_tmp_9  =   A_B1R11[bn1_ma_reg2];
			 data_tmp_10 =   A_B1R12[bn1_ma_reg2];
			 data_tmp_11 =   A_B1R13[bn1_ma_reg2];
			 data_tmp_12 =   A_B1R14[bn1_ma_reg2];
			 data_tmp_13 =   A_B1R15[bn1_ma_reg2];
			 A_B1R3[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg2];
			 A_B1R4[bn1_ma_reg2]  =  A_B1R2[bn1_ma_reg3];
			 A_B1R5[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg3];
			 A_B1R6[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg4];
			 A_B1R7[bn1_ma_reg2]  =  A_B1R2[bn1_ma_reg4];
			 A_B1R8[bn1_ma_reg2]  =  A_B1R2[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg2] =  A_B0R2[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg2] =  A_B1R2[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg2] =  A_B0R2[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg2] =  A_B1R2[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg2] =  A_B1R2[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg2] =  A_B0R2[bn0_ma_reg8];
			 A_B0R2[bn0_ma_reg2]  =  data_tmp_1; 
			 A_B1R2[bn1_ma_reg3]  =  data_tmp_2; 
			 A_B0R2[bn0_ma_reg3]  =  data_tmp_3; 
			 A_B0R2[bn0_ma_reg4]  =  data_tmp_4; 
			 A_B1R2[bn1_ma_reg4]  =  data_tmp_5; 
			 A_B1R2[bn1_ma_reg5]  =  data_tmp_6; 
			 A_B0R2[bn0_ma_reg5]  =  data_tmp_7; 
			 A_B0R2[bn0_ma_reg6]  =  data_tmp_8; 
			 A_B1R2[bn1_ma_reg6]  =  data_tmp_9; 
			 A_B0R2[bn0_ma_reg7]  =  data_tmp_10;
			 A_B1R2[bn1_ma_reg7]  =  data_tmp_11;
			 A_B1R2[bn1_ma_reg8]  =  data_tmp_12;
			 A_B0R2[bn0_ma_reg8]  =  data_tmp_13;
			 //******************************************************
			 data_tmp_1  = A_B0R4[bn0_ma_reg2];
			 data_tmp_2  = A_B0R5[bn0_ma_reg2];
			 data_tmp_3  = A_B0R6[bn0_ma_reg2];
			 data_tmp_4  = A_B0R7[bn0_ma_reg2];
			 data_tmp_5  = A_B0R8[bn0_ma_reg2];
			 data_tmp_6  = A_B0R9[bn0_ma_reg2];
			 data_tmp_7  = A_B0R10[bn0_ma_reg2];
			 data_tmp_8  = A_B0R11[bn0_ma_reg2];
			 data_tmp_9  = A_B0R12[bn0_ma_reg2];
			 data_tmp_10 = A_B0R13[bn0_ma_reg2];
			 data_tmp_11 = A_B0R14[bn0_ma_reg2];
			 data_tmp_12 = A_B0R15[bn0_ma_reg2];
			 A_B0R4[bn0_ma_reg2]  =  A_B1R3[bn1_ma_reg3];
			 A_B0R5[bn0_ma_reg2]  =  A_B0R3[bn0_ma_reg3];
			 A_B0R6[bn0_ma_reg2]  =  A_B0R3[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg2]  =  A_B1R3[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg2]  =  A_B1R3[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg2]  =  A_B0R3[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg2] =  A_B0R3[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg2] =  A_B1R3[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg2] =  A_B0R3[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg2] =  A_B1R3[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg2] =  A_B1R3[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg2] =  A_B0R3[bn0_ma_reg8];
			 A_B1R3[bn1_ma_reg3]  =  data_tmp_1;  
			 A_B0R3[bn0_ma_reg3]  =  data_tmp_2;  
			 A_B0R3[bn0_ma_reg4]  =  data_tmp_3;  
			 A_B1R3[bn1_ma_reg4]  =  data_tmp_4;  
			 A_B1R3[bn1_ma_reg5]  =  data_tmp_5;  
			 A_B0R3[bn0_ma_reg5]  =  data_tmp_6;  
			 A_B0R3[bn0_ma_reg6]  =  data_tmp_7;  
			 A_B1R3[bn1_ma_reg6]  =  data_tmp_8;  
			 A_B0R3[bn0_ma_reg7]  =  data_tmp_9;  
			 A_B1R3[bn1_ma_reg7]  =  data_tmp_10; 
			 A_B1R3[bn1_ma_reg8]  =  data_tmp_11; 
			 A_B0R3[bn0_ma_reg8]  =  data_tmp_12; 
			 //----------------------------------------------------------------------
             data_tmp_1  = A_B1R5[bn1_ma_reg3];
			 data_tmp_2  = A_B1R6[bn1_ma_reg3];
			 data_tmp_3  = A_B1R7[bn1_ma_reg3];
			 data_tmp_4  = A_B1R8[bn1_ma_reg3];
			 data_tmp_5  = A_B1R9[bn1_ma_reg3];
			 data_tmp_6  = A_B1R10[bn1_ma_reg3];
			 data_tmp_7  = A_B1R11[bn1_ma_reg3];
			 data_tmp_8  = A_B1R12[bn1_ma_reg3];
			 data_tmp_9  = A_B1R13[bn1_ma_reg3];
			 data_tmp_10 = A_B1R14[bn1_ma_reg3];
			 data_tmp_11 = A_B1R15[bn1_ma_reg3];
			 A_B1R5[bn1_ma_reg3]  =  A_B0R4[bn0_ma_reg3];
			 A_B1R6[bn1_ma_reg3]  =  A_B0R4[bn0_ma_reg4];
			 A_B1R7[bn1_ma_reg3]  =  A_B1R4[bn1_ma_reg4];
			 A_B1R8[bn1_ma_reg3]  =  A_B1R4[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg3]  =  A_B0R4[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg3] =  A_B0R4[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg3] =  A_B1R4[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg3] =  A_B0R4[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg3] =  A_B1R4[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg3] =  A_B1R4[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg3] =  A_B0R4[bn0_ma_reg8];
			 A_B0R4[bn0_ma_reg3]  =  data_tmp_1 ;
			 A_B0R4[bn0_ma_reg4]  =  data_tmp_2 ;
			 A_B1R4[bn1_ma_reg4]  =  data_tmp_3 ;
			 A_B1R4[bn1_ma_reg5]  =  data_tmp_4 ;
			 A_B0R4[bn0_ma_reg5]  =  data_tmp_5 ;
			 A_B0R4[bn0_ma_reg6]  =  data_tmp_6 ;
			 A_B1R4[bn1_ma_reg6]  =  data_tmp_7 ;
			 A_B0R4[bn0_ma_reg7]  =  data_tmp_8 ;
			 A_B1R4[bn1_ma_reg7]  =  data_tmp_9 ;
			 A_B1R4[bn1_ma_reg8]  =  data_tmp_10;
			 A_B0R4[bn0_ma_reg8]  =  data_tmp_11;
			 //************************************************************************
			 data_tmp_1  = A_B0R6[bn0_ma_reg3];
			 data_tmp_2  = A_B0R7[bn0_ma_reg3];
			 data_tmp_3  = A_B0R8[bn0_ma_reg3];
			 data_tmp_4  = A_B0R9[bn0_ma_reg3];
			 data_tmp_5  = A_B0R10[bn0_ma_reg3];
			 data_tmp_6  = A_B0R11[bn0_ma_reg3];
			 data_tmp_7  = A_B0R12[bn0_ma_reg3];
			 data_tmp_8  = A_B0R13[bn0_ma_reg3];
			 data_tmp_9  = A_B0R14[bn0_ma_reg3];
			 data_tmp_10 = A_B0R15[bn0_ma_reg3];
			 A_B0R6[bn0_ma_reg3]  =  A_B0R5[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg3]  =  A_B1R5[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg3]  =  A_B1R5[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg3]  =  A_B0R5[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg3] =  A_B0R5[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg3] =  A_B1R5[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg3] =  A_B0R5[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg3] =  A_B1R5[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg3] =  A_B1R5[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg3] =  A_B0R5[bn0_ma_reg8];
			 A_B0R5[bn0_ma_reg4]  =  data_tmp_1; 
			 A_B1R5[bn1_ma_reg4]  =  data_tmp_2; 
			 A_B1R5[bn1_ma_reg5]  =  data_tmp_3; 
			 A_B0R5[bn0_ma_reg5]  =  data_tmp_4; 
			 A_B0R5[bn0_ma_reg6]  =  data_tmp_5; 
			 A_B1R5[bn1_ma_reg6]  =  data_tmp_6; 
			 A_B0R5[bn0_ma_reg7]  =  data_tmp_7; 
			 A_B1R5[bn1_ma_reg7]  =  data_tmp_8; 
			 A_B1R5[bn1_ma_reg8]  =  data_tmp_9; 
			 A_B0R5[bn0_ma_reg8]  =  data_tmp_10;
			 //-----------------------------------------------------------------------
			 data_tmp_1  = A_B0R7[bn0_ma_reg4];
			 data_tmp_2  = A_B0R8[bn0_ma_reg4];
			 data_tmp_3  = A_B0R9[bn0_ma_reg4];
			 data_tmp_4  = A_B0R10[bn0_ma_reg4];
			 data_tmp_5  = A_B0R11[bn0_ma_reg4];
			 data_tmp_6  = A_B0R12[bn0_ma_reg4];
			 data_tmp_7  = A_B0R13[bn0_ma_reg4];
			 data_tmp_8  = A_B0R14[bn0_ma_reg4];
			 data_tmp_9  = A_B0R15[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg4]  =  A_B1R6[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg4]  =  A_B1R6[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg4]  =  A_B0R6[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg4] =  A_B0R6[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg4] =  A_B1R6[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg4] =  A_B0R6[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg4] =  A_B1R6[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg4] =  A_B1R6[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg4] =  A_B0R6[bn0_ma_reg8];
			 A_B1R6[bn1_ma_reg4]  =  data_tmp_1;
			 A_B1R6[bn1_ma_reg5]  =  data_tmp_2;
			 A_B0R6[bn0_ma_reg5]  =  data_tmp_3;
			 A_B0R6[bn0_ma_reg6]  =  data_tmp_4;
			 A_B1R6[bn1_ma_reg6]  =  data_tmp_5;
			 A_B0R6[bn0_ma_reg7]  =  data_tmp_6;
			 A_B1R6[bn1_ma_reg7]  =  data_tmp_7;
			 A_B1R6[bn1_ma_reg8]  =  data_tmp_8;
			 A_B0R6[bn0_ma_reg8]  =  data_tmp_9;
			 //----------------------------------------------------------------------
			 data_tmp_1  = A_B1R8[bn1_ma_reg4];
			 data_tmp_2  = A_B1R9[bn1_ma_reg4];
			 data_tmp_3  = A_B1R10[bn1_ma_reg4];
			 data_tmp_4  = A_B1R11[bn1_ma_reg4];
			 data_tmp_5  = A_B1R12[bn1_ma_reg4];
			 data_tmp_6  = A_B1R13[bn1_ma_reg4];
			 data_tmp_7  = A_B1R14[bn1_ma_reg4];
			 data_tmp_8  = A_B1R15[bn1_ma_reg4];
			 A_B1R8[bn1_ma_reg4]  =  A_B1R7[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg4]  =  A_B0R7[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg4] =  A_B0R7[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg4] =  A_B1R7[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg4] =  A_B0R7[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg4] =  A_B1R7[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg4] =  A_B1R7[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg4] =  A_B0R7[bn0_ma_reg8];
			 A_B1R7[bn1_ma_reg5]  =  data_tmp_1;
			 A_B0R7[bn0_ma_reg5]  =  data_tmp_2;
			 A_B0R7[bn0_ma_reg6]  =  data_tmp_3;
			 A_B1R7[bn1_ma_reg6]  =  data_tmp_4;
			 A_B0R7[bn0_ma_reg7]  =  data_tmp_5;
			 A_B1R7[bn1_ma_reg7]  =  data_tmp_6;
			 A_B1R7[bn1_ma_reg8]  =  data_tmp_7;
			 A_B0R7[bn0_ma_reg8]  =  data_tmp_8;
			 //----------------------------------------------------------------------
			 data_tmp_1 = A_B1R9[bn1_ma_reg5];
			 data_tmp_2 = A_B1R10[bn1_ma_reg5];
			 data_tmp_3 = A_B1R11[bn1_ma_reg5];
			 data_tmp_4 = A_B1R12[bn1_ma_reg5];
			 data_tmp_5 = A_B1R13[bn1_ma_reg5];
			 data_tmp_6 = A_B1R14[bn1_ma_reg5];
			 data_tmp_7 = A_B1R15[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg5]  = A_B0R8[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg5] = A_B0R8[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg5] = A_B1R8[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg5] = A_B0R8[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg5] = A_B1R8[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg5] = A_B1R8[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg5] = A_B0R8[bn0_ma_reg8];
			 A_B0R8[bn0_ma_reg5]  = data_tmp_1;
			 A_B0R8[bn0_ma_reg6]  = data_tmp_2;
			 A_B1R8[bn1_ma_reg6]  = data_tmp_3;
			 A_B0R8[bn0_ma_reg7]  = data_tmp_4;
			 A_B1R8[bn1_ma_reg7]  = data_tmp_5;
			 A_B1R8[bn1_ma_reg8]  = data_tmp_6;
			 A_B0R8[bn0_ma_reg8]  = data_tmp_7;
			 //---------------------------------------------------------------------
			 data_tmp_1  = A_B0R10[bn0_ma_reg5];
			 data_tmp_2  = A_B0R11[bn0_ma_reg5];
			 data_tmp_3  = A_B0R12[bn0_ma_reg5];
			 data_tmp_4  = A_B0R13[bn0_ma_reg5];
			 data_tmp_5  = A_B0R14[bn0_ma_reg5];
			 data_tmp_6  = A_B0R15[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg5] = A_B0R9[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg5] = A_B1R9[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg5] = A_B0R9[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg5] = A_B1R9[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg5] = A_B1R9[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg5] = A_B0R9[bn0_ma_reg8];
			 A_B0R9[bn0_ma_reg6]  = data_tmp_1;
			 A_B1R9[bn1_ma_reg6]  = data_tmp_2;
			 A_B0R9[bn0_ma_reg7]  = data_tmp_3;
			 A_B1R9[bn1_ma_reg7]  = data_tmp_4;
			 A_B1R9[bn1_ma_reg8]  = data_tmp_5;
			 A_B0R9[bn0_ma_reg8]  = data_tmp_6;
			 //--------------------------------------------------------------------
			 data_tmp_1  = A_B0R11[bn0_ma_reg6];
			 data_tmp_2  = A_B0R12[bn0_ma_reg6];
			 data_tmp_3  = A_B0R13[bn0_ma_reg6];
			 data_tmp_4  = A_B0R14[bn0_ma_reg6];
			 data_tmp_5  = A_B0R15[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg6] = A_B1R10[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg6] = A_B0R10[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg6] = A_B1R10[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg6] = A_B1R10[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg6] = A_B0R10[bn0_ma_reg8];
			 A_B1R10[bn1_ma_reg6] = data_tmp_1;
			 A_B0R10[bn0_ma_reg7] = data_tmp_2;
			 A_B1R10[bn1_ma_reg7] = data_tmp_3;
			 A_B1R10[bn1_ma_reg8] = data_tmp_4;
			 A_B0R10[bn0_ma_reg8] = data_tmp_5;
			 //--------------------------------------------------------------------
			 data_tmp_1  = A_B1R12[bn1_ma_reg6];
			 data_tmp_2  = A_B1R13[bn1_ma_reg6];
			 data_tmp_3  = A_B1R14[bn1_ma_reg6];
			 data_tmp_4  = A_B1R15[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg6] = A_B0R11[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg6] = A_B1R11[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg6] = A_B1R11[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg6] = A_B0R11[bn0_ma_reg8];
			 A_B0R11[bn0_ma_reg7] = data_tmp_1;
			 A_B1R11[bn1_ma_reg7] = data_tmp_2;
			 A_B1R11[bn1_ma_reg8] = data_tmp_3;
			 A_B0R11[bn0_ma_reg8] = data_tmp_4;
			 //--------------------------------------------------------------------
			 data_tmp_1 = A_B0R13[bn0_ma_reg7];
			 data_tmp_2 = A_B0R14[bn0_ma_reg7];
			 data_tmp_3 = A_B0R15[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg7] = A_B1R12[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg7] = A_B1R12[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg7] = A_B0R12[bn0_ma_reg8];
			 A_B1R12[bn1_ma_reg7] = data_tmp_1;
			 A_B1R12[bn1_ma_reg8] = data_tmp_2;
			 A_B0R12[bn0_ma_reg8] = data_tmp_3;
			 //--------------------------------------------------------------------
			 data_tmp_1 = A_B1R14[bn1_ma_reg7];
			 data_tmp_2 = A_B1R15[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg7] = A_B1R13[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg7] = A_B0R13[bn0_ma_reg8];
			 A_B1R13[bn1_ma_reg8] = data_tmp_1;
			 A_B0R13[bn0_ma_reg8] = data_tmp_2;
			 //--------------------------------------------------------------------
			 data_tmp_1 = A_B1R15[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg8] = A_B0R14[bn0_ma_reg8];
			 A_B0R14[bn0_ma_reg8] = data_tmp_1;
			 
		  }else
		  {	 
			 data_tmp_1  = A_B1R1[bn1_ma_reg1];
			 data_tmp_2  = A_B1R2[bn1_ma_reg1];
			 data_tmp_3  = A_B1R3[bn1_ma_reg1];
			 data_tmp_4  = A_B1R4[bn1_ma_reg1];
			 data_tmp_5  = A_B1R5[bn1_ma_reg1];
			 data_tmp_6  = A_B1R6[bn1_ma_reg1];
			 data_tmp_7  = A_B1R7[bn1_ma_reg1];
			 data_tmp_8  = A_B1R8[bn1_ma_reg1];
			 data_tmp_9  = A_B1R9[bn1_ma_reg1];
			 data_tmp_10 = A_B1R10[bn1_ma_reg1];
			 data_tmp_11 = A_B1R11[bn1_ma_reg1];
			 data_tmp_12 = A_B1R12[bn1_ma_reg1];
			 data_tmp_13 = A_B1R13[bn1_ma_reg1];
			 data_tmp_14 = A_B1R14[bn1_ma_reg1];
			 data_tmp_15 = A_B1R15[bn1_ma_reg1];
             A_B1R1[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg1];
             A_B1R2[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg2];
             A_B1R3[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg2];
             A_B1R4[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg3];
             A_B1R5[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg3];
             A_B1R6[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg4];
             A_B1R7[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg4];
             A_B1R8[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg5];
             A_B1R9[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg5];
             A_B1R10[bn1_ma_reg1] =  A_B1R0[bn1_ma_reg6];
             A_B1R11[bn1_ma_reg1] =  A_B0R0[bn0_ma_reg6];
             A_B1R12[bn1_ma_reg1] =  A_B1R0[bn1_ma_reg7];
             A_B1R13[bn1_ma_reg1] =  A_B0R0[bn0_ma_reg7];
             A_B1R14[bn1_ma_reg1] =  A_B0R0[bn0_ma_reg8];
             A_B1R15[bn1_ma_reg1] =  A_B1R0[bn1_ma_reg8];
             A_B0R0[bn0_ma_reg1]  =  data_tmp_1; 
             A_B0R0[bn0_ma_reg2]  =  data_tmp_2; 
             A_B1R0[bn1_ma_reg2]  =  data_tmp_3; 
             A_B0R0[bn0_ma_reg3]  =  data_tmp_4; 
             A_B1R0[bn1_ma_reg3]  =  data_tmp_5; 
             A_B1R0[bn1_ma_reg4]  =  data_tmp_6; 
             A_B0R0[bn0_ma_reg4]  =  data_tmp_7; 
             A_B0R0[bn0_ma_reg5]  =  data_tmp_8; 
             A_B1R0[bn1_ma_reg5]  =  data_tmp_9; 
             A_B1R0[bn1_ma_reg6]  =  data_tmp_10;
             A_B0R0[bn0_ma_reg6]  =  data_tmp_11;
             A_B1R0[bn1_ma_reg7]  =  data_tmp_12;
             A_B0R0[bn0_ma_reg7]  =  data_tmp_13;
             A_B0R0[bn0_ma_reg8]  =  data_tmp_14;
             A_B1R0[bn1_ma_reg8]  =  data_tmp_15;
             //-----------------------------------------------------------
             data_tmp_1  = A_B0R2[bn0_ma_reg1];
			 data_tmp_2  = A_B0R3[bn0_ma_reg1];
			 data_tmp_3  = A_B0R4[bn0_ma_reg1];
			 data_tmp_4  = A_B0R5[bn0_ma_reg1];
			 data_tmp_5  = A_B0R6[bn0_ma_reg1];
			 data_tmp_6  = A_B0R7[bn0_ma_reg1];
			 data_tmp_7  = A_B0R8[bn0_ma_reg1];
			 data_tmp_8  = A_B0R9[bn0_ma_reg1];
			 data_tmp_9  = A_B0R10[bn0_ma_reg1];
			 data_tmp_10 = A_B0R11[bn0_ma_reg1];
			 data_tmp_11 = A_B0R12[bn0_ma_reg1];
			 data_tmp_12 = A_B0R13[bn0_ma_reg1];
			 data_tmp_13 = A_B0R14[bn0_ma_reg1];
			 data_tmp_14 = A_B0R15[bn0_ma_reg1];
			 A_B0R2[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg2];
			 A_B0R3[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg2];
			 A_B0R4[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg3];
			 A_B0R5[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg3];
			 A_B0R6[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg4];
			 A_B0R7[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg1] =  A_B1R1[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg1] =  A_B0R1[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg1] =  A_B1R1[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg1] =  A_B0R1[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg1] =  A_B0R1[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg1] =  A_B1R1[bn1_ma_reg8];
             A_B0R1[bn0_ma_reg2]  =  data_tmp_1; 
             A_B1R1[bn1_ma_reg2]  =  data_tmp_2; 
             A_B0R1[bn0_ma_reg3]  =  data_tmp_3; 
             A_B1R1[bn1_ma_reg3]  =  data_tmp_4; 
             A_B1R1[bn1_ma_reg4]  =  data_tmp_5; 
             A_B0R1[bn0_ma_reg4]  =  data_tmp_6; 
             A_B0R1[bn0_ma_reg5]  =  data_tmp_7; 
             A_B1R1[bn1_ma_reg5]  =  data_tmp_8; 
             A_B1R1[bn1_ma_reg6]  =  data_tmp_9; 
             A_B0R1[bn0_ma_reg6]  =  data_tmp_10;
             A_B1R1[bn1_ma_reg7]  =  data_tmp_11;
             A_B0R1[bn0_ma_reg7]  =  data_tmp_12;
             A_B0R1[bn0_ma_reg8]  =  data_tmp_13;
             A_B1R1[bn1_ma_reg8]  =  data_tmp_14;
             //------------------------------------------------------------
             data_tmp_1  =  A_B0R3[bn0_ma_reg2];
			 data_tmp_2  =  A_B0R4[bn0_ma_reg2];
			 data_tmp_3  =  A_B0R5[bn0_ma_reg2];
			 data_tmp_4  =  A_B0R6[bn0_ma_reg2];
			 data_tmp_5  =  A_B0R7[bn0_ma_reg2];
			 data_tmp_6  =  A_B0R8[bn0_ma_reg2];
			 data_tmp_7  =  A_B0R9[bn0_ma_reg2];
			 data_tmp_8  =  A_B0R10[bn0_ma_reg2];
			 data_tmp_9  =  A_B0R11[bn0_ma_reg2];
			 data_tmp_10 =  A_B0R12[bn0_ma_reg2];
			 data_tmp_11 =  A_B0R13[bn0_ma_reg2];
			 data_tmp_12 =  A_B0R14[bn0_ma_reg2];
			 data_tmp_13 =  A_B0R15[bn0_ma_reg2];
			 A_B0R3[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg2];
			 A_B0R4[bn0_ma_reg2]  =  A_B0R2[bn0_ma_reg3];
			 A_B0R5[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg3];
			 A_B0R6[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg4];
			 A_B0R7[bn0_ma_reg2]  =  A_B0R2[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg2]  =  A_B0R2[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg2] =  A_B1R2[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg2] =  A_B0R2[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg2] =  A_B1R2[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg2] =  A_B0R2[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg2] =  A_B0R2[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg2] =  A_B1R2[bn1_ma_reg8];
			 A_B1R2[bn1_ma_reg2]  =  data_tmp_1; 
			 A_B0R2[bn0_ma_reg3]  =  data_tmp_2; 
			 A_B1R2[bn1_ma_reg3]  =  data_tmp_3; 
			 A_B1R2[bn1_ma_reg4]  =  data_tmp_4; 
			 A_B0R2[bn0_ma_reg4]  =  data_tmp_5; 
			 A_B0R2[bn0_ma_reg5]  =  data_tmp_6; 
			 A_B1R2[bn1_ma_reg5]  =  data_tmp_7; 
			 A_B1R2[bn1_ma_reg6]  =  data_tmp_8; 
			 A_B0R2[bn0_ma_reg6]  =  data_tmp_9; 
			 A_B1R2[bn1_ma_reg7]  =  data_tmp_10;
			 A_B0R2[bn0_ma_reg7]  =  data_tmp_11;
			 A_B0R2[bn0_ma_reg8]  =  data_tmp_12;
			 A_B1R2[bn1_ma_reg8]  =  data_tmp_13;
			 //-----------------------------------------------------------                        
			 data_tmp_1  = A_B1R4[bn1_ma_reg2];
			 data_tmp_2  = A_B1R5[bn1_ma_reg2];
			 data_tmp_3  = A_B1R6[bn1_ma_reg2];
			 data_tmp_4  = A_B1R7[bn1_ma_reg2];
			 data_tmp_5  = A_B1R8[bn1_ma_reg2];
			 data_tmp_6  = A_B1R9[bn1_ma_reg2];
			 data_tmp_7  = A_B1R10[bn1_ma_reg2];
			 data_tmp_8  = A_B1R11[bn1_ma_reg2];
			 data_tmp_9  = A_B1R12[bn1_ma_reg2];
			 data_tmp_10 = A_B1R13[bn1_ma_reg2];
			 data_tmp_11 = A_B1R14[bn1_ma_reg2];
			 data_tmp_12 = A_B1R15[bn1_ma_reg2];
			 A_B1R4[bn1_ma_reg2]  = A_B0R3[bn0_ma_reg3];
			 A_B1R5[bn1_ma_reg2]  = A_B1R3[bn1_ma_reg3];
			 A_B1R6[bn1_ma_reg2]  = A_B1R3[bn1_ma_reg4];
			 A_B1R7[bn1_ma_reg2]  = A_B0R3[bn0_ma_reg4];
			 A_B1R8[bn1_ma_reg2]  = A_B0R3[bn0_ma_reg5];
			 A_B1R9[bn1_ma_reg2]  = A_B1R3[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg2] = A_B1R3[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg2] = A_B0R3[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg2] = A_B1R3[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg2] = A_B0R3[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg2] = A_B0R3[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg2] = A_B1R3[bn1_ma_reg8];
			 A_B0R3[bn0_ma_reg3]  = data_tmp_1; 
			 A_B1R3[bn1_ma_reg3]  = data_tmp_2; 
			 A_B1R3[bn1_ma_reg4]  = data_tmp_3; 
			 A_B0R3[bn0_ma_reg4]  = data_tmp_4; 
			 A_B0R3[bn0_ma_reg5]  = data_tmp_5; 
			 A_B1R3[bn1_ma_reg5]  = data_tmp_6; 
			 A_B1R3[bn1_ma_reg6]  = data_tmp_7; 
			 A_B0R3[bn0_ma_reg6]  = data_tmp_8; 
			 A_B1R3[bn1_ma_reg7]  = data_tmp_9; 
			 A_B0R3[bn0_ma_reg7]  = data_tmp_10;
			 A_B0R3[bn0_ma_reg8]  = data_tmp_11;
			 A_B1R3[bn1_ma_reg8]  = data_tmp_12;
			 //------------------------------------------------------------
			 data_tmp_1  =  A_B0R5[bn0_ma_reg3];
			 data_tmp_2  =  A_B0R6[bn0_ma_reg3];
			 data_tmp_3  =  A_B0R7[bn0_ma_reg3];
			 data_tmp_4  =  A_B0R8[bn0_ma_reg3];
			 data_tmp_5  =  A_B0R9[bn0_ma_reg3];
			 data_tmp_6  =  A_B0R10[bn0_ma_reg3];
			 data_tmp_7  =  A_B0R11[bn0_ma_reg3];
			 data_tmp_8  =  A_B0R12[bn0_ma_reg3];
			 data_tmp_9  =  A_B0R13[bn0_ma_reg3];
			 data_tmp_10 =  A_B0R14[bn0_ma_reg3];
			 data_tmp_11 =  A_B0R15[bn0_ma_reg3];
			 A_B0R5[bn0_ma_reg3]  =  A_B1R4[bn1_ma_reg3];
			 A_B0R6[bn0_ma_reg3]  =  A_B1R4[bn1_ma_reg4];
			 A_B0R7[bn0_ma_reg3]  =  A_B0R4[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg3]  =  A_B0R4[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg3]  =  A_B1R4[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg3] =  A_B1R4[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg3] =  A_B0R4[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg3] =  A_B1R4[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg3] =  A_B0R4[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg3] =  A_B0R4[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg3] =  A_B1R4[bn1_ma_reg8];
			 A_B1R4[bn1_ma_reg3]  =  data_tmp_1; 
			 A_B1R4[bn1_ma_reg4]  =  data_tmp_2; 
			 A_B0R4[bn0_ma_reg4]  =  data_tmp_3; 
			 A_B0R4[bn0_ma_reg5]  =  data_tmp_4; 
			 A_B1R4[bn1_ma_reg5]  =  data_tmp_5; 
			 A_B1R4[bn1_ma_reg6]  =  data_tmp_6; 
			 A_B0R4[bn0_ma_reg6]  =  data_tmp_7; 
			 A_B1R4[bn1_ma_reg7]  =  data_tmp_8; 
			 A_B0R4[bn0_ma_reg7]  =  data_tmp_9; 
			 A_B0R4[bn0_ma_reg8]  =  data_tmp_10;
			 A_B1R4[bn1_ma_reg8]  =  data_tmp_11;
			 //-------------------------------------------------------------
			 data_tmp_1  =  A_B1R6[bn1_ma_reg3];
			 data_tmp_2  =  A_B1R7[bn1_ma_reg3];
			 data_tmp_3  =  A_B1R8[bn1_ma_reg3];
			 data_tmp_4  =  A_B1R9[bn1_ma_reg3];
			 data_tmp_5  =  A_B1R10[bn1_ma_reg3];
			 data_tmp_6  =  A_B1R11[bn1_ma_reg3];
			 data_tmp_7  =  A_B1R12[bn1_ma_reg3];
			 data_tmp_8  =  A_B1R13[bn1_ma_reg3];
			 data_tmp_9  =  A_B1R14[bn1_ma_reg3];
			 data_tmp_10 =  A_B1R15[bn1_ma_reg3];
			 A_B1R6[bn1_ma_reg3]  =  A_B1R5[bn1_ma_reg4];
			 A_B1R7[bn1_ma_reg3]  =  A_B0R5[bn0_ma_reg4];
			 A_B1R8[bn1_ma_reg3]  =  A_B0R5[bn0_ma_reg5];
			 A_B1R9[bn1_ma_reg3]  =  A_B1R5[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg3] =  A_B1R5[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg3] =  A_B0R5[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg3] =  A_B1R5[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg3] =  A_B0R5[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg3] =  A_B0R5[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg3] =  A_B1R5[bn1_ma_reg8];
			 A_B1R5[bn1_ma_reg4]  =  data_tmp_1; 
			 A_B0R5[bn0_ma_reg4]  =  data_tmp_2; 
			 A_B0R5[bn0_ma_reg5]  =  data_tmp_3; 
			 A_B1R5[bn1_ma_reg5]  =  data_tmp_4; 
			 A_B1R5[bn1_ma_reg6]  =  data_tmp_5; 
			 A_B0R5[bn0_ma_reg6]  =  data_tmp_6; 
			 A_B1R5[bn1_ma_reg7]  =  data_tmp_7; 
			 A_B0R5[bn0_ma_reg7]  =  data_tmp_8; 
			 A_B0R5[bn0_ma_reg8]  =  data_tmp_9; 
			 A_B1R5[bn1_ma_reg8]  =  data_tmp_10;
			 //----------------------------------------------------------------
			 data_tmp_1  =  A_B1R7[bn1_ma_reg4];
			 data_tmp_2  =  A_B1R8[bn1_ma_reg4];
			 data_tmp_3  =  A_B1R9[bn1_ma_reg4];
			 data_tmp_4  =  A_B1R10[bn1_ma_reg4];
			 data_tmp_5  =  A_B1R11[bn1_ma_reg4];
			 data_tmp_6  =  A_B1R12[bn1_ma_reg4];
			 data_tmp_7  =  A_B1R13[bn1_ma_reg4];
			 data_tmp_8  =  A_B1R14[bn1_ma_reg4];
			 data_tmp_9  =  A_B1R15[bn1_ma_reg4];
			 A_B1R7[bn1_ma_reg4]  =  A_B0R6[bn0_ma_reg4];
			 A_B1R8[bn1_ma_reg4]  =  A_B0R6[bn0_ma_reg5];
			 A_B1R9[bn1_ma_reg4]  =  A_B1R6[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg4] =  A_B1R6[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg4] =  A_B0R6[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg4] =  A_B1R6[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg4] =  A_B0R6[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg4] =  A_B0R6[bn0_ma_reg8];
             A_B1R15[bn1_ma_reg4] =  A_B1R6[bn1_ma_reg8];
             A_B0R6[bn0_ma_reg4]  =  data_tmp_1;
             A_B0R6[bn0_ma_reg5]  =  data_tmp_2;
             A_B1R6[bn1_ma_reg5]  =  data_tmp_3;
             A_B1R6[bn1_ma_reg6]  =  data_tmp_4;
             A_B0R6[bn0_ma_reg6]  =  data_tmp_5;
             A_B1R6[bn1_ma_reg7]  =  data_tmp_6;
             A_B0R6[bn0_ma_reg7]  =  data_tmp_7;
             A_B0R6[bn0_ma_reg8]  =  data_tmp_8;
             A_B1R6[bn1_ma_reg8]  =  data_tmp_9;
             //------------------------------------------------------------------
			 data_tmp_1  = A_B0R8[bn0_ma_reg4];
			 data_tmp_2  = A_B0R9[bn0_ma_reg4];
			 data_tmp_3  = A_B0R10[bn0_ma_reg4];
			 data_tmp_4  = A_B0R11[bn0_ma_reg4];
			 data_tmp_5  = A_B0R12[bn0_ma_reg4];
			 data_tmp_6  = A_B0R13[bn0_ma_reg4];
			 data_tmp_7  = A_B0R14[bn0_ma_reg4];
			 data_tmp_8  = A_B0R15[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg4]  =  A_B0R7[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg4]  =  A_B1R7[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg4] =  A_B1R7[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg4] =  A_B0R7[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg4] =  A_B1R7[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg4] =  A_B0R7[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg4] =  A_B0R7[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg4] =  A_B1R7[bn1_ma_reg8];
			 A_B0R7[bn0_ma_reg5]  =  data_tmp_1;
			 A_B1R7[bn1_ma_reg5]  =  data_tmp_2;
			 A_B1R7[bn1_ma_reg6]  =  data_tmp_3;
			 A_B0R7[bn0_ma_reg6]  =  data_tmp_4;
			 A_B1R7[bn1_ma_reg7]  =  data_tmp_5;
			 A_B0R7[bn0_ma_reg7]  =  data_tmp_6;
			 A_B0R7[bn0_ma_reg8]  =  data_tmp_7;
			 A_B1R7[bn1_ma_reg8]  =  data_tmp_8;
			 //----------------------------------------------------------------
			 data_tmp_1  = A_B0R9[bn0_ma_reg5];
			 data_tmp_2  = A_B0R10[bn0_ma_reg5];
			 data_tmp_3  = A_B0R11[bn0_ma_reg5];
			 data_tmp_4  = A_B0R12[bn0_ma_reg5];
			 data_tmp_5  = A_B0R13[bn0_ma_reg5];
			 data_tmp_6  = A_B0R14[bn0_ma_reg5];
			 data_tmp_7  = A_B0R15[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg5]  =  A_B1R8[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg8];
			 A_B1R8[bn1_ma_reg5]  =  data_tmp_1;
			 A_B1R8[bn1_ma_reg6]  =  data_tmp_2;
			 A_B0R8[bn0_ma_reg6]  =  data_tmp_3;
			 A_B1R8[bn1_ma_reg7]  =  data_tmp_4;
			 A_B0R8[bn0_ma_reg7]  =  data_tmp_5;
			 A_B0R8[bn0_ma_reg8]  =  data_tmp_6;
			 A_B1R8[bn1_ma_reg8]  =  data_tmp_7;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B1R10[bn1_ma_reg5];
			 data_tmp_2 = A_B1R11[bn1_ma_reg5];
			 data_tmp_3 = A_B1R12[bn1_ma_reg5];
			 data_tmp_4 = A_B1R13[bn1_ma_reg5];
			 data_tmp_5 = A_B1R14[bn1_ma_reg5];
			 data_tmp_6 = A_B1R15[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg5] = A_B1R9[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg5] = A_B0R9[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg5] = A_B1R9[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg5] = A_B0R9[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg5] = A_B0R9[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg5] = A_B1R9[bn1_ma_reg8];
			 A_B1R9[bn1_ma_reg6]  = data_tmp_1;
			 A_B0R9[bn0_ma_reg6]  = data_tmp_2;
			 A_B1R9[bn1_ma_reg7]  = data_tmp_3;
			 A_B0R9[bn0_ma_reg7]  = data_tmp_4;
			 A_B0R9[bn0_ma_reg8]  = data_tmp_5;
			 A_B1R9[bn1_ma_reg8]  = data_tmp_6;
			 //----------------------------------------------------------------
			 data_tmp_1 =  A_B1R11[bn1_ma_reg6];
			 data_tmp_2 =  A_B1R12[bn1_ma_reg6];
			 data_tmp_3 =  A_B1R13[bn1_ma_reg6];
			 data_tmp_4 =  A_B1R14[bn1_ma_reg6];
			 data_tmp_5 =  A_B1R15[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg6]  = A_B0R10[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg6]  = A_B1R10[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg6]  = A_B0R10[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg6]  = A_B0R10[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg6]  = A_B1R10[bn1_ma_reg8];
			 A_B0R10[bn0_ma_reg6]  = data_tmp_1;
			 A_B1R10[bn1_ma_reg7]  = data_tmp_2;
			 A_B0R10[bn0_ma_reg7]  = data_tmp_3;
			 A_B0R10[bn0_ma_reg8]  = data_tmp_4;
			 A_B1R10[bn1_ma_reg8]  = data_tmp_5;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B0R12[bn0_ma_reg6];
			 data_tmp_2 = A_B0R13[bn0_ma_reg6];
			 data_tmp_3 = A_B0R14[bn0_ma_reg6];
			 data_tmp_4 = A_B0R15[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg6] = A_B1R11[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg6] = A_B0R11[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg6] = A_B0R11[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg6] = A_B1R11[bn1_ma_reg8];
			 A_B1R11[bn1_ma_reg7] = data_tmp_1;
			 A_B0R11[bn0_ma_reg7] = data_tmp_2;
			 A_B0R11[bn0_ma_reg8] = data_tmp_3;
			 A_B1R11[bn1_ma_reg8] = data_tmp_4;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B1R13[bn1_ma_reg7];
			 data_tmp_2 = A_B1R14[bn1_ma_reg7];
			 data_tmp_3 = A_B1R15[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg7] = A_B0R12[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg7] = A_B0R12[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg7] = A_B1R12[bn1_ma_reg8];
			 A_B0R12[bn0_ma_reg7] = data_tmp_1;
			 A_B0R12[bn0_ma_reg8] = data_tmp_2;
			 A_B1R12[bn1_ma_reg8] = data_tmp_3;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B0R14[bn0_ma_reg7];
			 data_tmp_2 = A_B0R15[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg7] = A_B0R13[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg7] = A_B1R13[bn1_ma_reg8];
			 A_B0R13[bn0_ma_reg8] = data_tmp_1;
			 A_B1R13[bn1_ma_reg8] = data_tmp_2;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B0R15[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg8] = A_B1R14[bn1_ma_reg8];
			 A_B1R14[bn1_ma_reg8] = data_tmp_1;
			 //-----------------------------------------------------------------
		  }
		 }
		 //std::cout <<"A_B0R0["<<bn0_ma_reg<<"]: "<<A_B0R0[bn0_ma_reg]<<"\n";
		 //std::cout <<"A_B0R1["<<bn0_ma_reg<<"]: "<<A_B0R1[bn0_ma_reg]<<"\n";
		 //std::cout <<"A_B1R0["<<bn1_ma_reg<<"]: "<<A_B1R0[bn1_ma_reg]<<"\n";
		 //std::cout <<"A_B1R1["<<bn1_ma_reg<<"]: "<<A_B1R1[bn1_ma_reg]<<"\n";
         //std::cout <<"------------------------------------------\n";		 
		}
	}
	
	int index0;
    int index1;
    int index2;
    int index3;
    int index4;
    int index5;
    int index6;
    int index7;
    int index8;
    int index9;
    int index10;
    int index11;
    int index12;
    int index13;
    int index14;
    int index15;
	//data output
	//bit reverse
	for(int i = 0;i < group;i++){
		for(int j = 0;j < radix;j++){
			gray_i  = Gray(i,group);
			BC_tmp  = j * group + gray_i;
			RR_R16(BC_tmp,Stage-1,BC);
			AGU_R16(BC,bn_tmp,ma_tmp);
			if(bn_tmp == 0){
			   BR_R16(16 * BC_tmp    ,index0);
			   BR_R16(16 * BC_tmp + 1,index1);
			   BR_R16(16 * BC_tmp + 2,index2);
			   BR_R16(16 * BC_tmp + 3,index3);
			   BR_R16(16 * BC_tmp + 4,index4);
			   BR_R16(16 * BC_tmp + 5,index5);
			   BR_R16(16 * BC_tmp + 6,index6);
			   BR_R16(16 * BC_tmp + 7,index7);
			   BR_R16(16 * BC_tmp + 8,index8);
			   BR_R16(16 * BC_tmp + 9,index9);
			   BR_R16(16 * BC_tmp + 10,index10);
			   BR_R16(16 * BC_tmp + 11,index11);
			   BR_R16(16 * BC_tmp + 12,index12);
			   BR_R16(16 * BC_tmp + 13,index13);
			   BR_R16(16 * BC_tmp + 14,index14);
			   BR_R16(16 * BC_tmp + 15,index15);
               A[index0] = A_B0R0[ma_tmp];
			   A[index1] = A_B0R1[ma_tmp];
			   A[index2] = A_B0R2[ma_tmp];
			   A[index3] = A_B0R3[ma_tmp];
			   A[index4] = A_B0R4[ma_tmp];
			   A[index5] = A_B0R5[ma_tmp];
			   A[index6] = A_B0R6[ma_tmp];
			   A[index7] = A_B0R7[ma_tmp];
			   A[index8] = A_B0R8[ma_tmp];
			   A[index9] = A_B0R9[ma_tmp];
			   A[index10] = A_B0R10[ma_tmp];
			   A[index11] = A_B0R11[ma_tmp];
			   A[index12] = A_B0R12[ma_tmp];
			   A[index13] = A_B0R13[ma_tmp];
			   A[index14] = A_B0R14[ma_tmp];
			   A[index15] = A_B0R15[ma_tmp];
			}
			else {
			   //std::cout <<"BC_tmp*2: "<< BC_tmp*2<<"\n";
			   //std::cout <<"Bit reverse: "<< index_r0<<"\n";
			   BR_R16(16 * BC_tmp    ,index0);
			   BR_R16(16 * BC_tmp + 1,index1);
			   BR_R16(16 * BC_tmp + 2,index2);
			   BR_R16(16 * BC_tmp + 3,index3);
			   BR_R16(16 * BC_tmp + 4,index4);
			   BR_R16(16 * BC_tmp + 5,index5);
			   BR_R16(16 * BC_tmp + 6,index6);
			   BR_R16(16 * BC_tmp + 7,index7);
			   BR_R16(16 * BC_tmp + 8,index8);
			   BR_R16(16 * BC_tmp + 9,index9);
			   BR_R16(16 * BC_tmp + 10,index10);
			   BR_R16(16 * BC_tmp + 11,index11);
			   BR_R16(16 * BC_tmp + 12,index12);
			   BR_R16(16 * BC_tmp + 13,index13);
			   BR_R16(16 * BC_tmp + 14,index14);
			   BR_R16(16 * BC_tmp + 15,index15);
               A[index0]     = A_B1R0[ma_tmp];
               A[index1]     = A_B1R1[ma_tmp];
               A[index2]     = A_B1R2[ma_tmp];
               A[index3]     = A_B1R3[ma_tmp];
               A[index4]     = A_B1R4[ma_tmp];
               A[index5]     = A_B1R5[ma_tmp];
               A[index6]     = A_B1R6[ma_tmp];
               A[index7]     = A_B1R7[ma_tmp];
               A[index8]     = A_B1R8[ma_tmp];
               A[index9]     = A_B1R9[ma_tmp];
               A[index10]    = A_B1R10[ma_tmp];
               A[index11]    = A_B1R11[ma_tmp];
               A[index12]    = A_B1R12[ma_tmp];
               A[index13]    = A_B1R13[ma_tmp];
               A[index14]    = A_B1R14[ma_tmp];
               A[index15]    = A_B1R15[ma_tmp];
			}			
		}		
	}
	
}
//----------------------------------------------
//mixed radix fft
//radix-4 and raidx-2
void NTTSPMB::NTT_r4_r2(std::vector<ZZ> &A,
std::vector<ZZ> &B0R0,std::vector<ZZ> &B0R1,std::vector<ZZ> &B0R2,std::vector<ZZ> &B0R3,
std::vector<ZZ> &B1R0,std::vector<ZZ> &B1R1,std::vector<ZZ> &B1R2,std::vector<ZZ> &B1R3){
	int            Stage;//stage	
	int            word_size;
	int            offset;
	int            BC;
	int            BC_tmp;
    int            group;
    int            ma_tmp;
	int            bn_tmp;
    int            bit_tmp;
	int            bn0_bc_tmp;//frist in bc data
	int            bn1_bc_tmp;//frist in bc data
	int            bn0_ma_reg1;
	int            bn0_ma_reg2;
	int            bn1_ma_reg1;
	int            bn1_ma_reg2;
	int            gray_i;
    int            BC_WIDTH;
    int 	       tw_modulus;
    int 	       tw_modulus_tmp;
    std::vector<int> bit_array_tmp;
	std::cout << "****************************************\n";
	//std::cout << "NTT_R4_R2 Start calculate!\n";
	
	std::ofstream DATARECORD("./NTT_R4_R2_SPMB.txt");
	std::ofstream spmb_r4_r2("./SPMB_tw/spmb_r4_r2.txt");
	std::ofstream DTFAG_golden_st0("./SPMB_tw/DTFAG_golden_st0_128.txt");
	std::ofstream DTFAG_golden_st1("./SPMB_tw/DTFAG_golden_st1_128.txt");
	std::ofstream DTFAG_golden_st2("./SPMB_tw/DTFAG_golden_st2_128.txt");
	//-------------Bit Reverse--------------
	BitOperate RR, IntToVec, VecToInt;
	long long bit_width = log2(N);
	long long bit_width_s = log(radix)/log(2);
	vector<long long> bit_array;
	int debug = 0;
	//----------------------------------------

    Stage = (unsigned long)ceil(log2(N));
	if((Stage % 2) == 1) Stage = Stage - 1;
	Stage = (unsigned long)Stage / 2;
	BC_WIDTH    = (int)ceil(log2(N/4));
	offset      =  (int)N / 4 ;
	word_size   =  (int)N / 8; // 2 * 4
	group       =  (int)N / 16; // 4 * 4
	tw_modulus  =  (int)N / 4 ;
	bit_array_tmp.resize(BC_WIDTH);
	DATARECORD <<"Stage: "<< Stage<<"\n";
	
	//SRAM
	ZZ               data_tmp;
	ZZ               data_tmp_1;
	ZZ               data_tmp_2;
	std::vector<ZZ>  A_B0R0;
	std::vector<ZZ>  A_B0R1;
	std::vector<ZZ>  A_B0R2;
	std::vector<ZZ>  A_B0R3;
	std::vector<ZZ>  A_B1R0;
	std::vector<ZZ>  A_B1R1;
	std::vector<ZZ>  A_B1R2;
	std::vector<ZZ>  A_B1R3;
	
	A_B0R0.resize(word_size);
	A_B0R1.resize(word_size);
	A_B0R2.resize(word_size);
	A_B0R3.resize(word_size);
	A_B1R0.resize(word_size);
	A_B1R1.resize(word_size);
	A_B1R2.resize(word_size);
	A_B1R3.resize(word_size);
	/*******************************************/
	//2020/08/12 modify
	/******************************************/
	B0R0.resize(word_size);
	B0R1.resize(word_size);
	B0R2.resize(word_size);
	B0R3.resize(word_size);
	B1R0.resize(word_size);
	B1R1.resize(word_size);
	B1R2.resize(word_size);
	B1R3.resize(word_size);
	//-----------------------------------------
	int length;
	ZZ  factor;   //base factor
	ZZ  factor_t; //acctually mul factor
	ZZ  factor_2t;
	ZZ  factor_3t;
    
	//init load data
    for(int i = 0; i < group; i++){
		for(int j = 0 ; j < radix ; j++){
			bn_tmp = 0;
			ma_tmp = 0;
			gray_i  = Gray(i,group);
			BC = j * group + gray_i;
            //bit calculate
            for(int j=0; j < BC_WIDTH;j++){
                bit_tmp = BC % 2;
                BC = BC >> 1;
                bit_array_tmp[j] = bit_tmp;
            } 
			gray_i  = Gray(i,group);
			BC = j * group + gray_i;
            for(int rs = 0; rs < BC_WIDTH; rs++){
                if((bit_array_tmp[rs] == 1) && (rs != 0)) ma_tmp = ma_tmp + exp2((rs-1)); 
                bn_tmp = bn_tmp ^ bit_array_tmp[rs];
            }
			if(bn_tmp == 0){
				A_B0R0[ma_tmp] = A[BC];
				A_B0R1[ma_tmp] = A[BC + offset];
				A_B0R2[ma_tmp] = A[BC + 2 * offset];
				A_B0R3[ma_tmp] = A[BC + 3 * offset];
			}else {
				A_B1R0[ma_tmp] = A[BC];
				A_B1R1[ma_tmp] = A[BC + offset];
				A_B1R2[ma_tmp] = A[BC + 2 * offset];
				A_B1R3[ma_tmp] = A[BC + 3 * offset];
			}
		}
	}

	ma_tmp = 0;
	bn_tmp = 0;
	BC     = 0;
	DATARECORD << "init load over! \n";
	int tw_degree = 1; // siang
	//need modify to mult by twiddle factor
	for(int s = 0; s < Stage; s++){
		DATARECORD <<"----------------------------"<<"\n";
		DATARECORD <<"Stage: "<< s<<"\n";
		std::cout <<"----------------------------"<<"\n";
		std::cout <<"Stage: "<< s<<"\n";
		spmb_r4_r2 <<"----------------------------"<<"\n";
		spmb_r4_r2 <<"Stage: "<< s<<"\n";
		if(s == 0)factor = W;
		else {
			SqrMod(factor,factor,p);
			SqrMod(factor,factor,p);
			tw_degree = tw_degree * 4;
		}
		tw_modulus_tmp = tw_modulus >> (2 * s);
		for(int i = 0 ;i < group;i++){
			bn0_bc_tmp  = 0;
			bn1_bc_tmp  = 0;
			spmb_r4_r2 <<"---------------i = " << i << "-----------------"<<"\n";
			for(int j = 0;j < radix;j++){
				DATARECORD <<"----------------------------"<<"\n";
				DATARECORD <<"i: "<< i << ", j = " << j <<"\n";
				std::cout <<"----------------------------"<<"\n";
				std::cout <<"i: "<< i << ", j = " << j <<"\n";
				DATARECORD <<"twiddle factor =  " << factor << std::endl;
				gray_i  = Gray(i,group);
				BC_tmp  = j * group + gray_i;
				if(s == Stage-1)RR_R4_R2(BC_tmp,(2*s - 1),BC);
				else RR_R4(BC_tmp,s,BC);
				length = BC % tw_modulus_tmp;
				PowerMod(factor_t,factor,length,p);
				AGU_R4(BC,bn_tmp,ma_tmp);

				//-----------compute data idx-------------
				spmb_r4_r2 /*<< "BC = " << BC */<< "BC_tmp = " << BC_tmp << std::endl;
				
				IntToVec.IntToVec(BC_tmp, N, bit_array);
				rotate(bit_array.begin(), bit_array.begin()+ bit_width_s*s , bit_array.end());
				long long Data = VecToInt.VecToInt(bit_array, N);
				spmb_r4_r2 << "Data_index = ";
                spmb_r4_r2 << "( " ;					
				for(int k = 0; k < radix ; k++ ){
					spmb_r4_r2 << Data + k*(1<<(bit_width-bit_width_s-bit_width_s*s)) <<" ";	
				}
				spmb_r4_r2 << ") " ;
				spmb_r4_r2 << ", (w^" << 0 << ", w^" << tw_degree * length 
				<< ", w^" << tw_degree * length * 2 << ", w^" << tw_degree * length * 3
				<< ")" << std::endl;
				//-----------------------------------------

				DATARECORD << "BC_tmp = " << BC_tmp << ", length: " << length <<
				", tw_dg: " <<  tw_degree * length << " = " << tw_degree << " * " << length <<"\n";	// TF_record
				//DATARECORD <<"length =  " << length <<"\n";
				DATARECORD <<"factor_t =  " << factor_t << ", p = " <<  p <<"\n";

				if(bn_tmp == 0){
					
					if(j < 2)bn0_bc_tmp = BC_tmp;
					PowerMod(factor_2t,factor_t,2,p);
					PowerMod(factor_3t,factor_t,3,p);
					DATARECORD << "*****before BU*****" << std::endl;
					DATARECORD <<"A_B0R0[" << ma_tmp << "] = " << A_B0R0[ma_tmp] << std::endl;
					DATARECORD <<"A_B0R1[" << ma_tmp << "] = " << A_B0R1[ma_tmp] << std::endl;
					DATARECORD <<"A_B0R2[" << ma_tmp << "] = " << A_B0R2[ma_tmp] << std::endl;
					DATARECORD <<"A_B0R3[" << ma_tmp << "] = " << A_B0R3[ma_tmp] << std::endl;
					if(!debug) Radix4_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp],A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
					DATARECORD << "*****after BU*****" << std::endl;
					DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp] << ", factor_0 : w^" << tw_degree*0<<"\n";
					DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp] << ", facotr_t : w^" << tw_degree*length  <<"\n";
					DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<<A_B0R2[ma_tmp] << ", facotr_2t : w^" << tw_degree*length*2 <<"\n";
					DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<<A_B0R3[ma_tmp] << ", facotr_3t : w^" << tw_degree*length*3 <<"\n";
					if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],factor_t,p);
					if(!debug) MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],factor_2t,p);
					if(!debug) MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],factor_3t,p);
					DATARECORD << "*****after Mul*****" << std::endl;
					DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp] << std::endl;
					DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp] << ", facotr_t : " << factor_t  << ", length = " << length * 1 << std::endl;
					DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp] << ", facotr_2t : " << factor_2t << ", length = " << length * 2 << std::endl;
					DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp] << ", facotr_3t : " << factor_3t << ", length = " << length * 3 << std::endl;

					if(j <  2)bn0_ma_reg1 = ma_tmp;
					if(j >= 2)bn0_ma_reg2 = ma_tmp;
					//----------------DTFAG golden pattern------------------
					switch(s){
						case 0:
							DTFAG_golden_st0 << 1 		  << " \n ";
							DTFAG_golden_st0 << factor_t  << " \n ";
							DTFAG_golden_st0 << factor_2t << " \n ";
							DTFAG_golden_st0 << factor_3t << " \n ";
							break;
						case 1:
							DTFAG_golden_st1 << 1 		  << " \n ";
							DTFAG_golden_st1 << factor_t  << " \n ";
							DTFAG_golden_st1 << factor_2t << " \n ";
							DTFAG_golden_st1 << factor_3t << " \n ";
							break;
						case 2:
							DTFAG_golden_st2 << 1 		  << " \n ";
							DTFAG_golden_st2 << factor_t  << " \n ";
							DTFAG_golden_st2 << factor_2t << " \n ";
							DTFAG_golden_st2 << factor_3t << " \n ";
							break;
						default:
							break;
					}
					//--------------------------------------------------
				}
			    else {
					if(j < 2)bn1_bc_tmp = BC_tmp;
					PowerMod(factor_2t,factor_t,2,p);
					PowerMod(factor_3t,factor_t,3,p);
					DATARECORD << "*****before BU*****" << std::endl;
					DATARECORD <<"A_B1R0[" << ma_tmp << "] = " << A_B1R0[ma_tmp] << std::endl;
					DATARECORD <<"A_B1R1[" << ma_tmp << "] = " << A_B1R1[ma_tmp] << std::endl;
					DATARECORD <<"A_B1R2[" << ma_tmp << "] = " << A_B1R2[ma_tmp] << std::endl;
					DATARECORD <<"A_B1R3[" << ma_tmp << "] = " << A_B1R3[ma_tmp] << std::endl;
					if(!debug) Radix4_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp],A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
					DATARECORD << "*****after BU*****" << std::endl;
					DATARECORD << "A_B1R0[" << ma_tmp << "] = " << A_B1R0[ma_tmp] << ", factor_0 : w^" << tw_degree*0<<"\n";
					DATARECORD << "A_B1R1[" << ma_tmp << "] = " << A_B1R1[ma_tmp] << ", facotr_t : w^" << tw_degree*length  <<"\n";
					DATARECORD << "A_B1R2[" << ma_tmp << "] = " << A_B1R2[ma_tmp] << ", facotr_2t : w^" << tw_degree*length*2 <<"\n";
					DATARECORD << "A_B1R3[" << ma_tmp << "] = " << A_B1R3[ma_tmp] << ", facotr_3t : w^" << tw_degree*length*3 <<"\n";
					if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],factor_t,p);
					if(!debug) MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],factor_2t,p);
					if(!debug) MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],factor_3t,p);
					DATARECORD << "*****after mul*****" << std::endl;
					DATARECORD << "A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp] << std::endl;
					DATARECORD << "A_B1R1[" << ma_tmp << "] = " << A_B1R1[ma_tmp] <<", facotr_t : " << factor_t  << ", length = " << length * 1 << std::endl;
					DATARECORD << "A_B1R2[" << ma_tmp << "] = " << A_B1R2[ma_tmp] <<", facotr_2t : " << factor_2t << ", length = " << length * 2 << std::endl;
					DATARECORD << "A_B1R3[" << ma_tmp << "] = " << A_B1R3[ma_tmp] <<", facotr_3t : " << factor_3t << ", length = " << length * 3 << std::endl;
					if(j < 2) bn1_ma_reg1 = ma_tmp;
					if(j >= 2)bn1_ma_reg2 = ma_tmp;

					//----------------DTFAG golden pattern------------------
					switch(s){
						case 0:
							DTFAG_golden_st0 << 1 		  << " \n ";
							DTFAG_golden_st0 << factor_t  << " \n ";
							DTFAG_golden_st0 << factor_2t << " \n ";
							DTFAG_golden_st0 << factor_3t << " \n ";
							break;
						case 1:
							DTFAG_golden_st1 << 1 		  << " \n ";
							DTFAG_golden_st1 << factor_t  << " \n ";
							DTFAG_golden_st1 << factor_2t << " \n ";
							DTFAG_golden_st1 << factor_3t << " \n ";
							break;
						case 2:
							DTFAG_golden_st2 << 1 		  << " \n ";
							DTFAG_golden_st2 << factor_t  << " \n ";
							DTFAG_golden_st2 << factor_2t << " \n ";
							DTFAG_golden_st2 << factor_3t << " \n ";
							break;
						default:
							break;
					}
					//--------------------------------------------------
				}
			}
		//data relocation
		if(s != Stage-1){	
		  if(bn1_bc_tmp > bn0_bc_tmp){
             if(bn0_ma_reg1 > bn0_ma_reg2){
				 ma_tmp = bn0_ma_reg1;
				 bn0_ma_reg1 = bn0_ma_reg2;
				 bn0_ma_reg2 = ma_tmp;
			 }
			 if(bn1_ma_reg1 > bn1_ma_reg2){
				 ma_tmp = bn1_ma_reg1;
				 bn1_ma_reg1 = bn1_ma_reg2;
				 bn1_ma_reg2 = ma_tmp;
			 } 
			 data_tmp   = A_B0R1[bn0_ma_reg1];
			 data_tmp_1 = A_B0R2[bn0_ma_reg1];
			 data_tmp_2 = A_B0R3[bn0_ma_reg1];
			 A_B0R1[bn0_ma_reg1] = A_B1R0[bn1_ma_reg1];
			 A_B0R2[bn0_ma_reg1] = A_B1R0[bn1_ma_reg2];
			 A_B0R3[bn0_ma_reg1] = A_B0R0[bn0_ma_reg2];
			 A_B1R0[bn1_ma_reg1] = data_tmp;
			 A_B1R0[bn1_ma_reg2] = data_tmp_1;
			 A_B0R0[bn0_ma_reg2] = data_tmp_2;
			 data_tmp   = A_B1R2[bn1_ma_reg1];
			 data_tmp_1 = A_B1R3[bn1_ma_reg1];
			 A_B1R2[bn1_ma_reg1] = A_B1R1[bn1_ma_reg2];
			 A_B1R3[bn1_ma_reg1] = A_B0R1[bn0_ma_reg2];
			 A_B1R1[bn1_ma_reg2] = data_tmp;
			 A_B0R1[bn0_ma_reg2] = data_tmp_1;
			 data_tmp = A_B1R3[bn1_ma_reg2];
			 A_B1R3[bn1_ma_reg2] = A_B0R2[bn0_ma_reg2];
			 A_B0R2[bn0_ma_reg2] = data_tmp;
		  }else {
             if(bn0_ma_reg1 > bn0_ma_reg2){
				 ma_tmp = bn0_ma_reg1;
				 bn0_ma_reg1 = bn0_ma_reg2;
				 bn0_ma_reg2 = ma_tmp;
			 }
			 if(bn1_ma_reg1 > bn1_ma_reg2){
				 ma_tmp = bn1_ma_reg1;
				 bn1_ma_reg1 = bn1_ma_reg2;
				 bn1_ma_reg2 = ma_tmp;
			 } 	 	 
		 	 data_tmp   = A_B1R1[bn1_ma_reg1];
		 	 data_tmp_1 = A_B1R2[bn1_ma_reg1];
		 	 data_tmp_2 = A_B1R3[bn1_ma_reg1];
			 A_B1R1[bn1_ma_reg1] = A_B0R0[bn0_ma_reg1];
			 A_B1R2[bn1_ma_reg1] = A_B0R0[bn0_ma_reg2];
			 A_B1R3[bn1_ma_reg1] = A_B1R0[bn1_ma_reg2];
			 A_B0R0[bn0_ma_reg1] = data_tmp;
			 A_B0R0[bn0_ma_reg2] = data_tmp_1;
			 A_B1R0[bn1_ma_reg2] = data_tmp_2;
			 data_tmp   = A_B0R2[bn0_ma_reg1];
			 data_tmp_1 = A_B0R3[bn0_ma_reg1];
			 A_B0R2[bn0_ma_reg1]  = A_B0R1[bn0_ma_reg2];
			 A_B0R3[bn0_ma_reg1]  = A_B1R1[bn1_ma_reg2];
			 A_B0R1[bn0_ma_reg2]  = data_tmp;
		 	 A_B1R1[bn1_ma_reg2]  = data_tmp_1;
			 data_tmp = A_B0R3[bn0_ma_reg2];
             A_B0R3[bn0_ma_reg2] = A_B1R2[bn1_ma_reg2];
			 A_B1R2[bn1_ma_reg2] = data_tmp;
			 
		  }
		}else {
		     if(bn1_bc_tmp > bn0_bc_tmp){
                if(bn0_ma_reg1 > bn0_ma_reg2){
			   	 ma_tmp = bn0_ma_reg1;
			   	 bn0_ma_reg1 = bn0_ma_reg2;
			   	 bn0_ma_reg2 = ma_tmp;
			    }
			    if(bn1_ma_reg1 > bn1_ma_reg2){
			   	 ma_tmp = bn1_ma_reg1;
			   	 bn1_ma_reg1 = bn1_ma_reg2;
			   	 bn1_ma_reg2 = ma_tmp;
			    } 
			    //bn0 r0 r2 r4 r6
			    //bn1 r1 r3 r5 r7
			    //relocation 
			    //bn0 r0 r1 r2 r3
			    //bn1 r4 r5 r6 r7
			    data_tmp   = A_B0R1[bn0_ma_reg1]; //r2
			    data_tmp_1 = A_B0R2[bn0_ma_reg1]; //r4
			    data_tmp_2 = A_B0R3[bn0_ma_reg1]; //r6
                A_B0R1[bn0_ma_reg1] = A_B1R0[bn1_ma_reg1]; // r1
			    A_B0R2[bn0_ma_reg1] = data_tmp; //r2
			    A_B0R3[bn0_ma_reg1] = A_B1R1[bn1_ma_reg1]; //r3
			    A_B1R0[bn1_ma_reg1] = data_tmp_1; //r4
                A_B1R1[bn1_ma_reg1] = A_B1R2[bn1_ma_reg1];
                A_B1R2[bn1_ma_reg1] = data_tmp_2;
			    //bn1 r8 r10 r12 r14
			    //bn0 r9 r11 r13 r15
			    //after 
			    //bn1 r8 r9 r10 r11
			    //bn0 r12 r13 r14 r15
                data_tmp   = A_B1R1[bn1_ma_reg2]; //r10
		        data_tmp_1 = A_B1R2[bn1_ma_reg2]; //r12
		        data_tmp_2 = A_B1R3[bn1_ma_reg2]; //r14
		        A_B1R1[bn1_ma_reg2] = A_B0R0[bn0_ma_reg2]; //r9
		        A_B1R2[bn1_ma_reg2] = data_tmp; //r10
		        A_B1R3[bn1_ma_reg2] = A_B0R1[bn0_ma_reg2];//r11
                A_B0R0[bn0_ma_reg2] = data_tmp_1;
                A_B0R1[bn0_ma_reg2] = A_B0R2[bn0_ma_reg2];
			    A_B0R2[bn0_ma_reg2] = data_tmp_2;
		     }else {
                if(bn0_ma_reg1 > bn0_ma_reg2){
			     ma_tmp = bn0_ma_reg1;
			     bn0_ma_reg1 = bn0_ma_reg2;
			     bn0_ma_reg2 = ma_tmp;
			    }
			    if(bn1_ma_reg1 > bn1_ma_reg2){
			     ma_tmp = bn1_ma_reg1;
			     bn1_ma_reg1 = bn1_ma_reg2;
			     bn1_ma_reg2 = ma_tmp;
			    } 				  
			    //bn1 r16 r18 r20 r22 
			    //bn0 r17 r19 r21 r23
                //after relocation
				//bn1 r16 r17 r18 r19
				//bn0 r20 r21 r22 r23
			     data_tmp   = A_B1R1[bn1_ma_reg1];
			     data_tmp_1 = A_B1R2[bn1_ma_reg1];
			     data_tmp_2 = A_B1R3[bn1_ma_reg1];
			     A_B1R1[bn1_ma_reg1] = A_B0R0[bn0_ma_reg1];
				 A_B1R2[bn1_ma_reg1] = data_tmp;
				 A_B1R3[bn1_ma_reg1] = A_B0R1[bn0_ma_reg1];
				 A_B0R0[bn0_ma_reg1] = data_tmp_1;
				 A_B0R1[bn0_ma_reg1] = A_B0R2[bn0_ma_reg1];
				 A_B0R2[bn0_ma_reg1] = data_tmp_2;
				 //bn0 r24 r26 r28 r30
				 //bn1 r25 r27 r29 r31
				 //after relocation
				 //bn0 r24 r25 r26 r27
				 //bn1 r28 r29 r30 r31
				 data_tmp   = A_B0R1[bn0_ma_reg2];
			     data_tmp_1 = A_B0R2[bn0_ma_reg2];
			     data_tmp_2 = A_B0R3[bn0_ma_reg2];				 
				 A_B0R1[bn0_ma_reg2] = A_B1R0[bn1_ma_reg2];
				 A_B0R2[bn0_ma_reg2] = data_tmp;
				 A_B0R3[bn0_ma_reg2] = A_B1R1[bn1_ma_reg2];				 
				 A_B1R0[bn1_ma_reg2] = data_tmp_1;
				 A_B1R1[bn1_ma_reg2] = A_B1R2[bn1_ma_reg2];
				 A_B1R2[bn1_ma_reg2] = data_tmp_2; 
		     }
		    } 
	    }
	}
	DATARECORD <<"------------radix_2 part----------------"<< std::endl;	
	// radix-2 FFT compute
	for(int i = 0; i < group ; i++){
		for(int j=0; j < radix; j++){
			DATARECORD <<"------------------"<< std::endl;	
			gray_i  = Gray(i,group);
			BC_tmp  = j * group + gray_i;
			BC = BC_tmp;
			AGU_R4(BC,bn_tmp,ma_tmp);
			if(bn_tmp == 0){
				if(j < 2)bn0_bc_tmp = BC_tmp;
				DATARECORD <<"A_B0R0[" << ma_tmp << "] = " << A_B0R0[ma_tmp] <<", A_B0R1[" << ma_tmp << "] = " << A_B0R1[ma_tmp] << std::endl;
				DATARECORD <<"A_B0R2[" << ma_tmp << "] = " << A_B0R2[ma_tmp] <<", A_B0R3[" << ma_tmp << "] = " << A_B0R3[ma_tmp] << std::endl;
				Radix2_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp]);
				Radix2_BU(A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
				if(j <  2)bn0_ma_reg1 = ma_tmp;
				if(j >= 2)bn0_ma_reg2 = ma_tmp;
			}else {
			    if(j < 2)bn1_bc_tmp = BC_tmp;
				DATARECORD <<"A_B1R0[" << ma_tmp << "] = " << A_B1R0[ma_tmp] <<", A_B1R1[" << ma_tmp << "] = " << A_B1R1[ma_tmp] << std::endl;
				DATARECORD <<"A_B1R2[" << ma_tmp << "] = " << A_B1R2[ma_tmp] <<", A_B1R3[" << ma_tmp << "] = " << A_B1R3[ma_tmp] << std::endl;
			    Radix2_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp]);
			    Radix2_BU(A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
				if(j < 2) bn1_ma_reg1 = ma_tmp;
				if(j >= 2)bn1_ma_reg2 = ma_tmp;
			}			
		}
	}
	
	//*****************************************************************
	//data relocation , let order for INTT 
	//*****************************************************************
    for(int i = 0; i < group;i ++){
		for(int j = 0 ; j < radix; j++){
		   gray_i = Gray(i,group);	
		   BC_tmp = j * group + gray_i;
		   RR_R4_R2(BC_tmp,(2*(Stage-1) - 1),BC);
		   AGU_R4(BC,bn_tmp,ma_tmp);
		   if(bn_tmp == 0){
		   	if(j < 2)bn0_bc_tmp = BC_tmp;
		   	if(j <  2)bn0_ma_reg1 = ma_tmp;
		   	if(j >= 2)bn0_ma_reg2 = ma_tmp;
		   }else {
		    if(j < 2)bn1_bc_tmp = BC_tmp;   
		   	if(j < 2) bn1_ma_reg1 = ma_tmp;
		   	if(j >= 2)bn1_ma_reg2 = ma_tmp;
		   }			
		}
    if(bn1_bc_tmp > bn0_bc_tmp){
        if(bn0_ma_reg1 > bn0_ma_reg2){
	   	 ma_tmp = bn0_ma_reg1;
	   	 bn0_ma_reg1 = bn0_ma_reg2;
	   	 bn0_ma_reg2 = ma_tmp;
	    }
	    if(bn1_ma_reg1 > bn1_ma_reg2){
	   	 ma_tmp = bn1_ma_reg1;
	   	 bn1_ma_reg1 = bn1_ma_reg2;
	   	 bn1_ma_reg2 = ma_tmp;
	    } 
		//-----------------------------------------
		//bn: 0-------->1
        data_tmp   = A_B0R1[bn0_ma_reg1];
        data_tmp_1 = A_B0R2[bn0_ma_reg1];
        data_tmp_2 = A_B0R3[bn0_ma_reg1];
		A_B0R1[bn0_ma_reg1] = A_B1R0[bn1_ma_reg1];
		A_B0R2[bn0_ma_reg1] = data_tmp;
		A_B0R3[bn0_ma_reg1] = A_B1R1[bn1_ma_reg1];
		A_B1R0[bn1_ma_reg1] = data_tmp_1;
		A_B1R1[bn1_ma_reg1] = A_B1R2[bn1_ma_reg1];
		A_B1R2[bn1_ma_reg1] = data_tmp_2;
		//--------------------------------------
		//bn: 1-----------> 0
		data_tmp   = A_B1R1[bn1_ma_reg2];
		data_tmp_1 = A_B1R2[bn1_ma_reg2];
		data_tmp_2 = A_B1R3[bn1_ma_reg2];
		A_B1R1[bn1_ma_reg2] = A_B0R0[bn0_ma_reg2];
		A_B1R2[bn1_ma_reg2] = data_tmp;
		A_B1R3[bn1_ma_reg2] = A_B0R1[bn0_ma_reg2];
		A_B0R0[bn0_ma_reg2] = data_tmp_1;
		A_B0R1[bn0_ma_reg2] = A_B0R2[bn0_ma_reg2];
		A_B0R2[bn0_ma_reg2] = data_tmp_2;
    }else {
        if(bn0_ma_reg1 > bn0_ma_reg2){
	     ma_tmp = bn0_ma_reg1;
	     bn0_ma_reg1 = bn0_ma_reg2;
	     bn0_ma_reg2 = ma_tmp;
	    }
	    if(bn1_ma_reg1 > bn1_ma_reg2){
	     ma_tmp = bn1_ma_reg1;
	     bn1_ma_reg1 = bn1_ma_reg2;
	     bn1_ma_reg2 = ma_tmp;
	    } 
        //-------------------------------------
        //bn: 1 ------> 0
		data_tmp   = A_B1R1[bn1_ma_reg1];
		data_tmp_1 = A_B1R2[bn1_ma_reg1];
		data_tmp_2 = A_B1R3[bn1_ma_reg1];
		A_B1R1[bn1_ma_reg1] = A_B0R0[bn0_ma_reg1];
		A_B1R2[bn1_ma_reg1] = data_tmp;
		A_B1R3[bn1_ma_reg1] = A_B0R1[bn0_ma_reg1];
		A_B0R0[bn0_ma_reg1] = data_tmp_1;
		A_B0R1[bn0_ma_reg1] = A_B0R2[bn0_ma_reg1];
		A_B0R2[bn0_ma_reg1] = data_tmp_2;
        //------------------------------------
		//bn: 0------->1
		data_tmp   = A_B0R1[bn0_ma_reg2];
		data_tmp_1 = A_B0R2[bn0_ma_reg2];
		data_tmp_2 = A_B0R3[bn0_ma_reg2];
		A_B0R1[bn0_ma_reg2] = A_B1R0[bn1_ma_reg2];
		A_B0R2[bn0_ma_reg2] = data_tmp;
		A_B0R3[bn0_ma_reg2] = A_B1R1[bn1_ma_reg2];
		A_B1R0[bn1_ma_reg2] = data_tmp_1;
		A_B1R1[bn1_ma_reg2] = A_B1R2[bn1_ma_reg2];
		A_B1R2[bn1_ma_reg2] = data_tmp_2;
		
     }		
	}

    DATARECORD << "---------------------------------------------------------\n";
	DATARECORD << "*********************************************************\n";
	DATARECORD << "                  FINAL!!!!                              \n";
	DATARECORD << "*********************************************************\n";
	DATARECORD << "---------------------------------------------------------\n";
	
    int BC_bit_size;
	BC_bit_size = (int)ceil(log2(N/4));
	DATARECORD << "BC_bit_size: "<< BC_bit_size <<"\n";
	int index_out;
	for(int i = 0; i < group; i++){
    	for(int j = 0;j < radix;j++){
			DATARECORD << "-------------------------" << std::endl;
			gray_i  = Gray(i,group);
    		BC_tmp  = j * group + gray_i;
			REORDERBC_R4_R2_OUT(BC_tmp,BC);
			BC = BC_tmp;
			DATARECORD << "BC_tmp: " << BC_tmp <<"\n";
			DATARECORD << "BC: " << BC <<"\n";
    		AGU_R4(BC,bn_tmp,ma_tmp);
			DATARECORD << "bn_tmp: " << bn_tmp <<"\n";
			DATARECORD << "ma_tmp: " << ma_tmp <<"\n";
			index_out = (int)i * ( N / group) + 4 * j ;
			DATARECORD << "index_out: " << index_out << std::endl;
    		if(bn_tmp == 0){
     		   A[index_out+0] = A_B0R0[ma_tmp];
    		   A[index_out+1] = A_B0R1[ma_tmp];
    		   A[index_out+2] = A_B0R2[ma_tmp];
    		   A[index_out+3] = A_B0R3[ma_tmp];
			   DATARECORD << "BN: 0 \n";
			   DATARECORD << "A["<< index_out      << "]: " << A[index_out+0] << " \n";
			   DATARECORD << "A["<< index_out + 1  << "]: " << A[index_out+1] << " \n";
			   DATARECORD << "A["<< index_out + 2  << "]: " << A[index_out+2] << " \n";
			   DATARECORD << "A["<< index_out + 3  << "]: " << A[index_out+3] << " \n";
    		}
    	    else {
    		   A[index_out+0] = A_B1R0[ma_tmp];
    		   A[index_out+1] = A_B1R1[ma_tmp];
    		   A[index_out+2] = A_B1R2[ma_tmp];
    		   A[index_out+3] = A_B1R3[ma_tmp];
			   DATARECORD << "BN: 1 \n";
			   DATARECORD << "A["<< index_out      << "]: " << A[index_out+0] << " \n";
			   DATARECORD << "A["<< index_out + 1  << "]: " << A[index_out+1] << " \n";
			   DATARECORD << "A["<< index_out + 2  << "]: " << A[index_out+2] << " \n";
			   DATARECORD << "A["<< index_out + 3  << "]: " << A[index_out+3] << " \n";
    		}
    	}	   
	}
	
	for(int i = 0;i < word_size;i++){
		B0R0[i] = A_B0R0[i];
        B0R1[i] = A_B0R1[i];
        B0R2[i] = A_B0R2[i];
        B0R3[i] = A_B0R3[i];
        B1R0[i] = A_B1R0[i];
        B1R1[i] = A_B1R1[i];
        B1R2[i] = A_B1R2[i];
        B1R3[i] = A_B1R3[i];
	}
	DATARECORD << "---------------------------------------------------------\n";
	DATARECORD << "----                  OVER!!!!!!!                   -----\n";
	DATARECORD << "---------------------------------------------------------\n";
	DATARECORD << "---------------------------------------------------------\n";
	DATARECORD.close();
}
void NTTSPMB::INTT_r4_r2(std::vector<ZZ> &A,
std::vector<ZZ> &B0R0,std::vector<ZZ> &B0R1,std::vector<ZZ> &B0R2,std::vector<ZZ> &B0R3,
std::vector<ZZ> &B1R0,std::vector<ZZ> &B1R1,std::vector<ZZ> &B1R2,std::vector<ZZ> &B1R3){
	int            Stage;//stage	
	int            word_size;
	int            offset;
	int            BC_RR;
	int            BC_Reorder; //butterfly counter Reorder
	int            BC_tmp;
    int            group;
    int            ma_tmp;
	int            bn_tmp;
    int            bit_tmp;
	int            bn0_bc_tmp;//frist in bc data
	int            bn1_bc_tmp;//frist in bc data
	int            bn0_ma_reg1;
	int            bn0_ma_reg2;
	int            bn1_ma_reg1;
	int            bn1_ma_reg2;
	int            gray_i;
    int            BC_WIDTH;
    int 	       tw_modulus;
    int 	       tw_modulus_tmp;
    std::vector<int> bit_array_tmp;
	std::cout << "****************************************\n";
	//std::cout << "NTT_R4_R2 Start calculate!\n";
	
	std::ofstream DATARECORD("./INTT_R4_R2_SPMB.txt");
    Stage = (unsigned long)ceil(log2(N));
	if((Stage % 2) == 1) Stage = Stage - 1;
	Stage = (unsigned long)Stage / 2;
	BC_WIDTH    = (int)ceil(log2(N/4));
	offset      =  (int)N / 4 ;
	word_size   =  (int)N / 8; // 2 * 4
	group       =  (int)N / 16; // 4 * 4
	tw_modulus  =  (int)N / 4 ;
	bit_array_tmp.resize(BC_WIDTH);
	DATARECORD <<"Stage: "<< Stage<<"\n";
	
	//SRAM
	ZZ               data_tmp;
	ZZ               data_tmp_1;
	ZZ               data_tmp_2;
	std::vector<ZZ>  A_B0R0;
	std::vector<ZZ>  A_B0R1;
	std::vector<ZZ>  A_B0R2;
	std::vector<ZZ>  A_B0R3;
	std::vector<ZZ>  A_B1R0;
	std::vector<ZZ>  A_B1R1;
	std::vector<ZZ>  A_B1R2;
	std::vector<ZZ>  A_B1R3;
	
	A_B0R0.resize(word_size);
	A_B0R1.resize(word_size);
	A_B0R2.resize(word_size);
	A_B0R3.resize(word_size);
	A_B1R0.resize(word_size);
	A_B1R1.resize(word_size);
	A_B1R2.resize(word_size);
	A_B1R3.resize(word_size);
	int length;
	ZZ  factor;   //base factor
	ZZ  factor_t; //acctually mul factor
	ZZ  factor_2t;
	ZZ  factor_3t;

	std::cout <<"----------------------------------\n";
	std::cout <<"loading start!!\n";
	for(int i=0; i < word_size;i++){
		//std::cout << "i: " << i <<"\n";
		A_B0R0[i] = B0R0[i];
		A_B0R1[i] = B0R1[i];
		A_B0R2[i] = B0R2[i];
		A_B0R3[i] = B0R3[i];
		A_B1R0[i] = B1R0[i];
		A_B1R1[i] = B1R1[i];
		A_B1R2[i] = B1R2[i];
		A_B1R3[i] = B1R3[i];
	}
	
	
	ma_tmp = 0;
	bn_tmp = 0;
	BC_RR  = 0;
	std::cout << "init load over! \n";
	//need modify to mult by twiddle factor
	for(int s = 0; s < Stage; s++){
		DATARECORD << "*************************************************************\n";
		DATARECORD << "Stage: "<< s <<"\n";
		if(s == 0)factor = IW;
		else {
			SqrMod(factor,factor,p);
			SqrMod(factor,factor,p);
		}
		tw_modulus_tmp = tw_modulus >> (2 * s);
		for(int i = 0 ;i < group;i++){
			bn0_bc_tmp  = 0;
			bn1_bc_tmp  = 0;
			for(int j = 0;j < radix;j++){
				DATARECORD << "---------------------------------------------------\n";
				DATARECORD << "Stage: "<< s <<"\n";
				gray_i  = Gray(i,group);
				BC_tmp  = j * group + gray_i;
				DATARECORD << "BC_tmp: "<< BC_tmp <<"\n";
				if(s == Stage-1)RR_R4_R2(BC_tmp,(2*s - 1),BC_RR);
				else RR_R4(BC_tmp,s,BC_RR);
				DATARECORD << "BC_RR: "<< BC_RR <<"\n";
				BC_IFFT_Reorder_R4_R2(BC_RR,BC_Reorder);
				DATARECORD << "BC_Reorder: "<< BC_Reorder <<"\n";
				length = BC_RR % tw_modulus_tmp;
				DATARECORD << "factor: "<< factor <<"\n";
				DATARECORD << "twiddle length: "<< length <<"\n";
				PowerMod(factor_t,factor,length,p);
				AGU_R4(BC_Reorder,bn_tmp,ma_tmp);
				if(s==0)DATARECORD << "bn_tmp: "<< bn_tmp <<"\n";
				if(s==0)DATARECORD << "ma_tmp: "<< ma_tmp <<"\n";
				if(bn_tmp == 0){
					if(j < 2)bn0_bc_tmp = BC_tmp;
					PowerMod(factor_2t,factor_t,2,p);
					PowerMod(factor_3t,factor_t,3,p);
					DATARECORD << "Before multiplication and Radix-4 BU operation!!!\n";
					DATARECORD << "A_B0R0["<< ma_tmp  << "]: " << A_B0R0[ma_tmp] << " \n";
					DATARECORD << "A_B0R1["<< ma_tmp  << "]: " << A_B0R1[ma_tmp] << " \n";
					DATARECORD << "A_B0R2["<< ma_tmp  << "]: " << A_B0R2[ma_tmp] << " \n";
					DATARECORD << "A_B0R3["<< ma_tmp  << "]: " << A_B0R3[ma_tmp] << " \n";
					Radix4_BU_INTT(A_B0R0[ma_tmp],A_B0R1[ma_tmp],A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
					DATARECORD << "After Radix-4 BU operation!!!\n";
					DATARECORD << "A_B0R0["<< ma_tmp  << "]: " << A_B0R0[ma_tmp] << " \n";
					DATARECORD << "A_B0R1["<< ma_tmp  << "]: " << A_B0R1[ma_tmp] << " \n";
					DATARECORD << "A_B0R2["<< ma_tmp  << "]: " << A_B0R2[ma_tmp] << " \n";
					DATARECORD << "A_B0R3["<< ma_tmp  << "]: " << A_B0R3[ma_tmp] << " \n";					
					MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],factor_t,p);
					MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],factor_2t,p);
					MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],factor_3t,p);
					DATARECORD << "factor_t:  " << factor_t  << " \n";
					DATARECORD << "factor_2t: " << factor_2t << " \n";
					DATARECORD << "factor_3t: " << factor_3t << " \n";
					DATARECORD << "After multiplication  operation!!!\n";
					DATARECORD << "A_B0R0["<< ma_tmp  << "]: " << A_B0R0[ma_tmp] << " \n";
					DATARECORD << "A_B0R1["<< ma_tmp  << "]: " << A_B0R1[ma_tmp] << " \n";
					DATARECORD << "A_B0R2["<< ma_tmp  << "]: " << A_B0R2[ma_tmp] << " \n";
					DATARECORD << "A_B0R3["<< ma_tmp  << "]: " << A_B0R3[ma_tmp] << " \n";
					DATARECORD << "---------------------------------------------------\n";
					if(j <  2)bn0_ma_reg1 = ma_tmp;
					if(j >= 2)bn0_ma_reg2 = ma_tmp;
				}
			    else {
					if(j < 2)bn1_bc_tmp = BC_tmp;
					PowerMod(factor_2t,factor_t,2,p);
					PowerMod(factor_3t,factor_t,3,p);
					DATARECORD << "Before multiplication and Radix-4 BU operation!!!\n";
					DATARECORD << "A_B1R0["<< ma_tmp  << "]: " << A_B1R0[ma_tmp] << " \n";
					DATARECORD << "A_B1R1["<< ma_tmp  << "]: " << A_B1R1[ma_tmp] << " \n";
					DATARECORD << "A_B1R2["<< ma_tmp  << "]: " << A_B1R2[ma_tmp] << " \n";
					DATARECORD << "A_B1R3["<< ma_tmp  << "]: " << A_B1R3[ma_tmp] << " \n";
					Radix4_BU_INTT(A_B1R0[ma_tmp],A_B1R1[ma_tmp],A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
					DATARECORD << "After Radix-4 BU operation!!!\n";
					DATARECORD << "A_B1R0["<< ma_tmp  << "]: " << A_B1R0[ma_tmp] << " \n";
					DATARECORD << "A_B1R1["<< ma_tmp  << "]: " << A_B1R1[ma_tmp] << " \n";
					DATARECORD << "A_B1R2["<< ma_tmp  << "]: " << A_B1R2[ma_tmp] << " \n";
					DATARECORD << "A_B1R3["<< ma_tmp  << "]: " << A_B1R3[ma_tmp] << " \n";	
					MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],factor_t,p);
					MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],factor_2t,p);
					MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],factor_3t,p);
					DATARECORD << "factor_t:  " << factor_t << " \n";
					DATARECORD << "factor_2t: " << factor_2t << " \n";
					DATARECORD << "factor_3t: " << factor_3t<< " \n";
					DATARECORD << "After multiplication and Radix-4 BU operation!!!\n";
					DATARECORD << "A_B1R0["<< ma_tmp  << "]: " << A_B1R0[ma_tmp] << " \n";
					DATARECORD << "A_B1R1["<< ma_tmp  << "]: " << A_B1R1[ma_tmp] << " \n";
					DATARECORD << "A_B1R2["<< ma_tmp  << "]: " << A_B1R2[ma_tmp] << " \n";
					DATARECORD << "A_B1R3["<< ma_tmp  << "]: " << A_B1R3[ma_tmp] << " \n";
					DATARECORD << "---------------------------------------------------\n";
					if(j < 2) bn1_ma_reg1 = ma_tmp;
					if(j >= 2)bn1_ma_reg2 = ma_tmp;
				}
			}
		//data relocation
		if(s != Stage-1){	
		  if(bn1_bc_tmp > bn0_bc_tmp){
			DATARECORD << "****************************************************\n";
			DATARECORD << "BN: 0 1 1 0 !!\n";
			DATARECORD << "Before radix-4 data relocation!!!\n";
			DATARECORD << "A_B0R0["<< bn0_ma_reg1  << "]: " << A_B0R0[bn0_ma_reg1] << " \n";
			DATARECORD << "A_B0R1["<< bn0_ma_reg1  << "]: " << A_B0R1[bn0_ma_reg1] << " \n";
			DATARECORD << "A_B0R2["<< bn0_ma_reg1  << "]: " << A_B0R2[bn0_ma_reg1] << " \n";
			DATARECORD << "A_B0R3["<< bn0_ma_reg1  << "]: " << A_B0R3[bn0_ma_reg1] << " \n";
			DATARECORD << "A_B0R0["<< bn0_ma_reg2  << "]: " << A_B0R0[bn0_ma_reg2] << " \n";
			DATARECORD << "A_B0R1["<< bn0_ma_reg2  << "]: " << A_B0R1[bn0_ma_reg2] << " \n";
			DATARECORD << "A_B0R2["<< bn0_ma_reg2  << "]: " << A_B0R2[bn0_ma_reg2] << " \n";
			DATARECORD << "A_B0R3["<< bn0_ma_reg2  << "]: " << A_B0R3[bn0_ma_reg2] << " \n";
			DATARECORD << "A_B1R0["<< bn1_ma_reg1  << "]: " << A_B1R0[bn1_ma_reg1] << " \n";
			DATARECORD << "A_B1R1["<< bn1_ma_reg1  << "]: " << A_B1R1[bn1_ma_reg1] << " \n";
			DATARECORD << "A_B1R2["<< bn1_ma_reg1  << "]: " << A_B1R2[bn1_ma_reg1] << " \n";
			DATARECORD << "A_B1R3["<< bn1_ma_reg1  << "]: " << A_B1R3[bn1_ma_reg1] << " \n";				
			DATARECORD << "A_B1R0["<< bn1_ma_reg2  << "]: " << A_B1R0[bn1_ma_reg2] << " \n";
			DATARECORD << "A_B1R1["<< bn1_ma_reg2  << "]: " << A_B1R1[bn1_ma_reg2] << " \n";
			DATARECORD << "A_B1R2["<< bn1_ma_reg2  << "]: " << A_B1R2[bn1_ma_reg2] << " \n";
			DATARECORD << "A_B1R3["<< bn1_ma_reg2  << "]: " << A_B1R3[bn1_ma_reg2] << " \n";

			 data_tmp   = A_B0R1[bn0_ma_reg1];
			 data_tmp_1 = A_B0R2[bn0_ma_reg1];
			 data_tmp_2 = A_B0R3[bn0_ma_reg1];
			 A_B0R1[bn0_ma_reg1] = A_B1R0[bn1_ma_reg1];
			 A_B0R2[bn0_ma_reg1] = A_B1R0[bn1_ma_reg2];
			 A_B0R3[bn0_ma_reg1] = A_B0R0[bn0_ma_reg2];
			 A_B1R0[bn1_ma_reg1] = data_tmp;
			 A_B1R0[bn1_ma_reg2] = data_tmp_1;
			 A_B0R0[bn0_ma_reg2] = data_tmp_2;
			 data_tmp   = A_B1R2[bn1_ma_reg1];
			 data_tmp_1 = A_B1R3[bn1_ma_reg1];
			 A_B1R2[bn1_ma_reg1] = A_B1R1[bn1_ma_reg2];
			 A_B1R3[bn1_ma_reg1] = A_B0R1[bn0_ma_reg2];
			 A_B1R1[bn1_ma_reg2] = data_tmp;
			 A_B0R1[bn0_ma_reg2] = data_tmp_1;
			 data_tmp = A_B1R3[bn1_ma_reg2];
			 A_B1R3[bn1_ma_reg2] = A_B0R2[bn0_ma_reg2];
			 A_B0R2[bn0_ma_reg2] = data_tmp;
			 DATARECORD << "after radix-4 data relocation!!!\n";
			 DATARECORD << "A_B0R0["<< bn0_ma_reg1  << "]: " << A_B0R0[bn0_ma_reg1] << " \n";
			 DATARECORD << "A_B0R1["<< bn0_ma_reg1  << "]: " << A_B0R1[bn0_ma_reg1] << " \n";
			 DATARECORD << "A_B0R2["<< bn0_ma_reg1  << "]: " << A_B0R2[bn0_ma_reg1] << " \n";
			 DATARECORD << "A_B0R3["<< bn0_ma_reg1  << "]: " << A_B0R3[bn0_ma_reg1] << " \n";
			 DATARECORD << "A_B0R0["<< bn0_ma_reg2  << "]: " << A_B0R0[bn0_ma_reg2] << " \n";
			 DATARECORD << "A_B0R1["<< bn0_ma_reg2  << "]: " << A_B0R1[bn0_ma_reg2] << " \n";
			 DATARECORD << "A_B0R2["<< bn0_ma_reg2  << "]: " << A_B0R2[bn0_ma_reg2] << " \n";
			 DATARECORD << "A_B0R3["<< bn0_ma_reg2  << "]: " << A_B0R3[bn0_ma_reg2] << " \n";
			 DATARECORD << "A_B1R0["<< bn1_ma_reg1  << "]: " << A_B1R0[bn1_ma_reg1] << " \n";
			 DATARECORD << "A_B1R1["<< bn1_ma_reg1  << "]: " << A_B1R1[bn1_ma_reg1] << " \n";
			 DATARECORD << "A_B1R2["<< bn1_ma_reg1  << "]: " << A_B1R2[bn1_ma_reg1] << " \n";
			 DATARECORD << "A_B1R3["<< bn1_ma_reg1  << "]: " << A_B1R3[bn1_ma_reg1] << " \n";				
			 DATARECORD << "A_B1R0["<< bn1_ma_reg2  << "]: " << A_B1R0[bn1_ma_reg2] << " \n";
			 DATARECORD << "A_B1R1["<< bn1_ma_reg2  << "]: " << A_B1R1[bn1_ma_reg2] << " \n";
			 DATARECORD << "A_B1R2["<< bn1_ma_reg2  << "]: " << A_B1R2[bn1_ma_reg2] << " \n";
			 DATARECORD << "A_B1R3["<< bn1_ma_reg2  << "]: " << A_B1R3[bn1_ma_reg2] << " \n";
			 DATARECORD << "****************************************************\n";
			 
			 
		  }else {
			DATARECORD << "****************************************************\n";
			DATARECORD << "BN: 1 0 0 1 !!\n";
			DATARECORD << "Before radix-4 data relocation!!!\n";
			DATARECORD << "A_B0R0["<< bn0_ma_reg1  << "]: " << A_B0R0[bn0_ma_reg1] << " \n";
			DATARECORD << "A_B0R1["<< bn0_ma_reg1  << "]: " << A_B0R1[bn0_ma_reg1] << " \n";
			DATARECORD << "A_B0R2["<< bn0_ma_reg1  << "]: " << A_B0R2[bn0_ma_reg1] << " \n";
			DATARECORD << "A_B0R3["<< bn0_ma_reg1  << "]: " << A_B0R3[bn0_ma_reg1] << " \n";
			DATARECORD << "A_B0R0["<< bn0_ma_reg2  << "]: " << A_B0R0[bn0_ma_reg2] << " \n";
			DATARECORD << "A_B0R1["<< bn0_ma_reg2  << "]: " << A_B0R1[bn0_ma_reg2] << " \n";
			DATARECORD << "A_B0R2["<< bn0_ma_reg2  << "]: " << A_B0R2[bn0_ma_reg2] << " \n";
			DATARECORD << "A_B0R3["<< bn0_ma_reg2  << "]: " << A_B0R3[bn0_ma_reg2] << " \n";
			DATARECORD << "A_B1R0["<< bn1_ma_reg1  << "]: " << A_B1R0[bn1_ma_reg1] << " \n";
			DATARECORD << "A_B1R1["<< bn1_ma_reg1  << "]: " << A_B1R1[bn1_ma_reg1] << " \n";
			DATARECORD << "A_B1R2["<< bn1_ma_reg1  << "]: " << A_B1R2[bn1_ma_reg1] << " \n";
			DATARECORD << "A_B1R3["<< bn1_ma_reg1  << "]: " << A_B1R3[bn1_ma_reg1] << " \n";				
			DATARECORD << "A_B1R0["<< bn1_ma_reg2  << "]: " << A_B1R0[bn1_ma_reg2] << " \n";
			DATARECORD << "A_B1R1["<< bn1_ma_reg2  << "]: " << A_B1R1[bn1_ma_reg2] << " \n";
			DATARECORD << "A_B1R2["<< bn1_ma_reg2  << "]: " << A_B1R2[bn1_ma_reg2] << " \n";
			DATARECORD << "A_B1R3["<< bn1_ma_reg2  << "]: " << A_B1R3[bn1_ma_reg2] << " \n";
			 
		 	 data_tmp   = A_B1R1[bn1_ma_reg1];
		 	 data_tmp_1 = A_B1R2[bn1_ma_reg1];
		 	 data_tmp_2 = A_B1R3[bn1_ma_reg1];
			 A_B1R1[bn1_ma_reg1] = A_B0R0[bn0_ma_reg1];
			 A_B1R2[bn1_ma_reg1] = A_B0R0[bn0_ma_reg2];
			 A_B1R3[bn1_ma_reg1] = A_B1R0[bn1_ma_reg2];
			 A_B0R0[bn0_ma_reg1] = data_tmp;
			 A_B0R0[bn0_ma_reg2] = data_tmp_1;
			 A_B1R0[bn1_ma_reg2] = data_tmp_2;
			 data_tmp   = A_B0R2[bn0_ma_reg1];
			 data_tmp_1 = A_B0R3[bn0_ma_reg1];
			 A_B0R2[bn0_ma_reg1]  = A_B0R1[bn0_ma_reg2];
			 A_B0R3[bn0_ma_reg1]  = A_B1R1[bn1_ma_reg2];
			 A_B0R1[bn0_ma_reg2]  = data_tmp;
		 	 A_B1R1[bn1_ma_reg2]  = data_tmp_1;
			 data_tmp = A_B0R3[bn0_ma_reg2];
             A_B0R3[bn0_ma_reg2] = A_B1R2[bn1_ma_reg2];
			 A_B1R2[bn1_ma_reg2] = data_tmp;
			 DATARECORD << "after radix-4 data relocation!!!\n";
			 DATARECORD << "A_B0R0["<< bn0_ma_reg1  << "]: " << A_B0R0[bn0_ma_reg1] << " \n";
			 DATARECORD << "A_B0R1["<< bn0_ma_reg1  << "]: " << A_B0R1[bn0_ma_reg1] << " \n";
			 DATARECORD << "A_B0R2["<< bn0_ma_reg1  << "]: " << A_B0R2[bn0_ma_reg1] << " \n";
			 DATARECORD << "A_B0R3["<< bn0_ma_reg1  << "]: " << A_B0R3[bn0_ma_reg1] << " \n";
			 DATARECORD << "A_B0R0["<< bn0_ma_reg2  << "]: " << A_B0R0[bn0_ma_reg2] << " \n";
			 DATARECORD << "A_B0R1["<< bn0_ma_reg2  << "]: " << A_B0R1[bn0_ma_reg2] << " \n";
			 DATARECORD << "A_B0R2["<< bn0_ma_reg2  << "]: " << A_B0R2[bn0_ma_reg2] << " \n";
			 DATARECORD << "A_B0R3["<< bn0_ma_reg2  << "]: " << A_B0R3[bn0_ma_reg2] << " \n";
			 DATARECORD << "A_B1R0["<< bn1_ma_reg1  << "]: " << A_B1R0[bn1_ma_reg1] << " \n";
			 DATARECORD << "A_B1R1["<< bn1_ma_reg1  << "]: " << A_B1R1[bn1_ma_reg1] << " \n";
			 DATARECORD << "A_B1R2["<< bn1_ma_reg1  << "]: " << A_B1R2[bn1_ma_reg1] << " \n";
			 DATARECORD << "A_B1R3["<< bn1_ma_reg1  << "]: " << A_B1R3[bn1_ma_reg1] << " \n";				
			 DATARECORD << "A_B1R0["<< bn1_ma_reg2  << "]: " << A_B1R0[bn1_ma_reg2] << " \n";
			 DATARECORD << "A_B1R1["<< bn1_ma_reg2  << "]: " << A_B1R1[bn1_ma_reg2] << " \n";
			 DATARECORD << "A_B1R2["<< bn1_ma_reg2  << "]: " << A_B1R2[bn1_ma_reg2] << " \n";
			 DATARECORD << "A_B1R3["<< bn1_ma_reg2  << "]: " << A_B1R3[bn1_ma_reg2] << " \n";
			 DATARECORD << "****************************************************\n";
		  }
		}else {
		     if(bn1_bc_tmp > bn0_bc_tmp){
                if(bn0_ma_reg1 > bn0_ma_reg2){
			   	 ma_tmp = bn0_ma_reg1;
			   	 bn0_ma_reg1 = bn0_ma_reg2;
			   	 bn0_ma_reg2 = ma_tmp;
			    }
			    if(bn1_ma_reg1 > bn1_ma_reg2){
			   	 ma_tmp = bn1_ma_reg1;
			   	 bn1_ma_reg1 = bn1_ma_reg2;
			   	 bn1_ma_reg2 = ma_tmp;
			    } 
			    //bn0 r0 r2 r4 r6
			    //bn1 r1 r3 r5 r7
			    //relocation 
			    //bn0 r0 r1 r2 r3
			    //bn1 r4 r5 r6 r7
				//-------------------------------------------------------------------
				DATARECORD << "****************************************************\n";
				DATARECORD << "Before radix-2 data relocation!!!\n";
				DATARECORD << "A_B0R0["<< bn0_ma_reg1  << "]: " << A_B0R0[bn0_ma_reg1] << " \n";
				DATARECORD << "A_B0R1["<< bn0_ma_reg1  << "]: " << A_B0R1[bn0_ma_reg1] << " \n";
				DATARECORD << "A_B0R2["<< bn0_ma_reg1  << "]: " << A_B0R2[bn0_ma_reg1] << " \n";
				DATARECORD << "A_B0R3["<< bn0_ma_reg1  << "]: " << A_B0R3[bn0_ma_reg1] << " \n";
				DATARECORD << "A_B0R0["<< bn0_ma_reg2  << "]: " << A_B0R0[bn0_ma_reg2] << " \n";
				DATARECORD << "A_B0R1["<< bn0_ma_reg2  << "]: " << A_B0R1[bn0_ma_reg2] << " \n";
				DATARECORD << "A_B0R2["<< bn0_ma_reg2  << "]: " << A_B0R2[bn0_ma_reg2] << " \n";
				DATARECORD << "A_B0R3["<< bn0_ma_reg2  << "]: " << A_B0R3[bn0_ma_reg2] << " \n";
				DATARECORD << "A_B1R0["<< bn1_ma_reg1  << "]: " << A_B1R0[bn1_ma_reg1] << " \n";
				DATARECORD << "A_B1R1["<< bn1_ma_reg1  << "]: " << A_B1R1[bn1_ma_reg1] << " \n";
				DATARECORD << "A_B1R2["<< bn1_ma_reg1  << "]: " << A_B1R2[bn1_ma_reg1] << " \n";
				DATARECORD << "A_B1R3["<< bn1_ma_reg1  << "]: " << A_B1R3[bn1_ma_reg1] << " \n";				
				DATARECORD << "A_B1R0["<< bn1_ma_reg2  << "]: " << A_B1R0[bn1_ma_reg2] << " \n";
				DATARECORD << "A_B1R1["<< bn1_ma_reg2  << "]: " << A_B1R1[bn1_ma_reg2] << " \n";
				DATARECORD << "A_B1R2["<< bn1_ma_reg2  << "]: " << A_B1R2[bn1_ma_reg2] << " \n";
				DATARECORD << "A_B1R3["<< bn1_ma_reg2  << "]: " << A_B1R3[bn1_ma_reg2] << " \n";				
			
			    data_tmp   = A_B0R1[bn0_ma_reg1]; //r2
			    data_tmp_1 = A_B0R2[bn0_ma_reg1]; //r4
			    data_tmp_2 = A_B0R3[bn0_ma_reg1]; //r6
                A_B0R1[bn0_ma_reg1] = A_B1R0[bn1_ma_reg1]; // r1
			    A_B0R2[bn0_ma_reg1] = data_tmp; //r2
			    A_B0R3[bn0_ma_reg1] = A_B1R1[bn1_ma_reg1]; //r3
			    A_B1R0[bn1_ma_reg1] = data_tmp_1; //r4
                A_B1R1[bn1_ma_reg1] = A_B1R2[bn1_ma_reg1];
                A_B1R2[bn1_ma_reg1] = data_tmp_2;
			    //bn1 r8 r10 r12 r14
			    //bn0 r9 r11 r13 r15
			    //after 
			    //bn1 r8 r9 r10 r11
			    //bn0 r12 r13 r14 r15
                data_tmp   = A_B1R1[bn1_ma_reg2]; //r10
		        data_tmp_1 = A_B1R2[bn1_ma_reg2]; //r12
		        data_tmp_2 = A_B1R3[bn1_ma_reg2]; //r14
		        A_B1R1[bn1_ma_reg2] = A_B0R0[bn0_ma_reg2]; //r9
		        A_B1R2[bn1_ma_reg2] = data_tmp; //r10
		        A_B1R3[bn1_ma_reg2] = A_B0R1[bn0_ma_reg2];//r11
                A_B0R0[bn0_ma_reg2] = data_tmp_1;
                A_B0R1[bn0_ma_reg2] = A_B0R2[bn0_ma_reg2];
			    A_B0R2[bn0_ma_reg2] = data_tmp_2;
				
				DATARECORD << "After radix-2 data relocation!!!!!                 \n";
                DATARECORD << "A_B0R0["<< bn0_ma_reg1  << "]: " << A_B0R0[bn0_ma_reg1] << " \n";
                DATARECORD << "A_B0R1["<< bn0_ma_reg1  << "]: " << A_B0R1[bn0_ma_reg1] << " \n";
                DATARECORD << "A_B0R2["<< bn0_ma_reg1  << "]: " << A_B0R2[bn0_ma_reg1] << " \n";
                DATARECORD << "A_B0R3["<< bn0_ma_reg1  << "]: " << A_B0R3[bn0_ma_reg1] << " \n";
                DATARECORD << "A_B0R0["<< bn0_ma_reg2  << "]: " << A_B0R0[bn0_ma_reg2] << " \n";
                DATARECORD << "A_B0R1["<< bn0_ma_reg2  << "]: " << A_B0R1[bn0_ma_reg2] << " \n";
                DATARECORD << "A_B0R2["<< bn0_ma_reg2  << "]: " << A_B0R2[bn0_ma_reg2] << " \n";
                DATARECORD << "A_B0R3["<< bn0_ma_reg2  << "]: " << A_B0R3[bn0_ma_reg2] << " \n";
                DATARECORD << "A_B1R0["<< bn1_ma_reg1  << "]: " << A_B1R0[bn1_ma_reg1] << " \n";
                DATARECORD << "A_B1R1["<< bn1_ma_reg1  << "]: " << A_B1R1[bn1_ma_reg1] << " \n";
                DATARECORD << "A_B1R2["<< bn1_ma_reg1  << "]: " << A_B1R2[bn1_ma_reg1] << " \n";
                DATARECORD << "A_B1R3["<< bn1_ma_reg1  << "]: " << A_B1R3[bn1_ma_reg1] << " \n";				
                DATARECORD << "A_B1R0["<< bn1_ma_reg2  << "]: " << A_B1R0[bn1_ma_reg2] << " \n";
                DATARECORD << "A_B1R1["<< bn1_ma_reg2  << "]: " << A_B1R1[bn1_ma_reg2] << " \n";
                DATARECORD << "A_B1R2["<< bn1_ma_reg2  << "]: " << A_B1R2[bn1_ma_reg2] << " \n";
                DATARECORD << "A_B1R3["<< bn1_ma_reg2  << "]: " << A_B1R3[bn1_ma_reg2] << " \n";		
		     }else {
                if(bn0_ma_reg1 > bn0_ma_reg2){
			     ma_tmp = bn0_ma_reg1;
			     bn0_ma_reg1 = bn0_ma_reg2;
			     bn0_ma_reg2 = ma_tmp;
			    }
			    if(bn1_ma_reg1 > bn1_ma_reg2){
			     ma_tmp = bn1_ma_reg1;
			     bn1_ma_reg1 = bn1_ma_reg2;
			     bn1_ma_reg2 = ma_tmp;
			    } 				  
			    //bn1 r16 r18 r20 r22 
			    //bn0 r17 r19 r21 r23
                //after relocation
				//bn1 r16 r17 r18 r19
				//bn0 r20 r21 r22 r23
				//-------------------------------------------------------------------
				DATARECORD << "****************************************************\n";
				DATARECORD << "Before radix-2 data relocation!!!\n";
				DATARECORD << "A_B0R0["<< bn0_ma_reg1  << "]: " << A_B0R0[bn0_ma_reg1] << " \n";
				DATARECORD << "A_B0R1["<< bn0_ma_reg1  << "]: " << A_B0R1[bn0_ma_reg1] << " \n";
				DATARECORD << "A_B0R2["<< bn0_ma_reg1  << "]: " << A_B0R2[bn0_ma_reg1] << " \n";
				DATARECORD << "A_B0R3["<< bn0_ma_reg1  << "]: " << A_B0R3[bn0_ma_reg1] << " \n";
				DATARECORD << "A_B0R0["<< bn0_ma_reg2  << "]: " << A_B0R0[bn0_ma_reg2] << " \n";
				DATARECORD << "A_B0R1["<< bn0_ma_reg2  << "]: " << A_B0R1[bn0_ma_reg2] << " \n";
				DATARECORD << "A_B0R2["<< bn0_ma_reg2  << "]: " << A_B0R2[bn0_ma_reg2] << " \n";
				DATARECORD << "A_B0R3["<< bn0_ma_reg2  << "]: " << A_B0R3[bn0_ma_reg2] << " \n";
				DATARECORD << "A_B1R0["<< bn1_ma_reg1  << "]: " << A_B1R0[bn1_ma_reg1] << " \n";
				DATARECORD << "A_B1R1["<< bn1_ma_reg1  << "]: " << A_B1R1[bn1_ma_reg1] << " \n";
				DATARECORD << "A_B1R2["<< bn1_ma_reg1  << "]: " << A_B1R2[bn1_ma_reg1] << " \n";
				DATARECORD << "A_B1R3["<< bn1_ma_reg1  << "]: " << A_B1R3[bn1_ma_reg1] << " \n";				
				DATARECORD << "A_B1R0["<< bn1_ma_reg2  << "]: " << A_B1R0[bn1_ma_reg2] << " \n";
				DATARECORD << "A_B1R1["<< bn1_ma_reg2  << "]: " << A_B1R1[bn1_ma_reg2] << " \n";
				DATARECORD << "A_B1R2["<< bn1_ma_reg2  << "]: " << A_B1R2[bn1_ma_reg2] << " \n";
				DATARECORD << "A_B1R3["<< bn1_ma_reg2  << "]: " << A_B1R3[bn1_ma_reg2] << " \n";				
				
			     data_tmp   = A_B1R1[bn1_ma_reg1];
			     data_tmp_1 = A_B1R2[bn1_ma_reg1];
			     data_tmp_2 = A_B1R3[bn1_ma_reg1];
			     A_B1R1[bn1_ma_reg1] = A_B0R0[bn0_ma_reg1];
				 A_B1R2[bn1_ma_reg1] = data_tmp;
				 A_B1R3[bn1_ma_reg1] = A_B0R1[bn0_ma_reg1];
				 A_B0R0[bn0_ma_reg1] = data_tmp_1;
				 A_B0R1[bn0_ma_reg1] = A_B0R2[bn0_ma_reg1];
				 A_B0R2[bn0_ma_reg1] = data_tmp_2;
				 //bn0 r24 r26 r28 r30
				 //bn1 r25 r27 r29 r31
				 //after relocation
				 //bn0 r24 r25 r26 r27
				 //bn1 r28 r29 r30 r31
				 data_tmp   = A_B0R1[bn0_ma_reg2];
			     data_tmp_1 = A_B0R2[bn0_ma_reg2];
			     data_tmp_2 = A_B0R3[bn0_ma_reg2];				 
				 A_B0R1[bn0_ma_reg2] = A_B1R0[bn1_ma_reg2];
				 A_B0R2[bn0_ma_reg2] = data_tmp;
				 A_B0R3[bn0_ma_reg2] = A_B1R1[bn1_ma_reg2];				 
				 A_B1R0[bn1_ma_reg2] = data_tmp_1;
				 A_B1R1[bn1_ma_reg2] = A_B1R2[bn1_ma_reg2];
				 A_B1R2[bn1_ma_reg2] = data_tmp_2; 
				 
				 DATARECORD << "After radix-2 data relocation!!!!!                 \n";
                 DATARECORD << "A_B0R0["<< bn0_ma_reg1  << "]: " << A_B0R0[bn0_ma_reg1] << " \n";
                 DATARECORD << "A_B0R1["<< bn0_ma_reg1  << "]: " << A_B0R1[bn0_ma_reg1] << " \n";
                 DATARECORD << "A_B0R2["<< bn0_ma_reg1  << "]: " << A_B0R2[bn0_ma_reg1] << " \n";
                 DATARECORD << "A_B0R3["<< bn0_ma_reg1  << "]: " << A_B0R3[bn0_ma_reg1] << " \n";
                 DATARECORD << "A_B0R0["<< bn0_ma_reg2  << "]: " << A_B0R0[bn0_ma_reg2] << " \n";
                 DATARECORD << "A_B0R1["<< bn0_ma_reg2  << "]: " << A_B0R1[bn0_ma_reg2] << " \n";
                 DATARECORD << "A_B0R2["<< bn0_ma_reg2  << "]: " << A_B0R2[bn0_ma_reg2] << " \n";
                 DATARECORD << "A_B0R3["<< bn0_ma_reg2  << "]: " << A_B0R3[bn0_ma_reg2] << " \n";
                 DATARECORD << "A_B1R0["<< bn1_ma_reg1  << "]: " << A_B1R0[bn1_ma_reg1] << " \n";
                 DATARECORD << "A_B1R1["<< bn1_ma_reg1  << "]: " << A_B1R1[bn1_ma_reg1] << " \n";
                 DATARECORD << "A_B1R2["<< bn1_ma_reg1  << "]: " << A_B1R2[bn1_ma_reg1] << " \n";
                 DATARECORD << "A_B1R3["<< bn1_ma_reg1  << "]: " << A_B1R3[bn1_ma_reg1] << " \n";				
                 DATARECORD << "A_B1R0["<< bn1_ma_reg2  << "]: " << A_B1R0[bn1_ma_reg2] << " \n";
                 DATARECORD << "A_B1R1["<< bn1_ma_reg2  << "]: " << A_B1R1[bn1_ma_reg2] << " \n";
                 DATARECORD << "A_B1R2["<< bn1_ma_reg2  << "]: " << A_B1R2[bn1_ma_reg2] << " \n";
                 DATARECORD << "A_B1R3["<< bn1_ma_reg2  << "]: " << A_B1R3[bn1_ma_reg2] << " \n";							
		     }
		    } 
 	    }
	}
	DATARECORD << "---------------------------------------------------\n";
	DATARECORD << "Radix-2 fft!!!!!!!!!!                              \n";
	DATARECORD << "---------------------------------------------------\n";
	DATARECORD << "---------------------------------------------------\n";
	
	// radix-2 FFT compute
	for(int i = 0; i < group ; i++){
		for(int j=0; j < radix; j++){
			gray_i  = Gray(i,group);
			BC_tmp  = j * group + gray_i;
			DATARECORD << "BC_tmp:"<< BC_tmp <<"\n";
			BC_RR = BC_tmp;
			DATARECORD << "BC_RR:"<< BC_RR <<"\n";
			BC_IFFT_Reorder_R4_R2(BC_RR,BC_Reorder);
			DATARECORD << "BC_Reorder:"<< BC_Reorder <<"\n";
			AGU_R4(BC_Reorder,bn_tmp,ma_tmp);
		    DATARECORD << "bn_tmp:"<< bn_tmp <<"\n";
			DATARECORD << "ma_tmp:"<< ma_tmp <<"\n";
			if(bn_tmp == 0){
				if(j < 2)bn0_bc_tmp = BC_tmp;
				DATARECORD << "---------------------------------------------------\n";
				DATARECORD << "Before multiplication and Radix-2 BU operation!!!\n";
				DATARECORD << "A_B0R0["<< ma_tmp  << "]: " << A_B0R0[ma_tmp] << " \n";
				DATARECORD << "A_B0R1["<< ma_tmp  << "]: " << A_B0R1[ma_tmp] << " \n";
				DATARECORD << "A_B0R2["<< ma_tmp  << "]: " << A_B0R2[ma_tmp] << " \n";
				DATARECORD << "A_B0R3["<< ma_tmp  << "]: " << A_B0R3[ma_tmp] << " \n";
				Radix2_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp]);
				Radix2_BU(A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
				DATARECORD << "---------------------------------------------------\n";
				DATARECORD << "After multiplication and Radix-2 BU operation!!!\n";
				DATARECORD << "A_B0R0["<< ma_tmp  << "]: " << A_B0R0[ma_tmp] << " \n";
				DATARECORD << "A_B0R1["<< ma_tmp  << "]: " << A_B0R1[ma_tmp] << " \n";
				DATARECORD << "A_B0R2["<< ma_tmp  << "]: " << A_B0R2[ma_tmp] << " \n";
				DATARECORD << "A_B0R3["<< ma_tmp  << "]: " << A_B0R3[ma_tmp] << " \n";
				if(j <  2)bn0_ma_reg1 = ma_tmp;
				if(j >= 2)bn0_ma_reg2 = ma_tmp;
			}else {
			    if(j < 2)bn1_bc_tmp = BC_tmp;
				DATARECORD << "---------------------------------------------------\n";
				DATARECORD << "Before multiplication and Radix-2 BU operation!!!\n";
				DATARECORD << "A_B1R0["<< ma_tmp  << "]: " << A_B1R0[ma_tmp] << " \n";
				DATARECORD << "A_B1R1["<< ma_tmp  << "]: " << A_B1R1[ma_tmp] << " \n";
				DATARECORD << "A_B1R2["<< ma_tmp  << "]: " << A_B1R2[ma_tmp] << " \n";
				DATARECORD << "A_B1R3["<< ma_tmp  << "]: " << A_B1R3[ma_tmp] << " \n";
			    Radix2_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp]);
			    Radix2_BU(A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
				DATARECORD << "---------------------------------------------------\n";
				DATARECORD << "After multiplication and Radix-2 BU operation!!!\n";
				DATARECORD << "A_B1R0["<< ma_tmp  << "]: " << A_B1R0[ma_tmp] << " \n";
				DATARECORD << "A_B1R1["<< ma_tmp  << "]: " << A_B1R1[ma_tmp] << " \n";
				DATARECORD << "A_B1R2["<< ma_tmp  << "]: " << A_B1R2[ma_tmp] << " \n";
				DATARECORD << "A_B1R3["<< ma_tmp  << "]: " << A_B1R3[ma_tmp] << " \n";	
				if(j < 2) bn1_ma_reg1 = ma_tmp;
				if(j >= 2)bn1_ma_reg2 = ma_tmp;
			}			
		}
	/*
    if(bn1_bc_tmp > bn0_bc_tmp){
	    //bn: 0 1 1 0
		//For example N = 32 point!
		//--------------------------------
		//BN0 , MA: 0 , r0 r16 r4 r20
		//BN1 , MA: 1 , r1 r17 r5 r21
		//BN1 , MA: 2 , r2 r18 r6 r22
		//BN0 , MA: 3 , r3 r19 r7 r23
	    //--------------------------------
		//-----After Relocation!----------
		//--------------------------------
		//BN0 , MA: 0 ,  r0   r1   r2   r3
		//BN1 , MA: 1 , r16  r17  r18  r19
		//BN1 , MA: 2 , r4   r5   r6   r7
		//BN0 , MA: 3 , r20  r21  r22  r23
		//--------------------------------
	    data_tmp   = A_B0R1[bn0_ma_reg1]; //r2
	    data_tmp_1 = A_B0R2[bn0_ma_reg1]; //r4
	    data_tmp_2 = A_B0R3[bn0_ma_reg1]; //r6
		A_B0R1[bn0_ma_reg1] = A_B1R0[bn1_ma_reg1];
		A_B0R2[bn0_ma_reg1] = A_B1R0[bn1_ma_reg2];
		A_B0R3[bn0_ma_reg1] = A_B0R0[bn0_ma_reg2];
		A_B1R0[bn1_ma_reg1] = data_tmp;
		A_B1R0[bn1_ma_reg2] = data_tmp_1;
		A_B0R0[bn0_ma_reg2] = data_tmp_2;
		//--------------------------------
		//BN0 , MA: 0 ,  r0   r1   r2   r3
		//BN1 , MA: 1 , r16   r17  r5  r21
		//BN1 , MA: 2 , r4    r18  r6  r22
		//BN0 , MA: 3 , r20   r19  r7  r23
		//--------------------------------
		data_tmp    = A_B1R1[bn1_ma_reg2];
		data_tmp_1  = A_B0R1[bn0_ma_reg2];
		A_B1R1[bn1_ma_reg2] = A_B1R2[bn1_ma_reg1];
		A_B0R1[bn0_ma_reg2] = A_B1R3[bn1_ma_reg1];
		A_B1R2[bn1_ma_reg1] = data_tmp;
		A_B1R3[bn1_ma_reg1] = data_tmp_1;
		//--------------------------------
		//BN0 , MA: 0 ,  r0   r1   r2   r3
        //BN1 , MA: 1 , r16   r17  r18 r19
        //BN1 , MA: 2 , r4     r5  r6  r22
        //BN0 , MA: 3 , r20   r21  r7  r23
        //--------------------------------
		data_tmp   =  A_B0R2[bn0_ma_reg2];
		A_B0R2[bn0_ma_reg2]  = A_B1R3[bn1_ma_reg2];
		A_B1R3[bn1_ma_reg2]  = data_tmp;
		//--------------------------------
		//BN0 , MA: 0 ,  r0   r1   r2   r3
		//BN1 , MA: 1 , r16   r17  r18 r19
		//BN1 , MA: 2 , r4     r5  r6   r7
		//BN0 , MA: 3 , r20   r21  r22 r23
        //--------------------------------
    }else {
	    //bn: 1 0 0 1
		//For example N = 32 point!
		//--------------------------------
		//BN1 , MA: 0 , r8  r24 r12 r28
		//BN0 , MA: 1 , r9  r25 r13 r29
		//BN0 , MA: 2 , r10 r26 r14 r30
		//BN1 , MA: 3 , r11 r27 r15 r31
	    //--------------------------------
		//-----After Relocation!----------
		//--------------------------------
		//BN1 , MA: 0 ,  r8  r9   r10  r11
		//BN0 , MA: 1 , r24  r25  r26  r27
		//BN0 , MA: 2 , r12  r13  r14  r15
		//BN1 , MA: 3 , r28  r29  r30  r31
		//--------------------------------
        data_tmp   = A_B1R1[bn1_ma_reg1];
        data_tmp_1 = A_B1R2[bn1_ma_reg1];
        data_tmp_2 = A_B1R3[bn1_ma_reg1];
		A_B1R1[bn1_ma_reg1] = A_B0R0[bn0_ma_reg1];
		A_B1R2[bn1_ma_reg1] = A_B0R0[bn0_ma_reg2];
        A_B1R3[bn1_ma_reg1] = A_B1R0[bn1_ma_reg2];
        A_B0R0[bn0_ma_reg1] = data_tmp;
        A_B0R0[bn0_ma_reg2] = data_tmp_1;
        A_B1R0[bn1_ma_reg2] = data_tmp_2;
		//--------------------------------
        //BN1 , MA: 0 , r8  r9  r10 r11
        //BN0 , MA: 1 , r24 r25 r13 r29
        //BN0 , MA: 2 , r12 r26 r14 r30
        //BN1 , MA: 3 , r28 r27 r15 r31
        //--------------------------------
        data_tmp   = A_B0R1[bn0_ma_reg2];
        data_tmp_1 = A_B1R1[bn1_ma_reg2];
        A_B0R1[bn0_ma_reg2] = A_B0R2[bn0_ma_reg1];
        A_B1R1[bn1_ma_reg2] = A_B0R3[bn0_ma_reg1];
        A_B0R2[bn0_ma_reg1] = data_tmp;
        A_B0R3[bn0_ma_reg1] = data_tmp_1;
        //--------------------------------
        //BN1 , MA: 0 , r8  r9  r10 r11
        //BN0 , MA: 1 , r24 r25 r26 r27
        //BN0 , MA: 2 , r12 r13 r14 r30
        //BN1 , MA: 3 , r28 r29 r15 r31
        //--------------------------------
        data_tmp = A_B0R3[bn0_ma_reg2];
        A_B0R3[bn0_ma_reg2] = A_B1R2[bn1_ma_reg2];
        A_B1R2[bn1_ma_reg2] = data_tmp;
        //--------------------------------
        //BN1 , MA: 0 , r8  r9  r10 r11
        //BN0 , MA: 1 , r24 r25 r26 r27
        //BN0 , MA: 2 , r12 r13 r14 r15
        //BN1 , MA: 3 , r28 r29 r30 r31
        //--------------------------------  
     }
	 */
	}
	
    DATARECORD << "---------------------------------------------------------\n";
	DATARECORD << "*********************************************************\n";
	DATARECORD << "                  INTT OUTPUT !!!!                       \n";
	DATARECORD << "*********************************************************\n";
	DATARECORD << "---------------------------------------------------------\n";
	
    int BC_bit_size;
	BC_bit_size = (int)ceil(log2(N/4));
	//DATARECORD << "BC_bit_size: "<< BC_bit_size <<"\n";
	int index_out;
	for(int i = 0; i < group; i++){
    	for(int j = 0;j < radix;j++){
			DATARECORD << "------------------------------------------------------\n";
    		BC_tmp  = j * group + i;
			DATARECORD << "BC_tmp: " << BC_tmp <<"\n";
			INTT_REORDERBC_R4_R2_OUT(BC_tmp,BC_RR);
			DATARECORD << "BC_RR: " << BC_RR <<"\n";
    		AGU_R4(BC_RR,bn_tmp,ma_tmp);
			DATARECORD << "bn_tmp: " << bn_tmp <<"\n";
			DATARECORD << "ma_tmp: " << ma_tmp <<"\n";
			index_out = (int)i * ( N / group) + 4 * j ;
    		if(bn_tmp == 0){
			   MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],IN,p);
			   MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],IN,p);
			   MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],IN,p);
			   MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],IN,p);
     		   A[index_out+0] = A_B0R0[ma_tmp];
    		   A[index_out+1] = A_B0R1[ma_tmp];
    		   A[index_out+2] = A_B0R2[ma_tmp];
    		   A[index_out+3] = A_B0R3[ma_tmp];
			   DATARECORD << "BN: 0 \n";
			   DATARECORD << "A["<< index_out      << "]: " << A[index_out+0] << " \n";
			   DATARECORD << "A["<< index_out + 1  << "]: " << A[index_out+1] << " \n";
			   DATARECORD << "A["<< index_out + 2  << "]: " << A[index_out+2] << " \n";
			   DATARECORD << "A["<< index_out + 3  << "]: " << A[index_out+3] << " \n";
    		}
    	    else {
			   MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],IN,p);
			   MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],IN,p);
			   MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],IN,p);
			   MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],IN,p);
    		   A[index_out+0] = A_B1R0[ma_tmp];
    		   A[index_out+1] = A_B1R1[ma_tmp];
    		   A[index_out+2] = A_B1R2[ma_tmp];
    		   A[index_out+3] = A_B1R3[ma_tmp];
			   DATARECORD << "BN: 1 \n";
			   DATARECORD << "A["<< index_out      << "]: " << A[index_out+0] << " \n";
			   DATARECORD << "A["<< index_out + 1  << "]: " << A[index_out+1] << " \n";
			   DATARECORD << "A["<< index_out + 2  << "]: " << A[index_out+2] << " \n";
			   DATARECORD << "A["<< index_out + 3  << "]: " << A[index_out+3] << " \n";
    		}
    	}	   
	}
	
	DATARECORD << "---------------------------------------------------------\n";
	DATARECORD << "----                  OVER!!!!!!!                   -----\n";
	DATARECORD << "---------------------------------------------------------\n";
	DATARECORD << "---------------------------------------------------------\n";
	DATARECORD.close();
}
//----------------------------------------------
//mixed radix fft
//radix 16 and raidx-2
void NTTSPMB::NTT_r16_r2(std::vector<ZZ> &A,std::vector<ZZ> &B0R0,std::vector<ZZ> &B0R1,std::vector<ZZ> &B0R2,std::vector<ZZ> &B0R3,
std::vector<ZZ> &B0R4,std::vector<ZZ> &B0R5,std::vector<ZZ> &B0R6,std::vector<ZZ> &B0R7,std::vector<ZZ> &B0R8,std::vector<ZZ> &B0R9,
std::vector<ZZ> &B0R10,std::vector<ZZ> &B0R11,std::vector<ZZ> &B0R12,std::vector<ZZ> &B0R13,std::vector<ZZ> &B0R14,std::vector<ZZ> &B0R15,
std::vector<ZZ> &B1R0,std::vector<ZZ> &B1R1,std::vector<ZZ> &B1R2,std::vector<ZZ> &B1R3,std::vector<ZZ> &B1R4,std::vector<ZZ> &B1R5,
std::vector<ZZ> &B1R6,std::vector<ZZ> &B1R7,std::vector<ZZ> &B1R8,std::vector<ZZ> &B1R9,std::vector<ZZ> &B1R10,std::vector<ZZ> &B1R11,
std::vector<ZZ> &B1R12,std::vector<ZZ> &B1R13,std::vector<ZZ> &B1R14,std::vector<ZZ> &B1R15){
	
	int            Stage;//stage	
	int            word_size;
	int            offset;
	int            BC;
	int            BC_tmp;
    int            group;
    int            ma_tmp;
	int            bn_tmp;
    int            bit_tmp;
	int            bn0_bc_tmp;//frist in bc data
	int            bn1_bc_tmp;//frist in bc data
	int            bn0_ma_reg1;
	int            bn0_ma_reg2;
	int            bn0_ma_reg3;
	int            bn0_ma_reg4;
	int            bn0_ma_reg5;
	int            bn0_ma_reg6;
	int            bn0_ma_reg7;
	int            bn0_ma_reg8;
	int            bn1_ma_reg1;
	int            bn1_ma_reg2;
	int            bn1_ma_reg3;
	int            bn1_ma_reg4;
	int            bn1_ma_reg5;
	int            bn1_ma_reg6;
	int            bn1_ma_reg7;
	int            bn1_ma_reg8;
	int            gray_i;
    int            BC_WIDTH;
	int            tw_modulus;
    int            tw_modulus_tmp;
    double         Stage_double;	
    std::vector<int> bit_array_tmp;
	
	std::ofstream DATARECORD("./NTT_R16_R2_SPMB.txt");
	std::ofstream r16_r2_SPMB_TF_dg("./my_print_data/r16_r2_SPMB_TF_dg.txt");
	std::ofstream spmb_r16_r2("./SPMB_tw/spmb_r16_r2.txt");
	std::ofstream DTFAG_golden_st0("./SPMB_tw/DTFAG_golden_st0_8192.txt");
	std::ofstream DTFAG_golden_st1("./SPMB_tw/DTFAG_golden_st1_8192.txt");
	std::ofstream DTFAG_golden_st2("./SPMB_tw/DTFAG_golden_st2_8192.txt");

	//-------------Bit Reverse--------------
	BitOperate RR, IntToVec, VecToInt;
	long long bit_width = log2(N);
	long long bit_width_s = log(radix)/log(2);
	vector<long long> bit_array;
	//----------------------------------------

	//radix-16 Stage
    Stage_double  = log2(N);
    Stage       =  (int)floor(Stage_double/4);
	BC_WIDTH    =  (int)ceil(log2(N/16));
	offset      =  (int)N /  16;
	word_size   =  (int)N / (2 * 16);
	group       =  (int)N / (256);
	tw_modulus  =  (int)N /  16;
	bit_array_tmp.resize(BC_WIDTH);

    DATARECORD << "group: "    << group << "\n";
    DATARECORD << "BC_WIDTH: " << BC_WIDTH << "\n";
	
	
	ZZ               data_tmp_1;
	ZZ               data_tmp_2;
	ZZ               data_tmp_3;
	ZZ               data_tmp_4;
	ZZ               data_tmp_5;
	ZZ               data_tmp_6;
	ZZ               data_tmp_7;
	ZZ               data_tmp_8;
	ZZ               data_tmp_9;
	ZZ               data_tmp_10;
	ZZ               data_tmp_11;
	ZZ               data_tmp_12;
	ZZ               data_tmp_13;
	ZZ               data_tmp_14;
	ZZ               data_tmp_15;
	std::vector<ZZ>  A_B0R0;
	std::vector<ZZ>  A_B0R1;
	std::vector<ZZ>  A_B0R2;
	std::vector<ZZ>  A_B0R3;
	std::vector<ZZ>  A_B0R4;
	std::vector<ZZ>  A_B0R5;
	std::vector<ZZ>  A_B0R6;
	std::vector<ZZ>  A_B0R7;
	std::vector<ZZ>  A_B0R8;
	std::vector<ZZ>  A_B0R9;
	std::vector<ZZ>  A_B0R10;
	std::vector<ZZ>  A_B0R11;
	std::vector<ZZ>  A_B0R12;
	std::vector<ZZ>  A_B0R13;
	std::vector<ZZ>  A_B0R14;
	std::vector<ZZ>  A_B0R15;
	std::vector<ZZ>  A_B1R0;
	std::vector<ZZ>  A_B1R1;
	std::vector<ZZ>  A_B1R2;
	std::vector<ZZ>  A_B1R3;
	std::vector<ZZ>  A_B1R4;
	std::vector<ZZ>  A_B1R5;
	std::vector<ZZ>  A_B1R6;
	std::vector<ZZ>  A_B1R7;
	std::vector<ZZ>  A_B1R8;
	std::vector<ZZ>  A_B1R9;
	std::vector<ZZ>  A_B1R10;
	std::vector<ZZ>  A_B1R11;
	std::vector<ZZ>  A_B1R12;
	std::vector<ZZ>  A_B1R13;
	std::vector<ZZ>  A_B1R14;
	std::vector<ZZ>  A_B1R15;
	
	A_B0R0.resize(word_size);
	A_B0R1.resize(word_size);
	A_B0R2.resize(word_size);
	A_B0R3.resize(word_size);
	A_B0R4.resize(word_size);
	A_B0R5.resize(word_size);
	A_B0R6.resize(word_size);
	A_B0R7.resize(word_size);
	A_B0R8.resize(word_size);
	A_B0R9.resize(word_size);
	A_B0R10.resize(word_size);
	A_B0R11.resize(word_size);
	A_B0R12.resize(word_size);
	A_B0R13.resize(word_size);
	A_B0R14.resize(word_size);
	A_B0R15.resize(word_size);
	A_B1R0.resize(word_size);
	A_B1R1.resize(word_size);
	A_B1R2.resize(word_size);
	A_B1R3.resize(word_size);
	A_B1R4.resize(word_size);
	A_B1R5.resize(word_size);
	A_B1R6.resize(word_size);
	A_B1R7.resize(word_size);
	A_B1R8.resize(word_size);
	A_B1R9.resize(word_size);
	A_B1R10.resize(word_size);
	A_B1R11.resize(word_size);
	A_B1R12.resize(word_size);
	A_B1R13.resize(word_size);
	A_B1R14.resize(word_size);
	A_B1R15.resize(word_size);
	//----------------------------------------------------
	B0R0.resize(word_size);
	B0R1.resize(word_size);
	B0R2.resize(word_size);
	B0R3.resize(word_size);
	B0R4.resize(word_size);
	B0R5.resize(word_size);
	B0R6.resize(word_size);
	B0R7.resize(word_size);
	B0R8.resize(word_size);
	B0R9.resize(word_size);
	B0R10.resize(word_size);
	B0R11.resize(word_size);
	B0R12.resize(word_size);
	B0R13.resize(word_size);
	B0R14.resize(word_size);
	B0R15.resize(word_size);
	B1R0.resize(word_size);
	B1R1.resize(word_size);
	B1R2.resize(word_size);
	B1R3.resize(word_size);
	B1R4.resize(word_size);
	B1R5.resize(word_size);
	B1R6.resize(word_size);
	B1R7.resize(word_size);
	B1R8.resize(word_size);
	B1R9.resize(word_size);
	B1R10.resize(word_size);
	B1R11.resize(word_size);
	B1R12.resize(word_size);
	B1R13.resize(word_size);
	B1R14.resize(word_size);
	B1R15.resize(word_size);
	//----------------------------------------------------
	int length;
	ZZ  factor;   //base factor
	ZZ  factor_t; //acctually mul factor
	ZZ  factor_2t;
	ZZ  factor_3t;
	ZZ  factor_4t;
	ZZ  factor_5t;
	ZZ  factor_6t;
	ZZ  factor_7t;
	ZZ  factor_8t;
	ZZ  factor_9t;
	ZZ  factor_10t;
	ZZ  factor_11t;
	ZZ  factor_12t;
	ZZ  factor_13t;
	ZZ  factor_14t;
	ZZ  factor_15t;
    
	//init load data
    for(int i = 0; i < group; i++){
		for(int j = 0 ; j < radix ; j++){
			bn_tmp = 0;
			ma_tmp = 0;
			gray_i  = Gray(i,group);
			BC = j * group + gray_i;
            //bit calculate
            for(int j=0; j < BC_WIDTH;j++){
                bit_tmp = BC % 2;
                BC = BC >> 1;
                bit_array_tmp[j] = bit_tmp;
            } 
			gray_i  = Gray(i,group);
			BC = j * group + gray_i;
            for(int rs = 0; rs < BC_WIDTH; rs++){
                if((bit_array_tmp[rs] == 1) && (rs != 0)) ma_tmp = ma_tmp + exp2((rs-1)); 
                bn_tmp = bn_tmp ^ bit_array_tmp[rs];
            }
			if(bn_tmp == 0){
				A_B0R0[ma_tmp]  = A[BC];
				A_B0R1[ma_tmp]  = A[BC +      offset];
				A_B0R2[ma_tmp]  = A[BC + 2  * offset];
				A_B0R3[ma_tmp]  = A[BC + 3  * offset];
				A_B0R4[ma_tmp]  = A[BC + 4  * offset];
				A_B0R5[ma_tmp]  = A[BC + 5  * offset];
				A_B0R6[ma_tmp]  = A[BC + 6  * offset];
				A_B0R7[ma_tmp]  = A[BC + 7  * offset];
				A_B0R8[ma_tmp]  = A[BC + 8  * offset];
				A_B0R9[ma_tmp]  = A[BC + 9  * offset];
				A_B0R10[ma_tmp] = A[BC + 10 * offset];
				A_B0R11[ma_tmp] = A[BC + 11 * offset];
				A_B0R12[ma_tmp] = A[BC + 12 * offset];
				A_B0R13[ma_tmp] = A[BC + 13 * offset];
				A_B0R14[ma_tmp] = A[BC + 14 * offset];
				A_B0R15[ma_tmp] = A[BC + 15 * offset];
			}else {
				A_B1R0[ma_tmp]  = A[BC];
				A_B1R1[ma_tmp]  = A[BC +     offset];
				A_B1R2[ma_tmp]  = A[BC + 2 * offset];
				A_B1R3[ma_tmp]  = A[BC + 3 * offset];
				A_B1R4[ma_tmp]  = A[BC + 4 * offset];
				A_B1R5[ma_tmp]  = A[BC + 5 * offset];
				A_B1R6[ma_tmp]  = A[BC + 6 * offset];
				A_B1R7[ma_tmp]  = A[BC + 7 * offset];
				A_B1R8[ma_tmp]  = A[BC + 8 * offset];
				A_B1R9[ma_tmp]  = A[BC + 9 * offset];
				A_B1R10[ma_tmp] = A[BC + 10 * offset];
				A_B1R11[ma_tmp] = A[BC + 11 * offset];
				A_B1R12[ma_tmp] = A[BC + 12 * offset];
				A_B1R13[ma_tmp] = A[BC + 13 * offset];
				A_B1R14[ma_tmp] = A[BC + 14 * offset];
				A_B1R15[ma_tmp] = A[BC + 15 * offset];
			}
		}
	}
	ma_tmp = 0;
	bn_tmp = 0;
	BC     = 0;
	int tw_degree = 1; // siang
	std::cout << "init load over! \n";
	DATARECORD <<"radix-16 computing stage:  "<< Stage <<"\n";
	//need modify to mult by twiddle factor
	for(int s = 0; s < Stage; s++){
		if(s == 0)factor = W;
		else {
			SqrMod(factor,factor,p);
			SqrMod(factor,factor,p);
			SqrMod(factor,factor,p);
			SqrMod(factor,factor,p);
			tw_degree = tw_degree * 16;
		}
		DATARECORD <<"---------------------------------\n";
		DATARECORD <<"Now Stage: "<< s <<"\n";
		r16_r2_SPMB_TF_dg << "Now Stage: "<< s <<"\n";
		spmb_r16_r2 <<"Now Stage: "<< s <<"\n";
		spmb_r16_r2 <<"twiddle factor : "<< factor <<"\n";
		tw_modulus_tmp  = tw_modulus >> ( 4 * s);
		for(int i = 0 ;i < group;i++){
			bn0_bc_tmp  = 0;
			bn1_bc_tmp  = 0;
			spmb_r16_r2 << "-----------------i = " << i << "-------------------" << std::endl;
			r16_r2_SPMB_TF_dg << "-------------------------------\n";
			for(int j = 0;j < radix;j++){
				DATARECORD << "-------------------------------------"<< std::endl;
				DATARECORD << "i = " << i << ", j = " << j << std::endl;
				gray_i  = Gray(i,group);
				BC_tmp  = j * group + gray_i;
				DATARECORD << "BC_tmp: " << BC_tmp << "\n";
				if(s == Stage - 1) RR_R16_R2(BC_tmp,(4 * s - 3),BC);
				else RR_R16(BC_tmp,s,BC);
				DATARECORD << "After RR_R16 , BC : " << BC << "\n";
				length = BC % tw_modulus_tmp;
				r16_r2_SPMB_TF_dg << "length: " << length << "\n";
				PowerMod(factor_t,factor,length,p);
				AGU_R16(BC,bn_tmp,ma_tmp);

				//-----------compute data idx-------------
				spmb_r16_r2 /*<< "BC = " << BC */<< "BC_tmp = " << BC_tmp << std::endl;
				
				IntToVec.IntToVec(BC_tmp, N, bit_array);
				rotate(bit_array.begin(), bit_array.begin()+ bit_width_s*s , bit_array.end());
				long long Data = VecToInt.VecToInt(bit_array, N);
				spmb_r16_r2 << "Data_index = ";
                spmb_r16_r2 << "( " ;					
				for(int k = 0; k < radix ; k++ ){
					spmb_r16_r2 << Data + k*(1<<(bit_width-bit_width_s-bit_width_s*s)) <<" ";	
				}
				spmb_r16_r2 << ") " << std::endl;
				spmb_r16_r2 << ", (w^" << 0 << ", w^" << tw_degree * length << ", w^" << tw_degree * length * 2 
				<< ", w^" << tw_degree * length * 3  << ", w^" << tw_degree * length * 4 << ", w^" << tw_degree * length * 5 
				<< ", w^" << tw_degree * length * 6  << ", w^" << tw_degree * length * 7 << ", w^" << tw_degree * length * 8 
				<< ", w^" << tw_degree * length * 9  << ", w^" << tw_degree * length * 10 << ", w^" << tw_degree * length * 11 
				<< ", w^" << tw_degree * length * 12  << ", w^" << tw_degree * length * 13 << ", w^" << tw_degree * length * 14 
				<< ", w^" << tw_degree * length * 15
				<< ")" <<std::endl;
				//-----------------------------------------

				DATARECORD << "BN : " << bn_tmp << "\n";
				DATARECORD << "MA : " << ma_tmp << "\n";
				if(bn_tmp == 0){
					if(j < 2)bn0_bc_tmp = BC_tmp;
					PowerMod(factor_2t,factor_t,2,p);
					PowerMod(factor_3t,factor_t,3,p);
					PowerMod(factor_4t,factor_t,4,p);
					PowerMod(factor_5t,factor_t,5,p);
					PowerMod(factor_6t,factor_t,6,p);
					PowerMod(factor_7t,factor_t,7,p);
					PowerMod(factor_8t,factor_t,8,p);
					PowerMod(factor_9t,factor_t,9,p);
					PowerMod(factor_10t,factor_t,10,p);
					PowerMod(factor_11t,factor_t,11,p);
					PowerMod(factor_12t,factor_t,12,p);
					PowerMod(factor_13t,factor_t,13,p);
					PowerMod(factor_14t,factor_t,14,p);
					PowerMod(factor_15t,factor_t,15,p);
					DATARECORD <<"facotr_t : " << factor_t  << ", length = " << length * 1 << " \n";
					DATARECORD <<"facotr_2t : " << factor_2t << ", length = " << length * 2 << " \n";
					DATARECORD <<"facotr_3t : " << factor_3t << ", length = " << length * 3 << " \n";
					DATARECORD <<"facotr_4t : " << factor_4t << ", length = " << length * 4 << " \n";
					DATARECORD <<"facotr_5t : " << factor_5t << ", length = " << length * 5 << " \n";
					DATARECORD <<"facotr_6t : " << factor_6t << ", length = " << length * 6 << " \n";
					DATARECORD <<"facotr_7t : " << factor_7t << ", length = " << length * 7 << " \n";
					DATARECORD <<"facotr_8t : " << factor_8t << ", length = " << length * 8 << " \n";
					DATARECORD <<"facotr_9t : " << factor_9t << ", length = " << length * 9 << " \n";
					DATARECORD <<"facotr_10t : " << factor_10t << ", length = " << length * 10 << " \n";
					DATARECORD <<"facotr_11t : " << factor_11t << ", length = " << length * 11 << " \n";
					DATARECORD <<"facotr_12t : " << factor_12t << ", length = " << length * 12 << " \n";
					DATARECORD <<"facotr_13t : " << factor_13t << ", length = " << length * 13 << " \n";
					DATARECORD <<"facotr_14t : " << factor_14t << ", length = " << length * 14 << " \n";
					DATARECORD <<"facotr_15t : " << factor_15t << ", length = " << length * 15 << " \n";
					DATARECORD << "Before Radix-16 butterfly unit operation!!! \n";
				    //DATARECORD << "A_B0R0["<< ma_tmp <<"]: " << A_B0R0[ma_tmp]<<"\n";
				    //DATARECORD << "A_B0R1["<< ma_tmp <<"]: " << A_B0R1[ma_tmp]<<"\n";
				    //DATARECORD << "A_B0R2["<< ma_tmp <<"]: " << A_B0R2[ma_tmp]<<"\n";
				    //DATARECORD << "A_B0R3["<< ma_tmp <<"]: " << A_B0R3[ma_tmp]<<"\n";
				    //DATARECORD << "A_B0R4["<< ma_tmp <<"]: " << A_B0R4[ma_tmp]<<"\n";
				    //DATARECORD << "A_B0R5["<< ma_tmp <<"]: " << A_B0R5[ma_tmp]<<"\n";
				    //DATARECORD << "A_B0R6["<< ma_tmp <<"]: " << A_B0R6[ma_tmp]<<"\n";
				    //DATARECORD << "A_B0R7["<< ma_tmp <<"]: " << A_B0R7[ma_tmp]<<"\n";
				    //DATARECORD << "A_B0R8["<< ma_tmp <<"]: " << A_B0R8[ma_tmp]<<"\n";
				    //DATARECORD << "A_B0R9["<< ma_tmp <<"]: " << A_B0R9[ma_tmp]<<"\n";
				    //DATARECORD << "A_B0R10["<< ma_tmp <<"]: "<< A_B0R10[ma_tmp]<<"\n";
				    //DATARECORD << "A_B0R11["<< ma_tmp <<"]: "<< A_B0R11[ma_tmp]<<"\n";
				    //DATARECORD << "A_B0R12["<< ma_tmp <<"]: "<< A_B0R12[ma_tmp]<<"\n";
				    //DATARECORD << "A_B0R13["<< ma_tmp <<"]: "<< A_B0R13[ma_tmp]<<"\n";
				    //DATARECORD << "A_B0R14["<< ma_tmp <<"]: "<< A_B0R14[ma_tmp]<<"\n";
				    //DATARECORD << "A_B0R15["<< ma_tmp <<"]: "<< A_B0R15[ma_tmp]<<"\n";
					Radix16_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp], A_B0R2[ma_tmp], A_B0R3[ma_tmp],A_B0R4[ma_tmp],
							   A_B0R5[ma_tmp],A_B0R6[ma_tmp], A_B0R7[ma_tmp], A_B0R8[ma_tmp],A_B0R9[ma_tmp],
							   A_B0R10[ma_tmp],A_B0R11[ma_tmp],A_B0R12[ma_tmp],A_B0R13[ma_tmp],A_B0R14[ma_tmp],
							   A_B0R15[ma_tmp]);
					DATARECORD << "After Radix-16 butterfly unit operation!!! \n";
				    //DATARECORD << "A_B0R0["<< ma_tmp <<"]: " << A_B0R0[ma_tmp]<<"\n";
				    //DATARECORD << "A_B0R1["<< ma_tmp <<"]: " << A_B0R1[ma_tmp]<<"\n";
				    //DATARECORD << "A_B0R2["<< ma_tmp <<"]: " << A_B0R2[ma_tmp]<<"\n";
				    //DATARECORD << "A_B0R3["<< ma_tmp <<"]: " << A_B0R3[ma_tmp]<<"\n";
				    //DATARECORD << "A_B0R4["<< ma_tmp <<"]: " << A_B0R4[ma_tmp]<<"\n";
				    //DATARECORD << "A_B0R5["<< ma_tmp <<"]: " << A_B0R5[ma_tmp]<<"\n";
				    //DATARECORD << "A_B0R6["<< ma_tmp <<"]: " << A_B0R6[ma_tmp]<<"\n";
				    //DATARECORD << "A_B0R7["<< ma_tmp <<"]: " << A_B0R7[ma_tmp]<<"\n";
				    //DATARECORD << "A_B0R8["<< ma_tmp <<"]: " << A_B0R8[ma_tmp]<<"\n";
				    //DATARECORD << "A_B0R9["<< ma_tmp <<"]: " << A_B0R9[ma_tmp]<<"\n";
				    //DATARECORD << "A_B0R10["<< ma_tmp <<"]: "<< A_B0R10[ma_tmp]<<"\n";
				    //DATARECORD << "A_B0R11["<< ma_tmp <<"]: "<< A_B0R11[ma_tmp]<<"\n";
				    //DATARECORD << "A_B0R12["<< ma_tmp <<"]: "<< A_B0R12[ma_tmp]<<"\n";
				    //DATARECORD << "A_B0R13["<< ma_tmp <<"]: "<< A_B0R13[ma_tmp]<<"\n";
				    //DATARECORD << "A_B0R14["<< ma_tmp <<"]: "<< A_B0R14[ma_tmp]<<"\n";
				    //DATARECORD << "A_B0R15["<< ma_tmp <<"]: "<< A_B0R15[ma_tmp]<<"\n";
                    MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],factor_t,p);
					MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],factor_2t,p);
					MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],factor_3t,p);
					MulMod(A_B0R4[ma_tmp],A_B0R4[ma_tmp],factor_4t,p);
					MulMod(A_B0R5[ma_tmp],A_B0R5[ma_tmp],factor_5t,p);
					MulMod(A_B0R6[ma_tmp],A_B0R6[ma_tmp],factor_6t,p);
					MulMod(A_B0R7[ma_tmp],A_B0R7[ma_tmp],factor_7t,p);
					MulMod(A_B0R8[ma_tmp],A_B0R8[ma_tmp],factor_8t,p);
					MulMod(A_B0R9[ma_tmp],A_B0R9[ma_tmp],factor_9t,p);
					MulMod(A_B0R10[ma_tmp],A_B0R10[ma_tmp],factor_10t,p);
					MulMod(A_B0R11[ma_tmp],A_B0R11[ma_tmp],factor_11t,p);
					MulMod(A_B0R12[ma_tmp],A_B0R12[ma_tmp],factor_12t,p);
					MulMod(A_B0R13[ma_tmp],A_B0R13[ma_tmp],factor_13t,p);
					MulMod(A_B0R14[ma_tmp],A_B0R14[ma_tmp],factor_14t,p);
					MulMod(A_B0R15[ma_tmp],A_B0R15[ma_tmp],factor_15t,p);
					if(j <  2)bn0_ma_reg1 = ma_tmp;
					if((j >= 2)  && (j < 4))bn0_ma_reg2 = ma_tmp;
					if((j >= 4)  && (j < 6))bn0_ma_reg3 = ma_tmp;
					if((j >= 6)  && (j < 8))bn0_ma_reg4 = ma_tmp;
					if((j >= 8)  && (j < 10))bn0_ma_reg5 = ma_tmp;
					if((j >= 10) && (j < 12))bn0_ma_reg6 = ma_tmp;
					if((j >= 12) && (j < 14))bn0_ma_reg7 = ma_tmp;
					if((j >= 14) && (j < 16))bn0_ma_reg8 = ma_tmp;
					//----------------DTFAG golden pattern------------------
					switch(s){
						case 0:
							DTFAG_golden_st0 << 1 		  << " \n ";
							DTFAG_golden_st0 << factor_t  << " \n ";
							DTFAG_golden_st0 << factor_2t << " \n ";
							DTFAG_golden_st0 << factor_3t << " \n ";
							DTFAG_golden_st0 << factor_4t << " \n ";
							DTFAG_golden_st0 << factor_5t  << " \n ";
							DTFAG_golden_st0 << factor_6t << " \n ";
							DTFAG_golden_st0 << factor_7t << " \n ";
							DTFAG_golden_st0 << factor_8t << " \n ";
							DTFAG_golden_st0 << factor_9t  << " \n ";
							DTFAG_golden_st0 << factor_10t << " \n ";
							DTFAG_golden_st0 << factor_11t << " \n ";
							DTFAG_golden_st0 << factor_12t << " \n ";
							DTFAG_golden_st0 << factor_13t  << " \n ";
							DTFAG_golden_st0 << factor_14t << " \n ";
							DTFAG_golden_st0 << factor_15t << " \n ";
							break;
						case 1:
							DTFAG_golden_st1 << 1 		  << " \n ";
							DTFAG_golden_st1 << factor_t  << " \n ";
							DTFAG_golden_st1 << factor_2t << " \n ";
							DTFAG_golden_st1 << factor_3t << " \n ";
							DTFAG_golden_st1 << factor_4t << " \n ";
							DTFAG_golden_st1 << factor_5t  << " \n ";
							DTFAG_golden_st1 << factor_6t << " \n ";
							DTFAG_golden_st1 << factor_7t << " \n ";
							DTFAG_golden_st1 << factor_8t << " \n ";
							DTFAG_golden_st1 << factor_9t  << " \n ";
							DTFAG_golden_st1 << factor_10t << " \n ";
							DTFAG_golden_st1 << factor_11t << " \n ";
							DTFAG_golden_st1 << factor_12t << " \n ";
							DTFAG_golden_st1 << factor_13t  << " \n ";
							DTFAG_golden_st1 << factor_14t << " \n ";
							DTFAG_golden_st1 << factor_15t << " \n ";
							break;
						case 2:
							DTFAG_golden_st2 << 1 		  << " \n ";
							DTFAG_golden_st2 << factor_t  << " \n ";
							DTFAG_golden_st2 << factor_2t << " \n ";
							DTFAG_golden_st2 << factor_3t << " \n ";
							DTFAG_golden_st2 << factor_4t << " \n ";
							DTFAG_golden_st2 << factor_5t  << " \n ";
							DTFAG_golden_st2 << factor_6t << " \n ";
							DTFAG_golden_st2 << factor_7t << " \n ";
							DTFAG_golden_st2 << factor_8t << " \n ";
							DTFAG_golden_st2 << factor_9t  << " \n ";
							DTFAG_golden_st2 << factor_10t << " \n ";
							DTFAG_golden_st2 << factor_11t << " \n ";
							DTFAG_golden_st2 << factor_12t << " \n ";
							DTFAG_golden_st2 << factor_13t  << " \n ";
							DTFAG_golden_st2 << factor_14t << " \n ";
							DTFAG_golden_st2 << factor_15t << " \n ";
							break;
						default:
							break;
					}
					//--------------------------------------------------
				}
			    else {
					if(j < 2)bn1_bc_tmp = BC_tmp;
					PowerMod(factor_2t,factor_t,2,p);
					PowerMod(factor_3t,factor_t,3,p);
					PowerMod(factor_4t,factor_t,4,p);
					PowerMod(factor_5t,factor_t,5,p);
					PowerMod(factor_6t,factor_t,6,p);
					PowerMod(factor_7t,factor_t,7,p);
					PowerMod(factor_8t,factor_t,8,p);
					PowerMod(factor_9t,factor_t,9,p);
					PowerMod(factor_10t,factor_t,10,p);
					PowerMod(factor_11t,factor_t,11,p);
					PowerMod(factor_12t,factor_t,12,p);
					PowerMod(factor_13t,factor_t,13,p);
					PowerMod(factor_14t,factor_t,14,p);
					PowerMod(factor_15t,factor_t,15,p);
					DATARECORD <<"facotr_t : " << factor_t  << ", length = " << length * 1 << " \n";
					DATARECORD <<"facotr_2t : " << factor_2t << ", length = " << length * 2 << " \n";
					DATARECORD <<"facotr_3t : " << factor_3t << ", length = " << length * 3 << " \n";
					DATARECORD <<"facotr_4t : " << factor_4t  << ", length = " << length * 4 << " \n";
					DATARECORD <<"facotr_5t : " << factor_5t << ", length = " << length * 5 << " \n";
					DATARECORD <<"facotr_6t : " << factor_6t << ", length = " << length * 6 << " \n";
					DATARECORD <<"facotr_7t : " << factor_7t  << ", length = " << length * 7 << " \n";
					DATARECORD <<"facotr_8t : " << factor_8t << ", length = " << length * 8 << " \n";
					DATARECORD <<"facotr_9t : " << factor_9t << ", length = " << length * 9 << " \n";
					DATARECORD <<"facotr_10t : " << factor_10t  << ", length = " << length * 10 << " \n";
					DATARECORD <<"facotr_11t : " << factor_11t << ", length = " << length * 11 << " \n";
					DATARECORD <<"facotr_12t : " << factor_12t << ", length = " << length * 12 << " \n";
					DATARECORD <<"facotr_13t : " << factor_13t  << ", length = " << length * 13 << " \n";
					DATARECORD <<"facotr_14t : " << factor_14t << ", length = " << length * 14 << " \n";
					DATARECORD <<"facotr_15t : " << factor_15t << ", length = " << length * 15 << " \n";
					DATARECORD << "Before Radix-16 butterfly unit operation!!! \n";
				    //DATARECORD << "A_B1R0["<< ma_tmp <<"]: " << A_B1R0[ma_tmp]<<"\n";
				    //DATARECORD << "A_B1R1["<< ma_tmp <<"]: " << A_B1R1[ma_tmp]<<"\n";
				    //DATARECORD << "A_B1R2["<< ma_tmp <<"]: " << A_B1R2[ma_tmp]<<"\n";
				    //DATARECORD << "A_B1R3["<< ma_tmp <<"]: " << A_B1R3[ma_tmp]<<"\n";
				    //DATARECORD << "A_B1R4["<< ma_tmp <<"]: " << A_B1R4[ma_tmp]<<"\n";
				    //DATARECORD << "A_B1R5["<< ma_tmp <<"]: " << A_B1R5[ma_tmp]<<"\n";
				    //DATARECORD << "A_B1R6["<< ma_tmp <<"]: " << A_B1R6[ma_tmp]<<"\n";
				    //DATARECORD << "A_B1R7["<< ma_tmp <<"]: " << A_B1R7[ma_tmp]<<"\n";
				    //DATARECORD << "A_B1R8["<< ma_tmp <<"]: " << A_B1R8[ma_tmp]<<"\n";
				    //DATARECORD << "A_B1R9["<< ma_tmp <<"]: " << A_B1R9[ma_tmp]<<"\n";
				    //DATARECORD << "A_B1R10["<< ma_tmp <<"]: "<< A_B1R10[ma_tmp]<<"\n";
				    //DATARECORD << "A_B1R11["<< ma_tmp <<"]: "<< A_B1R11[ma_tmp]<<"\n";
				    //DATARECORD << "A_B1R12["<< ma_tmp <<"]: "<< A_B1R12[ma_tmp]<<"\n";
				    //DATARECORD << "A_B1R13["<< ma_tmp <<"]: "<< A_B1R13[ma_tmp]<<"\n";
				    //DATARECORD << "A_B1R14["<< ma_tmp <<"]: "<< A_B1R14[ma_tmp]<<"\n";
				    //DATARECORD << "A_B1R15["<< ma_tmp <<"]: "<< A_B1R15[ma_tmp]<<"\n";					
					Radix16_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp], A_B1R2[ma_tmp], A_B1R3[ma_tmp],A_B1R4[ma_tmp],
							   A_B1R5[ma_tmp],A_B1R6[ma_tmp], A_B1R7[ma_tmp], A_B1R8[ma_tmp],A_B1R9[ma_tmp],
							   A_B1R10[ma_tmp],A_B1R11[ma_tmp],A_B1R12[ma_tmp],A_B1R13[ma_tmp],A_B1R14[ma_tmp],
							   A_B1R15[ma_tmp]);
					DATARECORD << "After Radix-16 butterfly unit operation!!! \n";
				    //DATARECORD << "A_B1R0["<< ma_tmp <<"]: " << A_B1R0[ma_tmp]<<"\n";
				    //DATARECORD << "A_B1R1["<< ma_tmp <<"]: " << A_B1R1[ma_tmp]<<"\n";
				    //DATARECORD << "A_B1R2["<< ma_tmp <<"]: " << A_B1R2[ma_tmp]<<"\n";
				    //DATARECORD << "A_B1R3["<< ma_tmp <<"]: " << A_B1R3[ma_tmp]<<"\n";
				    //DATARECORD << "A_B1R4["<< ma_tmp <<"]: " << A_B1R4[ma_tmp]<<"\n";
				    //DATARECORD << "A_B1R5["<< ma_tmp <<"]: " << A_B1R5[ma_tmp]<<"\n";
				    //DATARECORD << "A_B1R6["<< ma_tmp <<"]: " << A_B1R6[ma_tmp]<<"\n";
				    //DATARECORD << "A_B1R7["<< ma_tmp <<"]: " << A_B1R7[ma_tmp]<<"\n";
				    //DATARECORD << "A_B1R8["<< ma_tmp <<"]: " << A_B1R8[ma_tmp]<<"\n";
				    //DATARECORD << "A_B1R9["<< ma_tmp <<"]: " << A_B1R9[ma_tmp]<<"\n";
				    //DATARECORD << "A_B1R10["<< ma_tmp <<"]: "<< A_B1R10[ma_tmp]<<"\n";
				    //DATARECORD << "A_B1R11["<< ma_tmp <<"]: "<< A_B1R11[ma_tmp]<<"\n";
				    //DATARECORD << "A_B1R12["<< ma_tmp <<"]: "<< A_B1R12[ma_tmp]<<"\n";
				    //DATARECORD << "A_B1R13["<< ma_tmp <<"]: "<< A_B1R13[ma_tmp]<<"\n";
				    //DATARECORD << "A_B1R14["<< ma_tmp <<"]: "<< A_B1R14[ma_tmp]<<"\n";
				    //DATARECORD << "A_B1R15["<< ma_tmp <<"]: "<< A_B1R15[ma_tmp]<<"\n";							   
					MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],factor_t,p);
					MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],factor_2t,p);
					MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],factor_3t,p);
					MulMod(A_B1R4[ma_tmp],A_B1R4[ma_tmp],factor_4t,p);
					MulMod(A_B1R5[ma_tmp],A_B1R5[ma_tmp],factor_5t,p);
					MulMod(A_B1R6[ma_tmp],A_B1R6[ma_tmp],factor_6t,p);
					MulMod(A_B1R7[ma_tmp],A_B1R7[ma_tmp],factor_7t,p);
					MulMod(A_B1R8[ma_tmp],A_B1R8[ma_tmp],factor_8t,p);
					MulMod(A_B1R9[ma_tmp],A_B1R9[ma_tmp],factor_9t,p);
					MulMod(A_B1R10[ma_tmp],A_B1R10[ma_tmp],factor_10t,p);
					MulMod(A_B1R11[ma_tmp],A_B1R11[ma_tmp],factor_11t,p);
					MulMod(A_B1R12[ma_tmp],A_B1R12[ma_tmp],factor_12t,p);
					MulMod(A_B1R13[ma_tmp],A_B1R13[ma_tmp],factor_13t,p);
					MulMod(A_B1R14[ma_tmp],A_B1R14[ma_tmp],factor_14t,p);
					MulMod(A_B1R15[ma_tmp],A_B1R15[ma_tmp],factor_15t,p);
                    if(j <  2)bn1_ma_reg1 = ma_tmp;					
                    if((j >= 2)  && (j < 4))bn1_ma_reg2 = ma_tmp;
                    if((j >= 4)  && (j < 6))bn1_ma_reg3 = ma_tmp;
                    if((j >= 6)  && (j < 8))bn1_ma_reg4 = ma_tmp;
                    if((j >= 8)  && (j < 10))bn1_ma_reg5 = ma_tmp;
                    if((j >= 10) && (j < 12))bn1_ma_reg6 = ma_tmp;
                    if((j >= 12) && (j < 14))bn1_ma_reg7 = ma_tmp;
                    if((j >= 14) && (j < 16))bn1_ma_reg8 = ma_tmp;	
					//----------------DTFAG golden pattern------------------
					switch(s){
						case 0:
							DTFAG_golden_st0 << 1 		  << " \n ";
							DTFAG_golden_st0 << factor_t  << " \n ";
							DTFAG_golden_st0 << factor_2t << " \n ";
							DTFAG_golden_st0 << factor_3t << " \n ";
							DTFAG_golden_st0 << factor_4t << " \n ";
							DTFAG_golden_st0 << factor_5t  << " \n ";
							DTFAG_golden_st0 << factor_6t << " \n ";
							DTFAG_golden_st0 << factor_7t << " \n ";
							DTFAG_golden_st0 << factor_8t << " \n ";
							DTFAG_golden_st0 << factor_9t  << " \n ";
							DTFAG_golden_st0 << factor_10t << " \n ";
							DTFAG_golden_st0 << factor_11t << " \n ";
							DTFAG_golden_st0 << factor_12t << " \n ";
							DTFAG_golden_st0 << factor_13t  << " \n ";
							DTFAG_golden_st0 << factor_14t << " \n ";
							DTFAG_golden_st0 << factor_15t << " \n ";
							break;
						case 1:
							DTFAG_golden_st1 << 1 		  << " \n ";
							DTFAG_golden_st1 << factor_t  << " \n ";
							DTFAG_golden_st1 << factor_2t << " \n ";
							DTFAG_golden_st1 << factor_3t << " \n ";
							DTFAG_golden_st1 << factor_4t << " \n ";
							DTFAG_golden_st1 << factor_5t  << " \n ";
							DTFAG_golden_st1 << factor_6t << " \n ";
							DTFAG_golden_st1 << factor_7t << " \n ";
							DTFAG_golden_st1 << factor_8t << " \n ";
							DTFAG_golden_st1 << factor_9t  << " \n ";
							DTFAG_golden_st1 << factor_10t << " \n ";
							DTFAG_golden_st1 << factor_11t << " \n ";
							DTFAG_golden_st1 << factor_12t << " \n ";
							DTFAG_golden_st1 << factor_13t  << " \n ";
							DTFAG_golden_st1 << factor_14t << " \n ";
							DTFAG_golden_st1 << factor_15t << " \n ";
							break;
						case 2:
							DTFAG_golden_st2 << 1 		  << " \n ";
							DTFAG_golden_st2 << factor_t  << " \n ";
							DTFAG_golden_st2 << factor_2t << " \n ";
							DTFAG_golden_st2 << factor_3t << " \n ";
							DTFAG_golden_st2 << factor_4t << " \n ";
							DTFAG_golden_st2 << factor_5t  << " \n ";
							DTFAG_golden_st2 << factor_6t << " \n ";
							DTFAG_golden_st2 << factor_7t << " \n ";
							DTFAG_golden_st2 << factor_8t << " \n ";
							DTFAG_golden_st2 << factor_9t  << " \n ";
							DTFAG_golden_st2 << factor_10t << " \n ";
							DTFAG_golden_st2 << factor_11t << " \n ";
							DTFAG_golden_st2 << factor_12t << " \n ";
							DTFAG_golden_st2 << factor_13t  << " \n ";
							DTFAG_golden_st2 << factor_14t << " \n ";
							DTFAG_golden_st2 << factor_15t << " \n ";
							break;
						default:
							break;
					}
					//--------------------------------------------------				
				}
			}
		//data relocation
		 if(s < Stage-1){
		  if(bn1_bc_tmp > bn0_bc_tmp){
			 data_tmp_1  = A_B0R1[bn0_ma_reg1];
			 data_tmp_2  = A_B0R2[bn0_ma_reg1];
			 data_tmp_3  = A_B0R3[bn0_ma_reg1];
			 data_tmp_4  = A_B0R4[bn0_ma_reg1];
			 data_tmp_5  = A_B0R5[bn0_ma_reg1];
			 data_tmp_6  = A_B0R6[bn0_ma_reg1];
			 data_tmp_7  = A_B0R7[bn0_ma_reg1];
			 data_tmp_8  = A_B0R8[bn0_ma_reg1];
			 data_tmp_9  = A_B0R9[bn0_ma_reg1];
			 data_tmp_10 = A_B0R10[bn0_ma_reg1];
			 data_tmp_11 = A_B0R11[bn0_ma_reg1];
			 data_tmp_12 = A_B0R12[bn0_ma_reg1];
			 data_tmp_13 = A_B0R13[bn0_ma_reg1];
			 data_tmp_14 = A_B0R14[bn0_ma_reg1];
			 data_tmp_15 = A_B0R15[bn0_ma_reg1];
			 A_B0R1[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg1];
			 A_B0R2[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg2];
			 A_B0R3[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg2];
			 A_B0R4[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg3];
			 A_B0R5[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg3];
			 A_B0R6[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg1] = A_B0R0[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg1] = A_B1R0[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg1] = A_B0R0[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg1] = A_B1R0[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg1] = A_B1R0[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg1] = A_B0R0[bn0_ma_reg8];
			 A_B1R0[bn1_ma_reg1]  = data_tmp_1; 
			 A_B1R0[bn1_ma_reg2]  = data_tmp_2; 
			 A_B0R0[bn0_ma_reg2]  = data_tmp_3; 
			 A_B1R0[bn1_ma_reg3]  = data_tmp_4; 
			 A_B0R0[bn0_ma_reg3]  = data_tmp_5; 
			 A_B0R0[bn0_ma_reg4]  = data_tmp_6; 
			 A_B1R0[bn1_ma_reg4]  = data_tmp_7; 
			 A_B1R0[bn1_ma_reg5]  = data_tmp_8; 
			 A_B0R0[bn0_ma_reg5]  = data_tmp_9; 
			 A_B0R0[bn0_ma_reg6]  = data_tmp_10;
			 A_B1R0[bn1_ma_reg6]  = data_tmp_11;
			 A_B0R0[bn0_ma_reg7]  = data_tmp_12;
			 A_B1R0[bn1_ma_reg7]  = data_tmp_13;
			 A_B1R0[bn1_ma_reg8]  = data_tmp_14;
			 A_B0R0[bn0_ma_reg8]  = data_tmp_15;
			 /*********************************************/
			 data_tmp_1  = A_B1R2[bn1_ma_reg1];
			 data_tmp_2  = A_B1R3[bn1_ma_reg1];
			 data_tmp_3  = A_B1R4[bn1_ma_reg1];
			 data_tmp_4  = A_B1R5[bn1_ma_reg1];
			 data_tmp_5  = A_B1R6[bn1_ma_reg1];
			 data_tmp_6  = A_B1R7[bn1_ma_reg1];
			 data_tmp_7  = A_B1R8[bn1_ma_reg1];
			 data_tmp_8  = A_B1R9[bn1_ma_reg1];
			 data_tmp_9  = A_B1R10[bn1_ma_reg1];
			 data_tmp_10 = A_B1R11[bn1_ma_reg1];
			 data_tmp_11 = A_B1R12[bn1_ma_reg1];
			 data_tmp_12 = A_B1R13[bn1_ma_reg1];
			 data_tmp_13 = A_B1R14[bn1_ma_reg1];
			 data_tmp_14 = A_B1R15[bn1_ma_reg1];
             A_B1R2[bn1_ma_reg1]  =	 A_B1R1[bn1_ma_reg2];
             A_B1R3[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg2];
             A_B1R4[bn1_ma_reg1]  =  A_B1R1[bn1_ma_reg3];
             A_B1R5[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg3];
             A_B1R6[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg4];
             A_B1R7[bn1_ma_reg1]  =  A_B1R1[bn1_ma_reg4];
             A_B1R8[bn1_ma_reg1]  =  A_B1R1[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg1] =  A_B0R1[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg1] =  A_B1R1[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg1] =  A_B0R1[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg1] =  A_B1R1[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg1] =  A_B1R1[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg1] =  A_B0R1[bn0_ma_reg8];
			 A_B1R1[bn1_ma_reg2]  =  data_tmp_1; 
			 A_B0R1[bn0_ma_reg2]  =  data_tmp_2; 
			 A_B1R1[bn1_ma_reg3]  =  data_tmp_3; 
			 A_B0R1[bn0_ma_reg3]  =  data_tmp_4; 
			 A_B0R1[bn0_ma_reg4]  =  data_tmp_5; 
			 A_B1R1[bn1_ma_reg4]  =  data_tmp_6; 
			 A_B1R1[bn1_ma_reg5]  =  data_tmp_7; 
			 A_B0R1[bn0_ma_reg5]  =  data_tmp_8; 
			 A_B0R1[bn0_ma_reg6]  =  data_tmp_9; 
			 A_B1R1[bn1_ma_reg6]  =  data_tmp_10;
			 A_B0R1[bn0_ma_reg7]  =  data_tmp_11;
			 A_B1R1[bn1_ma_reg7]  =  data_tmp_12;
			 A_B1R1[bn1_ma_reg8]  =  data_tmp_13;
			 A_B0R1[bn0_ma_reg8]  =  data_tmp_14;
			/************************************************************/ 
			 data_tmp_1  =   A_B1R3[bn1_ma_reg2];
			 data_tmp_2  =   A_B1R4[bn1_ma_reg2];
			 data_tmp_3  =   A_B1R5[bn1_ma_reg2];
			 data_tmp_4  =   A_B1R6[bn1_ma_reg2];
			 data_tmp_5  =   A_B1R7[bn1_ma_reg2];
			 data_tmp_6  =   A_B1R8[bn1_ma_reg2];
			 data_tmp_7  =   A_B1R9[bn1_ma_reg2];
			 data_tmp_8  =   A_B1R10[bn1_ma_reg2];
			 data_tmp_9  =   A_B1R11[bn1_ma_reg2];
			 data_tmp_10 =   A_B1R12[bn1_ma_reg2];
			 data_tmp_11 =   A_B1R13[bn1_ma_reg2];
			 data_tmp_12 =   A_B1R14[bn1_ma_reg2];
			 data_tmp_13 =   A_B1R15[bn1_ma_reg2];
			 A_B1R3[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg2];
			 A_B1R4[bn1_ma_reg2]  =  A_B1R2[bn1_ma_reg3];
			 A_B1R5[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg3];
			 A_B1R6[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg4];
			 A_B1R7[bn1_ma_reg2]  =  A_B1R2[bn1_ma_reg4];
			 A_B1R8[bn1_ma_reg2]  =  A_B1R2[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg2] =  A_B0R2[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg2] =  A_B1R2[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg2] =  A_B0R2[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg2] =  A_B1R2[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg2] =  A_B1R2[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg2] =  A_B0R2[bn0_ma_reg8];
			 A_B0R2[bn0_ma_reg2]  =  data_tmp_1; 
			 A_B1R2[bn1_ma_reg3]  =  data_tmp_2; 
			 A_B0R2[bn0_ma_reg3]  =  data_tmp_3; 
			 A_B0R2[bn0_ma_reg4]  =  data_tmp_4; 
			 A_B1R2[bn1_ma_reg4]  =  data_tmp_5; 
			 A_B1R2[bn1_ma_reg5]  =  data_tmp_6; 
			 A_B0R2[bn0_ma_reg5]  =  data_tmp_7; 
			 A_B0R2[bn0_ma_reg6]  =  data_tmp_8; 
			 A_B1R2[bn1_ma_reg6]  =  data_tmp_9; 
			 A_B0R2[bn0_ma_reg7]  =  data_tmp_10;
			 A_B1R2[bn1_ma_reg7]  =  data_tmp_11;
			 A_B1R2[bn1_ma_reg8]  =  data_tmp_12;
			 A_B0R2[bn0_ma_reg8]  =  data_tmp_13;
			 //******************************************************
			 data_tmp_1  = A_B0R4[bn0_ma_reg2];
			 data_tmp_2  = A_B0R5[bn0_ma_reg2];
			 data_tmp_3  = A_B0R6[bn0_ma_reg2];
			 data_tmp_4  = A_B0R7[bn0_ma_reg2];
			 data_tmp_5  = A_B0R8[bn0_ma_reg2];
			 data_tmp_6  = A_B0R9[bn0_ma_reg2];
			 data_tmp_7  = A_B0R10[bn0_ma_reg2];
			 data_tmp_8  = A_B0R11[bn0_ma_reg2];
			 data_tmp_9  = A_B0R12[bn0_ma_reg2];
			 data_tmp_10 = A_B0R13[bn0_ma_reg2];
			 data_tmp_11 = A_B0R14[bn0_ma_reg2];
			 data_tmp_12 = A_B0R15[bn0_ma_reg2];
			 A_B0R4[bn0_ma_reg2]  =  A_B1R3[bn1_ma_reg3];
			 A_B0R5[bn0_ma_reg2]  =  A_B0R3[bn0_ma_reg3];
			 A_B0R6[bn0_ma_reg2]  =  A_B0R3[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg2]  =  A_B1R3[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg2]  =  A_B1R3[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg2]  =  A_B0R3[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg2] =  A_B0R3[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg2] =  A_B1R3[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg2] =  A_B0R3[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg2] =  A_B1R3[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg2] =  A_B1R3[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg2] =  A_B0R3[bn0_ma_reg8];
			 A_B1R3[bn1_ma_reg3]  =  data_tmp_1;  
			 A_B0R3[bn0_ma_reg3]  =  data_tmp_2;  
			 A_B0R3[bn0_ma_reg4]  =  data_tmp_3;  
			 A_B1R3[bn1_ma_reg4]  =  data_tmp_4;  
			 A_B1R3[bn1_ma_reg5]  =  data_tmp_5;  
			 A_B0R3[bn0_ma_reg5]  =  data_tmp_6;  
			 A_B0R3[bn0_ma_reg6]  =  data_tmp_7;  
			 A_B1R3[bn1_ma_reg6]  =  data_tmp_8;  
			 A_B0R3[bn0_ma_reg7]  =  data_tmp_9;  
			 A_B1R3[bn1_ma_reg7]  =  data_tmp_10; 
			 A_B1R3[bn1_ma_reg8]  =  data_tmp_11; 
			 A_B0R3[bn0_ma_reg8]  =  data_tmp_12; 
			 //----------------------------------------------------------------------
             data_tmp_1  = A_B1R5[bn1_ma_reg3];
			 data_tmp_2  = A_B1R6[bn1_ma_reg3];
			 data_tmp_3  = A_B1R7[bn1_ma_reg3];
			 data_tmp_4  = A_B1R8[bn1_ma_reg3];
			 data_tmp_5  = A_B1R9[bn1_ma_reg3];
			 data_tmp_6  = A_B1R10[bn1_ma_reg3];
			 data_tmp_7  = A_B1R11[bn1_ma_reg3];
			 data_tmp_8  = A_B1R12[bn1_ma_reg3];
			 data_tmp_9  = A_B1R13[bn1_ma_reg3];
			 data_tmp_10 = A_B1R14[bn1_ma_reg3];
			 data_tmp_11 = A_B1R15[bn1_ma_reg3];
			 A_B1R5[bn1_ma_reg3]  =  A_B0R4[bn0_ma_reg3];
			 A_B1R6[bn1_ma_reg3]  =  A_B0R4[bn0_ma_reg4];
			 A_B1R7[bn1_ma_reg3]  =  A_B1R4[bn1_ma_reg4];
			 A_B1R8[bn1_ma_reg3]  =  A_B1R4[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg3]  =  A_B0R4[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg3] =  A_B0R4[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg3] =  A_B1R4[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg3] =  A_B0R4[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg3] =  A_B1R4[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg3] =  A_B1R4[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg3] =  A_B0R4[bn0_ma_reg8];
			 A_B0R4[bn0_ma_reg3]  =  data_tmp_1 ;
			 A_B0R4[bn0_ma_reg4]  =  data_tmp_2 ;
			 A_B1R4[bn1_ma_reg4]  =  data_tmp_3 ;
			 A_B1R4[bn1_ma_reg5]  =  data_tmp_4 ;
			 A_B0R4[bn0_ma_reg5]  =  data_tmp_5 ;
			 A_B0R4[bn0_ma_reg6]  =  data_tmp_6 ;
			 A_B1R4[bn1_ma_reg6]  =  data_tmp_7 ;
			 A_B0R4[bn0_ma_reg7]  =  data_tmp_8 ;
			 A_B1R4[bn1_ma_reg7]  =  data_tmp_9 ;
			 A_B1R4[bn1_ma_reg8]  =  data_tmp_10;
			 A_B0R4[bn0_ma_reg8]  =  data_tmp_11;
			 //************************************************************************
			 data_tmp_1  = A_B0R6[bn0_ma_reg3];
			 data_tmp_2  = A_B0R7[bn0_ma_reg3];
			 data_tmp_3  = A_B0R8[bn0_ma_reg3];
			 data_tmp_4  = A_B0R9[bn0_ma_reg3];
			 data_tmp_5  = A_B0R10[bn0_ma_reg3];
			 data_tmp_6  = A_B0R11[bn0_ma_reg3];
			 data_tmp_7  = A_B0R12[bn0_ma_reg3];
			 data_tmp_8  = A_B0R13[bn0_ma_reg3];
			 data_tmp_9  = A_B0R14[bn0_ma_reg3];
			 data_tmp_10 = A_B0R15[bn0_ma_reg3];
			 A_B0R6[bn0_ma_reg3]  =  A_B0R5[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg3]  =  A_B1R5[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg3]  =  A_B1R5[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg3]  =  A_B0R5[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg3] =  A_B0R5[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg3] =  A_B1R5[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg3] =  A_B0R5[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg3] =  A_B1R5[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg3] =  A_B1R5[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg3] =  A_B0R5[bn0_ma_reg8];
			 A_B0R5[bn0_ma_reg4]  =  data_tmp_1; 
			 A_B1R5[bn1_ma_reg4]  =  data_tmp_2; 
			 A_B1R5[bn1_ma_reg5]  =  data_tmp_3; 
			 A_B0R5[bn0_ma_reg5]  =  data_tmp_4; 
			 A_B0R5[bn0_ma_reg6]  =  data_tmp_5; 
			 A_B1R5[bn1_ma_reg6]  =  data_tmp_6; 
			 A_B0R5[bn0_ma_reg7]  =  data_tmp_7; 
			 A_B1R5[bn1_ma_reg7]  =  data_tmp_8; 
			 A_B1R5[bn1_ma_reg8]  =  data_tmp_9; 
			 A_B0R5[bn0_ma_reg8]  =  data_tmp_10;
			 //-----------------------------------------------------------------------
			 data_tmp_1  = A_B0R7[bn0_ma_reg4];
			 data_tmp_2  = A_B0R8[bn0_ma_reg4];
			 data_tmp_3  = A_B0R9[bn0_ma_reg4];
			 data_tmp_4  = A_B0R10[bn0_ma_reg4];
			 data_tmp_5  = A_B0R11[bn0_ma_reg4];
			 data_tmp_6  = A_B0R12[bn0_ma_reg4];
			 data_tmp_7  = A_B0R13[bn0_ma_reg4];
			 data_tmp_8  = A_B0R14[bn0_ma_reg4];
			 data_tmp_9  = A_B0R15[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg4]  =  A_B1R6[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg4]  =  A_B1R6[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg4]  =  A_B0R6[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg4] =  A_B0R6[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg4] =  A_B1R6[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg4] =  A_B0R6[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg4] =  A_B1R6[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg4] =  A_B1R6[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg4] =  A_B0R6[bn0_ma_reg8];
			 A_B1R6[bn1_ma_reg4]  =  data_tmp_1;
			 A_B1R6[bn1_ma_reg5]  =  data_tmp_2;
			 A_B0R6[bn0_ma_reg5]  =  data_tmp_3;
			 A_B0R6[bn0_ma_reg6]  =  data_tmp_4;
			 A_B1R6[bn1_ma_reg6]  =  data_tmp_5;
			 A_B0R6[bn0_ma_reg7]  =  data_tmp_6;
			 A_B1R6[bn1_ma_reg7]  =  data_tmp_7;
			 A_B1R6[bn1_ma_reg8]  =  data_tmp_8;
			 A_B0R6[bn0_ma_reg8]  =  data_tmp_9;
			 //----------------------------------------------------------------------
			 data_tmp_1  = A_B1R8[bn1_ma_reg4];
			 data_tmp_2  = A_B1R9[bn1_ma_reg4];
			 data_tmp_3  = A_B1R10[bn1_ma_reg4];
			 data_tmp_4  = A_B1R11[bn1_ma_reg4];
			 data_tmp_5  = A_B1R12[bn1_ma_reg4];
			 data_tmp_6  = A_B1R13[bn1_ma_reg4];
			 data_tmp_7  = A_B1R14[bn1_ma_reg4];
			 data_tmp_8  = A_B1R15[bn1_ma_reg4];
			 A_B1R8[bn1_ma_reg4]  =  A_B1R7[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg4]  =  A_B0R7[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg4] =  A_B0R7[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg4] =  A_B1R7[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg4] =  A_B0R7[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg4] =  A_B1R7[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg4] =  A_B1R7[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg4] =  A_B0R7[bn0_ma_reg8];
			 A_B1R7[bn1_ma_reg5]  =  data_tmp_1;
			 A_B0R7[bn0_ma_reg5]  =  data_tmp_2;
			 A_B0R7[bn0_ma_reg6]  =  data_tmp_3;
			 A_B1R7[bn1_ma_reg6]  =  data_tmp_4;
			 A_B0R7[bn0_ma_reg7]  =  data_tmp_5;
			 A_B1R7[bn1_ma_reg7]  =  data_tmp_6;
			 A_B1R7[bn1_ma_reg8]  =  data_tmp_7;
			 A_B0R7[bn0_ma_reg8]  =  data_tmp_8;
			 //----------------------------------------------------------------------
			 data_tmp_1 = A_B1R9[bn1_ma_reg5];
			 data_tmp_2 = A_B1R10[bn1_ma_reg5];
			 data_tmp_3 = A_B1R11[bn1_ma_reg5];
			 data_tmp_4 = A_B1R12[bn1_ma_reg5];
			 data_tmp_5 = A_B1R13[bn1_ma_reg5];
			 data_tmp_6 = A_B1R14[bn1_ma_reg5];
			 data_tmp_7 = A_B1R15[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg5]  = A_B0R8[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg5] = A_B0R8[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg5] = A_B1R8[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg5] = A_B0R8[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg5] = A_B1R8[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg5] = A_B1R8[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg5] = A_B0R8[bn0_ma_reg8];
			 A_B0R8[bn0_ma_reg5]  = data_tmp_1;
			 A_B0R8[bn0_ma_reg6]  = data_tmp_2;
			 A_B1R8[bn1_ma_reg6]  = data_tmp_3;
			 A_B0R8[bn0_ma_reg7]  = data_tmp_4;
			 A_B1R8[bn1_ma_reg7]  = data_tmp_5;
			 A_B1R8[bn1_ma_reg8]  = data_tmp_6;
			 A_B0R8[bn0_ma_reg8]  = data_tmp_7;
			 //---------------------------------------------------------------------
			 data_tmp_1  = A_B0R10[bn0_ma_reg5];
			 data_tmp_2  = A_B0R11[bn0_ma_reg5];
			 data_tmp_3  = A_B0R12[bn0_ma_reg5];
			 data_tmp_4  = A_B0R13[bn0_ma_reg5];
			 data_tmp_5  = A_B0R14[bn0_ma_reg5];
			 data_tmp_6  = A_B0R15[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg5] = A_B0R9[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg5] = A_B1R9[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg5] = A_B0R9[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg5] = A_B1R9[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg5] = A_B1R9[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg5] = A_B0R9[bn0_ma_reg8];
			 A_B0R9[bn0_ma_reg6]  = data_tmp_1;
			 A_B1R9[bn1_ma_reg6]  = data_tmp_2;
			 A_B0R9[bn0_ma_reg7]  = data_tmp_3;
			 A_B1R9[bn1_ma_reg7]  = data_tmp_4;
			 A_B1R9[bn1_ma_reg8]  = data_tmp_5;
			 A_B0R9[bn0_ma_reg8]  = data_tmp_6;
			 //--------------------------------------------------------------------
			 data_tmp_1  = A_B0R11[bn0_ma_reg6];
			 data_tmp_2  = A_B0R12[bn0_ma_reg6];
			 data_tmp_3  = A_B0R13[bn0_ma_reg6];
			 data_tmp_4  = A_B0R14[bn0_ma_reg6];
			 data_tmp_5  = A_B0R15[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg6] = A_B1R10[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg6] = A_B0R10[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg6] = A_B1R10[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg6] = A_B1R10[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg6] = A_B0R10[bn0_ma_reg8];
			 A_B1R10[bn1_ma_reg6] = data_tmp_1;
			 A_B0R10[bn0_ma_reg7] = data_tmp_2;
			 A_B1R10[bn1_ma_reg7] = data_tmp_3;
			 A_B1R10[bn1_ma_reg8] = data_tmp_4;
			 A_B0R10[bn0_ma_reg8] = data_tmp_5;
			 //--------------------------------------------------------------------
			 data_tmp_1  = A_B1R12[bn1_ma_reg6];
			 data_tmp_2  = A_B1R13[bn1_ma_reg6];
			 data_tmp_3  = A_B1R14[bn1_ma_reg6];
			 data_tmp_4  = A_B1R15[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg6] = A_B0R11[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg6] = A_B1R11[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg6] = A_B1R11[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg6] = A_B0R11[bn0_ma_reg8];
			 A_B0R11[bn0_ma_reg7] = data_tmp_1;
			 A_B1R11[bn1_ma_reg7] = data_tmp_2;
			 A_B1R11[bn1_ma_reg8] = data_tmp_3;
			 A_B0R11[bn0_ma_reg8] = data_tmp_4;
			 //--------------------------------------------------------------------
			 data_tmp_1 = A_B0R13[bn0_ma_reg7];
			 data_tmp_2 = A_B0R14[bn0_ma_reg7];
			 data_tmp_3 = A_B0R15[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg7] = A_B1R12[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg7] = A_B1R12[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg7] = A_B0R12[bn0_ma_reg8];
			 A_B1R12[bn1_ma_reg7] = data_tmp_1;
			 A_B1R12[bn1_ma_reg8] = data_tmp_2;
			 A_B0R12[bn0_ma_reg8] = data_tmp_3;
			 //--------------------------------------------------------------------
			 data_tmp_1 = A_B1R14[bn1_ma_reg7];
			 data_tmp_2 = A_B1R15[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg7] = A_B1R13[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg7] = A_B0R13[bn0_ma_reg8];
			 A_B1R13[bn1_ma_reg8] = data_tmp_1;
			 A_B0R13[bn0_ma_reg8] = data_tmp_2;
			 //--------------------------------------------------------------------
			 data_tmp_1 = A_B1R15[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg8] = A_B0R14[bn0_ma_reg8];
			 A_B0R14[bn0_ma_reg8] = data_tmp_1; 
		  }else{	 
			 data_tmp_1  = A_B1R1[bn1_ma_reg1];
			 data_tmp_2  = A_B1R2[bn1_ma_reg1];
			 data_tmp_3  = A_B1R3[bn1_ma_reg1];
			 data_tmp_4  = A_B1R4[bn1_ma_reg1];
			 data_tmp_5  = A_B1R5[bn1_ma_reg1];
			 data_tmp_6  = A_B1R6[bn1_ma_reg1];
			 data_tmp_7  = A_B1R7[bn1_ma_reg1];
			 data_tmp_8  = A_B1R8[bn1_ma_reg1];
			 data_tmp_9  = A_B1R9[bn1_ma_reg1];
			 data_tmp_10 = A_B1R10[bn1_ma_reg1];
			 data_tmp_11 = A_B1R11[bn1_ma_reg1];
			 data_tmp_12 = A_B1R12[bn1_ma_reg1];
			 data_tmp_13 = A_B1R13[bn1_ma_reg1];
			 data_tmp_14 = A_B1R14[bn1_ma_reg1];
			 data_tmp_15 = A_B1R15[bn1_ma_reg1];
             A_B1R1[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg1];
             A_B1R2[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg2];
             A_B1R3[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg2];
             A_B1R4[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg3];
             A_B1R5[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg3];
             A_B1R6[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg4];
             A_B1R7[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg4];
             A_B1R8[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg5];
             A_B1R9[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg5];
             A_B1R10[bn1_ma_reg1] =  A_B1R0[bn1_ma_reg6];
             A_B1R11[bn1_ma_reg1] =  A_B0R0[bn0_ma_reg6];
             A_B1R12[bn1_ma_reg1] =  A_B1R0[bn1_ma_reg7];
             A_B1R13[bn1_ma_reg1] =  A_B0R0[bn0_ma_reg7];
             A_B1R14[bn1_ma_reg1] =  A_B0R0[bn0_ma_reg8];
             A_B1R15[bn1_ma_reg1] =  A_B1R0[bn1_ma_reg8];
             A_B0R0[bn0_ma_reg1]  =  data_tmp_1; 
             A_B0R0[bn0_ma_reg2]  =  data_tmp_2; 
             A_B1R0[bn1_ma_reg2]  =  data_tmp_3; 
             A_B0R0[bn0_ma_reg3]  =  data_tmp_4; 
             A_B1R0[bn1_ma_reg3]  =  data_tmp_5; 
             A_B1R0[bn1_ma_reg4]  =  data_tmp_6; 
             A_B0R0[bn0_ma_reg4]  =  data_tmp_7; 
             A_B0R0[bn0_ma_reg5]  =  data_tmp_8; 
             A_B1R0[bn1_ma_reg5]  =  data_tmp_9; 
             A_B1R0[bn1_ma_reg6]  =  data_tmp_10;
             A_B0R0[bn0_ma_reg6]  =  data_tmp_11;
             A_B1R0[bn1_ma_reg7]  =  data_tmp_12;
             A_B0R0[bn0_ma_reg7]  =  data_tmp_13;
             A_B0R0[bn0_ma_reg8]  =  data_tmp_14;
             A_B1R0[bn1_ma_reg8]  =  data_tmp_15;
             //-----------------------------------------------------------
             data_tmp_1  = A_B0R2[bn0_ma_reg1];
			 data_tmp_2  = A_B0R3[bn0_ma_reg1];
			 data_tmp_3  = A_B0R4[bn0_ma_reg1];
			 data_tmp_4  = A_B0R5[bn0_ma_reg1];
			 data_tmp_5  = A_B0R6[bn0_ma_reg1];
			 data_tmp_6  = A_B0R7[bn0_ma_reg1];
			 data_tmp_7  = A_B0R8[bn0_ma_reg1];
			 data_tmp_8  = A_B0R9[bn0_ma_reg1];
			 data_tmp_9  = A_B0R10[bn0_ma_reg1];
			 data_tmp_10 = A_B0R11[bn0_ma_reg1];
			 data_tmp_11 = A_B0R12[bn0_ma_reg1];
			 data_tmp_12 = A_B0R13[bn0_ma_reg1];
			 data_tmp_13 = A_B0R14[bn0_ma_reg1];
			 data_tmp_14 = A_B0R15[bn0_ma_reg1];
			 A_B0R2[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg2];
			 A_B0R3[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg2];
			 A_B0R4[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg3];
			 A_B0R5[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg3];
			 A_B0R6[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg4];
			 A_B0R7[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg1] =  A_B1R1[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg1] =  A_B0R1[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg1] =  A_B1R1[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg1] =  A_B0R1[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg1] =  A_B0R1[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg1] =  A_B1R1[bn1_ma_reg8];
             A_B0R1[bn0_ma_reg2]  =  data_tmp_1; 
             A_B1R1[bn1_ma_reg2]  =  data_tmp_2; 
             A_B0R1[bn0_ma_reg3]  =  data_tmp_3; 
             A_B1R1[bn1_ma_reg3]  =  data_tmp_4; 
             A_B1R1[bn1_ma_reg4]  =  data_tmp_5; 
             A_B0R1[bn0_ma_reg4]  =  data_tmp_6; 
             A_B0R1[bn0_ma_reg5]  =  data_tmp_7; 
             A_B1R1[bn1_ma_reg5]  =  data_tmp_8; 
             A_B1R1[bn1_ma_reg6]  =  data_tmp_9; 
             A_B0R1[bn0_ma_reg6]  =  data_tmp_10;
             A_B1R1[bn1_ma_reg7]  =  data_tmp_11;
             A_B0R1[bn0_ma_reg7]  =  data_tmp_12;
             A_B0R1[bn0_ma_reg8]  =  data_tmp_13;
             A_B1R1[bn1_ma_reg8]  =  data_tmp_14;
             //------------------------------------------------------------
             data_tmp_1  =  A_B0R3[bn0_ma_reg2];
			 data_tmp_2  =  A_B0R4[bn0_ma_reg2];
			 data_tmp_3  =  A_B0R5[bn0_ma_reg2];
			 data_tmp_4  =  A_B0R6[bn0_ma_reg2];
			 data_tmp_5  =  A_B0R7[bn0_ma_reg2];
			 data_tmp_6  =  A_B0R8[bn0_ma_reg2];
			 data_tmp_7  =  A_B0R9[bn0_ma_reg2];
			 data_tmp_8  =  A_B0R10[bn0_ma_reg2];
			 data_tmp_9  =  A_B0R11[bn0_ma_reg2];
			 data_tmp_10 =  A_B0R12[bn0_ma_reg2];
			 data_tmp_11 =  A_B0R13[bn0_ma_reg2];
			 data_tmp_12 =  A_B0R14[bn0_ma_reg2];
			 data_tmp_13 =  A_B0R15[bn0_ma_reg2];
			 A_B0R3[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg2];
			 A_B0R4[bn0_ma_reg2]  =  A_B0R2[bn0_ma_reg3];
			 A_B0R5[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg3];
			 A_B0R6[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg4];
			 A_B0R7[bn0_ma_reg2]  =  A_B0R2[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg2]  =  A_B0R2[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg2] =  A_B1R2[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg2] =  A_B0R2[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg2] =  A_B1R2[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg2] =  A_B0R2[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg2] =  A_B0R2[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg2] =  A_B1R2[bn1_ma_reg8];
			 A_B1R2[bn1_ma_reg2]  =  data_tmp_1; 
			 A_B0R2[bn0_ma_reg3]  =  data_tmp_2; 
			 A_B1R2[bn1_ma_reg3]  =  data_tmp_3; 
			 A_B1R2[bn1_ma_reg4]  =  data_tmp_4; 
			 A_B0R2[bn0_ma_reg4]  =  data_tmp_5; 
			 A_B0R2[bn0_ma_reg5]  =  data_tmp_6; 
			 A_B1R2[bn1_ma_reg5]  =  data_tmp_7; 
			 A_B1R2[bn1_ma_reg6]  =  data_tmp_8; 
			 A_B0R2[bn0_ma_reg6]  =  data_tmp_9; 
			 A_B1R2[bn1_ma_reg7]  =  data_tmp_10;
			 A_B0R2[bn0_ma_reg7]  =  data_tmp_11;
			 A_B0R2[bn0_ma_reg8]  =  data_tmp_12;
			 A_B1R2[bn1_ma_reg8]  =  data_tmp_13;
			 //-----------------------------------------------------------                        
			 data_tmp_1  = A_B1R4[bn1_ma_reg2];
			 data_tmp_2  = A_B1R5[bn1_ma_reg2];
			 data_tmp_3  = A_B1R6[bn1_ma_reg2];
			 data_tmp_4  = A_B1R7[bn1_ma_reg2];
			 data_tmp_5  = A_B1R8[bn1_ma_reg2];
			 data_tmp_6  = A_B1R9[bn1_ma_reg2];
			 data_tmp_7  = A_B1R10[bn1_ma_reg2];
			 data_tmp_8  = A_B1R11[bn1_ma_reg2];
			 data_tmp_9  = A_B1R12[bn1_ma_reg2];
			 data_tmp_10 = A_B1R13[bn1_ma_reg2];
			 data_tmp_11 = A_B1R14[bn1_ma_reg2];
			 data_tmp_12 = A_B1R15[bn1_ma_reg2];
			 A_B1R4[bn1_ma_reg2]  = A_B0R3[bn0_ma_reg3];
			 A_B1R5[bn1_ma_reg2]  = A_B1R3[bn1_ma_reg3];
			 A_B1R6[bn1_ma_reg2]  = A_B1R3[bn1_ma_reg4];
			 A_B1R7[bn1_ma_reg2]  = A_B0R3[bn0_ma_reg4];
			 A_B1R8[bn1_ma_reg2]  = A_B0R3[bn0_ma_reg5];
			 A_B1R9[bn1_ma_reg2]  = A_B1R3[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg2] = A_B1R3[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg2] = A_B0R3[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg2] = A_B1R3[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg2] = A_B0R3[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg2] = A_B0R3[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg2] = A_B1R3[bn1_ma_reg8];
			 A_B0R3[bn0_ma_reg3]  = data_tmp_1; 
			 A_B1R3[bn1_ma_reg3]  = data_tmp_2; 
			 A_B1R3[bn1_ma_reg4]  = data_tmp_3; 
			 A_B0R3[bn0_ma_reg4]  = data_tmp_4; 
			 A_B0R3[bn0_ma_reg5]  = data_tmp_5; 
			 A_B1R3[bn1_ma_reg5]  = data_tmp_6; 
			 A_B1R3[bn1_ma_reg6]  = data_tmp_7; 
			 A_B0R3[bn0_ma_reg6]  = data_tmp_8; 
			 A_B1R3[bn1_ma_reg7]  = data_tmp_9; 
			 A_B0R3[bn0_ma_reg7]  = data_tmp_10;
			 A_B0R3[bn0_ma_reg8]  = data_tmp_11;
			 A_B1R3[bn1_ma_reg8]  = data_tmp_12;
			 //------------------------------------------------------------
			 data_tmp_1  =  A_B0R5[bn0_ma_reg3];
			 data_tmp_2  =  A_B0R6[bn0_ma_reg3];
			 data_tmp_3  =  A_B0R7[bn0_ma_reg3];
			 data_tmp_4  =  A_B0R8[bn0_ma_reg3];
			 data_tmp_5  =  A_B0R9[bn0_ma_reg3];
			 data_tmp_6  =  A_B0R10[bn0_ma_reg3];
			 data_tmp_7  =  A_B0R11[bn0_ma_reg3];
			 data_tmp_8  =  A_B0R12[bn0_ma_reg3];
			 data_tmp_9  =  A_B0R13[bn0_ma_reg3];
			 data_tmp_10 =  A_B0R14[bn0_ma_reg3];
			 data_tmp_11 =  A_B0R15[bn0_ma_reg3];
			 A_B0R5[bn0_ma_reg3]  =  A_B1R4[bn1_ma_reg3];
			 A_B0R6[bn0_ma_reg3]  =  A_B1R4[bn1_ma_reg4];
			 A_B0R7[bn0_ma_reg3]  =  A_B0R4[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg3]  =  A_B0R4[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg3]  =  A_B1R4[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg3] =  A_B1R4[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg3] =  A_B0R4[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg3] =  A_B1R4[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg3] =  A_B0R4[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg3] =  A_B0R4[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg3] =  A_B1R4[bn1_ma_reg8];
			 A_B1R4[bn1_ma_reg3]  =  data_tmp_1; 
			 A_B1R4[bn1_ma_reg4]  =  data_tmp_2; 
			 A_B0R4[bn0_ma_reg4]  =  data_tmp_3; 
			 A_B0R4[bn0_ma_reg5]  =  data_tmp_4; 
			 A_B1R4[bn1_ma_reg5]  =  data_tmp_5; 
			 A_B1R4[bn1_ma_reg6]  =  data_tmp_6; 
			 A_B0R4[bn0_ma_reg6]  =  data_tmp_7; 
			 A_B1R4[bn1_ma_reg7]  =  data_tmp_8; 
			 A_B0R4[bn0_ma_reg7]  =  data_tmp_9; 
			 A_B0R4[bn0_ma_reg8]  =  data_tmp_10;
			 A_B1R4[bn1_ma_reg8]  =  data_tmp_11;
			 //-------------------------------------------------------------
			 data_tmp_1  =  A_B1R6[bn1_ma_reg3];
			 data_tmp_2  =  A_B1R7[bn1_ma_reg3];
			 data_tmp_3  =  A_B1R8[bn1_ma_reg3];
			 data_tmp_4  =  A_B1R9[bn1_ma_reg3];
			 data_tmp_5  =  A_B1R10[bn1_ma_reg3];
			 data_tmp_6  =  A_B1R11[bn1_ma_reg3];
			 data_tmp_7  =  A_B1R12[bn1_ma_reg3];
			 data_tmp_8  =  A_B1R13[bn1_ma_reg3];
			 data_tmp_9  =  A_B1R14[bn1_ma_reg3];
			 data_tmp_10 =  A_B1R15[bn1_ma_reg3];
			 A_B1R6[bn1_ma_reg3]  =  A_B1R5[bn1_ma_reg4];
			 A_B1R7[bn1_ma_reg3]  =  A_B0R5[bn0_ma_reg4];
			 A_B1R8[bn1_ma_reg3]  =  A_B0R5[bn0_ma_reg5];
			 A_B1R9[bn1_ma_reg3]  =  A_B1R5[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg3] =  A_B1R5[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg3] =  A_B0R5[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg3] =  A_B1R5[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg3] =  A_B0R5[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg3] =  A_B0R5[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg3] =  A_B1R5[bn1_ma_reg8];
			 A_B1R5[bn1_ma_reg4]  =  data_tmp_1; 
			 A_B0R5[bn0_ma_reg4]  =  data_tmp_2; 
			 A_B0R5[bn0_ma_reg5]  =  data_tmp_3; 
			 A_B1R5[bn1_ma_reg5]  =  data_tmp_4; 
			 A_B1R5[bn1_ma_reg6]  =  data_tmp_5; 
			 A_B0R5[bn0_ma_reg6]  =  data_tmp_6; 
			 A_B1R5[bn1_ma_reg7]  =  data_tmp_7; 
			 A_B0R5[bn0_ma_reg7]  =  data_tmp_8; 
			 A_B0R5[bn0_ma_reg8]  =  data_tmp_9; 
			 A_B1R5[bn1_ma_reg8]  =  data_tmp_10;
			 //----------------------------------------------------------------
			 data_tmp_1  =  A_B1R7[bn1_ma_reg4];
			 data_tmp_2  =  A_B1R8[bn1_ma_reg4];
			 data_tmp_3  =  A_B1R9[bn1_ma_reg4];
			 data_tmp_4  =  A_B1R10[bn1_ma_reg4];
			 data_tmp_5  =  A_B1R11[bn1_ma_reg4];
			 data_tmp_6  =  A_B1R12[bn1_ma_reg4];
			 data_tmp_7  =  A_B1R13[bn1_ma_reg4];
			 data_tmp_8  =  A_B1R14[bn1_ma_reg4];
			 data_tmp_9  =  A_B1R15[bn1_ma_reg4];
			 A_B1R7[bn1_ma_reg4]  =  A_B0R6[bn0_ma_reg4];
			 A_B1R8[bn1_ma_reg4]  =  A_B0R6[bn0_ma_reg5];
			 A_B1R9[bn1_ma_reg4]  =  A_B1R6[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg4] =  A_B1R6[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg4] =  A_B0R6[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg4] =  A_B1R6[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg4] =  A_B0R6[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg4] =  A_B0R6[bn0_ma_reg8];
             A_B1R15[bn1_ma_reg4] =  A_B1R6[bn1_ma_reg8];
             A_B0R6[bn0_ma_reg4]  =  data_tmp_1;
             A_B0R6[bn0_ma_reg5]  =  data_tmp_2;
             A_B1R6[bn1_ma_reg5]  =  data_tmp_3;
             A_B1R6[bn1_ma_reg6]  =  data_tmp_4;
             A_B0R6[bn0_ma_reg6]  =  data_tmp_5;
             A_B1R6[bn1_ma_reg7]  =  data_tmp_6;
             A_B0R6[bn0_ma_reg7]  =  data_tmp_7;
             A_B0R6[bn0_ma_reg8]  =  data_tmp_8;
             A_B1R6[bn1_ma_reg8]  =  data_tmp_9;
             //------------------------------------------------------------------
			 data_tmp_1  = A_B0R8[bn0_ma_reg4];
			 data_tmp_2  = A_B0R9[bn0_ma_reg4];
			 data_tmp_3  = A_B0R10[bn0_ma_reg4];
			 data_tmp_4  = A_B0R11[bn0_ma_reg4];
			 data_tmp_5  = A_B0R12[bn0_ma_reg4];
			 data_tmp_6  = A_B0R13[bn0_ma_reg4];
			 data_tmp_7  = A_B0R14[bn0_ma_reg4];
			 data_tmp_8  = A_B0R15[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg4]  =  A_B0R7[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg4]  =  A_B1R7[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg4] =  A_B1R7[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg4] =  A_B0R7[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg4] =  A_B1R7[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg4] =  A_B0R7[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg4] =  A_B0R7[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg4] =  A_B1R7[bn1_ma_reg8];
			 A_B0R7[bn0_ma_reg5]  =  data_tmp_1;
			 A_B1R7[bn1_ma_reg5]  =  data_tmp_2;
			 A_B1R7[bn1_ma_reg6]  =  data_tmp_3;
			 A_B0R7[bn0_ma_reg6]  =  data_tmp_4;
			 A_B1R7[bn1_ma_reg7]  =  data_tmp_5;
			 A_B0R7[bn0_ma_reg7]  =  data_tmp_6;
			 A_B0R7[bn0_ma_reg8]  =  data_tmp_7;
			 A_B1R7[bn1_ma_reg8]  =  data_tmp_8;
			 //----------------------------------------------------------------
			 data_tmp_1  = A_B0R9[bn0_ma_reg5];
			 data_tmp_2  = A_B0R10[bn0_ma_reg5];
			 data_tmp_3  = A_B0R11[bn0_ma_reg5];
			 data_tmp_4  = A_B0R12[bn0_ma_reg5];
			 data_tmp_5  = A_B0R13[bn0_ma_reg5];
			 data_tmp_6  = A_B0R14[bn0_ma_reg5];
			 data_tmp_7  = A_B0R15[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg5]  =  A_B1R8[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg8];
			 A_B1R8[bn1_ma_reg5]  =  data_tmp_1;
			 A_B1R8[bn1_ma_reg6]  =  data_tmp_2;
			 A_B0R8[bn0_ma_reg6]  =  data_tmp_3;
			 A_B1R8[bn1_ma_reg7]  =  data_tmp_4;
			 A_B0R8[bn0_ma_reg7]  =  data_tmp_5;
			 A_B0R8[bn0_ma_reg8]  =  data_tmp_6;
			 A_B1R8[bn1_ma_reg8]  =  data_tmp_7;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B1R10[bn1_ma_reg5];
			 data_tmp_2 = A_B1R11[bn1_ma_reg5];
			 data_tmp_3 = A_B1R12[bn1_ma_reg5];
			 data_tmp_4 = A_B1R13[bn1_ma_reg5];
			 data_tmp_5 = A_B1R14[bn1_ma_reg5];
			 data_tmp_6 = A_B1R15[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg5] = A_B1R9[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg5] = A_B0R9[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg5] = A_B1R9[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg5] = A_B0R9[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg5] = A_B0R9[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg5] = A_B1R9[bn1_ma_reg8];
			 A_B1R9[bn1_ma_reg6]  = data_tmp_1;
			 A_B0R9[bn0_ma_reg6]  = data_tmp_2;
			 A_B1R9[bn1_ma_reg7]  = data_tmp_3;
			 A_B0R9[bn0_ma_reg7]  = data_tmp_4;
			 A_B0R9[bn0_ma_reg8]  = data_tmp_5;
			 A_B1R9[bn1_ma_reg8]  = data_tmp_6;
			 //----------------------------------------------------------------
			 data_tmp_1 =  A_B1R11[bn1_ma_reg6];
			 data_tmp_2 =  A_B1R12[bn1_ma_reg6];
			 data_tmp_3 =  A_B1R13[bn1_ma_reg6];
			 data_tmp_4 =  A_B1R14[bn1_ma_reg6];
			 data_tmp_5 =  A_B1R15[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg6]  = A_B0R10[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg6]  = A_B1R10[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg6]  = A_B0R10[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg6]  = A_B0R10[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg6]  = A_B1R10[bn1_ma_reg8];
			 A_B0R10[bn0_ma_reg6]  = data_tmp_1;
			 A_B1R10[bn1_ma_reg7]  = data_tmp_2;
			 A_B0R10[bn0_ma_reg7]  = data_tmp_3;
			 A_B0R10[bn0_ma_reg8]  = data_tmp_4;
			 A_B1R10[bn1_ma_reg8]  = data_tmp_5;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B0R12[bn0_ma_reg6];
			 data_tmp_2 = A_B0R13[bn0_ma_reg6];
			 data_tmp_3 = A_B0R14[bn0_ma_reg6];
			 data_tmp_4 = A_B0R15[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg6] = A_B1R11[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg6] = A_B0R11[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg6] = A_B0R11[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg6] = A_B1R11[bn1_ma_reg8];
			 A_B1R11[bn1_ma_reg7] = data_tmp_1;
			 A_B0R11[bn0_ma_reg7] = data_tmp_2;
			 A_B0R11[bn0_ma_reg8] = data_tmp_3;
			 A_B1R11[bn1_ma_reg8] = data_tmp_4;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B1R13[bn1_ma_reg7];
			 data_tmp_2 = A_B1R14[bn1_ma_reg7];
			 data_tmp_3 = A_B1R15[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg7] = A_B0R12[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg7] = A_B0R12[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg7] = A_B1R12[bn1_ma_reg8];
			 A_B0R12[bn0_ma_reg7] = data_tmp_1;
			 A_B0R12[bn0_ma_reg8] = data_tmp_2;
			 A_B1R12[bn1_ma_reg8] = data_tmp_3;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B0R14[bn0_ma_reg7];
			 data_tmp_2 = A_B0R15[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg7] = A_B0R13[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg7] = A_B1R13[bn1_ma_reg8];
			 A_B0R13[bn0_ma_reg8] = data_tmp_1;
			 A_B1R13[bn1_ma_reg8] = data_tmp_2;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B0R15[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg8] = A_B1R14[bn1_ma_reg8];
			 A_B1R14[bn1_ma_reg8] = data_tmp_1;
			 //-----------------------------------------------------------------
		  }
		 }else{
		  if(bn1_bc_tmp > bn0_bc_tmp){
			 //Ex 0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30
			 //Ex 1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31
			 data_tmp_1  = A_B0R1[bn0_ma_reg1];
			 A_B0R1[bn0_ma_reg1] = A_B1R0[bn1_ma_reg1];
             A_B1R0[bn1_ma_reg1] = data_tmp_1;
			 data_tmp_1  = A_B0R3[bn0_ma_reg1];
			 A_B0R3[bn0_ma_reg1] = A_B1R2[bn1_ma_reg1];
             A_B1R2[bn1_ma_reg1] = data_tmp_1;
			 data_tmp_1  = A_B0R5[bn0_ma_reg1];
			 A_B0R5[bn0_ma_reg1] = A_B1R4[bn1_ma_reg1];
             A_B1R4[bn1_ma_reg1] = data_tmp_1;        
			 data_tmp_1  = A_B0R7[bn0_ma_reg1];
			 A_B0R7[bn0_ma_reg1] = A_B1R6[bn1_ma_reg1];
             A_B1R6[bn1_ma_reg1] = data_tmp_1;
			 data_tmp_1  = A_B0R9[bn0_ma_reg1];
			 A_B0R9[bn0_ma_reg1] = A_B1R8[bn1_ma_reg1];
             A_B1R8[bn1_ma_reg1] = data_tmp_1;			 
			 data_tmp_1  = A_B0R11[bn0_ma_reg1];
			 A_B0R11[bn0_ma_reg1] = A_B1R10[bn1_ma_reg1];
             A_B1R10[bn1_ma_reg1] = data_tmp_1;
			 data_tmp_1  = A_B0R13[bn0_ma_reg1];
			 A_B0R13[bn0_ma_reg1] = A_B1R12[bn1_ma_reg1];
             A_B1R12[bn1_ma_reg1] = data_tmp_1;
			 data_tmp_1  = A_B0R15[bn0_ma_reg1];
			 A_B0R15[bn0_ma_reg1] = A_B1R14[bn1_ma_reg1];
             A_B1R14[bn1_ma_reg1] = data_tmp_1;
			 /*********************************************/
			 data_tmp_1  = A_B1R1[bn1_ma_reg2];
			 A_B1R1[bn1_ma_reg2] = A_B0R0[bn0_ma_reg2];
             A_B0R0[bn0_ma_reg2] = data_tmp_1;			 
			 data_tmp_1  = A_B1R3[bn1_ma_reg2];
			 A_B1R3[bn1_ma_reg2] = A_B0R2[bn0_ma_reg2];
             A_B0R2[bn0_ma_reg2] = data_tmp_1;				 
			 data_tmp_1  = A_B1R5[bn1_ma_reg2];
			 A_B1R5[bn1_ma_reg2] = A_B0R4[bn0_ma_reg2];
             A_B0R4[bn0_ma_reg2] = data_tmp_1;	
			 data_tmp_1  = A_B1R7[bn1_ma_reg2];
			 A_B1R7[bn1_ma_reg2] = A_B0R6[bn0_ma_reg2];
             A_B0R6[bn0_ma_reg2] = data_tmp_1;	
			 data_tmp_1  = A_B1R9[bn1_ma_reg2];
			 A_B1R9[bn1_ma_reg2] = A_B0R8[bn0_ma_reg2];
             A_B0R8[bn0_ma_reg2] = data_tmp_1;	
			 data_tmp_1  = A_B1R11[bn1_ma_reg2];
			 A_B1R11[bn1_ma_reg2] = A_B0R10[bn0_ma_reg2];
             A_B0R10[bn0_ma_reg2] = data_tmp_1;	
			 data_tmp_1  = A_B1R13[bn1_ma_reg2];
			 A_B1R13[bn1_ma_reg2] = A_B0R12[bn0_ma_reg2];
             A_B0R12[bn0_ma_reg2] = data_tmp_1;	
			 data_tmp_1  = A_B1R15[bn1_ma_reg2];
			 A_B1R15[bn1_ma_reg2] = A_B0R14[bn0_ma_reg2];
             A_B0R14[bn0_ma_reg2] = data_tmp_1;	
			 /*********************************************/			 
			 data_tmp_1  = A_B1R1[bn1_ma_reg3];
			 A_B1R1[bn1_ma_reg3] = A_B0R0[bn0_ma_reg3];
             A_B0R0[bn0_ma_reg3] = data_tmp_1;			 
			 data_tmp_1  = A_B1R3[bn1_ma_reg3];
			 A_B1R3[bn1_ma_reg3] = A_B0R2[bn0_ma_reg3];
             A_B0R2[bn0_ma_reg3] = data_tmp_1;				 
			 data_tmp_1  = A_B1R5[bn1_ma_reg3];
			 A_B1R5[bn1_ma_reg3] = A_B0R4[bn0_ma_reg3];
             A_B0R4[bn0_ma_reg3] = data_tmp_1;	
			 data_tmp_1  = A_B1R7[bn1_ma_reg3];
			 A_B1R7[bn1_ma_reg3] = A_B0R6[bn0_ma_reg3];
             A_B0R6[bn0_ma_reg3] = data_tmp_1;	
			 data_tmp_1  = A_B1R9[bn1_ma_reg3];
			 A_B1R9[bn1_ma_reg3] = A_B0R8[bn0_ma_reg3];
             A_B0R8[bn0_ma_reg3] = data_tmp_1;	
			 data_tmp_1  = A_B1R11[bn1_ma_reg3];
			 A_B1R11[bn1_ma_reg3] = A_B0R10[bn0_ma_reg3];
             A_B0R10[bn0_ma_reg3] = data_tmp_1;	
			 data_tmp_1  = A_B1R13[bn1_ma_reg3];
			 A_B1R13[bn1_ma_reg3] = A_B0R12[bn0_ma_reg3];
             A_B0R12[bn0_ma_reg3] = data_tmp_1;	
			 data_tmp_1  = A_B1R15[bn1_ma_reg3];
			 A_B1R15[bn1_ma_reg3] = A_B0R14[bn0_ma_reg3];
             A_B0R14[bn0_ma_reg3] = data_tmp_1;	
			 /*********************************************/				 
			 data_tmp_1  = A_B0R1[bn0_ma_reg4];
			 A_B0R1[bn0_ma_reg4] = A_B1R0[bn1_ma_reg4];
             A_B1R0[bn1_ma_reg4] = data_tmp_1;
			 data_tmp_1  = A_B0R3[bn0_ma_reg4];
			 A_B0R3[bn0_ma_reg4] = A_B1R2[bn1_ma_reg4];
             A_B1R2[bn1_ma_reg4] = data_tmp_1;
			 data_tmp_1  = A_B0R5[bn0_ma_reg4];
			 A_B0R5[bn0_ma_reg4] = A_B1R4[bn1_ma_reg4];
             A_B1R4[bn1_ma_reg4] = data_tmp_1;        
			 data_tmp_1  = A_B0R7[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg4] = A_B1R6[bn1_ma_reg4];
             A_B1R6[bn1_ma_reg4] = data_tmp_1;
			 data_tmp_1  = A_B0R9[bn0_ma_reg4];
			 A_B0R9[bn0_ma_reg4] = A_B1R8[bn1_ma_reg4];
             A_B1R8[bn1_ma_reg4] = data_tmp_1;			 
			 data_tmp_1  = A_B0R11[bn0_ma_reg4];
			 A_B0R11[bn0_ma_reg4] = A_B1R10[bn1_ma_reg4];
             A_B1R10[bn1_ma_reg4] = data_tmp_1;
			 data_tmp_1  = A_B0R13[bn0_ma_reg4];
			 A_B0R13[bn0_ma_reg4] = A_B1R12[bn1_ma_reg4];
             A_B1R12[bn1_ma_reg4] = data_tmp_1;
			 data_tmp_1  = A_B0R15[bn0_ma_reg4];
			 A_B0R15[bn0_ma_reg4] = A_B1R14[bn1_ma_reg4];
             A_B1R14[bn1_ma_reg4] = data_tmp_1;
			 /*********************************************/
			 data_tmp_1  = A_B1R1[bn1_ma_reg5];
			 A_B1R1[bn1_ma_reg5] = A_B0R0[bn0_ma_reg5];
             A_B0R0[bn0_ma_reg5] = data_tmp_1;			 
			 data_tmp_1  = A_B1R3[bn1_ma_reg5];
			 A_B1R3[bn1_ma_reg5] = A_B0R2[bn0_ma_reg5];
             A_B0R2[bn0_ma_reg5] = data_tmp_1;				 
			 data_tmp_1  = A_B1R5[bn1_ma_reg5];
			 A_B1R5[bn1_ma_reg5] = A_B0R4[bn0_ma_reg5];
             A_B0R4[bn0_ma_reg5] = data_tmp_1;	
			 data_tmp_1  = A_B1R7[bn1_ma_reg5];
			 A_B1R7[bn1_ma_reg5] = A_B0R6[bn0_ma_reg5];
             A_B0R6[bn0_ma_reg5] = data_tmp_1;	
			 data_tmp_1  = A_B1R9[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg5] = A_B0R8[bn0_ma_reg5];
             A_B0R8[bn0_ma_reg5] = data_tmp_1;	
			 data_tmp_1  = A_B1R11[bn1_ma_reg5];
			 A_B1R11[bn1_ma_reg5] = A_B0R10[bn0_ma_reg5];
             A_B0R10[bn0_ma_reg5] = data_tmp_1;	
			 data_tmp_1  = A_B1R13[bn1_ma_reg5];
			 A_B1R13[bn1_ma_reg5] = A_B0R12[bn0_ma_reg5];
             A_B0R12[bn0_ma_reg5] = data_tmp_1;	
			 data_tmp_1  = A_B1R15[bn1_ma_reg5];
			 A_B1R15[bn1_ma_reg5] = A_B0R14[bn0_ma_reg5];
             A_B0R14[bn0_ma_reg5] = data_tmp_1;	
			 /*********************************************/
			 data_tmp_1  = A_B0R1[bn0_ma_reg6];
			 A_B0R1[bn0_ma_reg6] = A_B1R0[bn1_ma_reg6];
             A_B1R0[bn1_ma_reg6] = data_tmp_1;
			 data_tmp_1  = A_B0R3[bn0_ma_reg6];
			 A_B0R3[bn0_ma_reg6] = A_B1R2[bn1_ma_reg6];
             A_B1R2[bn1_ma_reg6] = data_tmp_1;
			 data_tmp_1  = A_B0R5[bn0_ma_reg6];
			 A_B0R5[bn0_ma_reg6] = A_B1R4[bn1_ma_reg6];
             A_B1R4[bn1_ma_reg6] = data_tmp_1;        
			 data_tmp_1  = A_B0R7[bn0_ma_reg6];
			 A_B0R7[bn0_ma_reg6] = A_B1R6[bn1_ma_reg6];
             A_B1R6[bn1_ma_reg6] = data_tmp_1;
			 data_tmp_1  = A_B0R9[bn0_ma_reg6];
			 A_B0R9[bn0_ma_reg6] = A_B1R8[bn1_ma_reg6];
             A_B1R8[bn1_ma_reg6] = data_tmp_1;			 
			 data_tmp_1  = A_B0R11[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg6] = A_B1R10[bn1_ma_reg6];
             A_B1R10[bn1_ma_reg6] = data_tmp_1;
			 data_tmp_1  = A_B0R13[bn0_ma_reg6];
			 A_B0R13[bn0_ma_reg6] = A_B1R12[bn1_ma_reg6];
             A_B1R12[bn1_ma_reg6] = data_tmp_1;
			 data_tmp_1  = A_B0R15[bn0_ma_reg6];
			 A_B0R15[bn0_ma_reg6] = A_B1R14[bn1_ma_reg6];
             A_B1R14[bn1_ma_reg6] = data_tmp_1;
			 /*********************************************/	
			 data_tmp_1  = A_B0R1[bn0_ma_reg7];
			 A_B0R1[bn0_ma_reg7] = A_B1R0[bn1_ma_reg7];
             A_B1R0[bn1_ma_reg7] = data_tmp_1;
			 data_tmp_1  = A_B0R3[bn0_ma_reg7];
			 A_B0R3[bn0_ma_reg7] = A_B1R2[bn1_ma_reg7];
             A_B1R2[bn1_ma_reg7] = data_tmp_1;
			 data_tmp_1  = A_B0R5[bn0_ma_reg7];
			 A_B0R5[bn0_ma_reg7] = A_B1R4[bn1_ma_reg7];
             A_B1R4[bn1_ma_reg7] = data_tmp_1;        
			 data_tmp_1  = A_B0R7[bn0_ma_reg7];
			 A_B0R7[bn0_ma_reg7] = A_B1R6[bn1_ma_reg7];
             A_B1R6[bn1_ma_reg7] = data_tmp_1;
			 data_tmp_1  = A_B0R9[bn0_ma_reg7];
			 A_B0R9[bn0_ma_reg7] = A_B1R8[bn1_ma_reg7];
             A_B1R8[bn1_ma_reg7] = data_tmp_1;			 
			 data_tmp_1  = A_B0R11[bn0_ma_reg7];
			 A_B0R11[bn0_ma_reg7] = A_B1R10[bn1_ma_reg7];
             A_B1R10[bn1_ma_reg7] = data_tmp_1;
			 data_tmp_1  = A_B0R13[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg7] = A_B1R12[bn1_ma_reg7];
             A_B1R12[bn1_ma_reg7] = data_tmp_1;
			 data_tmp_1  = A_B0R15[bn0_ma_reg7];
			 A_B0R15[bn0_ma_reg7] = A_B1R14[bn1_ma_reg7];
             A_B1R14[bn1_ma_reg7] = data_tmp_1;
			 /*********************************************/
			 data_tmp_1  = A_B1R1[bn1_ma_reg8];
			 A_B1R1[bn1_ma_reg8] = A_B0R0[bn0_ma_reg8];
             A_B0R0[bn0_ma_reg8] = data_tmp_1;			 
			 data_tmp_1  = A_B1R3[bn1_ma_reg8];
			 A_B1R3[bn1_ma_reg8] = A_B0R2[bn0_ma_reg8];
             A_B0R2[bn0_ma_reg8] = data_tmp_1;				 
			 data_tmp_1  = A_B1R5[bn1_ma_reg8];
			 A_B1R5[bn1_ma_reg8] = A_B0R4[bn0_ma_reg8];
             A_B0R4[bn0_ma_reg8] = data_tmp_1;	
			 data_tmp_1  = A_B1R7[bn1_ma_reg8];
			 A_B1R7[bn1_ma_reg8] = A_B0R6[bn0_ma_reg8];
             A_B0R6[bn0_ma_reg8] = data_tmp_1;	
			 data_tmp_1  = A_B1R9[bn1_ma_reg8];
			 A_B1R9[bn1_ma_reg8] = A_B0R8[bn0_ma_reg8];
             A_B0R8[bn0_ma_reg8] = data_tmp_1;	
			 data_tmp_1  = A_B1R11[bn1_ma_reg8];
			 A_B1R11[bn1_ma_reg8] = A_B0R10[bn0_ma_reg8];
             A_B0R10[bn0_ma_reg8] = data_tmp_1;	
			 data_tmp_1  = A_B1R13[bn1_ma_reg8];
			 A_B1R13[bn1_ma_reg8] = A_B0R12[bn0_ma_reg8];
             A_B0R12[bn0_ma_reg8] = data_tmp_1;	
			 data_tmp_1  = A_B1R15[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg8] = A_B0R14[bn0_ma_reg8];
             A_B0R14[bn0_ma_reg8] = data_tmp_1;	
			 /*********************************************/			 
		  }else{	 
			 data_tmp_1  = A_B1R1[bn1_ma_reg1];
			 A_B1R1[bn1_ma_reg1] = A_B0R0[bn0_ma_reg1];
             A_B0R0[bn0_ma_reg1] = data_tmp_1;			 
			 data_tmp_1  = A_B1R3[bn1_ma_reg1];
			 A_B1R3[bn1_ma_reg1] = A_B0R2[bn0_ma_reg1];
             A_B0R2[bn0_ma_reg1] = data_tmp_1;				 
			 data_tmp_1  = A_B1R5[bn1_ma_reg1];
			 A_B1R5[bn1_ma_reg1] = A_B0R4[bn0_ma_reg1];
             A_B0R4[bn0_ma_reg1] = data_tmp_1;	
			 data_tmp_1  = A_B1R7[bn1_ma_reg1];
			 A_B1R7[bn1_ma_reg1] = A_B0R6[bn0_ma_reg1];
             A_B0R6[bn0_ma_reg1] = data_tmp_1;	
			 data_tmp_1  = A_B1R9[bn1_ma_reg1];
			 A_B1R9[bn1_ma_reg1] = A_B0R8[bn0_ma_reg1];
             A_B0R8[bn0_ma_reg1] = data_tmp_1;	
			 data_tmp_1  = A_B1R11[bn1_ma_reg1];
			 A_B1R11[bn1_ma_reg1] = A_B0R10[bn0_ma_reg1];
             A_B0R10[bn0_ma_reg1] = data_tmp_1;	
			 data_tmp_1  = A_B1R13[bn1_ma_reg1];
			 A_B1R13[bn1_ma_reg1] = A_B0R12[bn0_ma_reg1];
             A_B0R12[bn0_ma_reg1] = data_tmp_1;	
			 data_tmp_1  = A_B1R15[bn1_ma_reg1];
			 A_B1R15[bn1_ma_reg1] = A_B0R14[bn0_ma_reg1];
             A_B0R14[bn0_ma_reg1] = data_tmp_1;	
			 /*********************************************/
			 data_tmp_1  = A_B0R1[bn0_ma_reg2];
			 A_B0R1[bn0_ma_reg2] = A_B1R0[bn1_ma_reg2];
             A_B1R0[bn1_ma_reg2] = data_tmp_1;
			 data_tmp_1  = A_B0R3[bn0_ma_reg2];
			 A_B0R3[bn0_ma_reg2] = A_B1R2[bn1_ma_reg2];
             A_B1R2[bn1_ma_reg2] = data_tmp_1;
			 data_tmp_1  = A_B0R5[bn0_ma_reg2];
			 A_B0R5[bn0_ma_reg2] = A_B1R4[bn1_ma_reg2];
             A_B1R4[bn1_ma_reg2] = data_tmp_1;        
			 data_tmp_1  = A_B0R7[bn0_ma_reg2];
			 A_B0R7[bn0_ma_reg2] = A_B1R6[bn1_ma_reg2];
             A_B1R6[bn1_ma_reg2] = data_tmp_1;
			 data_tmp_1  = A_B0R9[bn0_ma_reg2];
			 A_B0R9[bn0_ma_reg2] = A_B1R8[bn1_ma_reg2];
             A_B1R8[bn1_ma_reg2] = data_tmp_1;			 
			 data_tmp_1  = A_B0R11[bn0_ma_reg2];
			 A_B0R11[bn0_ma_reg2] = A_B1R10[bn1_ma_reg2];
             A_B1R10[bn1_ma_reg2] = data_tmp_1;
			 data_tmp_1  = A_B0R13[bn0_ma_reg2];
			 A_B0R13[bn0_ma_reg2] = A_B1R12[bn1_ma_reg2];
             A_B1R12[bn1_ma_reg2] = data_tmp_1;
			 data_tmp_1  = A_B0R15[bn0_ma_reg2];
			 A_B0R15[bn0_ma_reg2] = A_B1R14[bn1_ma_reg2];
             A_B1R14[bn1_ma_reg2] = data_tmp_1;
			 /*********************************************/	
			 data_tmp_1  = A_B0R1[bn0_ma_reg3];
			 A_B0R1[bn0_ma_reg3] = A_B1R0[bn1_ma_reg3];
             A_B1R0[bn1_ma_reg3] = data_tmp_1;
			 data_tmp_1  = A_B0R3[bn0_ma_reg3];
			 A_B0R3[bn0_ma_reg3] = A_B1R2[bn1_ma_reg3];
             A_B1R2[bn1_ma_reg3] = data_tmp_1;
			 data_tmp_1  = A_B0R5[bn0_ma_reg3];
			 A_B0R5[bn0_ma_reg3] = A_B1R4[bn1_ma_reg3];
             A_B1R4[bn1_ma_reg3] = data_tmp_1;        
			 data_tmp_1  = A_B0R7[bn0_ma_reg3];
			 A_B0R7[bn0_ma_reg3] = A_B1R6[bn1_ma_reg3];
             A_B1R6[bn1_ma_reg3] = data_tmp_1;
			 data_tmp_1  = A_B0R9[bn0_ma_reg3];
			 A_B0R9[bn0_ma_reg3] = A_B1R8[bn1_ma_reg3];
             A_B1R8[bn1_ma_reg3] = data_tmp_1;			 
			 data_tmp_1  = A_B0R11[bn0_ma_reg3];
			 A_B0R11[bn0_ma_reg3] = A_B1R10[bn1_ma_reg3];
             A_B1R10[bn1_ma_reg3] = data_tmp_1;
			 data_tmp_1  = A_B0R13[bn0_ma_reg3];
			 A_B0R13[bn0_ma_reg3] = A_B1R12[bn1_ma_reg3];
             A_B1R12[bn1_ma_reg3] = data_tmp_1;
			 data_tmp_1  = A_B0R15[bn0_ma_reg3];
			 A_B0R15[bn0_ma_reg3] = A_B1R14[bn1_ma_reg3];
             A_B1R14[bn1_ma_reg3] = data_tmp_1;
			 /*********************************************/
			 data_tmp_1  = A_B1R1[bn1_ma_reg4];
			 A_B1R1[bn1_ma_reg4] = A_B0R0[bn0_ma_reg4];
             A_B0R0[bn0_ma_reg4] = data_tmp_1;			 
			 data_tmp_1  = A_B1R3[bn1_ma_reg4];
			 A_B1R3[bn1_ma_reg4] = A_B0R2[bn0_ma_reg4];
             A_B0R2[bn0_ma_reg4] = data_tmp_1;				 
			 data_tmp_1  = A_B1R5[bn1_ma_reg4];
			 A_B1R5[bn1_ma_reg4] = A_B0R4[bn0_ma_reg4];
             A_B0R4[bn0_ma_reg4] = data_tmp_1;	
			 data_tmp_1  = A_B1R7[bn1_ma_reg4];
			 A_B1R7[bn1_ma_reg4] = A_B0R6[bn0_ma_reg4];
             A_B0R6[bn0_ma_reg4] = data_tmp_1;	
			 data_tmp_1  = A_B1R9[bn1_ma_reg4];
			 A_B1R9[bn1_ma_reg4] = A_B0R8[bn0_ma_reg4];
             A_B0R8[bn0_ma_reg4] = data_tmp_1;	
			 data_tmp_1  = A_B1R11[bn1_ma_reg4];
			 A_B1R11[bn1_ma_reg4] = A_B0R10[bn0_ma_reg4];
             A_B0R10[bn0_ma_reg4] = data_tmp_1;	
			 data_tmp_1  = A_B1R13[bn1_ma_reg4];
			 A_B1R13[bn1_ma_reg4] = A_B0R12[bn0_ma_reg4];
             A_B0R12[bn0_ma_reg4] = data_tmp_1;	
			 data_tmp_1  = A_B1R15[bn1_ma_reg4];
			 A_B1R15[bn1_ma_reg4] = A_B0R14[bn0_ma_reg4];
             A_B0R14[bn0_ma_reg4] = data_tmp_1;	
			 /*********************************************/	
			 data_tmp_1  = A_B0R1[bn0_ma_reg5];
			 A_B0R1[bn0_ma_reg5] = A_B1R0[bn1_ma_reg5];
             A_B1R0[bn1_ma_reg5] = data_tmp_1;
			 data_tmp_1  = A_B0R3[bn0_ma_reg5];
			 A_B0R3[bn0_ma_reg5] = A_B1R2[bn1_ma_reg5];
             A_B1R2[bn1_ma_reg5] = data_tmp_1;
			 data_tmp_1  = A_B0R5[bn0_ma_reg5];
			 A_B0R5[bn0_ma_reg5] = A_B1R4[bn1_ma_reg5];
             A_B1R4[bn1_ma_reg5] = data_tmp_1;        
			 data_tmp_1  = A_B0R7[bn0_ma_reg5];
			 A_B0R7[bn0_ma_reg5] = A_B1R6[bn1_ma_reg5];
             A_B1R6[bn1_ma_reg5] = data_tmp_1;
			 data_tmp_1  = A_B0R9[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg5] = A_B1R8[bn1_ma_reg5];
             A_B1R8[bn1_ma_reg5] = data_tmp_1;			 
			 data_tmp_1  = A_B0R11[bn0_ma_reg5];
			 A_B0R11[bn0_ma_reg5] = A_B1R10[bn1_ma_reg5];
             A_B1R10[bn1_ma_reg5] = data_tmp_1;
			 data_tmp_1  = A_B0R13[bn0_ma_reg5];
			 A_B0R13[bn0_ma_reg5] = A_B1R12[bn1_ma_reg5];
             A_B1R12[bn1_ma_reg5] = data_tmp_1;
			 data_tmp_1  = A_B0R15[bn0_ma_reg5];
			 A_B0R15[bn0_ma_reg5] = A_B1R14[bn1_ma_reg5];
             A_B1R14[bn1_ma_reg5] = data_tmp_1;
			 /*********************************************/
			 data_tmp_1  = A_B1R1[bn1_ma_reg6];
			 A_B1R1[bn1_ma_reg6] = A_B0R0[bn0_ma_reg6];
             A_B0R0[bn0_ma_reg6] = data_tmp_1;			 
			 data_tmp_1  = A_B1R3[bn1_ma_reg6];
			 A_B1R3[bn1_ma_reg6] = A_B0R2[bn0_ma_reg6];
             A_B0R2[bn0_ma_reg6] = data_tmp_1;				 
			 data_tmp_1  = A_B1R5[bn1_ma_reg6];
			 A_B1R5[bn1_ma_reg6] = A_B0R4[bn0_ma_reg6];
             A_B0R4[bn0_ma_reg6] = data_tmp_1;	
			 data_tmp_1  = A_B1R7[bn1_ma_reg6];
			 A_B1R7[bn1_ma_reg6] = A_B0R6[bn0_ma_reg6];
             A_B0R6[bn0_ma_reg6] = data_tmp_1;	
			 data_tmp_1  = A_B1R9[bn1_ma_reg6];
			 A_B1R9[bn1_ma_reg6] = A_B0R8[bn0_ma_reg6];
             A_B0R8[bn0_ma_reg6] = data_tmp_1;	
			 data_tmp_1  = A_B1R11[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg6] = A_B0R10[bn0_ma_reg6];
             A_B0R10[bn0_ma_reg6] = data_tmp_1;	
			 data_tmp_1  = A_B1R13[bn1_ma_reg6];
			 A_B1R13[bn1_ma_reg6] = A_B0R12[bn0_ma_reg6];
             A_B0R12[bn0_ma_reg6] = data_tmp_1;	
			 data_tmp_1  = A_B1R15[bn1_ma_reg6];
			 A_B1R15[bn1_ma_reg6] = A_B0R14[bn0_ma_reg6];
             A_B0R14[bn0_ma_reg6] = data_tmp_1;	
			 /*********************************************/	
			 data_tmp_1  = A_B1R1[bn1_ma_reg7];
			 A_B1R1[bn1_ma_reg7] = A_B0R0[bn0_ma_reg7];
             A_B0R0[bn0_ma_reg7] = data_tmp_1;			 
			 data_tmp_1  = A_B1R3[bn1_ma_reg7];
			 A_B1R3[bn1_ma_reg7] = A_B0R2[bn0_ma_reg7];
             A_B0R2[bn0_ma_reg7] = data_tmp_1;				 
			 data_tmp_1  = A_B1R5[bn1_ma_reg7];
			 A_B1R5[bn1_ma_reg7] = A_B0R4[bn0_ma_reg7];
             A_B0R4[bn0_ma_reg7] = data_tmp_1;	
			 data_tmp_1  = A_B1R7[bn1_ma_reg7];
			 A_B1R7[bn1_ma_reg7] = A_B0R6[bn0_ma_reg7];
             A_B0R6[bn0_ma_reg7] = data_tmp_1;	
			 data_tmp_1  = A_B1R9[bn1_ma_reg7];
			 A_B1R9[bn1_ma_reg7] = A_B0R8[bn0_ma_reg7];
             A_B0R8[bn0_ma_reg7] = data_tmp_1;	
			 data_tmp_1  = A_B1R11[bn1_ma_reg7];
			 A_B1R11[bn1_ma_reg7] = A_B0R10[bn0_ma_reg7];
             A_B0R10[bn0_ma_reg7] = data_tmp_1;	
			 data_tmp_1  = A_B1R13[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg7] = A_B0R12[bn0_ma_reg7];
             A_B0R12[bn0_ma_reg7] = data_tmp_1;	
			 data_tmp_1  = A_B1R15[bn1_ma_reg7];
			 A_B1R15[bn1_ma_reg7] = A_B0R14[bn0_ma_reg7];
             A_B0R14[bn0_ma_reg7] = data_tmp_1;	
			 /*********************************************/	
			 data_tmp_1  = A_B0R1[bn0_ma_reg8];
			 A_B0R1[bn0_ma_reg8] = A_B1R0[bn1_ma_reg8];
             A_B1R0[bn1_ma_reg8] = data_tmp_1;
			 data_tmp_1  = A_B0R3[bn0_ma_reg8];
			 A_B0R3[bn0_ma_reg8] = A_B1R2[bn1_ma_reg8];
             A_B1R2[bn1_ma_reg8] = data_tmp_1;
			 data_tmp_1  = A_B0R5[bn0_ma_reg8];
			 A_B0R5[bn0_ma_reg8] = A_B1R4[bn1_ma_reg8];
             A_B1R4[bn1_ma_reg8] = data_tmp_1;        
			 data_tmp_1  = A_B0R7[bn0_ma_reg8];
			 A_B0R7[bn0_ma_reg8] = A_B1R6[bn1_ma_reg8];
             A_B1R6[bn1_ma_reg8] = data_tmp_1;
			 data_tmp_1  = A_B0R9[bn0_ma_reg8];
			 A_B0R9[bn0_ma_reg8] = A_B1R8[bn1_ma_reg8];
             A_B1R8[bn1_ma_reg8] = data_tmp_1;			 
			 data_tmp_1  = A_B0R11[bn0_ma_reg8];
			 A_B0R11[bn0_ma_reg8] = A_B1R10[bn1_ma_reg8];
             A_B1R10[bn1_ma_reg8] = data_tmp_1;
			 data_tmp_1  = A_B0R13[bn0_ma_reg8];
			 A_B0R13[bn0_ma_reg8] = A_B1R12[bn1_ma_reg8];
             A_B1R12[bn1_ma_reg8] = data_tmp_1;
			 data_tmp_1  = A_B0R15[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg8] = A_B1R14[bn1_ma_reg8];
             A_B1R14[bn1_ma_reg8] = data_tmp_1;			 
             //-----------------------------------------------------------

		    }	
		 } 
		}
	}
	
	std::cout << "radix-16 FFT computing over!!\n";
	
	//DATARECORD <<"----------------------------------------------------------------------------- \n";
	//DATARECORD <<" Radix-2 FFT computing start!!! \n";
	
	
	for(int i = 0;i < group; i++){
        for(int j=0; j < radix; j++){
        	gray_i  = Gray(i,group);
        	BC_tmp  = j * group + gray_i;
			//DATARECORD << "BC_tmp: " << BC_tmp << "\n";
        	BC = BC_tmp;
			//DATARECORD << "After RR_R16 , BC : " << BC << "\n";
        	AGU_R16(BC,bn_tmp,ma_tmp);
		    //DATARECORD << "BN : " << bn_tmp << "\n";
		    //DATARECORD << "MA : " << ma_tmp << "\n";			
        	if(bn_tmp == 0){
        		if(j < 2)bn0_bc_tmp = BC_tmp;
				//DATARECORD << "Before Radix-2 butterfly unit operation!!! \n";
				//DATARECORD << "A_B0R0["<< ma_tmp <<"]: " << A_B0R0[ma_tmp]<<"\n";
				//DATARECORD << "A_B0R1["<< ma_tmp <<"]: " << A_B0R1[ma_tmp]<<"\n";
				//DATARECORD << "A_B0R2["<< ma_tmp <<"]: " << A_B0R2[ma_tmp]<<"\n";
				//DATARECORD << "A_B0R3["<< ma_tmp <<"]: " << A_B0R3[ma_tmp]<<"\n";
				//DATARECORD << "A_B0R4["<< ma_tmp <<"]: " << A_B0R4[ma_tmp]<<"\n";
				//DATARECORD << "A_B0R5["<< ma_tmp <<"]: " << A_B0R5[ma_tmp]<<"\n";
				//DATARECORD << "A_B0R6["<< ma_tmp <<"]: " << A_B0R6[ma_tmp]<<"\n";
				//DATARECORD << "A_B0R7["<< ma_tmp <<"]: " << A_B0R7[ma_tmp]<<"\n";
				//DATARECORD << "A_B0R8["<< ma_tmp <<"]: " << A_B0R8[ma_tmp]<<"\n";
				//DATARECORD << "A_B0R9["<< ma_tmp <<"]: " << A_B0R9[ma_tmp]<<"\n";
				//DATARECORD << "A_B0R10["<< ma_tmp <<"]: "<< A_B0R10[ma_tmp]<<"\n";
				//DATARECORD << "A_B0R11["<< ma_tmp <<"]: "<< A_B0R11[ma_tmp]<<"\n";
				//DATARECORD << "A_B0R12["<< ma_tmp <<"]: "<< A_B0R12[ma_tmp]<<"\n";
				//DATARECORD << "A_B0R13["<< ma_tmp <<"]: "<< A_B0R13[ma_tmp]<<"\n";
				//DATARECORD << "A_B0R14["<< ma_tmp <<"]: "<< A_B0R14[ma_tmp]<<"\n";
				//DATARECORD << "A_B0R15["<< ma_tmp <<"]: "<< A_B0R15[ma_tmp]<<"\n";					
        		Radix2_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp]);
        		Radix2_BU(A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
        		Radix2_BU(A_B0R4[ma_tmp],A_B0R5[ma_tmp]);
        		Radix2_BU(A_B0R6[ma_tmp],A_B0R7[ma_tmp]);
        		Radix2_BU(A_B0R8[ma_tmp],A_B0R9[ma_tmp]);
        		Radix2_BU(A_B0R10[ma_tmp],A_B0R11[ma_tmp]);
        		Radix2_BU(A_B0R12[ma_tmp],A_B0R13[ma_tmp]);
        		Radix2_BU(A_B0R14[ma_tmp],A_B0R15[ma_tmp]);
				//DATARECORD << "After Radix-2 butterfly unit operation!!! \n";
				//DATARECORD << "A_B0R0["<< ma_tmp <<"]: " << A_B0R0[ma_tmp]<<"\n";
				//DATARECORD << "A_B0R1["<< ma_tmp <<"]: " << A_B0R1[ma_tmp]<<"\n";
				//DATARECORD << "A_B0R2["<< ma_tmp <<"]: " << A_B0R2[ma_tmp]<<"\n";
				//DATARECORD << "A_B0R3["<< ma_tmp <<"]: " << A_B0R3[ma_tmp]<<"\n";
				//DATARECORD << "A_B0R4["<< ma_tmp <<"]: " << A_B0R4[ma_tmp]<<"\n";
				//DATARECORD << "A_B0R5["<< ma_tmp <<"]: " << A_B0R5[ma_tmp]<<"\n";
				//DATARECORD << "A_B0R6["<< ma_tmp <<"]: " << A_B0R6[ma_tmp]<<"\n";
				//DATARECORD << "A_B0R7["<< ma_tmp <<"]: " << A_B0R7[ma_tmp]<<"\n";
				//DATARECORD << "A_B0R8["<< ma_tmp <<"]: " << A_B0R8[ma_tmp]<<"\n";
				//DATARECORD << "A_B0R9["<< ma_tmp <<"]: " << A_B0R9[ma_tmp]<<"\n";
				//DATARECORD << "A_B0R10["<< ma_tmp <<"]: "<< A_B0R10[ma_tmp]<<"\n";
				//DATARECORD << "A_B0R11["<< ma_tmp <<"]: "<< A_B0R11[ma_tmp]<<"\n";
				//DATARECORD << "A_B0R12["<< ma_tmp <<"]: "<< A_B0R12[ma_tmp]<<"\n";
				//DATARECORD << "A_B0R13["<< ma_tmp <<"]: "<< A_B0R13[ma_tmp]<<"\n";
				//DATARECORD << "A_B0R14["<< ma_tmp <<"]: "<< A_B0R14[ma_tmp]<<"\n";
				//DATARECORD << "A_B0R15["<< ma_tmp <<"]: "<< A_B0R15[ma_tmp]<<"\n";						
        	}else {
				//DATARECORD << "Before Radix-2 butterfly unit operation!!! \n";
				//DATARECORD << "A_B1R0["<< ma_tmp <<"]: " << A_B1R0[ma_tmp]<<"\n";
				//DATARECORD << "A_B1R1["<< ma_tmp <<"]: " << A_B1R1[ma_tmp]<<"\n";
				//DATARECORD << "A_B1R2["<< ma_tmp <<"]: " << A_B1R2[ma_tmp]<<"\n";
				//DATARECORD << "A_B1R3["<< ma_tmp <<"]: " << A_B1R3[ma_tmp]<<"\n";
				//DATARECORD << "A_B1R4["<< ma_tmp <<"]: " << A_B1R4[ma_tmp]<<"\n";
				//DATARECORD << "A_B1R5["<< ma_tmp <<"]: " << A_B1R5[ma_tmp]<<"\n";
				//DATARECORD << "A_B1R6["<< ma_tmp <<"]: " << A_B1R6[ma_tmp]<<"\n";
				//DATARECORD << "A_B1R7["<< ma_tmp <<"]: " << A_B1R7[ma_tmp]<<"\n";
				//DATARECORD << "A_B1R8["<< ma_tmp <<"]: " << A_B1R8[ma_tmp]<<"\n";
				//DATARECORD << "A_B1R9["<< ma_tmp <<"]: " << A_B1R9[ma_tmp]<<"\n";
				//DATARECORD << "A_B1R10["<< ma_tmp <<"]: "<< A_B1R10[ma_tmp]<<"\n";
				//DATARECORD << "A_B1R11["<< ma_tmp <<"]: "<< A_B1R11[ma_tmp]<<"\n";
				//DATARECORD << "A_B1R12["<< ma_tmp <<"]: "<< A_B1R12[ma_tmp]<<"\n";
				//DATARECORD << "A_B1R13["<< ma_tmp <<"]: "<< A_B1R13[ma_tmp]<<"\n";
				//DATARECORD << "A_B1R14["<< ma_tmp <<"]: "<< A_B1R14[ma_tmp]<<"\n";
				//DATARECORD << "A_B1R15["<< ma_tmp <<"]: "<< A_B1R15[ma_tmp]<<"\n";				
        	    Radix2_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp]);
        	    Radix2_BU(A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
        	    Radix2_BU(A_B1R4[ma_tmp],A_B1R5[ma_tmp]);
        	    Radix2_BU(A_B1R6[ma_tmp],A_B1R7[ma_tmp]);
        	    Radix2_BU(A_B1R8[ma_tmp],A_B1R9[ma_tmp]);
        	    Radix2_BU(A_B1R10[ma_tmp],A_B1R11[ma_tmp]);
        	    Radix2_BU(A_B1R12[ma_tmp],A_B1R13[ma_tmp]);
        	    Radix2_BU(A_B1R14[ma_tmp],A_B1R15[ma_tmp]);
				//DATARECORD << "After Radix-2 butterfly unit operation!!! \n";
				//DATARECORD << "A_B1R0["<< ma_tmp <<"]: " << A_B1R0[ma_tmp]<<"\n";
				//DATARECORD << "A_B1R1["<< ma_tmp <<"]: " << A_B1R1[ma_tmp]<<"\n";
				//DATARECORD << "A_B1R2["<< ma_tmp <<"]: " << A_B1R2[ma_tmp]<<"\n";
				//DATARECORD << "A_B1R3["<< ma_tmp <<"]: " << A_B1R3[ma_tmp]<<"\n";
				//DATARECORD << "A_B1R4["<< ma_tmp <<"]: " << A_B1R4[ma_tmp]<<"\n";
				//DATARECORD << "A_B1R5["<< ma_tmp <<"]: " << A_B1R5[ma_tmp]<<"\n";
				//DATARECORD << "A_B1R6["<< ma_tmp <<"]: " << A_B1R6[ma_tmp]<<"\n";
				//DATARECORD << "A_B1R7["<< ma_tmp <<"]: " << A_B1R7[ma_tmp]<<"\n";
				//DATARECORD << "A_B1R8["<< ma_tmp <<"]: " << A_B1R8[ma_tmp]<<"\n";
				//DATARECORD << "A_B1R9["<< ma_tmp <<"]: " << A_B1R9[ma_tmp]<<"\n";
				//DATARECORD << "A_B1R10["<< ma_tmp <<"]: "<< A_B1R10[ma_tmp]<<"\n";
				//DATARECORD << "A_B1R11["<< ma_tmp <<"]: "<< A_B1R11[ma_tmp]<<"\n";
				//DATARECORD << "A_B1R12["<< ma_tmp <<"]: "<< A_B1R12[ma_tmp]<<"\n";
				//DATARECORD << "A_B1R13["<< ma_tmp <<"]: "<< A_B1R13[ma_tmp]<<"\n";
				//DATARECORD << "A_B1R14["<< ma_tmp <<"]: "<< A_B1R14[ma_tmp]<<"\n";
				//DATARECORD << "A_B1R15["<< ma_tmp <<"]: "<< A_B1R15[ma_tmp]<<"\n";					
        	}			
        }		
	}
	
	int index0_i;
	int index1_i;
	int index2_i;
	int index3_i;
	int index4_i;
	int index5_i;
	int index6_i;
	int index7_i;
	int index8_i;
	int index9_i;
	int index10_i;
	int index11_i;
	int index12_i;
	int index13_i;
	int index14_i;
	int index15_i;
	
	int index0;
    int index1;
    int index2;
    int index3;
    int index4;
    int index5;
    int index6;
    int index7;
    int index8;
    int index9;
    int index10;
    int index11;
    int index12;
    int index13;
    int index14;
    int index15;
	
	int BC_INDEX_32FLAG;
	int BC_INDEX_2FLAG;
	
	//DATARECORD <<"***************************************************************\n";
	//DATARECORD <<"***** DATA OUTPUT!!                                         ***\n";
	//DATARECORD <<"***************************************************************\n";
	//data output
	// SPMB data output , output function is "NTT_REORDERINDEX_R16_R2_OUT".
	for(int i = 0;i < group;i++){
		for(int j = 0;j < radix;j++){
			gray_i  = Gray(i,group);
			BC_tmp  = j * group + gray_i;
			//DATARECORD <<"---------------------------------------------------\n";
			//DATARECORD <<" BC_tmp: "<< BC_tmp <<"\n";
			RR_R16_R2(BC_tmp,((4 * (Stage-1)) - 3),BC);
            BC_INDEX_32FLAG = BC >> 1;
			BC_INDEX_2FLAG  = BC % 2;
			//DATARECORD <<" BC: "<< BC <<"\n";
			//DATARECORD <<" BC_INDEX_32FLAG: "<< BC_INDEX_32FLAG <<"\n";
			//DATARECORD <<" BC_INDEX_2FLAG:  "<< BC_INDEX_2FLAG  <<"\n";
			AGU_R16(BC,bn_tmp,ma_tmp);
			//DATARECORD <<" bn_tmp:  "<< bn_tmp <<"\n";
			//DATARECORD <<" ma_tmp:  "<< ma_tmp <<"\n";
			index0_i  = BC_INDEX_32FLAG * 32  + BC_INDEX_2FLAG * 2 + 0;
			index1_i  = BC_INDEX_32FLAG * 32  + BC_INDEX_2FLAG * 2 + 1;
			index2_i  = BC_INDEX_32FLAG * 32  + BC_INDEX_2FLAG * 2 + 4;
			index3_i  = BC_INDEX_32FLAG * 32  + BC_INDEX_2FLAG * 2 + 5;
			index4_i  = BC_INDEX_32FLAG * 32  + BC_INDEX_2FLAG * 2 + 8;
			index5_i  = BC_INDEX_32FLAG * 32  + BC_INDEX_2FLAG * 2 + 9;
			index6_i  = BC_INDEX_32FLAG * 32  + BC_INDEX_2FLAG * 2 + 12;
			index7_i  = BC_INDEX_32FLAG * 32  + BC_INDEX_2FLAG * 2 + 13;
			index8_i  = BC_INDEX_32FLAG * 32  + BC_INDEX_2FLAG * 2 + 16;
			index9_i  = BC_INDEX_32FLAG * 32  + BC_INDEX_2FLAG * 2 + 17;
			index10_i = BC_INDEX_32FLAG * 32  + BC_INDEX_2FLAG * 2 + 20;
			index11_i = BC_INDEX_32FLAG * 32  + BC_INDEX_2FLAG * 2 + 21;
			index12_i = BC_INDEX_32FLAG * 32  + BC_INDEX_2FLAG * 2 + 24;
			index13_i = BC_INDEX_32FLAG * 32  + BC_INDEX_2FLAG * 2 + 25;
			index14_i = BC_INDEX_32FLAG * 32  + BC_INDEX_2FLAG * 2 + 28;
			index15_i = BC_INDEX_32FLAG * 32  + BC_INDEX_2FLAG * 2 + 29;
			
			//DATARECORD <<" index0_i:  "<< index0_i <<"\n";
			//DATARECORD <<" index1_i:  "<< index1_i <<"\n";
			//DATARECORD <<" index2_i:  "<< index2_i <<"\n";
			//DATARECORD <<" index3_i:  "<< index3_i <<"\n";
			//DATARECORD <<" index4_i:  "<< index4_i <<"\n";
			//DATARECORD <<" index5_i:  "<< index5_i <<"\n";
			//DATARECORD <<" index6_i:  "<< index6_i <<"\n";
			//DATARECORD <<" index7_i:  "<< index7_i <<"\n";
			//DATARECORD <<" index8_i:  "<< index8_i <<"\n";
			//DATARECORD <<" index9_i:  "<< index9_i <<"\n";
			//DATARECORD <<" index10_i: "<< index10_i <<"\n";
			//DATARECORD <<" index11_i: "<< index11_i <<"\n";
			//DATARECORD <<" index12_i: "<< index12_i <<"\n";
			//DATARECORD <<" index13_i: "<< index13_i <<"\n";
			//DATARECORD <<" index14_i: "<< index14_i <<"\n";
			//DATARECORD <<" index15_i: "<< index15_i <<"\n";
					
			if(bn_tmp == 0){
			   NTT_REORDERINDEX_R16_R2_OUT(index0_i,index0);
			   NTT_REORDERINDEX_R16_R2_OUT(index1_i,index1);
			   NTT_REORDERINDEX_R16_R2_OUT(index2_i,index2);
			   NTT_REORDERINDEX_R16_R2_OUT(index3_i,index3);
			   NTT_REORDERINDEX_R16_R2_OUT(index4_i,index4);
			   NTT_REORDERINDEX_R16_R2_OUT(index5_i,index5);
			   NTT_REORDERINDEX_R16_R2_OUT(index6_i,index6);
			   NTT_REORDERINDEX_R16_R2_OUT(index7_i,index7);
			   NTT_REORDERINDEX_R16_R2_OUT(index8_i,index8);
			   NTT_REORDERINDEX_R16_R2_OUT(index9_i,index9);
			   NTT_REORDERINDEX_R16_R2_OUT(index10_i,index10);
			   NTT_REORDERINDEX_R16_R2_OUT(index11_i,index11);
			   NTT_REORDERINDEX_R16_R2_OUT(index12_i,index12);
			   NTT_REORDERINDEX_R16_R2_OUT(index13_i,index13);
			   NTT_REORDERINDEX_R16_R2_OUT(index14_i,index14);
			   NTT_REORDERINDEX_R16_R2_OUT(index15_i,index15);
			   //DATARECORD <<" index0:  "<< index0 <<"\n";
			   //DATARECORD <<" index1:  "<< index1 <<"\n";
			   //DATARECORD <<" index2:  "<< index2 <<"\n";
			   //DATARECORD <<" index3:  "<< index3 <<"\n";
			   //DATARECORD <<" index4:  "<< index4 <<"\n";
			   //DATARECORD <<" index5:  "<< index5 <<"\n";
			   //DATARECORD <<" index6:  "<< index6 <<"\n";
			   //DATARECORD <<" index7:  "<< index7 <<"\n";
			   //DATARECORD <<" index8:  "<< index8 <<"\n";
			   //DATARECORD <<" index9:  "<< index9 <<"\n";
			   //DATARECORD <<" index10: "<< index10 <<"\n";
			   //DATARECORD <<" index11: "<< index11 <<"\n";
			   //DATARECORD <<" index12: "<< index12 <<"\n";
			   //DATARECORD <<" index13: "<< index13 <<"\n";
			   //DATARECORD <<" index14: "<< index14 <<"\n";
			   //DATARECORD <<" index15: "<< index15 <<"\n";
               A[index0] = A_B0R0[ma_tmp];
			   A[index1] = A_B0R1[ma_tmp];
			   A[index2] = A_B0R2[ma_tmp];
			   A[index3] = A_B0R3[ma_tmp];
			   A[index4] = A_B0R4[ma_tmp];
			   A[index5] = A_B0R5[ma_tmp];
			   A[index6] = A_B0R6[ma_tmp];
			   A[index7] = A_B0R7[ma_tmp];
			   A[index8] = A_B0R8[ma_tmp];
			   A[index9] = A_B0R9[ma_tmp];
			   A[index10] = A_B0R10[ma_tmp];
			   A[index11] = A_B0R11[ma_tmp];
			   A[index12] = A_B0R12[ma_tmp];
			   A[index13] = A_B0R13[ma_tmp];
			   A[index14] = A_B0R14[ma_tmp];
			   A[index15] = A_B0R15[ma_tmp];
			}
			else {
			   NTT_REORDERINDEX_R16_R2_OUT(index0_i,index0);
			   NTT_REORDERINDEX_R16_R2_OUT(index1_i,index1);
			   NTT_REORDERINDEX_R16_R2_OUT(index2_i,index2);
			   NTT_REORDERINDEX_R16_R2_OUT(index3_i,index3);
			   NTT_REORDERINDEX_R16_R2_OUT(index4_i,index4);
			   NTT_REORDERINDEX_R16_R2_OUT(index5_i,index5);
			   NTT_REORDERINDEX_R16_R2_OUT(index6_i,index6);
			   NTT_REORDERINDEX_R16_R2_OUT(index7_i,index7);
			   NTT_REORDERINDEX_R16_R2_OUT(index8_i,index8);
			   NTT_REORDERINDEX_R16_R2_OUT(index9_i,index9);
			   NTT_REORDERINDEX_R16_R2_OUT(index10_i,index10);
			   NTT_REORDERINDEX_R16_R2_OUT(index11_i,index11);
			   NTT_REORDERINDEX_R16_R2_OUT(index12_i,index12);
			   NTT_REORDERINDEX_R16_R2_OUT(index13_i,index13);
			   NTT_REORDERINDEX_R16_R2_OUT(index14_i,index14);
			   NTT_REORDERINDEX_R16_R2_OUT(index15_i,index15);
			   //DATARECORD <<" index0:  "<< index0 <<"\n";
			   //DATARECORD <<" index1:  "<< index1 <<"\n";
			   //DATARECORD <<" index2:  "<< index2 <<"\n";
			   //DATARECORD <<" index3:  "<< index3 <<"\n";
			   //DATARECORD <<" index4:  "<< index4 <<"\n";
			   //DATARECORD <<" index5:  "<< index5 <<"\n";
			   //DATARECORD <<" index6:  "<< index6 <<"\n";
			   //DATARECORD <<" index7:  "<< index7 <<"\n";
			   //DATARECORD <<" index8:  "<< index8 <<"\n";
			   //DATARECORD <<" index9:  "<< index9 <<"\n";
			   //DATARECORD <<" index10: "<< index10 <<"\n";
			   //DATARECORD <<" index11: "<< index11 <<"\n";
			   //DATARECORD <<" index12: "<< index12 <<"\n";
			   //DATARECORD <<" index13: "<< index13 <<"\n";
			   //DATARECORD <<" index14: "<< index14 <<"\n";
			   //DATARECORD <<" index15: "<< index15 <<"\n";
               A[index0]     = A_B1R0[ma_tmp];
               A[index1]     = A_B1R1[ma_tmp];
               A[index2]     = A_B1R2[ma_tmp];
               A[index3]     = A_B1R3[ma_tmp];
               A[index4]     = A_B1R4[ma_tmp];
               A[index5]     = A_B1R5[ma_tmp];
               A[index6]     = A_B1R6[ma_tmp];
               A[index7]     = A_B1R7[ma_tmp];
               A[index8]     = A_B1R8[ma_tmp];
               A[index9]     = A_B1R9[ma_tmp];
               A[index10]    = A_B1R10[ma_tmp];
               A[index11]    = A_B1R11[ma_tmp];
               A[index12]    = A_B1R12[ma_tmp];
               A[index13]    = A_B1R13[ma_tmp];
               A[index14]    = A_B1R14[ma_tmp];
               A[index15]    = A_B1R15[ma_tmp];
			}			
		}		
	}
	
	//data relocation for INTT input order
	for(int i = 0;i < group;i++){
		for(int j = 0;j < radix;j++){
	        gray_i  = Gray(i,group);
			BC_tmp  = j * group + gray_i;
	        RR_R16_R2(BC_tmp,((4 * (Stage-1)) - 3),BC);
			AGU_R16(BC,bn_tmp,ma_tmp);
			if(bn_tmp == 0){
		       data_tmp_1  = A_B0R1[ma_tmp]; 
		       data_tmp_2  = A_B0R2[ma_tmp]; 
		       data_tmp_3  = A_B0R3[ma_tmp];
		       data_tmp_4  = A_B0R4[ma_tmp];
		       data_tmp_5  = A_B0R5[ma_tmp];
		       data_tmp_6  = A_B0R6[ma_tmp];
		       data_tmp_7  = A_B0R7[ma_tmp];
		       data_tmp_8  = A_B0R8[ma_tmp];
		       data_tmp_9  = A_B0R9[ma_tmp];
		       data_tmp_10 = A_B0R10[ma_tmp];
		       data_tmp_11 = A_B0R11[ma_tmp];
		       data_tmp_12 = A_B0R12[ma_tmp];
		       data_tmp_13 = A_B0R13[ma_tmp];
		       data_tmp_14 = A_B0R14[ma_tmp];
			   A_B0R1[ma_tmp]  = data_tmp_2; 
			   A_B0R2[ma_tmp]  = data_tmp_4; 
			   A_B0R3[ma_tmp]  = data_tmp_6;
			   A_B0R4[ma_tmp]  = data_tmp_8;
			   A_B0R5[ma_tmp]  = data_tmp_10;
			   A_B0R6[ma_tmp]  = data_tmp_12;
			   A_B0R7[ma_tmp]  = data_tmp_14;
			   A_B0R8[ma_tmp]  = data_tmp_1;
			   A_B0R9[ma_tmp]  = data_tmp_3;
			   A_B0R10[ma_tmp] = data_tmp_5;
			   A_B0R11[ma_tmp] = data_tmp_7;
			   A_B0R12[ma_tmp] = data_tmp_9;
			   A_B0R13[ma_tmp] = data_tmp_11;
			   A_B0R14[ma_tmp] = data_tmp_13;
			}else {
		       data_tmp_1  = A_B1R1[ma_tmp]; 
		       data_tmp_2  = A_B1R2[ma_tmp]; 
		       data_tmp_3  = A_B1R3[ma_tmp];
		       data_tmp_4  = A_B1R4[ma_tmp];
		       data_tmp_5  = A_B1R5[ma_tmp];
		       data_tmp_6  = A_B1R6[ma_tmp];
		       data_tmp_7  = A_B1R7[ma_tmp];
		       data_tmp_8  = A_B1R8[ma_tmp];
		       data_tmp_9  = A_B1R9[ma_tmp];
		       data_tmp_10 = A_B1R10[ma_tmp];
		       data_tmp_11 = A_B1R11[ma_tmp];
		       data_tmp_12 = A_B1R12[ma_tmp];
		       data_tmp_13 = A_B1R13[ma_tmp];
		       data_tmp_14 = A_B1R14[ma_tmp];
			   A_B1R1[ma_tmp]  = data_tmp_2; 
			   A_B1R2[ma_tmp]  = data_tmp_4; 
			   A_B1R3[ma_tmp]  = data_tmp_6;
			   A_B1R4[ma_tmp]  = data_tmp_8;
			   A_B1R5[ma_tmp]  = data_tmp_10;
			   A_B1R6[ma_tmp]  = data_tmp_12;
			   A_B1R7[ma_tmp]  = data_tmp_14;
			   A_B1R8[ma_tmp]  = data_tmp_1;
			   A_B1R9[ma_tmp]  = data_tmp_3;
			   A_B1R10[ma_tmp] = data_tmp_5;
			   A_B1R11[ma_tmp] = data_tmp_7;
			   A_B1R12[ma_tmp] = data_tmp_9;
			   A_B1R13[ma_tmp] = data_tmp_11;
			   A_B1R14[ma_tmp] = data_tmp_13;
			}
		}
	}
	
	for(int ss = 0; ss < word_size; ss++){
		//bn0
		B0R0[ss]  = A_B0R0[ss];
		B0R1[ss]  = A_B0R1[ss];
		B0R2[ss]  = A_B0R2[ss];
		B0R3[ss]  = A_B0R3[ss];
		B0R4[ss]  = A_B0R4[ss];
		B0R5[ss]  = A_B0R5[ss];
		B0R6[ss]  = A_B0R6[ss];
		B0R7[ss]  = A_B0R7[ss];
		B0R8[ss]  = A_B0R8[ss];
		B0R9[ss]  = A_B0R9[ss];
		B0R10[ss] = A_B0R10[ss];
		B0R11[ss] = A_B0R11[ss];
		B0R12[ss] = A_B0R12[ss];
		B0R13[ss] = A_B0R13[ss];
		B0R14[ss] = A_B0R14[ss];
		B0R15[ss] = A_B0R15[ss];
		//----------------------------
		//bn1
		B1R0[ss]  = A_B1R0[ss];
		B1R1[ss]  = A_B1R1[ss];
		B1R2[ss]  = A_B1R2[ss];
		B1R3[ss]  = A_B1R3[ss];
		B1R4[ss]  = A_B1R4[ss];
		B1R5[ss]  = A_B1R5[ss];
		B1R6[ss]  = A_B1R6[ss];
		B1R7[ss]  = A_B1R7[ss];
		B1R8[ss]  = A_B1R8[ss];
		B1R9[ss]  = A_B1R9[ss];
		B1R10[ss] = A_B1R10[ss];
		B1R11[ss] = A_B1R11[ss];
		B1R12[ss] = A_B1R12[ss];
		B1R13[ss] = A_B1R13[ss];
		B1R14[ss] = A_B1R14[ss];
		B1R15[ss] = A_B1R15[ss];		
	}	
}
//need to modify - 2020/09/24
void NTTSPMB::INTT_r16_r2(std::vector<ZZ> &A,std::vector<ZZ> &B0R0,std::vector<ZZ> &B0R1,std::vector<ZZ> &B0R2,std::vector<ZZ> &B0R3,
std::vector<ZZ> &B0R4,std::vector<ZZ> &B0R5,std::vector<ZZ> &B0R6,std::vector<ZZ> &B0R7,std::vector<ZZ> &B0R8,std::vector<ZZ> &B0R9,
std::vector<ZZ> &B0R10,std::vector<ZZ> &B0R11,std::vector<ZZ> &B0R12,std::vector<ZZ> &B0R13,std::vector<ZZ> &B0R14,std::vector<ZZ> &B0R15,
std::vector<ZZ> &B1R0,std::vector<ZZ> &B1R1,std::vector<ZZ> &B1R2,std::vector<ZZ> &B1R3,std::vector<ZZ> &B1R4,std::vector<ZZ> &B1R5,
std::vector<ZZ> &B1R6,std::vector<ZZ> &B1R7,std::vector<ZZ> &B1R8,std::vector<ZZ> &B1R9,std::vector<ZZ> &B1R10,std::vector<ZZ> &B1R11,
std::vector<ZZ> &B1R12,std::vector<ZZ> &B1R13,std::vector<ZZ> &B1R14,std::vector<ZZ> &B1R15){
	
	int            Stage;//stage	
	int            word_size;
	int            offset;
	int            BC;
	int            BC_RR;
	int            BC_Reorder;
	int            BC_tmp;
    int            group;
    int            ma_tmp;
	int            bn_tmp;
    int            bit_tmp;
	int            bn0_bc_tmp;//frist in bc data
	int            bn1_bc_tmp;//frist in bc data
	int            bn0_ma_reg1;
	int            bn0_ma_reg2;
	int            bn0_ma_reg3;
	int            bn0_ma_reg4;
	int            bn0_ma_reg5;
	int            bn0_ma_reg6;
	int            bn0_ma_reg7;
	int            bn0_ma_reg8;
	int            bn1_ma_reg1;
	int            bn1_ma_reg2;
	int            bn1_ma_reg3;
	int            bn1_ma_reg4;
	int            bn1_ma_reg5;
	int            bn1_ma_reg6;
	int            bn1_ma_reg7;
	int            bn1_ma_reg8;
	int            gray_i;
    int            BC_WIDTH;
	int            tw_modulus;
    int            tw_modulus_tmp;
	//display      
	int            display;
    double         Stage_double;	
    std::vector<int> bit_array_tmp;
	
	std::ofstream DATARECORD("./INTT_R16_R2_SPMB.txt");
	//*********************************************************************
	//****************************display setting *************************
	display  = 1;
	//*********************************************************************
	
	//radix-16 Stage
    Stage_double  = log2(N);
    Stage       =  (int)floor(Stage_double/4);
	BC_WIDTH    =  (int)ceil(log2(N/16));
	offset      =  (int)N /  16;
	word_size   =  (int)N / (2 * 16);
	group       =  (int)N / (256);
	tw_modulus  =  (int)N /  16;
	bit_array_tmp.resize(BC_WIDTH);

    DATARECORD << "group: "    << group << "\n";
    DATARECORD << "BC_WIDTH: " << BC_WIDTH << "\n";
	
	ZZ               data_tmp_1;
	ZZ               data_tmp_2;
	ZZ               data_tmp_3;
	ZZ               data_tmp_4;
	ZZ               data_tmp_5;
	ZZ               data_tmp_6;
	ZZ               data_tmp_7;
	ZZ               data_tmp_8;
	ZZ               data_tmp_9;
	ZZ               data_tmp_10;
	ZZ               data_tmp_11;
	ZZ               data_tmp_12;
	ZZ               data_tmp_13;
	ZZ               data_tmp_14;
	ZZ               data_tmp_15;
	std::vector<ZZ>  A_B0R0;
	std::vector<ZZ>  A_B0R1;
	std::vector<ZZ>  A_B0R2;
	std::vector<ZZ>  A_B0R3;
	std::vector<ZZ>  A_B0R4;
	std::vector<ZZ>  A_B0R5;
	std::vector<ZZ>  A_B0R6;
	std::vector<ZZ>  A_B0R7;
	std::vector<ZZ>  A_B0R8;
	std::vector<ZZ>  A_B0R9;
	std::vector<ZZ>  A_B0R10;
	std::vector<ZZ>  A_B0R11;
	std::vector<ZZ>  A_B0R12;
	std::vector<ZZ>  A_B0R13;
	std::vector<ZZ>  A_B0R14;
	std::vector<ZZ>  A_B0R15;
	std::vector<ZZ>  A_B1R0;
	std::vector<ZZ>  A_B1R1;
	std::vector<ZZ>  A_B1R2;
	std::vector<ZZ>  A_B1R3;
	std::vector<ZZ>  A_B1R4;
	std::vector<ZZ>  A_B1R5;
	std::vector<ZZ>  A_B1R6;
	std::vector<ZZ>  A_B1R7;
	std::vector<ZZ>  A_B1R8;
	std::vector<ZZ>  A_B1R9;
	std::vector<ZZ>  A_B1R10;
	std::vector<ZZ>  A_B1R11;
	std::vector<ZZ>  A_B1R12;
	std::vector<ZZ>  A_B1R13;
	std::vector<ZZ>  A_B1R14;
	std::vector<ZZ>  A_B1R15;
	
	A_B0R0.resize(word_size);
	A_B0R1.resize(word_size);
	A_B0R2.resize(word_size);
	A_B0R3.resize(word_size);
	A_B0R4.resize(word_size);
	A_B0R5.resize(word_size);
	A_B0R6.resize(word_size);
	A_B0R7.resize(word_size);
	A_B0R8.resize(word_size);
	A_B0R9.resize(word_size);
	A_B0R10.resize(word_size);
	A_B0R11.resize(word_size);
	A_B0R12.resize(word_size);
	A_B0R13.resize(word_size);
	A_B0R14.resize(word_size);
	A_B0R15.resize(word_size);
	A_B1R0.resize(word_size);
	A_B1R1.resize(word_size);
	A_B1R2.resize(word_size);
	A_B1R3.resize(word_size);
	A_B1R4.resize(word_size);
	A_B1R5.resize(word_size);
	A_B1R6.resize(word_size);
	A_B1R7.resize(word_size);
	A_B1R8.resize(word_size);
	A_B1R9.resize(word_size);
	A_B1R10.resize(word_size);
	A_B1R11.resize(word_size);
	A_B1R12.resize(word_size);
	A_B1R13.resize(word_size);
	A_B1R14.resize(word_size);
	A_B1R15.resize(word_size);
	
	//----------------------------------------------------
	int length;
	ZZ  factor;   //base factor
	ZZ  factor_t; //acctually mul factor
	ZZ  factor_2t;
	ZZ  factor_3t;
	ZZ  factor_4t;
	ZZ  factor_5t;
	ZZ  factor_6t;
	ZZ  factor_7t;
	ZZ  factor_8t;
	ZZ  factor_9t;
	ZZ  factor_10t;
	ZZ  factor_11t;
	ZZ  factor_12t;
	ZZ  factor_13t;
	ZZ  factor_14t;
	ZZ  factor_15t;
    
	//init load data
    for(int i = 0; i < word_size; i++){
		A_B0R0[i]   =  B0R0[i];
        A_B0R1[i]   =  B0R1[i];
        A_B0R2[i]   =  B0R2[i];
        A_B0R3[i]   =  B0R3[i];
        A_B0R4[i]   =  B0R4[i];
        A_B0R5[i]   =  B0R5[i];
        A_B0R6[i]   =  B0R6[i];
        A_B0R7[i]   =  B0R7[i];
        A_B0R8[i]   =  B0R8[i];
        A_B0R9[i]   =  B0R9[i];
        A_B0R10[i]  =  B0R10[i];
        A_B0R11[i]  =  B0R11[i];
        A_B0R12[i]  =  B0R12[i];
        A_B0R13[i]  =  B0R13[i];
        A_B0R14[i]  =  B0R14[i];
        A_B0R15[i]  =  B0R15[i];
        A_B1R0[i]   =  B1R0[i];
        A_B1R1[i]   =  B1R1[i];
        A_B1R2[i]   =  B1R2[i];
        A_B1R3[i]   =  B1R3[i];
        A_B1R4[i]   =  B1R4[i];
        A_B1R5[i]   =  B1R5[i];
        A_B1R6[i]   =  B1R6[i];
        A_B1R7[i]   =  B1R7[i];
        A_B1R8[i]   =  B1R8[i];
        A_B1R9[i]   =  B1R9[i];
        A_B1R10[i]  =  B1R10[i];
        A_B1R11[i]  =  B1R11[i];
        A_B1R12[i]  =  B1R12[i];
        A_B1R13[i]  =  B1R13[i];
        A_B1R14[i]  =  B1R14[i];
        A_B1R15[i]  =  B1R15[i];
	}
	ma_tmp       = 0;
	bn_tmp       = 0;
	BC           = 0;
	BC_RR        = 0;
	BC_Reorder   = 0;
	std::cout << "init load over! \n";
	//need modify to mult by twiddle factor
	for(int s = 0; s < Stage; s++){
		if(s == 0)factor = IW;
		else {
			SqrMod(factor,factor,p);
			SqrMod(factor,factor,p);
			SqrMod(factor,factor,p);
			SqrMod(factor,factor,p);
		}
		tw_modulus_tmp  = tw_modulus >> ( 4 * s);
		for(int i = 0 ;i < group;i++){
			bn0_bc_tmp  = 0;
			bn1_bc_tmp  = 0;
			for(int j = 0;j < radix;j++){
				if(display == 1)DATARECORD <<"---------------------------------\n";
		        if(display == 1)DATARECORD <<"Now Stage: "<< s <<"\n";
				gray_i  = Gray(i,group);
				BC_tmp  = j * group + gray_i;
				if(display == 1)DATARECORD << "BC_tmp: "<< BC_tmp <<"\n";
				if(s == Stage - 1) RR_R16_R2(BC_tmp,(4 * s - 3),BC_RR);
				else RR_R16(BC_tmp,s,BC_RR);
				if(display == 1)DATARECORD <<"BC_RR: "<< BC_RR <<"\n";
				BC_IFFT_Reorder_R16_R2(BC_RR,BC_Reorder);
				if(display == 1)DATARECORD <<"BC_Reorder: "<< BC_Reorder <<"\n";
				length = BC_RR % tw_modulus_tmp;
				if(display == 1)DATARECORD <<"factor: "<< factor <<"\n";
				if(display == 1)DATARECORD <<"twiddle length: "<< length <<"\n";
				PowerMod(factor_t,factor,length,p);
				AGU_R16(BC_Reorder,bn_tmp,ma_tmp);
				if(display == 1)DATARECORD <<"bn_tmp: "<< bn_tmp <<"\n";
				if(display == 1)DATARECORD <<"ma_tmp: "<< ma_tmp <<"\n";
				if(bn_tmp == 0){
					if(j < 2)bn0_bc_tmp = BC_tmp;
					if(display == 1)DATARECORD <<" Before Radix-16 Butterfly unit computing!!! \n";
					if(display == 1)DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";
					PowerMod(factor_2t,factor_t,2,p);
					PowerMod(factor_3t,factor_t,3,p);
					PowerMod(factor_4t,factor_t,4,p);
					PowerMod(factor_5t,factor_t,5,p);
					PowerMod(factor_6t,factor_t,6,p);
					PowerMod(factor_7t,factor_t,7,p);
					PowerMod(factor_8t,factor_t,8,p);
					PowerMod(factor_9t,factor_t,9,p);
					PowerMod(factor_10t,factor_t,10,p);
					PowerMod(factor_11t,factor_t,11,p);
					PowerMod(factor_12t,factor_t,12,p);
					PowerMod(factor_13t,factor_t,13,p);
					PowerMod(factor_14t,factor_t,14,p);
					PowerMod(factor_15t,factor_t,15,p);
					if(display == 1)DATARECORD <<"factor_t: "<< factor_t <<"\n";
					if(display == 1)DATARECORD <<"factor_2t: "<< factor_2t <<"\n";
					if(display == 1)DATARECORD <<"factor_3t: "<< factor_3t <<"\n";
					if(display == 1)DATARECORD <<"factor_4t: "<< factor_4t <<"\n";
					if(display == 1)DATARECORD <<"factor_5t: "<< factor_5t <<"\n";
					if(display == 1)DATARECORD <<"factor_6t: "<< factor_6t <<"\n";
					if(display == 1)DATARECORD <<"factor_7t: "<< factor_7t <<"\n";
					if(display == 1)DATARECORD <<"factor_8t: "<< factor_8t <<"\n";
					if(display == 1)DATARECORD <<"factor_9t: "<< factor_9t <<"\n";
					if(display == 1)DATARECORD <<"factor_10t: "<< factor_10t <<"\n";
					if(display == 1)DATARECORD <<"factor_11t: "<< factor_11t <<"\n";
					if(display == 1)DATARECORD <<"factor_12t: "<< factor_12t <<"\n";
					if(display == 1)DATARECORD <<"factor_13t: "<< factor_13t <<"\n";
					if(display == 1)DATARECORD <<"factor_14t: "<< factor_14t <<"\n";
					if(display == 1)DATARECORD <<"factor_15t: "<< factor_15t <<"\n";
					Radix16_BU_INTT(A_B0R0[ma_tmp],A_B0R1[ma_tmp], A_B0R2[ma_tmp], A_B0R3[ma_tmp],A_B0R4[ma_tmp],
							   A_B0R5[ma_tmp],A_B0R6[ma_tmp], A_B0R7[ma_tmp], A_B0R8[ma_tmp],A_B0R9[ma_tmp],
							   A_B0R10[ma_tmp],A_B0R11[ma_tmp],A_B0R12[ma_tmp],A_B0R13[ma_tmp],A_B0R14[ma_tmp],
							   A_B0R15[ma_tmp]);
					if(display == 1)DATARECORD <<" After Radix-16 Butterfly unit computing!!! \n";
					if(display == 1)DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";		   
                    MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],factor_t,p);
					MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],factor_2t,p);
					MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],factor_3t,p);
					MulMod(A_B0R4[ma_tmp],A_B0R4[ma_tmp],factor_4t,p);
					MulMod(A_B0R5[ma_tmp],A_B0R5[ma_tmp],factor_5t,p);
					MulMod(A_B0R6[ma_tmp],A_B0R6[ma_tmp],factor_6t,p);
					MulMod(A_B0R7[ma_tmp],A_B0R7[ma_tmp],factor_7t,p);
					MulMod(A_B0R8[ma_tmp],A_B0R8[ma_tmp],factor_8t,p);
					MulMod(A_B0R9[ma_tmp],A_B0R9[ma_tmp],factor_9t,p);
					MulMod(A_B0R10[ma_tmp],A_B0R10[ma_tmp],factor_10t,p);
					MulMod(A_B0R11[ma_tmp],A_B0R11[ma_tmp],factor_11t,p);
					MulMod(A_B0R12[ma_tmp],A_B0R12[ma_tmp],factor_12t,p);
					MulMod(A_B0R13[ma_tmp],A_B0R13[ma_tmp],factor_13t,p);
					MulMod(A_B0R14[ma_tmp],A_B0R14[ma_tmp],factor_14t,p);
					MulMod(A_B0R15[ma_tmp],A_B0R15[ma_tmp],factor_15t,p);
					if(display == 1)DATARECORD <<" After mult by twiddle factor!!! \n";
					if(display == 1)DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";
					if(j <  2)bn0_ma_reg1 = ma_tmp;
					if((j >= 2)  && (j < 4))bn0_ma_reg2 = ma_tmp;
					if((j >= 4)  && (j < 6))bn0_ma_reg3 = ma_tmp;
					if((j >= 6)  && (j < 8))bn0_ma_reg4 = ma_tmp;
					if((j >= 8)  && (j < 10))bn0_ma_reg5 = ma_tmp;
					if((j >= 10) && (j < 12))bn0_ma_reg6 = ma_tmp;
					if((j >= 12) && (j < 14))bn0_ma_reg7 = ma_tmp;
					if((j >= 14) && (j < 16))bn0_ma_reg8 = ma_tmp;
				}
			    else {
					if(j < 2)bn1_bc_tmp = BC_tmp;
					if(display == 1)DATARECORD <<" Before Radix-16 Butterfly unit computing!!! \n";
					if(display == 1)DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";
					PowerMod(factor_2t,factor_t,2,p);
					PowerMod(factor_3t,factor_t,3,p);
					PowerMod(factor_4t,factor_t,4,p);
					PowerMod(factor_5t,factor_t,5,p);
					PowerMod(factor_6t,factor_t,6,p);
					PowerMod(factor_7t,factor_t,7,p);
					PowerMod(factor_8t,factor_t,8,p);
					PowerMod(factor_9t,factor_t,9,p);
					PowerMod(factor_10t,factor_t,10,p);
					PowerMod(factor_11t,factor_t,11,p);
					PowerMod(factor_12t,factor_t,12,p);
					PowerMod(factor_13t,factor_t,13,p);
					PowerMod(factor_14t,factor_t,14,p);
					PowerMod(factor_15t,factor_t,15,p);
					if(display == 1)DATARECORD <<"factor_t: "<< factor_t <<"\n";
					if(display == 1)DATARECORD <<"factor_2t: "<< factor_2t <<"\n";
					if(display == 1)DATARECORD <<"factor_3t: "<< factor_3t <<"\n";
					if(display == 1)DATARECORD <<"factor_4t: "<< factor_4t <<"\n";
					if(display == 1)DATARECORD <<"factor_5t: "<< factor_5t <<"\n";
					if(display == 1)DATARECORD <<"factor_6t: "<< factor_6t <<"\n";
					if(display == 1)DATARECORD <<"factor_7t: "<< factor_7t <<"\n";
					if(display == 1)DATARECORD <<"factor_8t: "<< factor_8t <<"\n";
					if(display == 1)DATARECORD <<"factor_9t: "<< factor_9t <<"\n";
					if(display == 1)DATARECORD <<"factor_10t: "<< factor_10t <<"\n";
					if(display == 1)DATARECORD <<"factor_11t: "<< factor_11t <<"\n";
					if(display == 1)DATARECORD <<"factor_12t: "<< factor_12t <<"\n";
					if(display == 1)DATARECORD <<"factor_13t: "<< factor_13t <<"\n";
					if(display == 1)DATARECORD <<"factor_14t: "<< factor_14t <<"\n";
					if(display == 1)DATARECORD <<"factor_15t: "<< factor_15t <<"\n";					
					Radix16_BU_INTT(A_B1R0[ma_tmp],A_B1R1[ma_tmp], A_B1R2[ma_tmp], A_B1R3[ma_tmp],A_B1R4[ma_tmp],
							   A_B1R5[ma_tmp],A_B1R6[ma_tmp], A_B1R7[ma_tmp], A_B1R8[ma_tmp],A_B1R9[ma_tmp],
							   A_B1R10[ma_tmp],A_B1R11[ma_tmp],A_B1R12[ma_tmp],A_B1R13[ma_tmp],A_B1R14[ma_tmp],
							   A_B1R15[ma_tmp]);
					if(display == 1)DATARECORD <<" After Radix-16 Butterfly unit computing!!! \n";
					if(display == 1)DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";									   
					MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],factor_t,p);
					MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],factor_2t,p);
					MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],factor_3t,p);
					MulMod(A_B1R4[ma_tmp],A_B1R4[ma_tmp],factor_4t,p);
					MulMod(A_B1R5[ma_tmp],A_B1R5[ma_tmp],factor_5t,p);
					MulMod(A_B1R6[ma_tmp],A_B1R6[ma_tmp],factor_6t,p);
					MulMod(A_B1R7[ma_tmp],A_B1R7[ma_tmp],factor_7t,p);
					MulMod(A_B1R8[ma_tmp],A_B1R8[ma_tmp],factor_8t,p);
					MulMod(A_B1R9[ma_tmp],A_B1R9[ma_tmp],factor_9t,p);
					MulMod(A_B1R10[ma_tmp],A_B1R10[ma_tmp],factor_10t,p);
					MulMod(A_B1R11[ma_tmp],A_B1R11[ma_tmp],factor_11t,p);
					MulMod(A_B1R12[ma_tmp],A_B1R12[ma_tmp],factor_12t,p);
					MulMod(A_B1R13[ma_tmp],A_B1R13[ma_tmp],factor_13t,p);
					MulMod(A_B1R14[ma_tmp],A_B1R14[ma_tmp],factor_14t,p);
					MulMod(A_B1R15[ma_tmp],A_B1R15[ma_tmp],factor_15t,p);
					if(display == 1)DATARECORD <<" After mult by twiddle factor!!! \n";
					if(display == 1)DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";					
                    if(j <  2)bn1_ma_reg1 = ma_tmp;					
                    if((j >= 2)  && (j < 4))bn1_ma_reg2 = ma_tmp;
                    if((j >= 4)  && (j < 6))bn1_ma_reg3 = ma_tmp;
                    if((j >= 6)  && (j < 8))bn1_ma_reg4 = ma_tmp;
                    if((j >= 8)  && (j < 10))bn1_ma_reg5 = ma_tmp;
                    if((j >= 10) && (j < 12))bn1_ma_reg6 = ma_tmp;
                    if((j >= 12) && (j < 14))bn1_ma_reg7 = ma_tmp;
                    if((j >= 14) && (j < 16))bn1_ma_reg8 = ma_tmp;					
				}
			}
		//data relocation
		 if(s < Stage-1){
		  if(bn1_bc_tmp > bn0_bc_tmp){
			 data_tmp_1  = A_B0R1[bn0_ma_reg1];
			 data_tmp_2  = A_B0R2[bn0_ma_reg1];
			 data_tmp_3  = A_B0R3[bn0_ma_reg1];
			 data_tmp_4  = A_B0R4[bn0_ma_reg1];
			 data_tmp_5  = A_B0R5[bn0_ma_reg1];
			 data_tmp_6  = A_B0R6[bn0_ma_reg1];
			 data_tmp_7  = A_B0R7[bn0_ma_reg1];
			 data_tmp_8  = A_B0R8[bn0_ma_reg1];
			 data_tmp_9  = A_B0R9[bn0_ma_reg1];
			 data_tmp_10 = A_B0R10[bn0_ma_reg1];
			 data_tmp_11 = A_B0R11[bn0_ma_reg1];
			 data_tmp_12 = A_B0R12[bn0_ma_reg1];
			 data_tmp_13 = A_B0R13[bn0_ma_reg1];
			 data_tmp_14 = A_B0R14[bn0_ma_reg1];
			 data_tmp_15 = A_B0R15[bn0_ma_reg1];
			 A_B0R1[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg1];
			 A_B0R2[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg2];
			 A_B0R3[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg2];
			 A_B0R4[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg3];
			 A_B0R5[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg3];
			 A_B0R6[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg1] = A_B0R0[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg1] = A_B1R0[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg1] = A_B0R0[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg1] = A_B1R0[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg1] = A_B1R0[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg1] = A_B0R0[bn0_ma_reg8];
			 A_B1R0[bn1_ma_reg1]  = data_tmp_1; 
			 A_B1R0[bn1_ma_reg2]  = data_tmp_2; 
			 A_B0R0[bn0_ma_reg2]  = data_tmp_3; 
			 A_B1R0[bn1_ma_reg3]  = data_tmp_4; 
			 A_B0R0[bn0_ma_reg3]  = data_tmp_5; 
			 A_B0R0[bn0_ma_reg4]  = data_tmp_6; 
			 A_B1R0[bn1_ma_reg4]  = data_tmp_7; 
			 A_B1R0[bn1_ma_reg5]  = data_tmp_8; 
			 A_B0R0[bn0_ma_reg5]  = data_tmp_9; 
			 A_B0R0[bn0_ma_reg6]  = data_tmp_10;
			 A_B1R0[bn1_ma_reg6]  = data_tmp_11;
			 A_B0R0[bn0_ma_reg7]  = data_tmp_12;
			 A_B1R0[bn1_ma_reg7]  = data_tmp_13;
			 A_B1R0[bn1_ma_reg8]  = data_tmp_14;
			 A_B0R0[bn0_ma_reg8]  = data_tmp_15;
			 data_tmp_1  = A_B1R2[bn1_ma_reg1];
			 data_tmp_2  = A_B1R3[bn1_ma_reg1];
			 data_tmp_3  = A_B1R4[bn1_ma_reg1];
			 data_tmp_4  = A_B1R5[bn1_ma_reg1];
			 data_tmp_5  = A_B1R6[bn1_ma_reg1];
			 data_tmp_6  = A_B1R7[bn1_ma_reg1];
			 data_tmp_7  = A_B1R8[bn1_ma_reg1];
			 data_tmp_8  = A_B1R9[bn1_ma_reg1];
			 data_tmp_9  = A_B1R10[bn1_ma_reg1];
			 data_tmp_10 = A_B1R11[bn1_ma_reg1];
			 data_tmp_11 = A_B1R12[bn1_ma_reg1];
			 data_tmp_12 = A_B1R13[bn1_ma_reg1];
			 data_tmp_13 = A_B1R14[bn1_ma_reg1];
			 data_tmp_14 = A_B1R15[bn1_ma_reg1];
             A_B1R2[bn1_ma_reg1]  =	 A_B1R1[bn1_ma_reg2];
             A_B1R3[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg2];
             A_B1R4[bn1_ma_reg1]  =  A_B1R1[bn1_ma_reg3];
             A_B1R5[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg3];
             A_B1R6[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg4];
             A_B1R7[bn1_ma_reg1]  =  A_B1R1[bn1_ma_reg4];
             A_B1R8[bn1_ma_reg1]  =  A_B1R1[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg1] =  A_B0R1[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg1] =  A_B1R1[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg1] =  A_B0R1[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg1] =  A_B1R1[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg1] =  A_B1R1[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg1] =  A_B0R1[bn0_ma_reg8];
			 A_B1R1[bn1_ma_reg2]  =  data_tmp_1; 
			 A_B0R1[bn0_ma_reg2]  =  data_tmp_2; 
			 A_B1R1[bn1_ma_reg3]  =  data_tmp_3; 
			 A_B0R1[bn0_ma_reg3]  =  data_tmp_4; 
			 A_B0R1[bn0_ma_reg4]  =  data_tmp_5; 
			 A_B1R1[bn1_ma_reg4]  =  data_tmp_6; 
			 A_B1R1[bn1_ma_reg5]  =  data_tmp_7; 
			 A_B0R1[bn0_ma_reg5]  =  data_tmp_8; 
			 A_B0R1[bn0_ma_reg6]  =  data_tmp_9; 
			 A_B1R1[bn1_ma_reg6]  =  data_tmp_10;
			 A_B0R1[bn0_ma_reg7]  =  data_tmp_11;
			 A_B1R1[bn1_ma_reg7]  =  data_tmp_12;
			 A_B1R1[bn1_ma_reg8]  =  data_tmp_13;
			 A_B0R1[bn0_ma_reg8]  =  data_tmp_14;
			/************************************************************/ 
			 data_tmp_1  =   A_B1R3[bn1_ma_reg2];
			 data_tmp_2  =   A_B1R4[bn1_ma_reg2];
			 data_tmp_3  =   A_B1R5[bn1_ma_reg2];
			 data_tmp_4  =   A_B1R6[bn1_ma_reg2];
			 data_tmp_5  =   A_B1R7[bn1_ma_reg2];
			 data_tmp_6  =   A_B1R8[bn1_ma_reg2];
			 data_tmp_7  =   A_B1R9[bn1_ma_reg2];
			 data_tmp_8  =   A_B1R10[bn1_ma_reg2];
			 data_tmp_9  =   A_B1R11[bn1_ma_reg2];
			 data_tmp_10 =   A_B1R12[bn1_ma_reg2];
			 data_tmp_11 =   A_B1R13[bn1_ma_reg2];
			 data_tmp_12 =   A_B1R14[bn1_ma_reg2];
			 data_tmp_13 =   A_B1R15[bn1_ma_reg2];
			 A_B1R3[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg2];
			 A_B1R4[bn1_ma_reg2]  =  A_B1R2[bn1_ma_reg3];
			 A_B1R5[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg3];
			 A_B1R6[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg4];
			 A_B1R7[bn1_ma_reg2]  =  A_B1R2[bn1_ma_reg4];
			 A_B1R8[bn1_ma_reg2]  =  A_B1R2[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg2] =  A_B0R2[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg2] =  A_B1R2[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg2] =  A_B0R2[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg2] =  A_B1R2[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg2] =  A_B1R2[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg2] =  A_B0R2[bn0_ma_reg8];
			 A_B0R2[bn0_ma_reg2]  =  data_tmp_1; 
			 A_B1R2[bn1_ma_reg3]  =  data_tmp_2; 
			 A_B0R2[bn0_ma_reg3]  =  data_tmp_3; 
			 A_B0R2[bn0_ma_reg4]  =  data_tmp_4; 
			 A_B1R2[bn1_ma_reg4]  =  data_tmp_5; 
			 A_B1R2[bn1_ma_reg5]  =  data_tmp_6; 
			 A_B0R2[bn0_ma_reg5]  =  data_tmp_7; 
			 A_B0R2[bn0_ma_reg6]  =  data_tmp_8; 
			 A_B1R2[bn1_ma_reg6]  =  data_tmp_9; 
			 A_B0R2[bn0_ma_reg7]  =  data_tmp_10;
			 A_B1R2[bn1_ma_reg7]  =  data_tmp_11;
			 A_B1R2[bn1_ma_reg8]  =  data_tmp_12;
			 A_B0R2[bn0_ma_reg8]  =  data_tmp_13;
			 //******************************************************
			 data_tmp_1  = A_B0R4[bn0_ma_reg2];
			 data_tmp_2  = A_B0R5[bn0_ma_reg2];
			 data_tmp_3  = A_B0R6[bn0_ma_reg2];
			 data_tmp_4  = A_B0R7[bn0_ma_reg2];
			 data_tmp_5  = A_B0R8[bn0_ma_reg2];
			 data_tmp_6  = A_B0R9[bn0_ma_reg2];
			 data_tmp_7  = A_B0R10[bn0_ma_reg2];
			 data_tmp_8  = A_B0R11[bn0_ma_reg2];
			 data_tmp_9  = A_B0R12[bn0_ma_reg2];
			 data_tmp_10 = A_B0R13[bn0_ma_reg2];
			 data_tmp_11 = A_B0R14[bn0_ma_reg2];
			 data_tmp_12 = A_B0R15[bn0_ma_reg2];
			 A_B0R4[bn0_ma_reg2]  =  A_B1R3[bn1_ma_reg3];
			 A_B0R5[bn0_ma_reg2]  =  A_B0R3[bn0_ma_reg3];
			 A_B0R6[bn0_ma_reg2]  =  A_B0R3[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg2]  =  A_B1R3[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg2]  =  A_B1R3[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg2]  =  A_B0R3[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg2] =  A_B0R3[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg2] =  A_B1R3[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg2] =  A_B0R3[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg2] =  A_B1R3[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg2] =  A_B1R3[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg2] =  A_B0R3[bn0_ma_reg8];
			 A_B1R3[bn1_ma_reg3]  =  data_tmp_1;  
			 A_B0R3[bn0_ma_reg3]  =  data_tmp_2;  
			 A_B0R3[bn0_ma_reg4]  =  data_tmp_3;  
			 A_B1R3[bn1_ma_reg4]  =  data_tmp_4;  
			 A_B1R3[bn1_ma_reg5]  =  data_tmp_5;  
			 A_B0R3[bn0_ma_reg5]  =  data_tmp_6;  
			 A_B0R3[bn0_ma_reg6]  =  data_tmp_7;  
			 A_B1R3[bn1_ma_reg6]  =  data_tmp_8;  
			 A_B0R3[bn0_ma_reg7]  =  data_tmp_9;  
			 A_B1R3[bn1_ma_reg7]  =  data_tmp_10; 
			 A_B1R3[bn1_ma_reg8]  =  data_tmp_11; 
			 A_B0R3[bn0_ma_reg8]  =  data_tmp_12; 
			 //----------------------------------------------------------------------
             data_tmp_1  = A_B1R5[bn1_ma_reg3];
			 data_tmp_2  = A_B1R6[bn1_ma_reg3];
			 data_tmp_3  = A_B1R7[bn1_ma_reg3];
			 data_tmp_4  = A_B1R8[bn1_ma_reg3];
			 data_tmp_5  = A_B1R9[bn1_ma_reg3];
			 data_tmp_6  = A_B1R10[bn1_ma_reg3];
			 data_tmp_7  = A_B1R11[bn1_ma_reg3];
			 data_tmp_8  = A_B1R12[bn1_ma_reg3];
			 data_tmp_9  = A_B1R13[bn1_ma_reg3];
			 data_tmp_10 = A_B1R14[bn1_ma_reg3];
			 data_tmp_11 = A_B1R15[bn1_ma_reg3];
			 A_B1R5[bn1_ma_reg3]  =  A_B0R4[bn0_ma_reg3];
			 A_B1R6[bn1_ma_reg3]  =  A_B0R4[bn0_ma_reg4];
			 A_B1R7[bn1_ma_reg3]  =  A_B1R4[bn1_ma_reg4];
			 A_B1R8[bn1_ma_reg3]  =  A_B1R4[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg3]  =  A_B0R4[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg3] =  A_B0R4[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg3] =  A_B1R4[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg3] =  A_B0R4[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg3] =  A_B1R4[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg3] =  A_B1R4[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg3] =  A_B0R4[bn0_ma_reg8];
			 A_B0R4[bn0_ma_reg3]  =  data_tmp_1 ;
			 A_B0R4[bn0_ma_reg4]  =  data_tmp_2 ;
			 A_B1R4[bn1_ma_reg4]  =  data_tmp_3 ;
			 A_B1R4[bn1_ma_reg5]  =  data_tmp_4 ;
			 A_B0R4[bn0_ma_reg5]  =  data_tmp_5 ;
			 A_B0R4[bn0_ma_reg6]  =  data_tmp_6 ;
			 A_B1R4[bn1_ma_reg6]  =  data_tmp_7 ;
			 A_B0R4[bn0_ma_reg7]  =  data_tmp_8 ;
			 A_B1R4[bn1_ma_reg7]  =  data_tmp_9 ;
			 A_B1R4[bn1_ma_reg8]  =  data_tmp_10;
			 A_B0R4[bn0_ma_reg8]  =  data_tmp_11;
			 //************************************************************************
			 data_tmp_1  = A_B0R6[bn0_ma_reg3];
			 data_tmp_2  = A_B0R7[bn0_ma_reg3];
			 data_tmp_3  = A_B0R8[bn0_ma_reg3];
			 data_tmp_4  = A_B0R9[bn0_ma_reg3];
			 data_tmp_5  = A_B0R10[bn0_ma_reg3];
			 data_tmp_6  = A_B0R11[bn0_ma_reg3];
			 data_tmp_7  = A_B0R12[bn0_ma_reg3];
			 data_tmp_8  = A_B0R13[bn0_ma_reg3];
			 data_tmp_9  = A_B0R14[bn0_ma_reg3];
			 data_tmp_10 = A_B0R15[bn0_ma_reg3];
			 A_B0R6[bn0_ma_reg3]  =  A_B0R5[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg3]  =  A_B1R5[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg3]  =  A_B1R5[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg3]  =  A_B0R5[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg3] =  A_B0R5[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg3] =  A_B1R5[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg3] =  A_B0R5[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg3] =  A_B1R5[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg3] =  A_B1R5[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg3] =  A_B0R5[bn0_ma_reg8];
			 A_B0R5[bn0_ma_reg4]  =  data_tmp_1; 
			 A_B1R5[bn1_ma_reg4]  =  data_tmp_2; 
			 A_B1R5[bn1_ma_reg5]  =  data_tmp_3; 
			 A_B0R5[bn0_ma_reg5]  =  data_tmp_4; 
			 A_B0R5[bn0_ma_reg6]  =  data_tmp_5; 
			 A_B1R5[bn1_ma_reg6]  =  data_tmp_6; 
			 A_B0R5[bn0_ma_reg7]  =  data_tmp_7; 
			 A_B1R5[bn1_ma_reg7]  =  data_tmp_8; 
			 A_B1R5[bn1_ma_reg8]  =  data_tmp_9; 
			 A_B0R5[bn0_ma_reg8]  =  data_tmp_10;
			 //-----------------------------------------------------------------------
			 data_tmp_1  = A_B0R7[bn0_ma_reg4];
			 data_tmp_2  = A_B0R8[bn0_ma_reg4];
			 data_tmp_3  = A_B0R9[bn0_ma_reg4];
			 data_tmp_4  = A_B0R10[bn0_ma_reg4];
			 data_tmp_5  = A_B0R11[bn0_ma_reg4];
			 data_tmp_6  = A_B0R12[bn0_ma_reg4];
			 data_tmp_7  = A_B0R13[bn0_ma_reg4];
			 data_tmp_8  = A_B0R14[bn0_ma_reg4];
			 data_tmp_9  = A_B0R15[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg4]  =  A_B1R6[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg4]  =  A_B1R6[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg4]  =  A_B0R6[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg4] =  A_B0R6[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg4] =  A_B1R6[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg4] =  A_B0R6[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg4] =  A_B1R6[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg4] =  A_B1R6[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg4] =  A_B0R6[bn0_ma_reg8];
			 A_B1R6[bn1_ma_reg4]  =  data_tmp_1;
			 A_B1R6[bn1_ma_reg5]  =  data_tmp_2;
			 A_B0R6[bn0_ma_reg5]  =  data_tmp_3;
			 A_B0R6[bn0_ma_reg6]  =  data_tmp_4;
			 A_B1R6[bn1_ma_reg6]  =  data_tmp_5;
			 A_B0R6[bn0_ma_reg7]  =  data_tmp_6;
			 A_B1R6[bn1_ma_reg7]  =  data_tmp_7;
			 A_B1R6[bn1_ma_reg8]  =  data_tmp_8;
			 A_B0R6[bn0_ma_reg8]  =  data_tmp_9;
			 //----------------------------------------------------------------------
			 data_tmp_1  = A_B1R8[bn1_ma_reg4];
			 data_tmp_2  = A_B1R9[bn1_ma_reg4];
			 data_tmp_3  = A_B1R10[bn1_ma_reg4];
			 data_tmp_4  = A_B1R11[bn1_ma_reg4];
			 data_tmp_5  = A_B1R12[bn1_ma_reg4];
			 data_tmp_6  = A_B1R13[bn1_ma_reg4];
			 data_tmp_7  = A_B1R14[bn1_ma_reg4];
			 data_tmp_8  = A_B1R15[bn1_ma_reg4];
			 A_B1R8[bn1_ma_reg4]  =  A_B1R7[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg4]  =  A_B0R7[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg4] =  A_B0R7[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg4] =  A_B1R7[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg4] =  A_B0R7[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg4] =  A_B1R7[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg4] =  A_B1R7[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg4] =  A_B0R7[bn0_ma_reg8];
			 A_B1R7[bn1_ma_reg5]  =  data_tmp_1;
			 A_B0R7[bn0_ma_reg5]  =  data_tmp_2;
			 A_B0R7[bn0_ma_reg6]  =  data_tmp_3;
			 A_B1R7[bn1_ma_reg6]  =  data_tmp_4;
			 A_B0R7[bn0_ma_reg7]  =  data_tmp_5;
			 A_B1R7[bn1_ma_reg7]  =  data_tmp_6;
			 A_B1R7[bn1_ma_reg8]  =  data_tmp_7;
			 A_B0R7[bn0_ma_reg8]  =  data_tmp_8;
			 //----------------------------------------------------------------------
			 data_tmp_1 = A_B1R9[bn1_ma_reg5];
			 data_tmp_2 = A_B1R10[bn1_ma_reg5];
			 data_tmp_3 = A_B1R11[bn1_ma_reg5];
			 data_tmp_4 = A_B1R12[bn1_ma_reg5];
			 data_tmp_5 = A_B1R13[bn1_ma_reg5];
			 data_tmp_6 = A_B1R14[bn1_ma_reg5];
			 data_tmp_7 = A_B1R15[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg5]  = A_B0R8[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg5] = A_B0R8[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg5] = A_B1R8[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg5] = A_B0R8[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg5] = A_B1R8[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg5] = A_B1R8[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg5] = A_B0R8[bn0_ma_reg8];
			 A_B0R8[bn0_ma_reg5]  = data_tmp_1;
			 A_B0R8[bn0_ma_reg6]  = data_tmp_2;
			 A_B1R8[bn1_ma_reg6]  = data_tmp_3;
			 A_B0R8[bn0_ma_reg7]  = data_tmp_4;
			 A_B1R8[bn1_ma_reg7]  = data_tmp_5;
			 A_B1R8[bn1_ma_reg8]  = data_tmp_6;
			 A_B0R8[bn0_ma_reg8]  = data_tmp_7;
			 //---------------------------------------------------------------------
			 data_tmp_1  = A_B0R10[bn0_ma_reg5];
			 data_tmp_2  = A_B0R11[bn0_ma_reg5];
			 data_tmp_3  = A_B0R12[bn0_ma_reg5];
			 data_tmp_4  = A_B0R13[bn0_ma_reg5];
			 data_tmp_5  = A_B0R14[bn0_ma_reg5];
			 data_tmp_6  = A_B0R15[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg5] = A_B0R9[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg5] = A_B1R9[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg5] = A_B0R9[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg5] = A_B1R9[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg5] = A_B1R9[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg5] = A_B0R9[bn0_ma_reg8];
			 A_B0R9[bn0_ma_reg6]  = data_tmp_1;
			 A_B1R9[bn1_ma_reg6]  = data_tmp_2;
			 A_B0R9[bn0_ma_reg7]  = data_tmp_3;
			 A_B1R9[bn1_ma_reg7]  = data_tmp_4;
			 A_B1R9[bn1_ma_reg8]  = data_tmp_5;
			 A_B0R9[bn0_ma_reg8]  = data_tmp_6;
			 //--------------------------------------------------------------------
			 data_tmp_1  = A_B0R11[bn0_ma_reg6];
			 data_tmp_2  = A_B0R12[bn0_ma_reg6];
			 data_tmp_3  = A_B0R13[bn0_ma_reg6];
			 data_tmp_4  = A_B0R14[bn0_ma_reg6];
			 data_tmp_5  = A_B0R15[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg6] = A_B1R10[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg6] = A_B0R10[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg6] = A_B1R10[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg6] = A_B1R10[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg6] = A_B0R10[bn0_ma_reg8];
			 A_B1R10[bn1_ma_reg6] = data_tmp_1;
			 A_B0R10[bn0_ma_reg7] = data_tmp_2;
			 A_B1R10[bn1_ma_reg7] = data_tmp_3;
			 A_B1R10[bn1_ma_reg8] = data_tmp_4;
			 A_B0R10[bn0_ma_reg8] = data_tmp_5;
			 //--------------------------------------------------------------------
			 data_tmp_1  = A_B1R12[bn1_ma_reg6];
			 data_tmp_2  = A_B1R13[bn1_ma_reg6];
			 data_tmp_3  = A_B1R14[bn1_ma_reg6];
			 data_tmp_4  = A_B1R15[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg6] = A_B0R11[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg6] = A_B1R11[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg6] = A_B1R11[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg6] = A_B0R11[bn0_ma_reg8];
			 A_B0R11[bn0_ma_reg7] = data_tmp_1;
			 A_B1R11[bn1_ma_reg7] = data_tmp_2;
			 A_B1R11[bn1_ma_reg8] = data_tmp_3;
			 A_B0R11[bn0_ma_reg8] = data_tmp_4;
			 //--------------------------------------------------------------------
			 data_tmp_1 = A_B0R13[bn0_ma_reg7];
			 data_tmp_2 = A_B0R14[bn0_ma_reg7];
			 data_tmp_3 = A_B0R15[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg7] = A_B1R12[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg7] = A_B1R12[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg7] = A_B0R12[bn0_ma_reg8];
			 A_B1R12[bn1_ma_reg7] = data_tmp_1;
			 A_B1R12[bn1_ma_reg8] = data_tmp_2;
			 A_B0R12[bn0_ma_reg8] = data_tmp_3;
			 //--------------------------------------------------------------------
			 data_tmp_1 = A_B1R14[bn1_ma_reg7];
			 data_tmp_2 = A_B1R15[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg7] = A_B1R13[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg7] = A_B0R13[bn0_ma_reg8];
			 A_B1R13[bn1_ma_reg8] = data_tmp_1;
			 A_B0R13[bn0_ma_reg8] = data_tmp_2;
			 //--------------------------------------------------------------------
			 data_tmp_1 = A_B1R15[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg8] = A_B0R14[bn0_ma_reg8];
			 A_B0R14[bn0_ma_reg8] = data_tmp_1; 
		  }else{	 
			 data_tmp_1  = A_B1R1[bn1_ma_reg1];
			 data_tmp_2  = A_B1R2[bn1_ma_reg1];
			 data_tmp_3  = A_B1R3[bn1_ma_reg1];
			 data_tmp_4  = A_B1R4[bn1_ma_reg1];
			 data_tmp_5  = A_B1R5[bn1_ma_reg1];
			 data_tmp_6  = A_B1R6[bn1_ma_reg1];
			 data_tmp_7  = A_B1R7[bn1_ma_reg1];
			 data_tmp_8  = A_B1R8[bn1_ma_reg1];
			 data_tmp_9  = A_B1R9[bn1_ma_reg1];
			 data_tmp_10 = A_B1R10[bn1_ma_reg1];
			 data_tmp_11 = A_B1R11[bn1_ma_reg1];
			 data_tmp_12 = A_B1R12[bn1_ma_reg1];
			 data_tmp_13 = A_B1R13[bn1_ma_reg1];
			 data_tmp_14 = A_B1R14[bn1_ma_reg1];
			 data_tmp_15 = A_B1R15[bn1_ma_reg1];
             A_B1R1[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg1];
             A_B1R2[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg2];
             A_B1R3[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg2];
             A_B1R4[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg3];
             A_B1R5[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg3];
             A_B1R6[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg4];
             A_B1R7[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg4];
             A_B1R8[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg5];
             A_B1R9[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg5];
             A_B1R10[bn1_ma_reg1] =  A_B1R0[bn1_ma_reg6];
             A_B1R11[bn1_ma_reg1] =  A_B0R0[bn0_ma_reg6];
             A_B1R12[bn1_ma_reg1] =  A_B1R0[bn1_ma_reg7];
             A_B1R13[bn1_ma_reg1] =  A_B0R0[bn0_ma_reg7];
             A_B1R14[bn1_ma_reg1] =  A_B0R0[bn0_ma_reg8];
             A_B1R15[bn1_ma_reg1] =  A_B1R0[bn1_ma_reg8];
             A_B0R0[bn0_ma_reg1]  =  data_tmp_1; 
             A_B0R0[bn0_ma_reg2]  =  data_tmp_2; 
             A_B1R0[bn1_ma_reg2]  =  data_tmp_3; 
             A_B0R0[bn0_ma_reg3]  =  data_tmp_4; 
             A_B1R0[bn1_ma_reg3]  =  data_tmp_5; 
             A_B1R0[bn1_ma_reg4]  =  data_tmp_6; 
             A_B0R0[bn0_ma_reg4]  =  data_tmp_7; 
             A_B0R0[bn0_ma_reg5]  =  data_tmp_8; 
             A_B1R0[bn1_ma_reg5]  =  data_tmp_9; 
             A_B1R0[bn1_ma_reg6]  =  data_tmp_10;
             A_B0R0[bn0_ma_reg6]  =  data_tmp_11;
             A_B1R0[bn1_ma_reg7]  =  data_tmp_12;
             A_B0R0[bn0_ma_reg7]  =  data_tmp_13;
             A_B0R0[bn0_ma_reg8]  =  data_tmp_14;
             A_B1R0[bn1_ma_reg8]  =  data_tmp_15;
             //-----------------------------------------------------------
             data_tmp_1  = A_B0R2[bn0_ma_reg1];
			 data_tmp_2  = A_B0R3[bn0_ma_reg1];
			 data_tmp_3  = A_B0R4[bn0_ma_reg1];
			 data_tmp_4  = A_B0R5[bn0_ma_reg1];
			 data_tmp_5  = A_B0R6[bn0_ma_reg1];
			 data_tmp_6  = A_B0R7[bn0_ma_reg1];
			 data_tmp_7  = A_B0R8[bn0_ma_reg1];
			 data_tmp_8  = A_B0R9[bn0_ma_reg1];
			 data_tmp_9  = A_B0R10[bn0_ma_reg1];
			 data_tmp_10 = A_B0R11[bn0_ma_reg1];
			 data_tmp_11 = A_B0R12[bn0_ma_reg1];
			 data_tmp_12 = A_B0R13[bn0_ma_reg1];
			 data_tmp_13 = A_B0R14[bn0_ma_reg1];
			 data_tmp_14 = A_B0R15[bn0_ma_reg1];
			 A_B0R2[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg2];
			 A_B0R3[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg2];
			 A_B0R4[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg3];
			 A_B0R5[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg3];
			 A_B0R6[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg4];
			 A_B0R7[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg1] =  A_B1R1[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg1] =  A_B0R1[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg1] =  A_B1R1[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg1] =  A_B0R1[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg1] =  A_B0R1[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg1] =  A_B1R1[bn1_ma_reg8];
             A_B0R1[bn0_ma_reg2]  =  data_tmp_1; 
             A_B1R1[bn1_ma_reg2]  =  data_tmp_2; 
             A_B0R1[bn0_ma_reg3]  =  data_tmp_3; 
             A_B1R1[bn1_ma_reg3]  =  data_tmp_4; 
             A_B1R1[bn1_ma_reg4]  =  data_tmp_5; 
             A_B0R1[bn0_ma_reg4]  =  data_tmp_6; 
             A_B0R1[bn0_ma_reg5]  =  data_tmp_7; 
             A_B1R1[bn1_ma_reg5]  =  data_tmp_8; 
             A_B1R1[bn1_ma_reg6]  =  data_tmp_9; 
             A_B0R1[bn0_ma_reg6]  =  data_tmp_10;
             A_B1R1[bn1_ma_reg7]  =  data_tmp_11;
             A_B0R1[bn0_ma_reg7]  =  data_tmp_12;
             A_B0R1[bn0_ma_reg8]  =  data_tmp_13;
             A_B1R1[bn1_ma_reg8]  =  data_tmp_14;
             //------------------------------------------------------------
             data_tmp_1  =  A_B0R3[bn0_ma_reg2];
			 data_tmp_2  =  A_B0R4[bn0_ma_reg2];
			 data_tmp_3  =  A_B0R5[bn0_ma_reg2];
			 data_tmp_4  =  A_B0R6[bn0_ma_reg2];
			 data_tmp_5  =  A_B0R7[bn0_ma_reg2];
			 data_tmp_6  =  A_B0R8[bn0_ma_reg2];
			 data_tmp_7  =  A_B0R9[bn0_ma_reg2];
			 data_tmp_8  =  A_B0R10[bn0_ma_reg2];
			 data_tmp_9  =  A_B0R11[bn0_ma_reg2];
			 data_tmp_10 =  A_B0R12[bn0_ma_reg2];
			 data_tmp_11 =  A_B0R13[bn0_ma_reg2];
			 data_tmp_12 =  A_B0R14[bn0_ma_reg2];
			 data_tmp_13 =  A_B0R15[bn0_ma_reg2];
			 A_B0R3[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg2];
			 A_B0R4[bn0_ma_reg2]  =  A_B0R2[bn0_ma_reg3];
			 A_B0R5[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg3];
			 A_B0R6[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg4];
			 A_B0R7[bn0_ma_reg2]  =  A_B0R2[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg2]  =  A_B0R2[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg2] =  A_B1R2[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg2] =  A_B0R2[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg2] =  A_B1R2[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg2] =  A_B0R2[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg2] =  A_B0R2[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg2] =  A_B1R2[bn1_ma_reg8];
			 A_B1R2[bn1_ma_reg2]  =  data_tmp_1; 
			 A_B0R2[bn0_ma_reg3]  =  data_tmp_2; 
			 A_B1R2[bn1_ma_reg3]  =  data_tmp_3; 
			 A_B1R2[bn1_ma_reg4]  =  data_tmp_4; 
			 A_B0R2[bn0_ma_reg4]  =  data_tmp_5; 
			 A_B0R2[bn0_ma_reg5]  =  data_tmp_6; 
			 A_B1R2[bn1_ma_reg5]  =  data_tmp_7; 
			 A_B1R2[bn1_ma_reg6]  =  data_tmp_8; 
			 A_B0R2[bn0_ma_reg6]  =  data_tmp_9; 
			 A_B1R2[bn1_ma_reg7]  =  data_tmp_10;
			 A_B0R2[bn0_ma_reg7]  =  data_tmp_11;
			 A_B0R2[bn0_ma_reg8]  =  data_tmp_12;
			 A_B1R2[bn1_ma_reg8]  =  data_tmp_13;
			 //-----------------------------------------------------------                        
			 data_tmp_1  = A_B1R4[bn1_ma_reg2];
			 data_tmp_2  = A_B1R5[bn1_ma_reg2];
			 data_tmp_3  = A_B1R6[bn1_ma_reg2];
			 data_tmp_4  = A_B1R7[bn1_ma_reg2];
			 data_tmp_5  = A_B1R8[bn1_ma_reg2];
			 data_tmp_6  = A_B1R9[bn1_ma_reg2];
			 data_tmp_7  = A_B1R10[bn1_ma_reg2];
			 data_tmp_8  = A_B1R11[bn1_ma_reg2];
			 data_tmp_9  = A_B1R12[bn1_ma_reg2];
			 data_tmp_10 = A_B1R13[bn1_ma_reg2];
			 data_tmp_11 = A_B1R14[bn1_ma_reg2];
			 data_tmp_12 = A_B1R15[bn1_ma_reg2];
			 A_B1R4[bn1_ma_reg2]  = A_B0R3[bn0_ma_reg3];
			 A_B1R5[bn1_ma_reg2]  = A_B1R3[bn1_ma_reg3];
			 A_B1R6[bn1_ma_reg2]  = A_B1R3[bn1_ma_reg4];
			 A_B1R7[bn1_ma_reg2]  = A_B0R3[bn0_ma_reg4];
			 A_B1R8[bn1_ma_reg2]  = A_B0R3[bn0_ma_reg5];
			 A_B1R9[bn1_ma_reg2]  = A_B1R3[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg2] = A_B1R3[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg2] = A_B0R3[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg2] = A_B1R3[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg2] = A_B0R3[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg2] = A_B0R3[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg2] = A_B1R3[bn1_ma_reg8];
			 A_B0R3[bn0_ma_reg3]  = data_tmp_1; 
			 A_B1R3[bn1_ma_reg3]  = data_tmp_2; 
			 A_B1R3[bn1_ma_reg4]  = data_tmp_3; 
			 A_B0R3[bn0_ma_reg4]  = data_tmp_4; 
			 A_B0R3[bn0_ma_reg5]  = data_tmp_5; 
			 A_B1R3[bn1_ma_reg5]  = data_tmp_6; 
			 A_B1R3[bn1_ma_reg6]  = data_tmp_7; 
			 A_B0R3[bn0_ma_reg6]  = data_tmp_8; 
			 A_B1R3[bn1_ma_reg7]  = data_tmp_9; 
			 A_B0R3[bn0_ma_reg7]  = data_tmp_10;
			 A_B0R3[bn0_ma_reg8]  = data_tmp_11;
			 A_B1R3[bn1_ma_reg8]  = data_tmp_12;
			 //------------------------------------------------------------
			 data_tmp_1  =  A_B0R5[bn0_ma_reg3];
			 data_tmp_2  =  A_B0R6[bn0_ma_reg3];
			 data_tmp_3  =  A_B0R7[bn0_ma_reg3];
			 data_tmp_4  =  A_B0R8[bn0_ma_reg3];
			 data_tmp_5  =  A_B0R9[bn0_ma_reg3];
			 data_tmp_6  =  A_B0R10[bn0_ma_reg3];
			 data_tmp_7  =  A_B0R11[bn0_ma_reg3];
			 data_tmp_8  =  A_B0R12[bn0_ma_reg3];
			 data_tmp_9  =  A_B0R13[bn0_ma_reg3];
			 data_tmp_10 =  A_B0R14[bn0_ma_reg3];
			 data_tmp_11 =  A_B0R15[bn0_ma_reg3];
			 A_B0R5[bn0_ma_reg3]  =  A_B1R4[bn1_ma_reg3];
			 A_B0R6[bn0_ma_reg3]  =  A_B1R4[bn1_ma_reg4];
			 A_B0R7[bn0_ma_reg3]  =  A_B0R4[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg3]  =  A_B0R4[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg3]  =  A_B1R4[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg3] =  A_B1R4[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg3] =  A_B0R4[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg3] =  A_B1R4[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg3] =  A_B0R4[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg3] =  A_B0R4[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg3] =  A_B1R4[bn1_ma_reg8];
			 A_B1R4[bn1_ma_reg3]  =  data_tmp_1; 
			 A_B1R4[bn1_ma_reg4]  =  data_tmp_2; 
			 A_B0R4[bn0_ma_reg4]  =  data_tmp_3; 
			 A_B0R4[bn0_ma_reg5]  =  data_tmp_4; 
			 A_B1R4[bn1_ma_reg5]  =  data_tmp_5; 
			 A_B1R4[bn1_ma_reg6]  =  data_tmp_6; 
			 A_B0R4[bn0_ma_reg6]  =  data_tmp_7; 
			 A_B1R4[bn1_ma_reg7]  =  data_tmp_8; 
			 A_B0R4[bn0_ma_reg7]  =  data_tmp_9; 
			 A_B0R4[bn0_ma_reg8]  =  data_tmp_10;
			 A_B1R4[bn1_ma_reg8]  =  data_tmp_11;
			 //-------------------------------------------------------------
			 data_tmp_1  =  A_B1R6[bn1_ma_reg3];
			 data_tmp_2  =  A_B1R7[bn1_ma_reg3];
			 data_tmp_3  =  A_B1R8[bn1_ma_reg3];
			 data_tmp_4  =  A_B1R9[bn1_ma_reg3];
			 data_tmp_5  =  A_B1R10[bn1_ma_reg3];
			 data_tmp_6  =  A_B1R11[bn1_ma_reg3];
			 data_tmp_7  =  A_B1R12[bn1_ma_reg3];
			 data_tmp_8  =  A_B1R13[bn1_ma_reg3];
			 data_tmp_9  =  A_B1R14[bn1_ma_reg3];
			 data_tmp_10 =  A_B1R15[bn1_ma_reg3];
			 A_B1R6[bn1_ma_reg3]  =  A_B1R5[bn1_ma_reg4];
			 A_B1R7[bn1_ma_reg3]  =  A_B0R5[bn0_ma_reg4];
			 A_B1R8[bn1_ma_reg3]  =  A_B0R5[bn0_ma_reg5];
			 A_B1R9[bn1_ma_reg3]  =  A_B1R5[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg3] =  A_B1R5[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg3] =  A_B0R5[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg3] =  A_B1R5[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg3] =  A_B0R5[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg3] =  A_B0R5[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg3] =  A_B1R5[bn1_ma_reg8];
			 A_B1R5[bn1_ma_reg4]  =  data_tmp_1; 
			 A_B0R5[bn0_ma_reg4]  =  data_tmp_2; 
			 A_B0R5[bn0_ma_reg5]  =  data_tmp_3; 
			 A_B1R5[bn1_ma_reg5]  =  data_tmp_4; 
			 A_B1R5[bn1_ma_reg6]  =  data_tmp_5; 
			 A_B0R5[bn0_ma_reg6]  =  data_tmp_6; 
			 A_B1R5[bn1_ma_reg7]  =  data_tmp_7; 
			 A_B0R5[bn0_ma_reg7]  =  data_tmp_8; 
			 A_B0R5[bn0_ma_reg8]  =  data_tmp_9; 
			 A_B1R5[bn1_ma_reg8]  =  data_tmp_10;
			 //----------------------------------------------------------------
			 data_tmp_1  =  A_B1R7[bn1_ma_reg4];
			 data_tmp_2  =  A_B1R8[bn1_ma_reg4];
			 data_tmp_3  =  A_B1R9[bn1_ma_reg4];
			 data_tmp_4  =  A_B1R10[bn1_ma_reg4];
			 data_tmp_5  =  A_B1R11[bn1_ma_reg4];
			 data_tmp_6  =  A_B1R12[bn1_ma_reg4];
			 data_tmp_7  =  A_B1R13[bn1_ma_reg4];
			 data_tmp_8  =  A_B1R14[bn1_ma_reg4];
			 data_tmp_9  =  A_B1R15[bn1_ma_reg4];
			 A_B1R7[bn1_ma_reg4]  =  A_B0R6[bn0_ma_reg4];
			 A_B1R8[bn1_ma_reg4]  =  A_B0R6[bn0_ma_reg5];
			 A_B1R9[bn1_ma_reg4]  =  A_B1R6[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg4] =  A_B1R6[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg4] =  A_B0R6[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg4] =  A_B1R6[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg4] =  A_B0R6[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg4] =  A_B0R6[bn0_ma_reg8];
             A_B1R15[bn1_ma_reg4] =  A_B1R6[bn1_ma_reg8];
             A_B0R6[bn0_ma_reg4]  =  data_tmp_1;
             A_B0R6[bn0_ma_reg5]  =  data_tmp_2;
             A_B1R6[bn1_ma_reg5]  =  data_tmp_3;
             A_B1R6[bn1_ma_reg6]  =  data_tmp_4;
             A_B0R6[bn0_ma_reg6]  =  data_tmp_5;
             A_B1R6[bn1_ma_reg7]  =  data_tmp_6;
             A_B0R6[bn0_ma_reg7]  =  data_tmp_7;
             A_B0R6[bn0_ma_reg8]  =  data_tmp_8;
             A_B1R6[bn1_ma_reg8]  =  data_tmp_9;
             //------------------------------------------------------------------
			 data_tmp_1  = A_B0R8[bn0_ma_reg4];
			 data_tmp_2  = A_B0R9[bn0_ma_reg4];
			 data_tmp_3  = A_B0R10[bn0_ma_reg4];
			 data_tmp_4  = A_B0R11[bn0_ma_reg4];
			 data_tmp_5  = A_B0R12[bn0_ma_reg4];
			 data_tmp_6  = A_B0R13[bn0_ma_reg4];
			 data_tmp_7  = A_B0R14[bn0_ma_reg4];
			 data_tmp_8  = A_B0R15[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg4]  =  A_B0R7[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg4]  =  A_B1R7[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg4] =  A_B1R7[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg4] =  A_B0R7[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg4] =  A_B1R7[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg4] =  A_B0R7[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg4] =  A_B0R7[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg4] =  A_B1R7[bn1_ma_reg8];
			 A_B0R7[bn0_ma_reg5]  =  data_tmp_1;
			 A_B1R7[bn1_ma_reg5]  =  data_tmp_2;
			 A_B1R7[bn1_ma_reg6]  =  data_tmp_3;
			 A_B0R7[bn0_ma_reg6]  =  data_tmp_4;
			 A_B1R7[bn1_ma_reg7]  =  data_tmp_5;
			 A_B0R7[bn0_ma_reg7]  =  data_tmp_6;
			 A_B0R7[bn0_ma_reg8]  =  data_tmp_7;
			 A_B1R7[bn1_ma_reg8]  =  data_tmp_8;
			 //----------------------------------------------------------------
			 data_tmp_1  = A_B0R9[bn0_ma_reg5];
			 data_tmp_2  = A_B0R10[bn0_ma_reg5];
			 data_tmp_3  = A_B0R11[bn0_ma_reg5];
			 data_tmp_4  = A_B0R12[bn0_ma_reg5];
			 data_tmp_5  = A_B0R13[bn0_ma_reg5];
			 data_tmp_6  = A_B0R14[bn0_ma_reg5];
			 data_tmp_7  = A_B0R15[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg5]  =  A_B1R8[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg8];
			 A_B1R8[bn1_ma_reg5]  =  data_tmp_1;
			 A_B1R8[bn1_ma_reg6]  =  data_tmp_2;
			 A_B0R8[bn0_ma_reg6]  =  data_tmp_3;
			 A_B1R8[bn1_ma_reg7]  =  data_tmp_4;
			 A_B0R8[bn0_ma_reg7]  =  data_tmp_5;
			 A_B0R8[bn0_ma_reg8]  =  data_tmp_6;
			 A_B1R8[bn1_ma_reg8]  =  data_tmp_7;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B1R10[bn1_ma_reg5];
			 data_tmp_2 = A_B1R11[bn1_ma_reg5];
			 data_tmp_3 = A_B1R12[bn1_ma_reg5];
			 data_tmp_4 = A_B1R13[bn1_ma_reg5];
			 data_tmp_5 = A_B1R14[bn1_ma_reg5];
			 data_tmp_6 = A_B1R15[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg5] = A_B1R9[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg5] = A_B0R9[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg5] = A_B1R9[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg5] = A_B0R9[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg5] = A_B0R9[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg5] = A_B1R9[bn1_ma_reg8];
			 A_B1R9[bn1_ma_reg6]  = data_tmp_1;
			 A_B0R9[bn0_ma_reg6]  = data_tmp_2;
			 A_B1R9[bn1_ma_reg7]  = data_tmp_3;
			 A_B0R9[bn0_ma_reg7]  = data_tmp_4;
			 A_B0R9[bn0_ma_reg8]  = data_tmp_5;
			 A_B1R9[bn1_ma_reg8]  = data_tmp_6;
			 //----------------------------------------------------------------
			 data_tmp_1 =  A_B1R11[bn1_ma_reg6];
			 data_tmp_2 =  A_B1R12[bn1_ma_reg6];
			 data_tmp_3 =  A_B1R13[bn1_ma_reg6];
			 data_tmp_4 =  A_B1R14[bn1_ma_reg6];
			 data_tmp_5 =  A_B1R15[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg6]  = A_B0R10[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg6]  = A_B1R10[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg6]  = A_B0R10[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg6]  = A_B0R10[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg6]  = A_B1R10[bn1_ma_reg8];
			 A_B0R10[bn0_ma_reg6]  = data_tmp_1;
			 A_B1R10[bn1_ma_reg7]  = data_tmp_2;
			 A_B0R10[bn0_ma_reg7]  = data_tmp_3;
			 A_B0R10[bn0_ma_reg8]  = data_tmp_4;
			 A_B1R10[bn1_ma_reg8]  = data_tmp_5;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B0R12[bn0_ma_reg6];
			 data_tmp_2 = A_B0R13[bn0_ma_reg6];
			 data_tmp_3 = A_B0R14[bn0_ma_reg6];
			 data_tmp_4 = A_B0R15[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg6] = A_B1R11[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg6] = A_B0R11[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg6] = A_B0R11[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg6] = A_B1R11[bn1_ma_reg8];
			 A_B1R11[bn1_ma_reg7] = data_tmp_1;
			 A_B0R11[bn0_ma_reg7] = data_tmp_2;
			 A_B0R11[bn0_ma_reg8] = data_tmp_3;
			 A_B1R11[bn1_ma_reg8] = data_tmp_4;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B1R13[bn1_ma_reg7];
			 data_tmp_2 = A_B1R14[bn1_ma_reg7];
			 data_tmp_3 = A_B1R15[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg7] = A_B0R12[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg7] = A_B0R12[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg7] = A_B1R12[bn1_ma_reg8];
			 A_B0R12[bn0_ma_reg7] = data_tmp_1;
			 A_B0R12[bn0_ma_reg8] = data_tmp_2;
			 A_B1R12[bn1_ma_reg8] = data_tmp_3;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B0R14[bn0_ma_reg7];
			 data_tmp_2 = A_B0R15[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg7] = A_B0R13[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg7] = A_B1R13[bn1_ma_reg8];
			 A_B0R13[bn0_ma_reg8] = data_tmp_1;
			 A_B1R13[bn1_ma_reg8] = data_tmp_2;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B0R15[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg8] = A_B1R14[bn1_ma_reg8];
			 A_B1R14[bn1_ma_reg8] = data_tmp_1;
			 //-----------------------------------------------------------------
		  }
		 }else{
		  if(bn1_bc_tmp > bn0_bc_tmp){
			 //Ex 0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30
			 //Ex 1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31
			 data_tmp_1  = A_B0R1[bn0_ma_reg1];
			 A_B0R1[bn0_ma_reg1] = A_B1R0[bn1_ma_reg1];
             A_B1R0[bn1_ma_reg1] = data_tmp_1;
			 data_tmp_1  = A_B0R3[bn0_ma_reg1];
			 A_B0R3[bn0_ma_reg1] = A_B1R2[bn1_ma_reg1];
             A_B1R2[bn1_ma_reg1] = data_tmp_1;
			 data_tmp_1  = A_B0R5[bn0_ma_reg1];
			 A_B0R5[bn0_ma_reg1] = A_B1R4[bn1_ma_reg1];
             A_B1R4[bn1_ma_reg1] = data_tmp_1;        
			 data_tmp_1  = A_B0R7[bn0_ma_reg1];
			 A_B0R7[bn0_ma_reg1] = A_B1R6[bn1_ma_reg1];
             A_B1R6[bn1_ma_reg1] = data_tmp_1;
			 data_tmp_1  = A_B0R9[bn0_ma_reg1];
			 A_B0R9[bn0_ma_reg1] = A_B1R8[bn1_ma_reg1];
             A_B1R8[bn1_ma_reg1] = data_tmp_1;			 
			 data_tmp_1  = A_B0R11[bn0_ma_reg1];
			 A_B0R11[bn0_ma_reg1] = A_B1R10[bn1_ma_reg1];
             A_B1R10[bn1_ma_reg1] = data_tmp_1;
			 data_tmp_1  = A_B0R13[bn0_ma_reg1];
			 A_B0R13[bn0_ma_reg1] = A_B1R12[bn1_ma_reg1];
             A_B1R12[bn1_ma_reg1] = data_tmp_1;
			 data_tmp_1  = A_B0R15[bn0_ma_reg1];
			 A_B0R15[bn0_ma_reg1] = A_B1R14[bn1_ma_reg1];
             A_B1R14[bn1_ma_reg1] = data_tmp_1;
			 /*********************************************/
			 data_tmp_1  = A_B1R1[bn1_ma_reg2];
			 A_B1R1[bn1_ma_reg2] = A_B0R0[bn0_ma_reg2];
             A_B0R0[bn0_ma_reg2] = data_tmp_1;			 
			 data_tmp_1  = A_B1R3[bn1_ma_reg2];
			 A_B1R3[bn1_ma_reg2] = A_B0R2[bn0_ma_reg2];
             A_B0R2[bn0_ma_reg2] = data_tmp_1;				 
			 data_tmp_1  = A_B1R5[bn1_ma_reg2];
			 A_B1R5[bn1_ma_reg2] = A_B0R4[bn0_ma_reg2];
             A_B0R4[bn0_ma_reg2] = data_tmp_1;	
			 data_tmp_1  = A_B1R7[bn1_ma_reg2];
			 A_B1R7[bn1_ma_reg2] = A_B0R6[bn0_ma_reg2];
             A_B0R6[bn0_ma_reg2] = data_tmp_1;	
			 data_tmp_1  = A_B1R9[bn1_ma_reg2];
			 A_B1R9[bn1_ma_reg2] = A_B0R8[bn0_ma_reg2];
             A_B0R8[bn0_ma_reg2] = data_tmp_1;	
			 data_tmp_1  = A_B1R11[bn1_ma_reg2];
			 A_B1R11[bn1_ma_reg2] = A_B0R10[bn0_ma_reg2];
             A_B0R10[bn0_ma_reg2] = data_tmp_1;	
			 data_tmp_1  = A_B1R13[bn1_ma_reg2];
			 A_B1R13[bn1_ma_reg2] = A_B0R12[bn0_ma_reg2];
             A_B0R12[bn0_ma_reg2] = data_tmp_1;	
			 data_tmp_1  = A_B1R15[bn1_ma_reg2];
			 A_B1R15[bn1_ma_reg2] = A_B0R14[bn0_ma_reg2];
             A_B0R14[bn0_ma_reg2] = data_tmp_1;	
			 /*********************************************/			 
			 data_tmp_1  = A_B1R1[bn1_ma_reg3];
			 A_B1R1[bn1_ma_reg3] = A_B0R0[bn0_ma_reg3];
             A_B0R0[bn0_ma_reg3] = data_tmp_1;			 
			 data_tmp_1  = A_B1R3[bn1_ma_reg3];
			 A_B1R3[bn1_ma_reg3] = A_B0R2[bn0_ma_reg3];
             A_B0R2[bn0_ma_reg3] = data_tmp_1;				 
			 data_tmp_1  = A_B1R5[bn1_ma_reg3];
			 A_B1R5[bn1_ma_reg3] = A_B0R4[bn0_ma_reg3];
             A_B0R4[bn0_ma_reg3] = data_tmp_1;	
			 data_tmp_1  = A_B1R7[bn1_ma_reg3];
			 A_B1R7[bn1_ma_reg3] = A_B0R6[bn0_ma_reg3];
             A_B0R6[bn0_ma_reg3] = data_tmp_1;	
			 data_tmp_1  = A_B1R9[bn1_ma_reg3];
			 A_B1R9[bn1_ma_reg3] = A_B0R8[bn0_ma_reg3];
             A_B0R8[bn0_ma_reg3] = data_tmp_1;	
			 data_tmp_1  = A_B1R11[bn1_ma_reg3];
			 A_B1R11[bn1_ma_reg3] = A_B0R10[bn0_ma_reg3];
             A_B0R10[bn0_ma_reg3] = data_tmp_1;	
			 data_tmp_1  = A_B1R13[bn1_ma_reg3];
			 A_B1R13[bn1_ma_reg3] = A_B0R12[bn0_ma_reg3];
             A_B0R12[bn0_ma_reg3] = data_tmp_1;	
			 data_tmp_1  = A_B1R15[bn1_ma_reg3];
			 A_B1R15[bn1_ma_reg3] = A_B0R14[bn0_ma_reg3];
             A_B0R14[bn0_ma_reg3] = data_tmp_1;	
			 /*********************************************/				 
			 data_tmp_1  = A_B0R1[bn0_ma_reg4];
			 A_B0R1[bn0_ma_reg4] = A_B1R0[bn1_ma_reg4];
             A_B1R0[bn1_ma_reg4] = data_tmp_1;
			 data_tmp_1  = A_B0R3[bn0_ma_reg4];
			 A_B0R3[bn0_ma_reg4] = A_B1R2[bn1_ma_reg4];
             A_B1R2[bn1_ma_reg4] = data_tmp_1;
			 data_tmp_1  = A_B0R5[bn0_ma_reg4];
			 A_B0R5[bn0_ma_reg4] = A_B1R4[bn1_ma_reg4];
             A_B1R4[bn1_ma_reg4] = data_tmp_1;        
			 data_tmp_1  = A_B0R7[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg4] = A_B1R6[bn1_ma_reg4];
             A_B1R6[bn1_ma_reg4] = data_tmp_1;
			 data_tmp_1  = A_B0R9[bn0_ma_reg4];
			 A_B0R9[bn0_ma_reg4] = A_B1R8[bn1_ma_reg4];
             A_B1R8[bn1_ma_reg4] = data_tmp_1;			 
			 data_tmp_1  = A_B0R11[bn0_ma_reg4];
			 A_B0R11[bn0_ma_reg4] = A_B1R10[bn1_ma_reg4];
             A_B1R10[bn1_ma_reg4] = data_tmp_1;
			 data_tmp_1  = A_B0R13[bn0_ma_reg4];
			 A_B0R13[bn0_ma_reg4] = A_B1R12[bn1_ma_reg4];
             A_B1R12[bn1_ma_reg4] = data_tmp_1;
			 data_tmp_1  = A_B0R15[bn0_ma_reg4];
			 A_B0R15[bn0_ma_reg4] = A_B1R14[bn1_ma_reg4];
             A_B1R14[bn1_ma_reg4] = data_tmp_1;
			 /*********************************************/
			 data_tmp_1  = A_B1R1[bn1_ma_reg5];
			 A_B1R1[bn1_ma_reg5] = A_B0R0[bn0_ma_reg5];
             A_B0R0[bn0_ma_reg5] = data_tmp_1;			 
			 data_tmp_1  = A_B1R3[bn1_ma_reg5];
			 A_B1R3[bn1_ma_reg5] = A_B0R2[bn0_ma_reg5];
             A_B0R2[bn0_ma_reg5] = data_tmp_1;				 
			 data_tmp_1  = A_B1R5[bn1_ma_reg5];
			 A_B1R5[bn1_ma_reg5] = A_B0R4[bn0_ma_reg5];
             A_B0R4[bn0_ma_reg5] = data_tmp_1;	
			 data_tmp_1  = A_B1R7[bn1_ma_reg5];
			 A_B1R7[bn1_ma_reg5] = A_B0R6[bn0_ma_reg5];
             A_B0R6[bn0_ma_reg5] = data_tmp_1;	
			 data_tmp_1  = A_B1R9[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg5] = A_B0R8[bn0_ma_reg5];
             A_B0R8[bn0_ma_reg5] = data_tmp_1;	
			 data_tmp_1  = A_B1R11[bn1_ma_reg5];
			 A_B1R11[bn1_ma_reg5] = A_B0R10[bn0_ma_reg5];
             A_B0R10[bn0_ma_reg5] = data_tmp_1;	
			 data_tmp_1  = A_B1R13[bn1_ma_reg5];
			 A_B1R13[bn1_ma_reg5] = A_B0R12[bn0_ma_reg5];
             A_B0R12[bn0_ma_reg5] = data_tmp_1;	
			 data_tmp_1  = A_B1R15[bn1_ma_reg5];
			 A_B1R15[bn1_ma_reg5] = A_B0R14[bn0_ma_reg5];
             A_B0R14[bn0_ma_reg5] = data_tmp_1;	
			 /*********************************************/
			 data_tmp_1  = A_B0R1[bn0_ma_reg6];
			 A_B0R1[bn0_ma_reg6] = A_B1R0[bn1_ma_reg6];
             A_B1R0[bn1_ma_reg6] = data_tmp_1;
			 data_tmp_1  = A_B0R3[bn0_ma_reg6];
			 A_B0R3[bn0_ma_reg6] = A_B1R2[bn1_ma_reg6];
             A_B1R2[bn1_ma_reg6] = data_tmp_1;
			 data_tmp_1  = A_B0R5[bn0_ma_reg6];
			 A_B0R5[bn0_ma_reg6] = A_B1R4[bn1_ma_reg6];
             A_B1R4[bn1_ma_reg6] = data_tmp_1;        
			 data_tmp_1  = A_B0R7[bn0_ma_reg6];
			 A_B0R7[bn0_ma_reg6] = A_B1R6[bn1_ma_reg6];
             A_B1R6[bn1_ma_reg6] = data_tmp_1;
			 data_tmp_1  = A_B0R9[bn0_ma_reg6];
			 A_B0R9[bn0_ma_reg6] = A_B1R8[bn1_ma_reg6];
             A_B1R8[bn1_ma_reg6] = data_tmp_1;			 
			 data_tmp_1  = A_B0R11[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg6] = A_B1R10[bn1_ma_reg6];
             A_B1R10[bn1_ma_reg6] = data_tmp_1;
			 data_tmp_1  = A_B0R13[bn0_ma_reg6];
			 A_B0R13[bn0_ma_reg6] = A_B1R12[bn1_ma_reg6];
             A_B1R12[bn1_ma_reg6] = data_tmp_1;
			 data_tmp_1  = A_B0R15[bn0_ma_reg6];
			 A_B0R15[bn0_ma_reg6] = A_B1R14[bn1_ma_reg6];
             A_B1R14[bn1_ma_reg6] = data_tmp_1;
			 /*********************************************/	
			 data_tmp_1  = A_B0R1[bn0_ma_reg7];
			 A_B0R1[bn0_ma_reg7] = A_B1R0[bn1_ma_reg7];
             A_B1R0[bn1_ma_reg7] = data_tmp_1;
			 data_tmp_1  = A_B0R3[bn0_ma_reg7];
			 A_B0R3[bn0_ma_reg7] = A_B1R2[bn1_ma_reg7];
             A_B1R2[bn1_ma_reg7] = data_tmp_1;
			 data_tmp_1  = A_B0R5[bn0_ma_reg7];
			 A_B0R5[bn0_ma_reg7] = A_B1R4[bn1_ma_reg7];
             A_B1R4[bn1_ma_reg7] = data_tmp_1;        
			 data_tmp_1  = A_B0R7[bn0_ma_reg7];
			 A_B0R7[bn0_ma_reg7] = A_B1R6[bn1_ma_reg7];
             A_B1R6[bn1_ma_reg7] = data_tmp_1;
			 data_tmp_1  = A_B0R9[bn0_ma_reg7];
			 A_B0R9[bn0_ma_reg7] = A_B1R8[bn1_ma_reg7];
             A_B1R8[bn1_ma_reg7] = data_tmp_1;			 
			 data_tmp_1  = A_B0R11[bn0_ma_reg7];
			 A_B0R11[bn0_ma_reg7] = A_B1R10[bn1_ma_reg7];
             A_B1R10[bn1_ma_reg7] = data_tmp_1;
			 data_tmp_1  = A_B0R13[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg7] = A_B1R12[bn1_ma_reg7];
             A_B1R12[bn1_ma_reg7] = data_tmp_1;
			 data_tmp_1  = A_B0R15[bn0_ma_reg7];
			 A_B0R15[bn0_ma_reg7] = A_B1R14[bn1_ma_reg7];
             A_B1R14[bn1_ma_reg7] = data_tmp_1;
			 /*********************************************/
			 data_tmp_1  = A_B1R1[bn1_ma_reg8];
			 A_B1R1[bn1_ma_reg8] = A_B0R0[bn0_ma_reg8];
             A_B0R0[bn0_ma_reg8] = data_tmp_1;			 
			 data_tmp_1  = A_B1R3[bn1_ma_reg8];
			 A_B1R3[bn1_ma_reg8] = A_B0R2[bn0_ma_reg8];
             A_B0R2[bn0_ma_reg8] = data_tmp_1;				 
			 data_tmp_1  = A_B1R5[bn1_ma_reg8];
			 A_B1R5[bn1_ma_reg8] = A_B0R4[bn0_ma_reg8];
             A_B0R4[bn0_ma_reg8] = data_tmp_1;	
			 data_tmp_1  = A_B1R7[bn1_ma_reg8];
			 A_B1R7[bn1_ma_reg8] = A_B0R6[bn0_ma_reg8];
             A_B0R6[bn0_ma_reg8] = data_tmp_1;	
			 data_tmp_1  = A_B1R9[bn1_ma_reg8];
			 A_B1R9[bn1_ma_reg8] = A_B0R8[bn0_ma_reg8];
             A_B0R8[bn0_ma_reg8] = data_tmp_1;	
			 data_tmp_1  = A_B1R11[bn1_ma_reg8];
			 A_B1R11[bn1_ma_reg8] = A_B0R10[bn0_ma_reg8];
             A_B0R10[bn0_ma_reg8] = data_tmp_1;	
			 data_tmp_1  = A_B1R13[bn1_ma_reg8];
			 A_B1R13[bn1_ma_reg8] = A_B0R12[bn0_ma_reg8];
             A_B0R12[bn0_ma_reg8] = data_tmp_1;	
			 data_tmp_1  = A_B1R15[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg8] = A_B0R14[bn0_ma_reg8];
             A_B0R14[bn0_ma_reg8] = data_tmp_1;	
			 /*********************************************/	
		  }else{	 
			 data_tmp_1  = A_B1R1[bn1_ma_reg1];
			 A_B1R1[bn1_ma_reg1] = A_B0R0[bn0_ma_reg1];
             A_B0R0[bn0_ma_reg1] = data_tmp_1;			 
			 data_tmp_1  = A_B1R3[bn1_ma_reg1];
			 A_B1R3[bn1_ma_reg1] = A_B0R2[bn0_ma_reg1];
             A_B0R2[bn0_ma_reg1] = data_tmp_1;				 
			 data_tmp_1  = A_B1R5[bn1_ma_reg1];
			 A_B1R5[bn1_ma_reg1] = A_B0R4[bn0_ma_reg1];
             A_B0R4[bn0_ma_reg1] = data_tmp_1;	
			 data_tmp_1  = A_B1R7[bn1_ma_reg1];
			 A_B1R7[bn1_ma_reg1] = A_B0R6[bn0_ma_reg1];
             A_B0R6[bn0_ma_reg1] = data_tmp_1;	
			 data_tmp_1  = A_B1R9[bn1_ma_reg1];
			 A_B1R9[bn1_ma_reg1] = A_B0R8[bn0_ma_reg1];
             A_B0R8[bn0_ma_reg1] = data_tmp_1;	
			 data_tmp_1  = A_B1R11[bn1_ma_reg1];
			 A_B1R11[bn1_ma_reg1] = A_B0R10[bn0_ma_reg1];
             A_B0R10[bn0_ma_reg1] = data_tmp_1;	
			 data_tmp_1  = A_B1R13[bn1_ma_reg1];
			 A_B1R13[bn1_ma_reg1] = A_B0R12[bn0_ma_reg1];
             A_B0R12[bn0_ma_reg1] = data_tmp_1;	
			 data_tmp_1  = A_B1R15[bn1_ma_reg1];
			 A_B1R15[bn1_ma_reg1] = A_B0R14[bn0_ma_reg1];
             A_B0R14[bn0_ma_reg1] = data_tmp_1;	
			 /*********************************************/
			 data_tmp_1  = A_B0R1[bn0_ma_reg2];
			 A_B0R1[bn0_ma_reg2] = A_B1R0[bn1_ma_reg2];
             A_B1R0[bn1_ma_reg2] = data_tmp_1;
			 data_tmp_1  = A_B0R3[bn0_ma_reg2];
			 A_B0R3[bn0_ma_reg2] = A_B1R2[bn1_ma_reg2];
             A_B1R2[bn1_ma_reg2] = data_tmp_1;
			 data_tmp_1  = A_B0R5[bn0_ma_reg2];
			 A_B0R5[bn0_ma_reg2] = A_B1R4[bn1_ma_reg2];
             A_B1R4[bn1_ma_reg2] = data_tmp_1;        
			 data_tmp_1  = A_B0R7[bn0_ma_reg2];
			 A_B0R7[bn0_ma_reg2] = A_B1R6[bn1_ma_reg2];
             A_B1R6[bn1_ma_reg2] = data_tmp_1;
			 data_tmp_1  = A_B0R9[bn0_ma_reg2];
			 A_B0R9[bn0_ma_reg2] = A_B1R8[bn1_ma_reg2];
             A_B1R8[bn1_ma_reg2] = data_tmp_1;			 
			 data_tmp_1  = A_B0R11[bn0_ma_reg2];
			 A_B0R11[bn0_ma_reg2] = A_B1R10[bn1_ma_reg2];
             A_B1R10[bn1_ma_reg2] = data_tmp_1;
			 data_tmp_1  = A_B0R13[bn0_ma_reg2];
			 A_B0R13[bn0_ma_reg2] = A_B1R12[bn1_ma_reg2];
             A_B1R12[bn1_ma_reg2] = data_tmp_1;
			 data_tmp_1  = A_B0R15[bn0_ma_reg2];
			 A_B0R15[bn0_ma_reg2] = A_B1R14[bn1_ma_reg2];
             A_B1R14[bn1_ma_reg2] = data_tmp_1;
			 /*********************************************/	
			 data_tmp_1  = A_B0R1[bn0_ma_reg3];
			 A_B0R1[bn0_ma_reg3] = A_B1R0[bn1_ma_reg3];
             A_B1R0[bn1_ma_reg3] = data_tmp_1;
			 data_tmp_1  = A_B0R3[bn0_ma_reg3];
			 A_B0R3[bn0_ma_reg3] = A_B1R2[bn1_ma_reg3];
             A_B1R2[bn1_ma_reg3] = data_tmp_1;
			 data_tmp_1  = A_B0R5[bn0_ma_reg3];
			 A_B0R5[bn0_ma_reg3] = A_B1R4[bn1_ma_reg3];
             A_B1R4[bn1_ma_reg3] = data_tmp_1;        
			 data_tmp_1  = A_B0R7[bn0_ma_reg3];
			 A_B0R7[bn0_ma_reg3] = A_B1R6[bn1_ma_reg3];
             A_B1R6[bn1_ma_reg3] = data_tmp_1;
			 data_tmp_1  = A_B0R9[bn0_ma_reg3];
			 A_B0R9[bn0_ma_reg3] = A_B1R8[bn1_ma_reg3];
             A_B1R8[bn1_ma_reg3] = data_tmp_1;			 
			 data_tmp_1  = A_B0R11[bn0_ma_reg3];
			 A_B0R11[bn0_ma_reg3] = A_B1R10[bn1_ma_reg3];
             A_B1R10[bn1_ma_reg3] = data_tmp_1;
			 data_tmp_1  = A_B0R13[bn0_ma_reg3];
			 A_B0R13[bn0_ma_reg3] = A_B1R12[bn1_ma_reg3];
             A_B1R12[bn1_ma_reg3] = data_tmp_1;
			 data_tmp_1  = A_B0R15[bn0_ma_reg3];
			 A_B0R15[bn0_ma_reg3] = A_B1R14[bn1_ma_reg3];
             A_B1R14[bn1_ma_reg3] = data_tmp_1;
			 /*********************************************/
			 data_tmp_1  = A_B1R1[bn1_ma_reg4];
			 A_B1R1[bn1_ma_reg4] = A_B0R0[bn0_ma_reg4];
             A_B0R0[bn0_ma_reg4] = data_tmp_1;			 
			 data_tmp_1  = A_B1R3[bn1_ma_reg4];
			 A_B1R3[bn1_ma_reg4] = A_B0R2[bn0_ma_reg4];
             A_B0R2[bn0_ma_reg4] = data_tmp_1;				 
			 data_tmp_1  = A_B1R5[bn1_ma_reg4];
			 A_B1R5[bn1_ma_reg4] = A_B0R4[bn0_ma_reg4];
             A_B0R4[bn0_ma_reg4] = data_tmp_1;	
			 data_tmp_1  = A_B1R7[bn1_ma_reg4];
			 A_B1R7[bn1_ma_reg4] = A_B0R6[bn0_ma_reg4];
             A_B0R6[bn0_ma_reg4] = data_tmp_1;	
			 data_tmp_1  = A_B1R9[bn1_ma_reg4];
			 A_B1R9[bn1_ma_reg4] = A_B0R8[bn0_ma_reg4];
             A_B0R8[bn0_ma_reg4] = data_tmp_1;	
			 data_tmp_1  = A_B1R11[bn1_ma_reg4];
			 A_B1R11[bn1_ma_reg4] = A_B0R10[bn0_ma_reg4];
             A_B0R10[bn0_ma_reg4] = data_tmp_1;	
			 data_tmp_1  = A_B1R13[bn1_ma_reg4];
			 A_B1R13[bn1_ma_reg4] = A_B0R12[bn0_ma_reg4];
             A_B0R12[bn0_ma_reg4] = data_tmp_1;	
			 data_tmp_1  = A_B1R15[bn1_ma_reg4];
			 A_B1R15[bn1_ma_reg4] = A_B0R14[bn0_ma_reg4];
             A_B0R14[bn0_ma_reg4] = data_tmp_1;	
			 /*********************************************/	
			 data_tmp_1  = A_B0R1[bn0_ma_reg5];
			 A_B0R1[bn0_ma_reg5] = A_B1R0[bn1_ma_reg5];
             A_B1R0[bn1_ma_reg5] = data_tmp_1;
			 data_tmp_1  = A_B0R3[bn0_ma_reg5];
			 A_B0R3[bn0_ma_reg5] = A_B1R2[bn1_ma_reg5];
             A_B1R2[bn1_ma_reg5] = data_tmp_1;
			 data_tmp_1  = A_B0R5[bn0_ma_reg5];
			 A_B0R5[bn0_ma_reg5] = A_B1R4[bn1_ma_reg5];
             A_B1R4[bn1_ma_reg5] = data_tmp_1;        
			 data_tmp_1  = A_B0R7[bn0_ma_reg5];
			 A_B0R7[bn0_ma_reg5] = A_B1R6[bn1_ma_reg5];
             A_B1R6[bn1_ma_reg5] = data_tmp_1;
			 data_tmp_1  = A_B0R9[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg5] = A_B1R8[bn1_ma_reg5];
             A_B1R8[bn1_ma_reg5] = data_tmp_1;			 
			 data_tmp_1  = A_B0R11[bn0_ma_reg5];
			 A_B0R11[bn0_ma_reg5] = A_B1R10[bn1_ma_reg5];
             A_B1R10[bn1_ma_reg5] = data_tmp_1;
			 data_tmp_1  = A_B0R13[bn0_ma_reg5];
			 A_B0R13[bn0_ma_reg5] = A_B1R12[bn1_ma_reg5];
             A_B1R12[bn1_ma_reg5] = data_tmp_1;
			 data_tmp_1  = A_B0R15[bn0_ma_reg5];
			 A_B0R15[bn0_ma_reg5] = A_B1R14[bn1_ma_reg5];
             A_B1R14[bn1_ma_reg5] = data_tmp_1;
			 /*********************************************/
			 data_tmp_1  = A_B1R1[bn1_ma_reg6];
			 A_B1R1[bn1_ma_reg6] = A_B0R0[bn0_ma_reg6];
             A_B0R0[bn0_ma_reg6] = data_tmp_1;			 
			 data_tmp_1  = A_B1R3[bn1_ma_reg6];
			 A_B1R3[bn1_ma_reg6] = A_B0R2[bn0_ma_reg6];
             A_B0R2[bn0_ma_reg6] = data_tmp_1;				 
			 data_tmp_1  = A_B1R5[bn1_ma_reg6];
			 A_B1R5[bn1_ma_reg6] = A_B0R4[bn0_ma_reg6];
             A_B0R4[bn0_ma_reg6] = data_tmp_1;	
			 data_tmp_1  = A_B1R7[bn1_ma_reg6];
			 A_B1R7[bn1_ma_reg6] = A_B0R6[bn0_ma_reg6];
             A_B0R6[bn0_ma_reg6] = data_tmp_1;	
			 data_tmp_1  = A_B1R9[bn1_ma_reg6];
			 A_B1R9[bn1_ma_reg6] = A_B0R8[bn0_ma_reg6];
             A_B0R8[bn0_ma_reg6] = data_tmp_1;	
			 data_tmp_1  = A_B1R11[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg6] = A_B0R10[bn0_ma_reg6];
             A_B0R10[bn0_ma_reg6] = data_tmp_1;	
			 data_tmp_1  = A_B1R13[bn1_ma_reg6];
			 A_B1R13[bn1_ma_reg6] = A_B0R12[bn0_ma_reg6];
             A_B0R12[bn0_ma_reg6] = data_tmp_1;	
			 data_tmp_1  = A_B1R15[bn1_ma_reg6];
			 A_B1R15[bn1_ma_reg6] = A_B0R14[bn0_ma_reg6];
             A_B0R14[bn0_ma_reg6] = data_tmp_1;	
			 /*********************************************/	
			 data_tmp_1  = A_B1R1[bn1_ma_reg7];
			 A_B1R1[bn1_ma_reg7] = A_B0R0[bn0_ma_reg7];
             A_B0R0[bn0_ma_reg7] = data_tmp_1;			 
			 data_tmp_1  = A_B1R3[bn1_ma_reg7];
			 A_B1R3[bn1_ma_reg7] = A_B0R2[bn0_ma_reg7];
             A_B0R2[bn0_ma_reg7] = data_tmp_1;				 
			 data_tmp_1  = A_B1R5[bn1_ma_reg7];
			 A_B1R5[bn1_ma_reg7] = A_B0R4[bn0_ma_reg7];
             A_B0R4[bn0_ma_reg7] = data_tmp_1;	
			 data_tmp_1  = A_B1R7[bn1_ma_reg7];
			 A_B1R7[bn1_ma_reg7] = A_B0R6[bn0_ma_reg7];
             A_B0R6[bn0_ma_reg7] = data_tmp_1;	
			 data_tmp_1  = A_B1R9[bn1_ma_reg7];
			 A_B1R9[bn1_ma_reg7] = A_B0R8[bn0_ma_reg7];
             A_B0R8[bn0_ma_reg7] = data_tmp_1;	
			 data_tmp_1  = A_B1R11[bn1_ma_reg7];
			 A_B1R11[bn1_ma_reg7] = A_B0R10[bn0_ma_reg7];
             A_B0R10[bn0_ma_reg7] = data_tmp_1;	
			 data_tmp_1  = A_B1R13[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg7] = A_B0R12[bn0_ma_reg7];
             A_B0R12[bn0_ma_reg7] = data_tmp_1;	
			 data_tmp_1  = A_B1R15[bn1_ma_reg7];
			 A_B1R15[bn1_ma_reg7] = A_B0R14[bn0_ma_reg7];
             A_B0R14[bn0_ma_reg7] = data_tmp_1;	
			 /*********************************************/	
			 data_tmp_1  = A_B0R1[bn0_ma_reg8];
			 A_B0R1[bn0_ma_reg8] = A_B1R0[bn1_ma_reg8];
             A_B1R0[bn1_ma_reg8] = data_tmp_1;
			 data_tmp_1  = A_B0R3[bn0_ma_reg8];
			 A_B0R3[bn0_ma_reg8] = A_B1R2[bn1_ma_reg8];
             A_B1R2[bn1_ma_reg8] = data_tmp_1;
			 data_tmp_1  = A_B0R5[bn0_ma_reg8];
			 A_B0R5[bn0_ma_reg8] = A_B1R4[bn1_ma_reg8];
             A_B1R4[bn1_ma_reg8] = data_tmp_1;        
			 data_tmp_1  = A_B0R7[bn0_ma_reg8];
			 A_B0R7[bn0_ma_reg8] = A_B1R6[bn1_ma_reg8];
             A_B1R6[bn1_ma_reg8] = data_tmp_1;
			 data_tmp_1  = A_B0R9[bn0_ma_reg8];
			 A_B0R9[bn0_ma_reg8] = A_B1R8[bn1_ma_reg8];
             A_B1R8[bn1_ma_reg8] = data_tmp_1;			 
			 data_tmp_1  = A_B0R11[bn0_ma_reg8];
			 A_B0R11[bn0_ma_reg8] = A_B1R10[bn1_ma_reg8];
             A_B1R10[bn1_ma_reg8] = data_tmp_1;
			 data_tmp_1  = A_B0R13[bn0_ma_reg8];
			 A_B0R13[bn0_ma_reg8] = A_B1R12[bn1_ma_reg8];
             A_B1R12[bn1_ma_reg8] = data_tmp_1;
			 data_tmp_1  = A_B0R15[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg8] = A_B1R14[bn1_ma_reg8];
             A_B1R14[bn1_ma_reg8] = data_tmp_1;			 
             //-----------------------------------------------------------
		    }	
		 } 
		}
	}
	
	std::cout << "radix-16 FFT computing over!!\n";
	
	if(display == 1)DATARECORD <<" ----------------------------------------------- \n";
	if(display == 1)DATARECORD <<" Radix-2 FFT computing start!!! \n";
	
	
	for(int i = 0;i < group; i++){
        for(int j=0; j < radix; j++){
        	gray_i  = Gray(i,group);
        	BC_tmp  = j * group + gray_i;
        	BC_RR = BC_tmp;
			BC_IFFT_Reorder_R16_R2(BC_RR,BC_Reorder);
        	AGU_R16(BC_Reorder,bn_tmp,ma_tmp);
			if(display == 1)DATARECORD <<"--------------------------------------------- \n";
			if(display == 1)DATARECORD <<"BC_tmp: "<< BC_tmp <<"\n";
			if(display == 1)DATARECORD <<"BC_RR: " << BC_RR  <<"\n";
			if(display == 1)DATARECORD <<"BC_Reorder: "<< BC_Reorder <<"\n";
			if(display == 1)DATARECORD <<"bn_tmp: "<< bn_tmp <<"\n";
			if(display == 1)DATARECORD <<"ma_tmp: "<< ma_tmp <<"\n";			
        	if(bn_tmp == 0){
				if(display == 1)DATARECORD <<" Before Radix-2 Butterfly unit computing!!! \n";
				if(display == 1)DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";				
        		Radix2_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp]);
        		Radix2_BU(A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
        		Radix2_BU(A_B0R4[ma_tmp],A_B0R5[ma_tmp]);
        		Radix2_BU(A_B0R6[ma_tmp],A_B0R7[ma_tmp]);
        		Radix2_BU(A_B0R8[ma_tmp],A_B0R9[ma_tmp]);
        		Radix2_BU(A_B0R10[ma_tmp],A_B0R11[ma_tmp]);
        		Radix2_BU(A_B0R12[ma_tmp],A_B0R13[ma_tmp]);
        		Radix2_BU(A_B0R14[ma_tmp],A_B0R15[ma_tmp]);
				if(display == 1)DATARECORD <<" After Radix-2 Butterfly unit computing!!! \n";
                if(display == 1)DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";			
        	}else {
				if(display == 1)DATARECORD <<" Before Radix-2 Butterfly unit computing!!! \n";
				if(display == 1)DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";					
        	    Radix2_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp]);
        	    Radix2_BU(A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
        	    Radix2_BU(A_B1R4[ma_tmp],A_B1R5[ma_tmp]);
        	    Radix2_BU(A_B1R6[ma_tmp],A_B1R7[ma_tmp]);
        	    Radix2_BU(A_B1R8[ma_tmp],A_B1R9[ma_tmp]);
        	    Radix2_BU(A_B1R10[ma_tmp],A_B1R11[ma_tmp]);
        	    Radix2_BU(A_B1R12[ma_tmp],A_B1R13[ma_tmp]);
        	    Radix2_BU(A_B1R14[ma_tmp],A_B1R15[ma_tmp]);
				if(display == 1)DATARECORD <<" After Radix-2 Butterfly unit computing!!! \n";
                if(display == 1)DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";				
        	}			
        }	
	}
	
	int index0_i;
	int index1_i;
	int index2_i;
	int index3_i;
	int index4_i;
	int index5_i;
	int index6_i;
	int index7_i;
	int index8_i;
	int index9_i;
	int index10_i;
	int index11_i;
	int index12_i;
	int index13_i;
	int index14_i;
	int index15_i;
	
	int index0;
    int index1;
    int index2;
    int index3;
    int index4;
    int index5;
    int index6;
    int index7;
    int index8;
    int index9;
    int index10;
    int index11;
    int index12;
    int index13;
    int index14;
    int index15;
	
	DATARECORD <<"***************************************************************\n";
	DATARECORD <<"***** INTT OUTPUT!!                                         ***\n";
	DATARECORD <<"***************************************************************\n";
	
	//--------------------------------------------------------------------------
    for(int i = 0; i < word_size;i++){
		MulMod(A_B0R0[i],A_B0R0[i],IN,p);
		MulMod(A_B0R1[i],A_B0R1[i],IN,p);
		MulMod(A_B0R2[i],A_B0R2[i],IN,p);
		MulMod(A_B0R3[i],A_B0R3[i],IN,p);
		MulMod(A_B0R4[i],A_B0R4[i],IN,p);
		MulMod(A_B0R5[i],A_B0R5[i],IN,p);
		MulMod(A_B0R6[i],A_B0R6[i],IN,p);
		MulMod(A_B0R7[i],A_B0R7[i],IN,p);
		MulMod(A_B0R8[i],A_B0R8[i],IN,p);
		MulMod(A_B0R9[i],A_B0R9[i],IN,p);
		MulMod(A_B0R10[i],A_B0R10[i],IN,p);
		MulMod(A_B0R11[i],A_B0R11[i],IN,p);
		MulMod(A_B0R12[i],A_B0R12[i],IN,p);
		MulMod(A_B0R13[i],A_B0R13[i],IN,p);
		MulMod(A_B0R14[i],A_B0R14[i],IN,p);
		MulMod(A_B0R15[i],A_B0R15[i],IN,p);
		//bank1 
		MulMod(A_B1R0[i],A_B1R0[i],IN,p);
		MulMod(A_B1R1[i],A_B1R1[i],IN,p);
		MulMod(A_B1R2[i],A_B1R2[i],IN,p);
		MulMod(A_B1R3[i],A_B1R3[i],IN,p);
		MulMod(A_B1R4[i],A_B1R4[i],IN,p);
		MulMod(A_B1R5[i],A_B1R5[i],IN,p);
		MulMod(A_B1R6[i],A_B1R6[i],IN,p);
		MulMod(A_B1R7[i],A_B1R7[i],IN,p);
		MulMod(A_B1R8[i],A_B1R8[i],IN,p);
		MulMod(A_B1R9[i],A_B1R9[i],IN,p);
		MulMod(A_B1R10[i],A_B1R10[i],IN,p);
		MulMod(A_B1R11[i],A_B1R11[i],IN,p);
		MulMod(A_B1R12[i],A_B1R12[i],IN,p);
		MulMod(A_B1R13[i],A_B1R13[i],IN,p);
		MulMod(A_B1R14[i],A_B1R14[i],IN,p);
		MulMod(A_B1R15[i],A_B1R15[i],IN,p);
	}
	
	std::ofstream BEFORE_DATA_RELOCATION("./INTT_R16_R2_Before_DATA_RELOCATION.txt");
	std::ofstream AFTER_DATA_RELOCATION("./INTT_R16_R2_After_DATA_RELOCATION.txt");
	
	//data relocation!!!!!!
	for(int i = 0;i < group;i++){
		for(int j = 0; j < radix; j++){
        	gray_i  = Gray(i,group);
        	BC_tmp  = j * group + gray_i;
        	BC_RR = BC_tmp;
			BC_IFFT_Reorder_R16_R2(BC_RR,BC_Reorder);
        	AGU_R16(BC_Reorder,bn_tmp,ma_tmp);
			BEFORE_DATA_RELOCATION <<"--------------------------------------------- \n";
			BEFORE_DATA_RELOCATION <<"i: "     << i <<"\n";
			BEFORE_DATA_RELOCATION <<"j: "     << j <<"\n";
			BEFORE_DATA_RELOCATION <<"gray_i: "<< gray_i <<"\n";
			BEFORE_DATA_RELOCATION <<"BC_tmp: "<< BC_tmp <<"\n";
			BEFORE_DATA_RELOCATION <<"BC_RR: " << BC_RR  <<"\n";
			BEFORE_DATA_RELOCATION <<"BC_Reorder: "<< BC_Reorder <<"\n";
			BEFORE_DATA_RELOCATION <<"bn_tmp: "<< bn_tmp <<"\n";
			BEFORE_DATA_RELOCATION <<"ma_tmp: "<< ma_tmp <<"\n";			
        	if(bn_tmp == 0){
				if(j < 2)bn0_bc_tmp = BC_tmp;
				BEFORE_DATA_RELOCATION <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]<<"\n";
				BEFORE_DATA_RELOCATION <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]<<"\n";
				BEFORE_DATA_RELOCATION <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]<<"\n";
				BEFORE_DATA_RELOCATION <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]<<"\n";
				BEFORE_DATA_RELOCATION <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]<<"\n";
				BEFORE_DATA_RELOCATION <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]<<"\n";
				BEFORE_DATA_RELOCATION <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]<<"\n";
				BEFORE_DATA_RELOCATION <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]<<"\n";
				BEFORE_DATA_RELOCATION <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]<<"\n";
				BEFORE_DATA_RELOCATION <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]<<"\n";
				BEFORE_DATA_RELOCATION <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";
				BEFORE_DATA_RELOCATION <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";
				BEFORE_DATA_RELOCATION <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";
				BEFORE_DATA_RELOCATION <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";
				BEFORE_DATA_RELOCATION <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";
				BEFORE_DATA_RELOCATION <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";	
				if(j <  2)bn0_ma_reg1 = ma_tmp;
				if((j >= 2)  && (j < 4))bn0_ma_reg2 = ma_tmp;
				if((j >= 4)  && (j < 6))bn0_ma_reg3 = ma_tmp;
				if((j >= 6)  && (j < 8))bn0_ma_reg4 = ma_tmp;
				if((j >= 8)  && (j < 10))bn0_ma_reg5 = ma_tmp;
				if((j >= 10) && (j < 12))bn0_ma_reg6 = ma_tmp;
				if((j >= 12) && (j < 14))bn0_ma_reg7 = ma_tmp;
				if((j >= 14) && (j < 16))bn0_ma_reg8 = ma_tmp;				
        	}else {
				if(j < 2)bn1_bc_tmp = BC_tmp;
				BEFORE_DATA_RELOCATION <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]<<"\n";
				BEFORE_DATA_RELOCATION <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]<<"\n";
				BEFORE_DATA_RELOCATION <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]<<"\n";
				BEFORE_DATA_RELOCATION <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]<<"\n";
				BEFORE_DATA_RELOCATION <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]<<"\n";
				BEFORE_DATA_RELOCATION <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]<<"\n";
				BEFORE_DATA_RELOCATION <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]<<"\n";
				BEFORE_DATA_RELOCATION <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]<<"\n";
				BEFORE_DATA_RELOCATION <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]<<"\n";
				BEFORE_DATA_RELOCATION <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]<<"\n";
				BEFORE_DATA_RELOCATION <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";
				BEFORE_DATA_RELOCATION <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";
				BEFORE_DATA_RELOCATION <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";
				BEFORE_DATA_RELOCATION <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";
				BEFORE_DATA_RELOCATION <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";
				BEFORE_DATA_RELOCATION <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";	
                if(j <  2)bn1_ma_reg1 = ma_tmp;					
                if((j >= 2)  && (j < 4))bn1_ma_reg2 = ma_tmp;
                if((j >= 4)  && (j < 6))bn1_ma_reg3 = ma_tmp;
                if((j >= 6)  && (j < 8))bn1_ma_reg4 = ma_tmp;
                if((j >= 8)  && (j < 10))bn1_ma_reg5 = ma_tmp;
                if((j >= 10) && (j < 12))bn1_ma_reg6 = ma_tmp;
                if((j >= 12) && (j < 14))bn1_ma_reg7 = ma_tmp;
                if((j >= 14) && (j < 16))bn1_ma_reg8 = ma_tmp;					
        	}
		}
		  if(bn1_bc_tmp > bn0_bc_tmp){
			 data_tmp_1  = A_B0R1[bn0_ma_reg1];
			 data_tmp_2  = A_B0R2[bn0_ma_reg1];
			 data_tmp_3  = A_B0R3[bn0_ma_reg1];
			 data_tmp_4  = A_B0R4[bn0_ma_reg1];
			 data_tmp_5  = A_B0R5[bn0_ma_reg1];
			 data_tmp_6  = A_B0R6[bn0_ma_reg1];
			 data_tmp_7  = A_B0R7[bn0_ma_reg1];
			 data_tmp_8  = A_B0R8[bn0_ma_reg1];
			 data_tmp_9  = A_B0R9[bn0_ma_reg1];
			 data_tmp_10 = A_B0R10[bn0_ma_reg1];
			 data_tmp_11 = A_B0R11[bn0_ma_reg1];
			 data_tmp_12 = A_B0R12[bn0_ma_reg1];
			 data_tmp_13 = A_B0R13[bn0_ma_reg1];
			 data_tmp_14 = A_B0R14[bn0_ma_reg1];
			 data_tmp_15 = A_B0R15[bn0_ma_reg1];
			 A_B0R1[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg1];
			 A_B0R2[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg2];
			 A_B0R3[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg2];
			 A_B0R4[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg3];
			 A_B0R5[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg3];
			 A_B0R6[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg1] = A_B0R0[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg1] = A_B1R0[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg1] = A_B0R0[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg1] = A_B1R0[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg1] = A_B1R0[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg1] = A_B0R0[bn0_ma_reg8];
			 A_B1R0[bn1_ma_reg1]  = data_tmp_1; 
			 A_B1R0[bn1_ma_reg2]  = data_tmp_2; 
			 A_B0R0[bn0_ma_reg2]  = data_tmp_3; 
			 A_B1R0[bn1_ma_reg3]  = data_tmp_4; 
			 A_B0R0[bn0_ma_reg3]  = data_tmp_5; 
			 A_B0R0[bn0_ma_reg4]  = data_tmp_6; 
			 A_B1R0[bn1_ma_reg4]  = data_tmp_7; 
			 A_B1R0[bn1_ma_reg5]  = data_tmp_8; 
			 A_B0R0[bn0_ma_reg5]  = data_tmp_9; 
			 A_B0R0[bn0_ma_reg6]  = data_tmp_10;
			 A_B1R0[bn1_ma_reg6]  = data_tmp_11;
			 A_B0R0[bn0_ma_reg7]  = data_tmp_12;
			 A_B1R0[bn1_ma_reg7]  = data_tmp_13;
			 A_B1R0[bn1_ma_reg8]  = data_tmp_14;
			 A_B0R0[bn0_ma_reg8]  = data_tmp_15;
			 data_tmp_1  = A_B1R2[bn1_ma_reg1];
			 data_tmp_2  = A_B1R3[bn1_ma_reg1];
			 data_tmp_3  = A_B1R4[bn1_ma_reg1];
			 data_tmp_4  = A_B1R5[bn1_ma_reg1];
			 data_tmp_5  = A_B1R6[bn1_ma_reg1];
			 data_tmp_6  = A_B1R7[bn1_ma_reg1];
			 data_tmp_7  = A_B1R8[bn1_ma_reg1];
			 data_tmp_8  = A_B1R9[bn1_ma_reg1];
			 data_tmp_9  = A_B1R10[bn1_ma_reg1];
			 data_tmp_10 = A_B1R11[bn1_ma_reg1];
			 data_tmp_11 = A_B1R12[bn1_ma_reg1];
			 data_tmp_12 = A_B1R13[bn1_ma_reg1];
			 data_tmp_13 = A_B1R14[bn1_ma_reg1];
			 data_tmp_14 = A_B1R15[bn1_ma_reg1];
             A_B1R2[bn1_ma_reg1]  =	 A_B1R1[bn1_ma_reg2];
             A_B1R3[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg2];
             A_B1R4[bn1_ma_reg1]  =  A_B1R1[bn1_ma_reg3];
             A_B1R5[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg3];
             A_B1R6[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg4];
             A_B1R7[bn1_ma_reg1]  =  A_B1R1[bn1_ma_reg4];
             A_B1R8[bn1_ma_reg1]  =  A_B1R1[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg1] =  A_B0R1[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg1] =  A_B1R1[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg1] =  A_B0R1[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg1] =  A_B1R1[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg1] =  A_B1R1[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg1] =  A_B0R1[bn0_ma_reg8];
			 A_B1R1[bn1_ma_reg2]  =  data_tmp_1; 
			 A_B0R1[bn0_ma_reg2]  =  data_tmp_2; 
			 A_B1R1[bn1_ma_reg3]  =  data_tmp_3; 
			 A_B0R1[bn0_ma_reg3]  =  data_tmp_4; 
			 A_B0R1[bn0_ma_reg4]  =  data_tmp_5; 
			 A_B1R1[bn1_ma_reg4]  =  data_tmp_6; 
			 A_B1R1[bn1_ma_reg5]  =  data_tmp_7; 
			 A_B0R1[bn0_ma_reg5]  =  data_tmp_8; 
			 A_B0R1[bn0_ma_reg6]  =  data_tmp_9; 
			 A_B1R1[bn1_ma_reg6]  =  data_tmp_10;
			 A_B0R1[bn0_ma_reg7]  =  data_tmp_11;
			 A_B1R1[bn1_ma_reg7]  =  data_tmp_12;
			 A_B1R1[bn1_ma_reg8]  =  data_tmp_13;
			 A_B0R1[bn0_ma_reg8]  =  data_tmp_14;
			/************************************************************/ 
			 data_tmp_1  =   A_B1R3[bn1_ma_reg2];
			 data_tmp_2  =   A_B1R4[bn1_ma_reg2];
			 data_tmp_3  =   A_B1R5[bn1_ma_reg2];
			 data_tmp_4  =   A_B1R6[bn1_ma_reg2];
			 data_tmp_5  =   A_B1R7[bn1_ma_reg2];
			 data_tmp_6  =   A_B1R8[bn1_ma_reg2];
			 data_tmp_7  =   A_B1R9[bn1_ma_reg2];
			 data_tmp_8  =   A_B1R10[bn1_ma_reg2];
			 data_tmp_9  =   A_B1R11[bn1_ma_reg2];
			 data_tmp_10 =   A_B1R12[bn1_ma_reg2];
			 data_tmp_11 =   A_B1R13[bn1_ma_reg2];
			 data_tmp_12 =   A_B1R14[bn1_ma_reg2];
			 data_tmp_13 =   A_B1R15[bn1_ma_reg2];
			 A_B1R3[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg2];
			 A_B1R4[bn1_ma_reg2]  =  A_B1R2[bn1_ma_reg3];
			 A_B1R5[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg3];
			 A_B1R6[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg4];
			 A_B1R7[bn1_ma_reg2]  =  A_B1R2[bn1_ma_reg4];
			 A_B1R8[bn1_ma_reg2]  =  A_B1R2[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg2] =  A_B0R2[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg2] =  A_B1R2[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg2] =  A_B0R2[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg2] =  A_B1R2[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg2] =  A_B1R2[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg2] =  A_B0R2[bn0_ma_reg8];
			 A_B0R2[bn0_ma_reg2]  =  data_tmp_1; 
			 A_B1R2[bn1_ma_reg3]  =  data_tmp_2; 
			 A_B0R2[bn0_ma_reg3]  =  data_tmp_3; 
			 A_B0R2[bn0_ma_reg4]  =  data_tmp_4; 
			 A_B1R2[bn1_ma_reg4]  =  data_tmp_5; 
			 A_B1R2[bn1_ma_reg5]  =  data_tmp_6; 
			 A_B0R2[bn0_ma_reg5]  =  data_tmp_7; 
			 A_B0R2[bn0_ma_reg6]  =  data_tmp_8; 
			 A_B1R2[bn1_ma_reg6]  =  data_tmp_9; 
			 A_B0R2[bn0_ma_reg7]  =  data_tmp_10;
			 A_B1R2[bn1_ma_reg7]  =  data_tmp_11;
			 A_B1R2[bn1_ma_reg8]  =  data_tmp_12;
			 A_B0R2[bn0_ma_reg8]  =  data_tmp_13;
			 //******************************************************
			 data_tmp_1  = A_B0R4[bn0_ma_reg2];
			 data_tmp_2  = A_B0R5[bn0_ma_reg2];
			 data_tmp_3  = A_B0R6[bn0_ma_reg2];
			 data_tmp_4  = A_B0R7[bn0_ma_reg2];
			 data_tmp_5  = A_B0R8[bn0_ma_reg2];
			 data_tmp_6  = A_B0R9[bn0_ma_reg2];
			 data_tmp_7  = A_B0R10[bn0_ma_reg2];
			 data_tmp_8  = A_B0R11[bn0_ma_reg2];
			 data_tmp_9  = A_B0R12[bn0_ma_reg2];
			 data_tmp_10 = A_B0R13[bn0_ma_reg2];
			 data_tmp_11 = A_B0R14[bn0_ma_reg2];
			 data_tmp_12 = A_B0R15[bn0_ma_reg2];
			 A_B0R4[bn0_ma_reg2]  =  A_B1R3[bn1_ma_reg3];
			 A_B0R5[bn0_ma_reg2]  =  A_B0R3[bn0_ma_reg3];
			 A_B0R6[bn0_ma_reg2]  =  A_B0R3[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg2]  =  A_B1R3[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg2]  =  A_B1R3[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg2]  =  A_B0R3[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg2] =  A_B0R3[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg2] =  A_B1R3[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg2] =  A_B0R3[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg2] =  A_B1R3[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg2] =  A_B1R3[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg2] =  A_B0R3[bn0_ma_reg8];
			 A_B1R3[bn1_ma_reg3]  =  data_tmp_1;  
			 A_B0R3[bn0_ma_reg3]  =  data_tmp_2;  
			 A_B0R3[bn0_ma_reg4]  =  data_tmp_3;  
			 A_B1R3[bn1_ma_reg4]  =  data_tmp_4;  
			 A_B1R3[bn1_ma_reg5]  =  data_tmp_5;  
			 A_B0R3[bn0_ma_reg5]  =  data_tmp_6;  
			 A_B0R3[bn0_ma_reg6]  =  data_tmp_7;  
			 A_B1R3[bn1_ma_reg6]  =  data_tmp_8;  
			 A_B0R3[bn0_ma_reg7]  =  data_tmp_9;  
			 A_B1R3[bn1_ma_reg7]  =  data_tmp_10; 
			 A_B1R3[bn1_ma_reg8]  =  data_tmp_11; 
			 A_B0R3[bn0_ma_reg8]  =  data_tmp_12; 
			 //----------------------------------------------------------------------
             data_tmp_1  = A_B1R5[bn1_ma_reg3];
			 data_tmp_2  = A_B1R6[bn1_ma_reg3];
			 data_tmp_3  = A_B1R7[bn1_ma_reg3];
			 data_tmp_4  = A_B1R8[bn1_ma_reg3];
			 data_tmp_5  = A_B1R9[bn1_ma_reg3];
			 data_tmp_6  = A_B1R10[bn1_ma_reg3];
			 data_tmp_7  = A_B1R11[bn1_ma_reg3];
			 data_tmp_8  = A_B1R12[bn1_ma_reg3];
			 data_tmp_9  = A_B1R13[bn1_ma_reg3];
			 data_tmp_10 = A_B1R14[bn1_ma_reg3];
			 data_tmp_11 = A_B1R15[bn1_ma_reg3];
			 A_B1R5[bn1_ma_reg3]  =  A_B0R4[bn0_ma_reg3];
			 A_B1R6[bn1_ma_reg3]  =  A_B0R4[bn0_ma_reg4];
			 A_B1R7[bn1_ma_reg3]  =  A_B1R4[bn1_ma_reg4];
			 A_B1R8[bn1_ma_reg3]  =  A_B1R4[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg3]  =  A_B0R4[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg3] =  A_B0R4[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg3] =  A_B1R4[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg3] =  A_B0R4[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg3] =  A_B1R4[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg3] =  A_B1R4[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg3] =  A_B0R4[bn0_ma_reg8];
			 A_B0R4[bn0_ma_reg3]  =  data_tmp_1 ;
			 A_B0R4[bn0_ma_reg4]  =  data_tmp_2 ;
			 A_B1R4[bn1_ma_reg4]  =  data_tmp_3 ;
			 A_B1R4[bn1_ma_reg5]  =  data_tmp_4 ;
			 A_B0R4[bn0_ma_reg5]  =  data_tmp_5 ;
			 A_B0R4[bn0_ma_reg6]  =  data_tmp_6 ;
			 A_B1R4[bn1_ma_reg6]  =  data_tmp_7 ;
			 A_B0R4[bn0_ma_reg7]  =  data_tmp_8 ;
			 A_B1R4[bn1_ma_reg7]  =  data_tmp_9 ;
			 A_B1R4[bn1_ma_reg8]  =  data_tmp_10;
			 A_B0R4[bn0_ma_reg8]  =  data_tmp_11;
			 //************************************************************************
			 data_tmp_1  = A_B0R6[bn0_ma_reg3];
			 data_tmp_2  = A_B0R7[bn0_ma_reg3];
			 data_tmp_3  = A_B0R8[bn0_ma_reg3];
			 data_tmp_4  = A_B0R9[bn0_ma_reg3];
			 data_tmp_5  = A_B0R10[bn0_ma_reg3];
			 data_tmp_6  = A_B0R11[bn0_ma_reg3];
			 data_tmp_7  = A_B0R12[bn0_ma_reg3];
			 data_tmp_8  = A_B0R13[bn0_ma_reg3];
			 data_tmp_9  = A_B0R14[bn0_ma_reg3];
			 data_tmp_10 = A_B0R15[bn0_ma_reg3];
			 A_B0R6[bn0_ma_reg3]  =  A_B0R5[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg3]  =  A_B1R5[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg3]  =  A_B1R5[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg3]  =  A_B0R5[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg3] =  A_B0R5[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg3] =  A_B1R5[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg3] =  A_B0R5[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg3] =  A_B1R5[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg3] =  A_B1R5[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg3] =  A_B0R5[bn0_ma_reg8];
			 A_B0R5[bn0_ma_reg4]  =  data_tmp_1; 
			 A_B1R5[bn1_ma_reg4]  =  data_tmp_2; 
			 A_B1R5[bn1_ma_reg5]  =  data_tmp_3; 
			 A_B0R5[bn0_ma_reg5]  =  data_tmp_4; 
			 A_B0R5[bn0_ma_reg6]  =  data_tmp_5; 
			 A_B1R5[bn1_ma_reg6]  =  data_tmp_6; 
			 A_B0R5[bn0_ma_reg7]  =  data_tmp_7; 
			 A_B1R5[bn1_ma_reg7]  =  data_tmp_8; 
			 A_B1R5[bn1_ma_reg8]  =  data_tmp_9; 
			 A_B0R5[bn0_ma_reg8]  =  data_tmp_10;
			 //-----------------------------------------------------------------------
			 data_tmp_1  = A_B0R7[bn0_ma_reg4];
			 data_tmp_2  = A_B0R8[bn0_ma_reg4];
			 data_tmp_3  = A_B0R9[bn0_ma_reg4];
			 data_tmp_4  = A_B0R10[bn0_ma_reg4];
			 data_tmp_5  = A_B0R11[bn0_ma_reg4];
			 data_tmp_6  = A_B0R12[bn0_ma_reg4];
			 data_tmp_7  = A_B0R13[bn0_ma_reg4];
			 data_tmp_8  = A_B0R14[bn0_ma_reg4];
			 data_tmp_9  = A_B0R15[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg4]  =  A_B1R6[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg4]  =  A_B1R6[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg4]  =  A_B0R6[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg4] =  A_B0R6[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg4] =  A_B1R6[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg4] =  A_B0R6[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg4] =  A_B1R6[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg4] =  A_B1R6[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg4] =  A_B0R6[bn0_ma_reg8];
			 A_B1R6[bn1_ma_reg4]  =  data_tmp_1;
			 A_B1R6[bn1_ma_reg5]  =  data_tmp_2;
			 A_B0R6[bn0_ma_reg5]  =  data_tmp_3;
			 A_B0R6[bn0_ma_reg6]  =  data_tmp_4;
			 A_B1R6[bn1_ma_reg6]  =  data_tmp_5;
			 A_B0R6[bn0_ma_reg7]  =  data_tmp_6;
			 A_B1R6[bn1_ma_reg7]  =  data_tmp_7;
			 A_B1R6[bn1_ma_reg8]  =  data_tmp_8;
			 A_B0R6[bn0_ma_reg8]  =  data_tmp_9;
			 //----------------------------------------------------------------------
			 data_tmp_1  = A_B1R8[bn1_ma_reg4];
			 data_tmp_2  = A_B1R9[bn1_ma_reg4];
			 data_tmp_3  = A_B1R10[bn1_ma_reg4];
			 data_tmp_4  = A_B1R11[bn1_ma_reg4];
			 data_tmp_5  = A_B1R12[bn1_ma_reg4];
			 data_tmp_6  = A_B1R13[bn1_ma_reg4];
			 data_tmp_7  = A_B1R14[bn1_ma_reg4];
			 data_tmp_8  = A_B1R15[bn1_ma_reg4];
			 A_B1R8[bn1_ma_reg4]  =  A_B1R7[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg4]  =  A_B0R7[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg4] =  A_B0R7[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg4] =  A_B1R7[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg4] =  A_B0R7[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg4] =  A_B1R7[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg4] =  A_B1R7[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg4] =  A_B0R7[bn0_ma_reg8];
			 A_B1R7[bn1_ma_reg5]  =  data_tmp_1;
			 A_B0R7[bn0_ma_reg5]  =  data_tmp_2;
			 A_B0R7[bn0_ma_reg6]  =  data_tmp_3;
			 A_B1R7[bn1_ma_reg6]  =  data_tmp_4;
			 A_B0R7[bn0_ma_reg7]  =  data_tmp_5;
			 A_B1R7[bn1_ma_reg7]  =  data_tmp_6;
			 A_B1R7[bn1_ma_reg8]  =  data_tmp_7;
			 A_B0R7[bn0_ma_reg8]  =  data_tmp_8;
			 //----------------------------------------------------------------------
			 data_tmp_1 = A_B1R9[bn1_ma_reg5];
			 data_tmp_2 = A_B1R10[bn1_ma_reg5];
			 data_tmp_3 = A_B1R11[bn1_ma_reg5];
			 data_tmp_4 = A_B1R12[bn1_ma_reg5];
			 data_tmp_5 = A_B1R13[bn1_ma_reg5];
			 data_tmp_6 = A_B1R14[bn1_ma_reg5];
			 data_tmp_7 = A_B1R15[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg5]  = A_B0R8[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg5] = A_B0R8[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg5] = A_B1R8[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg5] = A_B0R8[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg5] = A_B1R8[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg5] = A_B1R8[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg5] = A_B0R8[bn0_ma_reg8];
			 A_B0R8[bn0_ma_reg5]  = data_tmp_1;
			 A_B0R8[bn0_ma_reg6]  = data_tmp_2;
			 A_B1R8[bn1_ma_reg6]  = data_tmp_3;
			 A_B0R8[bn0_ma_reg7]  = data_tmp_4;
			 A_B1R8[bn1_ma_reg7]  = data_tmp_5;
			 A_B1R8[bn1_ma_reg8]  = data_tmp_6;
			 A_B0R8[bn0_ma_reg8]  = data_tmp_7;
			 //---------------------------------------------------------------------
			 data_tmp_1  = A_B0R10[bn0_ma_reg5];
			 data_tmp_2  = A_B0R11[bn0_ma_reg5];
			 data_tmp_3  = A_B0R12[bn0_ma_reg5];
			 data_tmp_4  = A_B0R13[bn0_ma_reg5];
			 data_tmp_5  = A_B0R14[bn0_ma_reg5];
			 data_tmp_6  = A_B0R15[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg5] = A_B0R9[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg5] = A_B1R9[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg5] = A_B0R9[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg5] = A_B1R9[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg5] = A_B1R9[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg5] = A_B0R9[bn0_ma_reg8];
			 A_B0R9[bn0_ma_reg6]  = data_tmp_1;
			 A_B1R9[bn1_ma_reg6]  = data_tmp_2;
			 A_B0R9[bn0_ma_reg7]  = data_tmp_3;
			 A_B1R9[bn1_ma_reg7]  = data_tmp_4;
			 A_B1R9[bn1_ma_reg8]  = data_tmp_5;
			 A_B0R9[bn0_ma_reg8]  = data_tmp_6;
			 //--------------------------------------------------------------------
			 data_tmp_1  = A_B0R11[bn0_ma_reg6];
			 data_tmp_2  = A_B0R12[bn0_ma_reg6];
			 data_tmp_3  = A_B0R13[bn0_ma_reg6];
			 data_tmp_4  = A_B0R14[bn0_ma_reg6];
			 data_tmp_5  = A_B0R15[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg6] = A_B1R10[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg6] = A_B0R10[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg6] = A_B1R10[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg6] = A_B1R10[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg6] = A_B0R10[bn0_ma_reg8];
			 A_B1R10[bn1_ma_reg6] = data_tmp_1;
			 A_B0R10[bn0_ma_reg7] = data_tmp_2;
			 A_B1R10[bn1_ma_reg7] = data_tmp_3;
			 A_B1R10[bn1_ma_reg8] = data_tmp_4;
			 A_B0R10[bn0_ma_reg8] = data_tmp_5;
			 //--------------------------------------------------------------------
			 data_tmp_1  = A_B1R12[bn1_ma_reg6];
			 data_tmp_2  = A_B1R13[bn1_ma_reg6];
			 data_tmp_3  = A_B1R14[bn1_ma_reg6];
			 data_tmp_4  = A_B1R15[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg6] = A_B0R11[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg6] = A_B1R11[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg6] = A_B1R11[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg6] = A_B0R11[bn0_ma_reg8];
			 A_B0R11[bn0_ma_reg7] = data_tmp_1;
			 A_B1R11[bn1_ma_reg7] = data_tmp_2;
			 A_B1R11[bn1_ma_reg8] = data_tmp_3;
			 A_B0R11[bn0_ma_reg8] = data_tmp_4;
			 //--------------------------------------------------------------------
			 data_tmp_1 = A_B0R13[bn0_ma_reg7];
			 data_tmp_2 = A_B0R14[bn0_ma_reg7];
			 data_tmp_3 = A_B0R15[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg7] = A_B1R12[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg7] = A_B1R12[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg7] = A_B0R12[bn0_ma_reg8];
			 A_B1R12[bn1_ma_reg7] = data_tmp_1;
			 A_B1R12[bn1_ma_reg8] = data_tmp_2;
			 A_B0R12[bn0_ma_reg8] = data_tmp_3;
			 //--------------------------------------------------------------------
			 data_tmp_1 = A_B1R14[bn1_ma_reg7];
			 data_tmp_2 = A_B1R15[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg7] = A_B1R13[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg7] = A_B0R13[bn0_ma_reg8];
			 A_B1R13[bn1_ma_reg8] = data_tmp_1;
			 A_B0R13[bn0_ma_reg8] = data_tmp_2;
			 //--------------------------------------------------------------------
			 data_tmp_1 = A_B1R15[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg8] = A_B0R14[bn0_ma_reg8];
			 A_B0R14[bn0_ma_reg8] = data_tmp_1; 
		  }else{	 
			 data_tmp_1  = A_B1R1[bn1_ma_reg1];
			 data_tmp_2  = A_B1R2[bn1_ma_reg1];
			 data_tmp_3  = A_B1R3[bn1_ma_reg1];
			 data_tmp_4  = A_B1R4[bn1_ma_reg1];
			 data_tmp_5  = A_B1R5[bn1_ma_reg1];
			 data_tmp_6  = A_B1R6[bn1_ma_reg1];
			 data_tmp_7  = A_B1R7[bn1_ma_reg1];
			 data_tmp_8  = A_B1R8[bn1_ma_reg1];
			 data_tmp_9  = A_B1R9[bn1_ma_reg1];
			 data_tmp_10 = A_B1R10[bn1_ma_reg1];
			 data_tmp_11 = A_B1R11[bn1_ma_reg1];
			 data_tmp_12 = A_B1R12[bn1_ma_reg1];
			 data_tmp_13 = A_B1R13[bn1_ma_reg1];
			 data_tmp_14 = A_B1R14[bn1_ma_reg1];
			 data_tmp_15 = A_B1R15[bn1_ma_reg1];
             A_B1R1[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg1];
             A_B1R2[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg2];
             A_B1R3[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg2];
             A_B1R4[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg3];
             A_B1R5[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg3];
             A_B1R6[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg4];
             A_B1R7[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg4];
             A_B1R8[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg5];
             A_B1R9[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg5];
             A_B1R10[bn1_ma_reg1] =  A_B1R0[bn1_ma_reg6];
             A_B1R11[bn1_ma_reg1] =  A_B0R0[bn0_ma_reg6];
             A_B1R12[bn1_ma_reg1] =  A_B1R0[bn1_ma_reg7];
             A_B1R13[bn1_ma_reg1] =  A_B0R0[bn0_ma_reg7];
             A_B1R14[bn1_ma_reg1] =  A_B0R0[bn0_ma_reg8];
             A_B1R15[bn1_ma_reg1] =  A_B1R0[bn1_ma_reg8];
             A_B0R0[bn0_ma_reg1]  =  data_tmp_1; 
             A_B0R0[bn0_ma_reg2]  =  data_tmp_2; 
             A_B1R0[bn1_ma_reg2]  =  data_tmp_3; 
             A_B0R0[bn0_ma_reg3]  =  data_tmp_4; 
             A_B1R0[bn1_ma_reg3]  =  data_tmp_5; 
             A_B1R0[bn1_ma_reg4]  =  data_tmp_6; 
             A_B0R0[bn0_ma_reg4]  =  data_tmp_7; 
             A_B0R0[bn0_ma_reg5]  =  data_tmp_8; 
             A_B1R0[bn1_ma_reg5]  =  data_tmp_9; 
             A_B1R0[bn1_ma_reg6]  =  data_tmp_10;
             A_B0R0[bn0_ma_reg6]  =  data_tmp_11;
             A_B1R0[bn1_ma_reg7]  =  data_tmp_12;
             A_B0R0[bn0_ma_reg7]  =  data_tmp_13;
             A_B0R0[bn0_ma_reg8]  =  data_tmp_14;
             A_B1R0[bn1_ma_reg8]  =  data_tmp_15;
             //-----------------------------------------------------------
             data_tmp_1  = A_B0R2[bn0_ma_reg1];
			 data_tmp_2  = A_B0R3[bn0_ma_reg1];
			 data_tmp_3  = A_B0R4[bn0_ma_reg1];
			 data_tmp_4  = A_B0R5[bn0_ma_reg1];
			 data_tmp_5  = A_B0R6[bn0_ma_reg1];
			 data_tmp_6  = A_B0R7[bn0_ma_reg1];
			 data_tmp_7  = A_B0R8[bn0_ma_reg1];
			 data_tmp_8  = A_B0R9[bn0_ma_reg1];
			 data_tmp_9  = A_B0R10[bn0_ma_reg1];
			 data_tmp_10 = A_B0R11[bn0_ma_reg1];
			 data_tmp_11 = A_B0R12[bn0_ma_reg1];
			 data_tmp_12 = A_B0R13[bn0_ma_reg1];
			 data_tmp_13 = A_B0R14[bn0_ma_reg1];
			 data_tmp_14 = A_B0R15[bn0_ma_reg1];
			 A_B0R2[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg2];
			 A_B0R3[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg2];
			 A_B0R4[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg3];
			 A_B0R5[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg3];
			 A_B0R6[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg4];
			 A_B0R7[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg1] =  A_B1R1[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg1] =  A_B0R1[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg1] =  A_B1R1[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg1] =  A_B0R1[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg1] =  A_B0R1[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg1] =  A_B1R1[bn1_ma_reg8];
             A_B0R1[bn0_ma_reg2]  =  data_tmp_1; 
             A_B1R1[bn1_ma_reg2]  =  data_tmp_2; 
             A_B0R1[bn0_ma_reg3]  =  data_tmp_3; 
             A_B1R1[bn1_ma_reg3]  =  data_tmp_4; 
             A_B1R1[bn1_ma_reg4]  =  data_tmp_5; 
             A_B0R1[bn0_ma_reg4]  =  data_tmp_6; 
             A_B0R1[bn0_ma_reg5]  =  data_tmp_7; 
             A_B1R1[bn1_ma_reg5]  =  data_tmp_8; 
             A_B1R1[bn1_ma_reg6]  =  data_tmp_9; 
             A_B0R1[bn0_ma_reg6]  =  data_tmp_10;
             A_B1R1[bn1_ma_reg7]  =  data_tmp_11;
             A_B0R1[bn0_ma_reg7]  =  data_tmp_12;
             A_B0R1[bn0_ma_reg8]  =  data_tmp_13;
             A_B1R1[bn1_ma_reg8]  =  data_tmp_14;
             //------------------------------------------------------------
             data_tmp_1  =  A_B0R3[bn0_ma_reg2];
			 data_tmp_2  =  A_B0R4[bn0_ma_reg2];
			 data_tmp_3  =  A_B0R5[bn0_ma_reg2];
			 data_tmp_4  =  A_B0R6[bn0_ma_reg2];
			 data_tmp_5  =  A_B0R7[bn0_ma_reg2];
			 data_tmp_6  =  A_B0R8[bn0_ma_reg2];
			 data_tmp_7  =  A_B0R9[bn0_ma_reg2];
			 data_tmp_8  =  A_B0R10[bn0_ma_reg2];
			 data_tmp_9  =  A_B0R11[bn0_ma_reg2];
			 data_tmp_10 =  A_B0R12[bn0_ma_reg2];
			 data_tmp_11 =  A_B0R13[bn0_ma_reg2];
			 data_tmp_12 =  A_B0R14[bn0_ma_reg2];
			 data_tmp_13 =  A_B0R15[bn0_ma_reg2];
			 A_B0R3[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg2];
			 A_B0R4[bn0_ma_reg2]  =  A_B0R2[bn0_ma_reg3];
			 A_B0R5[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg3];
			 A_B0R6[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg4];
			 A_B0R7[bn0_ma_reg2]  =  A_B0R2[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg2]  =  A_B0R2[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg2] =  A_B1R2[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg2] =  A_B0R2[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg2] =  A_B1R2[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg2] =  A_B0R2[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg2] =  A_B0R2[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg2] =  A_B1R2[bn1_ma_reg8];
			 A_B1R2[bn1_ma_reg2]  =  data_tmp_1; 
			 A_B0R2[bn0_ma_reg3]  =  data_tmp_2; 
			 A_B1R2[bn1_ma_reg3]  =  data_tmp_3; 
			 A_B1R2[bn1_ma_reg4]  =  data_tmp_4; 
			 A_B0R2[bn0_ma_reg4]  =  data_tmp_5; 
			 A_B0R2[bn0_ma_reg5]  =  data_tmp_6; 
			 A_B1R2[bn1_ma_reg5]  =  data_tmp_7; 
			 A_B1R2[bn1_ma_reg6]  =  data_tmp_8; 
			 A_B0R2[bn0_ma_reg6]  =  data_tmp_9; 
			 A_B1R2[bn1_ma_reg7]  =  data_tmp_10;
			 A_B0R2[bn0_ma_reg7]  =  data_tmp_11;
			 A_B0R2[bn0_ma_reg8]  =  data_tmp_12;
			 A_B1R2[bn1_ma_reg8]  =  data_tmp_13;
			 //-----------------------------------------------------------                        
			 data_tmp_1  = A_B1R4[bn1_ma_reg2];
			 data_tmp_2  = A_B1R5[bn1_ma_reg2];
			 data_tmp_3  = A_B1R6[bn1_ma_reg2];
			 data_tmp_4  = A_B1R7[bn1_ma_reg2];
			 data_tmp_5  = A_B1R8[bn1_ma_reg2];
			 data_tmp_6  = A_B1R9[bn1_ma_reg2];
			 data_tmp_7  = A_B1R10[bn1_ma_reg2];
			 data_tmp_8  = A_B1R11[bn1_ma_reg2];
			 data_tmp_9  = A_B1R12[bn1_ma_reg2];
			 data_tmp_10 = A_B1R13[bn1_ma_reg2];
			 data_tmp_11 = A_B1R14[bn1_ma_reg2];
			 data_tmp_12 = A_B1R15[bn1_ma_reg2];
			 A_B1R4[bn1_ma_reg2]  = A_B0R3[bn0_ma_reg3];
			 A_B1R5[bn1_ma_reg2]  = A_B1R3[bn1_ma_reg3];
			 A_B1R6[bn1_ma_reg2]  = A_B1R3[bn1_ma_reg4];
			 A_B1R7[bn1_ma_reg2]  = A_B0R3[bn0_ma_reg4];
			 A_B1R8[bn1_ma_reg2]  = A_B0R3[bn0_ma_reg5];
			 A_B1R9[bn1_ma_reg2]  = A_B1R3[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg2] = A_B1R3[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg2] = A_B0R3[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg2] = A_B1R3[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg2] = A_B0R3[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg2] = A_B0R3[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg2] = A_B1R3[bn1_ma_reg8];
			 A_B0R3[bn0_ma_reg3]  = data_tmp_1; 
			 A_B1R3[bn1_ma_reg3]  = data_tmp_2; 
			 A_B1R3[bn1_ma_reg4]  = data_tmp_3; 
			 A_B0R3[bn0_ma_reg4]  = data_tmp_4; 
			 A_B0R3[bn0_ma_reg5]  = data_tmp_5; 
			 A_B1R3[bn1_ma_reg5]  = data_tmp_6; 
			 A_B1R3[bn1_ma_reg6]  = data_tmp_7; 
			 A_B0R3[bn0_ma_reg6]  = data_tmp_8; 
			 A_B1R3[bn1_ma_reg7]  = data_tmp_9; 
			 A_B0R3[bn0_ma_reg7]  = data_tmp_10;
			 A_B0R3[bn0_ma_reg8]  = data_tmp_11;
			 A_B1R3[bn1_ma_reg8]  = data_tmp_12;
			 //------------------------------------------------------------
			 data_tmp_1  =  A_B0R5[bn0_ma_reg3];
			 data_tmp_2  =  A_B0R6[bn0_ma_reg3];
			 data_tmp_3  =  A_B0R7[bn0_ma_reg3];
			 data_tmp_4  =  A_B0R8[bn0_ma_reg3];
			 data_tmp_5  =  A_B0R9[bn0_ma_reg3];
			 data_tmp_6  =  A_B0R10[bn0_ma_reg3];
			 data_tmp_7  =  A_B0R11[bn0_ma_reg3];
			 data_tmp_8  =  A_B0R12[bn0_ma_reg3];
			 data_tmp_9  =  A_B0R13[bn0_ma_reg3];
			 data_tmp_10 =  A_B0R14[bn0_ma_reg3];
			 data_tmp_11 =  A_B0R15[bn0_ma_reg3];
			 A_B0R5[bn0_ma_reg3]  =  A_B1R4[bn1_ma_reg3];
			 A_B0R6[bn0_ma_reg3]  =  A_B1R4[bn1_ma_reg4];
			 A_B0R7[bn0_ma_reg3]  =  A_B0R4[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg3]  =  A_B0R4[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg3]  =  A_B1R4[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg3] =  A_B1R4[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg3] =  A_B0R4[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg3] =  A_B1R4[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg3] =  A_B0R4[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg3] =  A_B0R4[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg3] =  A_B1R4[bn1_ma_reg8];
			 A_B1R4[bn1_ma_reg3]  =  data_tmp_1; 
			 A_B1R4[bn1_ma_reg4]  =  data_tmp_2; 
			 A_B0R4[bn0_ma_reg4]  =  data_tmp_3; 
			 A_B0R4[bn0_ma_reg5]  =  data_tmp_4; 
			 A_B1R4[bn1_ma_reg5]  =  data_tmp_5; 
			 A_B1R4[bn1_ma_reg6]  =  data_tmp_6; 
			 A_B0R4[bn0_ma_reg6]  =  data_tmp_7; 
			 A_B1R4[bn1_ma_reg7]  =  data_tmp_8; 
			 A_B0R4[bn0_ma_reg7]  =  data_tmp_9; 
			 A_B0R4[bn0_ma_reg8]  =  data_tmp_10;
			 A_B1R4[bn1_ma_reg8]  =  data_tmp_11;
			 //-------------------------------------------------------------
			 data_tmp_1  =  A_B1R6[bn1_ma_reg3];
			 data_tmp_2  =  A_B1R7[bn1_ma_reg3];
			 data_tmp_3  =  A_B1R8[bn1_ma_reg3];
			 data_tmp_4  =  A_B1R9[bn1_ma_reg3];
			 data_tmp_5  =  A_B1R10[bn1_ma_reg3];
			 data_tmp_6  =  A_B1R11[bn1_ma_reg3];
			 data_tmp_7  =  A_B1R12[bn1_ma_reg3];
			 data_tmp_8  =  A_B1R13[bn1_ma_reg3];
			 data_tmp_9  =  A_B1R14[bn1_ma_reg3];
			 data_tmp_10 =  A_B1R15[bn1_ma_reg3];
			 A_B1R6[bn1_ma_reg3]  =  A_B1R5[bn1_ma_reg4];
			 A_B1R7[bn1_ma_reg3]  =  A_B0R5[bn0_ma_reg4];
			 A_B1R8[bn1_ma_reg3]  =  A_B0R5[bn0_ma_reg5];
			 A_B1R9[bn1_ma_reg3]  =  A_B1R5[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg3] =  A_B1R5[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg3] =  A_B0R5[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg3] =  A_B1R5[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg3] =  A_B0R5[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg3] =  A_B0R5[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg3] =  A_B1R5[bn1_ma_reg8];
			 A_B1R5[bn1_ma_reg4]  =  data_tmp_1; 
			 A_B0R5[bn0_ma_reg4]  =  data_tmp_2; 
			 A_B0R5[bn0_ma_reg5]  =  data_tmp_3; 
			 A_B1R5[bn1_ma_reg5]  =  data_tmp_4; 
			 A_B1R5[bn1_ma_reg6]  =  data_tmp_5; 
			 A_B0R5[bn0_ma_reg6]  =  data_tmp_6; 
			 A_B1R5[bn1_ma_reg7]  =  data_tmp_7; 
			 A_B0R5[bn0_ma_reg7]  =  data_tmp_8; 
			 A_B0R5[bn0_ma_reg8]  =  data_tmp_9; 
			 A_B1R5[bn1_ma_reg8]  =  data_tmp_10;
			 //----------------------------------------------------------------
			 data_tmp_1  =  A_B1R7[bn1_ma_reg4];
			 data_tmp_2  =  A_B1R8[bn1_ma_reg4];
			 data_tmp_3  =  A_B1R9[bn1_ma_reg4];
			 data_tmp_4  =  A_B1R10[bn1_ma_reg4];
			 data_tmp_5  =  A_B1R11[bn1_ma_reg4];
			 data_tmp_6  =  A_B1R12[bn1_ma_reg4];
			 data_tmp_7  =  A_B1R13[bn1_ma_reg4];
			 data_tmp_8  =  A_B1R14[bn1_ma_reg4];
			 data_tmp_9  =  A_B1R15[bn1_ma_reg4];
			 A_B1R7[bn1_ma_reg4]  =  A_B0R6[bn0_ma_reg4];
			 A_B1R8[bn1_ma_reg4]  =  A_B0R6[bn0_ma_reg5];
			 A_B1R9[bn1_ma_reg4]  =  A_B1R6[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg4] =  A_B1R6[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg4] =  A_B0R6[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg4] =  A_B1R6[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg4] =  A_B0R6[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg4] =  A_B0R6[bn0_ma_reg8];
             A_B1R15[bn1_ma_reg4] =  A_B1R6[bn1_ma_reg8];
             A_B0R6[bn0_ma_reg4]  =  data_tmp_1;
             A_B0R6[bn0_ma_reg5]  =  data_tmp_2;
             A_B1R6[bn1_ma_reg5]  =  data_tmp_3;
             A_B1R6[bn1_ma_reg6]  =  data_tmp_4;
             A_B0R6[bn0_ma_reg6]  =  data_tmp_5;
             A_B1R6[bn1_ma_reg7]  =  data_tmp_6;
             A_B0R6[bn0_ma_reg7]  =  data_tmp_7;
             A_B0R6[bn0_ma_reg8]  =  data_tmp_8;
             A_B1R6[bn1_ma_reg8]  =  data_tmp_9;
             //------------------------------------------------------------------
			 data_tmp_1  = A_B0R8[bn0_ma_reg4];
			 data_tmp_2  = A_B0R9[bn0_ma_reg4];
			 data_tmp_3  = A_B0R10[bn0_ma_reg4];
			 data_tmp_4  = A_B0R11[bn0_ma_reg4];
			 data_tmp_5  = A_B0R12[bn0_ma_reg4];
			 data_tmp_6  = A_B0R13[bn0_ma_reg4];
			 data_tmp_7  = A_B0R14[bn0_ma_reg4];
			 data_tmp_8  = A_B0R15[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg4]  =  A_B0R7[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg4]  =  A_B1R7[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg4] =  A_B1R7[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg4] =  A_B0R7[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg4] =  A_B1R7[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg4] =  A_B0R7[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg4] =  A_B0R7[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg4] =  A_B1R7[bn1_ma_reg8];
			 A_B0R7[bn0_ma_reg5]  =  data_tmp_1;
			 A_B1R7[bn1_ma_reg5]  =  data_tmp_2;
			 A_B1R7[bn1_ma_reg6]  =  data_tmp_3;
			 A_B0R7[bn0_ma_reg6]  =  data_tmp_4;
			 A_B1R7[bn1_ma_reg7]  =  data_tmp_5;
			 A_B0R7[bn0_ma_reg7]  =  data_tmp_6;
			 A_B0R7[bn0_ma_reg8]  =  data_tmp_7;
			 A_B1R7[bn1_ma_reg8]  =  data_tmp_8;
			 //----------------------------------------------------------------
			 data_tmp_1  = A_B0R9[bn0_ma_reg5];
			 data_tmp_2  = A_B0R10[bn0_ma_reg5];
			 data_tmp_3  = A_B0R11[bn0_ma_reg5];
			 data_tmp_4  = A_B0R12[bn0_ma_reg5];
			 data_tmp_5  = A_B0R13[bn0_ma_reg5];
			 data_tmp_6  = A_B0R14[bn0_ma_reg5];
			 data_tmp_7  = A_B0R15[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg5]  =  A_B1R8[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg8];
			 A_B1R8[bn1_ma_reg5]  =  data_tmp_1;
			 A_B1R8[bn1_ma_reg6]  =  data_tmp_2;
			 A_B0R8[bn0_ma_reg6]  =  data_tmp_3;
			 A_B1R8[bn1_ma_reg7]  =  data_tmp_4;
			 A_B0R8[bn0_ma_reg7]  =  data_tmp_5;
			 A_B0R8[bn0_ma_reg8]  =  data_tmp_6;
			 A_B1R8[bn1_ma_reg8]  =  data_tmp_7;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B1R10[bn1_ma_reg5];
			 data_tmp_2 = A_B1R11[bn1_ma_reg5];
			 data_tmp_3 = A_B1R12[bn1_ma_reg5];
			 data_tmp_4 = A_B1R13[bn1_ma_reg5];
			 data_tmp_5 = A_B1R14[bn1_ma_reg5];
			 data_tmp_6 = A_B1R15[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg5] = A_B1R9[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg5] = A_B0R9[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg5] = A_B1R9[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg5] = A_B0R9[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg5] = A_B0R9[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg5] = A_B1R9[bn1_ma_reg8];
			 A_B1R9[bn1_ma_reg6]  = data_tmp_1;
			 A_B0R9[bn0_ma_reg6]  = data_tmp_2;
			 A_B1R9[bn1_ma_reg7]  = data_tmp_3;
			 A_B0R9[bn0_ma_reg7]  = data_tmp_4;
			 A_B0R9[bn0_ma_reg8]  = data_tmp_5;
			 A_B1R9[bn1_ma_reg8]  = data_tmp_6;
			 //----------------------------------------------------------------
			 data_tmp_1 =  A_B1R11[bn1_ma_reg6];
			 data_tmp_2 =  A_B1R12[bn1_ma_reg6];
			 data_tmp_3 =  A_B1R13[bn1_ma_reg6];
			 data_tmp_4 =  A_B1R14[bn1_ma_reg6];
			 data_tmp_5 =  A_B1R15[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg6]  = A_B0R10[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg6]  = A_B1R10[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg6]  = A_B0R10[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg6]  = A_B0R10[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg6]  = A_B1R10[bn1_ma_reg8];
			 A_B0R10[bn0_ma_reg6]  = data_tmp_1;
			 A_B1R10[bn1_ma_reg7]  = data_tmp_2;
			 A_B0R10[bn0_ma_reg7]  = data_tmp_3;
			 A_B0R10[bn0_ma_reg8]  = data_tmp_4;
			 A_B1R10[bn1_ma_reg8]  = data_tmp_5;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B0R12[bn0_ma_reg6];
			 data_tmp_2 = A_B0R13[bn0_ma_reg6];
			 data_tmp_3 = A_B0R14[bn0_ma_reg6];
			 data_tmp_4 = A_B0R15[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg6] = A_B1R11[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg6] = A_B0R11[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg6] = A_B0R11[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg6] = A_B1R11[bn1_ma_reg8];
			 A_B1R11[bn1_ma_reg7] = data_tmp_1;
			 A_B0R11[bn0_ma_reg7] = data_tmp_2;
			 A_B0R11[bn0_ma_reg8] = data_tmp_3;
			 A_B1R11[bn1_ma_reg8] = data_tmp_4;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B1R13[bn1_ma_reg7];
			 data_tmp_2 = A_B1R14[bn1_ma_reg7];
			 data_tmp_3 = A_B1R15[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg7] = A_B0R12[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg7] = A_B0R12[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg7] = A_B1R12[bn1_ma_reg8];
			 A_B0R12[bn0_ma_reg7] = data_tmp_1;
			 A_B0R12[bn0_ma_reg8] = data_tmp_2;
			 A_B1R12[bn1_ma_reg8] = data_tmp_3;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B0R14[bn0_ma_reg7];
			 data_tmp_2 = A_B0R15[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg7] = A_B0R13[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg7] = A_B1R13[bn1_ma_reg8];
			 A_B0R13[bn0_ma_reg8] = data_tmp_1;
			 A_B1R13[bn1_ma_reg8] = data_tmp_2;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B0R15[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg8] = A_B1R14[bn1_ma_reg8];
			 A_B1R14[bn1_ma_reg8] = data_tmp_1;
			 //-----------------------------------------------------------------
		  }		
	}

    //--------------------------------------------------------------------------------------------
	//After data relocation for INTT output
	for(int i = 0;i < group;i++){
		for(int j = 0; j < radix; j++){
        	gray_i  = Gray(i,group);
        	BC_tmp  = j * group + gray_i;
        	BC_RR = BC_tmp;
			BC_IFFT_Reorder_R16_R2(BC_RR,BC_Reorder);
        	AGU_R16(BC_Reorder,bn_tmp,ma_tmp);
			AFTER_DATA_RELOCATION <<"--------------------------------------------- \n";
			AFTER_DATA_RELOCATION <<"i: "     << i <<"\n";
			AFTER_DATA_RELOCATION <<"j: "     << j <<"\n";
			AFTER_DATA_RELOCATION <<"gray_i: "<< gray_i <<"\n";			
			AFTER_DATA_RELOCATION <<"BC_tmp: "<< BC_tmp <<"\n";
			AFTER_DATA_RELOCATION <<"BC_RR: " << BC_RR  <<"\n";
			AFTER_DATA_RELOCATION <<"BC_Reorder: "<< BC_Reorder <<"\n";
			AFTER_DATA_RELOCATION <<"bn_tmp: "<< bn_tmp <<"\n";
			AFTER_DATA_RELOCATION <<"ma_tmp: "<< ma_tmp <<"\n";			
        	if(bn_tmp == 0){
				AFTER_DATA_RELOCATION <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";						
        	}else {
				AFTER_DATA_RELOCATION <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";								
        	}	
		}
	}
	
	//-------------------------------------------------------------------------
	
	std::ofstream INTT_DATA_OUTPUT("./INTT_R16_R2_DATA_OUTPUT.txt");
	
	int Num_of_data_In_Group;
	int BC_Reorder_R2_Index;
	Num_of_data_In_Group = (int)N / group;
	//data output
	// SPMB data output , output function is "NTT_REORDERINDEX_R16_R2_OUT".		
	for(int i = 0;i < group;i++){
		for(int j = 0;j < radix;j++){
			BC_tmp  = j * group + i;
			INTT_DATA_OUTPUT <<"---------------------------------------------------\n";
			INTT_DATA_OUTPUT <<" group: " << group <<"\n";
			INTT_DATA_OUTPUT <<" i: "     << i <<"\n";
			INTT_DATA_OUTPUT <<" j: "     << j <<"\n";
			INTT_DATA_OUTPUT <<" BC_tmp: "<< BC_tmp <<"\n";
			BC_IFFT_Reorder_R16_R2(BC_tmp,BC_Reorder_R2_Index);
			INTT_DATA_OUTPUT <<"BC_Reorder_R2_Index: "<< BC_Reorder_R2_Index <<"\n";
			BC_IFFT_Reorder_R16_R2_OUT(BC_tmp,BC_Reorder);
			INTT_DATA_OUTPUT <<"BC_Reorder_OUT_Index: "<< BC_Reorder <<"\n";
			AGU_R16(BC_Reorder,bn_tmp,ma_tmp);
			INTT_DATA_OUTPUT <<" bn_tmp:  "<< bn_tmp <<"\n";
			INTT_DATA_OUTPUT <<" ma_tmp:  "<< ma_tmp <<"\n";
			index0_i  = i * Num_of_data_In_Group+ 16 * j;
			index1_i  = i * Num_of_data_In_Group+ 16 * j  + 1;
			index2_i  = i * Num_of_data_In_Group+ 16 * j  + 2;
			index3_i  = i * Num_of_data_In_Group+ 16 * j  + 3;
			index4_i  = i * Num_of_data_In_Group+ 16 * j  + 4;
			index5_i  = i * Num_of_data_In_Group+ 16 * j  + 5;
			index6_i  = i * Num_of_data_In_Group+ 16 * j  + 6;
			index7_i  = i * Num_of_data_In_Group+ 16 * j  + 7;
			index8_i  = i * Num_of_data_In_Group+ 16 * j  + 8;
			index9_i  = i * Num_of_data_In_Group+ 16 * j  + 9;
			index10_i = i * Num_of_data_In_Group+ 16 * j + 10;
			index11_i = i * Num_of_data_In_Group+ 16 * j + 11;
			index12_i = i * Num_of_data_In_Group+ 16 * j + 12;
			index13_i = i * Num_of_data_In_Group+ 16 * j + 13;
			index14_i = i * Num_of_data_In_Group+ 16 * j + 14;
			index15_i = i * Num_of_data_In_Group+ 16 * j + 15;
			
			INTT_DATA_OUTPUT <<" index0_i:  "<< index0_i <<"\n";
			INTT_DATA_OUTPUT <<" index1_i:  "<< index1_i <<"\n";
			INTT_DATA_OUTPUT <<" index2_i:  "<< index2_i <<"\n";
			INTT_DATA_OUTPUT <<" index3_i:  "<< index3_i <<"\n";
			INTT_DATA_OUTPUT <<" index4_i:  "<< index4_i <<"\n";
			INTT_DATA_OUTPUT <<" index5_i:  "<< index5_i <<"\n";
			INTT_DATA_OUTPUT <<" index6_i:  "<< index6_i <<"\n";
			INTT_DATA_OUTPUT <<" index7_i:  "<< index7_i <<"\n";
			INTT_DATA_OUTPUT <<" index8_i:  "<< index8_i <<"\n";
			INTT_DATA_OUTPUT <<" index9_i:  "<< index9_i <<"\n";
			INTT_DATA_OUTPUT <<" index10_i: "<< index10_i <<"\n";
			INTT_DATA_OUTPUT <<" index11_i: "<< index11_i <<"\n";
			INTT_DATA_OUTPUT <<" index12_i: "<< index12_i <<"\n";
			INTT_DATA_OUTPUT <<" index13_i: "<< index13_i <<"\n";
			INTT_DATA_OUTPUT <<" index14_i: "<< index14_i <<"\n";
			INTT_DATA_OUTPUT <<" index15_i: "<< index15_i <<"\n";
					
			if(bn_tmp == 0){
			   index0  = index0_i;
			   index1  = index1_i;
			   index2  = index2_i;
			   index3  = index3_i;
			   index4  = index4_i;
			   index5  = index5_i;
			   index6  = index6_i;
			   index7  = index7_i;
			   index8  = index8_i;
			   index9  = index9_i;
			   index10 = index10_i;
			   index11 = index11_i;
			   index12 = index12_i;
			   index13 = index13_i;
			   index14 = index14_i;
			   index15 = index15_i;
               A[index0] = A_B0R0[ma_tmp];
			   A[index1] = A_B0R1[ma_tmp];
			   A[index2] = A_B0R2[ma_tmp];
			   A[index3] = A_B0R3[ma_tmp];
			   A[index4] = A_B0R4[ma_tmp];
			   A[index5] = A_B0R5[ma_tmp];
			   A[index6] = A_B0R6[ma_tmp];
			   A[index7] = A_B0R7[ma_tmp];
			   A[index8] = A_B0R8[ma_tmp];
			   A[index9] = A_B0R9[ma_tmp];
			   A[index10] = A_B0R10[ma_tmp];
			   A[index11] = A_B0R11[ma_tmp];
			   A[index12] = A_B0R12[ma_tmp];
			   A[index13] = A_B0R13[ma_tmp];
			   A[index14] = A_B0R14[ma_tmp];
			   A[index15] = A_B0R15[ma_tmp];
			   INTT_DATA_OUTPUT <<" A["<< index0  <<"]: "<< A_B0R0[ma_tmp] <<"\n";
			   INTT_DATA_OUTPUT <<" A["<< index1  <<"]: "<< A_B0R1[ma_tmp] <<"\n";
			   INTT_DATA_OUTPUT <<" A["<< index2  <<"]: "<< A_B0R2[ma_tmp] <<"\n";
			   INTT_DATA_OUTPUT <<" A["<< index3  <<"]: "<< A_B0R3[ma_tmp] <<"\n";
			   INTT_DATA_OUTPUT <<" A["<< index4  <<"]: "<< A_B0R4[ma_tmp] <<"\n";
			   INTT_DATA_OUTPUT <<" A["<< index5  <<"]: "<< A_B0R5[ma_tmp] <<"\n";
			   INTT_DATA_OUTPUT <<" A["<< index6  <<"]: "<< A_B0R6[ma_tmp] <<"\n";
			   INTT_DATA_OUTPUT <<" A["<< index7  <<"]: "<< A_B0R7[ma_tmp] <<"\n";
			   INTT_DATA_OUTPUT <<" A["<< index8  <<"]: "<< A_B0R8[ma_tmp] <<"\n";
			   INTT_DATA_OUTPUT <<" A["<< index9  <<"]: "<< A_B0R9[ma_tmp] <<"\n";
			   INTT_DATA_OUTPUT <<" A["<< index10 <<"]: "<< A_B0R10[ma_tmp] <<"\n";
			   INTT_DATA_OUTPUT <<" A["<< index11 <<"]: "<< A_B0R11[ma_tmp] <<"\n";
			   INTT_DATA_OUTPUT <<" A["<< index12 <<"]: "<< A_B0R12[ma_tmp] <<"\n";
			   INTT_DATA_OUTPUT <<" A["<< index13 <<"]: "<< A_B0R13[ma_tmp] <<"\n";
			   INTT_DATA_OUTPUT <<" A["<< index14 <<"]: "<< A_B0R14[ma_tmp] <<"\n";
			   INTT_DATA_OUTPUT <<" A["<< index15 <<"]: "<< A_B0R15[ma_tmp] <<"\n";			   
			}
			else {
			   index0  = index0_i;
			   index1  = index1_i;
			   index2  = index2_i;
			   index3  = index3_i;
			   index4  = index4_i;
			   index5  = index5_i;
			   index6  = index6_i;
			   index7  = index7_i;
			   index8  = index8_i;
			   index9  = index9_i;
			   index10 = index10_i;
			   index11 = index11_i;
			   index12 = index12_i;
			   index13 = index13_i;
			   index14 = index14_i;
			   index15 = index15_i;
               A[index0]     = A_B1R0[ma_tmp];
               A[index1]     = A_B1R1[ma_tmp];
               A[index2]     = A_B1R2[ma_tmp];
               A[index3]     = A_B1R3[ma_tmp];
               A[index4]     = A_B1R4[ma_tmp];
               A[index5]     = A_B1R5[ma_tmp];
               A[index6]     = A_B1R6[ma_tmp];
               A[index7]     = A_B1R7[ma_tmp];
               A[index8]     = A_B1R8[ma_tmp];
               A[index9]     = A_B1R9[ma_tmp];
               A[index10]    = A_B1R10[ma_tmp];
               A[index11]    = A_B1R11[ma_tmp];
               A[index12]    = A_B1R12[ma_tmp];
               A[index13]    = A_B1R13[ma_tmp];
               A[index14]    = A_B1R14[ma_tmp];
               A[index15]    = A_B1R15[ma_tmp];
			   INTT_DATA_OUTPUT <<" A["<< index0  <<"]: "<< A_B1R0[ma_tmp] <<"\n";
			   INTT_DATA_OUTPUT <<" A["<< index1  <<"]: "<< A_B1R1[ma_tmp] <<"\n";
			   INTT_DATA_OUTPUT <<" A["<< index2  <<"]: "<< A_B1R2[ma_tmp] <<"\n";
			   INTT_DATA_OUTPUT <<" A["<< index3  <<"]: "<< A_B1R3[ma_tmp] <<"\n";
			   INTT_DATA_OUTPUT <<" A["<< index4  <<"]: "<< A_B1R4[ma_tmp] <<"\n";
			   INTT_DATA_OUTPUT <<" A["<< index5  <<"]: "<< A_B1R5[ma_tmp] <<"\n";
			   INTT_DATA_OUTPUT <<" A["<< index6  <<"]: "<< A_B1R6[ma_tmp] <<"\n";
			   INTT_DATA_OUTPUT <<" A["<< index7  <<"]: "<< A_B1R7[ma_tmp] <<"\n";
			   INTT_DATA_OUTPUT <<" A["<< index8  <<"]: "<< A_B1R8[ma_tmp] <<"\n";
			   INTT_DATA_OUTPUT <<" A["<< index9  <<"]: "<< A_B1R9[ma_tmp] <<"\n";
			   INTT_DATA_OUTPUT <<" A["<< index10 <<"]: "<< A_B1R10[ma_tmp] <<"\n";
			   INTT_DATA_OUTPUT <<" A["<< index11 <<"]: "<< A_B1R11[ma_tmp] <<"\n";
			   INTT_DATA_OUTPUT <<" A["<< index12 <<"]: "<< A_B1R12[ma_tmp] <<"\n";
			   INTT_DATA_OUTPUT <<" A["<< index13 <<"]: "<< A_B1R13[ma_tmp] <<"\n";
			   INTT_DATA_OUTPUT <<" A["<< index14 <<"]: "<< A_B1R14[ma_tmp] <<"\n";
			   INTT_DATA_OUTPUT <<" A["<< index15 <<"]: "<< A_B1R15[ma_tmp] <<"\n";				   
			}			
		}		
	}
}
//-----------------
//mixed radix fft
//radix 16 and raidx-4
void NTTSPMB::NTT_r16_r4(std::vector<ZZ> &A,std::vector<ZZ> &B0R0,std::vector<ZZ> &B0R1,std::vector<ZZ> &B0R2,std::vector<ZZ> &B0R3,
std::vector<ZZ> &B0R4,std::vector<ZZ> &B0R5,std::vector<ZZ> &B0R6,std::vector<ZZ> &B0R7,std::vector<ZZ> &B0R8,std::vector<ZZ> &B0R9,
std::vector<ZZ> &B0R10,std::vector<ZZ> &B0R11,std::vector<ZZ> &B0R12,std::vector<ZZ> &B0R13,std::vector<ZZ> &B0R14,std::vector<ZZ> &B0R15,
std::vector<ZZ> &B1R0,std::vector<ZZ> &B1R1,std::vector<ZZ> &B1R2,std::vector<ZZ> &B1R3,std::vector<ZZ> &B1R4,std::vector<ZZ> &B1R5,
std::vector<ZZ> &B1R6,std::vector<ZZ> &B1R7,std::vector<ZZ> &B1R8,std::vector<ZZ> &B1R9,std::vector<ZZ> &B1R10,std::vector<ZZ> &B1R11,
std::vector<ZZ> &B1R12,std::vector<ZZ> &B1R13,std::vector<ZZ> &B1R14,std::vector<ZZ> &B1R15){
	
	int            Stage;//stage	
	int            word_size;
	int            offset;
	int            BC;
	int            BC_tmp;
    int            group;
    int            ma_tmp;
	int            bn_tmp;
    int            bit_tmp;
	int            bn0_bc_tmp;//frist in bc data
	int            bn1_bc_tmp;//frist in bc data
	int            bn0_ma_reg1;
	int            bn0_ma_reg2;
	int            bn0_ma_reg3;
	int            bn0_ma_reg4;
	int            bn0_ma_reg5;
	int            bn0_ma_reg6;
	int            bn0_ma_reg7;
	int            bn0_ma_reg8;
	int            bn1_ma_reg1;
	int            bn1_ma_reg2;
	int            bn1_ma_reg3;
	int            bn1_ma_reg4;
	int            bn1_ma_reg5;
	int            bn1_ma_reg6;
	int            bn1_ma_reg7;
	int            bn1_ma_reg8;
	int            gray_i;
    int            BC_WIDTH;
	int            tw_modulus;
    int            tw_modulus_tmp;
	//display
	int            display;
    double         Stage_double;	
    std::vector<int> bit_array_tmp;
	
	std::ofstream DATARECORD("./NTT_R16_R4_SPMB.txt");
	std::ofstream siang_record("./my_print_data/r16_r4_SPMB.txt");
	std::ofstream r16_r4_SPMB_TF_dg("./my_print_data/r16_r4_SPMB_TF_dg.txt");

	std::ofstream spmb_r16_r4("./SPMB_tw/spmb_r16_r4.txt");
	std::ofstream DTFAG_golden_st0("./SPMB_tw/DTFAG_golden_st0_16384.txt");
	std::ofstream DTFAG_golden_st1("./SPMB_tw/DTFAG_golden_st1_16384.txt");
	std::ofstream DTFAG_golden_st2("./SPMB_tw/DTFAG_golden_st2_16384.txt");
	//**********************************************************
	//**********************************************************
	display  = 1;
	//**********************************************************

	//-------------Bit Reverse--------------
	BitOperate RR, IntToVec, VecToInt;
	long long bit_width = log2(N);
	long long bit_width_s = log(radix)/log(2);
	vector<long long> bit_array;
	//----------------------------------------
	
	//radix-16 Stage
    Stage_double  = log2(N);
    Stage       =  (int)floor(Stage_double/4);
	BC_WIDTH    =  (int)ceil(log2(N/16));
	offset      =  (int)N /  16;
	word_size   =  (int)N / (2 * 16);
	group       =  (int)N / (256);
	tw_modulus  =  (int)N /  16;
	bit_array_tmp.resize(BC_WIDTH);

    DATARECORD << "group: "    << group << "\n";
    DATARECORD << "BC_WIDTH: " << BC_WIDTH << "\n";
	siang_record << "group: "    << group << "\n";	// siang record
    siang_record << "BC_WIDTH: " << BC_WIDTH << "\n";	// siagn record
	
	ZZ               data_tmp_1;
	ZZ               data_tmp_2;
	ZZ               data_tmp_3;
	ZZ               data_tmp_4;
	ZZ               data_tmp_5;
	ZZ               data_tmp_6;
	ZZ               data_tmp_7;
	ZZ               data_tmp_8;
	ZZ               data_tmp_9;
	ZZ               data_tmp_10;
	ZZ               data_tmp_11;
	ZZ               data_tmp_12;
	ZZ               data_tmp_13;
	ZZ               data_tmp_14;
	ZZ               data_tmp_15;
	std::vector<ZZ>  A_B0R0;
	std::vector<ZZ>  A_B0R1;
	std::vector<ZZ>  A_B0R2;
	std::vector<ZZ>  A_B0R3;
	std::vector<ZZ>  A_B0R4;
	std::vector<ZZ>  A_B0R5;
	std::vector<ZZ>  A_B0R6;
	std::vector<ZZ>  A_B0R7;
	std::vector<ZZ>  A_B0R8;
	std::vector<ZZ>  A_B0R9;
	std::vector<ZZ>  A_B0R10;
	std::vector<ZZ>  A_B0R11;
	std::vector<ZZ>  A_B0R12;
	std::vector<ZZ>  A_B0R13;
	std::vector<ZZ>  A_B0R14;
	std::vector<ZZ>  A_B0R15;
	std::vector<ZZ>  A_B1R0;
	std::vector<ZZ>  A_B1R1;
	std::vector<ZZ>  A_B1R2;
	std::vector<ZZ>  A_B1R3;
	std::vector<ZZ>  A_B1R4;
	std::vector<ZZ>  A_B1R5;
	std::vector<ZZ>  A_B1R6;
	std::vector<ZZ>  A_B1R7;
	std::vector<ZZ>  A_B1R8;
	std::vector<ZZ>  A_B1R9;
	std::vector<ZZ>  A_B1R10;
	std::vector<ZZ>  A_B1R11;
	std::vector<ZZ>  A_B1R12;
	std::vector<ZZ>  A_B1R13;
	std::vector<ZZ>  A_B1R14;
	std::vector<ZZ>  A_B1R15;
	
	A_B0R0.resize(word_size);
	A_B0R1.resize(word_size);
	A_B0R2.resize(word_size);
	A_B0R3.resize(word_size);
	A_B0R4.resize(word_size);
	A_B0R5.resize(word_size);
	A_B0R6.resize(word_size);
	A_B0R7.resize(word_size);
	A_B0R8.resize(word_size);
	A_B0R9.resize(word_size);
	A_B0R10.resize(word_size);
	A_B0R11.resize(word_size);
	A_B0R12.resize(word_size);
	A_B0R13.resize(word_size);
	A_B0R14.resize(word_size);
	A_B0R15.resize(word_size);
	A_B1R0.resize(word_size);
	A_B1R1.resize(word_size);
	A_B1R2.resize(word_size);
	A_B1R3.resize(word_size);
	A_B1R4.resize(word_size);
	A_B1R5.resize(word_size);
	A_B1R6.resize(word_size);
	A_B1R7.resize(word_size);
	A_B1R8.resize(word_size);
	A_B1R9.resize(word_size);
	A_B1R10.resize(word_size);
	A_B1R11.resize(word_size);
	A_B1R12.resize(word_size);
	A_B1R13.resize(word_size);
	A_B1R14.resize(word_size);
	A_B1R15.resize(word_size);
	//----------------------------------------------------
	B0R0.resize(word_size);
	B0R1.resize(word_size);
	B0R2.resize(word_size);
	B0R3.resize(word_size);
	B0R4.resize(word_size);
	B0R5.resize(word_size);
	B0R6.resize(word_size);
	B0R7.resize(word_size);
	B0R8.resize(word_size);
	B0R9.resize(word_size);
	B0R10.resize(word_size);
	B0R11.resize(word_size);
	B0R12.resize(word_size);
	B0R13.resize(word_size);
	B0R14.resize(word_size);
	B0R15.resize(word_size);
	B1R0.resize(word_size);
	B1R1.resize(word_size);
	B1R2.resize(word_size);
	B1R3.resize(word_size);
	B1R4.resize(word_size);
	B1R5.resize(word_size);
	B1R6.resize(word_size);
	B1R7.resize(word_size);
	B1R8.resize(word_size);
	B1R9.resize(word_size);
	B1R10.resize(word_size);
	B1R11.resize(word_size);
	B1R12.resize(word_size);
	B1R13.resize(word_size);
	B1R14.resize(word_size);
	B1R15.resize(word_size);
	//----------------------------------------------------
	int length;
	ZZ  factor;   //base factor
	ZZ  factor_t; //acctually mul factor
	ZZ  factor_2t;
	ZZ  factor_3t;
	ZZ  factor_4t;
	ZZ  factor_5t;
	ZZ  factor_6t;
	ZZ  factor_7t;
	ZZ  factor_8t;
	ZZ  factor_9t;
	ZZ  factor_10t;
	ZZ  factor_11t;
	ZZ  factor_12t;
	ZZ  factor_13t;
	ZZ  factor_14t;
	ZZ  factor_15t;
    
	//init load data
    for(int i = 0; i < group; i++){
		for(int j = 0 ; j < radix ; j++){
			bn_tmp = 0;
			ma_tmp = 0;
			gray_i  = Gray(i,group);
			BC = j * group + gray_i;
            //bit calculate
            for(int j=0; j < BC_WIDTH;j++){
                bit_tmp = BC % 2;
                BC = BC >> 1;
                bit_array_tmp[j] = bit_tmp;
            } 
			gray_i  = Gray(i,group);
			BC = j * group + gray_i;
            for(int rs = 0; rs < BC_WIDTH; rs++){
                if((bit_array_tmp[rs] == 1) && (rs != 0)) ma_tmp = ma_tmp + exp2((rs-1)); 
                bn_tmp = bn_tmp ^ bit_array_tmp[rs];
            }
			if(bn_tmp == 0){
				A_B0R0[ma_tmp]  = A[BC];
				A_B0R1[ma_tmp]  = A[BC +      offset];
				A_B0R2[ma_tmp]  = A[BC + 2  * offset];
				A_B0R3[ma_tmp]  = A[BC + 3  * offset];
				A_B0R4[ma_tmp]  = A[BC + 4  * offset];
				A_B0R5[ma_tmp]  = A[BC + 5  * offset];
				A_B0R6[ma_tmp]  = A[BC + 6  * offset];
				A_B0R7[ma_tmp]  = A[BC + 7  * offset];
				A_B0R8[ma_tmp]  = A[BC + 8  * offset];
				A_B0R9[ma_tmp]  = A[BC + 9  * offset];
				A_B0R10[ma_tmp] = A[BC + 10 * offset];
				A_B0R11[ma_tmp] = A[BC + 11 * offset];
				A_B0R12[ma_tmp] = A[BC + 12 * offset];
				A_B0R13[ma_tmp] = A[BC + 13 * offset];
				A_B0R14[ma_tmp] = A[BC + 14 * offset];
				A_B0R15[ma_tmp] = A[BC + 15 * offset];
			}else {
				A_B1R0[ma_tmp]  = A[BC];
				A_B1R1[ma_tmp]  = A[BC +     offset];
				A_B1R2[ma_tmp]  = A[BC + 2 * offset];
				A_B1R3[ma_tmp]  = A[BC + 3 * offset];
				A_B1R4[ma_tmp]  = A[BC + 4 * offset];
				A_B1R5[ma_tmp]  = A[BC + 5 * offset];
				A_B1R6[ma_tmp]  = A[BC + 6 * offset];
				A_B1R7[ma_tmp]  = A[BC + 7 * offset];
				A_B1R8[ma_tmp]  = A[BC + 8 * offset];
				A_B1R9[ma_tmp]  = A[BC + 9 * offset];
				A_B1R10[ma_tmp] = A[BC + 10 * offset];
				A_B1R11[ma_tmp] = A[BC + 11 * offset];
				A_B1R12[ma_tmp] = A[BC + 12 * offset];
				A_B1R13[ma_tmp] = A[BC + 13 * offset];
				A_B1R14[ma_tmp] = A[BC + 14 * offset];
				A_B1R15[ma_tmp] = A[BC + 15 * offset];
			}
		}
	}
	ma_tmp = 0;
	bn_tmp = 0;
	BC     = 0;
	int tw_degree = 1;
	std::cout << "init load over! \n";
	if(display == 1)DATARECORD <<"radix-16 computing stage:  "<< Stage <<"\n";
	if(display == 1)siang_record <<"radix-16 computing stage:  "<< Stage <<"\n"; // siang_record
	//need modify to mult by twiddle factor
	for(int s = 0; s < Stage; s++){
		if(s == 0)factor = W;
		else {
			SqrMod(factor,factor,p);
			SqrMod(factor,factor,p);
			SqrMod(factor,factor,p);
			SqrMod(factor,factor,p);
			tw_degree = tw_degree * 16;
		}
		std::cout << "factor = " << factor << std::endl;
		if(display == 1)DATARECORD <<"---------------------------------\n";
		if(display == 1)DATARECORD <<"Now Stage: "<< s <<"\n";		
		if(display == 1)siang_record <<"---------------------------------\n";	// siang_record
		if(display == 1)siang_record <<"Now Stage: "<< s <<"\n";	// siang_record
		if(display == 1)siang_record <<"factor: "<< factor <<"\n";	// siang_record
		spmb_r16_r4 <<"Now Stage: "<< s <<"\n";
		spmb_r16_r4 <<"twiddle factor : "<< factor <<"\n";
		r16_r4_SPMB_TF_dg << "Now Stage : " << s << "\n";
		tw_modulus_tmp  = tw_modulus >> ( 4 * s);
		int tw_record = 0; // siang_record
		int cnt = 0;//siang
		for(int i = 0 ;i < group;i++){
			bn0_bc_tmp  = 0;
			bn1_bc_tmp  = 0;
			r16_r4_SPMB_TF_dg << "-------------i = "<< i <<"-------------" << "\n";
			spmb_r16_r4 << "-----------------i = " << i << "-------------------" << std::endl;
			//r16_r4_SPMB_TF_dg << "i = " << i << std::endl;
			for(int j = 0;j < radix;j++){
				gray_i  = Gray(i,group);
			    BC_tmp  = j * group + gray_i;
				if(display == 1)DATARECORD <<"---------------------------------\n";
				if(display == 1)DATARECORD << "BC_tmp: " << BC_tmp << "\n";	
				if(display == 1)siang_record <<"---------------------------------\n"; // siang_record	
				if(display == 1)siang_record << "i = " << i << ", j: " << j << ", group: " << "\n";	 // siang_record		
				if(display == 1)siang_record << "BC_tmp: " << BC_tmp << "\n";	 // siang_record
				if(s == Stage - 1) RR_R16_R4(BC_tmp,(4 * s - 2),BC);
				else RR_R16(BC_tmp,s,BC);
				if(display == 1)DATARECORD << "After RR_R16 , BC : " << BC << "\n";		
				if(display == 1)siang_record << "After RR_R16 , BC : " << BC << "\n";	// siang_record	
				length = BC % tw_modulus_tmp;
				siang_record << "length : " << length << "\n" ;
				r16_r4_SPMB_TF_dg << "BC = " << BC << ", tw_modulus_tmp = " << tw_modulus_tmp << ", tw_modulus = " << tw_modulus
				<< ", length : " << length << ", tw_dg: " <<  tw_degree * length << " = " << tw_degree << " * " << length <<"\n";
				PowerMod(factor_t,factor,length,p);
				if(display == 1)siang_record << "factor = " << factor << ", length : " << length << "\n";	// siang_record	
				AGU_R16(BC,bn_tmp,ma_tmp);
				if(display == 1)DATARECORD << "BN : " << bn_tmp << "\n";
				if(display == 1)DATARECORD << "MA : " << ma_tmp << "\n";	

				//-----------compute data idx-------------
				spmb_r16_r4 << "BC_tmp = " << BC_tmp << std::endl;
				
				IntToVec.IntToVec(BC_tmp, N, bit_array);
				rotate(bit_array.begin(), bit_array.begin()+ bit_width_s*s , bit_array.end());
				long long Data = VecToInt.VecToInt(bit_array, N);
				spmb_r16_r4 << "Data_index = ";
                spmb_r16_r4 << "( " ;					
				for(int k = 0; k < radix ; k++ ){
					spmb_r16_r4 << Data + k*(1<<(bit_width-bit_width_s-bit_width_s*s)) <<" ";	
				}
				spmb_r16_r4 << ") " << std::endl;
				spmb_r16_r4 << ", (w^" << 0 << ", w^" << tw_degree * length << ", w^" << tw_degree * length * 2 
				<< ", w^" << tw_degree * length * 3  << ", w^" << tw_degree * length * 4 << ", w^" << tw_degree * length * 5 
				<< ", w^" << tw_degree * length * 6  << ", w^" << tw_degree * length * 7 << ", w^" << tw_degree * length * 8 
				<< ", w^" << tw_degree * length * 9  << ", w^" << tw_degree * length * 10 << ", w^" << tw_degree * length * 11 
				<< ", w^" << tw_degree * length * 12  << ", w^" << tw_degree * length * 13 << ", w^" << tw_degree * length * 14 
				<< ", w^" << tw_degree * length * 15
				<< ")" <<std::endl;
				//-----------------------------------------	


				if(bn_tmp == 0){
					if(j < 2)bn0_bc_tmp = BC_tmp;
					PowerMod(factor_2t,factor_t,2,p);
					PowerMod(factor_3t,factor_t,3,p);
					PowerMod(factor_4t,factor_t,4,p);
					PowerMod(factor_5t,factor_t,5,p);
					PowerMod(factor_6t,factor_t,6,p);
					PowerMod(factor_7t,factor_t,7,p);
					PowerMod(factor_8t,factor_t,8,p);
					PowerMod(factor_9t,factor_t,9,p);
					PowerMod(factor_10t,factor_t,10,p);
					PowerMod(factor_11t,factor_t,11,p);
					PowerMod(factor_12t,factor_t,12,p);
					PowerMod(factor_13t,factor_t,13,p);
					PowerMod(factor_14t,factor_t,14,p);
					PowerMod(factor_15t,factor_t,15,p);
					if(display == 1)DATARECORD <<"factor_t: "  << factor_t <<"\n";
					if(display == 1)DATARECORD <<"factor_2t: " << factor_2t <<"\n";
					if(display == 1)DATARECORD <<"factor_3t: " << factor_3t <<"\n";
					if(display == 1)DATARECORD <<"factor_4t: " << factor_4t <<"\n";
					if(display == 1)DATARECORD <<"factor_5t: " << factor_5t <<"\n";
					if(display == 1)DATARECORD <<"factor_6t: " << factor_6t <<"\n";
					if(display == 1)DATARECORD <<"factor_7t: " << factor_7t <<"\n";
					if(display == 1)DATARECORD <<"factor_8t: " << factor_8t <<"\n";
					if(display == 1)DATARECORD <<"factor_9t: " << factor_9t <<"\n";
					if(display == 1)DATARECORD <<"factor_10t: "<< factor_10t <<"\n";
					if(display == 1)DATARECORD <<"factor_11t: "<< factor_11t <<"\n";
					if(display == 1)DATARECORD <<"factor_12t: "<< factor_12t <<"\n";
					if(display == 1)DATARECORD <<"factor_13t: "<< factor_13t <<"\n";
					if(display == 1)DATARECORD <<"factor_14t: "<< factor_14t <<"\n";
					if(display == 1)DATARECORD <<"factor_15t: "<< factor_15t <<"\n";	
					//----------------siang record-----------------------------------------
					if(display == 1)siang_record << "p : " << p << "\n";	// siang_record
					if(display == 1)siang_record << "j : " << j << "\n";	// siang_record
					if(display == 1)siang_record << "i : " << i << "\n";	// siang_record
					if(display == 1)siang_record << "factor_t : "   		   << factor_t << "\n";	// siang_record	
					if(display == 1)siang_record << "factor_2t_factor_3t : "   << factor_2t  << "_" << factor_3t << "\n";	// siang_record	
					if(display == 1)siang_record << "factor_4t_factor_5t : "   << factor_4t  << "_" << factor_5t << "\n";	// siang_record	
					if(display == 1)siang_record << "factor_6t_factor_7t : "   << factor_6t  << "_" << factor_7t << "\n";	// siang_record	
					if(display == 1)siang_record << "factor_8t_factor_9t : "   << factor_8t  << "_" << factor_9t << "\n";	// siang_record	
					if(display == 1)siang_record << "factor_10t_factor_11t : " << factor_10t << "_" << factor_11t << "\n";	// siang_record	
					if(display == 1)siang_record << "factor_12t_factor_13t : " << factor_12t << "_" << factor_13t << "\n";	// siang_record	
					if(display == 1)siang_record << "factor_14t_factor_15t : " << factor_14t << "_" << factor_15t << "\n";	// siang_record	


					//-----------------------------------------------------------------------

					//if(display == 1)DATARECORD << "Before Radix-16 butterfly unit operation!!! \n";
				    /*if(display == 1)DATARECORD << "A_B0R0["<< ma_tmp <<"]: " << A_B0R0[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B0R1["<< ma_tmp <<"]: " << A_B0R1[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B0R2["<< ma_tmp <<"]: " << A_B0R2[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B0R3["<< ma_tmp <<"]: " << A_B0R3[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B0R4["<< ma_tmp <<"]: " << A_B0R4[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B0R5["<< ma_tmp <<"]: " << A_B0R5[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B0R6["<< ma_tmp <<"]: " << A_B0R6[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B0R7["<< ma_tmp <<"]: " << A_B0R7[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B0R8["<< ma_tmp <<"]: " << A_B0R8[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B0R9["<< ma_tmp <<"]: " << A_B0R9[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B0R10["<< ma_tmp <<"]: "<< A_B0R10[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B0R11["<< ma_tmp <<"]: "<< A_B0R11[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B0R12["<< ma_tmp <<"]: "<< A_B0R12[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B0R13["<< ma_tmp <<"]: "<< A_B0R13[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B0R14["<< ma_tmp <<"]: "<< A_B0R14[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B0R15["<< ma_tmp <<"]: "<< A_B0R15[ma_tmp]<<"\n";	*/				
					Radix16_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp], A_B0R2[ma_tmp], A_B0R3[ma_tmp],A_B0R4[ma_tmp],
							   A_B0R5[ma_tmp],A_B0R6[ma_tmp], A_B0R7[ma_tmp], A_B0R8[ma_tmp],A_B0R9[ma_tmp],
							   A_B0R10[ma_tmp],A_B0R11[ma_tmp],A_B0R12[ma_tmp],A_B0R13[ma_tmp],A_B0R14[ma_tmp],
							   A_B0R15[ma_tmp]);
					//if(display == 1)DATARECORD << "After Radix-16 butterfly unit operation!!! \n";
				   /* if(display == 1)DATARECORD << "A_B0R0["<< ma_tmp <<"]: " << A_B0R0[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B0R1["<< ma_tmp <<"]: " << A_B0R1[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B0R2["<< ma_tmp <<"]: " << A_B0R2[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B0R3["<< ma_tmp <<"]: " << A_B0R3[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B0R4["<< ma_tmp <<"]: " << A_B0R4[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B0R5["<< ma_tmp <<"]: " << A_B0R5[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B0R6["<< ma_tmp <<"]: " << A_B0R6[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B0R7["<< ma_tmp <<"]: " << A_B0R7[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B0R8["<< ma_tmp <<"]: " << A_B0R8[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B0R9["<< ma_tmp <<"]: " << A_B0R9[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B0R10["<< ma_tmp <<"]: "<< A_B0R10[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B0R11["<< ma_tmp <<"]: "<< A_B0R11[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B0R12["<< ma_tmp <<"]: "<< A_B0R12[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B0R13["<< ma_tmp <<"]: "<< A_B0R13[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B0R14["<< ma_tmp <<"]: "<< A_B0R14[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B0R15["<< ma_tmp <<"]: "<< A_B0R15[ma_tmp]<<"\n";		*/					   
					MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],factor_t,p);
					MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],factor_2t,p);
					MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],factor_3t,p);
					MulMod(A_B0R4[ma_tmp],A_B0R4[ma_tmp],factor_4t,p);
					MulMod(A_B0R5[ma_tmp],A_B0R5[ma_tmp],factor_5t,p);
					MulMod(A_B0R6[ma_tmp],A_B0R6[ma_tmp],factor_6t,p);
					MulMod(A_B0R7[ma_tmp],A_B0R7[ma_tmp],factor_7t,p);
					MulMod(A_B0R8[ma_tmp],A_B0R8[ma_tmp],factor_8t,p);
					MulMod(A_B0R9[ma_tmp],A_B0R9[ma_tmp],factor_9t,p);
					MulMod(A_B0R10[ma_tmp],A_B0R10[ma_tmp],factor_10t,p);
					MulMod(A_B0R11[ma_tmp],A_B0R11[ma_tmp],factor_11t,p);
					MulMod(A_B0R12[ma_tmp],A_B0R12[ma_tmp],factor_12t,p);
					MulMod(A_B0R13[ma_tmp],A_B0R13[ma_tmp],factor_13t,p);
					MulMod(A_B0R14[ma_tmp],A_B0R14[ma_tmp],factor_14t,p);
					MulMod(A_B0R15[ma_tmp],A_B0R15[ma_tmp],factor_15t,p);
					//if(display == 1)DATARECORD <<" After mult by twiddle factor!!! \n";
					/*if(display == 1)DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";		*/				
					if(j <  2)bn0_ma_reg1 = ma_tmp;
					if((j >= 2)  && (j < 4))bn0_ma_reg2 = ma_tmp;
					if((j >= 4)  && (j < 6))bn0_ma_reg3 = ma_tmp;
					if((j >= 6)  && (j < 8))bn0_ma_reg4 = ma_tmp;
					if((j >= 8)  && (j < 10))bn0_ma_reg5 = ma_tmp;
					if((j >= 10) && (j < 12))bn0_ma_reg6 = ma_tmp;
					if((j >= 12) && (j < 14))bn0_ma_reg7 = ma_tmp;
					if((j >= 14) && (j < 16))bn0_ma_reg8 = ma_tmp;
					//----------------DTFAG golden pattern------------------
					switch(s){
						case 0:
							DTFAG_golden_st0 << 1 		  << " \n ";
							DTFAG_golden_st0 << factor_t  << " \n ";
							DTFAG_golden_st0 << factor_2t << " \n ";
							DTFAG_golden_st0 << factor_3t << " \n ";
							DTFAG_golden_st0 << factor_4t << " \n ";
							DTFAG_golden_st0 << factor_5t  << " \n ";
							DTFAG_golden_st0 << factor_6t << " \n ";
							DTFAG_golden_st0 << factor_7t << " \n ";
							DTFAG_golden_st0 << factor_8t << " \n ";
							DTFAG_golden_st0 << factor_9t  << " \n ";
							DTFAG_golden_st0 << factor_10t << " \n ";
							DTFAG_golden_st0 << factor_11t << " \n ";
							DTFAG_golden_st0 << factor_12t << " \n ";
							DTFAG_golden_st0 << factor_13t  << " \n ";
							DTFAG_golden_st0 << factor_14t << " \n ";
							DTFAG_golden_st0 << factor_15t << " \n ";
							break;
						case 1:
							DTFAG_golden_st1 << 1 		  << " \n ";
							DTFAG_golden_st1 << factor_t  << " \n ";
							DTFAG_golden_st1 << factor_2t << " \n ";
							DTFAG_golden_st1 << factor_3t << " \n ";
							DTFAG_golden_st1 << factor_4t << " \n ";
							DTFAG_golden_st1 << factor_5t  << " \n ";
							DTFAG_golden_st1 << factor_6t << " \n ";
							DTFAG_golden_st1 << factor_7t << " \n ";
							DTFAG_golden_st1 << factor_8t << " \n ";
							DTFAG_golden_st1 << factor_9t  << " \n ";
							DTFAG_golden_st1 << factor_10t << " \n ";
							DTFAG_golden_st1 << factor_11t << " \n ";
							DTFAG_golden_st1 << factor_12t << " \n ";
							DTFAG_golden_st1 << factor_13t  << " \n ";
							DTFAG_golden_st1 << factor_14t << " \n ";
							DTFAG_golden_st1 << factor_15t << " \n ";
							break;
						case 2:
							DTFAG_golden_st2 << 1 		  << " \n ";
							DTFAG_golden_st2 << factor_t  << " \n ";
							DTFAG_golden_st2 << factor_2t << " \n ";
							DTFAG_golden_st2 << factor_3t << " \n ";
							DTFAG_golden_st2 << factor_4t << " \n ";
							DTFAG_golden_st2 << factor_5t  << " \n ";
							DTFAG_golden_st2 << factor_6t << " \n ";
							DTFAG_golden_st2 << factor_7t << " \n ";
							DTFAG_golden_st2 << factor_8t << " \n ";
							DTFAG_golden_st2 << factor_9t  << " \n ";
							DTFAG_golden_st2 << factor_10t << " \n ";
							DTFAG_golden_st2 << factor_11t << " \n ";
							DTFAG_golden_st2 << factor_12t << " \n ";
							DTFAG_golden_st2 << factor_13t  << " \n ";
							DTFAG_golden_st2 << factor_14t << " \n ";
							DTFAG_golden_st2 << factor_15t << " \n ";
							break;
						default:
							break;
					}
					//--------------------------------------------------
				}
			    else {
					if(j < 2)bn1_bc_tmp = BC_tmp;
					PowerMod(factor_2t,factor_t,2,p);
					PowerMod(factor_3t,factor_t,3,p);
					PowerMod(factor_4t,factor_t,4,p);
					PowerMod(factor_5t,factor_t,5,p);
					PowerMod(factor_6t,factor_t,6,p);
					PowerMod(factor_7t,factor_t,7,p);
					PowerMod(factor_8t,factor_t,8,p);
					PowerMod(factor_9t,factor_t,9,p);
					PowerMod(factor_10t,factor_t,10,p);
					PowerMod(factor_11t,factor_t,11,p);
					PowerMod(factor_12t,factor_t,12,p);
					PowerMod(factor_13t,factor_t,13,p);
					PowerMod(factor_14t,factor_t,14,p);
					PowerMod(factor_15t,factor_t,15,p);
					if(display == 1)DATARECORD <<"factor_t: "  << factor_t <<"\n";
					if(display == 1)DATARECORD <<"factor_2t: " << factor_2t <<"\n";
					if(display == 1)DATARECORD <<"factor_3t: " << factor_3t <<"\n";
					if(display == 1)DATARECORD <<"factor_4t: " << factor_4t <<"\n";
					if(display == 1)DATARECORD <<"factor_5t: " << factor_5t <<"\n";
					if(display == 1)DATARECORD <<"factor_6t: " << factor_6t <<"\n";
					if(display == 1)DATARECORD <<"factor_7t: " << factor_7t <<"\n";
					if(display == 1)DATARECORD <<"factor_8t: " << factor_8t <<"\n";
					if(display == 1)DATARECORD <<"factor_9t: " << factor_9t <<"\n";
					if(display == 1)DATARECORD <<"factor_10t: "<< factor_10t <<"\n";
					if(display == 1)DATARECORD <<"factor_11t: "<< factor_11t <<"\n";
					if(display == 1)DATARECORD <<"factor_12t: "<< factor_12t <<"\n";
					if(display == 1)DATARECORD <<"factor_13t: "<< factor_13t <<"\n";
					if(display == 1)DATARECORD <<"factor_14t: "<< factor_14t <<"\n";
					if(display == 1)DATARECORD <<"factor_15t: "<< factor_15t <<"\n";	
					//----------------siang record-----------------------------------------
					if(display == 1)siang_record << "p : " << p << "\n";	// siang_record
					if(display == 1)siang_record << "j : " << j << "\n";	// siang_record
					if(display == 1)siang_record << "i : " << i << "\n";	// siang_record
					if(display == 1)siang_record << "factor_t : "   		   << factor_t << "\n";	// siang_record	
					if(display == 1)siang_record << "factor_2t_factor_3t : "   << factor_2t  << "_"  << factor_3t << "\n";	// siang_record	
					if(display == 1)siang_record << "factor_4t_factor_5t : "   << factor_4t  << "_"  << factor_5t << "\n";	// siang_record	
					if(display == 1)siang_record << "factor_6t_factor_7t : "   << factor_6t  << "_"  << factor_7t << "\n";	// siang_record	
					if(display == 1)siang_record << "factor_8t_factor_9t : "   << factor_8t  << "_"  << factor_9t << "\n";	// siang_record	
					if(display == 1)siang_record << "factor_10t_factor_11t : " << factor_10t << "_"  << factor_11t << "\n";	// siang_record	
					if(display == 1)siang_record << "factor_12t_factor_13t : " << factor_12t << "_"  << factor_13t << "\n";	// siang_record	
					if(display == 1)siang_record << "factor_14t_factor_15t : " << factor_14t << "_"  << factor_15t << "\n";	// siang_record		
					//-----------------------------------------------------------------------

					//if(display == 1)DATARECORD << "Before Radix-16 butterfly unit operation!!! \n";
				   /* if(display == 1)DATARECORD << "A_B1R0["<< ma_tmp <<"]: " << A_B1R0[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B1R1["<< ma_tmp <<"]: " << A_B1R1[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B1R2["<< ma_tmp <<"]: " << A_B1R2[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B1R3["<< ma_tmp <<"]: " << A_B1R3[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B1R4["<< ma_tmp <<"]: " << A_B1R4[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B1R5["<< ma_tmp <<"]: " << A_B1R5[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B1R6["<< ma_tmp <<"]: " << A_B1R6[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B1R7["<< ma_tmp <<"]: " << A_B1R7[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B1R8["<< ma_tmp <<"]: " << A_B1R8[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B1R9["<< ma_tmp <<"]: " << A_B1R9[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B1R10["<< ma_tmp <<"]: "<< A_B1R10[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B1R11["<< ma_tmp <<"]: "<< A_B1R11[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B1R12["<< ma_tmp <<"]: "<< A_B1R12[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B1R13["<< ma_tmp <<"]: "<< A_B1R13[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B1R14["<< ma_tmp <<"]: "<< A_B1R14[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B1R15["<< ma_tmp <<"]: "<< A_B1R15[ma_tmp]<<"\n";		*/				
					Radix16_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp], A_B1R2[ma_tmp], A_B1R3[ma_tmp],A_B1R4[ma_tmp],
							   A_B1R5[ma_tmp],A_B1R6[ma_tmp], A_B1R7[ma_tmp], A_B1R8[ma_tmp],A_B1R9[ma_tmp],
							   A_B1R10[ma_tmp],A_B1R11[ma_tmp],A_B1R12[ma_tmp],A_B1R13[ma_tmp],A_B1R14[ma_tmp],
							   A_B1R15[ma_tmp]);
					//if(display == 1)DATARECORD << "After Radix-16 butterfly unit operation!!! \n";
				    /*if(display == 1)DATARECORD << "A_B1R0["<< ma_tmp <<"]: " << A_B1R0[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B1R1["<< ma_tmp <<"]: " << A_B1R1[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B1R2["<< ma_tmp <<"]: " << A_B1R2[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B1R3["<< ma_tmp <<"]: " << A_B1R3[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B1R4["<< ma_tmp <<"]: " << A_B1R4[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B1R5["<< ma_tmp <<"]: " << A_B1R5[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B1R6["<< ma_tmp <<"]: " << A_B1R6[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B1R7["<< ma_tmp <<"]: " << A_B1R7[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B1R8["<< ma_tmp <<"]: " << A_B1R8[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B1R9["<< ma_tmp <<"]: " << A_B1R9[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B1R10["<< ma_tmp <<"]: "<< A_B1R10[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B1R11["<< ma_tmp <<"]: "<< A_B1R11[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B1R12["<< ma_tmp <<"]: "<< A_B1R12[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B1R13["<< ma_tmp <<"]: "<< A_B1R13[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B1R14["<< ma_tmp <<"]: "<< A_B1R14[ma_tmp]<<"\n";
				    if(display == 1)DATARECORD << "A_B1R15["<< ma_tmp <<"]: "<< A_B1R15[ma_tmp]<<"\n";*/							   
					MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],factor_t,p);
					MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],factor_2t,p);
					MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],factor_3t,p);
					MulMod(A_B1R4[ma_tmp],A_B1R4[ma_tmp],factor_4t,p);
					MulMod(A_B1R5[ma_tmp],A_B1R5[ma_tmp],factor_5t,p);
					MulMod(A_B1R6[ma_tmp],A_B1R6[ma_tmp],factor_6t,p);
					MulMod(A_B1R7[ma_tmp],A_B1R7[ma_tmp],factor_7t,p);
					MulMod(A_B1R8[ma_tmp],A_B1R8[ma_tmp],factor_8t,p);
					MulMod(A_B1R9[ma_tmp],A_B1R9[ma_tmp],factor_9t,p);
					MulMod(A_B1R10[ma_tmp],A_B1R10[ma_tmp],factor_10t,p);
					MulMod(A_B1R11[ma_tmp],A_B1R11[ma_tmp],factor_11t,p);
					MulMod(A_B1R12[ma_tmp],A_B1R12[ma_tmp],factor_12t,p);
					MulMod(A_B1R13[ma_tmp],A_B1R13[ma_tmp],factor_13t,p);
					MulMod(A_B1R14[ma_tmp],A_B1R14[ma_tmp],factor_14t,p);
					MulMod(A_B1R15[ma_tmp],A_B1R15[ma_tmp],factor_15t,p);
					//if(display == 1)DATARECORD <<" After mult by twiddle factor!!! \n";
					/*if(display == 1)DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";	*/			
                    if(j <  2)bn1_ma_reg1 = ma_tmp;					
                    if((j >= 2)  && (j < 4))bn1_ma_reg2 = ma_tmp;
                    if((j >= 4)  && (j < 6))bn1_ma_reg3 = ma_tmp;
                    if((j >= 6)  && (j < 8))bn1_ma_reg4 = ma_tmp;
                    if((j >= 8)  && (j < 10))bn1_ma_reg5 = ma_tmp;
                    if((j >= 10) && (j < 12))bn1_ma_reg6 = ma_tmp;
                    if((j >= 12) && (j < 14))bn1_ma_reg7 = ma_tmp;
                    if((j >= 14) && (j < 16))bn1_ma_reg8 = ma_tmp;	
					//----------------DTFAG golden pattern------------------
					switch(s){
						case 0:
							DTFAG_golden_st0 << 1 		  << " \n ";
							DTFAG_golden_st0 << factor_t  << " \n ";
							DTFAG_golden_st0 << factor_2t << " \n ";
							DTFAG_golden_st0 << factor_3t << " \n ";
							DTFAG_golden_st0 << factor_4t << " \n ";
							DTFAG_golden_st0 << factor_5t  << " \n ";
							DTFAG_golden_st0 << factor_6t << " \n ";
							DTFAG_golden_st0 << factor_7t << " \n ";
							DTFAG_golden_st0 << factor_8t << " \n ";
							DTFAG_golden_st0 << factor_9t  << " \n ";
							DTFAG_golden_st0 << factor_10t << " \n ";
							DTFAG_golden_st0 << factor_11t << " \n ";
							DTFAG_golden_st0 << factor_12t << " \n ";
							DTFAG_golden_st0 << factor_13t  << " \n ";
							DTFAG_golden_st0 << factor_14t << " \n ";
							DTFAG_golden_st0 << factor_15t << " \n ";
							break;
						case 1:
							DTFAG_golden_st1 << 1 		  << " \n ";
							DTFAG_golden_st1 << factor_t  << " \n ";
							DTFAG_golden_st1 << factor_2t << " \n ";
							DTFAG_golden_st1 << factor_3t << " \n ";
							DTFAG_golden_st1 << factor_4t << " \n ";
							DTFAG_golden_st1 << factor_5t  << " \n ";
							DTFAG_golden_st1 << factor_6t << " \n ";
							DTFAG_golden_st1 << factor_7t << " \n ";
							DTFAG_golden_st1 << factor_8t << " \n ";
							DTFAG_golden_st1 << factor_9t  << " \n ";
							DTFAG_golden_st1 << factor_10t << " \n ";
							DTFAG_golden_st1 << factor_11t << " \n ";
							DTFAG_golden_st1 << factor_12t << " \n ";
							DTFAG_golden_st1 << factor_13t  << " \n ";
							DTFAG_golden_st1 << factor_14t << " \n ";
							DTFAG_golden_st1 << factor_15t << " \n ";
							break;
						case 2:
							DTFAG_golden_st2 << 1 		  << " \n ";
							DTFAG_golden_st2 << factor_t  << " \n ";
							DTFAG_golden_st2 << factor_2t << " \n ";
							DTFAG_golden_st2 << factor_3t << " \n ";
							DTFAG_golden_st2 << factor_4t << " \n ";
							DTFAG_golden_st2 << factor_5t  << " \n ";
							DTFAG_golden_st2 << factor_6t << " \n ";
							DTFAG_golden_st2 << factor_7t << " \n ";
							DTFAG_golden_st2 << factor_8t << " \n ";
							DTFAG_golden_st2 << factor_9t  << " \n ";
							DTFAG_golden_st2 << factor_10t << " \n ";
							DTFAG_golden_st2 << factor_11t << " \n ";
							DTFAG_golden_st2 << factor_12t << " \n ";
							DTFAG_golden_st2 << factor_13t  << " \n ";
							DTFAG_golden_st2 << factor_14t << " \n ";
							DTFAG_golden_st2 << factor_15t << " \n ";
							break;
						default:
							break;
					}
					//--------------------------------------------------				
				}
			}
		//data relocation
		 if(s < Stage-1){
		  if(bn1_bc_tmp > bn0_bc_tmp){
			 data_tmp_1  = A_B0R1[bn0_ma_reg1];
			 data_tmp_2  = A_B0R2[bn0_ma_reg1];
			 data_tmp_3  = A_B0R3[bn0_ma_reg1];
			 data_tmp_4  = A_B0R4[bn0_ma_reg1];
			 data_tmp_5  = A_B0R5[bn0_ma_reg1];
			 data_tmp_6  = A_B0R6[bn0_ma_reg1];
			 data_tmp_7  = A_B0R7[bn0_ma_reg1];
			 data_tmp_8  = A_B0R8[bn0_ma_reg1];
			 data_tmp_9  = A_B0R9[bn0_ma_reg1];
			 data_tmp_10 = A_B0R10[bn0_ma_reg1];
			 data_tmp_11 = A_B0R11[bn0_ma_reg1];
			 data_tmp_12 = A_B0R12[bn0_ma_reg1];
			 data_tmp_13 = A_B0R13[bn0_ma_reg1];
			 data_tmp_14 = A_B0R14[bn0_ma_reg1];
			 data_tmp_15 = A_B0R15[bn0_ma_reg1];
			 A_B0R1[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg1];
			 A_B0R2[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg2];
			 A_B0R3[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg2];
			 A_B0R4[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg3];
			 A_B0R5[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg3];
			 A_B0R6[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg1] = A_B0R0[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg1] = A_B1R0[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg1] = A_B0R0[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg1] = A_B1R0[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg1] = A_B1R0[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg1] = A_B0R0[bn0_ma_reg8];
			 A_B1R0[bn1_ma_reg1]  = data_tmp_1; 
			 A_B1R0[bn1_ma_reg2]  = data_tmp_2; 
			 A_B0R0[bn0_ma_reg2]  = data_tmp_3; 
			 A_B1R0[bn1_ma_reg3]  = data_tmp_4; 
			 A_B0R0[bn0_ma_reg3]  = data_tmp_5; 
			 A_B0R0[bn0_ma_reg4]  = data_tmp_6; 
			 A_B1R0[bn1_ma_reg4]  = data_tmp_7; 
			 A_B1R0[bn1_ma_reg5]  = data_tmp_8; 
			 A_B0R0[bn0_ma_reg5]  = data_tmp_9; 
			 A_B0R0[bn0_ma_reg6]  = data_tmp_10;
			 A_B1R0[bn1_ma_reg6]  = data_tmp_11;
			 A_B0R0[bn0_ma_reg7]  = data_tmp_12;
			 A_B1R0[bn1_ma_reg7]  = data_tmp_13;
			 A_B1R0[bn1_ma_reg8]  = data_tmp_14;
			 A_B0R0[bn0_ma_reg8]  = data_tmp_15;
			 /*********************************************/
			 data_tmp_1  = A_B1R2[bn1_ma_reg1];
			 data_tmp_2  = A_B1R3[bn1_ma_reg1];
			 data_tmp_3  = A_B1R4[bn1_ma_reg1];
			 data_tmp_4  = A_B1R5[bn1_ma_reg1];
			 data_tmp_5  = A_B1R6[bn1_ma_reg1];
			 data_tmp_6  = A_B1R7[bn1_ma_reg1];
			 data_tmp_7  = A_B1R8[bn1_ma_reg1];
			 data_tmp_8  = A_B1R9[bn1_ma_reg1];
			 data_tmp_9  = A_B1R10[bn1_ma_reg1];
			 data_tmp_10 = A_B1R11[bn1_ma_reg1];
			 data_tmp_11 = A_B1R12[bn1_ma_reg1];
			 data_tmp_12 = A_B1R13[bn1_ma_reg1];
			 data_tmp_13 = A_B1R14[bn1_ma_reg1];
			 data_tmp_14 = A_B1R15[bn1_ma_reg1];
             A_B1R2[bn1_ma_reg1]  =	 A_B1R1[bn1_ma_reg2];
             A_B1R3[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg2];
             A_B1R4[bn1_ma_reg1]  =  A_B1R1[bn1_ma_reg3];
             A_B1R5[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg3];
             A_B1R6[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg4];
             A_B1R7[bn1_ma_reg1]  =  A_B1R1[bn1_ma_reg4];
             A_B1R8[bn1_ma_reg1]  =  A_B1R1[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg1] =  A_B0R1[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg1] =  A_B1R1[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg1] =  A_B0R1[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg1] =  A_B1R1[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg1] =  A_B1R1[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg1] =  A_B0R1[bn0_ma_reg8];
			 A_B1R1[bn1_ma_reg2]  =  data_tmp_1; 
			 A_B0R1[bn0_ma_reg2]  =  data_tmp_2; 
			 A_B1R1[bn1_ma_reg3]  =  data_tmp_3; 
			 A_B0R1[bn0_ma_reg3]  =  data_tmp_4; 
			 A_B0R1[bn0_ma_reg4]  =  data_tmp_5; 
			 A_B1R1[bn1_ma_reg4]  =  data_tmp_6; 
			 A_B1R1[bn1_ma_reg5]  =  data_tmp_7; 
			 A_B0R1[bn0_ma_reg5]  =  data_tmp_8; 
			 A_B0R1[bn0_ma_reg6]  =  data_tmp_9; 
			 A_B1R1[bn1_ma_reg6]  =  data_tmp_10;
			 A_B0R1[bn0_ma_reg7]  =  data_tmp_11;
			 A_B1R1[bn1_ma_reg7]  =  data_tmp_12;
			 A_B1R1[bn1_ma_reg8]  =  data_tmp_13;
			 A_B0R1[bn0_ma_reg8]  =  data_tmp_14;
			/************************************************************/ 
			 data_tmp_1  =   A_B1R3[bn1_ma_reg2];
			 data_tmp_2  =   A_B1R4[bn1_ma_reg2];
			 data_tmp_3  =   A_B1R5[bn1_ma_reg2];
			 data_tmp_4  =   A_B1R6[bn1_ma_reg2];
			 data_tmp_5  =   A_B1R7[bn1_ma_reg2];
			 data_tmp_6  =   A_B1R8[bn1_ma_reg2];
			 data_tmp_7  =   A_B1R9[bn1_ma_reg2];
			 data_tmp_8  =   A_B1R10[bn1_ma_reg2];
			 data_tmp_9  =   A_B1R11[bn1_ma_reg2];
			 data_tmp_10 =   A_B1R12[bn1_ma_reg2];
			 data_tmp_11 =   A_B1R13[bn1_ma_reg2];
			 data_tmp_12 =   A_B1R14[bn1_ma_reg2];
			 data_tmp_13 =   A_B1R15[bn1_ma_reg2];
			 A_B1R3[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg2];
			 A_B1R4[bn1_ma_reg2]  =  A_B1R2[bn1_ma_reg3];
			 A_B1R5[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg3];
			 A_B1R6[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg4];
			 A_B1R7[bn1_ma_reg2]  =  A_B1R2[bn1_ma_reg4];
			 A_B1R8[bn1_ma_reg2]  =  A_B1R2[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg2] =  A_B0R2[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg2] =  A_B1R2[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg2] =  A_B0R2[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg2] =  A_B1R2[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg2] =  A_B1R2[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg2] =  A_B0R2[bn0_ma_reg8];
			 A_B0R2[bn0_ma_reg2]  =  data_tmp_1; 
			 A_B1R2[bn1_ma_reg3]  =  data_tmp_2; 
			 A_B0R2[bn0_ma_reg3]  =  data_tmp_3; 
			 A_B0R2[bn0_ma_reg4]  =  data_tmp_4; 
			 A_B1R2[bn1_ma_reg4]  =  data_tmp_5; 
			 A_B1R2[bn1_ma_reg5]  =  data_tmp_6; 
			 A_B0R2[bn0_ma_reg5]  =  data_tmp_7; 
			 A_B0R2[bn0_ma_reg6]  =  data_tmp_8; 
			 A_B1R2[bn1_ma_reg6]  =  data_tmp_9; 
			 A_B0R2[bn0_ma_reg7]  =  data_tmp_10;
			 A_B1R2[bn1_ma_reg7]  =  data_tmp_11;
			 A_B1R2[bn1_ma_reg8]  =  data_tmp_12;
			 A_B0R2[bn0_ma_reg8]  =  data_tmp_13;
			 //******************************************************
			 data_tmp_1  = A_B0R4[bn0_ma_reg2];
			 data_tmp_2  = A_B0R5[bn0_ma_reg2];
			 data_tmp_3  = A_B0R6[bn0_ma_reg2];
			 data_tmp_4  = A_B0R7[bn0_ma_reg2];
			 data_tmp_5  = A_B0R8[bn0_ma_reg2];
			 data_tmp_6  = A_B0R9[bn0_ma_reg2];
			 data_tmp_7  = A_B0R10[bn0_ma_reg2];
			 data_tmp_8  = A_B0R11[bn0_ma_reg2];
			 data_tmp_9  = A_B0R12[bn0_ma_reg2];
			 data_tmp_10 = A_B0R13[bn0_ma_reg2];
			 data_tmp_11 = A_B0R14[bn0_ma_reg2];
			 data_tmp_12 = A_B0R15[bn0_ma_reg2];
			 A_B0R4[bn0_ma_reg2]  =  A_B1R3[bn1_ma_reg3];
			 A_B0R5[bn0_ma_reg2]  =  A_B0R3[bn0_ma_reg3];
			 A_B0R6[bn0_ma_reg2]  =  A_B0R3[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg2]  =  A_B1R3[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg2]  =  A_B1R3[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg2]  =  A_B0R3[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg2] =  A_B0R3[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg2] =  A_B1R3[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg2] =  A_B0R3[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg2] =  A_B1R3[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg2] =  A_B1R3[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg2] =  A_B0R3[bn0_ma_reg8];
			 A_B1R3[bn1_ma_reg3]  =  data_tmp_1;  
			 A_B0R3[bn0_ma_reg3]  =  data_tmp_2;  
			 A_B0R3[bn0_ma_reg4]  =  data_tmp_3;  
			 A_B1R3[bn1_ma_reg4]  =  data_tmp_4;  
			 A_B1R3[bn1_ma_reg5]  =  data_tmp_5;  
			 A_B0R3[bn0_ma_reg5]  =  data_tmp_6;  
			 A_B0R3[bn0_ma_reg6]  =  data_tmp_7;  
			 A_B1R3[bn1_ma_reg6]  =  data_tmp_8;  
			 A_B0R3[bn0_ma_reg7]  =  data_tmp_9;  
			 A_B1R3[bn1_ma_reg7]  =  data_tmp_10; 
			 A_B1R3[bn1_ma_reg8]  =  data_tmp_11; 
			 A_B0R3[bn0_ma_reg8]  =  data_tmp_12; 
			 //----------------------------------------------------------------------
             data_tmp_1  = A_B1R5[bn1_ma_reg3];
			 data_tmp_2  = A_B1R6[bn1_ma_reg3];
			 data_tmp_3  = A_B1R7[bn1_ma_reg3];
			 data_tmp_4  = A_B1R8[bn1_ma_reg3];
			 data_tmp_5  = A_B1R9[bn1_ma_reg3];
			 data_tmp_6  = A_B1R10[bn1_ma_reg3];
			 data_tmp_7  = A_B1R11[bn1_ma_reg3];
			 data_tmp_8  = A_B1R12[bn1_ma_reg3];
			 data_tmp_9  = A_B1R13[bn1_ma_reg3];
			 data_tmp_10 = A_B1R14[bn1_ma_reg3];
			 data_tmp_11 = A_B1R15[bn1_ma_reg3];
			 A_B1R5[bn1_ma_reg3]  =  A_B0R4[bn0_ma_reg3];
			 A_B1R6[bn1_ma_reg3]  =  A_B0R4[bn0_ma_reg4];
			 A_B1R7[bn1_ma_reg3]  =  A_B1R4[bn1_ma_reg4];
			 A_B1R8[bn1_ma_reg3]  =  A_B1R4[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg3]  =  A_B0R4[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg3] =  A_B0R4[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg3] =  A_B1R4[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg3] =  A_B0R4[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg3] =  A_B1R4[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg3] =  A_B1R4[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg3] =  A_B0R4[bn0_ma_reg8];
			 A_B0R4[bn0_ma_reg3]  =  data_tmp_1 ;
			 A_B0R4[bn0_ma_reg4]  =  data_tmp_2 ;
			 A_B1R4[bn1_ma_reg4]  =  data_tmp_3 ;
			 A_B1R4[bn1_ma_reg5]  =  data_tmp_4 ;
			 A_B0R4[bn0_ma_reg5]  =  data_tmp_5 ;
			 A_B0R4[bn0_ma_reg6]  =  data_tmp_6 ;
			 A_B1R4[bn1_ma_reg6]  =  data_tmp_7 ;
			 A_B0R4[bn0_ma_reg7]  =  data_tmp_8 ;
			 A_B1R4[bn1_ma_reg7]  =  data_tmp_9 ;
			 A_B1R4[bn1_ma_reg8]  =  data_tmp_10;
			 A_B0R4[bn0_ma_reg8]  =  data_tmp_11;
			 //************************************************************************
			 data_tmp_1  = A_B0R6[bn0_ma_reg3];
			 data_tmp_2  = A_B0R7[bn0_ma_reg3];
			 data_tmp_3  = A_B0R8[bn0_ma_reg3];
			 data_tmp_4  = A_B0R9[bn0_ma_reg3];
			 data_tmp_5  = A_B0R10[bn0_ma_reg3];
			 data_tmp_6  = A_B0R11[bn0_ma_reg3];
			 data_tmp_7  = A_B0R12[bn0_ma_reg3];
			 data_tmp_8  = A_B0R13[bn0_ma_reg3];
			 data_tmp_9  = A_B0R14[bn0_ma_reg3];
			 data_tmp_10 = A_B0R15[bn0_ma_reg3];
			 A_B0R6[bn0_ma_reg3]  =  A_B0R5[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg3]  =  A_B1R5[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg3]  =  A_B1R5[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg3]  =  A_B0R5[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg3] =  A_B0R5[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg3] =  A_B1R5[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg3] =  A_B0R5[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg3] =  A_B1R5[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg3] =  A_B1R5[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg3] =  A_B0R5[bn0_ma_reg8];
			 A_B0R5[bn0_ma_reg4]  =  data_tmp_1; 
			 A_B1R5[bn1_ma_reg4]  =  data_tmp_2; 
			 A_B1R5[bn1_ma_reg5]  =  data_tmp_3; 
			 A_B0R5[bn0_ma_reg5]  =  data_tmp_4; 
			 A_B0R5[bn0_ma_reg6]  =  data_tmp_5; 
			 A_B1R5[bn1_ma_reg6]  =  data_tmp_6; 
			 A_B0R5[bn0_ma_reg7]  =  data_tmp_7; 
			 A_B1R5[bn1_ma_reg7]  =  data_tmp_8; 
			 A_B1R5[bn1_ma_reg8]  =  data_tmp_9; 
			 A_B0R5[bn0_ma_reg8]  =  data_tmp_10;
			 //-----------------------------------------------------------------------
			 data_tmp_1  = A_B0R7[bn0_ma_reg4];
			 data_tmp_2  = A_B0R8[bn0_ma_reg4];
			 data_tmp_3  = A_B0R9[bn0_ma_reg4];
			 data_tmp_4  = A_B0R10[bn0_ma_reg4];
			 data_tmp_5  = A_B0R11[bn0_ma_reg4];
			 data_tmp_6  = A_B0R12[bn0_ma_reg4];
			 data_tmp_7  = A_B0R13[bn0_ma_reg4];
			 data_tmp_8  = A_B0R14[bn0_ma_reg4];
			 data_tmp_9  = A_B0R15[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg4]  =  A_B1R6[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg4]  =  A_B1R6[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg4]  =  A_B0R6[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg4] =  A_B0R6[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg4] =  A_B1R6[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg4] =  A_B0R6[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg4] =  A_B1R6[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg4] =  A_B1R6[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg4] =  A_B0R6[bn0_ma_reg8];
			 A_B1R6[bn1_ma_reg4]  =  data_tmp_1;
			 A_B1R6[bn1_ma_reg5]  =  data_tmp_2;
			 A_B0R6[bn0_ma_reg5]  =  data_tmp_3;
			 A_B0R6[bn0_ma_reg6]  =  data_tmp_4;
			 A_B1R6[bn1_ma_reg6]  =  data_tmp_5;
			 A_B0R6[bn0_ma_reg7]  =  data_tmp_6;
			 A_B1R6[bn1_ma_reg7]  =  data_tmp_7;
			 A_B1R6[bn1_ma_reg8]  =  data_tmp_8;
			 A_B0R6[bn0_ma_reg8]  =  data_tmp_9;
			 //----------------------------------------------------------------------
			 data_tmp_1  = A_B1R8[bn1_ma_reg4];
			 data_tmp_2  = A_B1R9[bn1_ma_reg4];
			 data_tmp_3  = A_B1R10[bn1_ma_reg4];
			 data_tmp_4  = A_B1R11[bn1_ma_reg4];
			 data_tmp_5  = A_B1R12[bn1_ma_reg4];
			 data_tmp_6  = A_B1R13[bn1_ma_reg4];
			 data_tmp_7  = A_B1R14[bn1_ma_reg4];
			 data_tmp_8  = A_B1R15[bn1_ma_reg4];
			 A_B1R8[bn1_ma_reg4]  =  A_B1R7[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg4]  =  A_B0R7[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg4] =  A_B0R7[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg4] =  A_B1R7[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg4] =  A_B0R7[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg4] =  A_B1R7[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg4] =  A_B1R7[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg4] =  A_B0R7[bn0_ma_reg8];
			 A_B1R7[bn1_ma_reg5]  =  data_tmp_1;
			 A_B0R7[bn0_ma_reg5]  =  data_tmp_2;
			 A_B0R7[bn0_ma_reg6]  =  data_tmp_3;
			 A_B1R7[bn1_ma_reg6]  =  data_tmp_4;
			 A_B0R7[bn0_ma_reg7]  =  data_tmp_5;
			 A_B1R7[bn1_ma_reg7]  =  data_tmp_6;
			 A_B1R7[bn1_ma_reg8]  =  data_tmp_7;
			 A_B0R7[bn0_ma_reg8]  =  data_tmp_8;
			 //----------------------------------------------------------------------
			 data_tmp_1 = A_B1R9[bn1_ma_reg5];
			 data_tmp_2 = A_B1R10[bn1_ma_reg5];
			 data_tmp_3 = A_B1R11[bn1_ma_reg5];
			 data_tmp_4 = A_B1R12[bn1_ma_reg5];
			 data_tmp_5 = A_B1R13[bn1_ma_reg5];
			 data_tmp_6 = A_B1R14[bn1_ma_reg5];
			 data_tmp_7 = A_B1R15[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg5]  = A_B0R8[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg5] = A_B0R8[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg5] = A_B1R8[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg5] = A_B0R8[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg5] = A_B1R8[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg5] = A_B1R8[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg5] = A_B0R8[bn0_ma_reg8];
			 A_B0R8[bn0_ma_reg5]  = data_tmp_1;
			 A_B0R8[bn0_ma_reg6]  = data_tmp_2;
			 A_B1R8[bn1_ma_reg6]  = data_tmp_3;
			 A_B0R8[bn0_ma_reg7]  = data_tmp_4;
			 A_B1R8[bn1_ma_reg7]  = data_tmp_5;
			 A_B1R8[bn1_ma_reg8]  = data_tmp_6;
			 A_B0R8[bn0_ma_reg8]  = data_tmp_7;
			 //---------------------------------------------------------------------
			 data_tmp_1  = A_B0R10[bn0_ma_reg5];
			 data_tmp_2  = A_B0R11[bn0_ma_reg5];
			 data_tmp_3  = A_B0R12[bn0_ma_reg5];
			 data_tmp_4  = A_B0R13[bn0_ma_reg5];
			 data_tmp_5  = A_B0R14[bn0_ma_reg5];
			 data_tmp_6  = A_B0R15[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg5] = A_B0R9[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg5] = A_B1R9[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg5] = A_B0R9[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg5] = A_B1R9[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg5] = A_B1R9[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg5] = A_B0R9[bn0_ma_reg8];
			 A_B0R9[bn0_ma_reg6]  = data_tmp_1;
			 A_B1R9[bn1_ma_reg6]  = data_tmp_2;
			 A_B0R9[bn0_ma_reg7]  = data_tmp_3;
			 A_B1R9[bn1_ma_reg7]  = data_tmp_4;
			 A_B1R9[bn1_ma_reg8]  = data_tmp_5;
			 A_B0R9[bn0_ma_reg8]  = data_tmp_6;
			 //--------------------------------------------------------------------
			 data_tmp_1  = A_B0R11[bn0_ma_reg6];
			 data_tmp_2  = A_B0R12[bn0_ma_reg6];
			 data_tmp_3  = A_B0R13[bn0_ma_reg6];
			 data_tmp_4  = A_B0R14[bn0_ma_reg6];
			 data_tmp_5  = A_B0R15[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg6] = A_B1R10[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg6] = A_B0R10[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg6] = A_B1R10[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg6] = A_B1R10[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg6] = A_B0R10[bn0_ma_reg8];
			 A_B1R10[bn1_ma_reg6] = data_tmp_1;
			 A_B0R10[bn0_ma_reg7] = data_tmp_2;
			 A_B1R10[bn1_ma_reg7] = data_tmp_3;
			 A_B1R10[bn1_ma_reg8] = data_tmp_4;
			 A_B0R10[bn0_ma_reg8] = data_tmp_5;
			 //--------------------------------------------------------------------
			 data_tmp_1  = A_B1R12[bn1_ma_reg6];
			 data_tmp_2  = A_B1R13[bn1_ma_reg6];
			 data_tmp_3  = A_B1R14[bn1_ma_reg6];
			 data_tmp_4  = A_B1R15[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg6] = A_B0R11[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg6] = A_B1R11[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg6] = A_B1R11[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg6] = A_B0R11[bn0_ma_reg8];
			 A_B0R11[bn0_ma_reg7] = data_tmp_1;
			 A_B1R11[bn1_ma_reg7] = data_tmp_2;
			 A_B1R11[bn1_ma_reg8] = data_tmp_3;
			 A_B0R11[bn0_ma_reg8] = data_tmp_4;
			 //--------------------------------------------------------------------
			 data_tmp_1 = A_B0R13[bn0_ma_reg7];
			 data_tmp_2 = A_B0R14[bn0_ma_reg7];
			 data_tmp_3 = A_B0R15[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg7] = A_B1R12[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg7] = A_B1R12[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg7] = A_B0R12[bn0_ma_reg8];
			 A_B1R12[bn1_ma_reg7] = data_tmp_1;
			 A_B1R12[bn1_ma_reg8] = data_tmp_2;
			 A_B0R12[bn0_ma_reg8] = data_tmp_3;
			 //--------------------------------------------------------------------
			 data_tmp_1 = A_B1R14[bn1_ma_reg7];
			 data_tmp_2 = A_B1R15[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg7] = A_B1R13[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg7] = A_B0R13[bn0_ma_reg8];
			 A_B1R13[bn1_ma_reg8] = data_tmp_1;
			 A_B0R13[bn0_ma_reg8] = data_tmp_2;
			 //--------------------------------------------------------------------
			 data_tmp_1 = A_B1R15[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg8] = A_B0R14[bn0_ma_reg8];
			 A_B0R14[bn0_ma_reg8] = data_tmp_1; 
		  }else{	 
			 data_tmp_1  = A_B1R1[bn1_ma_reg1];
			 data_tmp_2  = A_B1R2[bn1_ma_reg1];
			 data_tmp_3  = A_B1R3[bn1_ma_reg1];
			 data_tmp_4  = A_B1R4[bn1_ma_reg1];
			 data_tmp_5  = A_B1R5[bn1_ma_reg1];
			 data_tmp_6  = A_B1R6[bn1_ma_reg1];
			 data_tmp_7  = A_B1R7[bn1_ma_reg1];
			 data_tmp_8  = A_B1R8[bn1_ma_reg1];
			 data_tmp_9  = A_B1R9[bn1_ma_reg1];
			 data_tmp_10 = A_B1R10[bn1_ma_reg1];
			 data_tmp_11 = A_B1R11[bn1_ma_reg1];
			 data_tmp_12 = A_B1R12[bn1_ma_reg1];
			 data_tmp_13 = A_B1R13[bn1_ma_reg1];
			 data_tmp_14 = A_B1R14[bn1_ma_reg1];
			 data_tmp_15 = A_B1R15[bn1_ma_reg1];
             A_B1R1[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg1];
             A_B1R2[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg2];
             A_B1R3[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg2];
             A_B1R4[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg3];
             A_B1R5[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg3];
             A_B1R6[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg4];
             A_B1R7[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg4];
             A_B1R8[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg5];
             A_B1R9[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg5];
             A_B1R10[bn1_ma_reg1] =  A_B1R0[bn1_ma_reg6];
             A_B1R11[bn1_ma_reg1] =  A_B0R0[bn0_ma_reg6];
             A_B1R12[bn1_ma_reg1] =  A_B1R0[bn1_ma_reg7];
             A_B1R13[bn1_ma_reg1] =  A_B0R0[bn0_ma_reg7];
             A_B1R14[bn1_ma_reg1] =  A_B0R0[bn0_ma_reg8];
             A_B1R15[bn1_ma_reg1] =  A_B1R0[bn1_ma_reg8];
             A_B0R0[bn0_ma_reg1]  =  data_tmp_1; 
             A_B0R0[bn0_ma_reg2]  =  data_tmp_2; 
             A_B1R0[bn1_ma_reg2]  =  data_tmp_3; 
             A_B0R0[bn0_ma_reg3]  =  data_tmp_4; 
             A_B1R0[bn1_ma_reg3]  =  data_tmp_5; 
             A_B1R0[bn1_ma_reg4]  =  data_tmp_6; 
             A_B0R0[bn0_ma_reg4]  =  data_tmp_7; 
             A_B0R0[bn0_ma_reg5]  =  data_tmp_8; 
             A_B1R0[bn1_ma_reg5]  =  data_tmp_9; 
             A_B1R0[bn1_ma_reg6]  =  data_tmp_10;
             A_B0R0[bn0_ma_reg6]  =  data_tmp_11;
             A_B1R0[bn1_ma_reg7]  =  data_tmp_12;
             A_B0R0[bn0_ma_reg7]  =  data_tmp_13;
             A_B0R0[bn0_ma_reg8]  =  data_tmp_14;
             A_B1R0[bn1_ma_reg8]  =  data_tmp_15;
             //-----------------------------------------------------------
             data_tmp_1  = A_B0R2[bn0_ma_reg1];
			 data_tmp_2  = A_B0R3[bn0_ma_reg1];
			 data_tmp_3  = A_B0R4[bn0_ma_reg1];
			 data_tmp_4  = A_B0R5[bn0_ma_reg1];
			 data_tmp_5  = A_B0R6[bn0_ma_reg1];
			 data_tmp_6  = A_B0R7[bn0_ma_reg1];
			 data_tmp_7  = A_B0R8[bn0_ma_reg1];
			 data_tmp_8  = A_B0R9[bn0_ma_reg1];
			 data_tmp_9  = A_B0R10[bn0_ma_reg1];
			 data_tmp_10 = A_B0R11[bn0_ma_reg1];
			 data_tmp_11 = A_B0R12[bn0_ma_reg1];
			 data_tmp_12 = A_B0R13[bn0_ma_reg1];
			 data_tmp_13 = A_B0R14[bn0_ma_reg1];
			 data_tmp_14 = A_B0R15[bn0_ma_reg1];
			 A_B0R2[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg2];
			 A_B0R3[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg2];
			 A_B0R4[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg3];
			 A_B0R5[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg3];
			 A_B0R6[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg4];
			 A_B0R7[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg1] =  A_B1R1[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg1] =  A_B0R1[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg1] =  A_B1R1[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg1] =  A_B0R1[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg1] =  A_B0R1[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg1] =  A_B1R1[bn1_ma_reg8];
             A_B0R1[bn0_ma_reg2]  =  data_tmp_1; 
             A_B1R1[bn1_ma_reg2]  =  data_tmp_2; 
             A_B0R1[bn0_ma_reg3]  =  data_tmp_3; 
             A_B1R1[bn1_ma_reg3]  =  data_tmp_4; 
             A_B1R1[bn1_ma_reg4]  =  data_tmp_5; 
             A_B0R1[bn0_ma_reg4]  =  data_tmp_6; 
             A_B0R1[bn0_ma_reg5]  =  data_tmp_7; 
             A_B1R1[bn1_ma_reg5]  =  data_tmp_8; 
             A_B1R1[bn1_ma_reg6]  =  data_tmp_9; 
             A_B0R1[bn0_ma_reg6]  =  data_tmp_10;
             A_B1R1[bn1_ma_reg7]  =  data_tmp_11;
             A_B0R1[bn0_ma_reg7]  =  data_tmp_12;
             A_B0R1[bn0_ma_reg8]  =  data_tmp_13;
             A_B1R1[bn1_ma_reg8]  =  data_tmp_14;
             //------------------------------------------------------------
             data_tmp_1  =  A_B0R3[bn0_ma_reg2];
			 data_tmp_2  =  A_B0R4[bn0_ma_reg2];
			 data_tmp_3  =  A_B0R5[bn0_ma_reg2];
			 data_tmp_4  =  A_B0R6[bn0_ma_reg2];
			 data_tmp_5  =  A_B0R7[bn0_ma_reg2];
			 data_tmp_6  =  A_B0R8[bn0_ma_reg2];
			 data_tmp_7  =  A_B0R9[bn0_ma_reg2];
			 data_tmp_8  =  A_B0R10[bn0_ma_reg2];
			 data_tmp_9  =  A_B0R11[bn0_ma_reg2];
			 data_tmp_10 =  A_B0R12[bn0_ma_reg2];
			 data_tmp_11 =  A_B0R13[bn0_ma_reg2];
			 data_tmp_12 =  A_B0R14[bn0_ma_reg2];
			 data_tmp_13 =  A_B0R15[bn0_ma_reg2];
			 A_B0R3[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg2];
			 A_B0R4[bn0_ma_reg2]  =  A_B0R2[bn0_ma_reg3];
			 A_B0R5[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg3];
			 A_B0R6[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg4];
			 A_B0R7[bn0_ma_reg2]  =  A_B0R2[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg2]  =  A_B0R2[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg2] =  A_B1R2[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg2] =  A_B0R2[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg2] =  A_B1R2[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg2] =  A_B0R2[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg2] =  A_B0R2[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg2] =  A_B1R2[bn1_ma_reg8];
			 A_B1R2[bn1_ma_reg2]  =  data_tmp_1; 
			 A_B0R2[bn0_ma_reg3]  =  data_tmp_2; 
			 A_B1R2[bn1_ma_reg3]  =  data_tmp_3; 
			 A_B1R2[bn1_ma_reg4]  =  data_tmp_4; 
			 A_B0R2[bn0_ma_reg4]  =  data_tmp_5; 
			 A_B0R2[bn0_ma_reg5]  =  data_tmp_6; 
			 A_B1R2[bn1_ma_reg5]  =  data_tmp_7; 
			 A_B1R2[bn1_ma_reg6]  =  data_tmp_8; 
			 A_B0R2[bn0_ma_reg6]  =  data_tmp_9; 
			 A_B1R2[bn1_ma_reg7]  =  data_tmp_10;
			 A_B0R2[bn0_ma_reg7]  =  data_tmp_11;
			 A_B0R2[bn0_ma_reg8]  =  data_tmp_12;
			 A_B1R2[bn1_ma_reg8]  =  data_tmp_13;
			 //-----------------------------------------------------------                        
			 data_tmp_1  = A_B1R4[bn1_ma_reg2];
			 data_tmp_2  = A_B1R5[bn1_ma_reg2];
			 data_tmp_3  = A_B1R6[bn1_ma_reg2];
			 data_tmp_4  = A_B1R7[bn1_ma_reg2];
			 data_tmp_5  = A_B1R8[bn1_ma_reg2];
			 data_tmp_6  = A_B1R9[bn1_ma_reg2];
			 data_tmp_7  = A_B1R10[bn1_ma_reg2];
			 data_tmp_8  = A_B1R11[bn1_ma_reg2];
			 data_tmp_9  = A_B1R12[bn1_ma_reg2];
			 data_tmp_10 = A_B1R13[bn1_ma_reg2];
			 data_tmp_11 = A_B1R14[bn1_ma_reg2];
			 data_tmp_12 = A_B1R15[bn1_ma_reg2];
			 A_B1R4[bn1_ma_reg2]  = A_B0R3[bn0_ma_reg3];
			 A_B1R5[bn1_ma_reg2]  = A_B1R3[bn1_ma_reg3];
			 A_B1R6[bn1_ma_reg2]  = A_B1R3[bn1_ma_reg4];
			 A_B1R7[bn1_ma_reg2]  = A_B0R3[bn0_ma_reg4];
			 A_B1R8[bn1_ma_reg2]  = A_B0R3[bn0_ma_reg5];
			 A_B1R9[bn1_ma_reg2]  = A_B1R3[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg2] = A_B1R3[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg2] = A_B0R3[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg2] = A_B1R3[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg2] = A_B0R3[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg2] = A_B0R3[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg2] = A_B1R3[bn1_ma_reg8];
			 A_B0R3[bn0_ma_reg3]  = data_tmp_1; 
			 A_B1R3[bn1_ma_reg3]  = data_tmp_2; 
			 A_B1R3[bn1_ma_reg4]  = data_tmp_3; 
			 A_B0R3[bn0_ma_reg4]  = data_tmp_4; 
			 A_B0R3[bn0_ma_reg5]  = data_tmp_5; 
			 A_B1R3[bn1_ma_reg5]  = data_tmp_6; 
			 A_B1R3[bn1_ma_reg6]  = data_tmp_7; 
			 A_B0R3[bn0_ma_reg6]  = data_tmp_8; 
			 A_B1R3[bn1_ma_reg7]  = data_tmp_9; 
			 A_B0R3[bn0_ma_reg7]  = data_tmp_10;
			 A_B0R3[bn0_ma_reg8]  = data_tmp_11;
			 A_B1R3[bn1_ma_reg8]  = data_tmp_12;
			 //------------------------------------------------------------
			 data_tmp_1  =  A_B0R5[bn0_ma_reg3];
			 data_tmp_2  =  A_B0R6[bn0_ma_reg3];
			 data_tmp_3  =  A_B0R7[bn0_ma_reg3];
			 data_tmp_4  =  A_B0R8[bn0_ma_reg3];
			 data_tmp_5  =  A_B0R9[bn0_ma_reg3];
			 data_tmp_6  =  A_B0R10[bn0_ma_reg3];
			 data_tmp_7  =  A_B0R11[bn0_ma_reg3];
			 data_tmp_8  =  A_B0R12[bn0_ma_reg3];
			 data_tmp_9  =  A_B0R13[bn0_ma_reg3];
			 data_tmp_10 =  A_B0R14[bn0_ma_reg3];
			 data_tmp_11 =  A_B0R15[bn0_ma_reg3];
			 A_B0R5[bn0_ma_reg3]  =  A_B1R4[bn1_ma_reg3];
			 A_B0R6[bn0_ma_reg3]  =  A_B1R4[bn1_ma_reg4];
			 A_B0R7[bn0_ma_reg3]  =  A_B0R4[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg3]  =  A_B0R4[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg3]  =  A_B1R4[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg3] =  A_B1R4[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg3] =  A_B0R4[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg3] =  A_B1R4[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg3] =  A_B0R4[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg3] =  A_B0R4[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg3] =  A_B1R4[bn1_ma_reg8];
			 A_B1R4[bn1_ma_reg3]  =  data_tmp_1; 
			 A_B1R4[bn1_ma_reg4]  =  data_tmp_2; 
			 A_B0R4[bn0_ma_reg4]  =  data_tmp_3; 
			 A_B0R4[bn0_ma_reg5]  =  data_tmp_4; 
			 A_B1R4[bn1_ma_reg5]  =  data_tmp_5; 
			 A_B1R4[bn1_ma_reg6]  =  data_tmp_6; 
			 A_B0R4[bn0_ma_reg6]  =  data_tmp_7; 
			 A_B1R4[bn1_ma_reg7]  =  data_tmp_8; 
			 A_B0R4[bn0_ma_reg7]  =  data_tmp_9; 
			 A_B0R4[bn0_ma_reg8]  =  data_tmp_10;
			 A_B1R4[bn1_ma_reg8]  =  data_tmp_11;
			 //-------------------------------------------------------------
			 data_tmp_1  =  A_B1R6[bn1_ma_reg3];
			 data_tmp_2  =  A_B1R7[bn1_ma_reg3];
			 data_tmp_3  =  A_B1R8[bn1_ma_reg3];
			 data_tmp_4  =  A_B1R9[bn1_ma_reg3];
			 data_tmp_5  =  A_B1R10[bn1_ma_reg3];
			 data_tmp_6  =  A_B1R11[bn1_ma_reg3];
			 data_tmp_7  =  A_B1R12[bn1_ma_reg3];
			 data_tmp_8  =  A_B1R13[bn1_ma_reg3];
			 data_tmp_9  =  A_B1R14[bn1_ma_reg3];
			 data_tmp_10 =  A_B1R15[bn1_ma_reg3];
			 A_B1R6[bn1_ma_reg3]  =  A_B1R5[bn1_ma_reg4];
			 A_B1R7[bn1_ma_reg3]  =  A_B0R5[bn0_ma_reg4];
			 A_B1R8[bn1_ma_reg3]  =  A_B0R5[bn0_ma_reg5];
			 A_B1R9[bn1_ma_reg3]  =  A_B1R5[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg3] =  A_B1R5[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg3] =  A_B0R5[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg3] =  A_B1R5[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg3] =  A_B0R5[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg3] =  A_B0R5[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg3] =  A_B1R5[bn1_ma_reg8];
			 A_B1R5[bn1_ma_reg4]  =  data_tmp_1; 
			 A_B0R5[bn0_ma_reg4]  =  data_tmp_2; 
			 A_B0R5[bn0_ma_reg5]  =  data_tmp_3; 
			 A_B1R5[bn1_ma_reg5]  =  data_tmp_4; 
			 A_B1R5[bn1_ma_reg6]  =  data_tmp_5; 
			 A_B0R5[bn0_ma_reg6]  =  data_tmp_6; 
			 A_B1R5[bn1_ma_reg7]  =  data_tmp_7; 
			 A_B0R5[bn0_ma_reg7]  =  data_tmp_8; 
			 A_B0R5[bn0_ma_reg8]  =  data_tmp_9; 
			 A_B1R5[bn1_ma_reg8]  =  data_tmp_10;
			 //----------------------------------------------------------------
			 data_tmp_1  =  A_B1R7[bn1_ma_reg4];
			 data_tmp_2  =  A_B1R8[bn1_ma_reg4];
			 data_tmp_3  =  A_B1R9[bn1_ma_reg4];
			 data_tmp_4  =  A_B1R10[bn1_ma_reg4];
			 data_tmp_5  =  A_B1R11[bn1_ma_reg4];
			 data_tmp_6  =  A_B1R12[bn1_ma_reg4];
			 data_tmp_7  =  A_B1R13[bn1_ma_reg4];
			 data_tmp_8  =  A_B1R14[bn1_ma_reg4];
			 data_tmp_9  =  A_B1R15[bn1_ma_reg4];
			 A_B1R7[bn1_ma_reg4]  =  A_B0R6[bn0_ma_reg4];
			 A_B1R8[bn1_ma_reg4]  =  A_B0R6[bn0_ma_reg5];
			 A_B1R9[bn1_ma_reg4]  =  A_B1R6[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg4] =  A_B1R6[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg4] =  A_B0R6[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg4] =  A_B1R6[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg4] =  A_B0R6[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg4] =  A_B0R6[bn0_ma_reg8];
             A_B1R15[bn1_ma_reg4] =  A_B1R6[bn1_ma_reg8];
             A_B0R6[bn0_ma_reg4]  =  data_tmp_1;
             A_B0R6[bn0_ma_reg5]  =  data_tmp_2;
             A_B1R6[bn1_ma_reg5]  =  data_tmp_3;
             A_B1R6[bn1_ma_reg6]  =  data_tmp_4;
             A_B0R6[bn0_ma_reg6]  =  data_tmp_5;
             A_B1R6[bn1_ma_reg7]  =  data_tmp_6;
             A_B0R6[bn0_ma_reg7]  =  data_tmp_7;
             A_B0R6[bn0_ma_reg8]  =  data_tmp_8;
             A_B1R6[bn1_ma_reg8]  =  data_tmp_9;
             //------------------------------------------------------------------
			 data_tmp_1  = A_B0R8[bn0_ma_reg4];
			 data_tmp_2  = A_B0R9[bn0_ma_reg4];
			 data_tmp_3  = A_B0R10[bn0_ma_reg4];
			 data_tmp_4  = A_B0R11[bn0_ma_reg4];
			 data_tmp_5  = A_B0R12[bn0_ma_reg4];
			 data_tmp_6  = A_B0R13[bn0_ma_reg4];
			 data_tmp_7  = A_B0R14[bn0_ma_reg4];
			 data_tmp_8  = A_B0R15[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg4]  =  A_B0R7[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg4]  =  A_B1R7[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg4] =  A_B1R7[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg4] =  A_B0R7[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg4] =  A_B1R7[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg4] =  A_B0R7[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg4] =  A_B0R7[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg4] =  A_B1R7[bn1_ma_reg8];
			 A_B0R7[bn0_ma_reg5]  =  data_tmp_1;
			 A_B1R7[bn1_ma_reg5]  =  data_tmp_2;
			 A_B1R7[bn1_ma_reg6]  =  data_tmp_3;
			 A_B0R7[bn0_ma_reg6]  =  data_tmp_4;
			 A_B1R7[bn1_ma_reg7]  =  data_tmp_5;
			 A_B0R7[bn0_ma_reg7]  =  data_tmp_6;
			 A_B0R7[bn0_ma_reg8]  =  data_tmp_7;
			 A_B1R7[bn1_ma_reg8]  =  data_tmp_8;
			 //----------------------------------------------------------------
			 data_tmp_1  = A_B0R9[bn0_ma_reg5];
			 data_tmp_2  = A_B0R10[bn0_ma_reg5];
			 data_tmp_3  = A_B0R11[bn0_ma_reg5];
			 data_tmp_4  = A_B0R12[bn0_ma_reg5];
			 data_tmp_5  = A_B0R13[bn0_ma_reg5];
			 data_tmp_6  = A_B0R14[bn0_ma_reg5];
			 data_tmp_7  = A_B0R15[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg5]  =  A_B1R8[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg8];
			 A_B1R8[bn1_ma_reg5]  =  data_tmp_1;
			 A_B1R8[bn1_ma_reg6]  =  data_tmp_2;
			 A_B0R8[bn0_ma_reg6]  =  data_tmp_3;
			 A_B1R8[bn1_ma_reg7]  =  data_tmp_4;
			 A_B0R8[bn0_ma_reg7]  =  data_tmp_5;
			 A_B0R8[bn0_ma_reg8]  =  data_tmp_6;
			 A_B1R8[bn1_ma_reg8]  =  data_tmp_7;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B1R10[bn1_ma_reg5];
			 data_tmp_2 = A_B1R11[bn1_ma_reg5];
			 data_tmp_3 = A_B1R12[bn1_ma_reg5];
			 data_tmp_4 = A_B1R13[bn1_ma_reg5];
			 data_tmp_5 = A_B1R14[bn1_ma_reg5];
			 data_tmp_6 = A_B1R15[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg5] = A_B1R9[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg5] = A_B0R9[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg5] = A_B1R9[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg5] = A_B0R9[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg5] = A_B0R9[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg5] = A_B1R9[bn1_ma_reg8];
			 A_B1R9[bn1_ma_reg6]  = data_tmp_1;
			 A_B0R9[bn0_ma_reg6]  = data_tmp_2;
			 A_B1R9[bn1_ma_reg7]  = data_tmp_3;
			 A_B0R9[bn0_ma_reg7]  = data_tmp_4;
			 A_B0R9[bn0_ma_reg8]  = data_tmp_5;
			 A_B1R9[bn1_ma_reg8]  = data_tmp_6;
			 //----------------------------------------------------------------
			 data_tmp_1 =  A_B1R11[bn1_ma_reg6];
			 data_tmp_2 =  A_B1R12[bn1_ma_reg6];
			 data_tmp_3 =  A_B1R13[bn1_ma_reg6];
			 data_tmp_4 =  A_B1R14[bn1_ma_reg6];
			 data_tmp_5 =  A_B1R15[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg6]  = A_B0R10[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg6]  = A_B1R10[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg6]  = A_B0R10[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg6]  = A_B0R10[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg6]  = A_B1R10[bn1_ma_reg8];
			 A_B0R10[bn0_ma_reg6]  = data_tmp_1;
			 A_B1R10[bn1_ma_reg7]  = data_tmp_2;
			 A_B0R10[bn0_ma_reg7]  = data_tmp_3;
			 A_B0R10[bn0_ma_reg8]  = data_tmp_4;
			 A_B1R10[bn1_ma_reg8]  = data_tmp_5;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B0R12[bn0_ma_reg6];
			 data_tmp_2 = A_B0R13[bn0_ma_reg6];
			 data_tmp_3 = A_B0R14[bn0_ma_reg6];
			 data_tmp_4 = A_B0R15[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg6] = A_B1R11[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg6] = A_B0R11[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg6] = A_B0R11[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg6] = A_B1R11[bn1_ma_reg8];
			 A_B1R11[bn1_ma_reg7] = data_tmp_1;
			 A_B0R11[bn0_ma_reg7] = data_tmp_2;
			 A_B0R11[bn0_ma_reg8] = data_tmp_3;
			 A_B1R11[bn1_ma_reg8] = data_tmp_4;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B1R13[bn1_ma_reg7];
			 data_tmp_2 = A_B1R14[bn1_ma_reg7];
			 data_tmp_3 = A_B1R15[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg7] = A_B0R12[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg7] = A_B0R12[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg7] = A_B1R12[bn1_ma_reg8];
			 A_B0R12[bn0_ma_reg7] = data_tmp_1;
			 A_B0R12[bn0_ma_reg8] = data_tmp_2;
			 A_B1R12[bn1_ma_reg8] = data_tmp_3;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B0R14[bn0_ma_reg7];
			 data_tmp_2 = A_B0R15[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg7] = A_B0R13[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg7] = A_B1R13[bn1_ma_reg8];
			 A_B0R13[bn0_ma_reg8] = data_tmp_1;
			 A_B1R13[bn1_ma_reg8] = data_tmp_2;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B0R15[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg8] = A_B1R14[bn1_ma_reg8];
			 A_B1R14[bn1_ma_reg8] = data_tmp_1;
			 //-----------------------------------------------------------------
		  }
		 }
		 else{
		  //radix-16 final stage data relocation
		  if(bn1_bc_tmp > bn0_bc_tmp){
			//bn: 0 1 1 0
		    data_tmp_1 = A_B0R1[bn0_ma_reg1];
		    data_tmp_2 = A_B0R2[bn0_ma_reg1];
		    data_tmp_3 = A_B0R3[bn0_ma_reg1];
			A_B0R1[bn0_ma_reg1] = A_B1R0[bn1_ma_reg1];
			A_B0R2[bn0_ma_reg1] = A_B1R0[bn1_ma_reg2];
			A_B0R3[bn0_ma_reg1] = A_B0R0[bn0_ma_reg2];
			A_B1R0[bn1_ma_reg1] = data_tmp_1;
			A_B1R0[bn1_ma_reg2] = data_tmp_2;
			A_B0R0[bn0_ma_reg2] = data_tmp_3;
			data_tmp_1 = A_B1R2[bn1_ma_reg1];
			data_tmp_2 = A_B1R3[bn1_ma_reg1];
			A_B1R2[bn1_ma_reg1] = A_B1R1[bn1_ma_reg2];
			A_B1R3[bn1_ma_reg1] = A_B0R1[bn0_ma_reg2];
			A_B1R1[bn1_ma_reg2] = data_tmp_1;
			A_B0R1[bn0_ma_reg2] = data_tmp_2;
			data_tmp_1 = A_B1R3[bn1_ma_reg2];
			A_B1R3[bn1_ma_reg2] = A_B0R2[bn0_ma_reg2];
			A_B0R2[bn0_ma_reg2] = data_tmp_1;
			//-----------------------------------------
			data_tmp_1 = A_B0R5[bn0_ma_reg1];
			data_tmp_2 = A_B0R6[bn0_ma_reg1];
			data_tmp_3 = A_B0R7[bn0_ma_reg1];
			A_B0R5[bn0_ma_reg1] = A_B1R4[bn1_ma_reg1];
			A_B0R6[bn0_ma_reg1] = A_B1R4[bn1_ma_reg2];
			A_B0R7[bn0_ma_reg1] = A_B0R4[bn0_ma_reg2];
			A_B1R4[bn1_ma_reg1] = data_tmp_1;
			A_B1R4[bn1_ma_reg2] = data_tmp_2;
			A_B0R4[bn0_ma_reg2] = data_tmp_3;
			data_tmp_1 = A_B1R6[bn1_ma_reg1];
			data_tmp_2 = A_B1R7[bn1_ma_reg1];
			A_B1R6[bn1_ma_reg1] = A_B1R5[bn1_ma_reg2];
			A_B1R7[bn1_ma_reg1] = A_B0R5[bn0_ma_reg2];
			A_B1R5[bn1_ma_reg2] = data_tmp_1;
			A_B0R5[bn0_ma_reg2] = data_tmp_2;
			data_tmp_1 = A_B1R7[bn1_ma_reg2];
			A_B1R7[bn1_ma_reg2] = A_B0R6[bn0_ma_reg2];
			A_B0R6[bn0_ma_reg2] = data_tmp_1;
			//----------------------------------------
			data_tmp_1 = A_B0R9[bn0_ma_reg1];
			data_tmp_2 = A_B0R10[bn0_ma_reg1];
			data_tmp_3 = A_B0R11[bn0_ma_reg1];
			A_B0R9[bn0_ma_reg1]  = A_B1R8[bn1_ma_reg1];
			A_B0R10[bn0_ma_reg1] = A_B1R8[bn1_ma_reg2];
			A_B0R11[bn0_ma_reg1] = A_B0R8[bn0_ma_reg2];
			A_B1R8[bn1_ma_reg1] = data_tmp_1;
			A_B1R8[bn1_ma_reg2] = data_tmp_2;
			A_B0R8[bn0_ma_reg2] = data_tmp_3;
			data_tmp_1 = A_B1R10[bn1_ma_reg1];
			data_tmp_2 = A_B1R11[bn1_ma_reg1];
			A_B1R10[bn1_ma_reg1] = A_B1R9[bn1_ma_reg2];
			A_B1R11[bn1_ma_reg1] = A_B0R9[bn0_ma_reg2];
			A_B1R9[bn1_ma_reg2] = data_tmp_1;
			A_B0R9[bn0_ma_reg2] = data_tmp_2;
			data_tmp_1 = A_B1R11[bn1_ma_reg2];
			A_B1R11[bn1_ma_reg2] = A_B0R10[bn0_ma_reg2];
			A_B0R10[bn0_ma_reg2] = data_tmp_1;			
			//----------------------------------------
			data_tmp_1 = A_B0R13[bn0_ma_reg1];
			data_tmp_2 = A_B0R14[bn0_ma_reg1];
			data_tmp_3 = A_B0R15[bn0_ma_reg1];
			A_B0R13[bn0_ma_reg1] = A_B1R12[bn1_ma_reg1];
			A_B0R14[bn0_ma_reg1] = A_B1R12[bn1_ma_reg2];
			A_B0R15[bn0_ma_reg1] = A_B0R12[bn0_ma_reg2];
			A_B1R12[bn1_ma_reg1] = data_tmp_1;
			A_B1R12[bn1_ma_reg2] = data_tmp_2;
			A_B0R12[bn0_ma_reg2] = data_tmp_3;
			data_tmp_1 = A_B1R14[bn1_ma_reg1];
			data_tmp_2 = A_B1R15[bn1_ma_reg1];
			A_B1R14[bn1_ma_reg1] = A_B1R13[bn1_ma_reg2];
			A_B1R15[bn1_ma_reg1] = A_B0R13[bn0_ma_reg2];
			A_B1R13[bn1_ma_reg2] = data_tmp_1;
			A_B0R13[bn0_ma_reg2] = data_tmp_2;
			data_tmp_1 = A_B1R15[bn1_ma_reg2];
			A_B1R15[bn1_ma_reg2] = A_B0R14[bn0_ma_reg2];
			A_B0R14[bn0_ma_reg2] = data_tmp_1;				
			//******************************************************
			//bn 1 0 0 1
		    data_tmp_1 = A_B1R1[bn1_ma_reg3];
		    data_tmp_2 = A_B1R2[bn1_ma_reg3];
		    data_tmp_3 = A_B1R3[bn1_ma_reg3];
			A_B1R1[bn1_ma_reg3] = A_B0R0[bn0_ma_reg3];
			A_B1R2[bn1_ma_reg3] = A_B0R0[bn0_ma_reg4];
			A_B1R3[bn1_ma_reg3] = A_B1R0[bn1_ma_reg4];
			A_B0R0[bn0_ma_reg3] = data_tmp_1;
			A_B0R0[bn0_ma_reg4] = data_tmp_2;
			A_B1R0[bn1_ma_reg4] = data_tmp_3;
			data_tmp_1 = A_B0R2[bn0_ma_reg3];
			data_tmp_2 = A_B0R3[bn0_ma_reg3];
			A_B0R2[bn0_ma_reg3] = A_B0R1[bn0_ma_reg4];
			A_B0R3[bn0_ma_reg3] = A_B1R1[bn1_ma_reg4];
			A_B0R1[bn0_ma_reg4] = data_tmp_1;
			A_B1R1[bn1_ma_reg4] = data_tmp_2;
			data_tmp_1 = A_B0R3[bn0_ma_reg4];
			A_B0R3[bn0_ma_reg4] = A_B1R2[bn1_ma_reg4];
			A_B1R2[bn1_ma_reg4] = data_tmp_1;
            //-----------------------------------------------------
		    data_tmp_1 = A_B1R5[bn1_ma_reg3];
		    data_tmp_2 = A_B1R6[bn1_ma_reg3];
		    data_tmp_3 = A_B1R7[bn1_ma_reg3];
			A_B1R5[bn1_ma_reg3] = A_B0R4[bn0_ma_reg3];
			A_B1R6[bn1_ma_reg3] = A_B0R4[bn0_ma_reg4];
			A_B1R7[bn1_ma_reg3] = A_B1R4[bn1_ma_reg4];
			A_B0R4[bn0_ma_reg3] = data_tmp_1;
			A_B0R4[bn0_ma_reg4] = data_tmp_2;
			A_B1R4[bn1_ma_reg4] = data_tmp_3;
			data_tmp_1 = A_B0R6[bn0_ma_reg3];
			data_tmp_2 = A_B0R7[bn0_ma_reg3];
			A_B0R6[bn0_ma_reg3] = A_B0R5[bn0_ma_reg4];
			A_B0R7[bn0_ma_reg3] = A_B1R5[bn1_ma_reg4];
			A_B0R5[bn0_ma_reg4] = data_tmp_1;
			A_B1R5[bn1_ma_reg4] = data_tmp_2;
			data_tmp_1 = A_B0R7[bn0_ma_reg4];
			A_B0R7[bn0_ma_reg4] = A_B1R6[bn1_ma_reg4];
			A_B1R6[bn1_ma_reg4] = data_tmp_1;
            //-----------------------------------------------------
		    data_tmp_1 = A_B1R9[bn1_ma_reg3];
		    data_tmp_2 = A_B1R10[bn1_ma_reg3];
		    data_tmp_3 = A_B1R11[bn1_ma_reg3];
			A_B1R9[bn1_ma_reg3]  = A_B0R8[bn0_ma_reg3];
			A_B1R10[bn1_ma_reg3] = A_B0R8[bn0_ma_reg4];
			A_B1R11[bn1_ma_reg3] = A_B1R8[bn1_ma_reg4];
			A_B0R8[bn0_ma_reg3] = data_tmp_1;
			A_B0R8[bn0_ma_reg4] = data_tmp_2;
			A_B1R8[bn1_ma_reg4] = data_tmp_3;
			data_tmp_1 = A_B0R10[bn0_ma_reg3];
			data_tmp_2 = A_B0R11[bn0_ma_reg3];
			A_B0R10[bn0_ma_reg3] = A_B0R9[bn0_ma_reg4];
			A_B0R11[bn0_ma_reg3] = A_B1R9[bn1_ma_reg4];
			A_B0R9[bn0_ma_reg4] = data_tmp_1;
			A_B1R9[bn1_ma_reg4] = data_tmp_2;
			data_tmp_1 = A_B0R11[bn0_ma_reg4];
			A_B0R11[bn0_ma_reg4] = A_B1R10[bn1_ma_reg4];
			A_B1R10[bn1_ma_reg4] = data_tmp_1;
            //-----------------------------------------------------
		    data_tmp_1 = A_B1R13[bn1_ma_reg3];
		    data_tmp_2 = A_B1R14[bn1_ma_reg3];
		    data_tmp_3 = A_B1R15[bn1_ma_reg3];
			A_B1R13[bn1_ma_reg3] = A_B0R12[bn0_ma_reg3];
			A_B1R14[bn1_ma_reg3] = A_B0R12[bn0_ma_reg4];
			A_B1R15[bn1_ma_reg3] = A_B1R12[bn1_ma_reg4];
			A_B0R12[bn0_ma_reg3] = data_tmp_1;
			A_B0R12[bn0_ma_reg4] = data_tmp_2;
			A_B1R12[bn1_ma_reg4] = data_tmp_3;
			data_tmp_1 = A_B0R14[bn0_ma_reg3];
			data_tmp_2 = A_B0R15[bn0_ma_reg3];
			A_B0R14[bn0_ma_reg3] = A_B0R13[bn0_ma_reg4];
			A_B0R15[bn0_ma_reg3] = A_B1R13[bn1_ma_reg4];
			A_B0R13[bn0_ma_reg4] = data_tmp_1;
			A_B1R13[bn1_ma_reg4] = data_tmp_2;
			data_tmp_1 = A_B0R15[bn0_ma_reg4];
			A_B0R15[bn0_ma_reg4] = A_B1R14[bn1_ma_reg4];
			A_B1R14[bn1_ma_reg4] = data_tmp_1;
            //*********************************************************	
			//bn 1 0 0 1
		    data_tmp_1 = A_B1R1[bn1_ma_reg5];
		    data_tmp_2 = A_B1R2[bn1_ma_reg5];
		    data_tmp_3 = A_B1R3[bn1_ma_reg5];
			A_B1R1[bn1_ma_reg5] = A_B0R0[bn0_ma_reg5];
			A_B1R2[bn1_ma_reg5] = A_B0R0[bn0_ma_reg6];
			A_B1R3[bn1_ma_reg5] = A_B1R0[bn1_ma_reg6];
			A_B0R0[bn0_ma_reg5] = data_tmp_1;
			A_B0R0[bn0_ma_reg6] = data_tmp_2;
			A_B1R0[bn1_ma_reg6] = data_tmp_3;
			data_tmp_1 = A_B0R2[bn0_ma_reg5];
			data_tmp_2 = A_B0R3[bn0_ma_reg5];
			A_B0R2[bn0_ma_reg5] = A_B0R1[bn0_ma_reg6];
			A_B0R3[bn0_ma_reg5] = A_B1R1[bn1_ma_reg6];
			A_B0R1[bn0_ma_reg6] = data_tmp_1;
			A_B1R1[bn1_ma_reg6] = data_tmp_2;
			data_tmp_1 = A_B0R3[bn0_ma_reg6];
			A_B0R3[bn0_ma_reg6] = A_B1R2[bn1_ma_reg6];
			A_B1R2[bn1_ma_reg6] = data_tmp_1;
            //-----------------------------------------------------
		    data_tmp_1 = A_B1R5[bn1_ma_reg5];
		    data_tmp_2 = A_B1R6[bn1_ma_reg5];
		    data_tmp_3 = A_B1R7[bn1_ma_reg5];
			A_B1R5[bn1_ma_reg5] = A_B0R4[bn0_ma_reg5];
			A_B1R6[bn1_ma_reg5] = A_B0R4[bn0_ma_reg6];
			A_B1R7[bn1_ma_reg5] = A_B1R4[bn1_ma_reg6];
			A_B0R4[bn0_ma_reg5] = data_tmp_1;
			A_B0R4[bn0_ma_reg6] = data_tmp_2;
			A_B1R4[bn1_ma_reg6] = data_tmp_3;
			data_tmp_1 = A_B0R6[bn0_ma_reg5];
			data_tmp_2 = A_B0R7[bn0_ma_reg5];
			A_B0R6[bn0_ma_reg5] = A_B0R5[bn0_ma_reg6];
			A_B0R7[bn0_ma_reg5] = A_B1R5[bn1_ma_reg6];
			A_B0R5[bn0_ma_reg6] = data_tmp_1;
			A_B1R5[bn1_ma_reg6] = data_tmp_2;
			data_tmp_1 = A_B0R7[bn0_ma_reg6];
			A_B0R7[bn0_ma_reg6] = A_B1R6[bn1_ma_reg6];
			A_B1R6[bn1_ma_reg6] = data_tmp_1;
            //-----------------------------------------------------
		    data_tmp_1 = A_B1R9[bn1_ma_reg5];
		    data_tmp_2 = A_B1R10[bn1_ma_reg5];
		    data_tmp_3 = A_B1R11[bn1_ma_reg5];
			A_B1R9[bn1_ma_reg5]  = A_B0R8[bn0_ma_reg5];
			A_B1R10[bn1_ma_reg5] = A_B0R8[bn0_ma_reg6];
			A_B1R11[bn1_ma_reg5] = A_B1R8[bn1_ma_reg6];
			A_B0R8[bn0_ma_reg5] = data_tmp_1;
			A_B0R8[bn0_ma_reg6] = data_tmp_2;
			A_B1R8[bn1_ma_reg6] = data_tmp_3;
			data_tmp_1 = A_B0R10[bn0_ma_reg5];
			data_tmp_2 = A_B0R11[bn0_ma_reg5];
			A_B0R10[bn0_ma_reg5] = A_B0R9[bn0_ma_reg6];
			A_B0R11[bn0_ma_reg5] = A_B1R9[bn1_ma_reg6];
			A_B0R9[bn0_ma_reg6] = data_tmp_1;
			A_B1R9[bn1_ma_reg6] = data_tmp_2;
			data_tmp_1 = A_B0R11[bn0_ma_reg6];
			A_B0R11[bn0_ma_reg6] = A_B1R10[bn1_ma_reg6];
			A_B1R10[bn1_ma_reg6] = data_tmp_1;
            //-----------------------------------------------------
		    data_tmp_1 = A_B1R13[bn1_ma_reg5];
		    data_tmp_2 = A_B1R14[bn1_ma_reg5];
		    data_tmp_3 = A_B1R15[bn1_ma_reg5];
			A_B1R13[bn1_ma_reg5] = A_B0R12[bn0_ma_reg5];
			A_B1R14[bn1_ma_reg5] = A_B0R12[bn0_ma_reg6];
			A_B1R15[bn1_ma_reg5] = A_B1R12[bn1_ma_reg6];
			A_B0R12[bn0_ma_reg5] = data_tmp_1;
			A_B0R12[bn0_ma_reg6] = data_tmp_2;
			A_B1R12[bn1_ma_reg6] = data_tmp_3;
			data_tmp_1 = A_B0R14[bn0_ma_reg5];
			data_tmp_2 = A_B0R15[bn0_ma_reg5];
			A_B0R14[bn0_ma_reg5] = A_B0R13[bn0_ma_reg6];
			A_B0R15[bn0_ma_reg5] = A_B1R13[bn1_ma_reg6];
			A_B0R13[bn0_ma_reg6] = data_tmp_1;
			A_B1R13[bn1_ma_reg6] = data_tmp_2;
			data_tmp_1 = A_B0R15[bn0_ma_reg6];
			A_B0R15[bn0_ma_reg6] = A_B1R14[bn1_ma_reg6];
			A_B1R14[bn1_ma_reg6] = data_tmp_1;
            //****************************************
            //bn: 0 1 1 0
		    data_tmp_1 = A_B0R1[bn0_ma_reg7];
		    data_tmp_2 = A_B0R2[bn0_ma_reg7];
		    data_tmp_3 = A_B0R3[bn0_ma_reg7];
			A_B0R1[bn0_ma_reg7] = A_B1R0[bn1_ma_reg7];
			A_B0R2[bn0_ma_reg7] = A_B1R0[bn1_ma_reg8];
			A_B0R3[bn0_ma_reg7] = A_B0R0[bn0_ma_reg8];
			A_B1R0[bn1_ma_reg7] = data_tmp_1;
			A_B1R0[bn1_ma_reg8] = data_tmp_2;
			A_B0R0[bn0_ma_reg8] = data_tmp_3;
			data_tmp_1 = A_B1R2[bn1_ma_reg7];
			data_tmp_2 = A_B1R3[bn1_ma_reg7];
			A_B1R2[bn1_ma_reg7] = A_B1R1[bn1_ma_reg8];
			A_B1R3[bn1_ma_reg7] = A_B0R1[bn0_ma_reg8];
			A_B1R1[bn1_ma_reg8] = data_tmp_1;
			A_B0R1[bn0_ma_reg8] = data_tmp_2;
			data_tmp_1 = A_B1R3[bn1_ma_reg8];
			A_B1R3[bn1_ma_reg8] = A_B0R2[bn0_ma_reg8];
			A_B0R2[bn0_ma_reg8] = data_tmp_1;
			//-----------------------------------------
			data_tmp_1 = A_B0R5[bn0_ma_reg7];
			data_tmp_2 = A_B0R6[bn0_ma_reg7];
			data_tmp_3 = A_B0R7[bn0_ma_reg7];
			A_B0R5[bn0_ma_reg7] = A_B1R4[bn1_ma_reg7];
			A_B0R6[bn0_ma_reg7] = A_B1R4[bn1_ma_reg8];
			A_B0R7[bn0_ma_reg7] = A_B0R4[bn0_ma_reg8];
			A_B1R4[bn1_ma_reg7] = data_tmp_1;
			A_B1R4[bn1_ma_reg8] = data_tmp_2;
			A_B0R4[bn0_ma_reg8] = data_tmp_3;
			data_tmp_1 = A_B1R6[bn1_ma_reg7];
			data_tmp_2 = A_B1R7[bn1_ma_reg7];
			A_B1R6[bn1_ma_reg7] = A_B1R5[bn1_ma_reg8];
			A_B1R7[bn1_ma_reg7] = A_B0R5[bn0_ma_reg8];
			A_B1R5[bn1_ma_reg8] = data_tmp_1;
			A_B0R5[bn0_ma_reg8] = data_tmp_2;
			data_tmp_1 = A_B1R7[bn1_ma_reg8];
			A_B1R7[bn1_ma_reg8] = A_B0R6[bn0_ma_reg8];
			A_B0R6[bn0_ma_reg8] = data_tmp_1;
			//----------------------------------------
			data_tmp_1 = A_B0R9[bn0_ma_reg7];
			data_tmp_2 = A_B0R10[bn0_ma_reg7];
			data_tmp_3 = A_B0R11[bn0_ma_reg7];
			A_B0R9[bn0_ma_reg7]  = A_B1R8[bn1_ma_reg7];
			A_B0R10[bn0_ma_reg7] = A_B1R8[bn1_ma_reg8];
			A_B0R11[bn0_ma_reg7] = A_B0R8[bn0_ma_reg8];
			A_B1R8[bn1_ma_reg7] = data_tmp_1;
			A_B1R8[bn1_ma_reg8] = data_tmp_2;
			A_B0R8[bn0_ma_reg8] = data_tmp_3;
			data_tmp_1 = A_B1R10[bn1_ma_reg7];
			data_tmp_2 = A_B1R11[bn1_ma_reg7];
			A_B1R10[bn1_ma_reg7] = A_B1R9[bn1_ma_reg8];
			A_B1R11[bn1_ma_reg7] = A_B0R9[bn0_ma_reg8];
			A_B1R9[bn1_ma_reg8] = data_tmp_1;
			A_B0R9[bn0_ma_reg8] = data_tmp_2;
			data_tmp_1 = A_B1R11[bn1_ma_reg8];
			A_B1R11[bn1_ma_reg8] = A_B0R10[bn0_ma_reg8];
			A_B0R10[bn0_ma_reg8] = data_tmp_1;			
			//----------------------------------------
			data_tmp_1 = A_B0R13[bn0_ma_reg7];
			data_tmp_2 = A_B0R14[bn0_ma_reg7];
			data_tmp_3 = A_B0R15[bn0_ma_reg7];
			A_B0R13[bn0_ma_reg7] = A_B1R12[bn1_ma_reg7];
			A_B0R14[bn0_ma_reg7] = A_B1R12[bn1_ma_reg8];
			A_B0R15[bn0_ma_reg7] = A_B0R12[bn0_ma_reg8];
			A_B1R12[bn1_ma_reg7] = data_tmp_1;
			A_B1R12[bn1_ma_reg8] = data_tmp_2;
			A_B0R12[bn0_ma_reg8] = data_tmp_3;
			data_tmp_1 = A_B1R14[bn1_ma_reg7];
			data_tmp_2 = A_B1R15[bn1_ma_reg7];
			A_B1R14[bn1_ma_reg7] = A_B1R13[bn1_ma_reg8];
			A_B1R15[bn1_ma_reg7] = A_B0R13[bn0_ma_reg8];
			A_B1R13[bn1_ma_reg8] = data_tmp_1;
			A_B0R13[bn0_ma_reg8] = data_tmp_2;
			data_tmp_1 = A_B1R15[bn1_ma_reg8];
			A_B1R15[bn1_ma_reg8] = A_B0R14[bn0_ma_reg8];
			A_B0R14[bn0_ma_reg8] = data_tmp_1;		
            //------------------------------------------------------			
		  }else{	 
			//bn 1 0 0 1
		    data_tmp_1 = A_B1R1[bn1_ma_reg1];
		    data_tmp_2 = A_B1R2[bn1_ma_reg1];
		    data_tmp_3 = A_B1R3[bn1_ma_reg1];
			A_B1R1[bn1_ma_reg1] = A_B0R0[bn0_ma_reg1];
			A_B1R2[bn1_ma_reg1] = A_B0R0[bn0_ma_reg2];
			A_B1R3[bn1_ma_reg1] = A_B1R0[bn1_ma_reg2];
			A_B0R0[bn0_ma_reg1] = data_tmp_1;
			A_B0R0[bn0_ma_reg2] = data_tmp_2;
			A_B1R0[bn1_ma_reg2] = data_tmp_3;
			data_tmp_1 = A_B0R2[bn0_ma_reg1];
			data_tmp_2 = A_B0R3[bn0_ma_reg1];
			A_B0R2[bn0_ma_reg1] = A_B0R1[bn0_ma_reg2];
			A_B0R3[bn0_ma_reg1] = A_B1R1[bn1_ma_reg2];
			A_B0R1[bn0_ma_reg2] = data_tmp_1;
			A_B1R1[bn1_ma_reg2] = data_tmp_2;
			data_tmp_1 = A_B0R3[bn0_ma_reg2];
			A_B0R3[bn0_ma_reg2] = A_B1R2[bn1_ma_reg2];
			A_B1R2[bn1_ma_reg2] = data_tmp_1;
            //-----------------------------------------------------
		    data_tmp_1 = A_B1R5[bn1_ma_reg1];
		    data_tmp_2 = A_B1R6[bn1_ma_reg1];
		    data_tmp_3 = A_B1R7[bn1_ma_reg1];
			A_B1R5[bn1_ma_reg1] = A_B0R4[bn0_ma_reg1];
			A_B1R6[bn1_ma_reg1] = A_B0R4[bn0_ma_reg2];
			A_B1R7[bn1_ma_reg1] = A_B1R4[bn1_ma_reg2];
			A_B0R4[bn0_ma_reg1] = data_tmp_1;
			A_B0R4[bn0_ma_reg2] = data_tmp_2;
			A_B1R4[bn1_ma_reg2] = data_tmp_3;
			data_tmp_1 = A_B0R6[bn0_ma_reg1];
			data_tmp_2 = A_B0R7[bn0_ma_reg1];
			A_B0R6[bn0_ma_reg1] = A_B0R5[bn0_ma_reg2];
			A_B0R7[bn0_ma_reg1] = A_B1R5[bn1_ma_reg2];
			A_B0R5[bn0_ma_reg2] = data_tmp_1;
			A_B1R5[bn1_ma_reg2] = data_tmp_2;
			data_tmp_1 = A_B0R7[bn0_ma_reg2];
			A_B0R7[bn0_ma_reg2] = A_B1R6[bn1_ma_reg2];
			A_B1R6[bn1_ma_reg2] = data_tmp_1;
            //-----------------------------------------------------
		    data_tmp_1 = A_B1R9[bn1_ma_reg1];
		    data_tmp_2 = A_B1R10[bn1_ma_reg1];
		    data_tmp_3 = A_B1R11[bn1_ma_reg1];
			A_B1R9[bn1_ma_reg1]  = A_B0R8[bn0_ma_reg1];
			A_B1R10[bn1_ma_reg1] = A_B0R8[bn0_ma_reg2];
			A_B1R11[bn1_ma_reg1] = A_B1R8[bn1_ma_reg2];
			A_B0R8[bn0_ma_reg1] = data_tmp_1;
			A_B0R8[bn0_ma_reg2] = data_tmp_2;
			A_B1R8[bn1_ma_reg2] = data_tmp_3;
			data_tmp_1 = A_B0R10[bn0_ma_reg1];
			data_tmp_2 = A_B0R11[bn0_ma_reg1];
			A_B0R10[bn0_ma_reg1] = A_B0R9[bn0_ma_reg2];
			A_B0R11[bn0_ma_reg1] = A_B1R9[bn1_ma_reg2];
			A_B0R9[bn0_ma_reg2] = data_tmp_1;
			A_B1R9[bn1_ma_reg2] = data_tmp_2;
			data_tmp_1 = A_B0R11[bn0_ma_reg2];
			A_B0R11[bn0_ma_reg2] = A_B1R10[bn1_ma_reg2];
			A_B1R10[bn1_ma_reg2] = data_tmp_1;
            //-----------------------------------------------------
		    data_tmp_1 = A_B1R13[bn1_ma_reg1];
		    data_tmp_2 = A_B1R14[bn1_ma_reg1];
		    data_tmp_3 = A_B1R15[bn1_ma_reg1];
			A_B1R13[bn1_ma_reg1] = A_B0R12[bn0_ma_reg1];
			A_B1R14[bn1_ma_reg1] = A_B0R12[bn0_ma_reg2];
			A_B1R15[bn1_ma_reg1] = A_B1R12[bn1_ma_reg2];
			A_B0R12[bn0_ma_reg1] = data_tmp_1;
			A_B0R12[bn0_ma_reg2] = data_tmp_2;
			A_B1R12[bn1_ma_reg2] = data_tmp_3;
			data_tmp_1 = A_B0R14[bn0_ma_reg1];
			data_tmp_2 = A_B0R15[bn0_ma_reg1];
			A_B0R14[bn0_ma_reg1] = A_B0R13[bn0_ma_reg2];
			A_B0R15[bn0_ma_reg1] = A_B1R13[bn1_ma_reg2];
			A_B0R13[bn0_ma_reg2] = data_tmp_1;
			A_B1R13[bn1_ma_reg2] = data_tmp_2;
			data_tmp_1 = A_B0R15[bn0_ma_reg2];
			A_B0R15[bn0_ma_reg2] = A_B1R14[bn1_ma_reg2];
			A_B1R14[bn1_ma_reg2] = data_tmp_1;
            //******************************************************
            //bn: 0 1 1 0
		    data_tmp_1 = A_B0R1[bn0_ma_reg3];
		    data_tmp_2 = A_B0R2[bn0_ma_reg3];
		    data_tmp_3 = A_B0R3[bn0_ma_reg3];
			A_B0R1[bn0_ma_reg3] = A_B1R0[bn1_ma_reg3];
			A_B0R2[bn0_ma_reg3] = A_B1R0[bn1_ma_reg4];
			A_B0R3[bn0_ma_reg3] = A_B0R0[bn0_ma_reg4];
			A_B1R0[bn1_ma_reg3] = data_tmp_1;
			A_B1R0[bn1_ma_reg4] = data_tmp_2;
			A_B0R0[bn0_ma_reg4] = data_tmp_3;
			data_tmp_1 = A_B1R2[bn1_ma_reg3];
			data_tmp_2 = A_B1R3[bn1_ma_reg3];
			A_B1R2[bn1_ma_reg3] = A_B1R1[bn1_ma_reg4];
			A_B1R3[bn1_ma_reg3] = A_B0R1[bn0_ma_reg4];
			A_B1R1[bn1_ma_reg4] = data_tmp_1;
			A_B0R1[bn0_ma_reg4] = data_tmp_2;
			data_tmp_1 = A_B1R3[bn1_ma_reg4];
			A_B1R3[bn1_ma_reg4] = A_B0R2[bn0_ma_reg4];
			A_B0R2[bn0_ma_reg4] = data_tmp_1;
			//-----------------------------------------
			data_tmp_1 = A_B0R5[bn0_ma_reg3];
			data_tmp_2 = A_B0R6[bn0_ma_reg3];
			data_tmp_3 = A_B0R7[bn0_ma_reg3];
			A_B0R5[bn0_ma_reg3] = A_B1R4[bn1_ma_reg3];
			A_B0R6[bn0_ma_reg3] = A_B1R4[bn1_ma_reg4];
			A_B0R7[bn0_ma_reg3] = A_B0R4[bn0_ma_reg4];
			A_B1R4[bn1_ma_reg3] = data_tmp_1;
			A_B1R4[bn1_ma_reg4] = data_tmp_2;
			A_B0R4[bn0_ma_reg4] = data_tmp_3;
			data_tmp_1 = A_B1R6[bn1_ma_reg3];
			data_tmp_2 = A_B1R7[bn1_ma_reg3];
			A_B1R6[bn1_ma_reg3] = A_B1R5[bn1_ma_reg4];
			A_B1R7[bn1_ma_reg3] = A_B0R5[bn0_ma_reg4];
			A_B1R5[bn1_ma_reg4] = data_tmp_1;
			A_B0R5[bn0_ma_reg4] = data_tmp_2;
			data_tmp_1 = A_B1R7[bn1_ma_reg4];
			A_B1R7[bn1_ma_reg4] = A_B0R6[bn0_ma_reg4];
			A_B0R6[bn0_ma_reg4] = data_tmp_1;
			//----------------------------------------
			data_tmp_1 = A_B0R9[bn0_ma_reg3];
			data_tmp_2 = A_B0R10[bn0_ma_reg3];
			data_tmp_3 = A_B0R11[bn0_ma_reg3];
			A_B0R9[bn0_ma_reg3]  = A_B1R8[bn1_ma_reg3];
			A_B0R10[bn0_ma_reg3] = A_B1R8[bn1_ma_reg4];
			A_B0R11[bn0_ma_reg3] = A_B0R8[bn0_ma_reg4];
			A_B1R8[bn1_ma_reg3] = data_tmp_1;
			A_B1R8[bn1_ma_reg4] = data_tmp_2;
			A_B0R8[bn0_ma_reg4] = data_tmp_3;
			data_tmp_1 = A_B1R10[bn1_ma_reg3];
			data_tmp_2 = A_B1R11[bn1_ma_reg3];
			A_B1R10[bn1_ma_reg3] = A_B1R9[bn1_ma_reg4];
			A_B1R11[bn1_ma_reg3] = A_B0R9[bn0_ma_reg4];
			A_B1R9[bn1_ma_reg4] = data_tmp_1;
			A_B0R9[bn0_ma_reg4] = data_tmp_2;
			data_tmp_1 = A_B1R11[bn1_ma_reg4];
			A_B1R11[bn1_ma_reg4] = A_B0R10[bn0_ma_reg4];
			A_B0R10[bn0_ma_reg4] = data_tmp_1;			
			//----------------------------------------
			data_tmp_1 = A_B0R13[bn0_ma_reg3];
			data_tmp_2 = A_B0R14[bn0_ma_reg3];
			data_tmp_3 = A_B0R15[bn0_ma_reg3];
			A_B0R13[bn0_ma_reg3] = A_B1R12[bn1_ma_reg3];
			A_B0R14[bn0_ma_reg3] = A_B1R12[bn1_ma_reg4];
			A_B0R15[bn0_ma_reg3] = A_B0R12[bn0_ma_reg4];
			A_B1R12[bn1_ma_reg3] = data_tmp_1;
			A_B1R12[bn1_ma_reg4] = data_tmp_2;
			A_B0R12[bn0_ma_reg4] = data_tmp_3;
			data_tmp_1 = A_B1R14[bn1_ma_reg3];
			data_tmp_2 = A_B1R15[bn1_ma_reg3];
			A_B1R14[bn1_ma_reg3] = A_B1R13[bn1_ma_reg4];
			A_B1R15[bn1_ma_reg3] = A_B0R13[bn0_ma_reg4];
			A_B1R13[bn1_ma_reg4] = data_tmp_1;
			A_B0R13[bn0_ma_reg4] = data_tmp_2;
			data_tmp_1 = A_B1R15[bn1_ma_reg4];
			A_B1R15[bn1_ma_reg4] = A_B0R14[bn0_ma_reg4];
			A_B0R14[bn0_ma_reg4] = data_tmp_1;
            //******************************************************
            //bn: 0 1 1 0
		    data_tmp_1 = A_B0R1[bn0_ma_reg5];
		    data_tmp_2 = A_B0R2[bn0_ma_reg5];
		    data_tmp_3 = A_B0R3[bn0_ma_reg5];
			A_B0R1[bn0_ma_reg5] = A_B1R0[bn1_ma_reg5];
			A_B0R2[bn0_ma_reg5] = A_B1R0[bn1_ma_reg6];
			A_B0R3[bn0_ma_reg5] = A_B0R0[bn0_ma_reg6];
			A_B1R0[bn1_ma_reg5] = data_tmp_1;
			A_B1R0[bn1_ma_reg6] = data_tmp_2;
			A_B0R0[bn0_ma_reg6] = data_tmp_3;
			data_tmp_1 = A_B1R2[bn1_ma_reg5];
			data_tmp_2 = A_B1R3[bn1_ma_reg5];
			A_B1R2[bn1_ma_reg5] = A_B1R1[bn1_ma_reg6];
			A_B1R3[bn1_ma_reg5] = A_B0R1[bn0_ma_reg6];
			A_B1R1[bn1_ma_reg6] = data_tmp_1;
			A_B0R1[bn0_ma_reg6] = data_tmp_2;
			data_tmp_1 = A_B1R3[bn1_ma_reg6];
			A_B1R3[bn1_ma_reg6] = A_B0R2[bn0_ma_reg6];
			A_B0R2[bn0_ma_reg6] = data_tmp_1;
			//-----------------------------------------
			data_tmp_1 = A_B0R5[bn0_ma_reg5];
			data_tmp_2 = A_B0R6[bn0_ma_reg5];
			data_tmp_3 = A_B0R7[bn0_ma_reg5];
			A_B0R5[bn0_ma_reg5] = A_B1R4[bn1_ma_reg5];
			A_B0R6[bn0_ma_reg5] = A_B1R4[bn1_ma_reg6];
			A_B0R7[bn0_ma_reg5] = A_B0R4[bn0_ma_reg6];
			A_B1R4[bn1_ma_reg5] = data_tmp_1;
			A_B1R4[bn1_ma_reg6] = data_tmp_2;
			A_B0R4[bn0_ma_reg6] = data_tmp_3;
			data_tmp_1 = A_B1R6[bn1_ma_reg5];
			data_tmp_2 = A_B1R7[bn1_ma_reg5];
			A_B1R6[bn1_ma_reg5] = A_B1R5[bn1_ma_reg6];
			A_B1R7[bn1_ma_reg5] = A_B0R5[bn0_ma_reg6];
			A_B1R5[bn1_ma_reg6] = data_tmp_1;
			A_B0R5[bn0_ma_reg6] = data_tmp_2;
			data_tmp_1 = A_B1R7[bn1_ma_reg6];
			A_B1R7[bn1_ma_reg6] = A_B0R6[bn0_ma_reg6];
			A_B0R6[bn0_ma_reg6] = data_tmp_1;
			//----------------------------------------
			data_tmp_1 = A_B0R9[bn0_ma_reg5];
			data_tmp_2 = A_B0R10[bn0_ma_reg5];
			data_tmp_3 = A_B0R11[bn0_ma_reg5];
			A_B0R9[bn0_ma_reg5]  = A_B1R8[bn1_ma_reg5];
			A_B0R10[bn0_ma_reg5] = A_B1R8[bn1_ma_reg6];
			A_B0R11[bn0_ma_reg5] = A_B0R8[bn0_ma_reg6];
			A_B1R8[bn1_ma_reg5] = data_tmp_1;
			A_B1R8[bn1_ma_reg6] = data_tmp_2;
			A_B0R8[bn0_ma_reg6] = data_tmp_3;
			data_tmp_1 = A_B1R10[bn1_ma_reg5];
			data_tmp_2 = A_B1R11[bn1_ma_reg5];
			A_B1R10[bn1_ma_reg5] = A_B1R9[bn1_ma_reg6];
			A_B1R11[bn1_ma_reg5] = A_B0R9[bn0_ma_reg6];
			A_B1R9[bn1_ma_reg6] = data_tmp_1;
			A_B0R9[bn0_ma_reg6] = data_tmp_2;
			data_tmp_1 = A_B1R11[bn1_ma_reg6];
			A_B1R11[bn1_ma_reg6] = A_B0R10[bn0_ma_reg6];
			A_B0R10[bn0_ma_reg6] = data_tmp_1;			
			//----------------------------------------
			data_tmp_1 = A_B0R13[bn0_ma_reg5];
			data_tmp_2 = A_B0R14[bn0_ma_reg5];
			data_tmp_3 = A_B0R15[bn0_ma_reg5];
			A_B0R13[bn0_ma_reg5] = A_B1R12[bn1_ma_reg5];
			A_B0R14[bn0_ma_reg5] = A_B1R12[bn1_ma_reg6];
			A_B0R15[bn0_ma_reg5] = A_B0R12[bn0_ma_reg6];
			A_B1R12[bn1_ma_reg5] = data_tmp_1;
			A_B1R12[bn1_ma_reg6] = data_tmp_2;
			A_B0R12[bn0_ma_reg6] = data_tmp_3;
			data_tmp_1 = A_B1R14[bn1_ma_reg5];
			data_tmp_2 = A_B1R15[bn1_ma_reg5];
			A_B1R14[bn1_ma_reg5] = A_B1R13[bn1_ma_reg6];
			A_B1R15[bn1_ma_reg5] = A_B0R13[bn0_ma_reg6];
			A_B1R13[bn1_ma_reg6] = data_tmp_1;
			A_B0R13[bn0_ma_reg6] = data_tmp_2;
			data_tmp_1 = A_B1R15[bn1_ma_reg6];
			A_B1R15[bn1_ma_reg6] = A_B0R14[bn0_ma_reg6];
			A_B0R14[bn0_ma_reg6] = data_tmp_1;			
			//***************************************************
			//bn 1 0 0 1
		    data_tmp_1 = A_B1R1[bn1_ma_reg7];
		    data_tmp_2 = A_B1R2[bn1_ma_reg7];
		    data_tmp_3 = A_B1R3[bn1_ma_reg7];
			A_B1R1[bn1_ma_reg7] = A_B0R0[bn0_ma_reg7];
			A_B1R2[bn1_ma_reg7] = A_B0R0[bn0_ma_reg8];
			A_B1R3[bn1_ma_reg7] = A_B1R0[bn1_ma_reg8];
			A_B0R0[bn0_ma_reg7] = data_tmp_1;
			A_B0R0[bn0_ma_reg8] = data_tmp_2;
			A_B1R0[bn1_ma_reg8] = data_tmp_3;
			data_tmp_1 = A_B0R2[bn0_ma_reg7];
			data_tmp_2 = A_B0R3[bn0_ma_reg7];
			A_B0R2[bn0_ma_reg7] = A_B0R1[bn0_ma_reg8];
			A_B0R3[bn0_ma_reg7] = A_B1R1[bn1_ma_reg8];
			A_B0R1[bn0_ma_reg8] = data_tmp_1;
			A_B1R1[bn1_ma_reg8] = data_tmp_2;
			data_tmp_1 = A_B0R3[bn0_ma_reg8];
			A_B0R3[bn0_ma_reg8] = A_B1R2[bn1_ma_reg8];
			A_B1R2[bn1_ma_reg8] = data_tmp_1;
            //-----------------------------------------------------
		    data_tmp_1 = A_B1R5[bn1_ma_reg7];
		    data_tmp_2 = A_B1R6[bn1_ma_reg7];
		    data_tmp_3 = A_B1R7[bn1_ma_reg7];
			A_B1R5[bn1_ma_reg7] = A_B0R4[bn0_ma_reg7];
			A_B1R6[bn1_ma_reg7] = A_B0R4[bn0_ma_reg8];
			A_B1R7[bn1_ma_reg7] = A_B1R4[bn1_ma_reg8];
			A_B0R4[bn0_ma_reg7] = data_tmp_1;
			A_B0R4[bn0_ma_reg8] = data_tmp_2;
			A_B1R4[bn1_ma_reg8] = data_tmp_3;
			data_tmp_1 = A_B0R6[bn0_ma_reg7];
			data_tmp_2 = A_B0R7[bn0_ma_reg7];
			A_B0R6[bn0_ma_reg7] = A_B0R5[bn0_ma_reg8];
			A_B0R7[bn0_ma_reg7] = A_B1R5[bn1_ma_reg8];
			A_B0R5[bn0_ma_reg8] = data_tmp_1;
			A_B1R5[bn1_ma_reg8] = data_tmp_2;
			data_tmp_1 = A_B0R7[bn0_ma_reg8];
			A_B0R7[bn0_ma_reg8] = A_B1R6[bn1_ma_reg8];
			A_B1R6[bn1_ma_reg8] = data_tmp_1;
            //-----------------------------------------------------
		    data_tmp_1 = A_B1R9[bn1_ma_reg7];
		    data_tmp_2 = A_B1R10[bn1_ma_reg7];
		    data_tmp_3 = A_B1R11[bn1_ma_reg7];
			A_B1R9[bn1_ma_reg7]  = A_B0R8[bn0_ma_reg7];
			A_B1R10[bn1_ma_reg7] = A_B0R8[bn0_ma_reg8];
			A_B1R11[bn1_ma_reg7] = A_B1R8[bn1_ma_reg8];
			A_B0R8[bn0_ma_reg7] = data_tmp_1;
			A_B0R8[bn0_ma_reg8] = data_tmp_2;
			A_B1R8[bn1_ma_reg8] = data_tmp_3;
			data_tmp_1 = A_B0R10[bn0_ma_reg7];
			data_tmp_2 = A_B0R11[bn0_ma_reg7];
			A_B0R10[bn0_ma_reg7] = A_B0R9[bn0_ma_reg8];
			A_B0R11[bn0_ma_reg7] = A_B1R9[bn1_ma_reg8];
			A_B0R9[bn0_ma_reg8] = data_tmp_1;
			A_B1R9[bn1_ma_reg8] = data_tmp_2;
			data_tmp_1 = A_B0R11[bn0_ma_reg8];
			A_B0R11[bn0_ma_reg8] = A_B1R10[bn1_ma_reg8];
			A_B1R10[bn1_ma_reg8] = data_tmp_1;
            //-----------------------------------------------------
		    data_tmp_1 = A_B1R13[bn1_ma_reg7];
		    data_tmp_2 = A_B1R14[bn1_ma_reg7];
		    data_tmp_3 = A_B1R15[bn1_ma_reg7];
			A_B1R13[bn1_ma_reg7] = A_B0R12[bn0_ma_reg7];
			A_B1R14[bn1_ma_reg7] = A_B0R12[bn0_ma_reg8];
			A_B1R15[bn1_ma_reg7] = A_B1R12[bn1_ma_reg8];
			A_B0R12[bn0_ma_reg7] = data_tmp_1;
			A_B0R12[bn0_ma_reg8] = data_tmp_2;
			A_B1R12[bn1_ma_reg8] = data_tmp_3;
			data_tmp_1 = A_B0R14[bn0_ma_reg7];
			data_tmp_2 = A_B0R15[bn0_ma_reg7];
			A_B0R14[bn0_ma_reg7] = A_B0R13[bn0_ma_reg8];
			A_B0R15[bn0_ma_reg7] = A_B1R13[bn1_ma_reg8];
			A_B0R13[bn0_ma_reg8] = data_tmp_1;
			A_B1R13[bn1_ma_reg8] = data_tmp_2;
			data_tmp_1 = A_B0R15[bn0_ma_reg8];
			A_B0R15[bn0_ma_reg8] = A_B1R14[bn1_ma_reg8];
			A_B1R14[bn1_ma_reg8] = data_tmp_1;			
		  }	
		 } 
		}
	}
	
	std::cout << "radix-16 FFT computing over!!\n";
	
	//if(display == 1)DATARECORD <<" ----------------------------------------------- \n";
	//if(display == 1)DATARECORD <<" Radix-4 FFT computing start!!! \n";
	//if(display == 1)siang_record <<" ----------------------------------------------- \n";
	//if(display == 1)siang_record <<" Radix-4 FFT computing start!!! \n";
	
	for(int i = 0;i < group; i++){
        for(int j=0; j < radix; j++){
        	gray_i  = Gray(i,group);
        	BC_tmp  = j * group + gray_i;
        	RR_R16_R4(BC_tmp,((4 * (Stage-1)) - 2),BC);
        	AGU_R16(BC,bn_tmp,ma_tmp);
			//if(display == 1)DATARECORD <<"--------------------------------------------- \n";
			//if(display == 1)DATARECORD <<"BC_tmp: "<< BC_tmp <<"\n";
			//if(display == 1)DATARECORD <<"BC   : " << BC     <<"\n";
			//if(display == 1)DATARECORD <<"bn_tmp: "<< bn_tmp <<"\n";
			//if(display == 1)DATARECORD <<"ma_tmp: "<< ma_tmp <<"\n";	
			if(display == 1)siang_record <<"--------------------------------------------- \n";	// siang_record
			if(display == 1)siang_record <<"BC_tmp: "<< BC_tmp <<"\n"; // siang_record
			if(display == 1)siang_record <<"BC   : " << BC     <<"\n"; // siang_record				
        	if(bn_tmp == 0){
				//if(display == 1)DATARECORD <<" Before Radix-4 Butterfly unit computing!!! \n";
				/*if(display == 1)DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";	*/				
                Radix4_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp],A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
        		Radix4_BU(A_B0R4[ma_tmp],A_B0R5[ma_tmp],A_B0R6[ma_tmp],A_B0R7[ma_tmp]);
        		Radix4_BU(A_B0R8[ma_tmp],A_B0R9[ma_tmp],A_B0R10[ma_tmp],A_B0R11[ma_tmp]);
        		Radix4_BU(A_B0R12[ma_tmp],A_B0R13[ma_tmp],A_B0R14[ma_tmp],A_B0R15[ma_tmp]);
				//if(display == 1)DATARECORD <<" After Radix-4 Butterfly unit computing!!! \n";
				/*if(display == 1)DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";		*/			
        	}else {
				//if(display == 1)DATARECORD <<" Before Radix-4 Butterfly unit computing!!! \n";
				/*if(display == 1)DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";	*/				
        	    Radix4_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp],A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
        	    Radix4_BU(A_B1R4[ma_tmp],A_B1R5[ma_tmp],A_B1R6[ma_tmp],A_B1R7[ma_tmp]);
        	    Radix4_BU(A_B1R8[ma_tmp],A_B1R9[ma_tmp],A_B1R10[ma_tmp],A_B1R11[ma_tmp]);
        	    Radix4_BU(A_B1R12[ma_tmp],A_B1R13[ma_tmp],A_B1R14[ma_tmp],A_B1R15[ma_tmp]);
				//if(display == 1)DATARECORD <<" After Radix-4 Butterfly unit computing!!! \n";
				/*if(display == 1)DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";	*/				
        	}			
        }		
	}
	
	int index0_i;
	int index1_i;
	int index2_i;
	int index3_i;
	int index4_i;
	int index5_i;
	int index6_i;
	int index7_i;
	int index8_i;
	int index9_i;
	int index10_i;
	int index11_i;
	int index12_i;
	int index13_i;
	int index14_i;
	int index15_i;
	
	int index0;
    int index1;
    int index2;
    int index3;
    int index4;
    int index5;
    int index6;
    int index7;
    int index8;
    int index9;
    int index10;
    int index11;
    int index12;
    int index13;
    int index14;
    int index15;
	
	std::ofstream NTT_DATA_OUT("./NTT_R16_R4_DATA_OUT.txt");
	
	NTT_DATA_OUT <<"***************************************************************\n";
	NTT_DATA_OUT <<"***** DATA OUTPUT!!                                         ***\n";
	NTT_DATA_OUT <<"***************************************************************\n";
	//data output
	// SPMB data output , output function is "NTT_REORDERINDEX_R16_R4_OUT".

	for(int i = 0;i < group;i++){
		for(int j = 0;j < radix;j++){
			//gray_i  = Gray(i,group);
			BC_tmp  = j * group + i;
			NTT_DATA_OUT <<"---------------------------------------------------\n";
			NTT_DATA_OUT <<" BC_tmp: "<< BC_tmp <<"\n";
			RR_R16_R4(BC_tmp,((4 * (Stage-1)) - 2),BC);
			//DATARECORD <<" BC: "<< BC <<"\n";
			AGU_R16(BC,bn_tmp,ma_tmp);
			NTT_DATA_OUT <<" bn_tmp:  "<< bn_tmp <<"\n";
			NTT_DATA_OUT <<" ma_tmp:  "<< ma_tmp <<"\n";
			index0_i  = 16 * BC;
			index1_i  = 16 * BC  + 1;
			index2_i  = 16 * BC  + 2;
			index3_i  = 16 * BC  + 3;
			index4_i  = 16 * BC  + 4;
			index5_i  = 16 * BC  + 5;
			index6_i  = 16 * BC  + 6;
			index7_i  = 16 * BC  + 7;
			index8_i  = 16 * BC  + 8;
			index9_i  = 16 * BC  + 9;
			index10_i = 16 * BC + 10;
			index11_i = 16 * BC + 11;
			index12_i = 16 * BC + 12;
			index13_i = 16 * BC + 13;
			index14_i = 16 * BC + 14;
			index15_i = 16 * BC + 15;

			NTT_DATA_OUT <<" index0_i:  "<< index0_i <<"\n";
			NTT_DATA_OUT <<" index1_i:  "<< index1_i <<"\n";
			NTT_DATA_OUT <<" index2_i:  "<< index2_i <<"\n";
			NTT_DATA_OUT <<" index3_i:  "<< index3_i <<"\n";
			NTT_DATA_OUT <<" index4_i:  "<< index4_i <<"\n";
			NTT_DATA_OUT <<" index5_i:  "<< index5_i <<"\n";
			NTT_DATA_OUT <<" index6_i:  "<< index6_i <<"\n";
			NTT_DATA_OUT <<" index7_i:  "<< index7_i <<"\n";
			NTT_DATA_OUT <<" index8_i:  "<< index8_i <<"\n";
			NTT_DATA_OUT <<" index9_i:  "<< index9_i <<"\n";
			NTT_DATA_OUT <<" index10_i: "<< index10_i <<"\n";
			NTT_DATA_OUT <<" index11_i: "<< index11_i <<"\n";
			NTT_DATA_OUT <<" index12_i: "<< index12_i <<"\n";
			NTT_DATA_OUT <<" index13_i: "<< index13_i <<"\n";
			NTT_DATA_OUT <<" index14_i: "<< index14_i <<"\n";
			NTT_DATA_OUT <<" index15_i: "<< index15_i <<"\n";
					
			if(bn_tmp == 0){
			   NTT_REORDERINDEX_R16_R4_OUT(index0_i,index0);
			   NTT_REORDERINDEX_R16_R4_OUT(index1_i,index1);
			   NTT_REORDERINDEX_R16_R4_OUT(index2_i,index2);
			   NTT_REORDERINDEX_R16_R4_OUT(index3_i,index3);
			   NTT_REORDERINDEX_R16_R4_OUT(index4_i,index4);
			   NTT_REORDERINDEX_R16_R4_OUT(index5_i,index5);
			   NTT_REORDERINDEX_R16_R4_OUT(index6_i,index6);
			   NTT_REORDERINDEX_R16_R4_OUT(index7_i,index7);
			   NTT_REORDERINDEX_R16_R4_OUT(index8_i,index8);
			   NTT_REORDERINDEX_R16_R4_OUT(index9_i,index9);
			   NTT_REORDERINDEX_R16_R4_OUT(index10_i,index10);
			   NTT_REORDERINDEX_R16_R4_OUT(index11_i,index11);
			   NTT_REORDERINDEX_R16_R4_OUT(index12_i,index12);
			   NTT_REORDERINDEX_R16_R4_OUT(index13_i,index13);
			   NTT_REORDERINDEX_R16_R4_OUT(index14_i,index14);
			   NTT_REORDERINDEX_R16_R4_OUT(index15_i,index15);
			   NTT_DATA_OUT <<" index0:  "<< index0 <<"\n";
			   NTT_DATA_OUT <<" index1:  "<< index1 <<"\n";
			   NTT_DATA_OUT <<" index2:  "<< index2 <<"\n";
			   NTT_DATA_OUT <<" index3:  "<< index3 <<"\n";
			   NTT_DATA_OUT <<" index4:  "<< index4 <<"\n";
			   NTT_DATA_OUT <<" index5:  "<< index5 <<"\n";
			   NTT_DATA_OUT <<" index6:  "<< index6 <<"\n";
			   NTT_DATA_OUT <<" index7:  "<< index7 <<"\n";
			   NTT_DATA_OUT <<" index8:  "<< index8 <<"\n";
			   NTT_DATA_OUT <<" index9:  "<< index9 <<"\n";
			   NTT_DATA_OUT <<" index10: "<< index10 <<"\n";
			   NTT_DATA_OUT <<" index11: "<< index11 <<"\n";
			   NTT_DATA_OUT <<" index12: "<< index12 <<"\n";
			   NTT_DATA_OUT <<" index13: "<< index13 <<"\n";
			   NTT_DATA_OUT <<" index14: "<< index14 <<"\n";
			   NTT_DATA_OUT <<" index15: "<< index15 <<"\n";
               A[index0] = A_B0R0[ma_tmp];
			   A[index1] = A_B0R1[ma_tmp];
			   A[index2] = A_B0R2[ma_tmp];
			   A[index3] = A_B0R3[ma_tmp];
			   A[index4] = A_B0R4[ma_tmp];
			   A[index5] = A_B0R5[ma_tmp];
			   A[index6] = A_B0R6[ma_tmp];
			   A[index7] = A_B0R7[ma_tmp];
			   A[index8] = A_B0R8[ma_tmp];
			   A[index9] = A_B0R9[ma_tmp];
			   A[index10] = A_B0R10[ma_tmp];
			   A[index11] = A_B0R11[ma_tmp];
			   A[index12] = A_B0R12[ma_tmp];
			   A[index13] = A_B0R13[ma_tmp];
			   A[index14] = A_B0R14[ma_tmp];
			   A[index15] = A_B0R15[ma_tmp];
			   NTT_DATA_OUT << "A_B0R0[" << ma_tmp <<"]: "<< A_B0R0[ma_tmp] << "\n";
			   NTT_DATA_OUT << "A_B0R1[" << ma_tmp <<"]: "<< A_B0R1[ma_tmp] << "\n";
			   NTT_DATA_OUT << "A_B0R2[" << ma_tmp <<"]: "<< A_B0R2[ma_tmp] << "\n";
			   NTT_DATA_OUT << "A_B0R3[" << ma_tmp <<"]: "<< A_B0R3[ma_tmp] << "\n";
			   NTT_DATA_OUT << "A_B0R4[" << ma_tmp <<"]: "<< A_B0R4[ma_tmp] << "\n";
			   NTT_DATA_OUT << "A_B0R5[" << ma_tmp <<"]: "<< A_B0R5[ma_tmp] << "\n";
			   NTT_DATA_OUT << "A_B0R6[" << ma_tmp <<"]: "<< A_B0R6[ma_tmp] << "\n";
			   NTT_DATA_OUT << "A_B0R7[" << ma_tmp <<"]: "<< A_B0R7[ma_tmp] << "\n";
			   NTT_DATA_OUT << "A_B0R8[" << ma_tmp <<"]: "<< A_B0R8[ma_tmp] << "\n";
			   NTT_DATA_OUT << "A_B0R9[" << ma_tmp <<"]: "<< A_B0R9[ma_tmp] << "\n";
			   NTT_DATA_OUT << "A_B0R10[" << ma_tmp <<"]: "<< A_B0R10[ma_tmp] << "\n";					
			   NTT_DATA_OUT << "A_B0R11[" << ma_tmp <<"]: "<< A_B0R11[ma_tmp] << "\n";					
			   NTT_DATA_OUT << "A_B0R12[" << ma_tmp <<"]: "<< A_B0R12[ma_tmp] << "\n";					
			   NTT_DATA_OUT << "A_B0R13[" << ma_tmp <<"]: "<< A_B0R13[ma_tmp] << "\n";					
			   NTT_DATA_OUT << "A_B0R14[" << ma_tmp <<"]: "<< A_B0R14[ma_tmp] << "\n";					
			   NTT_DATA_OUT << "A_B0R15[" << ma_tmp <<"]: "<< A_B0R15[ma_tmp] <<"\n";
			}
			else {
			   NTT_REORDERINDEX_R16_R4_OUT(index0_i,index0);
			   NTT_REORDERINDEX_R16_R4_OUT(index1_i,index1);
			   NTT_REORDERINDEX_R16_R4_OUT(index2_i,index2);
			   NTT_REORDERINDEX_R16_R4_OUT(index3_i,index3);
			   NTT_REORDERINDEX_R16_R4_OUT(index4_i,index4);
			   NTT_REORDERINDEX_R16_R4_OUT(index5_i,index5);
			   NTT_REORDERINDEX_R16_R4_OUT(index6_i,index6);
			   NTT_REORDERINDEX_R16_R4_OUT(index7_i,index7);
			   NTT_REORDERINDEX_R16_R4_OUT(index8_i,index8);
			   NTT_REORDERINDEX_R16_R4_OUT(index9_i,index9);
			   NTT_REORDERINDEX_R16_R4_OUT(index10_i,index10);
			   NTT_REORDERINDEX_R16_R4_OUT(index11_i,index11);
			   NTT_REORDERINDEX_R16_R4_OUT(index12_i,index12);
			   NTT_REORDERINDEX_R16_R4_OUT(index13_i,index13);
			   NTT_REORDERINDEX_R16_R4_OUT(index14_i,index14);
			   NTT_REORDERINDEX_R16_R4_OUT(index15_i,index15);
			   NTT_DATA_OUT <<" index0:  "<< index0 <<"\n";
			   NTT_DATA_OUT <<" index1:  "<< index1 <<"\n";
			   NTT_DATA_OUT <<" index2:  "<< index2 <<"\n";
			   NTT_DATA_OUT <<" index3:  "<< index3 <<"\n";
			   NTT_DATA_OUT <<" index4:  "<< index4 <<"\n";
			   NTT_DATA_OUT <<" index5:  "<< index5 <<"\n";
			   NTT_DATA_OUT <<" index6:  "<< index6 <<"\n";
			   NTT_DATA_OUT <<" index7:  "<< index7 <<"\n";
			   NTT_DATA_OUT <<" index8:  "<< index8 <<"\n";
			   NTT_DATA_OUT <<" index9:  "<< index9 <<"\n";
			   NTT_DATA_OUT <<" index10: "<< index10 <<"\n";
			   NTT_DATA_OUT <<" index11: "<< index11 <<"\n";
			   NTT_DATA_OUT <<" index12: "<< index12 <<"\n";
			   NTT_DATA_OUT <<" index13: "<< index13 <<"\n";
			   NTT_DATA_OUT <<" index14: "<< index14 <<"\n";
			   NTT_DATA_OUT <<" index15: "<< index15 <<"\n";
               A[index0]     = A_B1R0[ma_tmp];
               A[index1]     = A_B1R1[ma_tmp];
               A[index2]     = A_B1R2[ma_tmp];
               A[index3]     = A_B1R3[ma_tmp];
               A[index4]     = A_B1R4[ma_tmp];
               A[index5]     = A_B1R5[ma_tmp];
               A[index6]     = A_B1R6[ma_tmp];
               A[index7]     = A_B1R7[ma_tmp];
               A[index8]     = A_B1R8[ma_tmp];
               A[index9]     = A_B1R9[ma_tmp];
               A[index10]    = A_B1R10[ma_tmp];
               A[index11]    = A_B1R11[ma_tmp];
               A[index12]    = A_B1R12[ma_tmp];
               A[index13]    = A_B1R13[ma_tmp];
               A[index14]    = A_B1R14[ma_tmp];
               A[index15]    = A_B1R15[ma_tmp];
			   NTT_DATA_OUT << "A_B1R0[" << ma_tmp <<"]: "<< A_B1R0[ma_tmp] << "\n";
			   NTT_DATA_OUT << "A_B1R1[" << ma_tmp <<"]: "<< A_B1R1[ma_tmp] << "\n";
			   NTT_DATA_OUT << "A_B1R2[" << ma_tmp <<"]: "<< A_B1R2[ma_tmp] << "\n";
			   NTT_DATA_OUT << "A_B1R3[" << ma_tmp <<"]: "<< A_B1R3[ma_tmp] << "\n";
			   NTT_DATA_OUT << "A_B1R4[" << ma_tmp <<"]: "<< A_B1R4[ma_tmp] << "\n";
			   NTT_DATA_OUT << "A_B1R5[" << ma_tmp <<"]: "<< A_B1R5[ma_tmp] << "\n";
			   NTT_DATA_OUT << "A_B1R6[" << ma_tmp <<"]: "<< A_B1R6[ma_tmp] << "\n";
			   NTT_DATA_OUT << "A_B1R7[" << ma_tmp <<"]: "<< A_B1R7[ma_tmp] << "\n";
			   NTT_DATA_OUT << "A_B1R8[" << ma_tmp <<"]: "<< A_B1R8[ma_tmp] << "\n";
			   NTT_DATA_OUT << "A_B1R9[" << ma_tmp <<"]: "<< A_B1R9[ma_tmp] << "\n";
			   NTT_DATA_OUT << "A_B1R10[" << ma_tmp <<"]: "<< A_B1R10[ma_tmp] << "\n";					
			   NTT_DATA_OUT << "A_B1R11[" << ma_tmp <<"]: "<< A_B1R11[ma_tmp] << "\n";					
			   NTT_DATA_OUT << "A_B1R12[" << ma_tmp <<"]: "<< A_B1R12[ma_tmp] << "\n";					
			   NTT_DATA_OUT << "A_B1R13[" << ma_tmp <<"]: "<< A_B1R13[ma_tmp] << "\n";					
			   NTT_DATA_OUT << "A_B1R14[" << ma_tmp <<"]: "<< A_B1R14[ma_tmp] << "\n";					
			   NTT_DATA_OUT << "A_B1R15[" << ma_tmp <<"]: "<< A_B1R15[ma_tmp] <<"\n";			   
			}			
		}		
	}
	
	//************************************************
	//point : 1024P
	//          [0] [1]  [2]  [3]  [4]  [5]  [6]  [7]  [8]  [9]  [10] [11] [12] [13] [14] [15]
	//Original: 0 , 256, 512, 768, 64,  320, 576, 832, 128, 384, 640, 896, 192, 448, 704, 960
	//          [0] [1]  [2]  [3]  [4]  [5]  [6]  [7]  [8]  [9]  [10] [11] [12] [13] [14] [15]
	//--------> 0 , 64,  128, 192, 256, 320, 384, 448, 512, 576, 640, 704, 768, 832, 896, 960
	//data relocation for INTT input order
	std::ofstream BEFORE_DATA_RELOCATION("./NTT_R16_R4_Before_DATA_RELOCATION.txt");
	std::ofstream AFTER_DATA_RELOCATION("./NTT_R16_R4_After_DATA_RELOCATION.txt");
	
	for(int i = 0;i < group;i++){
		for(int j = 0;j < radix;j++){
	        gray_i  = Gray(i,group);
			BC_tmp  = j * group + gray_i;
	        RR_R16_R4(BC_tmp,((4 * (Stage-1)) - 2),BC);
			AGU_R16(BC,bn_tmp,ma_tmp);
			BEFORE_DATA_RELOCATION <<"--------------------------------------------- \n";
			BEFORE_DATA_RELOCATION <<"i: "     << i <<"\n";
			BEFORE_DATA_RELOCATION <<"j: "     << j <<"\n";
			BEFORE_DATA_RELOCATION <<"gray_i: "<< gray_i <<"\n";
			BEFORE_DATA_RELOCATION <<"BC_tmp: "<< BC_tmp <<"\n";
			BEFORE_DATA_RELOCATION <<"BC:    " << BC     <<"\n";
			BEFORE_DATA_RELOCATION <<"bn_tmp: "<< bn_tmp <<"\n";
			BEFORE_DATA_RELOCATION <<"ma_tmp: "<< ma_tmp <<"\n";
			if(bn_tmp == 0){
			   BEFORE_DATA_RELOCATION <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]<<"\n";
			   BEFORE_DATA_RELOCATION <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]<<"\n";
			   BEFORE_DATA_RELOCATION <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]<<"\n";
			   BEFORE_DATA_RELOCATION <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]<<"\n";
			   BEFORE_DATA_RELOCATION <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]<<"\n";
			   BEFORE_DATA_RELOCATION <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]<<"\n";
			   BEFORE_DATA_RELOCATION <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]<<"\n";
			   BEFORE_DATA_RELOCATION <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]<<"\n";
			   BEFORE_DATA_RELOCATION <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]<<"\n";
			   BEFORE_DATA_RELOCATION <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]<<"\n";
			   BEFORE_DATA_RELOCATION <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";
			   BEFORE_DATA_RELOCATION <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";
			   BEFORE_DATA_RELOCATION <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";
			   BEFORE_DATA_RELOCATION <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";
			   BEFORE_DATA_RELOCATION <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";
			   BEFORE_DATA_RELOCATION <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";					
               data_tmp_1  = A_B0R1[ma_tmp];
               data_tmp_2  = A_B0R2[ma_tmp];
               data_tmp_3  = A_B0R3[ma_tmp];
               data_tmp_4  = A_B0R4[ma_tmp];
               data_tmp_5  = A_B0R5[ma_tmp];
               data_tmp_6  = A_B0R6[ma_tmp];
               data_tmp_7  = A_B0R7[ma_tmp];
               data_tmp_8  = A_B0R8[ma_tmp];
               data_tmp_9  = A_B0R9[ma_tmp];
               data_tmp_10 = A_B0R10[ma_tmp];
               data_tmp_11 = A_B0R11[ma_tmp];
               data_tmp_12 = A_B0R12[ma_tmp];
               data_tmp_13 = A_B0R13[ma_tmp];
               data_tmp_14 = A_B0R14[ma_tmp];
               data_tmp_15 = A_B0R15[ma_tmp];
			   A_B0R1[ma_tmp]  = data_tmp_4;
			   A_B0R2[ma_tmp]  = data_tmp_8;
			   A_B0R3[ma_tmp]  = data_tmp_12;
			   A_B0R4[ma_tmp]  = data_tmp_1;
			   A_B0R6[ma_tmp]  = data_tmp_9;
			   A_B0R7[ma_tmp]  = data_tmp_13;
			   A_B0R8[ma_tmp]  = data_tmp_2;
			   A_B0R9[ma_tmp]  = data_tmp_6;
			   A_B0R11[ma_tmp] = data_tmp_14;
			   A_B0R12[ma_tmp] = data_tmp_3;
			   A_B0R13[ma_tmp] = data_tmp_7;
			   A_B0R14[ma_tmp] = data_tmp_11;
			}else {
			   BEFORE_DATA_RELOCATION <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]<<"\n";
			   BEFORE_DATA_RELOCATION <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]<<"\n";
			   BEFORE_DATA_RELOCATION <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]<<"\n";
			   BEFORE_DATA_RELOCATION <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]<<"\n";
			   BEFORE_DATA_RELOCATION <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]<<"\n";
			   BEFORE_DATA_RELOCATION <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]<<"\n";
			   BEFORE_DATA_RELOCATION <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]<<"\n";
			   BEFORE_DATA_RELOCATION <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]<<"\n";
			   BEFORE_DATA_RELOCATION <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]<<"\n";
			   BEFORE_DATA_RELOCATION <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]<<"\n";
			   BEFORE_DATA_RELOCATION <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";
			   BEFORE_DATA_RELOCATION <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";
			   BEFORE_DATA_RELOCATION <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";
			   BEFORE_DATA_RELOCATION <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";
			   BEFORE_DATA_RELOCATION <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";
			   BEFORE_DATA_RELOCATION <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";					
               data_tmp_1  = A_B1R1[ma_tmp];
               data_tmp_2  = A_B1R2[ma_tmp];
               data_tmp_3  = A_B1R3[ma_tmp];
               data_tmp_4  = A_B1R4[ma_tmp];
               data_tmp_5  = A_B1R5[ma_tmp];
               data_tmp_6  = A_B1R6[ma_tmp];
               data_tmp_7  = A_B1R7[ma_tmp];
               data_tmp_8  = A_B1R8[ma_tmp];
               data_tmp_9  = A_B1R9[ma_tmp];
               data_tmp_10 = A_B1R10[ma_tmp];
               data_tmp_11 = A_B1R11[ma_tmp];
               data_tmp_12 = A_B1R12[ma_tmp];
               data_tmp_13 = A_B1R13[ma_tmp];
               data_tmp_14 = A_B1R14[ma_tmp];
               data_tmp_15 = A_B1R15[ma_tmp];
			   A_B1R1[ma_tmp]  = data_tmp_4;
			   A_B1R2[ma_tmp]  = data_tmp_8;
			   A_B1R3[ma_tmp]  = data_tmp_12;
			   A_B1R4[ma_tmp]  = data_tmp_1;
			   A_B1R6[ma_tmp]  = data_tmp_9;
			   A_B1R7[ma_tmp]  = data_tmp_13;
			   A_B1R8[ma_tmp]  = data_tmp_2;
			   A_B1R9[ma_tmp]  = data_tmp_6;
			   A_B1R11[ma_tmp] = data_tmp_14;
			   A_B1R12[ma_tmp] = data_tmp_3;
			   A_B1R13[ma_tmp] = data_tmp_7;
			   A_B1R14[ma_tmp] = data_tmp_11;
			}
	    }
	}
	
	//After data relocation for INTT output
	for(int i = 0;i < group;i++){
		for(int j = 0; j < radix; j++){
        	gray_i  = Gray(i,group);
        	BC_tmp  = j * group + gray_i;
        	BC = BC_tmp;
        	AGU_R16(BC,bn_tmp,ma_tmp);
			AFTER_DATA_RELOCATION <<"--------------------------------------------- \n";
			AFTER_DATA_RELOCATION <<"i: "     << i <<"\n";
			AFTER_DATA_RELOCATION <<"j: "     << j <<"\n";
			AFTER_DATA_RELOCATION <<"gray_i: "<< gray_i <<"\n";			
			AFTER_DATA_RELOCATION <<"BC_tmp: "<< BC_tmp <<"\n";
			AFTER_DATA_RELOCATION <<"BC: "    << BC     <<"\n";
			AFTER_DATA_RELOCATION <<"bn_tmp: "<< bn_tmp <<"\n";
			AFTER_DATA_RELOCATION <<"ma_tmp: "<< ma_tmp <<"\n";			
        	if(bn_tmp == 0){
				AFTER_DATA_RELOCATION <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";						
        	}else {
				AFTER_DATA_RELOCATION <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";								
        	}	
		}
	}
	
	std::ofstream  B0R0_SPMB_O("./B0R0_R16_R4_SPMB.txt");
	
	for(int ss = 0; ss < word_size; ss++){
		//bn0
		B0R0[ss]  = A_B0R0[ss];
		B0R1[ss]  = A_B0R1[ss];
		B0R2[ss]  = A_B0R2[ss];
		B0R3[ss]  = A_B0R3[ss];
		B0R4[ss]  = A_B0R4[ss];
		B0R5[ss]  = A_B0R5[ss];
		B0R6[ss]  = A_B0R6[ss];
		B0R7[ss]  = A_B0R7[ss];
		B0R8[ss]  = A_B0R8[ss];
		B0R9[ss]  = A_B0R9[ss];
		B0R10[ss] = A_B0R10[ss];
		B0R11[ss] = A_B0R11[ss];
		B0R12[ss] = A_B0R12[ss];
		B0R13[ss] = A_B0R13[ss];
		B0R14[ss] = A_B0R14[ss];
		B0R15[ss] = A_B0R15[ss];
		if(ss==0)B0R0_SPMB_O << B0R0[ss] <<"\n";
		if(ss==0)B0R0_SPMB_O << B0R1[ss] <<"\n";
		if(ss==0)B0R0_SPMB_O << B0R2[ss] <<"\n";
		if(ss==0)B0R0_SPMB_O << B0R3[ss] <<"\n";
		if(ss==0)B0R0_SPMB_O << B0R4[ss] <<"\n";
		if(ss==0)B0R0_SPMB_O << B0R5[ss] <<"\n";
		if(ss==0)B0R0_SPMB_O << B0R6[ss] <<"\n";
		if(ss==0)B0R0_SPMB_O << B0R7[ss] <<"\n";
		if(ss==0)B0R0_SPMB_O << B0R8[ss] <<"\n";
		if(ss==0)B0R0_SPMB_O << B0R9[ss] <<"\n";
		if(ss==0)B0R0_SPMB_O << B0R10[ss] <<"\n";
		if(ss==0)B0R0_SPMB_O << B0R11[ss] <<"\n";
		if(ss==0)B0R0_SPMB_O << B0R12[ss] <<"\n";
		if(ss==0)B0R0_SPMB_O << B0R13[ss] <<"\n";
		if(ss==0)B0R0_SPMB_O << B0R14[ss] <<"\n";
		if(ss==0)B0R0_SPMB_O << B0R15[ss] <<"\n";
		//----------------------------
		//bn1
		B1R0[ss]  = A_B1R0[ss];
		B1R1[ss]  = A_B1R1[ss];
		B1R2[ss]  = A_B1R2[ss];
		B1R3[ss]  = A_B1R3[ss];
		B1R4[ss]  = A_B1R4[ss];
		B1R5[ss]  = A_B1R5[ss];
		B1R6[ss]  = A_B1R6[ss];
		B1R7[ss]  = A_B1R7[ss];
		B1R8[ss]  = A_B1R8[ss];
		B1R9[ss]  = A_B1R9[ss];
		B1R10[ss] = A_B1R10[ss];
		B1R11[ss] = A_B1R11[ss];
		B1R12[ss] = A_B1R12[ss];
		B1R13[ss] = A_B1R13[ss];
		B1R14[ss] = A_B1R14[ss];
		B1R15[ss] = A_B1R15[ss];		
	}	
}
//INTT Mixed radix-16 and radix-4
void NTTSPMB::INTT_r16_r4(std::vector<ZZ> &A,std::vector<ZZ> &B0R0,std::vector<ZZ> &B0R1,std::vector<ZZ> &B0R2,std::vector<ZZ> &B0R3,
std::vector<ZZ> &B0R4,std::vector<ZZ> &B0R5,std::vector<ZZ> &B0R6,std::vector<ZZ> &B0R7,std::vector<ZZ> &B0R8,std::vector<ZZ> &B0R9,
std::vector<ZZ> &B0R10,std::vector<ZZ> &B0R11,std::vector<ZZ> &B0R12,std::vector<ZZ> &B0R13,std::vector<ZZ> &B0R14,std::vector<ZZ> &B0R15,
std::vector<ZZ> &B1R0,std::vector<ZZ> &B1R1,std::vector<ZZ> &B1R2,std::vector<ZZ> &B1R3,std::vector<ZZ> &B1R4,std::vector<ZZ> &B1R5,
std::vector<ZZ> &B1R6,std::vector<ZZ> &B1R7,std::vector<ZZ> &B1R8,std::vector<ZZ> &B1R9,std::vector<ZZ> &B1R10,std::vector<ZZ> &B1R11,
std::vector<ZZ> &B1R12,std::vector<ZZ> &B1R13,std::vector<ZZ> &B1R14,std::vector<ZZ> &B1R15){
	
	int            Stage;//stage	
	int            word_size;
	int            offset;
	int            BC;
	int            BC_RR;
	int            BC_tmp;
    int            BC_Reorder;
	int            group;
    int            ma_tmp;
	int            bn_tmp;
    int            bit_tmp;
	int            bn0_bc_tmp;//frist in bc data
	int            bn1_bc_tmp;//frist in bc data
	int            bn0_ma_reg1;
	int            bn0_ma_reg2;
	int            bn0_ma_reg3;
	int            bn0_ma_reg4;
	int            bn0_ma_reg5;
	int            bn0_ma_reg6;
	int            bn0_ma_reg7;
	int            bn0_ma_reg8;
	int            bn1_ma_reg1;
	int            bn1_ma_reg2;
	int            bn1_ma_reg3;
	int            bn1_ma_reg4;
	int            bn1_ma_reg5;
	int            bn1_ma_reg6;
	int            bn1_ma_reg7;
	int            bn1_ma_reg8;
	int            gray_i;
    int            BC_WIDTH;
	int            tw_modulus;
    int            tw_modulus_tmp;
	//display
	int            display;
    double         Stage_double;	
    std::vector<int> bit_array_tmp;
	
	std::ofstream DATARECORD("./INTT_R16_R4_SPMB.txt");
	//**************************************************************
	//**************************************************************
	display = 1;
	//**************************************************************
	//radix-16 Stage
    Stage_double  = log2(N);
    Stage       =  (int)floor(Stage_double/4);
	BC_WIDTH    =  (int)ceil(log2(N/16));
	offset      =  (int)N /  16;
	word_size   =  (int)N / (2 * 16);
	group       =  (int)N / (256);
	tw_modulus  =  (int)N /  16;
	bit_array_tmp.resize(BC_WIDTH);

    DATARECORD << "group: "    << group << "\n";
    DATARECORD << "BC_WIDTH: " << BC_WIDTH << "\n";
	
	ZZ               data_tmp_1;
	ZZ               data_tmp_2;
	ZZ               data_tmp_3;
	ZZ               data_tmp_4;
	ZZ               data_tmp_5;
	ZZ               data_tmp_6;
	ZZ               data_tmp_7;
	ZZ               data_tmp_8;
	ZZ               data_tmp_9;
	ZZ               data_tmp_10;
	ZZ               data_tmp_11;
	ZZ               data_tmp_12;
	ZZ               data_tmp_13;
	ZZ               data_tmp_14;
	ZZ               data_tmp_15;
	std::vector<ZZ>  A_B0R0;
	std::vector<ZZ>  A_B0R1;
	std::vector<ZZ>  A_B0R2;
	std::vector<ZZ>  A_B0R3;
	std::vector<ZZ>  A_B0R4;
	std::vector<ZZ>  A_B0R5;
	std::vector<ZZ>  A_B0R6;
	std::vector<ZZ>  A_B0R7;
	std::vector<ZZ>  A_B0R8;
	std::vector<ZZ>  A_B0R9;
	std::vector<ZZ>  A_B0R10;
	std::vector<ZZ>  A_B0R11;
	std::vector<ZZ>  A_B0R12;
	std::vector<ZZ>  A_B0R13;
	std::vector<ZZ>  A_B0R14;
	std::vector<ZZ>  A_B0R15;
	std::vector<ZZ>  A_B1R0;
	std::vector<ZZ>  A_B1R1;
	std::vector<ZZ>  A_B1R2;
	std::vector<ZZ>  A_B1R3;
	std::vector<ZZ>  A_B1R4;
	std::vector<ZZ>  A_B1R5;
	std::vector<ZZ>  A_B1R6;
	std::vector<ZZ>  A_B1R7;
	std::vector<ZZ>  A_B1R8;
	std::vector<ZZ>  A_B1R9;
	std::vector<ZZ>  A_B1R10;
	std::vector<ZZ>  A_B1R11;
	std::vector<ZZ>  A_B1R12;
	std::vector<ZZ>  A_B1R13;
	std::vector<ZZ>  A_B1R14;
	std::vector<ZZ>  A_B1R15;
	
	A_B0R0.resize(word_size);
	A_B0R1.resize(word_size);
	A_B0R2.resize(word_size);
	A_B0R3.resize(word_size);
	A_B0R4.resize(word_size);
	A_B0R5.resize(word_size);
	A_B0R6.resize(word_size);
	A_B0R7.resize(word_size);
	A_B0R8.resize(word_size);
	A_B0R9.resize(word_size);
	A_B0R10.resize(word_size);
	A_B0R11.resize(word_size);
	A_B0R12.resize(word_size);
	A_B0R13.resize(word_size);
	A_B0R14.resize(word_size);
	A_B0R15.resize(word_size);
	A_B1R0.resize(word_size);
	A_B1R1.resize(word_size);
	A_B1R2.resize(word_size);
	A_B1R3.resize(word_size);
	A_B1R4.resize(word_size);
	A_B1R5.resize(word_size);
	A_B1R6.resize(word_size);
	A_B1R7.resize(word_size);
	A_B1R8.resize(word_size);
	A_B1R9.resize(word_size);
	A_B1R10.resize(word_size);
	A_B1R11.resize(word_size);
	A_B1R12.resize(word_size);
	A_B1R13.resize(word_size);
	A_B1R14.resize(word_size);
	A_B1R15.resize(word_size);

	//----------------------------------------------------
	int length;
	ZZ  factor;   //base factor
	ZZ  factor_t; //acctually mul factor
	ZZ  factor_2t;
	ZZ  factor_3t;
	ZZ  factor_4t;
	ZZ  factor_5t;
	ZZ  factor_6t;
	ZZ  factor_7t;
	ZZ  factor_8t;
	ZZ  factor_9t;
	ZZ  factor_10t;
	ZZ  factor_11t;
	ZZ  factor_12t;
	ZZ  factor_13t;
	ZZ  factor_14t;
	ZZ  factor_15t;

    //--------------------------------------
	std::cout << "init load start! \n";
	for(int i=0; i < word_size;i++){	
		A_B0R0[i]   = B0R0[i];
		A_B0R1[i]   = B0R1[i];
		A_B0R2[i]   = B0R2[i];
		A_B0R3[i]   = B0R3[i];
		A_B0R4[i]   = B0R4[i];
		A_B0R5[i]   = B0R5[i];
		A_B0R6[i]   = B0R6[i];
		A_B0R7[i]   = B0R7[i];
		A_B0R8[i]   = B0R8[i];
		A_B0R9[i]   = B0R9[i];
		A_B0R10[i]  = B0R10[i];
		A_B0R11[i]  = B0R11[i];
		A_B0R12[i]  = B0R12[i];
		A_B0R13[i]  = B0R13[i];
		A_B0R14[i]  = B0R14[i];
		A_B0R15[i]  = B0R15[i];
		A_B1R0[i]   = B1R0[i];
		A_B1R1[i]   = B1R1[i];
		A_B1R2[i]   = B1R2[i];
		A_B1R3[i]   = B1R3[i];
		A_B1R4[i]   = B1R4[i];
		A_B1R5[i]   = B1R5[i];
		A_B1R6[i]   = B1R6[i];
		A_B1R7[i]   = B1R7[i];
		A_B1R8[i]   = B1R8[i];
		A_B1R9[i]   = B1R9[i];
		A_B1R10[i]  = B1R10[i];
		A_B1R11[i]  = B1R11[i];
		A_B1R12[i]  = B1R12[i];
		A_B1R13[i]  = B1R13[i];
		A_B1R14[i]  = B1R14[i];
		A_B1R15[i]  = B1R15[i];
	}    

	ma_tmp      = 0;
	bn_tmp      = 0;
	BC          = 0;
	BC_RR       = 0;
	BC_Reorder  = 0;
	std::cout << "init load over! \n";
	//need modify to mult by twiddle factor
	for(int s = 0; s < Stage; s++){
		if(s == 0)factor = IW;
		else {
			SqrMod(factor,factor,p);
			SqrMod(factor,factor,p);
			SqrMod(factor,factor,p);
			SqrMod(factor,factor,p);
		}
		tw_modulus_tmp  = tw_modulus >> ( 4 * s);
		if(display == 1)DATARECORD << "*********************************\n";
		if(display == 1)DATARECORD << "NOW stage :"<< s <<"\n";
		for(int i = 0 ;i < group;i++){
			bn0_bc_tmp  = 0;
			bn1_bc_tmp  = 0;
			for(int j = 0;j < radix;j++){
				gray_i  = Gray(i,group);
			    BC_tmp  = j * group + gray_i;
				if(display == 1)DATARECORD << "-----------------------------------\n";
				if(display == 1)DATARECORD << "Now stage :"<< s <<"\n";
				if(display == 1)DATARECORD << "i: "<< i <<"\n";
				if(display == 1)DATARECORD << "j: "<< j <<"\n";
				if(display == 1)DATARECORD << "gray_i: "<< gray_i <<"\n";
				if(display == 1)DATARECORD << "BC_tmp: "<< BC_tmp <<"\n";
				if(s == Stage - 1) RR_R16_R4(BC_tmp,(4 * s - 2),BC_RR);
				else RR_R16(BC_tmp,s,BC_RR);
				if(display == 1)DATARECORD << "BC_RR: "<< BC_RR <<"\n";
				length = BC_RR % tw_modulus_tmp;
				PowerMod(factor_t,factor,length,p);
				BC_IFFT_Reorder_R16_R4(BC_RR,BC_Reorder);
				if(display == 1)DATARECORD << "BC_Reorder: "<< BC_Reorder <<"\n";
				AGU_R16(BC_Reorder,bn_tmp,ma_tmp);
				if(display == 1)DATARECORD << "bn_tmp: "<< bn_tmp <<"\n";
				if(display == 1)DATARECORD << "ma_tmp: "<< ma_tmp <<"\n";
				if(bn_tmp == 0){
					if(j < 2)bn0_bc_tmp = BC_tmp;
					PowerMod(factor_2t,factor_t,2,p);
					PowerMod(factor_3t,factor_t,3,p);
					PowerMod(factor_4t,factor_t,4,p);
					PowerMod(factor_5t,factor_t,5,p);
					PowerMod(factor_6t,factor_t,6,p);
					PowerMod(factor_7t,factor_t,7,p);
					PowerMod(factor_8t,factor_t,8,p);
					PowerMod(factor_9t,factor_t,9,p);
					PowerMod(factor_10t,factor_t,10,p);
					PowerMod(factor_11t,factor_t,11,p);
					PowerMod(factor_12t,factor_t,12,p);
					PowerMod(factor_13t,factor_t,13,p);
					PowerMod(factor_14t,factor_t,14,p);
					PowerMod(factor_15t,factor_t,15,p);
					if(display == 1)DATARECORD <<"factor_t: "<< factor_t <<"\n";
					if(display == 1)DATARECORD <<"factor_2t: "<< factor_2t <<"\n";
					if(display == 1)DATARECORD <<"factor_3t: "<< factor_3t <<"\n";
					if(display == 1)DATARECORD <<"factor_4t: "<< factor_4t <<"\n";
					if(display == 1)DATARECORD <<"factor_5t: "<< factor_5t <<"\n";
					if(display == 1)DATARECORD <<"factor_6t: "<< factor_6t <<"\n";
					if(display == 1)DATARECORD <<"factor_7t: "<< factor_7t <<"\n";
					if(display == 1)DATARECORD <<"factor_8t: "<< factor_8t <<"\n";
					if(display == 1)DATARECORD <<"factor_9t: "<< factor_9t <<"\n";
					if(display == 1)DATARECORD <<"factor_10t: "<< factor_10t <<"\n";
					if(display == 1)DATARECORD <<"factor_11t: "<< factor_11t <<"\n";
					if(display == 1)DATARECORD <<"factor_12t: "<< factor_12t <<"\n";
					if(display == 1)DATARECORD <<"factor_13t: "<< factor_13t <<"\n";
					if(display == 1)DATARECORD <<"factor_14t: "<< factor_14t <<"\n";
					if(display == 1)DATARECORD <<"factor_15t: "<< factor_15t <<"\n";					
					if(display == 1)DATARECORD <<" Before Radix-16 Butterfly unit computing!!! \n";
					if(display == 1)DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";					
					Radix16_BU_INTT(A_B0R0[ma_tmp],A_B0R1[ma_tmp], A_B0R2[ma_tmp], A_B0R3[ma_tmp],A_B0R4[ma_tmp],
							   A_B0R5[ma_tmp],A_B0R6[ma_tmp], A_B0R7[ma_tmp], A_B0R8[ma_tmp],A_B0R9[ma_tmp],
							   A_B0R10[ma_tmp],A_B0R11[ma_tmp],A_B0R12[ma_tmp],A_B0R13[ma_tmp],A_B0R14[ma_tmp],
							   A_B0R15[ma_tmp]);
					if(display == 1)DATARECORD <<" After Radix-16 Butterfly unit computing!!! \n";
					if(display == 1)DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";					
					MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],factor_t,p);
					MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],factor_2t,p);
					MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],factor_3t,p);
					MulMod(A_B0R4[ma_tmp],A_B0R4[ma_tmp],factor_4t,p);
					MulMod(A_B0R5[ma_tmp],A_B0R5[ma_tmp],factor_5t,p);
					MulMod(A_B0R6[ma_tmp],A_B0R6[ma_tmp],factor_6t,p);
					MulMod(A_B0R7[ma_tmp],A_B0R7[ma_tmp],factor_7t,p);
					MulMod(A_B0R8[ma_tmp],A_B0R8[ma_tmp],factor_8t,p);
					MulMod(A_B0R9[ma_tmp],A_B0R9[ma_tmp],factor_9t,p);
					MulMod(A_B0R10[ma_tmp],A_B0R10[ma_tmp],factor_10t,p);
					MulMod(A_B0R11[ma_tmp],A_B0R11[ma_tmp],factor_11t,p);
					MulMod(A_B0R12[ma_tmp],A_B0R12[ma_tmp],factor_12t,p);
					MulMod(A_B0R13[ma_tmp],A_B0R13[ma_tmp],factor_13t,p);
					MulMod(A_B0R14[ma_tmp],A_B0R14[ma_tmp],factor_14t,p);
					MulMod(A_B0R15[ma_tmp],A_B0R15[ma_tmp],factor_15t,p);
					if(display == 1)DATARECORD <<" After mult by twiddle factor!!! \n";
					if(display == 1)DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";					
					if(j <  2)bn0_ma_reg1 = ma_tmp;
					if((j >= 2)  && (j < 4))bn0_ma_reg2 = ma_tmp;
					if((j >= 4)  && (j < 6))bn0_ma_reg3 = ma_tmp;
					if((j >= 6)  && (j < 8))bn0_ma_reg4 = ma_tmp;
					if((j >= 8)  && (j < 10))bn0_ma_reg5 = ma_tmp;
					if((j >= 10) && (j < 12))bn0_ma_reg6 = ma_tmp;
					if((j >= 12) && (j < 14))bn0_ma_reg7 = ma_tmp;
					if((j >= 14) && (j < 16))bn0_ma_reg8 = ma_tmp;
				}
			    else {
					if(j < 2)bn1_bc_tmp = BC_tmp;				
					PowerMod(factor_2t,factor_t,2,p);
					PowerMod(factor_3t,factor_t,3,p);
					PowerMod(factor_4t,factor_t,4,p);
					PowerMod(factor_5t,factor_t,5,p);
					PowerMod(factor_6t,factor_t,6,p);
					PowerMod(factor_7t,factor_t,7,p);
					PowerMod(factor_8t,factor_t,8,p);
					PowerMod(factor_9t,factor_t,9,p);
					PowerMod(factor_10t,factor_t,10,p);
					PowerMod(factor_11t,factor_t,11,p);
					PowerMod(factor_12t,factor_t,12,p);
					PowerMod(factor_13t,factor_t,13,p);
					PowerMod(factor_14t,factor_t,14,p);
					PowerMod(factor_15t,factor_t,15,p);
					if(display == 1)DATARECORD <<"factor_t: "<< factor_t <<"\n";
					if(display == 1)DATARECORD <<"factor_2t: "<< factor_2t <<"\n";
					if(display == 1)DATARECORD <<"factor_3t: "<< factor_3t <<"\n";
					if(display == 1)DATARECORD <<"factor_4t: "<< factor_4t <<"\n";
					if(display == 1)DATARECORD <<"factor_5t: "<< factor_5t <<"\n";
					if(display == 1)DATARECORD <<"factor_6t: "<< factor_6t <<"\n";
					if(display == 1)DATARECORD <<"factor_7t: "<< factor_7t <<"\n";
					if(display == 1)DATARECORD <<"factor_8t: "<< factor_8t <<"\n";
					if(display == 1)DATARECORD <<"factor_9t: "<< factor_9t <<"\n";
					if(display == 1)DATARECORD <<"factor_10t: "<< factor_10t <<"\n";
					if(display == 1)DATARECORD <<"factor_11t: "<< factor_11t <<"\n";
					if(display == 1)DATARECORD <<"factor_12t: "<< factor_12t <<"\n";
					if(display == 1)DATARECORD <<"factor_13t: "<< factor_13t <<"\n";
					if(display == 1)DATARECORD <<"factor_14t: "<< factor_14t <<"\n";
					if(display == 1)DATARECORD <<"factor_15t: "<< factor_15t <<"\n";					
					if(display == 1)DATARECORD <<" Before Radix-16 Butterfly unit computing!!! \n";
					if(display == 1)DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";						
					Radix16_BU_INTT(A_B1R0[ma_tmp],A_B1R1[ma_tmp], A_B1R2[ma_tmp], A_B1R3[ma_tmp],A_B1R4[ma_tmp],
							   A_B1R5[ma_tmp],A_B1R6[ma_tmp], A_B1R7[ma_tmp], A_B1R8[ma_tmp],A_B1R9[ma_tmp],
							   A_B1R10[ma_tmp],A_B1R11[ma_tmp],A_B1R12[ma_tmp],A_B1R13[ma_tmp],A_B1R14[ma_tmp],
							   A_B1R15[ma_tmp]);
					if(display == 1)DATARECORD <<" After Radix-16 Butterfly unit computing!!! \n";
					if(display == 1)DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";								   
					MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],factor_t,p);
					MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],factor_2t,p);
					MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],factor_3t,p);
					MulMod(A_B1R4[ma_tmp],A_B1R4[ma_tmp],factor_4t,p);
					MulMod(A_B1R5[ma_tmp],A_B1R5[ma_tmp],factor_5t,p);
					MulMod(A_B1R6[ma_tmp],A_B1R6[ma_tmp],factor_6t,p);
					MulMod(A_B1R7[ma_tmp],A_B1R7[ma_tmp],factor_7t,p);
					MulMod(A_B1R8[ma_tmp],A_B1R8[ma_tmp],factor_8t,p);
					MulMod(A_B1R9[ma_tmp],A_B1R9[ma_tmp],factor_9t,p);
					MulMod(A_B1R10[ma_tmp],A_B1R10[ma_tmp],factor_10t,p);
					MulMod(A_B1R11[ma_tmp],A_B1R11[ma_tmp],factor_11t,p);
					MulMod(A_B1R12[ma_tmp],A_B1R12[ma_tmp],factor_12t,p);
					MulMod(A_B1R13[ma_tmp],A_B1R13[ma_tmp],factor_13t,p);
					MulMod(A_B1R14[ma_tmp],A_B1R14[ma_tmp],factor_14t,p);
					MulMod(A_B1R15[ma_tmp],A_B1R15[ma_tmp],factor_15t,p);
					if(display == 1)DATARECORD <<" After mult by twiddle factor!!! \n";
					if(display == 1)DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";
					if(display == 1)DATARECORD <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";						
                    if(j <  2)bn1_ma_reg1 = ma_tmp;					
                    if((j >= 2)  && (j < 4))bn1_ma_reg2 = ma_tmp;
                    if((j >= 4)  && (j < 6))bn1_ma_reg3 = ma_tmp;
                    if((j >= 6)  && (j < 8))bn1_ma_reg4 = ma_tmp;
                    if((j >= 8)  && (j < 10))bn1_ma_reg5 = ma_tmp;
                    if((j >= 10) && (j < 12))bn1_ma_reg6 = ma_tmp;
                    if((j >= 12) && (j < 14))bn1_ma_reg7 = ma_tmp;
                    if((j >= 14) && (j < 16))bn1_ma_reg8 = ma_tmp;					
				}
			}
		//data relocation
		 if(s < Stage-1){
		  if(bn1_bc_tmp > bn0_bc_tmp){
			 data_tmp_1  = A_B0R1[bn0_ma_reg1];
			 data_tmp_2  = A_B0R2[bn0_ma_reg1];
			 data_tmp_3  = A_B0R3[bn0_ma_reg1];
			 data_tmp_4  = A_B0R4[bn0_ma_reg1];
			 data_tmp_5  = A_B0R5[bn0_ma_reg1];
			 data_tmp_6  = A_B0R6[bn0_ma_reg1];
			 data_tmp_7  = A_B0R7[bn0_ma_reg1];
			 data_tmp_8  = A_B0R8[bn0_ma_reg1];
			 data_tmp_9  = A_B0R9[bn0_ma_reg1];
			 data_tmp_10 = A_B0R10[bn0_ma_reg1];
			 data_tmp_11 = A_B0R11[bn0_ma_reg1];
			 data_tmp_12 = A_B0R12[bn0_ma_reg1];
			 data_tmp_13 = A_B0R13[bn0_ma_reg1];
			 data_tmp_14 = A_B0R14[bn0_ma_reg1];
			 data_tmp_15 = A_B0R15[bn0_ma_reg1];
			 A_B0R1[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg1];
			 A_B0R2[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg2];
			 A_B0R3[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg2];
			 A_B0R4[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg3];
			 A_B0R5[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg3];
			 A_B0R6[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg1] = A_B0R0[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg1] = A_B1R0[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg1] = A_B0R0[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg1] = A_B1R0[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg1] = A_B1R0[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg1] = A_B0R0[bn0_ma_reg8];
			 A_B1R0[bn1_ma_reg1]  = data_tmp_1; 
			 A_B1R0[bn1_ma_reg2]  = data_tmp_2; 
			 A_B0R0[bn0_ma_reg2]  = data_tmp_3; 
			 A_B1R0[bn1_ma_reg3]  = data_tmp_4; 
			 A_B0R0[bn0_ma_reg3]  = data_tmp_5; 
			 A_B0R0[bn0_ma_reg4]  = data_tmp_6; 
			 A_B1R0[bn1_ma_reg4]  = data_tmp_7; 
			 A_B1R0[bn1_ma_reg5]  = data_tmp_8; 
			 A_B0R0[bn0_ma_reg5]  = data_tmp_9; 
			 A_B0R0[bn0_ma_reg6]  = data_tmp_10;
			 A_B1R0[bn1_ma_reg6]  = data_tmp_11;
			 A_B0R0[bn0_ma_reg7]  = data_tmp_12;
			 A_B1R0[bn1_ma_reg7]  = data_tmp_13;
			 A_B1R0[bn1_ma_reg8]  = data_tmp_14;
			 A_B0R0[bn0_ma_reg8]  = data_tmp_15;
			 /*********************************************/
			 data_tmp_1  = A_B1R2[bn1_ma_reg1];
			 data_tmp_2  = A_B1R3[bn1_ma_reg1];
			 data_tmp_3  = A_B1R4[bn1_ma_reg1];
			 data_tmp_4  = A_B1R5[bn1_ma_reg1];
			 data_tmp_5  = A_B1R6[bn1_ma_reg1];
			 data_tmp_6  = A_B1R7[bn1_ma_reg1];
			 data_tmp_7  = A_B1R8[bn1_ma_reg1];
			 data_tmp_8  = A_B1R9[bn1_ma_reg1];
			 data_tmp_9  = A_B1R10[bn1_ma_reg1];
			 data_tmp_10 = A_B1R11[bn1_ma_reg1];
			 data_tmp_11 = A_B1R12[bn1_ma_reg1];
			 data_tmp_12 = A_B1R13[bn1_ma_reg1];
			 data_tmp_13 = A_B1R14[bn1_ma_reg1];
			 data_tmp_14 = A_B1R15[bn1_ma_reg1];
             A_B1R2[bn1_ma_reg1]  =	 A_B1R1[bn1_ma_reg2];
             A_B1R3[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg2];
             A_B1R4[bn1_ma_reg1]  =  A_B1R1[bn1_ma_reg3];
             A_B1R5[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg3];
             A_B1R6[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg4];
             A_B1R7[bn1_ma_reg1]  =  A_B1R1[bn1_ma_reg4];
             A_B1R8[bn1_ma_reg1]  =  A_B1R1[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg1] =  A_B0R1[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg1] =  A_B1R1[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg1] =  A_B0R1[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg1] =  A_B1R1[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg1] =  A_B1R1[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg1] =  A_B0R1[bn0_ma_reg8];
			 A_B1R1[bn1_ma_reg2]  =  data_tmp_1; 
			 A_B0R1[bn0_ma_reg2]  =  data_tmp_2; 
			 A_B1R1[bn1_ma_reg3]  =  data_tmp_3; 
			 A_B0R1[bn0_ma_reg3]  =  data_tmp_4; 
			 A_B0R1[bn0_ma_reg4]  =  data_tmp_5; 
			 A_B1R1[bn1_ma_reg4]  =  data_tmp_6; 
			 A_B1R1[bn1_ma_reg5]  =  data_tmp_7; 
			 A_B0R1[bn0_ma_reg5]  =  data_tmp_8; 
			 A_B0R1[bn0_ma_reg6]  =  data_tmp_9; 
			 A_B1R1[bn1_ma_reg6]  =  data_tmp_10;
			 A_B0R1[bn0_ma_reg7]  =  data_tmp_11;
			 A_B1R1[bn1_ma_reg7]  =  data_tmp_12;
			 A_B1R1[bn1_ma_reg8]  =  data_tmp_13;
			 A_B0R1[bn0_ma_reg8]  =  data_tmp_14;
			/************************************************************/ 
			 data_tmp_1  =   A_B1R3[bn1_ma_reg2];
			 data_tmp_2  =   A_B1R4[bn1_ma_reg2];
			 data_tmp_3  =   A_B1R5[bn1_ma_reg2];
			 data_tmp_4  =   A_B1R6[bn1_ma_reg2];
			 data_tmp_5  =   A_B1R7[bn1_ma_reg2];
			 data_tmp_6  =   A_B1R8[bn1_ma_reg2];
			 data_tmp_7  =   A_B1R9[bn1_ma_reg2];
			 data_tmp_8  =   A_B1R10[bn1_ma_reg2];
			 data_tmp_9  =   A_B1R11[bn1_ma_reg2];
			 data_tmp_10 =   A_B1R12[bn1_ma_reg2];
			 data_tmp_11 =   A_B1R13[bn1_ma_reg2];
			 data_tmp_12 =   A_B1R14[bn1_ma_reg2];
			 data_tmp_13 =   A_B1R15[bn1_ma_reg2];
			 A_B1R3[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg2];
			 A_B1R4[bn1_ma_reg2]  =  A_B1R2[bn1_ma_reg3];
			 A_B1R5[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg3];
			 A_B1R6[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg4];
			 A_B1R7[bn1_ma_reg2]  =  A_B1R2[bn1_ma_reg4];
			 A_B1R8[bn1_ma_reg2]  =  A_B1R2[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg2] =  A_B0R2[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg2] =  A_B1R2[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg2] =  A_B0R2[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg2] =  A_B1R2[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg2] =  A_B1R2[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg2] =  A_B0R2[bn0_ma_reg8];
			 A_B0R2[bn0_ma_reg2]  =  data_tmp_1; 
			 A_B1R2[bn1_ma_reg3]  =  data_tmp_2; 
			 A_B0R2[bn0_ma_reg3]  =  data_tmp_3; 
			 A_B0R2[bn0_ma_reg4]  =  data_tmp_4; 
			 A_B1R2[bn1_ma_reg4]  =  data_tmp_5; 
			 A_B1R2[bn1_ma_reg5]  =  data_tmp_6; 
			 A_B0R2[bn0_ma_reg5]  =  data_tmp_7; 
			 A_B0R2[bn0_ma_reg6]  =  data_tmp_8; 
			 A_B1R2[bn1_ma_reg6]  =  data_tmp_9; 
			 A_B0R2[bn0_ma_reg7]  =  data_tmp_10;
			 A_B1R2[bn1_ma_reg7]  =  data_tmp_11;
			 A_B1R2[bn1_ma_reg8]  =  data_tmp_12;
			 A_B0R2[bn0_ma_reg8]  =  data_tmp_13;
			 //******************************************************
			 data_tmp_1  = A_B0R4[bn0_ma_reg2];
			 data_tmp_2  = A_B0R5[bn0_ma_reg2];
			 data_tmp_3  = A_B0R6[bn0_ma_reg2];
			 data_tmp_4  = A_B0R7[bn0_ma_reg2];
			 data_tmp_5  = A_B0R8[bn0_ma_reg2];
			 data_tmp_6  = A_B0R9[bn0_ma_reg2];
			 data_tmp_7  = A_B0R10[bn0_ma_reg2];
			 data_tmp_8  = A_B0R11[bn0_ma_reg2];
			 data_tmp_9  = A_B0R12[bn0_ma_reg2];
			 data_tmp_10 = A_B0R13[bn0_ma_reg2];
			 data_tmp_11 = A_B0R14[bn0_ma_reg2];
			 data_tmp_12 = A_B0R15[bn0_ma_reg2];
			 A_B0R4[bn0_ma_reg2]  =  A_B1R3[bn1_ma_reg3];
			 A_B0R5[bn0_ma_reg2]  =  A_B0R3[bn0_ma_reg3];
			 A_B0R6[bn0_ma_reg2]  =  A_B0R3[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg2]  =  A_B1R3[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg2]  =  A_B1R3[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg2]  =  A_B0R3[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg2] =  A_B0R3[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg2] =  A_B1R3[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg2] =  A_B0R3[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg2] =  A_B1R3[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg2] =  A_B1R3[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg2] =  A_B0R3[bn0_ma_reg8];
			 A_B1R3[bn1_ma_reg3]  =  data_tmp_1;  
			 A_B0R3[bn0_ma_reg3]  =  data_tmp_2;  
			 A_B0R3[bn0_ma_reg4]  =  data_tmp_3;  
			 A_B1R3[bn1_ma_reg4]  =  data_tmp_4;  
			 A_B1R3[bn1_ma_reg5]  =  data_tmp_5;  
			 A_B0R3[bn0_ma_reg5]  =  data_tmp_6;  
			 A_B0R3[bn0_ma_reg6]  =  data_tmp_7;  
			 A_B1R3[bn1_ma_reg6]  =  data_tmp_8;  
			 A_B0R3[bn0_ma_reg7]  =  data_tmp_9;  
			 A_B1R3[bn1_ma_reg7]  =  data_tmp_10; 
			 A_B1R3[bn1_ma_reg8]  =  data_tmp_11; 
			 A_B0R3[bn0_ma_reg8]  =  data_tmp_12; 
			 //----------------------------------------------------------------------
             data_tmp_1  = A_B1R5[bn1_ma_reg3];
			 data_tmp_2  = A_B1R6[bn1_ma_reg3];
			 data_tmp_3  = A_B1R7[bn1_ma_reg3];
			 data_tmp_4  = A_B1R8[bn1_ma_reg3];
			 data_tmp_5  = A_B1R9[bn1_ma_reg3];
			 data_tmp_6  = A_B1R10[bn1_ma_reg3];
			 data_tmp_7  = A_B1R11[bn1_ma_reg3];
			 data_tmp_8  = A_B1R12[bn1_ma_reg3];
			 data_tmp_9  = A_B1R13[bn1_ma_reg3];
			 data_tmp_10 = A_B1R14[bn1_ma_reg3];
			 data_tmp_11 = A_B1R15[bn1_ma_reg3];
			 A_B1R5[bn1_ma_reg3]  =  A_B0R4[bn0_ma_reg3];
			 A_B1R6[bn1_ma_reg3]  =  A_B0R4[bn0_ma_reg4];
			 A_B1R7[bn1_ma_reg3]  =  A_B1R4[bn1_ma_reg4];
			 A_B1R8[bn1_ma_reg3]  =  A_B1R4[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg3]  =  A_B0R4[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg3] =  A_B0R4[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg3] =  A_B1R4[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg3] =  A_B0R4[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg3] =  A_B1R4[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg3] =  A_B1R4[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg3] =  A_B0R4[bn0_ma_reg8];
			 A_B0R4[bn0_ma_reg3]  =  data_tmp_1 ;
			 A_B0R4[bn0_ma_reg4]  =  data_tmp_2 ;
			 A_B1R4[bn1_ma_reg4]  =  data_tmp_3 ;
			 A_B1R4[bn1_ma_reg5]  =  data_tmp_4 ;
			 A_B0R4[bn0_ma_reg5]  =  data_tmp_5 ;
			 A_B0R4[bn0_ma_reg6]  =  data_tmp_6 ;
			 A_B1R4[bn1_ma_reg6]  =  data_tmp_7 ;
			 A_B0R4[bn0_ma_reg7]  =  data_tmp_8 ;
			 A_B1R4[bn1_ma_reg7]  =  data_tmp_9 ;
			 A_B1R4[bn1_ma_reg8]  =  data_tmp_10;
			 A_B0R4[bn0_ma_reg8]  =  data_tmp_11;
			 //************************************************************************
			 data_tmp_1  = A_B0R6[bn0_ma_reg3];
			 data_tmp_2  = A_B0R7[bn0_ma_reg3];
			 data_tmp_3  = A_B0R8[bn0_ma_reg3];
			 data_tmp_4  = A_B0R9[bn0_ma_reg3];
			 data_tmp_5  = A_B0R10[bn0_ma_reg3];
			 data_tmp_6  = A_B0R11[bn0_ma_reg3];
			 data_tmp_7  = A_B0R12[bn0_ma_reg3];
			 data_tmp_8  = A_B0R13[bn0_ma_reg3];
			 data_tmp_9  = A_B0R14[bn0_ma_reg3];
			 data_tmp_10 = A_B0R15[bn0_ma_reg3];
			 A_B0R6[bn0_ma_reg3]  =  A_B0R5[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg3]  =  A_B1R5[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg3]  =  A_B1R5[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg3]  =  A_B0R5[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg3] =  A_B0R5[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg3] =  A_B1R5[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg3] =  A_B0R5[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg3] =  A_B1R5[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg3] =  A_B1R5[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg3] =  A_B0R5[bn0_ma_reg8];
			 A_B0R5[bn0_ma_reg4]  =  data_tmp_1; 
			 A_B1R5[bn1_ma_reg4]  =  data_tmp_2; 
			 A_B1R5[bn1_ma_reg5]  =  data_tmp_3; 
			 A_B0R5[bn0_ma_reg5]  =  data_tmp_4; 
			 A_B0R5[bn0_ma_reg6]  =  data_tmp_5; 
			 A_B1R5[bn1_ma_reg6]  =  data_tmp_6; 
			 A_B0R5[bn0_ma_reg7]  =  data_tmp_7; 
			 A_B1R5[bn1_ma_reg7]  =  data_tmp_8; 
			 A_B1R5[bn1_ma_reg8]  =  data_tmp_9; 
			 A_B0R5[bn0_ma_reg8]  =  data_tmp_10;
			 //-----------------------------------------------------------------------
			 data_tmp_1  = A_B0R7[bn0_ma_reg4];
			 data_tmp_2  = A_B0R8[bn0_ma_reg4];
			 data_tmp_3  = A_B0R9[bn0_ma_reg4];
			 data_tmp_4  = A_B0R10[bn0_ma_reg4];
			 data_tmp_5  = A_B0R11[bn0_ma_reg4];
			 data_tmp_6  = A_B0R12[bn0_ma_reg4];
			 data_tmp_7  = A_B0R13[bn0_ma_reg4];
			 data_tmp_8  = A_B0R14[bn0_ma_reg4];
			 data_tmp_9  = A_B0R15[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg4]  =  A_B1R6[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg4]  =  A_B1R6[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg4]  =  A_B0R6[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg4] =  A_B0R6[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg4] =  A_B1R6[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg4] =  A_B0R6[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg4] =  A_B1R6[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg4] =  A_B1R6[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg4] =  A_B0R6[bn0_ma_reg8];
			 A_B1R6[bn1_ma_reg4]  =  data_tmp_1;
			 A_B1R6[bn1_ma_reg5]  =  data_tmp_2;
			 A_B0R6[bn0_ma_reg5]  =  data_tmp_3;
			 A_B0R6[bn0_ma_reg6]  =  data_tmp_4;
			 A_B1R6[bn1_ma_reg6]  =  data_tmp_5;
			 A_B0R6[bn0_ma_reg7]  =  data_tmp_6;
			 A_B1R6[bn1_ma_reg7]  =  data_tmp_7;
			 A_B1R6[bn1_ma_reg8]  =  data_tmp_8;
			 A_B0R6[bn0_ma_reg8]  =  data_tmp_9;
			 //----------------------------------------------------------------------
			 data_tmp_1  = A_B1R8[bn1_ma_reg4];
			 data_tmp_2  = A_B1R9[bn1_ma_reg4];
			 data_tmp_3  = A_B1R10[bn1_ma_reg4];
			 data_tmp_4  = A_B1R11[bn1_ma_reg4];
			 data_tmp_5  = A_B1R12[bn1_ma_reg4];
			 data_tmp_6  = A_B1R13[bn1_ma_reg4];
			 data_tmp_7  = A_B1R14[bn1_ma_reg4];
			 data_tmp_8  = A_B1R15[bn1_ma_reg4];
			 A_B1R8[bn1_ma_reg4]  =  A_B1R7[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg4]  =  A_B0R7[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg4] =  A_B0R7[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg4] =  A_B1R7[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg4] =  A_B0R7[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg4] =  A_B1R7[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg4] =  A_B1R7[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg4] =  A_B0R7[bn0_ma_reg8];
			 A_B1R7[bn1_ma_reg5]  =  data_tmp_1;
			 A_B0R7[bn0_ma_reg5]  =  data_tmp_2;
			 A_B0R7[bn0_ma_reg6]  =  data_tmp_3;
			 A_B1R7[bn1_ma_reg6]  =  data_tmp_4;
			 A_B0R7[bn0_ma_reg7]  =  data_tmp_5;
			 A_B1R7[bn1_ma_reg7]  =  data_tmp_6;
			 A_B1R7[bn1_ma_reg8]  =  data_tmp_7;
			 A_B0R7[bn0_ma_reg8]  =  data_tmp_8;
			 //----------------------------------------------------------------------
			 data_tmp_1 = A_B1R9[bn1_ma_reg5];
			 data_tmp_2 = A_B1R10[bn1_ma_reg5];
			 data_tmp_3 = A_B1R11[bn1_ma_reg5];
			 data_tmp_4 = A_B1R12[bn1_ma_reg5];
			 data_tmp_5 = A_B1R13[bn1_ma_reg5];
			 data_tmp_6 = A_B1R14[bn1_ma_reg5];
			 data_tmp_7 = A_B1R15[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg5]  = A_B0R8[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg5] = A_B0R8[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg5] = A_B1R8[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg5] = A_B0R8[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg5] = A_B1R8[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg5] = A_B1R8[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg5] = A_B0R8[bn0_ma_reg8];
			 A_B0R8[bn0_ma_reg5]  = data_tmp_1;
			 A_B0R8[bn0_ma_reg6]  = data_tmp_2;
			 A_B1R8[bn1_ma_reg6]  = data_tmp_3;
			 A_B0R8[bn0_ma_reg7]  = data_tmp_4;
			 A_B1R8[bn1_ma_reg7]  = data_tmp_5;
			 A_B1R8[bn1_ma_reg8]  = data_tmp_6;
			 A_B0R8[bn0_ma_reg8]  = data_tmp_7;
			 //---------------------------------------------------------------------
			 data_tmp_1  = A_B0R10[bn0_ma_reg5];
			 data_tmp_2  = A_B0R11[bn0_ma_reg5];
			 data_tmp_3  = A_B0R12[bn0_ma_reg5];
			 data_tmp_4  = A_B0R13[bn0_ma_reg5];
			 data_tmp_5  = A_B0R14[bn0_ma_reg5];
			 data_tmp_6  = A_B0R15[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg5] = A_B0R9[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg5] = A_B1R9[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg5] = A_B0R9[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg5] = A_B1R9[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg5] = A_B1R9[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg5] = A_B0R9[bn0_ma_reg8];
			 A_B0R9[bn0_ma_reg6]  = data_tmp_1;
			 A_B1R9[bn1_ma_reg6]  = data_tmp_2;
			 A_B0R9[bn0_ma_reg7]  = data_tmp_3;
			 A_B1R9[bn1_ma_reg7]  = data_tmp_4;
			 A_B1R9[bn1_ma_reg8]  = data_tmp_5;
			 A_B0R9[bn0_ma_reg8]  = data_tmp_6;
			 //--------------------------------------------------------------------
			 data_tmp_1  = A_B0R11[bn0_ma_reg6];
			 data_tmp_2  = A_B0R12[bn0_ma_reg6];
			 data_tmp_3  = A_B0R13[bn0_ma_reg6];
			 data_tmp_4  = A_B0R14[bn0_ma_reg6];
			 data_tmp_5  = A_B0R15[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg6] = A_B1R10[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg6] = A_B0R10[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg6] = A_B1R10[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg6] = A_B1R10[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg6] = A_B0R10[bn0_ma_reg8];
			 A_B1R10[bn1_ma_reg6] = data_tmp_1;
			 A_B0R10[bn0_ma_reg7] = data_tmp_2;
			 A_B1R10[bn1_ma_reg7] = data_tmp_3;
			 A_B1R10[bn1_ma_reg8] = data_tmp_4;
			 A_B0R10[bn0_ma_reg8] = data_tmp_5;
			 //--------------------------------------------------------------------
			 data_tmp_1  = A_B1R12[bn1_ma_reg6];
			 data_tmp_2  = A_B1R13[bn1_ma_reg6];
			 data_tmp_3  = A_B1R14[bn1_ma_reg6];
			 data_tmp_4  = A_B1R15[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg6] = A_B0R11[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg6] = A_B1R11[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg6] = A_B1R11[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg6] = A_B0R11[bn0_ma_reg8];
			 A_B0R11[bn0_ma_reg7] = data_tmp_1;
			 A_B1R11[bn1_ma_reg7] = data_tmp_2;
			 A_B1R11[bn1_ma_reg8] = data_tmp_3;
			 A_B0R11[bn0_ma_reg8] = data_tmp_4;
			 //--------------------------------------------------------------------
			 data_tmp_1 = A_B0R13[bn0_ma_reg7];
			 data_tmp_2 = A_B0R14[bn0_ma_reg7];
			 data_tmp_3 = A_B0R15[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg7] = A_B1R12[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg7] = A_B1R12[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg7] = A_B0R12[bn0_ma_reg8];
			 A_B1R12[bn1_ma_reg7] = data_tmp_1;
			 A_B1R12[bn1_ma_reg8] = data_tmp_2;
			 A_B0R12[bn0_ma_reg8] = data_tmp_3;
			 //--------------------------------------------------------------------
			 data_tmp_1 = A_B1R14[bn1_ma_reg7];
			 data_tmp_2 = A_B1R15[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg7] = A_B1R13[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg7] = A_B0R13[bn0_ma_reg8];
			 A_B1R13[bn1_ma_reg8] = data_tmp_1;
			 A_B0R13[bn0_ma_reg8] = data_tmp_2;
			 //--------------------------------------------------------------------
			 data_tmp_1 = A_B1R15[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg8] = A_B0R14[bn0_ma_reg8];
			 A_B0R14[bn0_ma_reg8] = data_tmp_1; 
		  }else{	 
			 data_tmp_1  = A_B1R1[bn1_ma_reg1];
			 data_tmp_2  = A_B1R2[bn1_ma_reg1];
			 data_tmp_3  = A_B1R3[bn1_ma_reg1];
			 data_tmp_4  = A_B1R4[bn1_ma_reg1];
			 data_tmp_5  = A_B1R5[bn1_ma_reg1];
			 data_tmp_6  = A_B1R6[bn1_ma_reg1];
			 data_tmp_7  = A_B1R7[bn1_ma_reg1];
			 data_tmp_8  = A_B1R8[bn1_ma_reg1];
			 data_tmp_9  = A_B1R9[bn1_ma_reg1];
			 data_tmp_10 = A_B1R10[bn1_ma_reg1];
			 data_tmp_11 = A_B1R11[bn1_ma_reg1];
			 data_tmp_12 = A_B1R12[bn1_ma_reg1];
			 data_tmp_13 = A_B1R13[bn1_ma_reg1];
			 data_tmp_14 = A_B1R14[bn1_ma_reg1];
			 data_tmp_15 = A_B1R15[bn1_ma_reg1];
             A_B1R1[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg1];
             A_B1R2[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg2];
             A_B1R3[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg2];
             A_B1R4[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg3];
             A_B1R5[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg3];
             A_B1R6[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg4];
             A_B1R7[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg4];
             A_B1R8[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg5];
             A_B1R9[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg5];
             A_B1R10[bn1_ma_reg1] =  A_B1R0[bn1_ma_reg6];
             A_B1R11[bn1_ma_reg1] =  A_B0R0[bn0_ma_reg6];
             A_B1R12[bn1_ma_reg1] =  A_B1R0[bn1_ma_reg7];
             A_B1R13[bn1_ma_reg1] =  A_B0R0[bn0_ma_reg7];
             A_B1R14[bn1_ma_reg1] =  A_B0R0[bn0_ma_reg8];
             A_B1R15[bn1_ma_reg1] =  A_B1R0[bn1_ma_reg8];
             A_B0R0[bn0_ma_reg1]  =  data_tmp_1; 
             A_B0R0[bn0_ma_reg2]  =  data_tmp_2; 
             A_B1R0[bn1_ma_reg2]  =  data_tmp_3; 
             A_B0R0[bn0_ma_reg3]  =  data_tmp_4; 
             A_B1R0[bn1_ma_reg3]  =  data_tmp_5; 
             A_B1R0[bn1_ma_reg4]  =  data_tmp_6; 
             A_B0R0[bn0_ma_reg4]  =  data_tmp_7; 
             A_B0R0[bn0_ma_reg5]  =  data_tmp_8; 
             A_B1R0[bn1_ma_reg5]  =  data_tmp_9; 
             A_B1R0[bn1_ma_reg6]  =  data_tmp_10;
             A_B0R0[bn0_ma_reg6]  =  data_tmp_11;
             A_B1R0[bn1_ma_reg7]  =  data_tmp_12;
             A_B0R0[bn0_ma_reg7]  =  data_tmp_13;
             A_B0R0[bn0_ma_reg8]  =  data_tmp_14;
             A_B1R0[bn1_ma_reg8]  =  data_tmp_15;
             //-----------------------------------------------------------
             data_tmp_1  = A_B0R2[bn0_ma_reg1];
			 data_tmp_2  = A_B0R3[bn0_ma_reg1];
			 data_tmp_3  = A_B0R4[bn0_ma_reg1];
			 data_tmp_4  = A_B0R5[bn0_ma_reg1];
			 data_tmp_5  = A_B0R6[bn0_ma_reg1];
			 data_tmp_6  = A_B0R7[bn0_ma_reg1];
			 data_tmp_7  = A_B0R8[bn0_ma_reg1];
			 data_tmp_8  = A_B0R9[bn0_ma_reg1];
			 data_tmp_9  = A_B0R10[bn0_ma_reg1];
			 data_tmp_10 = A_B0R11[bn0_ma_reg1];
			 data_tmp_11 = A_B0R12[bn0_ma_reg1];
			 data_tmp_12 = A_B0R13[bn0_ma_reg1];
			 data_tmp_13 = A_B0R14[bn0_ma_reg1];
			 data_tmp_14 = A_B0R15[bn0_ma_reg1];
			 A_B0R2[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg2];
			 A_B0R3[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg2];
			 A_B0R4[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg3];
			 A_B0R5[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg3];
			 A_B0R6[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg4];
			 A_B0R7[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg1] =  A_B1R1[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg1] =  A_B0R1[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg1] =  A_B1R1[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg1] =  A_B0R1[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg1] =  A_B0R1[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg1] =  A_B1R1[bn1_ma_reg8];
             A_B0R1[bn0_ma_reg2]  =  data_tmp_1; 
             A_B1R1[bn1_ma_reg2]  =  data_tmp_2; 
             A_B0R1[bn0_ma_reg3]  =  data_tmp_3; 
             A_B1R1[bn1_ma_reg3]  =  data_tmp_4; 
             A_B1R1[bn1_ma_reg4]  =  data_tmp_5; 
             A_B0R1[bn0_ma_reg4]  =  data_tmp_6; 
             A_B0R1[bn0_ma_reg5]  =  data_tmp_7; 
             A_B1R1[bn1_ma_reg5]  =  data_tmp_8; 
             A_B1R1[bn1_ma_reg6]  =  data_tmp_9; 
             A_B0R1[bn0_ma_reg6]  =  data_tmp_10;
             A_B1R1[bn1_ma_reg7]  =  data_tmp_11;
             A_B0R1[bn0_ma_reg7]  =  data_tmp_12;
             A_B0R1[bn0_ma_reg8]  =  data_tmp_13;
             A_B1R1[bn1_ma_reg8]  =  data_tmp_14;
             //------------------------------------------------------------
             data_tmp_1  =  A_B0R3[bn0_ma_reg2];
			 data_tmp_2  =  A_B0R4[bn0_ma_reg2];
			 data_tmp_3  =  A_B0R5[bn0_ma_reg2];
			 data_tmp_4  =  A_B0R6[bn0_ma_reg2];
			 data_tmp_5  =  A_B0R7[bn0_ma_reg2];
			 data_tmp_6  =  A_B0R8[bn0_ma_reg2];
			 data_tmp_7  =  A_B0R9[bn0_ma_reg2];
			 data_tmp_8  =  A_B0R10[bn0_ma_reg2];
			 data_tmp_9  =  A_B0R11[bn0_ma_reg2];
			 data_tmp_10 =  A_B0R12[bn0_ma_reg2];
			 data_tmp_11 =  A_B0R13[bn0_ma_reg2];
			 data_tmp_12 =  A_B0R14[bn0_ma_reg2];
			 data_tmp_13 =  A_B0R15[bn0_ma_reg2];
			 A_B0R3[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg2];
			 A_B0R4[bn0_ma_reg2]  =  A_B0R2[bn0_ma_reg3];
			 A_B0R5[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg3];
			 A_B0R6[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg4];
			 A_B0R7[bn0_ma_reg2]  =  A_B0R2[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg2]  =  A_B0R2[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg2] =  A_B1R2[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg2] =  A_B0R2[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg2] =  A_B1R2[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg2] =  A_B0R2[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg2] =  A_B0R2[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg2] =  A_B1R2[bn1_ma_reg8];
			 A_B1R2[bn1_ma_reg2]  =  data_tmp_1; 
			 A_B0R2[bn0_ma_reg3]  =  data_tmp_2; 
			 A_B1R2[bn1_ma_reg3]  =  data_tmp_3; 
			 A_B1R2[bn1_ma_reg4]  =  data_tmp_4; 
			 A_B0R2[bn0_ma_reg4]  =  data_tmp_5; 
			 A_B0R2[bn0_ma_reg5]  =  data_tmp_6; 
			 A_B1R2[bn1_ma_reg5]  =  data_tmp_7; 
			 A_B1R2[bn1_ma_reg6]  =  data_tmp_8; 
			 A_B0R2[bn0_ma_reg6]  =  data_tmp_9; 
			 A_B1R2[bn1_ma_reg7]  =  data_tmp_10;
			 A_B0R2[bn0_ma_reg7]  =  data_tmp_11;
			 A_B0R2[bn0_ma_reg8]  =  data_tmp_12;
			 A_B1R2[bn1_ma_reg8]  =  data_tmp_13;
			 //-----------------------------------------------------------                        
			 data_tmp_1  = A_B1R4[bn1_ma_reg2];
			 data_tmp_2  = A_B1R5[bn1_ma_reg2];
			 data_tmp_3  = A_B1R6[bn1_ma_reg2];
			 data_tmp_4  = A_B1R7[bn1_ma_reg2];
			 data_tmp_5  = A_B1R8[bn1_ma_reg2];
			 data_tmp_6  = A_B1R9[bn1_ma_reg2];
			 data_tmp_7  = A_B1R10[bn1_ma_reg2];
			 data_tmp_8  = A_B1R11[bn1_ma_reg2];
			 data_tmp_9  = A_B1R12[bn1_ma_reg2];
			 data_tmp_10 = A_B1R13[bn1_ma_reg2];
			 data_tmp_11 = A_B1R14[bn1_ma_reg2];
			 data_tmp_12 = A_B1R15[bn1_ma_reg2];
			 A_B1R4[bn1_ma_reg2]  = A_B0R3[bn0_ma_reg3];
			 A_B1R5[bn1_ma_reg2]  = A_B1R3[bn1_ma_reg3];
			 A_B1R6[bn1_ma_reg2]  = A_B1R3[bn1_ma_reg4];
			 A_B1R7[bn1_ma_reg2]  = A_B0R3[bn0_ma_reg4];
			 A_B1R8[bn1_ma_reg2]  = A_B0R3[bn0_ma_reg5];
			 A_B1R9[bn1_ma_reg2]  = A_B1R3[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg2] = A_B1R3[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg2] = A_B0R3[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg2] = A_B1R3[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg2] = A_B0R3[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg2] = A_B0R3[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg2] = A_B1R3[bn1_ma_reg8];
			 A_B0R3[bn0_ma_reg3]  = data_tmp_1; 
			 A_B1R3[bn1_ma_reg3]  = data_tmp_2; 
			 A_B1R3[bn1_ma_reg4]  = data_tmp_3; 
			 A_B0R3[bn0_ma_reg4]  = data_tmp_4; 
			 A_B0R3[bn0_ma_reg5]  = data_tmp_5; 
			 A_B1R3[bn1_ma_reg5]  = data_tmp_6; 
			 A_B1R3[bn1_ma_reg6]  = data_tmp_7; 
			 A_B0R3[bn0_ma_reg6]  = data_tmp_8; 
			 A_B1R3[bn1_ma_reg7]  = data_tmp_9; 
			 A_B0R3[bn0_ma_reg7]  = data_tmp_10;
			 A_B0R3[bn0_ma_reg8]  = data_tmp_11;
			 A_B1R3[bn1_ma_reg8]  = data_tmp_12;
			 //------------------------------------------------------------
			 data_tmp_1  =  A_B0R5[bn0_ma_reg3];
			 data_tmp_2  =  A_B0R6[bn0_ma_reg3];
			 data_tmp_3  =  A_B0R7[bn0_ma_reg3];
			 data_tmp_4  =  A_B0R8[bn0_ma_reg3];
			 data_tmp_5  =  A_B0R9[bn0_ma_reg3];
			 data_tmp_6  =  A_B0R10[bn0_ma_reg3];
			 data_tmp_7  =  A_B0R11[bn0_ma_reg3];
			 data_tmp_8  =  A_B0R12[bn0_ma_reg3];
			 data_tmp_9  =  A_B0R13[bn0_ma_reg3];
			 data_tmp_10 =  A_B0R14[bn0_ma_reg3];
			 data_tmp_11 =  A_B0R15[bn0_ma_reg3];
			 A_B0R5[bn0_ma_reg3]  =  A_B1R4[bn1_ma_reg3];
			 A_B0R6[bn0_ma_reg3]  =  A_B1R4[bn1_ma_reg4];
			 A_B0R7[bn0_ma_reg3]  =  A_B0R4[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg3]  =  A_B0R4[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg3]  =  A_B1R4[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg3] =  A_B1R4[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg3] =  A_B0R4[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg3] =  A_B1R4[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg3] =  A_B0R4[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg3] =  A_B0R4[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg3] =  A_B1R4[bn1_ma_reg8];
			 A_B1R4[bn1_ma_reg3]  =  data_tmp_1; 
			 A_B1R4[bn1_ma_reg4]  =  data_tmp_2; 
			 A_B0R4[bn0_ma_reg4]  =  data_tmp_3; 
			 A_B0R4[bn0_ma_reg5]  =  data_tmp_4; 
			 A_B1R4[bn1_ma_reg5]  =  data_tmp_5; 
			 A_B1R4[bn1_ma_reg6]  =  data_tmp_6; 
			 A_B0R4[bn0_ma_reg6]  =  data_tmp_7; 
			 A_B1R4[bn1_ma_reg7]  =  data_tmp_8; 
			 A_B0R4[bn0_ma_reg7]  =  data_tmp_9; 
			 A_B0R4[bn0_ma_reg8]  =  data_tmp_10;
			 A_B1R4[bn1_ma_reg8]  =  data_tmp_11;
			 //-------------------------------------------------------------
			 data_tmp_1  =  A_B1R6[bn1_ma_reg3];
			 data_tmp_2  =  A_B1R7[bn1_ma_reg3];
			 data_tmp_3  =  A_B1R8[bn1_ma_reg3];
			 data_tmp_4  =  A_B1R9[bn1_ma_reg3];
			 data_tmp_5  =  A_B1R10[bn1_ma_reg3];
			 data_tmp_6  =  A_B1R11[bn1_ma_reg3];
			 data_tmp_7  =  A_B1R12[bn1_ma_reg3];
			 data_tmp_8  =  A_B1R13[bn1_ma_reg3];
			 data_tmp_9  =  A_B1R14[bn1_ma_reg3];
			 data_tmp_10 =  A_B1R15[bn1_ma_reg3];
			 A_B1R6[bn1_ma_reg3]  =  A_B1R5[bn1_ma_reg4];
			 A_B1R7[bn1_ma_reg3]  =  A_B0R5[bn0_ma_reg4];
			 A_B1R8[bn1_ma_reg3]  =  A_B0R5[bn0_ma_reg5];
			 A_B1R9[bn1_ma_reg3]  =  A_B1R5[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg3] =  A_B1R5[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg3] =  A_B0R5[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg3] =  A_B1R5[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg3] =  A_B0R5[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg3] =  A_B0R5[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg3] =  A_B1R5[bn1_ma_reg8];
			 A_B1R5[bn1_ma_reg4]  =  data_tmp_1; 
			 A_B0R5[bn0_ma_reg4]  =  data_tmp_2; 
			 A_B0R5[bn0_ma_reg5]  =  data_tmp_3; 
			 A_B1R5[bn1_ma_reg5]  =  data_tmp_4; 
			 A_B1R5[bn1_ma_reg6]  =  data_tmp_5; 
			 A_B0R5[bn0_ma_reg6]  =  data_tmp_6; 
			 A_B1R5[bn1_ma_reg7]  =  data_tmp_7; 
			 A_B0R5[bn0_ma_reg7]  =  data_tmp_8; 
			 A_B0R5[bn0_ma_reg8]  =  data_tmp_9; 
			 A_B1R5[bn1_ma_reg8]  =  data_tmp_10;
			 //----------------------------------------------------------------
			 data_tmp_1  =  A_B1R7[bn1_ma_reg4];
			 data_tmp_2  =  A_B1R8[bn1_ma_reg4];
			 data_tmp_3  =  A_B1R9[bn1_ma_reg4];
			 data_tmp_4  =  A_B1R10[bn1_ma_reg4];
			 data_tmp_5  =  A_B1R11[bn1_ma_reg4];
			 data_tmp_6  =  A_B1R12[bn1_ma_reg4];
			 data_tmp_7  =  A_B1R13[bn1_ma_reg4];
			 data_tmp_8  =  A_B1R14[bn1_ma_reg4];
			 data_tmp_9  =  A_B1R15[bn1_ma_reg4];
			 A_B1R7[bn1_ma_reg4]  =  A_B0R6[bn0_ma_reg4];
			 A_B1R8[bn1_ma_reg4]  =  A_B0R6[bn0_ma_reg5];
			 A_B1R9[bn1_ma_reg4]  =  A_B1R6[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg4] =  A_B1R6[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg4] =  A_B0R6[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg4] =  A_B1R6[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg4] =  A_B0R6[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg4] =  A_B0R6[bn0_ma_reg8];
             A_B1R15[bn1_ma_reg4] =  A_B1R6[bn1_ma_reg8];
             A_B0R6[bn0_ma_reg4]  =  data_tmp_1;
             A_B0R6[bn0_ma_reg5]  =  data_tmp_2;
             A_B1R6[bn1_ma_reg5]  =  data_tmp_3;
             A_B1R6[bn1_ma_reg6]  =  data_tmp_4;
             A_B0R6[bn0_ma_reg6]  =  data_tmp_5;
             A_B1R6[bn1_ma_reg7]  =  data_tmp_6;
             A_B0R6[bn0_ma_reg7]  =  data_tmp_7;
             A_B0R6[bn0_ma_reg8]  =  data_tmp_8;
             A_B1R6[bn1_ma_reg8]  =  data_tmp_9;
             //------------------------------------------------------------------
			 data_tmp_1  = A_B0R8[bn0_ma_reg4];
			 data_tmp_2  = A_B0R9[bn0_ma_reg4];
			 data_tmp_3  = A_B0R10[bn0_ma_reg4];
			 data_tmp_4  = A_B0R11[bn0_ma_reg4];
			 data_tmp_5  = A_B0R12[bn0_ma_reg4];
			 data_tmp_6  = A_B0R13[bn0_ma_reg4];
			 data_tmp_7  = A_B0R14[bn0_ma_reg4];
			 data_tmp_8  = A_B0R15[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg4]  =  A_B0R7[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg4]  =  A_B1R7[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg4] =  A_B1R7[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg4] =  A_B0R7[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg4] =  A_B1R7[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg4] =  A_B0R7[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg4] =  A_B0R7[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg4] =  A_B1R7[bn1_ma_reg8];
			 A_B0R7[bn0_ma_reg5]  =  data_tmp_1;
			 A_B1R7[bn1_ma_reg5]  =  data_tmp_2;
			 A_B1R7[bn1_ma_reg6]  =  data_tmp_3;
			 A_B0R7[bn0_ma_reg6]  =  data_tmp_4;
			 A_B1R7[bn1_ma_reg7]  =  data_tmp_5;
			 A_B0R7[bn0_ma_reg7]  =  data_tmp_6;
			 A_B0R7[bn0_ma_reg8]  =  data_tmp_7;
			 A_B1R7[bn1_ma_reg8]  =  data_tmp_8;
			 //----------------------------------------------------------------
			 data_tmp_1  = A_B0R9[bn0_ma_reg5];
			 data_tmp_2  = A_B0R10[bn0_ma_reg5];
			 data_tmp_3  = A_B0R11[bn0_ma_reg5];
			 data_tmp_4  = A_B0R12[bn0_ma_reg5];
			 data_tmp_5  = A_B0R13[bn0_ma_reg5];
			 data_tmp_6  = A_B0R14[bn0_ma_reg5];
			 data_tmp_7  = A_B0R15[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg5]  =  A_B1R8[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg8];
			 A_B1R8[bn1_ma_reg5]  =  data_tmp_1;
			 A_B1R8[bn1_ma_reg6]  =  data_tmp_2;
			 A_B0R8[bn0_ma_reg6]  =  data_tmp_3;
			 A_B1R8[bn1_ma_reg7]  =  data_tmp_4;
			 A_B0R8[bn0_ma_reg7]  =  data_tmp_5;
			 A_B0R8[bn0_ma_reg8]  =  data_tmp_6;
			 A_B1R8[bn1_ma_reg8]  =  data_tmp_7;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B1R10[bn1_ma_reg5];
			 data_tmp_2 = A_B1R11[bn1_ma_reg5];
			 data_tmp_3 = A_B1R12[bn1_ma_reg5];
			 data_tmp_4 = A_B1R13[bn1_ma_reg5];
			 data_tmp_5 = A_B1R14[bn1_ma_reg5];
			 data_tmp_6 = A_B1R15[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg5] = A_B1R9[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg5] = A_B0R9[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg5] = A_B1R9[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg5] = A_B0R9[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg5] = A_B0R9[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg5] = A_B1R9[bn1_ma_reg8];
			 A_B1R9[bn1_ma_reg6]  = data_tmp_1;
			 A_B0R9[bn0_ma_reg6]  = data_tmp_2;
			 A_B1R9[bn1_ma_reg7]  = data_tmp_3;
			 A_B0R9[bn0_ma_reg7]  = data_tmp_4;
			 A_B0R9[bn0_ma_reg8]  = data_tmp_5;
			 A_B1R9[bn1_ma_reg8]  = data_tmp_6;
			 //----------------------------------------------------------------
			 data_tmp_1 =  A_B1R11[bn1_ma_reg6];
			 data_tmp_2 =  A_B1R12[bn1_ma_reg6];
			 data_tmp_3 =  A_B1R13[bn1_ma_reg6];
			 data_tmp_4 =  A_B1R14[bn1_ma_reg6];
			 data_tmp_5 =  A_B1R15[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg6]  = A_B0R10[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg6]  = A_B1R10[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg6]  = A_B0R10[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg6]  = A_B0R10[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg6]  = A_B1R10[bn1_ma_reg8];
			 A_B0R10[bn0_ma_reg6]  = data_tmp_1;
			 A_B1R10[bn1_ma_reg7]  = data_tmp_2;
			 A_B0R10[bn0_ma_reg7]  = data_tmp_3;
			 A_B0R10[bn0_ma_reg8]  = data_tmp_4;
			 A_B1R10[bn1_ma_reg8]  = data_tmp_5;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B0R12[bn0_ma_reg6];
			 data_tmp_2 = A_B0R13[bn0_ma_reg6];
			 data_tmp_3 = A_B0R14[bn0_ma_reg6];
			 data_tmp_4 = A_B0R15[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg6] = A_B1R11[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg6] = A_B0R11[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg6] = A_B0R11[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg6] = A_B1R11[bn1_ma_reg8];
			 A_B1R11[bn1_ma_reg7] = data_tmp_1;
			 A_B0R11[bn0_ma_reg7] = data_tmp_2;
			 A_B0R11[bn0_ma_reg8] = data_tmp_3;
			 A_B1R11[bn1_ma_reg8] = data_tmp_4;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B1R13[bn1_ma_reg7];
			 data_tmp_2 = A_B1R14[bn1_ma_reg7];
			 data_tmp_3 = A_B1R15[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg7] = A_B0R12[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg7] = A_B0R12[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg7] = A_B1R12[bn1_ma_reg8];
			 A_B0R12[bn0_ma_reg7] = data_tmp_1;
			 A_B0R12[bn0_ma_reg8] = data_tmp_2;
			 A_B1R12[bn1_ma_reg8] = data_tmp_3;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B0R14[bn0_ma_reg7];
			 data_tmp_2 = A_B0R15[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg7] = A_B0R13[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg7] = A_B1R13[bn1_ma_reg8];
			 A_B0R13[bn0_ma_reg8] = data_tmp_1;
			 A_B1R13[bn1_ma_reg8] = data_tmp_2;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B0R15[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg8] = A_B1R14[bn1_ma_reg8];
			 A_B1R14[bn1_ma_reg8] = data_tmp_1;
			 //-----------------------------------------------------------------
		  }
		 }
		 else{
		  //radix-16 final stage data relocation
		  if(bn1_bc_tmp > bn0_bc_tmp){
			//bn: 0 1 1 0
		    data_tmp_1 = A_B0R1[bn0_ma_reg1];
		    data_tmp_2 = A_B0R2[bn0_ma_reg1];
		    data_tmp_3 = A_B0R3[bn0_ma_reg1];
			A_B0R1[bn0_ma_reg1] = A_B1R0[bn1_ma_reg1];
			A_B0R2[bn0_ma_reg1] = A_B1R0[bn1_ma_reg2];
			A_B0R3[bn0_ma_reg1] = A_B0R0[bn0_ma_reg2];
			A_B1R0[bn1_ma_reg1] = data_tmp_1;
			A_B1R0[bn1_ma_reg2] = data_tmp_2;
			A_B0R0[bn0_ma_reg2] = data_tmp_3;
			data_tmp_1 = A_B1R2[bn1_ma_reg1];
			data_tmp_2 = A_B1R3[bn1_ma_reg1];
			A_B1R2[bn1_ma_reg1] = A_B1R1[bn1_ma_reg2];
			A_B1R3[bn1_ma_reg1] = A_B0R1[bn0_ma_reg2];
			A_B1R1[bn1_ma_reg2] = data_tmp_1;
			A_B0R1[bn0_ma_reg2] = data_tmp_2;
			data_tmp_1 = A_B1R3[bn1_ma_reg2];
			A_B1R3[bn1_ma_reg2] = A_B0R2[bn0_ma_reg2];
			A_B0R2[bn0_ma_reg2] = data_tmp_1;
			//-----------------------------------------
			data_tmp_1 = A_B0R5[bn0_ma_reg1];
			data_tmp_2 = A_B0R6[bn0_ma_reg1];
			data_tmp_3 = A_B0R7[bn0_ma_reg1];
			A_B0R5[bn0_ma_reg1] = A_B1R4[bn1_ma_reg1];
			A_B0R6[bn0_ma_reg1] = A_B1R4[bn1_ma_reg2];
			A_B0R7[bn0_ma_reg1] = A_B0R4[bn0_ma_reg2];
			A_B1R4[bn1_ma_reg1] = data_tmp_1;
			A_B1R4[bn1_ma_reg2] = data_tmp_2;
			A_B0R4[bn0_ma_reg2] = data_tmp_3;
			data_tmp_1 = A_B1R6[bn1_ma_reg1];
			data_tmp_2 = A_B1R7[bn1_ma_reg1];
			A_B1R6[bn1_ma_reg1] = A_B1R5[bn1_ma_reg2];
			A_B1R7[bn1_ma_reg1] = A_B0R5[bn0_ma_reg2];
			A_B1R5[bn1_ma_reg2] = data_tmp_1;
			A_B0R5[bn0_ma_reg2] = data_tmp_2;
			data_tmp_1 = A_B1R7[bn1_ma_reg2];
			A_B1R7[bn1_ma_reg2] = A_B0R6[bn0_ma_reg2];
			A_B0R6[bn0_ma_reg2] = data_tmp_1;
			//----------------------------------------
			data_tmp_1 = A_B0R9[bn0_ma_reg1];
			data_tmp_2 = A_B0R10[bn0_ma_reg1];
			data_tmp_3 = A_B0R11[bn0_ma_reg1];
			A_B0R9[bn0_ma_reg1]  = A_B1R8[bn1_ma_reg1];
			A_B0R10[bn0_ma_reg1] = A_B1R8[bn1_ma_reg2];
			A_B0R11[bn0_ma_reg1] = A_B0R8[bn0_ma_reg2];
			A_B1R8[bn1_ma_reg1] = data_tmp_1;
			A_B1R8[bn1_ma_reg2] = data_tmp_2;
			A_B0R8[bn0_ma_reg2] = data_tmp_3;
			data_tmp_1 = A_B1R10[bn1_ma_reg1];
			data_tmp_2 = A_B1R11[bn1_ma_reg1];
			A_B1R10[bn1_ma_reg1] = A_B1R9[bn1_ma_reg2];
			A_B1R11[bn1_ma_reg1] = A_B0R9[bn0_ma_reg2];
			A_B1R9[bn1_ma_reg2] = data_tmp_1;
			A_B0R9[bn0_ma_reg2] = data_tmp_2;
			data_tmp_1 = A_B1R11[bn1_ma_reg2];
			A_B1R11[bn1_ma_reg2] = A_B0R10[bn0_ma_reg2];
			A_B0R10[bn0_ma_reg2] = data_tmp_1;			
			//----------------------------------------
			data_tmp_1 = A_B0R13[bn0_ma_reg1];
			data_tmp_2 = A_B0R14[bn0_ma_reg1];
			data_tmp_3 = A_B0R15[bn0_ma_reg1];
			A_B0R13[bn0_ma_reg1] = A_B1R12[bn1_ma_reg1];
			A_B0R14[bn0_ma_reg1] = A_B1R12[bn1_ma_reg2];
			A_B0R15[bn0_ma_reg1] = A_B0R12[bn0_ma_reg2];
			A_B1R12[bn1_ma_reg1] = data_tmp_1;
			A_B1R12[bn1_ma_reg2] = data_tmp_2;
			A_B0R12[bn0_ma_reg2] = data_tmp_3;
			data_tmp_1 = A_B1R14[bn1_ma_reg1];
			data_tmp_2 = A_B1R15[bn1_ma_reg1];
			A_B1R14[bn1_ma_reg1] = A_B1R13[bn1_ma_reg2];
			A_B1R15[bn1_ma_reg1] = A_B0R13[bn0_ma_reg2];
			A_B1R13[bn1_ma_reg2] = data_tmp_1;
			A_B0R13[bn0_ma_reg2] = data_tmp_2;
			data_tmp_1 = A_B1R15[bn1_ma_reg2];
			A_B1R15[bn1_ma_reg2] = A_B0R14[bn0_ma_reg2];
			A_B0R14[bn0_ma_reg2] = data_tmp_1;				
			//******************************************************
			//bn 1 0 0 1
		    data_tmp_1 = A_B1R1[bn1_ma_reg3];
		    data_tmp_2 = A_B1R2[bn1_ma_reg3];
		    data_tmp_3 = A_B1R3[bn1_ma_reg3];
			A_B1R1[bn1_ma_reg3] = A_B0R0[bn0_ma_reg3];
			A_B1R2[bn1_ma_reg3] = A_B0R0[bn0_ma_reg4];
			A_B1R3[bn1_ma_reg3] = A_B1R0[bn1_ma_reg4];
			A_B0R0[bn0_ma_reg3] = data_tmp_1;
			A_B0R0[bn0_ma_reg4] = data_tmp_2;
			A_B1R0[bn1_ma_reg4] = data_tmp_3;
			data_tmp_1 = A_B0R2[bn0_ma_reg3];
			data_tmp_2 = A_B0R3[bn0_ma_reg3];
			A_B0R2[bn0_ma_reg3] = A_B0R1[bn0_ma_reg4];
			A_B0R3[bn0_ma_reg3] = A_B1R1[bn1_ma_reg4];
			A_B0R1[bn0_ma_reg4] = data_tmp_1;
			A_B1R1[bn1_ma_reg4] = data_tmp_2;
			data_tmp_1 = A_B0R3[bn0_ma_reg4];
			A_B0R3[bn0_ma_reg4] = A_B1R2[bn1_ma_reg4];
			A_B1R2[bn1_ma_reg4] = data_tmp_1;
            //-----------------------------------------------------
		    data_tmp_1 = A_B1R5[bn1_ma_reg3];
		    data_tmp_2 = A_B1R6[bn1_ma_reg3];
		    data_tmp_3 = A_B1R7[bn1_ma_reg3];
			A_B1R5[bn1_ma_reg3] = A_B0R4[bn0_ma_reg3];
			A_B1R6[bn1_ma_reg3] = A_B0R4[bn0_ma_reg4];
			A_B1R7[bn1_ma_reg3] = A_B1R4[bn1_ma_reg4];
			A_B0R4[bn0_ma_reg3] = data_tmp_1;
			A_B0R4[bn0_ma_reg4] = data_tmp_2;
			A_B1R4[bn1_ma_reg4] = data_tmp_3;
			data_tmp_1 = A_B0R6[bn0_ma_reg3];
			data_tmp_2 = A_B0R7[bn0_ma_reg3];
			A_B0R6[bn0_ma_reg3] = A_B0R5[bn0_ma_reg4];
			A_B0R7[bn0_ma_reg3] = A_B1R5[bn1_ma_reg4];
			A_B0R5[bn0_ma_reg4] = data_tmp_1;
			A_B1R5[bn1_ma_reg4] = data_tmp_2;
			data_tmp_1 = A_B0R7[bn0_ma_reg4];
			A_B0R7[bn0_ma_reg4] = A_B1R6[bn1_ma_reg4];
			A_B1R6[bn1_ma_reg4] = data_tmp_1;
            //-----------------------------------------------------
		    data_tmp_1 = A_B1R9[bn1_ma_reg3];
		    data_tmp_2 = A_B1R10[bn1_ma_reg3];
		    data_tmp_3 = A_B1R11[bn1_ma_reg3];
			A_B1R9[bn1_ma_reg3]  = A_B0R8[bn0_ma_reg3];
			A_B1R10[bn1_ma_reg3] = A_B0R8[bn0_ma_reg4];
			A_B1R11[bn1_ma_reg3] = A_B1R8[bn1_ma_reg4];
			A_B0R8[bn0_ma_reg3] = data_tmp_1;
			A_B0R8[bn0_ma_reg4] = data_tmp_2;
			A_B1R8[bn1_ma_reg4] = data_tmp_3;
			data_tmp_1 = A_B0R10[bn0_ma_reg3];
			data_tmp_2 = A_B0R11[bn0_ma_reg3];
			A_B0R10[bn0_ma_reg3] = A_B0R9[bn0_ma_reg4];
			A_B0R11[bn0_ma_reg3] = A_B1R9[bn1_ma_reg4];
			A_B0R9[bn0_ma_reg4] = data_tmp_1;
			A_B1R9[bn1_ma_reg4] = data_tmp_2;
			data_tmp_1 = A_B0R11[bn0_ma_reg4];
			A_B0R11[bn0_ma_reg4] = A_B1R10[bn1_ma_reg4];
			A_B1R10[bn1_ma_reg4] = data_tmp_1;
            //-----------------------------------------------------
		    data_tmp_1 = A_B1R13[bn1_ma_reg3];
		    data_tmp_2 = A_B1R14[bn1_ma_reg3];
		    data_tmp_3 = A_B1R15[bn1_ma_reg3];
			A_B1R13[bn1_ma_reg3] = A_B0R12[bn0_ma_reg3];
			A_B1R14[bn1_ma_reg3] = A_B0R12[bn0_ma_reg4];
			A_B1R15[bn1_ma_reg3] = A_B1R12[bn1_ma_reg4];
			A_B0R12[bn0_ma_reg3] = data_tmp_1;
			A_B0R12[bn0_ma_reg4] = data_tmp_2;
			A_B1R12[bn1_ma_reg4] = data_tmp_3;
			data_tmp_1 = A_B0R14[bn0_ma_reg3];
			data_tmp_2 = A_B0R15[bn0_ma_reg3];
			A_B0R14[bn0_ma_reg3] = A_B0R13[bn0_ma_reg4];
			A_B0R15[bn0_ma_reg3] = A_B1R13[bn1_ma_reg4];
			A_B0R13[bn0_ma_reg4] = data_tmp_1;
			A_B1R13[bn1_ma_reg4] = data_tmp_2;
			data_tmp_1 = A_B0R15[bn0_ma_reg4];
			A_B0R15[bn0_ma_reg4] = A_B1R14[bn1_ma_reg4];
			A_B1R14[bn1_ma_reg4] = data_tmp_1;
            //*********************************************************	
			//bn 1 0 0 1
		    data_tmp_1 = A_B1R1[bn1_ma_reg5];
		    data_tmp_2 = A_B1R2[bn1_ma_reg5];
		    data_tmp_3 = A_B1R3[bn1_ma_reg5];
			A_B1R1[bn1_ma_reg5] = A_B0R0[bn0_ma_reg5];
			A_B1R2[bn1_ma_reg5] = A_B0R0[bn0_ma_reg6];
			A_B1R3[bn1_ma_reg5] = A_B1R0[bn1_ma_reg6];
			A_B0R0[bn0_ma_reg5] = data_tmp_1;
			A_B0R0[bn0_ma_reg6] = data_tmp_2;
			A_B1R0[bn1_ma_reg6] = data_tmp_3;
			data_tmp_1 = A_B0R2[bn0_ma_reg5];
			data_tmp_2 = A_B0R3[bn0_ma_reg5];
			A_B0R2[bn0_ma_reg5] = A_B0R1[bn0_ma_reg6];
			A_B0R3[bn0_ma_reg5] = A_B1R1[bn1_ma_reg6];
			A_B0R1[bn0_ma_reg6] = data_tmp_1;
			A_B1R1[bn1_ma_reg6] = data_tmp_2;
			data_tmp_1 = A_B0R3[bn0_ma_reg6];
			A_B0R3[bn0_ma_reg6] = A_B1R2[bn1_ma_reg6];
			A_B1R2[bn1_ma_reg6] = data_tmp_1;
            //-----------------------------------------------------
		    data_tmp_1 = A_B1R5[bn1_ma_reg5];
		    data_tmp_2 = A_B1R6[bn1_ma_reg5];
		    data_tmp_3 = A_B1R7[bn1_ma_reg5];
			A_B1R5[bn1_ma_reg5] = A_B0R4[bn0_ma_reg5];
			A_B1R6[bn1_ma_reg5] = A_B0R4[bn0_ma_reg6];
			A_B1R7[bn1_ma_reg5] = A_B1R4[bn1_ma_reg6];
			A_B0R4[bn0_ma_reg5] = data_tmp_1;
			A_B0R4[bn0_ma_reg6] = data_tmp_2;
			A_B1R4[bn1_ma_reg6] = data_tmp_3;
			data_tmp_1 = A_B0R6[bn0_ma_reg5];
			data_tmp_2 = A_B0R7[bn0_ma_reg5];
			A_B0R6[bn0_ma_reg5] = A_B0R5[bn0_ma_reg6];
			A_B0R7[bn0_ma_reg5] = A_B1R5[bn1_ma_reg6];
			A_B0R5[bn0_ma_reg6] = data_tmp_1;
			A_B1R5[bn1_ma_reg6] = data_tmp_2;
			data_tmp_1 = A_B0R7[bn0_ma_reg6];
			A_B0R7[bn0_ma_reg6] = A_B1R6[bn1_ma_reg6];
			A_B1R6[bn1_ma_reg6] = data_tmp_1;
            //-----------------------------------------------------
		    data_tmp_1 = A_B1R9[bn1_ma_reg5];
		    data_tmp_2 = A_B1R10[bn1_ma_reg5];
		    data_tmp_3 = A_B1R11[bn1_ma_reg5];
			A_B1R9[bn1_ma_reg5]  = A_B0R8[bn0_ma_reg5];
			A_B1R10[bn1_ma_reg5] = A_B0R8[bn0_ma_reg6];
			A_B1R11[bn1_ma_reg5] = A_B1R8[bn1_ma_reg6];
			A_B0R8[bn0_ma_reg5] = data_tmp_1;
			A_B0R8[bn0_ma_reg6] = data_tmp_2;
			A_B1R8[bn1_ma_reg6] = data_tmp_3;
			data_tmp_1 = A_B0R10[bn0_ma_reg5];
			data_tmp_2 = A_B0R11[bn0_ma_reg5];
			A_B0R10[bn0_ma_reg5] = A_B0R9[bn0_ma_reg6];
			A_B0R11[bn0_ma_reg5] = A_B1R9[bn1_ma_reg6];
			A_B0R9[bn0_ma_reg6] = data_tmp_1;
			A_B1R9[bn1_ma_reg6] = data_tmp_2;
			data_tmp_1 = A_B0R11[bn0_ma_reg6];
			A_B0R11[bn0_ma_reg6] = A_B1R10[bn1_ma_reg6];
			A_B1R10[bn1_ma_reg6] = data_tmp_1;
            //-----------------------------------------------------
		    data_tmp_1 = A_B1R13[bn1_ma_reg5];
		    data_tmp_2 = A_B1R14[bn1_ma_reg5];
		    data_tmp_3 = A_B1R15[bn1_ma_reg5];
			A_B1R13[bn1_ma_reg5] = A_B0R12[bn0_ma_reg5];
			A_B1R14[bn1_ma_reg5] = A_B0R12[bn0_ma_reg6];
			A_B1R15[bn1_ma_reg5] = A_B1R12[bn1_ma_reg6];
			A_B0R12[bn0_ma_reg5] = data_tmp_1;
			A_B0R12[bn0_ma_reg6] = data_tmp_2;
			A_B1R12[bn1_ma_reg6] = data_tmp_3;
			data_tmp_1 = A_B0R14[bn0_ma_reg5];
			data_tmp_2 = A_B0R15[bn0_ma_reg5];
			A_B0R14[bn0_ma_reg5] = A_B0R13[bn0_ma_reg6];
			A_B0R15[bn0_ma_reg5] = A_B1R13[bn1_ma_reg6];
			A_B0R13[bn0_ma_reg6] = data_tmp_1;
			A_B1R13[bn1_ma_reg6] = data_tmp_2;
			data_tmp_1 = A_B0R15[bn0_ma_reg6];
			A_B0R15[bn0_ma_reg6] = A_B1R14[bn1_ma_reg6];
			A_B1R14[bn1_ma_reg6] = data_tmp_1;
            //****************************************
            //bn: 0 1 1 0
		    data_tmp_1 = A_B0R1[bn0_ma_reg7];
		    data_tmp_2 = A_B0R2[bn0_ma_reg7];
		    data_tmp_3 = A_B0R3[bn0_ma_reg7];
			A_B0R1[bn0_ma_reg7] = A_B1R0[bn1_ma_reg7];
			A_B0R2[bn0_ma_reg7] = A_B1R0[bn1_ma_reg8];
			A_B0R3[bn0_ma_reg7] = A_B0R0[bn0_ma_reg8];
			A_B1R0[bn1_ma_reg7] = data_tmp_1;
			A_B1R0[bn1_ma_reg8] = data_tmp_2;
			A_B0R0[bn0_ma_reg8] = data_tmp_3;
			data_tmp_1 = A_B1R2[bn1_ma_reg7];
			data_tmp_2 = A_B1R3[bn1_ma_reg7];
			A_B1R2[bn1_ma_reg7] = A_B1R1[bn1_ma_reg8];
			A_B1R3[bn1_ma_reg7] = A_B0R1[bn0_ma_reg8];
			A_B1R1[bn1_ma_reg8] = data_tmp_1;
			A_B0R1[bn0_ma_reg8] = data_tmp_2;
			data_tmp_1 = A_B1R3[bn1_ma_reg8];
			A_B1R3[bn1_ma_reg8] = A_B0R2[bn0_ma_reg8];
			A_B0R2[bn0_ma_reg8] = data_tmp_1;
			//-----------------------------------------
			data_tmp_1 = A_B0R5[bn0_ma_reg7];
			data_tmp_2 = A_B0R6[bn0_ma_reg7];
			data_tmp_3 = A_B0R7[bn0_ma_reg7];
			A_B0R5[bn0_ma_reg7] = A_B1R4[bn1_ma_reg7];
			A_B0R6[bn0_ma_reg7] = A_B1R4[bn1_ma_reg8];
			A_B0R7[bn0_ma_reg7] = A_B0R4[bn0_ma_reg8];
			A_B1R4[bn1_ma_reg7] = data_tmp_1;
			A_B1R4[bn1_ma_reg8] = data_tmp_2;
			A_B0R4[bn0_ma_reg8] = data_tmp_3;
			data_tmp_1 = A_B1R6[bn1_ma_reg7];
			data_tmp_2 = A_B1R7[bn1_ma_reg7];
			A_B1R6[bn1_ma_reg7] = A_B1R5[bn1_ma_reg8];
			A_B1R7[bn1_ma_reg7] = A_B0R5[bn0_ma_reg8];
			A_B1R5[bn1_ma_reg8] = data_tmp_1;
			A_B0R5[bn0_ma_reg8] = data_tmp_2;
			data_tmp_1 = A_B1R7[bn1_ma_reg8];
			A_B1R7[bn1_ma_reg8] = A_B0R6[bn0_ma_reg8];
			A_B0R6[bn0_ma_reg8] = data_tmp_1;
			//----------------------------------------
			data_tmp_1 = A_B0R9[bn0_ma_reg7];
			data_tmp_2 = A_B0R10[bn0_ma_reg7];
			data_tmp_3 = A_B0R11[bn0_ma_reg7];
			A_B0R9[bn0_ma_reg7]  = A_B1R8[bn1_ma_reg7];
			A_B0R10[bn0_ma_reg7] = A_B1R8[bn1_ma_reg8];
			A_B0R11[bn0_ma_reg7] = A_B0R8[bn0_ma_reg8];
			A_B1R8[bn1_ma_reg7] = data_tmp_1;
			A_B1R8[bn1_ma_reg8] = data_tmp_2;
			A_B0R8[bn0_ma_reg8] = data_tmp_3;
			data_tmp_1 = A_B1R10[bn1_ma_reg7];
			data_tmp_2 = A_B1R11[bn1_ma_reg7];
			A_B1R10[bn1_ma_reg7] = A_B1R9[bn1_ma_reg8];
			A_B1R11[bn1_ma_reg7] = A_B0R9[bn0_ma_reg8];
			A_B1R9[bn1_ma_reg8] = data_tmp_1;
			A_B0R9[bn0_ma_reg8] = data_tmp_2;
			data_tmp_1 = A_B1R11[bn1_ma_reg8];
			A_B1R11[bn1_ma_reg8] = A_B0R10[bn0_ma_reg8];
			A_B0R10[bn0_ma_reg8] = data_tmp_1;			
			//----------------------------------------
			data_tmp_1 = A_B0R13[bn0_ma_reg7];
			data_tmp_2 = A_B0R14[bn0_ma_reg7];
			data_tmp_3 = A_B0R15[bn0_ma_reg7];
			A_B0R13[bn0_ma_reg7] = A_B1R12[bn1_ma_reg7];
			A_B0R14[bn0_ma_reg7] = A_B1R12[bn1_ma_reg8];
			A_B0R15[bn0_ma_reg7] = A_B0R12[bn0_ma_reg8];
			A_B1R12[bn1_ma_reg7] = data_tmp_1;
			A_B1R12[bn1_ma_reg8] = data_tmp_2;
			A_B0R12[bn0_ma_reg8] = data_tmp_3;
			data_tmp_1 = A_B1R14[bn1_ma_reg7];
			data_tmp_2 = A_B1R15[bn1_ma_reg7];
			A_B1R14[bn1_ma_reg7] = A_B1R13[bn1_ma_reg8];
			A_B1R15[bn1_ma_reg7] = A_B0R13[bn0_ma_reg8];
			A_B1R13[bn1_ma_reg8] = data_tmp_1;
			A_B0R13[bn0_ma_reg8] = data_tmp_2;
			data_tmp_1 = A_B1R15[bn1_ma_reg8];
			A_B1R15[bn1_ma_reg8] = A_B0R14[bn0_ma_reg8];
			A_B0R14[bn0_ma_reg8] = data_tmp_1;		
            //------------------------------------------------------			
		  }else{	 
			//bn 1 0 0 1
		    data_tmp_1 = A_B1R1[bn1_ma_reg1];
		    data_tmp_2 = A_B1R2[bn1_ma_reg1];
		    data_tmp_3 = A_B1R3[bn1_ma_reg1];
			A_B1R1[bn1_ma_reg1] = A_B0R0[bn0_ma_reg1];
			A_B1R2[bn1_ma_reg1] = A_B0R0[bn0_ma_reg2];
			A_B1R3[bn1_ma_reg1] = A_B1R0[bn1_ma_reg2];
			A_B0R0[bn0_ma_reg1] = data_tmp_1;
			A_B0R0[bn0_ma_reg2] = data_tmp_2;
			A_B1R0[bn1_ma_reg2] = data_tmp_3;
			data_tmp_1 = A_B0R2[bn0_ma_reg1];
			data_tmp_2 = A_B0R3[bn0_ma_reg1];
			A_B0R2[bn0_ma_reg1] = A_B0R1[bn0_ma_reg2];
			A_B0R3[bn0_ma_reg1] = A_B1R1[bn1_ma_reg2];
			A_B0R1[bn0_ma_reg2] = data_tmp_1;
			A_B1R1[bn1_ma_reg2] = data_tmp_2;
			data_tmp_1 = A_B0R3[bn0_ma_reg2];
			A_B0R3[bn0_ma_reg2] = A_B1R2[bn1_ma_reg2];
			A_B1R2[bn1_ma_reg2] = data_tmp_1;
            //-----------------------------------------------------
		    data_tmp_1 = A_B1R5[bn1_ma_reg1];
		    data_tmp_2 = A_B1R6[bn1_ma_reg1];
		    data_tmp_3 = A_B1R7[bn1_ma_reg1];
			A_B1R5[bn1_ma_reg1] = A_B0R4[bn0_ma_reg1];
			A_B1R6[bn1_ma_reg1] = A_B0R4[bn0_ma_reg2];
			A_B1R7[bn1_ma_reg1] = A_B1R4[bn1_ma_reg2];
			A_B0R4[bn0_ma_reg1] = data_tmp_1;
			A_B0R4[bn0_ma_reg2] = data_tmp_2;
			A_B1R4[bn1_ma_reg2] = data_tmp_3;
			data_tmp_1 = A_B0R6[bn0_ma_reg1];
			data_tmp_2 = A_B0R7[bn0_ma_reg1];
			A_B0R6[bn0_ma_reg1] = A_B0R5[bn0_ma_reg2];
			A_B0R7[bn0_ma_reg1] = A_B1R5[bn1_ma_reg2];
			A_B0R5[bn0_ma_reg2] = data_tmp_1;
			A_B1R5[bn1_ma_reg2] = data_tmp_2;
			data_tmp_1 = A_B0R7[bn0_ma_reg2];
			A_B0R7[bn0_ma_reg2] = A_B1R6[bn1_ma_reg2];
			A_B1R6[bn1_ma_reg2] = data_tmp_1;
            //-----------------------------------------------------
		    data_tmp_1 = A_B1R9[bn1_ma_reg1];
		    data_tmp_2 = A_B1R10[bn1_ma_reg1];
		    data_tmp_3 = A_B1R11[bn1_ma_reg1];
			A_B1R9[bn1_ma_reg1]  = A_B0R8[bn0_ma_reg1];
			A_B1R10[bn1_ma_reg1] = A_B0R8[bn0_ma_reg2];
			A_B1R11[bn1_ma_reg1] = A_B1R8[bn1_ma_reg2];
			A_B0R8[bn0_ma_reg1] = data_tmp_1;
			A_B0R8[bn0_ma_reg2] = data_tmp_2;
			A_B1R8[bn1_ma_reg2] = data_tmp_3;
			data_tmp_1 = A_B0R10[bn0_ma_reg1];
			data_tmp_2 = A_B0R11[bn0_ma_reg1];
			A_B0R10[bn0_ma_reg1] = A_B0R9[bn0_ma_reg2];
			A_B0R11[bn0_ma_reg1] = A_B1R9[bn1_ma_reg2];
			A_B0R9[bn0_ma_reg2] = data_tmp_1;
			A_B1R9[bn1_ma_reg2] = data_tmp_2;
			data_tmp_1 = A_B0R11[bn0_ma_reg2];
			A_B0R11[bn0_ma_reg2] = A_B1R10[bn1_ma_reg2];
			A_B1R10[bn1_ma_reg2] = data_tmp_1;
            //-----------------------------------------------------
		    data_tmp_1 = A_B1R13[bn1_ma_reg1];
		    data_tmp_2 = A_B1R14[bn1_ma_reg1];
		    data_tmp_3 = A_B1R15[bn1_ma_reg1];
			A_B1R13[bn1_ma_reg1] = A_B0R12[bn0_ma_reg1];
			A_B1R14[bn1_ma_reg1] = A_B0R12[bn0_ma_reg2];
			A_B1R15[bn1_ma_reg1] = A_B1R12[bn1_ma_reg2];
			A_B0R12[bn0_ma_reg1] = data_tmp_1;
			A_B0R12[bn0_ma_reg2] = data_tmp_2;
			A_B1R12[bn1_ma_reg2] = data_tmp_3;
			data_tmp_1 = A_B0R14[bn0_ma_reg1];
			data_tmp_2 = A_B0R15[bn0_ma_reg1];
			A_B0R14[bn0_ma_reg1] = A_B0R13[bn0_ma_reg2];
			A_B0R15[bn0_ma_reg1] = A_B1R13[bn1_ma_reg2];
			A_B0R13[bn0_ma_reg2] = data_tmp_1;
			A_B1R13[bn1_ma_reg2] = data_tmp_2;
			data_tmp_1 = A_B0R15[bn0_ma_reg2];
			A_B0R15[bn0_ma_reg2] = A_B1R14[bn1_ma_reg2];
			A_B1R14[bn1_ma_reg2] = data_tmp_1;
            //******************************************************
            //bn: 0 1 1 0
		    data_tmp_1 = A_B0R1[bn0_ma_reg3];
		    data_tmp_2 = A_B0R2[bn0_ma_reg3];
		    data_tmp_3 = A_B0R3[bn0_ma_reg3];
			A_B0R1[bn0_ma_reg3] = A_B1R0[bn1_ma_reg3];
			A_B0R2[bn0_ma_reg3] = A_B1R0[bn1_ma_reg4];
			A_B0R3[bn0_ma_reg3] = A_B0R0[bn0_ma_reg4];
			A_B1R0[bn1_ma_reg3] = data_tmp_1;
			A_B1R0[bn1_ma_reg4] = data_tmp_2;
			A_B0R0[bn0_ma_reg4] = data_tmp_3;
			data_tmp_1 = A_B1R2[bn1_ma_reg3];
			data_tmp_2 = A_B1R3[bn1_ma_reg3];
			A_B1R2[bn1_ma_reg3] = A_B1R1[bn1_ma_reg4];
			A_B1R3[bn1_ma_reg3] = A_B0R1[bn0_ma_reg4];
			A_B1R1[bn1_ma_reg4] = data_tmp_1;
			A_B0R1[bn0_ma_reg4] = data_tmp_2;
			data_tmp_1 = A_B1R3[bn1_ma_reg4];
			A_B1R3[bn1_ma_reg4] = A_B0R2[bn0_ma_reg4];
			A_B0R2[bn0_ma_reg4] = data_tmp_1;
			//-----------------------------------------
			data_tmp_1 = A_B0R5[bn0_ma_reg3];
			data_tmp_2 = A_B0R6[bn0_ma_reg3];
			data_tmp_3 = A_B0R7[bn0_ma_reg3];
			A_B0R5[bn0_ma_reg3] = A_B1R4[bn1_ma_reg3];
			A_B0R6[bn0_ma_reg3] = A_B1R4[bn1_ma_reg4];
			A_B0R7[bn0_ma_reg3] = A_B0R4[bn0_ma_reg4];
			A_B1R4[bn1_ma_reg3] = data_tmp_1;
			A_B1R4[bn1_ma_reg4] = data_tmp_2;
			A_B0R4[bn0_ma_reg4] = data_tmp_3;
			data_tmp_1 = A_B1R6[bn1_ma_reg3];
			data_tmp_2 = A_B1R7[bn1_ma_reg3];
			A_B1R6[bn1_ma_reg3] = A_B1R5[bn1_ma_reg4];
			A_B1R7[bn1_ma_reg3] = A_B0R5[bn0_ma_reg4];
			A_B1R5[bn1_ma_reg4] = data_tmp_1;
			A_B0R5[bn0_ma_reg4] = data_tmp_2;
			data_tmp_1 = A_B1R7[bn1_ma_reg4];
			A_B1R7[bn1_ma_reg4] = A_B0R6[bn0_ma_reg4];
			A_B0R6[bn0_ma_reg4] = data_tmp_1;
			//----------------------------------------
			data_tmp_1 = A_B0R9[bn0_ma_reg3];
			data_tmp_2 = A_B0R10[bn0_ma_reg3];
			data_tmp_3 = A_B0R11[bn0_ma_reg3];
			A_B0R9[bn0_ma_reg3]  = A_B1R8[bn1_ma_reg3];
			A_B0R10[bn0_ma_reg3] = A_B1R8[bn1_ma_reg4];
			A_B0R11[bn0_ma_reg3] = A_B0R8[bn0_ma_reg4];
			A_B1R8[bn1_ma_reg3] = data_tmp_1;
			A_B1R8[bn1_ma_reg4] = data_tmp_2;
			A_B0R8[bn0_ma_reg4] = data_tmp_3;
			data_tmp_1 = A_B1R10[bn1_ma_reg3];
			data_tmp_2 = A_B1R11[bn1_ma_reg3];
			A_B1R10[bn1_ma_reg3] = A_B1R9[bn1_ma_reg4];
			A_B1R11[bn1_ma_reg3] = A_B0R9[bn0_ma_reg4];
			A_B1R9[bn1_ma_reg4] = data_tmp_1;
			A_B0R9[bn0_ma_reg4] = data_tmp_2;
			data_tmp_1 = A_B1R11[bn1_ma_reg4];
			A_B1R11[bn1_ma_reg4] = A_B0R10[bn0_ma_reg4];
			A_B0R10[bn0_ma_reg4] = data_tmp_1;			
			//----------------------------------------
			data_tmp_1 = A_B0R13[bn0_ma_reg3];
			data_tmp_2 = A_B0R14[bn0_ma_reg3];
			data_tmp_3 = A_B0R15[bn0_ma_reg3];
			A_B0R13[bn0_ma_reg3] = A_B1R12[bn1_ma_reg3];
			A_B0R14[bn0_ma_reg3] = A_B1R12[bn1_ma_reg4];
			A_B0R15[bn0_ma_reg3] = A_B0R12[bn0_ma_reg4];
			A_B1R12[bn1_ma_reg3] = data_tmp_1;
			A_B1R12[bn1_ma_reg4] = data_tmp_2;
			A_B0R12[bn0_ma_reg4] = data_tmp_3;
			data_tmp_1 = A_B1R14[bn1_ma_reg3];
			data_tmp_2 = A_B1R15[bn1_ma_reg3];
			A_B1R14[bn1_ma_reg3] = A_B1R13[bn1_ma_reg4];
			A_B1R15[bn1_ma_reg3] = A_B0R13[bn0_ma_reg4];
			A_B1R13[bn1_ma_reg4] = data_tmp_1;
			A_B0R13[bn0_ma_reg4] = data_tmp_2;
			data_tmp_1 = A_B1R15[bn1_ma_reg4];
			A_B1R15[bn1_ma_reg4] = A_B0R14[bn0_ma_reg4];
			A_B0R14[bn0_ma_reg4] = data_tmp_1;
            //******************************************************
            //bn: 0 1 1 0
		    data_tmp_1 = A_B0R1[bn0_ma_reg5];
		    data_tmp_2 = A_B0R2[bn0_ma_reg5];
		    data_tmp_3 = A_B0R3[bn0_ma_reg5];
			A_B0R1[bn0_ma_reg5] = A_B1R0[bn1_ma_reg5];
			A_B0R2[bn0_ma_reg5] = A_B1R0[bn1_ma_reg6];
			A_B0R3[bn0_ma_reg5] = A_B0R0[bn0_ma_reg6];
			A_B1R0[bn1_ma_reg5] = data_tmp_1;
			A_B1R0[bn1_ma_reg6] = data_tmp_2;
			A_B0R0[bn0_ma_reg6] = data_tmp_3;
			data_tmp_1 = A_B1R2[bn1_ma_reg5];
			data_tmp_2 = A_B1R3[bn1_ma_reg5];
			A_B1R2[bn1_ma_reg5] = A_B1R1[bn1_ma_reg6];
			A_B1R3[bn1_ma_reg5] = A_B0R1[bn0_ma_reg6];
			A_B1R1[bn1_ma_reg6] = data_tmp_1;
			A_B0R1[bn0_ma_reg6] = data_tmp_2;
			data_tmp_1 = A_B1R3[bn1_ma_reg6];
			A_B1R3[bn1_ma_reg6] = A_B0R2[bn0_ma_reg6];
			A_B0R2[bn0_ma_reg6] = data_tmp_1;
			//-----------------------------------------
			data_tmp_1 = A_B0R5[bn0_ma_reg5];
			data_tmp_2 = A_B0R6[bn0_ma_reg5];
			data_tmp_3 = A_B0R7[bn0_ma_reg5];
			A_B0R5[bn0_ma_reg5] = A_B1R4[bn1_ma_reg5];
			A_B0R6[bn0_ma_reg5] = A_B1R4[bn1_ma_reg6];
			A_B0R7[bn0_ma_reg5] = A_B0R4[bn0_ma_reg6];
			A_B1R4[bn1_ma_reg5] = data_tmp_1;
			A_B1R4[bn1_ma_reg6] = data_tmp_2;
			A_B0R4[bn0_ma_reg6] = data_tmp_3;
			data_tmp_1 = A_B1R6[bn1_ma_reg5];
			data_tmp_2 = A_B1R7[bn1_ma_reg5];
			A_B1R6[bn1_ma_reg5] = A_B1R5[bn1_ma_reg6];
			A_B1R7[bn1_ma_reg5] = A_B0R5[bn0_ma_reg6];
			A_B1R5[bn1_ma_reg6] = data_tmp_1;
			A_B0R5[bn0_ma_reg6] = data_tmp_2;
			data_tmp_1 = A_B1R7[bn1_ma_reg6];
			A_B1R7[bn1_ma_reg6] = A_B0R6[bn0_ma_reg6];
			A_B0R6[bn0_ma_reg6] = data_tmp_1;
			//----------------------------------------
			data_tmp_1 = A_B0R9[bn0_ma_reg5];
			data_tmp_2 = A_B0R10[bn0_ma_reg5];
			data_tmp_3 = A_B0R11[bn0_ma_reg5];
			A_B0R9[bn0_ma_reg5]  = A_B1R8[bn1_ma_reg5];
			A_B0R10[bn0_ma_reg5] = A_B1R8[bn1_ma_reg6];
			A_B0R11[bn0_ma_reg5] = A_B0R8[bn0_ma_reg6];
			A_B1R8[bn1_ma_reg5] = data_tmp_1;
			A_B1R8[bn1_ma_reg6] = data_tmp_2;
			A_B0R8[bn0_ma_reg6] = data_tmp_3;
			data_tmp_1 = A_B1R10[bn1_ma_reg5];
			data_tmp_2 = A_B1R11[bn1_ma_reg5];
			A_B1R10[bn1_ma_reg5] = A_B1R9[bn1_ma_reg6];
			A_B1R11[bn1_ma_reg5] = A_B0R9[bn0_ma_reg6];
			A_B1R9[bn1_ma_reg6] = data_tmp_1;
			A_B0R9[bn0_ma_reg6] = data_tmp_2;
			data_tmp_1 = A_B1R11[bn1_ma_reg6];
			A_B1R11[bn1_ma_reg6] = A_B0R10[bn0_ma_reg6];
			A_B0R10[bn0_ma_reg6] = data_tmp_1;			
			//----------------------------------------
			data_tmp_1 = A_B0R13[bn0_ma_reg5];
			data_tmp_2 = A_B0R14[bn0_ma_reg5];
			data_tmp_3 = A_B0R15[bn0_ma_reg5];
			A_B0R13[bn0_ma_reg5] = A_B1R12[bn1_ma_reg5];
			A_B0R14[bn0_ma_reg5] = A_B1R12[bn1_ma_reg6];
			A_B0R15[bn0_ma_reg5] = A_B0R12[bn0_ma_reg6];
			A_B1R12[bn1_ma_reg5] = data_tmp_1;
			A_B1R12[bn1_ma_reg6] = data_tmp_2;
			A_B0R12[bn0_ma_reg6] = data_tmp_3;
			data_tmp_1 = A_B1R14[bn1_ma_reg5];
			data_tmp_2 = A_B1R15[bn1_ma_reg5];
			A_B1R14[bn1_ma_reg5] = A_B1R13[bn1_ma_reg6];
			A_B1R15[bn1_ma_reg5] = A_B0R13[bn0_ma_reg6];
			A_B1R13[bn1_ma_reg6] = data_tmp_1;
			A_B0R13[bn0_ma_reg6] = data_tmp_2;
			data_tmp_1 = A_B1R15[bn1_ma_reg6];
			A_B1R15[bn1_ma_reg6] = A_B0R14[bn0_ma_reg6];
			A_B0R14[bn0_ma_reg6] = data_tmp_1;			
			//***************************************************
			//bn 1 0 0 1
		    data_tmp_1 = A_B1R1[bn1_ma_reg7];
		    data_tmp_2 = A_B1R2[bn1_ma_reg7];
		    data_tmp_3 = A_B1R3[bn1_ma_reg7];
			A_B1R1[bn1_ma_reg7] = A_B0R0[bn0_ma_reg7];
			A_B1R2[bn1_ma_reg7] = A_B0R0[bn0_ma_reg8];
			A_B1R3[bn1_ma_reg7] = A_B1R0[bn1_ma_reg8];
			A_B0R0[bn0_ma_reg7] = data_tmp_1;
			A_B0R0[bn0_ma_reg8] = data_tmp_2;
			A_B1R0[bn1_ma_reg8] = data_tmp_3;
			data_tmp_1 = A_B0R2[bn0_ma_reg7];
			data_tmp_2 = A_B0R3[bn0_ma_reg7];
			A_B0R2[bn0_ma_reg7] = A_B0R1[bn0_ma_reg8];
			A_B0R3[bn0_ma_reg7] = A_B1R1[bn1_ma_reg8];
			A_B0R1[bn0_ma_reg8] = data_tmp_1;
			A_B1R1[bn1_ma_reg8] = data_tmp_2;
			data_tmp_1 = A_B0R3[bn0_ma_reg8];
			A_B0R3[bn0_ma_reg8] = A_B1R2[bn1_ma_reg8];
			A_B1R2[bn1_ma_reg8] = data_tmp_1;
            //-----------------------------------------------------
		    data_tmp_1 = A_B1R5[bn1_ma_reg7];
		    data_tmp_2 = A_B1R6[bn1_ma_reg7];
		    data_tmp_3 = A_B1R7[bn1_ma_reg7];
			A_B1R5[bn1_ma_reg7] = A_B0R4[bn0_ma_reg7];
			A_B1R6[bn1_ma_reg7] = A_B0R4[bn0_ma_reg8];
			A_B1R7[bn1_ma_reg7] = A_B1R4[bn1_ma_reg8];
			A_B0R4[bn0_ma_reg7] = data_tmp_1;
			A_B0R4[bn0_ma_reg8] = data_tmp_2;
			A_B1R4[bn1_ma_reg8] = data_tmp_3;
			data_tmp_1 = A_B0R6[bn0_ma_reg7];
			data_tmp_2 = A_B0R7[bn0_ma_reg7];
			A_B0R6[bn0_ma_reg7] = A_B0R5[bn0_ma_reg8];
			A_B0R7[bn0_ma_reg7] = A_B1R5[bn1_ma_reg8];
			A_B0R5[bn0_ma_reg8] = data_tmp_1;
			A_B1R5[bn1_ma_reg8] = data_tmp_2;
			data_tmp_1 = A_B0R7[bn0_ma_reg8];
			A_B0R7[bn0_ma_reg8] = A_B1R6[bn1_ma_reg8];
			A_B1R6[bn1_ma_reg8] = data_tmp_1;
            //-----------------------------------------------------
		    data_tmp_1 = A_B1R9[bn1_ma_reg7];
		    data_tmp_2 = A_B1R10[bn1_ma_reg7];
		    data_tmp_3 = A_B1R11[bn1_ma_reg7];
			A_B1R9[bn1_ma_reg7]  = A_B0R8[bn0_ma_reg7];
			A_B1R10[bn1_ma_reg7] = A_B0R8[bn0_ma_reg8];
			A_B1R11[bn1_ma_reg7] = A_B1R8[bn1_ma_reg8];
			A_B0R8[bn0_ma_reg7] = data_tmp_1;
			A_B0R8[bn0_ma_reg8] = data_tmp_2;
			A_B1R8[bn1_ma_reg8] = data_tmp_3;
			data_tmp_1 = A_B0R10[bn0_ma_reg7];
			data_tmp_2 = A_B0R11[bn0_ma_reg7];
			A_B0R10[bn0_ma_reg7] = A_B0R9[bn0_ma_reg8];
			A_B0R11[bn0_ma_reg7] = A_B1R9[bn1_ma_reg8];
			A_B0R9[bn0_ma_reg8] = data_tmp_1;
			A_B1R9[bn1_ma_reg8] = data_tmp_2;
			data_tmp_1 = A_B0R11[bn0_ma_reg8];
			A_B0R11[bn0_ma_reg8] = A_B1R10[bn1_ma_reg8];
			A_B1R10[bn1_ma_reg8] = data_tmp_1;
            //-----------------------------------------------------
		    data_tmp_1 = A_B1R13[bn1_ma_reg7];
		    data_tmp_2 = A_B1R14[bn1_ma_reg7];
		    data_tmp_3 = A_B1R15[bn1_ma_reg7];
			A_B1R13[bn1_ma_reg7] = A_B0R12[bn0_ma_reg7];
			A_B1R14[bn1_ma_reg7] = A_B0R12[bn0_ma_reg8];
			A_B1R15[bn1_ma_reg7] = A_B1R12[bn1_ma_reg8];
			A_B0R12[bn0_ma_reg7] = data_tmp_1;
			A_B0R12[bn0_ma_reg8] = data_tmp_2;
			A_B1R12[bn1_ma_reg8] = data_tmp_3;
			data_tmp_1 = A_B0R14[bn0_ma_reg7];
			data_tmp_2 = A_B0R15[bn0_ma_reg7];
			A_B0R14[bn0_ma_reg7] = A_B0R13[bn0_ma_reg8];
			A_B0R15[bn0_ma_reg7] = A_B1R13[bn1_ma_reg8];
			A_B0R13[bn0_ma_reg8] = data_tmp_1;
			A_B1R13[bn1_ma_reg8] = data_tmp_2;
			data_tmp_1 = A_B0R15[bn0_ma_reg8];
			A_B0R15[bn0_ma_reg8] = A_B1R14[bn1_ma_reg8];
			A_B1R14[bn1_ma_reg8] = data_tmp_1;			
		  }	
		 } 
		}
	}
	
	std::cout << "radix-16 IFFT computing over!!\n";
	
	DATARECORD <<" ----------------------------------------------- \n";
	DATARECORD <<" Radix-4 IFFT computing start!!! \n";
	
	for(int i = 0;i < group; i++){
        for(int j=0; j < radix; j++){
        	gray_i  = Gray(i,group);
        	BC_tmp  = j * group + gray_i;
			BC_RR = BC_tmp;
		    if(display ==1 )DATARECORD <<"---------------------------------------------------\n";
			if(display ==1 )DATARECORD <<" BC_tmp: "<< BC_tmp <<"\n";
			if(display ==1 )DATARECORD <<" ----------------------------------------\n";
			BC_IFFT_Reorder_R16_R4(BC_RR,BC_Reorder);
			if(display ==1 )DATARECORD <<" BC_Reorder:  "<< BC_Reorder <<"\n";
        	AGU_R16(BC_Reorder,bn_tmp,ma_tmp);
			if(display ==1 )DATARECORD <<" bn_tmp:  "<< bn_tmp <<"\n";
			if(display ==1 )DATARECORD <<" ma_tmp:  "<< ma_tmp <<"\n";			
        	if(bn_tmp == 0){
                if(j < 2)bn0_bc_tmp = BC_tmp;
				if(display == 1)DATARECORD <<" Before Radix-4 Butterfly unit computing!!! \n";
				if(display == 1)DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";										
                Radix4_BU_INTT(A_B0R0[ma_tmp],A_B0R1[ma_tmp],A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
        		Radix4_BU_INTT(A_B0R4[ma_tmp],A_B0R5[ma_tmp],A_B0R6[ma_tmp],A_B0R7[ma_tmp]);
        		Radix4_BU_INTT(A_B0R8[ma_tmp],A_B0R9[ma_tmp],A_B0R10[ma_tmp],A_B0R11[ma_tmp]);
        		Radix4_BU_INTT(A_B0R12[ma_tmp],A_B0R13[ma_tmp],A_B0R14[ma_tmp],A_B0R15[ma_tmp]);
				if(display == 1)DATARECORD <<" After Radix-4 Butterfly unit computing!!! \n";
				if(display == 1)DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";					
				if(j <  2)bn0_ma_reg1 = ma_tmp;
				if((j >= 2)  && (j < 4))bn0_ma_reg2 = ma_tmp;
				if((j >= 4)  && (j < 6))bn0_ma_reg3 = ma_tmp;
				if((j >= 6)  && (j < 8))bn0_ma_reg4 = ma_tmp;
				if((j >= 8)  && (j < 10))bn0_ma_reg5 = ma_tmp;
				if((j >= 10) && (j < 12))bn0_ma_reg6 = ma_tmp;
				if((j >= 12) && (j < 14))bn0_ma_reg7 = ma_tmp;
				if((j >= 14) && (j < 16))bn0_ma_reg8 = ma_tmp;				
        	}else {	
                if(j < 2)bn1_bc_tmp = BC_tmp;
				if(display == 1)DATARECORD <<" Before Radix-4 Butterfly unit computing!!! \n";
				if(display == 1)DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";					
        	    Radix4_BU_INTT(A_B1R0[ma_tmp],A_B1R1[ma_tmp],A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
        	    Radix4_BU_INTT(A_B1R4[ma_tmp],A_B1R5[ma_tmp],A_B1R6[ma_tmp],A_B1R7[ma_tmp]);
        	    Radix4_BU_INTT(A_B1R8[ma_tmp],A_B1R9[ma_tmp],A_B1R10[ma_tmp],A_B1R11[ma_tmp]);
        	    Radix4_BU_INTT(A_B1R12[ma_tmp],A_B1R13[ma_tmp],A_B1R14[ma_tmp],A_B1R15[ma_tmp]);
				if(display == 1)DATARECORD <<" After Radix-4 Butterfly unit computing!!! \n";
				if(display == 1)DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";
				if(display == 1)DATARECORD <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";					
                if(j <  2)bn1_ma_reg1 = ma_tmp;					
                if((j >= 2)  && (j < 4))bn1_ma_reg2 = ma_tmp;
                if((j >= 4)  && (j < 6))bn1_ma_reg3 = ma_tmp;
                if((j >= 6)  && (j < 8))bn1_ma_reg4 = ma_tmp;
                if((j >= 8)  && (j < 10))bn1_ma_reg5 = ma_tmp;
                if((j >= 10) && (j < 12))bn1_ma_reg6 = ma_tmp;
                if((j >= 12) && (j < 14))bn1_ma_reg7 = ma_tmp;
                if((j >= 14) && (j < 16))bn1_ma_reg8 = ma_tmp;				
        	}			
        }
		//data relocation for output 
		//0~ 15 ,16~31
		
		if(bn1_bc_tmp > bn0_bc_tmp){
			 data_tmp_1  = A_B0R1[bn0_ma_reg1];
			 data_tmp_2  = A_B0R2[bn0_ma_reg1];
			 data_tmp_3  = A_B0R3[bn0_ma_reg1];
			 data_tmp_4  = A_B0R4[bn0_ma_reg1];
			 data_tmp_5  = A_B0R5[bn0_ma_reg1];
			 data_tmp_6  = A_B0R6[bn0_ma_reg1];
			 data_tmp_7  = A_B0R7[bn0_ma_reg1];
			 data_tmp_8  = A_B0R8[bn0_ma_reg1];
			 data_tmp_9  = A_B0R9[bn0_ma_reg1];
			 data_tmp_10 = A_B0R10[bn0_ma_reg1];
			 data_tmp_11 = A_B0R11[bn0_ma_reg1];
			 data_tmp_12 = A_B0R12[bn0_ma_reg1];
			 data_tmp_13 = A_B0R13[bn0_ma_reg1];
			 data_tmp_14 = A_B0R14[bn0_ma_reg1];
			 data_tmp_15 = A_B0R15[bn0_ma_reg1];
			 A_B0R1[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg1];
			 A_B0R2[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg2];
			 A_B0R3[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg2];
			 A_B0R4[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg3];
			 A_B0R5[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg3];
			 A_B0R6[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg1] = A_B0R0[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg1] = A_B1R0[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg1] = A_B0R0[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg1] = A_B1R0[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg1] = A_B1R0[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg1] = A_B0R0[bn0_ma_reg8];
			 A_B1R0[bn1_ma_reg1]  = data_tmp_1; 
			 A_B1R0[bn1_ma_reg2]  = data_tmp_2; 
			 A_B0R0[bn0_ma_reg2]  = data_tmp_3; 
			 A_B1R0[bn1_ma_reg3]  = data_tmp_4; 
			 A_B0R0[bn0_ma_reg3]  = data_tmp_5; 
			 A_B0R0[bn0_ma_reg4]  = data_tmp_6; 
			 A_B1R0[bn1_ma_reg4]  = data_tmp_7; 
			 A_B1R0[bn1_ma_reg5]  = data_tmp_8; 
			 A_B0R0[bn0_ma_reg5]  = data_tmp_9; 
			 A_B0R0[bn0_ma_reg6]  = data_tmp_10;
			 A_B1R0[bn1_ma_reg6]  = data_tmp_11;
			 A_B0R0[bn0_ma_reg7]  = data_tmp_12;
			 A_B1R0[bn1_ma_reg7]  = data_tmp_13;
			 A_B1R0[bn1_ma_reg8]  = data_tmp_14;
			 A_B0R0[bn0_ma_reg8]  = data_tmp_15;
			 /*********************************************/
			 data_tmp_1  = A_B1R2[bn1_ma_reg1];
			 data_tmp_2  = A_B1R3[bn1_ma_reg1];
			 data_tmp_3  = A_B1R4[bn1_ma_reg1];
			 data_tmp_4  = A_B1R5[bn1_ma_reg1];
			 data_tmp_5  = A_B1R6[bn1_ma_reg1];
			 data_tmp_6  = A_B1R7[bn1_ma_reg1];
			 data_tmp_7  = A_B1R8[bn1_ma_reg1];
			 data_tmp_8  = A_B1R9[bn1_ma_reg1];
			 data_tmp_9  = A_B1R10[bn1_ma_reg1];
			 data_tmp_10 = A_B1R11[bn1_ma_reg1];
			 data_tmp_11 = A_B1R12[bn1_ma_reg1];
			 data_tmp_12 = A_B1R13[bn1_ma_reg1];
			 data_tmp_13 = A_B1R14[bn1_ma_reg1];
			 data_tmp_14 = A_B1R15[bn1_ma_reg1];
             A_B1R2[bn1_ma_reg1]  =	 A_B1R1[bn1_ma_reg2];
             A_B1R3[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg2];
             A_B1R4[bn1_ma_reg1]  =  A_B1R1[bn1_ma_reg3];
             A_B1R5[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg3];
             A_B1R6[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg4];
             A_B1R7[bn1_ma_reg1]  =  A_B1R1[bn1_ma_reg4];
             A_B1R8[bn1_ma_reg1]  =  A_B1R1[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg1] =  A_B0R1[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg1] =  A_B1R1[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg1] =  A_B0R1[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg1] =  A_B1R1[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg1] =  A_B1R1[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg1] =  A_B0R1[bn0_ma_reg8];
			 A_B1R1[bn1_ma_reg2]  =  data_tmp_1; 
			 A_B0R1[bn0_ma_reg2]  =  data_tmp_2; 
			 A_B1R1[bn1_ma_reg3]  =  data_tmp_3; 
			 A_B0R1[bn0_ma_reg3]  =  data_tmp_4; 
			 A_B0R1[bn0_ma_reg4]  =  data_tmp_5; 
			 A_B1R1[bn1_ma_reg4]  =  data_tmp_6; 
			 A_B1R1[bn1_ma_reg5]  =  data_tmp_7; 
			 A_B0R1[bn0_ma_reg5]  =  data_tmp_8; 
			 A_B0R1[bn0_ma_reg6]  =  data_tmp_9; 
			 A_B1R1[bn1_ma_reg6]  =  data_tmp_10;
			 A_B0R1[bn0_ma_reg7]  =  data_tmp_11;
			 A_B1R1[bn1_ma_reg7]  =  data_tmp_12;
			 A_B1R1[bn1_ma_reg8]  =  data_tmp_13;
			 A_B0R1[bn0_ma_reg8]  =  data_tmp_14;
			/************************************************************/ 
			 data_tmp_1  =   A_B1R3[bn1_ma_reg2];
			 data_tmp_2  =   A_B1R4[bn1_ma_reg2];
			 data_tmp_3  =   A_B1R5[bn1_ma_reg2];
			 data_tmp_4  =   A_B1R6[bn1_ma_reg2];
			 data_tmp_5  =   A_B1R7[bn1_ma_reg2];
			 data_tmp_6  =   A_B1R8[bn1_ma_reg2];
			 data_tmp_7  =   A_B1R9[bn1_ma_reg2];
			 data_tmp_8  =   A_B1R10[bn1_ma_reg2];
			 data_tmp_9  =   A_B1R11[bn1_ma_reg2];
			 data_tmp_10 =   A_B1R12[bn1_ma_reg2];
			 data_tmp_11 =   A_B1R13[bn1_ma_reg2];
			 data_tmp_12 =   A_B1R14[bn1_ma_reg2];
			 data_tmp_13 =   A_B1R15[bn1_ma_reg2];
			 A_B1R3[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg2];
			 A_B1R4[bn1_ma_reg2]  =  A_B1R2[bn1_ma_reg3];
			 A_B1R5[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg3];
			 A_B1R6[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg4];
			 A_B1R7[bn1_ma_reg2]  =  A_B1R2[bn1_ma_reg4];
			 A_B1R8[bn1_ma_reg2]  =  A_B1R2[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg2] =  A_B0R2[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg2] =  A_B1R2[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg2] =  A_B0R2[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg2] =  A_B1R2[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg2] =  A_B1R2[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg2] =  A_B0R2[bn0_ma_reg8];
			 A_B0R2[bn0_ma_reg2]  =  data_tmp_1; 
			 A_B1R2[bn1_ma_reg3]  =  data_tmp_2; 
			 A_B0R2[bn0_ma_reg3]  =  data_tmp_3; 
			 A_B0R2[bn0_ma_reg4]  =  data_tmp_4; 
			 A_B1R2[bn1_ma_reg4]  =  data_tmp_5; 
			 A_B1R2[bn1_ma_reg5]  =  data_tmp_6; 
			 A_B0R2[bn0_ma_reg5]  =  data_tmp_7; 
			 A_B0R2[bn0_ma_reg6]  =  data_tmp_8; 
			 A_B1R2[bn1_ma_reg6]  =  data_tmp_9; 
			 A_B0R2[bn0_ma_reg7]  =  data_tmp_10;
			 A_B1R2[bn1_ma_reg7]  =  data_tmp_11;
			 A_B1R2[bn1_ma_reg8]  =  data_tmp_12;
			 A_B0R2[bn0_ma_reg8]  =  data_tmp_13;
			 //******************************************************
			 data_tmp_1  = A_B0R4[bn0_ma_reg2];
			 data_tmp_2  = A_B0R5[bn0_ma_reg2];
			 data_tmp_3  = A_B0R6[bn0_ma_reg2];
			 data_tmp_4  = A_B0R7[bn0_ma_reg2];
			 data_tmp_5  = A_B0R8[bn0_ma_reg2];
			 data_tmp_6  = A_B0R9[bn0_ma_reg2];
			 data_tmp_7  = A_B0R10[bn0_ma_reg2];
			 data_tmp_8  = A_B0R11[bn0_ma_reg2];
			 data_tmp_9  = A_B0R12[bn0_ma_reg2];
			 data_tmp_10 = A_B0R13[bn0_ma_reg2];
			 data_tmp_11 = A_B0R14[bn0_ma_reg2];
			 data_tmp_12 = A_B0R15[bn0_ma_reg2];
			 A_B0R4[bn0_ma_reg2]  =  A_B1R3[bn1_ma_reg3];
			 A_B0R5[bn0_ma_reg2]  =  A_B0R3[bn0_ma_reg3];
			 A_B0R6[bn0_ma_reg2]  =  A_B0R3[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg2]  =  A_B1R3[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg2]  =  A_B1R3[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg2]  =  A_B0R3[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg2] =  A_B0R3[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg2] =  A_B1R3[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg2] =  A_B0R3[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg2] =  A_B1R3[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg2] =  A_B1R3[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg2] =  A_B0R3[bn0_ma_reg8];
			 A_B1R3[bn1_ma_reg3]  =  data_tmp_1;  
			 A_B0R3[bn0_ma_reg3]  =  data_tmp_2;  
			 A_B0R3[bn0_ma_reg4]  =  data_tmp_3;  
			 A_B1R3[bn1_ma_reg4]  =  data_tmp_4;  
			 A_B1R3[bn1_ma_reg5]  =  data_tmp_5;  
			 A_B0R3[bn0_ma_reg5]  =  data_tmp_6;  
			 A_B0R3[bn0_ma_reg6]  =  data_tmp_7;  
			 A_B1R3[bn1_ma_reg6]  =  data_tmp_8;  
			 A_B0R3[bn0_ma_reg7]  =  data_tmp_9;  
			 A_B1R3[bn1_ma_reg7]  =  data_tmp_10; 
			 A_B1R3[bn1_ma_reg8]  =  data_tmp_11; 
			 A_B0R3[bn0_ma_reg8]  =  data_tmp_12; 
			 //----------------------------------------------------------------------
             data_tmp_1  = A_B1R5[bn1_ma_reg3];
			 data_tmp_2  = A_B1R6[bn1_ma_reg3];
			 data_tmp_3  = A_B1R7[bn1_ma_reg3];
			 data_tmp_4  = A_B1R8[bn1_ma_reg3];
			 data_tmp_5  = A_B1R9[bn1_ma_reg3];
			 data_tmp_6  = A_B1R10[bn1_ma_reg3];
			 data_tmp_7  = A_B1R11[bn1_ma_reg3];
			 data_tmp_8  = A_B1R12[bn1_ma_reg3];
			 data_tmp_9  = A_B1R13[bn1_ma_reg3];
			 data_tmp_10 = A_B1R14[bn1_ma_reg3];
			 data_tmp_11 = A_B1R15[bn1_ma_reg3];
			 A_B1R5[bn1_ma_reg3]  =  A_B0R4[bn0_ma_reg3];
			 A_B1R6[bn1_ma_reg3]  =  A_B0R4[bn0_ma_reg4];
			 A_B1R7[bn1_ma_reg3]  =  A_B1R4[bn1_ma_reg4];
			 A_B1R8[bn1_ma_reg3]  =  A_B1R4[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg3]  =  A_B0R4[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg3] =  A_B0R4[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg3] =  A_B1R4[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg3] =  A_B0R4[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg3] =  A_B1R4[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg3] =  A_B1R4[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg3] =  A_B0R4[bn0_ma_reg8];
			 A_B0R4[bn0_ma_reg3]  =  data_tmp_1 ;
			 A_B0R4[bn0_ma_reg4]  =  data_tmp_2 ;
			 A_B1R4[bn1_ma_reg4]  =  data_tmp_3 ;
			 A_B1R4[bn1_ma_reg5]  =  data_tmp_4 ;
			 A_B0R4[bn0_ma_reg5]  =  data_tmp_5 ;
			 A_B0R4[bn0_ma_reg6]  =  data_tmp_6 ;
			 A_B1R4[bn1_ma_reg6]  =  data_tmp_7 ;
			 A_B0R4[bn0_ma_reg7]  =  data_tmp_8 ;
			 A_B1R4[bn1_ma_reg7]  =  data_tmp_9 ;
			 A_B1R4[bn1_ma_reg8]  =  data_tmp_10;
			 A_B0R4[bn0_ma_reg8]  =  data_tmp_11;
			 //************************************************************************
			 data_tmp_1  = A_B0R6[bn0_ma_reg3];
			 data_tmp_2  = A_B0R7[bn0_ma_reg3];
			 data_tmp_3  = A_B0R8[bn0_ma_reg3];
			 data_tmp_4  = A_B0R9[bn0_ma_reg3];
			 data_tmp_5  = A_B0R10[bn0_ma_reg3];
			 data_tmp_6  = A_B0R11[bn0_ma_reg3];
			 data_tmp_7  = A_B0R12[bn0_ma_reg3];
			 data_tmp_8  = A_B0R13[bn0_ma_reg3];
			 data_tmp_9  = A_B0R14[bn0_ma_reg3];
			 data_tmp_10 = A_B0R15[bn0_ma_reg3];
			 A_B0R6[bn0_ma_reg3]  =  A_B0R5[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg3]  =  A_B1R5[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg3]  =  A_B1R5[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg3]  =  A_B0R5[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg3] =  A_B0R5[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg3] =  A_B1R5[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg3] =  A_B0R5[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg3] =  A_B1R5[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg3] =  A_B1R5[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg3] =  A_B0R5[bn0_ma_reg8];
			 A_B0R5[bn0_ma_reg4]  =  data_tmp_1; 
			 A_B1R5[bn1_ma_reg4]  =  data_tmp_2; 
			 A_B1R5[bn1_ma_reg5]  =  data_tmp_3; 
			 A_B0R5[bn0_ma_reg5]  =  data_tmp_4; 
			 A_B0R5[bn0_ma_reg6]  =  data_tmp_5; 
			 A_B1R5[bn1_ma_reg6]  =  data_tmp_6; 
			 A_B0R5[bn0_ma_reg7]  =  data_tmp_7; 
			 A_B1R5[bn1_ma_reg7]  =  data_tmp_8; 
			 A_B1R5[bn1_ma_reg8]  =  data_tmp_9; 
			 A_B0R5[bn0_ma_reg8]  =  data_tmp_10;
			 //-----------------------------------------------------------------------
			 data_tmp_1  = A_B0R7[bn0_ma_reg4];
			 data_tmp_2  = A_B0R8[bn0_ma_reg4];
			 data_tmp_3  = A_B0R9[bn0_ma_reg4];
			 data_tmp_4  = A_B0R10[bn0_ma_reg4];
			 data_tmp_5  = A_B0R11[bn0_ma_reg4];
			 data_tmp_6  = A_B0R12[bn0_ma_reg4];
			 data_tmp_7  = A_B0R13[bn0_ma_reg4];
			 data_tmp_8  = A_B0R14[bn0_ma_reg4];
			 data_tmp_9  = A_B0R15[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg4]  =  A_B1R6[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg4]  =  A_B1R6[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg4]  =  A_B0R6[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg4] =  A_B0R6[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg4] =  A_B1R6[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg4] =  A_B0R6[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg4] =  A_B1R6[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg4] =  A_B1R6[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg4] =  A_B0R6[bn0_ma_reg8];
			 A_B1R6[bn1_ma_reg4]  =  data_tmp_1;
			 A_B1R6[bn1_ma_reg5]  =  data_tmp_2;
			 A_B0R6[bn0_ma_reg5]  =  data_tmp_3;
			 A_B0R6[bn0_ma_reg6]  =  data_tmp_4;
			 A_B1R6[bn1_ma_reg6]  =  data_tmp_5;
			 A_B0R6[bn0_ma_reg7]  =  data_tmp_6;
			 A_B1R6[bn1_ma_reg7]  =  data_tmp_7;
			 A_B1R6[bn1_ma_reg8]  =  data_tmp_8;
			 A_B0R6[bn0_ma_reg8]  =  data_tmp_9;
			 //----------------------------------------------------------------------
			 data_tmp_1  = A_B1R8[bn1_ma_reg4];
			 data_tmp_2  = A_B1R9[bn1_ma_reg4];
			 data_tmp_3  = A_B1R10[bn1_ma_reg4];
			 data_tmp_4  = A_B1R11[bn1_ma_reg4];
			 data_tmp_5  = A_B1R12[bn1_ma_reg4];
			 data_tmp_6  = A_B1R13[bn1_ma_reg4];
			 data_tmp_7  = A_B1R14[bn1_ma_reg4];
			 data_tmp_8  = A_B1R15[bn1_ma_reg4];
			 A_B1R8[bn1_ma_reg4]  =  A_B1R7[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg4]  =  A_B0R7[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg4] =  A_B0R7[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg4] =  A_B1R7[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg4] =  A_B0R7[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg4] =  A_B1R7[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg4] =  A_B1R7[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg4] =  A_B0R7[bn0_ma_reg8];
			 A_B1R7[bn1_ma_reg5]  =  data_tmp_1;
			 A_B0R7[bn0_ma_reg5]  =  data_tmp_2;
			 A_B0R7[bn0_ma_reg6]  =  data_tmp_3;
			 A_B1R7[bn1_ma_reg6]  =  data_tmp_4;
			 A_B0R7[bn0_ma_reg7]  =  data_tmp_5;
			 A_B1R7[bn1_ma_reg7]  =  data_tmp_6;
			 A_B1R7[bn1_ma_reg8]  =  data_tmp_7;
			 A_B0R7[bn0_ma_reg8]  =  data_tmp_8;
			 //----------------------------------------------------------------------
			 data_tmp_1 = A_B1R9[bn1_ma_reg5];
			 data_tmp_2 = A_B1R10[bn1_ma_reg5];
			 data_tmp_3 = A_B1R11[bn1_ma_reg5];
			 data_tmp_4 = A_B1R12[bn1_ma_reg5];
			 data_tmp_5 = A_B1R13[bn1_ma_reg5];
			 data_tmp_6 = A_B1R14[bn1_ma_reg5];
			 data_tmp_7 = A_B1R15[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg5]  = A_B0R8[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg5] = A_B0R8[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg5] = A_B1R8[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg5] = A_B0R8[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg5] = A_B1R8[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg5] = A_B1R8[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg5] = A_B0R8[bn0_ma_reg8];
			 A_B0R8[bn0_ma_reg5]  = data_tmp_1;
			 A_B0R8[bn0_ma_reg6]  = data_tmp_2;
			 A_B1R8[bn1_ma_reg6]  = data_tmp_3;
			 A_B0R8[bn0_ma_reg7]  = data_tmp_4;
			 A_B1R8[bn1_ma_reg7]  = data_tmp_5;
			 A_B1R8[bn1_ma_reg8]  = data_tmp_6;
			 A_B0R8[bn0_ma_reg8]  = data_tmp_7;
			 //---------------------------------------------------------------------
			 data_tmp_1  = A_B0R10[bn0_ma_reg5];
			 data_tmp_2  = A_B0R11[bn0_ma_reg5];
			 data_tmp_3  = A_B0R12[bn0_ma_reg5];
			 data_tmp_4  = A_B0R13[bn0_ma_reg5];
			 data_tmp_5  = A_B0R14[bn0_ma_reg5];
			 data_tmp_6  = A_B0R15[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg5] = A_B0R9[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg5] = A_B1R9[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg5] = A_B0R9[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg5] = A_B1R9[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg5] = A_B1R9[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg5] = A_B0R9[bn0_ma_reg8];
			 A_B0R9[bn0_ma_reg6]  = data_tmp_1;
			 A_B1R9[bn1_ma_reg6]  = data_tmp_2;
			 A_B0R9[bn0_ma_reg7]  = data_tmp_3;
			 A_B1R9[bn1_ma_reg7]  = data_tmp_4;
			 A_B1R9[bn1_ma_reg8]  = data_tmp_5;
			 A_B0R9[bn0_ma_reg8]  = data_tmp_6;
			 //--------------------------------------------------------------------
			 data_tmp_1  = A_B0R11[bn0_ma_reg6];
			 data_tmp_2  = A_B0R12[bn0_ma_reg6];
			 data_tmp_3  = A_B0R13[bn0_ma_reg6];
			 data_tmp_4  = A_B0R14[bn0_ma_reg6];
			 data_tmp_5  = A_B0R15[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg6] = A_B1R10[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg6] = A_B0R10[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg6] = A_B1R10[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg6] = A_B1R10[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg6] = A_B0R10[bn0_ma_reg8];
			 A_B1R10[bn1_ma_reg6] = data_tmp_1;
			 A_B0R10[bn0_ma_reg7] = data_tmp_2;
			 A_B1R10[bn1_ma_reg7] = data_tmp_3;
			 A_B1R10[bn1_ma_reg8] = data_tmp_4;
			 A_B0R10[bn0_ma_reg8] = data_tmp_5;
			 //--------------------------------------------------------------------
			 data_tmp_1  = A_B1R12[bn1_ma_reg6];
			 data_tmp_2  = A_B1R13[bn1_ma_reg6];
			 data_tmp_3  = A_B1R14[bn1_ma_reg6];
			 data_tmp_4  = A_B1R15[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg6] = A_B0R11[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg6] = A_B1R11[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg6] = A_B1R11[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg6] = A_B0R11[bn0_ma_reg8];
			 A_B0R11[bn0_ma_reg7] = data_tmp_1;
			 A_B1R11[bn1_ma_reg7] = data_tmp_2;
			 A_B1R11[bn1_ma_reg8] = data_tmp_3;
			 A_B0R11[bn0_ma_reg8] = data_tmp_4;
			 //--------------------------------------------------------------------
			 data_tmp_1 = A_B0R13[bn0_ma_reg7];
			 data_tmp_2 = A_B0R14[bn0_ma_reg7];
			 data_tmp_3 = A_B0R15[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg7] = A_B1R12[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg7] = A_B1R12[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg7] = A_B0R12[bn0_ma_reg8];
			 A_B1R12[bn1_ma_reg7] = data_tmp_1;
			 A_B1R12[bn1_ma_reg8] = data_tmp_2;
			 A_B0R12[bn0_ma_reg8] = data_tmp_3;
			 //--------------------------------------------------------------------
			 data_tmp_1 = A_B1R14[bn1_ma_reg7];
			 data_tmp_2 = A_B1R15[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg7] = A_B1R13[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg7] = A_B0R13[bn0_ma_reg8];
			 A_B1R13[bn1_ma_reg8] = data_tmp_1;
			 A_B0R13[bn0_ma_reg8] = data_tmp_2;
			 //--------------------------------------------------------------------
			 data_tmp_1 = A_B1R15[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg8] = A_B0R14[bn0_ma_reg8];
			 A_B0R14[bn0_ma_reg8] = data_tmp_1; 
		  }else{	 
			 data_tmp_1  = A_B1R1[bn1_ma_reg1];
			 data_tmp_2  = A_B1R2[bn1_ma_reg1];
			 data_tmp_3  = A_B1R3[bn1_ma_reg1];
			 data_tmp_4  = A_B1R4[bn1_ma_reg1];
			 data_tmp_5  = A_B1R5[bn1_ma_reg1];
			 data_tmp_6  = A_B1R6[bn1_ma_reg1];
			 data_tmp_7  = A_B1R7[bn1_ma_reg1];
			 data_tmp_8  = A_B1R8[bn1_ma_reg1];
			 data_tmp_9  = A_B1R9[bn1_ma_reg1];
			 data_tmp_10 = A_B1R10[bn1_ma_reg1];
			 data_tmp_11 = A_B1R11[bn1_ma_reg1];
			 data_tmp_12 = A_B1R12[bn1_ma_reg1];
			 data_tmp_13 = A_B1R13[bn1_ma_reg1];
			 data_tmp_14 = A_B1R14[bn1_ma_reg1];
			 data_tmp_15 = A_B1R15[bn1_ma_reg1];
             A_B1R1[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg1];
             A_B1R2[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg2];
             A_B1R3[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg2];
             A_B1R4[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg3];
             A_B1R5[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg3];
             A_B1R6[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg4];
             A_B1R7[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg4];
             A_B1R8[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg5];
             A_B1R9[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg5];
             A_B1R10[bn1_ma_reg1] =  A_B1R0[bn1_ma_reg6];
             A_B1R11[bn1_ma_reg1] =  A_B0R0[bn0_ma_reg6];
             A_B1R12[bn1_ma_reg1] =  A_B1R0[bn1_ma_reg7];
             A_B1R13[bn1_ma_reg1] =  A_B0R0[bn0_ma_reg7];
             A_B1R14[bn1_ma_reg1] =  A_B0R0[bn0_ma_reg8];
             A_B1R15[bn1_ma_reg1] =  A_B1R0[bn1_ma_reg8];
             A_B0R0[bn0_ma_reg1]  =  data_tmp_1; 
             A_B0R0[bn0_ma_reg2]  =  data_tmp_2; 
             A_B1R0[bn1_ma_reg2]  =  data_tmp_3; 
             A_B0R0[bn0_ma_reg3]  =  data_tmp_4; 
             A_B1R0[bn1_ma_reg3]  =  data_tmp_5; 
             A_B1R0[bn1_ma_reg4]  =  data_tmp_6; 
             A_B0R0[bn0_ma_reg4]  =  data_tmp_7; 
             A_B0R0[bn0_ma_reg5]  =  data_tmp_8; 
             A_B1R0[bn1_ma_reg5]  =  data_tmp_9; 
             A_B1R0[bn1_ma_reg6]  =  data_tmp_10;
             A_B0R0[bn0_ma_reg6]  =  data_tmp_11;
             A_B1R0[bn1_ma_reg7]  =  data_tmp_12;
             A_B0R0[bn0_ma_reg7]  =  data_tmp_13;
             A_B0R0[bn0_ma_reg8]  =  data_tmp_14;
             A_B1R0[bn1_ma_reg8]  =  data_tmp_15;
             //-----------------------------------------------------------
             data_tmp_1  = A_B0R2[bn0_ma_reg1];
			 data_tmp_2  = A_B0R3[bn0_ma_reg1];
			 data_tmp_3  = A_B0R4[bn0_ma_reg1];
			 data_tmp_4  = A_B0R5[bn0_ma_reg1];
			 data_tmp_5  = A_B0R6[bn0_ma_reg1];
			 data_tmp_6  = A_B0R7[bn0_ma_reg1];
			 data_tmp_7  = A_B0R8[bn0_ma_reg1];
			 data_tmp_8  = A_B0R9[bn0_ma_reg1];
			 data_tmp_9  = A_B0R10[bn0_ma_reg1];
			 data_tmp_10 = A_B0R11[bn0_ma_reg1];
			 data_tmp_11 = A_B0R12[bn0_ma_reg1];
			 data_tmp_12 = A_B0R13[bn0_ma_reg1];
			 data_tmp_13 = A_B0R14[bn0_ma_reg1];
			 data_tmp_14 = A_B0R15[bn0_ma_reg1];
			 A_B0R2[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg2];
			 A_B0R3[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg2];
			 A_B0R4[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg3];
			 A_B0R5[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg3];
			 A_B0R6[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg4];
			 A_B0R7[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg1] =  A_B1R1[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg1] =  A_B0R1[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg1] =  A_B1R1[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg1] =  A_B0R1[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg1] =  A_B0R1[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg1] =  A_B1R1[bn1_ma_reg8];
             A_B0R1[bn0_ma_reg2]  =  data_tmp_1; 
             A_B1R1[bn1_ma_reg2]  =  data_tmp_2; 
             A_B0R1[bn0_ma_reg3]  =  data_tmp_3; 
             A_B1R1[bn1_ma_reg3]  =  data_tmp_4; 
             A_B1R1[bn1_ma_reg4]  =  data_tmp_5; 
             A_B0R1[bn0_ma_reg4]  =  data_tmp_6; 
             A_B0R1[bn0_ma_reg5]  =  data_tmp_7; 
             A_B1R1[bn1_ma_reg5]  =  data_tmp_8; 
             A_B1R1[bn1_ma_reg6]  =  data_tmp_9; 
             A_B0R1[bn0_ma_reg6]  =  data_tmp_10;
             A_B1R1[bn1_ma_reg7]  =  data_tmp_11;
             A_B0R1[bn0_ma_reg7]  =  data_tmp_12;
             A_B0R1[bn0_ma_reg8]  =  data_tmp_13;
             A_B1R1[bn1_ma_reg8]  =  data_tmp_14;
             //------------------------------------------------------------
             data_tmp_1  =  A_B0R3[bn0_ma_reg2];
			 data_tmp_2  =  A_B0R4[bn0_ma_reg2];
			 data_tmp_3  =  A_B0R5[bn0_ma_reg2];
			 data_tmp_4  =  A_B0R6[bn0_ma_reg2];
			 data_tmp_5  =  A_B0R7[bn0_ma_reg2];
			 data_tmp_6  =  A_B0R8[bn0_ma_reg2];
			 data_tmp_7  =  A_B0R9[bn0_ma_reg2];
			 data_tmp_8  =  A_B0R10[bn0_ma_reg2];
			 data_tmp_9  =  A_B0R11[bn0_ma_reg2];
			 data_tmp_10 =  A_B0R12[bn0_ma_reg2];
			 data_tmp_11 =  A_B0R13[bn0_ma_reg2];
			 data_tmp_12 =  A_B0R14[bn0_ma_reg2];
			 data_tmp_13 =  A_B0R15[bn0_ma_reg2];
			 A_B0R3[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg2];
			 A_B0R4[bn0_ma_reg2]  =  A_B0R2[bn0_ma_reg3];
			 A_B0R5[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg3];
			 A_B0R6[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg4];
			 A_B0R7[bn0_ma_reg2]  =  A_B0R2[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg2]  =  A_B0R2[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg2] =  A_B1R2[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg2] =  A_B0R2[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg2] =  A_B1R2[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg2] =  A_B0R2[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg2] =  A_B0R2[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg2] =  A_B1R2[bn1_ma_reg8];
			 A_B1R2[bn1_ma_reg2]  =  data_tmp_1; 
			 A_B0R2[bn0_ma_reg3]  =  data_tmp_2; 
			 A_B1R2[bn1_ma_reg3]  =  data_tmp_3; 
			 A_B1R2[bn1_ma_reg4]  =  data_tmp_4; 
			 A_B0R2[bn0_ma_reg4]  =  data_tmp_5; 
			 A_B0R2[bn0_ma_reg5]  =  data_tmp_6; 
			 A_B1R2[bn1_ma_reg5]  =  data_tmp_7; 
			 A_B1R2[bn1_ma_reg6]  =  data_tmp_8; 
			 A_B0R2[bn0_ma_reg6]  =  data_tmp_9; 
			 A_B1R2[bn1_ma_reg7]  =  data_tmp_10;
			 A_B0R2[bn0_ma_reg7]  =  data_tmp_11;
			 A_B0R2[bn0_ma_reg8]  =  data_tmp_12;
			 A_B1R2[bn1_ma_reg8]  =  data_tmp_13;
			 //-----------------------------------------------------------                        
			 data_tmp_1  = A_B1R4[bn1_ma_reg2];
			 data_tmp_2  = A_B1R5[bn1_ma_reg2];
			 data_tmp_3  = A_B1R6[bn1_ma_reg2];
			 data_tmp_4  = A_B1R7[bn1_ma_reg2];
			 data_tmp_5  = A_B1R8[bn1_ma_reg2];
			 data_tmp_6  = A_B1R9[bn1_ma_reg2];
			 data_tmp_7  = A_B1R10[bn1_ma_reg2];
			 data_tmp_8  = A_B1R11[bn1_ma_reg2];
			 data_tmp_9  = A_B1R12[bn1_ma_reg2];
			 data_tmp_10 = A_B1R13[bn1_ma_reg2];
			 data_tmp_11 = A_B1R14[bn1_ma_reg2];
			 data_tmp_12 = A_B1R15[bn1_ma_reg2];
			 A_B1R4[bn1_ma_reg2]  = A_B0R3[bn0_ma_reg3];
			 A_B1R5[bn1_ma_reg2]  = A_B1R3[bn1_ma_reg3];
			 A_B1R6[bn1_ma_reg2]  = A_B1R3[bn1_ma_reg4];
			 A_B1R7[bn1_ma_reg2]  = A_B0R3[bn0_ma_reg4];
			 A_B1R8[bn1_ma_reg2]  = A_B0R3[bn0_ma_reg5];
			 A_B1R9[bn1_ma_reg2]  = A_B1R3[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg2] = A_B1R3[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg2] = A_B0R3[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg2] = A_B1R3[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg2] = A_B0R3[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg2] = A_B0R3[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg2] = A_B1R3[bn1_ma_reg8];
			 A_B0R3[bn0_ma_reg3]  = data_tmp_1; 
			 A_B1R3[bn1_ma_reg3]  = data_tmp_2; 
			 A_B1R3[bn1_ma_reg4]  = data_tmp_3; 
			 A_B0R3[bn0_ma_reg4]  = data_tmp_4; 
			 A_B0R3[bn0_ma_reg5]  = data_tmp_5; 
			 A_B1R3[bn1_ma_reg5]  = data_tmp_6; 
			 A_B1R3[bn1_ma_reg6]  = data_tmp_7; 
			 A_B0R3[bn0_ma_reg6]  = data_tmp_8; 
			 A_B1R3[bn1_ma_reg7]  = data_tmp_9; 
			 A_B0R3[bn0_ma_reg7]  = data_tmp_10;
			 A_B0R3[bn0_ma_reg8]  = data_tmp_11;
			 A_B1R3[bn1_ma_reg8]  = data_tmp_12;
			 //------------------------------------------------------------
			 data_tmp_1  =  A_B0R5[bn0_ma_reg3];
			 data_tmp_2  =  A_B0R6[bn0_ma_reg3];
			 data_tmp_3  =  A_B0R7[bn0_ma_reg3];
			 data_tmp_4  =  A_B0R8[bn0_ma_reg3];
			 data_tmp_5  =  A_B0R9[bn0_ma_reg3];
			 data_tmp_6  =  A_B0R10[bn0_ma_reg3];
			 data_tmp_7  =  A_B0R11[bn0_ma_reg3];
			 data_tmp_8  =  A_B0R12[bn0_ma_reg3];
			 data_tmp_9  =  A_B0R13[bn0_ma_reg3];
			 data_tmp_10 =  A_B0R14[bn0_ma_reg3];
			 data_tmp_11 =  A_B0R15[bn0_ma_reg3];
			 A_B0R5[bn0_ma_reg3]  =  A_B1R4[bn1_ma_reg3];
			 A_B0R6[bn0_ma_reg3]  =  A_B1R4[bn1_ma_reg4];
			 A_B0R7[bn0_ma_reg3]  =  A_B0R4[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg3]  =  A_B0R4[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg3]  =  A_B1R4[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg3] =  A_B1R4[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg3] =  A_B0R4[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg3] =  A_B1R4[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg3] =  A_B0R4[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg3] =  A_B0R4[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg3] =  A_B1R4[bn1_ma_reg8];
			 A_B1R4[bn1_ma_reg3]  =  data_tmp_1; 
			 A_B1R4[bn1_ma_reg4]  =  data_tmp_2; 
			 A_B0R4[bn0_ma_reg4]  =  data_tmp_3; 
			 A_B0R4[bn0_ma_reg5]  =  data_tmp_4; 
			 A_B1R4[bn1_ma_reg5]  =  data_tmp_5; 
			 A_B1R4[bn1_ma_reg6]  =  data_tmp_6; 
			 A_B0R4[bn0_ma_reg6]  =  data_tmp_7; 
			 A_B1R4[bn1_ma_reg7]  =  data_tmp_8; 
			 A_B0R4[bn0_ma_reg7]  =  data_tmp_9; 
			 A_B0R4[bn0_ma_reg8]  =  data_tmp_10;
			 A_B1R4[bn1_ma_reg8]  =  data_tmp_11;
			 //-------------------------------------------------------------
			 data_tmp_1  =  A_B1R6[bn1_ma_reg3];
			 data_tmp_2  =  A_B1R7[bn1_ma_reg3];
			 data_tmp_3  =  A_B1R8[bn1_ma_reg3];
			 data_tmp_4  =  A_B1R9[bn1_ma_reg3];
			 data_tmp_5  =  A_B1R10[bn1_ma_reg3];
			 data_tmp_6  =  A_B1R11[bn1_ma_reg3];
			 data_tmp_7  =  A_B1R12[bn1_ma_reg3];
			 data_tmp_8  =  A_B1R13[bn1_ma_reg3];
			 data_tmp_9  =  A_B1R14[bn1_ma_reg3];
			 data_tmp_10 =  A_B1R15[bn1_ma_reg3];
			 A_B1R6[bn1_ma_reg3]  =  A_B1R5[bn1_ma_reg4];
			 A_B1R7[bn1_ma_reg3]  =  A_B0R5[bn0_ma_reg4];
			 A_B1R8[bn1_ma_reg3]  =  A_B0R5[bn0_ma_reg5];
			 A_B1R9[bn1_ma_reg3]  =  A_B1R5[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg3] =  A_B1R5[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg3] =  A_B0R5[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg3] =  A_B1R5[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg3] =  A_B0R5[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg3] =  A_B0R5[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg3] =  A_B1R5[bn1_ma_reg8];
			 A_B1R5[bn1_ma_reg4]  =  data_tmp_1; 
			 A_B0R5[bn0_ma_reg4]  =  data_tmp_2; 
			 A_B0R5[bn0_ma_reg5]  =  data_tmp_3; 
			 A_B1R5[bn1_ma_reg5]  =  data_tmp_4; 
			 A_B1R5[bn1_ma_reg6]  =  data_tmp_5; 
			 A_B0R5[bn0_ma_reg6]  =  data_tmp_6; 
			 A_B1R5[bn1_ma_reg7]  =  data_tmp_7; 
			 A_B0R5[bn0_ma_reg7]  =  data_tmp_8; 
			 A_B0R5[bn0_ma_reg8]  =  data_tmp_9; 
			 A_B1R5[bn1_ma_reg8]  =  data_tmp_10;
			 //----------------------------------------------------------------
			 data_tmp_1  =  A_B1R7[bn1_ma_reg4];
			 data_tmp_2  =  A_B1R8[bn1_ma_reg4];
			 data_tmp_3  =  A_B1R9[bn1_ma_reg4];
			 data_tmp_4  =  A_B1R10[bn1_ma_reg4];
			 data_tmp_5  =  A_B1R11[bn1_ma_reg4];
			 data_tmp_6  =  A_B1R12[bn1_ma_reg4];
			 data_tmp_7  =  A_B1R13[bn1_ma_reg4];
			 data_tmp_8  =  A_B1R14[bn1_ma_reg4];
			 data_tmp_9  =  A_B1R15[bn1_ma_reg4];
			 A_B1R7[bn1_ma_reg4]  =  A_B0R6[bn0_ma_reg4];
			 A_B1R8[bn1_ma_reg4]  =  A_B0R6[bn0_ma_reg5];
			 A_B1R9[bn1_ma_reg4]  =  A_B1R6[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg4] =  A_B1R6[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg4] =  A_B0R6[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg4] =  A_B1R6[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg4] =  A_B0R6[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg4] =  A_B0R6[bn0_ma_reg8];
             A_B1R15[bn1_ma_reg4] =  A_B1R6[bn1_ma_reg8];
             A_B0R6[bn0_ma_reg4]  =  data_tmp_1;
             A_B0R6[bn0_ma_reg5]  =  data_tmp_2;
             A_B1R6[bn1_ma_reg5]  =  data_tmp_3;
             A_B1R6[bn1_ma_reg6]  =  data_tmp_4;
             A_B0R6[bn0_ma_reg6]  =  data_tmp_5;
             A_B1R6[bn1_ma_reg7]  =  data_tmp_6;
             A_B0R6[bn0_ma_reg7]  =  data_tmp_7;
             A_B0R6[bn0_ma_reg8]  =  data_tmp_8;
             A_B1R6[bn1_ma_reg8]  =  data_tmp_9;
             //------------------------------------------------------------------
			 data_tmp_1  = A_B0R8[bn0_ma_reg4];
			 data_tmp_2  = A_B0R9[bn0_ma_reg4];
			 data_tmp_3  = A_B0R10[bn0_ma_reg4];
			 data_tmp_4  = A_B0R11[bn0_ma_reg4];
			 data_tmp_5  = A_B0R12[bn0_ma_reg4];
			 data_tmp_6  = A_B0R13[bn0_ma_reg4];
			 data_tmp_7  = A_B0R14[bn0_ma_reg4];
			 data_tmp_8  = A_B0R15[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg4]  =  A_B0R7[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg4]  =  A_B1R7[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg4] =  A_B1R7[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg4] =  A_B0R7[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg4] =  A_B1R7[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg4] =  A_B0R7[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg4] =  A_B0R7[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg4] =  A_B1R7[bn1_ma_reg8];
			 A_B0R7[bn0_ma_reg5]  =  data_tmp_1;
			 A_B1R7[bn1_ma_reg5]  =  data_tmp_2;
			 A_B1R7[bn1_ma_reg6]  =  data_tmp_3;
			 A_B0R7[bn0_ma_reg6]  =  data_tmp_4;
			 A_B1R7[bn1_ma_reg7]  =  data_tmp_5;
			 A_B0R7[bn0_ma_reg7]  =  data_tmp_6;
			 A_B0R7[bn0_ma_reg8]  =  data_tmp_7;
			 A_B1R7[bn1_ma_reg8]  =  data_tmp_8;
			 //----------------------------------------------------------------
			 data_tmp_1  = A_B0R9[bn0_ma_reg5];
			 data_tmp_2  = A_B0R10[bn0_ma_reg5];
			 data_tmp_3  = A_B0R11[bn0_ma_reg5];
			 data_tmp_4  = A_B0R12[bn0_ma_reg5];
			 data_tmp_5  = A_B0R13[bn0_ma_reg5];
			 data_tmp_6  = A_B0R14[bn0_ma_reg5];
			 data_tmp_7  = A_B0R15[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg5]  =  A_B1R8[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg8];
			 A_B1R8[bn1_ma_reg5]  =  data_tmp_1;
			 A_B1R8[bn1_ma_reg6]  =  data_tmp_2;
			 A_B0R8[bn0_ma_reg6]  =  data_tmp_3;
			 A_B1R8[bn1_ma_reg7]  =  data_tmp_4;
			 A_B0R8[bn0_ma_reg7]  =  data_tmp_5;
			 A_B0R8[bn0_ma_reg8]  =  data_tmp_6;
			 A_B1R8[bn1_ma_reg8]  =  data_tmp_7;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B1R10[bn1_ma_reg5];
			 data_tmp_2 = A_B1R11[bn1_ma_reg5];
			 data_tmp_3 = A_B1R12[bn1_ma_reg5];
			 data_tmp_4 = A_B1R13[bn1_ma_reg5];
			 data_tmp_5 = A_B1R14[bn1_ma_reg5];
			 data_tmp_6 = A_B1R15[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg5] = A_B1R9[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg5] = A_B0R9[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg5] = A_B1R9[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg5] = A_B0R9[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg5] = A_B0R9[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg5] = A_B1R9[bn1_ma_reg8];
			 A_B1R9[bn1_ma_reg6]  = data_tmp_1;
			 A_B0R9[bn0_ma_reg6]  = data_tmp_2;
			 A_B1R9[bn1_ma_reg7]  = data_tmp_3;
			 A_B0R9[bn0_ma_reg7]  = data_tmp_4;
			 A_B0R9[bn0_ma_reg8]  = data_tmp_5;
			 A_B1R9[bn1_ma_reg8]  = data_tmp_6;
			 //----------------------------------------------------------------
			 data_tmp_1 =  A_B1R11[bn1_ma_reg6];
			 data_tmp_2 =  A_B1R12[bn1_ma_reg6];
			 data_tmp_3 =  A_B1R13[bn1_ma_reg6];
			 data_tmp_4 =  A_B1R14[bn1_ma_reg6];
			 data_tmp_5 =  A_B1R15[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg6]  = A_B0R10[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg6]  = A_B1R10[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg6]  = A_B0R10[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg6]  = A_B0R10[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg6]  = A_B1R10[bn1_ma_reg8];
			 A_B0R10[bn0_ma_reg6]  = data_tmp_1;
			 A_B1R10[bn1_ma_reg7]  = data_tmp_2;
			 A_B0R10[bn0_ma_reg7]  = data_tmp_3;
			 A_B0R10[bn0_ma_reg8]  = data_tmp_4;
			 A_B1R10[bn1_ma_reg8]  = data_tmp_5;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B0R12[bn0_ma_reg6];
			 data_tmp_2 = A_B0R13[bn0_ma_reg6];
			 data_tmp_3 = A_B0R14[bn0_ma_reg6];
			 data_tmp_4 = A_B0R15[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg6] = A_B1R11[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg6] = A_B0R11[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg6] = A_B0R11[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg6] = A_B1R11[bn1_ma_reg8];
			 A_B1R11[bn1_ma_reg7] = data_tmp_1;
			 A_B0R11[bn0_ma_reg7] = data_tmp_2;
			 A_B0R11[bn0_ma_reg8] = data_tmp_3;
			 A_B1R11[bn1_ma_reg8] = data_tmp_4;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B1R13[bn1_ma_reg7];
			 data_tmp_2 = A_B1R14[bn1_ma_reg7];
			 data_tmp_3 = A_B1R15[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg7] = A_B0R12[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg7] = A_B0R12[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg7] = A_B1R12[bn1_ma_reg8];
			 A_B0R12[bn0_ma_reg7] = data_tmp_1;
			 A_B0R12[bn0_ma_reg8] = data_tmp_2;
			 A_B1R12[bn1_ma_reg8] = data_tmp_3;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B0R14[bn0_ma_reg7];
			 data_tmp_2 = A_B0R15[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg7] = A_B0R13[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg7] = A_B1R13[bn1_ma_reg8];
			 A_B0R13[bn0_ma_reg8] = data_tmp_1;
			 A_B1R13[bn1_ma_reg8] = data_tmp_2;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B0R15[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg8] = A_B1R14[bn1_ma_reg8];
			 A_B1R14[bn1_ma_reg8] = data_tmp_1;
			 //-----------------------------------------------------------------
		}			
	}
	
	
	int index0_i;
	int index1_i;
	int index2_i;
	int index3_i;
	int index4_i;
	int index5_i;
	int index6_i;
	int index7_i;
	int index8_i;
	int index9_i;
	int index10_i;
	int index11_i;
	int index12_i;
	int index13_i;
	int index14_i;
	int index15_i;
	
	int index0;
    int index1;
    int index2;
    int index3;
    int index4;
    int index5;
    int index6;
    int index7;
    int index8;
    int index9;
    int index10;
    int index11;
    int index12;
    int index13;
    int index14;
    int index15;
	
	for(int i = 0; i < word_size; i++){
		MulMod(A_B0R0[i],A_B0R0[i],IN,p);
		MulMod(A_B0R1[i],A_B0R1[i],IN,p);
		MulMod(A_B0R2[i],A_B0R2[i],IN,p);
		MulMod(A_B0R3[i],A_B0R3[i],IN,p);
		MulMod(A_B0R4[i],A_B0R4[i],IN,p);
		MulMod(A_B0R5[i],A_B0R5[i],IN,p);
		MulMod(A_B0R6[i],A_B0R6[i],IN,p);
		MulMod(A_B0R7[i],A_B0R7[i],IN,p);
		MulMod(A_B0R8[i],A_B0R8[i],IN,p);
		MulMod(A_B0R9[i],A_B0R9[i],IN,p);
		MulMod(A_B0R10[i],A_B0R10[i],IN,p);
		MulMod(A_B0R11[i],A_B0R11[i],IN,p);
		MulMod(A_B0R12[i],A_B0R12[i],IN,p);
		MulMod(A_B0R13[i],A_B0R13[i],IN,p);
		MulMod(A_B0R14[i],A_B0R14[i],IN,p);
		MulMod(A_B0R15[i],A_B0R15[i],IN,p);
		//bank1 
		MulMod(A_B1R0[i],A_B1R0[i],IN,p);
		MulMod(A_B1R1[i],A_B1R1[i],IN,p);
		MulMod(A_B1R2[i],A_B1R2[i],IN,p);
		MulMod(A_B1R3[i],A_B1R3[i],IN,p);
		MulMod(A_B1R4[i],A_B1R4[i],IN,p);
		MulMod(A_B1R5[i],A_B1R5[i],IN,p);
		MulMod(A_B1R6[i],A_B1R6[i],IN,p);
		MulMod(A_B1R7[i],A_B1R7[i],IN,p);
		MulMod(A_B1R8[i],A_B1R8[i],IN,p);
		MulMod(A_B1R9[i],A_B1R9[i],IN,p);
		MulMod(A_B1R10[i],A_B1R10[i],IN,p);
		MulMod(A_B1R11[i],A_B1R11[i],IN,p);
		MulMod(A_B1R12[i],A_B1R12[i],IN,p);
		MulMod(A_B1R13[i],A_B1R13[i],IN,p);
		MulMod(A_B1R14[i],A_B1R14[i],IN,p);
		MulMod(A_B1R15[i],A_B1R15[i],IN,p);		
	}

	std::ofstream AFTER_DATA_RELOCATION("./INTT_R16_R4_After_DATA_RELOCATION.txt");
   //--------------------------------------------------------------------------------------------
	//After data relocation for INTT output
	for(int i = 0;i < group;i++){
		for(int j = 0; j < radix; j++){
        	gray_i  = Gray(i,group);
        	BC_tmp  = j * group + gray_i;
        	BC_RR = BC_tmp;
			BC_IFFT_Reorder_R16_R2(BC_RR,BC_Reorder);
        	AGU_R16(BC_Reorder,bn_tmp,ma_tmp);
			AFTER_DATA_RELOCATION <<"--------------------------------------------- \n";
			AFTER_DATA_RELOCATION <<"i: "     << i <<"\n";
			AFTER_DATA_RELOCATION <<"j: "     << j <<"\n";
			AFTER_DATA_RELOCATION <<"gray_i: "<< gray_i <<"\n";			
			AFTER_DATA_RELOCATION <<"BC_tmp: "<< BC_tmp <<"\n";
			AFTER_DATA_RELOCATION <<"BC_RR: " << BC_RR  <<"\n";
			AFTER_DATA_RELOCATION <<"BC_Reorder: "<< BC_Reorder <<"\n";
			AFTER_DATA_RELOCATION <<"bn_tmp: "<< bn_tmp <<"\n";
			AFTER_DATA_RELOCATION <<"ma_tmp: "<< ma_tmp <<"\n";			
        	if(bn_tmp == 0){
				AFTER_DATA_RELOCATION <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";						
        	}else {
				AFTER_DATA_RELOCATION <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";
				AFTER_DATA_RELOCATION <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";								
        	}	
		}
	}	

    //need to  data relocation! 
    std::ofstream INTT_DATA_OUTPUT("./INTT_R16_R4_DATA_OUTPUT.txt");	
	
	INTT_DATA_OUTPUT <<"***************************************************************\n";
	INTT_DATA_OUTPUT <<"***** DATA OUTPUT!!                                         ***\n";
	INTT_DATA_OUTPUT <<"***************************************************************\n";
	//data output
	for(int i = 0;i < group;i++){
		for(int j = 0;j < radix;j++){
			//gray_i  = Gray(i,group);
			BC_tmp  = j * group + i;
			INTT_DATA_OUTPUT <<"---------------------------------------------------\n";
			INTT_DATA_OUTPUT <<" i : "<< i <<"\n";
			INTT_DATA_OUTPUT <<" j : "<< j <<"\n";
			INTT_DATA_OUTPUT <<" BC_tmp: "<< BC_tmp <<"\n";
			BC_RR = BC_tmp;
			INTT_DATA_OUTPUT <<" BC_RR: "<< BC_RR <<"\n";
			BC_IFFT_Reorder_R16_R4_OUT(BC_RR,BC_Reorder);
			INTT_DATA_OUTPUT <<" BC_Reorder:  "<< BC_Reorder <<"\n";
			AGU_R16(BC_Reorder,bn_tmp,ma_tmp);
			INTT_DATA_OUTPUT <<" bn_tmp:  "<< bn_tmp <<"\n";
			INTT_DATA_OUTPUT <<" ma_tmp:  "<< ma_tmp <<"\n";
			index0_i  = (256 * i )+ (16 * j);
			index1_i  = (256 * i )+ (16 * j)+ 1;
			index2_i  = (256 * i )+ (16 * j)+ 2;
			index3_i  = (256 * i )+ (16 * j)+ 3;
			index4_i  = (256 * i )+ (16 * j)+ 4;
			index5_i  = (256 * i )+ (16 * j)+ 5;
			index6_i  = (256 * i )+ (16 * j)+ 6;
			index7_i  = (256 * i )+ (16 * j)+ 7;
			index8_i  = (256 * i )+ (16 * j)+ 8;
			index9_i  = (256 * i )+ (16 * j)+ 9;
			index10_i = (256 * i )+ (16 * j)+ 10;
			index11_i = (256 * i )+ (16 * j)+ 11;
			index12_i = (256 * i )+ (16 * j)+ 12;
			index13_i = (256 * i )+ (16 * j)+ 13;
			index14_i = (256 * i )+ (16 * j)+ 14;
			index15_i = (256 * i )+ (16 * j)+ 15;
			
			index0   =  index0_i;    
			index1   =  index1_i;  
			index2   =  index2_i;  
			index3   =  index3_i;  
			index4   =  index4_i;  
			index5   =  index5_i;  
			index6   =  index6_i;  
			index7   =  index7_i;  
			index8   =  index8_i;  
			index9   =  index9_i;  
			index10  =  index10_i; 
			index11  =  index11_i; 
			index12  =  index12_i; 
			index13  =  index13_i; 
			index14  =  index14_i; 
			index15  =  index15_i; 
			
			INTT_DATA_OUTPUT <<" index0_i:  "<< index0_i <<"\n";
			INTT_DATA_OUTPUT <<" index1_i:  "<< index1_i <<"\n";
			INTT_DATA_OUTPUT <<" index2_i:  "<< index2_i <<"\n";
			INTT_DATA_OUTPUT <<" index3_i:  "<< index3_i <<"\n";
			INTT_DATA_OUTPUT <<" index4_i:  "<< index4_i <<"\n";
			INTT_DATA_OUTPUT <<" index5_i:  "<< index5_i <<"\n";
			INTT_DATA_OUTPUT <<" index6_i:  "<< index6_i <<"\n";
			INTT_DATA_OUTPUT <<" index7_i:  "<< index7_i <<"\n";
			INTT_DATA_OUTPUT <<" index8_i:  "<< index8_i <<"\n";
			INTT_DATA_OUTPUT <<" index9_i:  "<< index9_i <<"\n";
			INTT_DATA_OUTPUT <<" index10_i: "<< index10_i <<"\n";
			INTT_DATA_OUTPUT <<" index11_i: "<< index11_i <<"\n";
			INTT_DATA_OUTPUT <<" index12_i: "<< index12_i <<"\n";
			INTT_DATA_OUTPUT <<" index13_i: "<< index13_i <<"\n";
			INTT_DATA_OUTPUT <<" index14_i: "<< index14_i <<"\n";
			INTT_DATA_OUTPUT <<" index15_i: "<< index15_i <<"\n";
					
			if(bn_tmp == 0){
               A[index0] = A_B0R0[ma_tmp];
			   A[index1] = A_B0R1[ma_tmp];
			   A[index2] = A_B0R2[ma_tmp];
			   A[index3] = A_B0R3[ma_tmp];
			   A[index4] = A_B0R4[ma_tmp];
			   A[index5] = A_B0R5[ma_tmp];
			   A[index6] = A_B0R6[ma_tmp];
			   A[index7] = A_B0R7[ma_tmp];
			   A[index8] = A_B0R8[ma_tmp];
			   A[index9] = A_B0R9[ma_tmp];
			   A[index10] = A_B0R10[ma_tmp];
			   A[index11] = A_B0R11[ma_tmp];
			   A[index12] = A_B0R12[ma_tmp];
			   A[index13] = A_B0R13[ma_tmp];
			   A[index14] = A_B0R14[ma_tmp];
			   A[index15] = A_B0R15[ma_tmp];
			   INTT_DATA_OUTPUT<< "A_B0R0[" << ma_tmp <<"]: "<< A_B0R0[ma_tmp] << "\n";
			   INTT_DATA_OUTPUT<< "A_B0R1[" << ma_tmp <<"]: "<< A_B0R1[ma_tmp] << "\n";
			   INTT_DATA_OUTPUT<< "A_B0R2[" << ma_tmp <<"]: "<< A_B0R2[ma_tmp] << "\n";
			   INTT_DATA_OUTPUT<< "A_B0R3[" << ma_tmp <<"]: "<< A_B0R3[ma_tmp] << "\n";
			   INTT_DATA_OUTPUT<< "A_B0R4[" << ma_tmp <<"]: "<< A_B0R4[ma_tmp] << "\n";
			   INTT_DATA_OUTPUT<< "A_B0R5[" << ma_tmp <<"]: "<< A_B0R5[ma_tmp] << "\n";
			   INTT_DATA_OUTPUT<< "A_B0R6[" << ma_tmp <<"]: "<< A_B0R6[ma_tmp] << "\n";
			   INTT_DATA_OUTPUT<< "A_B0R7[" << ma_tmp <<"]: "<< A_B0R7[ma_tmp] << "\n";
			   INTT_DATA_OUTPUT<< "A_B0R8[" << ma_tmp <<"]: "<< A_B0R8[ma_tmp] << "\n";
			   INTT_DATA_OUTPUT<< "A_B0R9[" << ma_tmp <<"]: "<< A_B0R9[ma_tmp] << "\n";
			   INTT_DATA_OUTPUT<< "A_B0R10[" << ma_tmp <<"]: "<< A_B0R10[ma_tmp] << "\n";					
			   INTT_DATA_OUTPUT<< "A_B0R11[" << ma_tmp <<"]: "<< A_B0R11[ma_tmp] << "\n";					
			   INTT_DATA_OUTPUT<< "A_B0R12[" << ma_tmp <<"]: "<< A_B0R12[ma_tmp] << "\n";					
			   INTT_DATA_OUTPUT<< "A_B0R13[" << ma_tmp <<"]: "<< A_B0R13[ma_tmp] << "\n";					
			   INTT_DATA_OUTPUT<< "A_B0R14[" << ma_tmp <<"]: "<< A_B0R14[ma_tmp] << "\n";					
			   INTT_DATA_OUTPUT<< "A_B0R15[" << ma_tmp <<"]: "<< A_B0R15[ma_tmp] <<"\n";
			}
			else {
               A[index0]     = A_B1R0[ma_tmp];
               A[index1]     = A_B1R1[ma_tmp];
               A[index2]     = A_B1R2[ma_tmp];
               A[index3]     = A_B1R3[ma_tmp];
               A[index4]     = A_B1R4[ma_tmp];
               A[index5]     = A_B1R5[ma_tmp];
               A[index6]     = A_B1R6[ma_tmp];
               A[index7]     = A_B1R7[ma_tmp];
               A[index8]     = A_B1R8[ma_tmp];
               A[index9]     = A_B1R9[ma_tmp];
               A[index10]    = A_B1R10[ma_tmp];
               A[index11]    = A_B1R11[ma_tmp];
               A[index12]    = A_B1R12[ma_tmp];
               A[index13]    = A_B1R13[ma_tmp];
               A[index14]    = A_B1R14[ma_tmp];
               A[index15]    = A_B1R15[ma_tmp];
			   INTT_DATA_OUTPUT << "A_B1R0[" << ma_tmp <<"]: "<< A_B1R0[ma_tmp] << "\n";
			   INTT_DATA_OUTPUT << "A_B1R1[" << ma_tmp <<"]: "<< A_B1R1[ma_tmp] << "\n";
			   INTT_DATA_OUTPUT << "A_B1R2[" << ma_tmp <<"]: "<< A_B1R2[ma_tmp] << "\n";
			   INTT_DATA_OUTPUT << "A_B1R3[" << ma_tmp <<"]: "<< A_B1R3[ma_tmp] << "\n";
			   INTT_DATA_OUTPUT << "A_B1R4[" << ma_tmp <<"]: "<< A_B1R4[ma_tmp] << "\n";
			   INTT_DATA_OUTPUT << "A_B1R5[" << ma_tmp <<"]: "<< A_B1R5[ma_tmp] << "\n";
			   INTT_DATA_OUTPUT << "A_B1R6[" << ma_tmp <<"]: "<< A_B1R6[ma_tmp] << "\n";
			   INTT_DATA_OUTPUT << "A_B1R7[" << ma_tmp <<"]: "<< A_B1R7[ma_tmp] << "\n";
			   INTT_DATA_OUTPUT << "A_B1R8[" << ma_tmp <<"]: "<< A_B1R8[ma_tmp] << "\n";
			   INTT_DATA_OUTPUT << "A_B1R9[" << ma_tmp <<"]: "<< A_B1R9[ma_tmp] << "\n";
			   INTT_DATA_OUTPUT << "A_B1R10[" << ma_tmp <<"]: "<< A_B1R10[ma_tmp] << "\n";					
			   INTT_DATA_OUTPUT << "A_B1R11[" << ma_tmp <<"]: "<< A_B1R11[ma_tmp] << "\n";					
			   INTT_DATA_OUTPUT << "A_B1R12[" << ma_tmp <<"]: "<< A_B1R12[ma_tmp] << "\n";					
			   INTT_DATA_OUTPUT << "A_B1R13[" << ma_tmp <<"]: "<< A_B1R13[ma_tmp] << "\n";					
			   INTT_DATA_OUTPUT << "A_B1R14[" << ma_tmp <<"]: "<< A_B1R14[ma_tmp] << "\n";					
			   INTT_DATA_OUTPUT << "A_B1R15[" << ma_tmp <<"]: "<< A_B1R15[ma_tmp] <<"\n";			   
			}			
		}		
	}
}
//mixed radix and radix-8
void NTTSPMB::NTT_r16_r8(std::vector<ZZ> &A,std::vector<ZZ> &B0R0,std::vector<ZZ> &B0R1,std::vector<ZZ> &B0R2,std::vector<ZZ> &B0R3,
std::vector<ZZ> &B0R4,std::vector<ZZ> &B0R5,std::vector<ZZ> &B0R6,std::vector<ZZ> &B0R7,std::vector<ZZ> &B0R8,std::vector<ZZ> &B0R9,
std::vector<ZZ> &B0R10,std::vector<ZZ> &B0R11,std::vector<ZZ> &B0R12,std::vector<ZZ> &B0R13,std::vector<ZZ> &B0R14,std::vector<ZZ> &B0R15,
std::vector<ZZ> &B1R0,std::vector<ZZ> &B1R1,std::vector<ZZ> &B1R2,std::vector<ZZ> &B1R3,std::vector<ZZ> &B1R4,std::vector<ZZ> &B1R5,
std::vector<ZZ> &B1R6,std::vector<ZZ> &B1R7,std::vector<ZZ> &B1R8,std::vector<ZZ> &B1R9,std::vector<ZZ> &B1R10,std::vector<ZZ> &B1R11,
std::vector<ZZ> &B1R12,std::vector<ZZ> &B1R13,std::vector<ZZ> &B1R14,std::vector<ZZ> &B1R15){
	
	int            Stage;//stage	
	int            word_size;
	int            offset;
	int            BC;
	int            BC_tmp;
    int            group;
    int            ma_tmp;
	int            bn_tmp;
    int            bit_tmp;
	int            bn0_bc_tmp;//frist in bc data
	int            bn1_bc_tmp;//frist in bc data
	int            bn0_ma_reg1;
	int            bn0_ma_reg2;
	int            bn0_ma_reg3;
	int            bn0_ma_reg4;
	int            bn0_ma_reg5;
	int            bn0_ma_reg6;
	int            bn0_ma_reg7;
	int            bn0_ma_reg8;
	int            bn1_ma_reg1;
	int            bn1_ma_reg2;
	int            bn1_ma_reg3;
	int            bn1_ma_reg4;
	int            bn1_ma_reg5;
	int            bn1_ma_reg6;
	int            bn1_ma_reg7;
	int            bn1_ma_reg8;
	int            gray_i;
    int            BC_WIDTH;
	int            tw_modulus;
    int            tw_modulus_tmp;
	//display 
	int            display;
    double         Stage_double;	
    std::vector<int> bit_array_tmp;
	//---------------------------------------------
	//display parameter , 1: display , 0:not display
	display = 1;
	//--------------------------------------------
	std::ofstream DATARECORD("./NTT_R16_R8_SPMB.txt");
	std::ofstream siang_record("./my_print_data/r16_r8_SPMB.txt");
	std::ofstream r16_r8_SPMB_TF_dg("./my_print_data/r16_r8_SPMB_TF_dg.txt");
	
	std::ofstream spmb_r16_r8("./SPMB_tw/spmb_r16_r8.txt");
	std::ofstream DTFAG_golden_st0("./SPMB_tw/DTFAG_golden_st0_32768.txt");
	std::ofstream DTFAG_golden_st1("./SPMB_tw/DTFAG_golden_st1_32768.txt");
	std::ofstream DTFAG_golden_st2("./SPMB_tw/DTFAG_golden_st2_32768.txt");
	
	//-------------Bit Reverse--------------
	BitOperate RR, IntToVec, VecToInt;
	long long bit_width = log2(N);
	long long bit_width_s = log(radix)/log(2);
	vector<long long> bit_array;
	//----------------------------------------
	
	//radix-16 Stage
    Stage_double  = log2(N);
    Stage       =  (int)floor(Stage_double/4);
	BC_WIDTH    =  (int)ceil(log2(N/16));
	offset      =  (int)N /  16;
	word_size   =  (int)N / (2 * 16);
	group       =  (int)N / (256);
	tw_modulus  =  (int)N /  16;
	bit_array_tmp.resize(BC_WIDTH);

    if(display == 1)DATARECORD << "group: "    << group << "\n";
    if(display == 1)DATARECORD << "BC_WIDTH: " << BC_WIDTH << "\n";
	siang_record << "group: "    << group << "\n";	// siang record
    siang_record << "BC_WIDTH: " << BC_WIDTH << "\n";	// siagn record
	
	ZZ               data_tmp_1;
	ZZ               data_tmp_2;
	ZZ               data_tmp_3;
	ZZ               data_tmp_4;
	ZZ               data_tmp_5;
	ZZ               data_tmp_6;
	ZZ               data_tmp_7;
	ZZ               data_tmp_8;
	ZZ               data_tmp_9;
	ZZ               data_tmp_10;
	ZZ               data_tmp_11;
	ZZ               data_tmp_12;
	ZZ               data_tmp_13;
	ZZ               data_tmp_14;
	ZZ               data_tmp_15;
	//data tmp for index 8 of each array
	//16 * 8 matrix
    std::vector<ZZ>  data_tmp_array_1;  
    std::vector<ZZ>  data_tmp_array_2;  
    std::vector<ZZ>  data_tmp_array_3;  
    std::vector<ZZ>  data_tmp_array_4;  
    std::vector<ZZ>  data_tmp_array_5;  
    std::vector<ZZ>  data_tmp_array_6;  
    std::vector<ZZ>  data_tmp_array_7;  
    std::vector<ZZ>  data_tmp_array_8;

    data_tmp_array_1.resize(16);
    data_tmp_array_2.resize(16);
    data_tmp_array_3.resize(16);
    data_tmp_array_4.resize(16);
    data_tmp_array_5.resize(16);
    data_tmp_array_6.resize(16);
	data_tmp_array_7.resize(16);
	data_tmp_array_8.resize(16);
	//---------------------------------
	std::vector<ZZ>  A_B0R0;
	std::vector<ZZ>  A_B0R1;
	std::vector<ZZ>  A_B0R2;
	std::vector<ZZ>  A_B0R3;
	std::vector<ZZ>  A_B0R4;
	std::vector<ZZ>  A_B0R5;
	std::vector<ZZ>  A_B0R6;
	std::vector<ZZ>  A_B0R7;
	std::vector<ZZ>  A_B0R8;
	std::vector<ZZ>  A_B0R9;
	std::vector<ZZ>  A_B0R10;
	std::vector<ZZ>  A_B0R11;
	std::vector<ZZ>  A_B0R12;
	std::vector<ZZ>  A_B0R13;
	std::vector<ZZ>  A_B0R14;
	std::vector<ZZ>  A_B0R15;
	std::vector<ZZ>  A_B1R0;
	std::vector<ZZ>  A_B1R1;
	std::vector<ZZ>  A_B1R2;
	std::vector<ZZ>  A_B1R3;
	std::vector<ZZ>  A_B1R4;
	std::vector<ZZ>  A_B1R5;
	std::vector<ZZ>  A_B1R6;
	std::vector<ZZ>  A_B1R7;
	std::vector<ZZ>  A_B1R8;
	std::vector<ZZ>  A_B1R9;
	std::vector<ZZ>  A_B1R10;
	std::vector<ZZ>  A_B1R11;
	std::vector<ZZ>  A_B1R12;
	std::vector<ZZ>  A_B1R13;
	std::vector<ZZ>  A_B1R14;
	std::vector<ZZ>  A_B1R15;
	
	A_B0R0.resize(word_size);
	A_B0R1.resize(word_size);
	A_B0R2.resize(word_size);
	A_B0R3.resize(word_size);
	A_B0R4.resize(word_size);
	A_B0R5.resize(word_size);
	A_B0R6.resize(word_size);
	A_B0R7.resize(word_size);
	A_B0R8.resize(word_size);
	A_B0R9.resize(word_size);
	A_B0R10.resize(word_size);
	A_B0R11.resize(word_size);
	A_B0R12.resize(word_size);
	A_B0R13.resize(word_size);
	A_B0R14.resize(word_size);
	A_B0R15.resize(word_size);
	A_B1R0.resize(word_size);
	A_B1R1.resize(word_size);
	A_B1R2.resize(word_size);
	A_B1R3.resize(word_size);
	A_B1R4.resize(word_size);
	A_B1R5.resize(word_size);
	A_B1R6.resize(word_size);
	A_B1R7.resize(word_size);
	A_B1R8.resize(word_size);
	A_B1R9.resize(word_size);
	A_B1R10.resize(word_size);
	A_B1R11.resize(word_size);
	A_B1R12.resize(word_size);
	A_B1R13.resize(word_size);
	A_B1R14.resize(word_size);
	A_B1R15.resize(word_size);
	//----------------------------------------------------
	B0R0.resize(word_size);
	B0R1.resize(word_size);
	B0R2.resize(word_size);
	B0R3.resize(word_size);
	B0R4.resize(word_size);
	B0R5.resize(word_size);
	B0R6.resize(word_size);
	B0R7.resize(word_size);
	B0R8.resize(word_size);
	B0R9.resize(word_size);
	B0R10.resize(word_size);
	B0R11.resize(word_size);
	B0R12.resize(word_size);
	B0R13.resize(word_size);
	B0R14.resize(word_size);
	B0R15.resize(word_size);
	B1R0.resize(word_size);
	B1R1.resize(word_size);
	B1R2.resize(word_size);
	B1R3.resize(word_size);
	B1R4.resize(word_size);
	B1R5.resize(word_size);
	B1R6.resize(word_size);
	B1R7.resize(word_size);
	B1R8.resize(word_size);
	B1R9.resize(word_size);
	B1R10.resize(word_size);
	B1R11.resize(word_size);
	B1R12.resize(word_size);
	B1R13.resize(word_size);
	B1R14.resize(word_size);
	B1R15.resize(word_size);
	//----------------------------------------------------
	int length;
	ZZ  factor;   //base factor
	ZZ  factor_t; //acctually mul factor
	ZZ  factor_2t;
	ZZ  factor_3t;
	ZZ  factor_4t;
	ZZ  factor_5t;
	ZZ  factor_6t;
	ZZ  factor_7t;
	ZZ  factor_8t;
	ZZ  factor_9t;
	ZZ  factor_10t;
	ZZ  factor_11t;
	ZZ  factor_12t;
	ZZ  factor_13t;
	ZZ  factor_14t;
	ZZ  factor_15t;
    
	//init load data
    for(int i = 0; i < group; i++){
		for(int j = 0 ; j < radix ; j++){
			bn_tmp = 0;
			ma_tmp = 0;
			gray_i  = Gray(i,group);
			BC = j * group + gray_i;
            //bit calculate
            for(int j=0; j < BC_WIDTH;j++){
                bit_tmp = BC % 2;
                BC = BC >> 1;
                bit_array_tmp[j] = bit_tmp;
            } 
			gray_i  = Gray(i,group);
			BC = j * group + gray_i;
            for(int rs = 0; rs < BC_WIDTH; rs++){
                if((bit_array_tmp[rs] == 1) && (rs != 0)) ma_tmp = ma_tmp + exp2((rs-1)); 
                bn_tmp = bn_tmp ^ bit_array_tmp[rs];
            }
			if(bn_tmp == 0){
				A_B0R0[ma_tmp]  = A[BC];
				A_B0R1[ma_tmp]  = A[BC +      offset];
				A_B0R2[ma_tmp]  = A[BC + 2  * offset];
				A_B0R3[ma_tmp]  = A[BC + 3  * offset];
				A_B0R4[ma_tmp]  = A[BC + 4  * offset];
				A_B0R5[ma_tmp]  = A[BC + 5  * offset];
				A_B0R6[ma_tmp]  = A[BC + 6  * offset];
				A_B0R7[ma_tmp]  = A[BC + 7  * offset];
				A_B0R8[ma_tmp]  = A[BC + 8  * offset];
				A_B0R9[ma_tmp]  = A[BC + 9  * offset];
				A_B0R10[ma_tmp] = A[BC + 10 * offset];
				A_B0R11[ma_tmp] = A[BC + 11 * offset];
				A_B0R12[ma_tmp] = A[BC + 12 * offset];
				A_B0R13[ma_tmp] = A[BC + 13 * offset];
				A_B0R14[ma_tmp] = A[BC + 14 * offset];
				A_B0R15[ma_tmp] = A[BC + 15 * offset];
			}else {
				A_B1R0[ma_tmp]  = A[BC];
				A_B1R1[ma_tmp]  = A[BC +     offset];
				A_B1R2[ma_tmp]  = A[BC + 2 * offset];
				A_B1R3[ma_tmp]  = A[BC + 3 * offset];
				A_B1R4[ma_tmp]  = A[BC + 4 * offset];
				A_B1R5[ma_tmp]  = A[BC + 5 * offset];
				A_B1R6[ma_tmp]  = A[BC + 6 * offset];
				A_B1R7[ma_tmp]  = A[BC + 7 * offset];
				A_B1R8[ma_tmp]  = A[BC + 8 * offset];
				A_B1R9[ma_tmp]  = A[BC + 9 * offset];
				A_B1R10[ma_tmp] = A[BC + 10 * offset];
				A_B1R11[ma_tmp] = A[BC + 11 * offset];
				A_B1R12[ma_tmp] = A[BC + 12 * offset];
				A_B1R13[ma_tmp] = A[BC + 13 * offset];
				A_B1R14[ma_tmp] = A[BC + 14 * offset];
				A_B1R15[ma_tmp] = A[BC + 15 * offset];
			}
		}
	}
	ma_tmp = 0;
	bn_tmp = 0;
	BC     = 0;
	std::cout << "init load over! \n";
	int tw_degree = 1;
	if(display == 1)DATARECORD <<"Stage: " << Stage << "\n";	
	if(display == 1)siang_record <<"radix-16 computing stage:  "<< Stage <<"\n"; // siang_record	
	//need modify to mult by twiddle factor
	for(int s = 0; s < Stage; s++){
		if(s == 0)factor = W;
		else {
			SqrMod(factor,factor,p);
			SqrMod(factor,factor,p);
			SqrMod(factor,factor,p);
			SqrMod(factor,factor,p);
			tw_degree = tw_degree * 16;
		}
		std::cout << "factor = " << factor << "\n";
		if(display == 1)DATARECORD <<"**************************************\n";					
		if(display == 1)siang_record <<"---------------------------------\n";	// siang_record
		if(display == 1)siang_record <<"Now Stage: "<< s <<"\n";	// siang_record
		if(display == 1)siang_record <<"factor: "<< factor <<"\n";	// siang_record	
		r16_r8_SPMB_TF_dg << "NOW stage: " << s << "!!!!\n";					
		spmb_r16_r8 <<"Now Stage: "<< s <<"\n";
		spmb_r16_r8 <<"twiddle factor : "<< factor <<"\n";
		tw_modulus_tmp  = tw_modulus >> ( 4 * s);
		for(int i = 0 ;i < group;i++){
			bn0_bc_tmp  = 0;
			bn1_bc_tmp  = 0;
			if(display == 1)DATARECORD <<"-------------------------------------\n";										
			if(display == 1)DATARECORD <<"NOW stage: " << s << "!!!!\n";	
			r16_r8_SPMB_TF_dg <<"-------------------------------------\n";							
			spmb_r16_r8 << "-----------------i = " << i << "-------------------" << std::endl;
			for(int j = 0;j < radix;j++){
				siang_record << "---------------------------------------" << std::endl;
				siang_record << "s = " << s << ", i = " << i << ", j = " << j << std::endl;
				gray_i  = Gray(i,group);
			    BC_tmp  = j * group + gray_i;
				if(display == 1)DATARECORD <<"---------------------------------\n";
				if(display == 1)DATARECORD << "i = " << i << ", j = " << j << std::endl;
				if(display == 1)DATARECORD << "BC_tmp: " << BC_tmp << "\n";					
				if(s == Stage - 1) RR_R16_R8(BC_tmp,(4 * s - 1),BC);
				else RR_R16(BC_tmp,s,BC);
				if(display == 1)DATARECORD << "After RR_R16 , BC : " << BC << "\n";		
				length = BC % tw_modulus_tmp;
				siang_record << "length: " << length << std::endl;
				r16_r8_SPMB_TF_dg << "BC = " << BC << ", tw_modulus_tmp = " << tw_modulus_tmp << ", tw_modulus = " << tw_modulus
				<< ", length : "<< length << ", tw_dg: " <<  tw_degree * length << " = " << tw_degree << " * " << length <<"\n";
				PowerMod(factor_t,factor,length,p);
				AGU_R16(BC,bn_tmp,ma_tmp);
				if(display == 1)DATARECORD << "BN : " << bn_tmp << "\n";
				if(display == 1)DATARECORD << "MA : " << ma_tmp << "\n";	
				
				//-----------compute data idx-------------
				spmb_r16_r8 /*<< "BC = " << BC */<< "BC_tmp = " << BC_tmp << std::endl;
				
				IntToVec.IntToVec(BC_tmp, N, bit_array);
				rotate(bit_array.begin(), bit_array.begin()+ bit_width_s*s , bit_array.end());
				long long Data = VecToInt.VecToInt(bit_array, N);
				spmb_r16_r8 << "Data_index = ";
                spmb_r16_r8 << "( " ;					
				for(int k = 0; k < radix ; k++ ){
					spmb_r16_r8 << Data + k*(1<<(bit_width-bit_width_s-bit_width_s*s)) <<" ";	
				}
				spmb_r16_r8 << ") " << std::endl;
				spmb_r16_r8 << ", (w^" << 0 << ", w^" << tw_degree * length << ", w^" << tw_degree * length * 2 
				<< ", w^" << tw_degree * length * 3  << ", w^" << tw_degree * length * 4 << ", w^" << tw_degree * length * 5 
				<< ", w^" << tw_degree * length * 6  << ", w^" << tw_degree * length * 7 << ", w^" << tw_degree * length * 8 
				<< ", w^" << tw_degree * length * 9  << ", w^" << tw_degree * length * 10 << ", w^" << tw_degree * length * 11 
				<< ", w^" << tw_degree * length * 12  << ", w^" << tw_degree * length * 13 << ", w^" << tw_degree * length * 14 
				<< ", w^" << tw_degree * length * 15
				<< ")" <<std::endl;
				//-----------------------------------------
				
				if(bn_tmp == 0){
					if(j < 2)bn0_bc_tmp = BC_tmp;
                    //if(display == 1)DATARECORD <<" Before radix-16 BU operation!! \n";					
                    //if(display == 1)DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";				
					PowerMod(factor_2t,factor_t,2,p);
					PowerMod(factor_3t,factor_t,3,p);
					PowerMod(factor_4t,factor_t,4,p);
					PowerMod(factor_5t,factor_t,5,p);
					PowerMod(factor_6t,factor_t,6,p);
					PowerMod(factor_7t,factor_t,7,p);
					PowerMod(factor_8t,factor_t,8,p);
					PowerMod(factor_9t,factor_t,9,p);
					PowerMod(factor_10t,factor_t,10,p);
					PowerMod(factor_11t,factor_t,11,p);
					PowerMod(factor_12t,factor_t,12,p);
					PowerMod(factor_13t,factor_t,13,p);
					PowerMod(factor_14t,factor_t,14,p);
					PowerMod(factor_15t,factor_t,15,p);

					if(display == 1)DATARECORD <<"factor_t: "  << factor_t <<"\n";
					if(display == 1)DATARECORD <<"factor_2t: " << factor_2t <<"\n";
					if(display == 1)DATARECORD <<"factor_3t: " << factor_3t <<"\n";
					if(display == 1)DATARECORD <<"factor_4t: " << factor_4t <<"\n";
					if(display == 1)DATARECORD <<"factor_5t: " << factor_5t <<"\n";
					if(display == 1)DATARECORD <<"factor_6t: " << factor_6t <<"\n";
					if(display == 1)DATARECORD <<"factor_7t: " << factor_7t <<"\n";
					if(display == 1)DATARECORD <<"factor_8t: " << factor_8t <<"\n";
					if(display == 1)DATARECORD <<"factor_9t: " << factor_9t <<"\n";
					if(display == 1)DATARECORD <<"factor_10t: "<< factor_10t <<"\n";
					if(display == 1)DATARECORD <<"factor_11t: "<< factor_11t <<"\n";
					if(display == 1)DATARECORD <<"factor_12t: "<< factor_12t <<"\n";
					if(display == 1)DATARECORD <<"factor_13t: "<< factor_13t <<"\n";
					if(display == 1)DATARECORD <<"factor_14t: "<< factor_14t <<"\n";
					if(display == 1)DATARECORD <<"factor_15t: "<< factor_15t <<"\n";	

					//----------------------siang-------------
					if(display == 1)siang_record << "p : " << p << "\n";	// siang_record
					if(display == 1)siang_record << "factor_t : "   		   << factor_t << "\n";	// siang_record	
					if(display == 1)siang_record << "factor_2t_factor_3t : "   << factor_2t  << "_" << factor_3t << "\n";	// siang_record	
					if(display == 1)siang_record << "factor_4t_factor_5t : "   << factor_4t  << "_" << factor_5t << "\n";	// siang_record	
					if(display == 1)siang_record << "factor_6t_factor_7t : "   << factor_6t  << "_" << factor_7t << "\n";	// siang_record	
					if(display == 1)siang_record << "factor_8t_factor_9t : "   << factor_8t  << "_" << factor_9t << "\n";	// siang_record	
					if(display == 1)siang_record << "factor_10t_factor_11t : " << factor_10t << "_" << factor_11t << "\n";	// siang_record	
					if(display == 1)siang_record << "factor_12t_factor_13t : " << factor_12t << "_" << factor_13t << "\n";	// siang_record	
					if(display == 1)siang_record << "factor_14t_factor_15t : " << factor_14t << "_" << factor_15t << "\n";	// siang_record	

					//----------------------------------------
					Radix16_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp], A_B0R2[ma_tmp], A_B0R3[ma_tmp],A_B0R4[ma_tmp],
							   A_B0R5[ma_tmp],A_B0R6[ma_tmp], A_B0R7[ma_tmp], A_B0R8[ma_tmp],A_B0R9[ma_tmp],
							   A_B0R10[ma_tmp],A_B0R11[ma_tmp],A_B0R12[ma_tmp],A_B0R13[ma_tmp],A_B0R14[ma_tmp],
							   A_B0R15[ma_tmp]);
                    //if(display == 1)DATARECORD <<" After radix-16 BU Operation!!!\n";												   
                    //if(display == 1)DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";												   
					MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],factor_t,p);
					MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],factor_2t,p);
					MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],factor_3t,p);
					MulMod(A_B0R4[ma_tmp],A_B0R4[ma_tmp],factor_4t,p);
					MulMod(A_B0R5[ma_tmp],A_B0R5[ma_tmp],factor_5t,p);
					MulMod(A_B0R6[ma_tmp],A_B0R6[ma_tmp],factor_6t,p);
					MulMod(A_B0R7[ma_tmp],A_B0R7[ma_tmp],factor_7t,p);
					MulMod(A_B0R8[ma_tmp],A_B0R8[ma_tmp],factor_8t,p);
					MulMod(A_B0R9[ma_tmp],A_B0R9[ma_tmp],factor_9t,p);
					MulMod(A_B0R10[ma_tmp],A_B0R10[ma_tmp],factor_10t,p);
					MulMod(A_B0R11[ma_tmp],A_B0R11[ma_tmp],factor_11t,p);
					MulMod(A_B0R12[ma_tmp],A_B0R12[ma_tmp],factor_12t,p);
					MulMod(A_B0R13[ma_tmp],A_B0R13[ma_tmp],factor_13t,p);
					MulMod(A_B0R14[ma_tmp],A_B0R14[ma_tmp],factor_14t,p);
					MulMod(A_B0R15[ma_tmp],A_B0R15[ma_tmp],factor_15t,p);
                    //if(display == 1)DATARECORD <<" After multiplication!!!\n";					
                    //if(display == 1)DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";											
					if(j <  2)bn0_ma_reg1 = ma_tmp;
					if((j >= 2)  && (j < 4))bn0_ma_reg2 = ma_tmp;
					if((j >= 4)  && (j < 6))bn0_ma_reg3 = ma_tmp;
					if((j >= 6)  && (j < 8))bn0_ma_reg4 = ma_tmp;
					if((j >= 8)  && (j < 10))bn0_ma_reg5 = ma_tmp;
					if((j >= 10) && (j < 12))bn0_ma_reg6 = ma_tmp;
					if((j >= 12) && (j < 14))bn0_ma_reg7 = ma_tmp;
					if((j >= 14) && (j < 16))bn0_ma_reg8 = ma_tmp;
					//----------------DTFAG golden pattern------------------
					switch(s){
						case 0:
							DTFAG_golden_st0 << 1 		  << " \n ";
							DTFAG_golden_st0 << factor_t  << " \n ";
							DTFAG_golden_st0 << factor_2t << " \n ";
							DTFAG_golden_st0 << factor_3t << " \n ";
							DTFAG_golden_st0 << factor_4t << " \n ";
							DTFAG_golden_st0 << factor_5t  << " \n ";
							DTFAG_golden_st0 << factor_6t << " \n ";
							DTFAG_golden_st0 << factor_7t << " \n ";
							DTFAG_golden_st0 << factor_8t << " \n ";
							DTFAG_golden_st0 << factor_9t  << " \n ";
							DTFAG_golden_st0 << factor_10t << " \n ";
							DTFAG_golden_st0 << factor_11t << " \n ";
							DTFAG_golden_st0 << factor_12t << " \n ";
							DTFAG_golden_st0 << factor_13t  << " \n ";
							DTFAG_golden_st0 << factor_14t << " \n ";
							DTFAG_golden_st0 << factor_15t << " \n ";
							break;
						case 1:
							DTFAG_golden_st1 << 1 		  << " \n ";
							DTFAG_golden_st1 << factor_t  << " \n ";
							DTFAG_golden_st1 << factor_2t << " \n ";
							DTFAG_golden_st1 << factor_3t << " \n ";
							DTFAG_golden_st1 << factor_4t << " \n ";
							DTFAG_golden_st1 << factor_5t  << " \n ";
							DTFAG_golden_st1 << factor_6t << " \n ";
							DTFAG_golden_st1 << factor_7t << " \n ";
							DTFAG_golden_st1 << factor_8t << " \n ";
							DTFAG_golden_st1 << factor_9t  << " \n ";
							DTFAG_golden_st1 << factor_10t << " \n ";
							DTFAG_golden_st1 << factor_11t << " \n ";
							DTFAG_golden_st1 << factor_12t << " \n ";
							DTFAG_golden_st1 << factor_13t  << " \n ";
							DTFAG_golden_st1 << factor_14t << " \n ";
							DTFAG_golden_st1 << factor_15t << " \n ";
							break;
						case 2:
							DTFAG_golden_st2 << 1 		  << " \n ";
							DTFAG_golden_st2 << factor_t  << " \n ";
							DTFAG_golden_st2 << factor_2t << " \n ";
							DTFAG_golden_st2 << factor_3t << " \n ";
							DTFAG_golden_st2 << factor_4t << " \n ";
							DTFAG_golden_st2 << factor_5t  << " \n ";
							DTFAG_golden_st2 << factor_6t << " \n ";
							DTFAG_golden_st2 << factor_7t << " \n ";
							DTFAG_golden_st2 << factor_8t << " \n ";
							DTFAG_golden_st2 << factor_9t  << " \n ";
							DTFAG_golden_st2 << factor_10t << " \n ";
							DTFAG_golden_st2 << factor_11t << " \n ";
							DTFAG_golden_st2 << factor_12t << " \n ";
							DTFAG_golden_st2 << factor_13t  << " \n ";
							DTFAG_golden_st2 << factor_14t << " \n ";
							DTFAG_golden_st2 << factor_15t << " \n ";
							break;
						default:
							break;
					}
					//--------------------------------------------------
				}
			    else {
					if(j < 2)bn1_bc_tmp = BC_tmp;
                    //if(display == 1)DATARECORD <<" Before radix-16 BU operation!! \n";					
                    //if(display == 1)DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";						
					PowerMod(factor_2t,factor_t,2,p);
					PowerMod(factor_3t,factor_t,3,p);
					PowerMod(factor_4t,factor_t,4,p);
					PowerMod(factor_5t,factor_t,5,p);
					PowerMod(factor_6t,factor_t,6,p);
					PowerMod(factor_7t,factor_t,7,p);
					PowerMod(factor_8t,factor_t,8,p);
					PowerMod(factor_9t,factor_t,9,p);
					PowerMod(factor_10t,factor_t,10,p);
					PowerMod(factor_11t,factor_t,11,p);
					PowerMod(factor_12t,factor_t,12,p);
					PowerMod(factor_13t,factor_t,13,p);
					PowerMod(factor_14t,factor_t,14,p);
					PowerMod(factor_15t,factor_t,15,p);	
					if(display == 1)DATARECORD <<"factor_t: "  << factor_t <<"\n";
					if(display == 1)DATARECORD <<"factor_2t: " << factor_2t <<"\n";
					if(display == 1)DATARECORD <<"factor_3t: " << factor_3t <<"\n";
					if(display == 1)DATARECORD <<"factor_4t: " << factor_4t <<"\n";
					if(display == 1)DATARECORD <<"factor_5t: " << factor_5t <<"\n";
					if(display == 1)DATARECORD <<"factor_6t: " << factor_6t <<"\n";
					if(display == 1)DATARECORD <<"factor_7t: " << factor_7t <<"\n";
					if(display == 1)DATARECORD <<"factor_8t: " << factor_8t <<"\n";
					if(display == 1)DATARECORD <<"factor_9t: " << factor_9t <<"\n";
					if(display == 1)DATARECORD <<"factor_10t: "<< factor_10t <<"\n";
					if(display == 1)DATARECORD <<"factor_11t: "<< factor_11t <<"\n";
					if(display == 1)DATARECORD <<"factor_12t: "<< factor_12t <<"\n";
					if(display == 1)DATARECORD <<"factor_13t: "<< factor_13t <<"\n";
					if(display == 1)DATARECORD <<"factor_14t: "<< factor_14t <<"\n";
					if(display == 1)DATARECORD <<"factor_15t: "<< factor_15t <<"\n";	
					//------------------siang record-----------------------
					if(display == 1)siang_record << "p : " << p << "\n";	// siang_record
					if(display == 1)siang_record << "factor_t : "   		   << factor_t << "\n";	// siang_record	
					if(display == 1)siang_record << "factor_2t_factor_3t : "   << factor_2t  << "_" << factor_3t << "\n";	// siang_record	
					if(display == 1)siang_record << "factor_4t_factor_5t : "   << factor_4t  << "_" << factor_5t << "\n";	// siang_record	
					if(display == 1)siang_record << "factor_6t_factor_7t : "   << factor_6t  << "_" << factor_7t << "\n";	// siang_record	
					if(display == 1)siang_record << "factor_8t_factor_9t : "   << factor_8t  << "_" << factor_9t << "\n";	// siang_record	
					if(display == 1)siang_record << "factor_10t_factor_11t : " << factor_10t << "_" << factor_11t << "\n";	// siang_record	
					if(display == 1)siang_record << "factor_12t_factor_13t : " << factor_12t << "_" << factor_13t << "\n";	// siang_record	
					if(display == 1)siang_record << "factor_14t_factor_15t : " << factor_14t << "_" << factor_15t << "\n";	// siang_record	
					//-----------------------------------------------------				
					Radix16_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp], A_B1R2[ma_tmp], A_B1R3[ma_tmp],A_B1R4[ma_tmp],
							   A_B1R5[ma_tmp],A_B1R6[ma_tmp], A_B1R7[ma_tmp], A_B1R8[ma_tmp],A_B1R9[ma_tmp],
							   A_B1R10[ma_tmp],A_B1R11[ma_tmp],A_B1R12[ma_tmp],A_B1R13[ma_tmp],A_B1R14[ma_tmp],
							   A_B1R15[ma_tmp]);
                    //if(display == 1)DATARECORD <<" After radix-16 BU Operation!!!\n";												   
                    //if(display == 1)DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";								   
					MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],factor_t,p);
					MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],factor_2t,p);
					MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],factor_3t,p);
					MulMod(A_B1R4[ma_tmp],A_B1R4[ma_tmp],factor_4t,p);
					MulMod(A_B1R5[ma_tmp],A_B1R5[ma_tmp],factor_5t,p);
					MulMod(A_B1R6[ma_tmp],A_B1R6[ma_tmp],factor_6t,p);
					MulMod(A_B1R7[ma_tmp],A_B1R7[ma_tmp],factor_7t,p);
					MulMod(A_B1R8[ma_tmp],A_B1R8[ma_tmp],factor_8t,p);
					MulMod(A_B1R9[ma_tmp],A_B1R9[ma_tmp],factor_9t,p);
					MulMod(A_B1R10[ma_tmp],A_B1R10[ma_tmp],factor_10t,p);
					MulMod(A_B1R11[ma_tmp],A_B1R11[ma_tmp],factor_11t,p);
					MulMod(A_B1R12[ma_tmp],A_B1R12[ma_tmp],factor_12t,p);
					MulMod(A_B1R13[ma_tmp],A_B1R13[ma_tmp],factor_13t,p);
					MulMod(A_B1R14[ma_tmp],A_B1R14[ma_tmp],factor_14t,p);
					MulMod(A_B1R15[ma_tmp],A_B1R15[ma_tmp],factor_15t,p);
                    //if(display == 1)DATARECORD <<" After multiplication!!!\n";					
                    //if(display == 1)DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";					
                    //if(display == 1)DATARECORD <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";						
                    if(j <  2)bn1_ma_reg1 = ma_tmp;					
                    if((j >= 2)  && (j < 4))bn1_ma_reg2 = ma_tmp;
                    if((j >= 4)  && (j < 6))bn1_ma_reg3 = ma_tmp;
                    if((j >= 6)  && (j < 8))bn1_ma_reg4 = ma_tmp;
                    if((j >= 8)  && (j < 10))bn1_ma_reg5 = ma_tmp;
                    if((j >= 10) && (j < 12))bn1_ma_reg6 = ma_tmp;
                    if((j >= 12) && (j < 14))bn1_ma_reg7 = ma_tmp;
                    if((j >= 14) && (j < 16))bn1_ma_reg8 = ma_tmp;
					//----------------DTFAG golden pattern------------------
					switch(s){
						case 0:
							DTFAG_golden_st0 << 1 		  << " \n ";
							DTFAG_golden_st0 << factor_t  << " \n ";
							DTFAG_golden_st0 << factor_2t << " \n ";
							DTFAG_golden_st0 << factor_3t << " \n ";
							DTFAG_golden_st0 << factor_4t << " \n ";
							DTFAG_golden_st0 << factor_5t  << " \n ";
							DTFAG_golden_st0 << factor_6t << " \n ";
							DTFAG_golden_st0 << factor_7t << " \n ";
							DTFAG_golden_st0 << factor_8t << " \n ";
							DTFAG_golden_st0 << factor_9t  << " \n ";
							DTFAG_golden_st0 << factor_10t << " \n ";
							DTFAG_golden_st0 << factor_11t << " \n ";
							DTFAG_golden_st0 << factor_12t << " \n ";
							DTFAG_golden_st0 << factor_13t  << " \n ";
							DTFAG_golden_st0 << factor_14t << " \n ";
							DTFAG_golden_st0 << factor_15t << " \n ";
							break;
						case 1:
							DTFAG_golden_st1 << 1 		  << " \n ";
							DTFAG_golden_st1 << factor_t  << " \n ";
							DTFAG_golden_st1 << factor_2t << " \n ";
							DTFAG_golden_st1 << factor_3t << " \n ";
							DTFAG_golden_st1 << factor_4t << " \n ";
							DTFAG_golden_st1 << factor_5t  << " \n ";
							DTFAG_golden_st1 << factor_6t << " \n ";
							DTFAG_golden_st1 << factor_7t << " \n ";
							DTFAG_golden_st1 << factor_8t << " \n ";
							DTFAG_golden_st1 << factor_9t  << " \n ";
							DTFAG_golden_st1 << factor_10t << " \n ";
							DTFAG_golden_st1 << factor_11t << " \n ";
							DTFAG_golden_st1 << factor_12t << " \n ";
							DTFAG_golden_st1 << factor_13t  << " \n ";
							DTFAG_golden_st1 << factor_14t << " \n ";
							DTFAG_golden_st1 << factor_15t << " \n ";
							break;
						case 2:
							DTFAG_golden_st2 << 1 		  << " \n ";
							DTFAG_golden_st2 << factor_t  << " \n ";
							DTFAG_golden_st2 << factor_2t << " \n ";
							DTFAG_golden_st2 << factor_3t << " \n ";
							DTFAG_golden_st2 << factor_4t << " \n ";
							DTFAG_golden_st2 << factor_5t  << " \n ";
							DTFAG_golden_st2 << factor_6t << " \n ";
							DTFAG_golden_st2 << factor_7t << " \n ";
							DTFAG_golden_st2 << factor_8t << " \n ";
							DTFAG_golden_st2 << factor_9t  << " \n ";
							DTFAG_golden_st2 << factor_10t << " \n ";
							DTFAG_golden_st2 << factor_11t << " \n ";
							DTFAG_golden_st2 << factor_12t << " \n ";
							DTFAG_golden_st2 << factor_13t  << " \n ";
							DTFAG_golden_st2 << factor_14t << " \n ";
							DTFAG_golden_st2 << factor_15t << " \n ";
							break;
						default:
							break;
					}
					//--------------------------------------------------					
				}
			}
		//data relocation
		 if(s < Stage-1){
		  if(bn1_bc_tmp > bn0_bc_tmp){
			 data_tmp_1  = A_B0R1[bn0_ma_reg1];
			 data_tmp_2  = A_B0R2[bn0_ma_reg1];
			 data_tmp_3  = A_B0R3[bn0_ma_reg1];
			 data_tmp_4  = A_B0R4[bn0_ma_reg1];
			 data_tmp_5  = A_B0R5[bn0_ma_reg1];
			 data_tmp_6  = A_B0R6[bn0_ma_reg1];
			 data_tmp_7  = A_B0R7[bn0_ma_reg1];
			 data_tmp_8  = A_B0R8[bn0_ma_reg1];
			 data_tmp_9  = A_B0R9[bn0_ma_reg1];
			 data_tmp_10 = A_B0R10[bn0_ma_reg1];
			 data_tmp_11 = A_B0R11[bn0_ma_reg1];
			 data_tmp_12 = A_B0R12[bn0_ma_reg1];
			 data_tmp_13 = A_B0R13[bn0_ma_reg1];
			 data_tmp_14 = A_B0R14[bn0_ma_reg1];
			 data_tmp_15 = A_B0R15[bn0_ma_reg1];
			 A_B0R1[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg1];
			 A_B0R2[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg2];
			 A_B0R3[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg2];
			 A_B0R4[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg3];
			 A_B0R5[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg3];
			 A_B0R6[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg1] = A_B0R0[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg1] = A_B1R0[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg1] = A_B0R0[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg1] = A_B1R0[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg1] = A_B1R0[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg1] = A_B0R0[bn0_ma_reg8];
			 A_B1R0[bn1_ma_reg1]  = data_tmp_1; 
			 A_B1R0[bn1_ma_reg2]  = data_tmp_2; 
			 A_B0R0[bn0_ma_reg2]  = data_tmp_3; 
			 A_B1R0[bn1_ma_reg3]  = data_tmp_4; 
			 A_B0R0[bn0_ma_reg3]  = data_tmp_5; 
			 A_B0R0[bn0_ma_reg4]  = data_tmp_6; 
			 A_B1R0[bn1_ma_reg4]  = data_tmp_7; 
			 A_B1R0[bn1_ma_reg5]  = data_tmp_8; 
			 A_B0R0[bn0_ma_reg5]  = data_tmp_9; 
			 A_B0R0[bn0_ma_reg6]  = data_tmp_10;
			 A_B1R0[bn1_ma_reg6]  = data_tmp_11;
			 A_B0R0[bn0_ma_reg7]  = data_tmp_12;
			 A_B1R0[bn1_ma_reg7]  = data_tmp_13;
			 A_B1R0[bn1_ma_reg8]  = data_tmp_14;
			 A_B0R0[bn0_ma_reg8]  = data_tmp_15;
			 /*********************************************/
			 data_tmp_1  = A_B1R2[bn1_ma_reg1];
			 data_tmp_2  = A_B1R3[bn1_ma_reg1];
			 data_tmp_3  = A_B1R4[bn1_ma_reg1];
			 data_tmp_4  = A_B1R5[bn1_ma_reg1];
			 data_tmp_5  = A_B1R6[bn1_ma_reg1];
			 data_tmp_6  = A_B1R7[bn1_ma_reg1];
			 data_tmp_7  = A_B1R8[bn1_ma_reg1];
			 data_tmp_8  = A_B1R9[bn1_ma_reg1];
			 data_tmp_9  = A_B1R10[bn1_ma_reg1];
			 data_tmp_10 = A_B1R11[bn1_ma_reg1];
			 data_tmp_11 = A_B1R12[bn1_ma_reg1];
			 data_tmp_12 = A_B1R13[bn1_ma_reg1];
			 data_tmp_13 = A_B1R14[bn1_ma_reg1];
			 data_tmp_14 = A_B1R15[bn1_ma_reg1];
             A_B1R2[bn1_ma_reg1]  =	 A_B1R1[bn1_ma_reg2];
             A_B1R3[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg2];
             A_B1R4[bn1_ma_reg1]  =  A_B1R1[bn1_ma_reg3];
             A_B1R5[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg3];
             A_B1R6[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg4];
             A_B1R7[bn1_ma_reg1]  =  A_B1R1[bn1_ma_reg4];
             A_B1R8[bn1_ma_reg1]  =  A_B1R1[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg1] =  A_B0R1[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg1] =  A_B1R1[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg1] =  A_B0R1[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg1] =  A_B1R1[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg1] =  A_B1R1[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg1] =  A_B0R1[bn0_ma_reg8];
			 A_B1R1[bn1_ma_reg2]  =  data_tmp_1; 
			 A_B0R1[bn0_ma_reg2]  =  data_tmp_2; 
			 A_B1R1[bn1_ma_reg3]  =  data_tmp_3; 
			 A_B0R1[bn0_ma_reg3]  =  data_tmp_4; 
			 A_B0R1[bn0_ma_reg4]  =  data_tmp_5; 
			 A_B1R1[bn1_ma_reg4]  =  data_tmp_6; 
			 A_B1R1[bn1_ma_reg5]  =  data_tmp_7; 
			 A_B0R1[bn0_ma_reg5]  =  data_tmp_8; 
			 A_B0R1[bn0_ma_reg6]  =  data_tmp_9; 
			 A_B1R1[bn1_ma_reg6]  =  data_tmp_10;
			 A_B0R1[bn0_ma_reg7]  =  data_tmp_11;
			 A_B1R1[bn1_ma_reg7]  =  data_tmp_12;
			 A_B1R1[bn1_ma_reg8]  =  data_tmp_13;
			 A_B0R1[bn0_ma_reg8]  =  data_tmp_14;
			/************************************************************/ 
			 data_tmp_1  =   A_B1R3[bn1_ma_reg2];
			 data_tmp_2  =   A_B1R4[bn1_ma_reg2];
			 data_tmp_3  =   A_B1R5[bn1_ma_reg2];
			 data_tmp_4  =   A_B1R6[bn1_ma_reg2];
			 data_tmp_5  =   A_B1R7[bn1_ma_reg2];
			 data_tmp_6  =   A_B1R8[bn1_ma_reg2];
			 data_tmp_7  =   A_B1R9[bn1_ma_reg2];
			 data_tmp_8  =   A_B1R10[bn1_ma_reg2];
			 data_tmp_9  =   A_B1R11[bn1_ma_reg2];
			 data_tmp_10 =   A_B1R12[bn1_ma_reg2];
			 data_tmp_11 =   A_B1R13[bn1_ma_reg2];
			 data_tmp_12 =   A_B1R14[bn1_ma_reg2];
			 data_tmp_13 =   A_B1R15[bn1_ma_reg2];
			 A_B1R3[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg2];
			 A_B1R4[bn1_ma_reg2]  =  A_B1R2[bn1_ma_reg3];
			 A_B1R5[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg3];
			 A_B1R6[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg4];
			 A_B1R7[bn1_ma_reg2]  =  A_B1R2[bn1_ma_reg4];
			 A_B1R8[bn1_ma_reg2]  =  A_B1R2[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg2] =  A_B0R2[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg2] =  A_B1R2[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg2] =  A_B0R2[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg2] =  A_B1R2[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg2] =  A_B1R2[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg2] =  A_B0R2[bn0_ma_reg8];
			 A_B0R2[bn0_ma_reg2]  =  data_tmp_1; 
			 A_B1R2[bn1_ma_reg3]  =  data_tmp_2; 
			 A_B0R2[bn0_ma_reg3]  =  data_tmp_3; 
			 A_B0R2[bn0_ma_reg4]  =  data_tmp_4; 
			 A_B1R2[bn1_ma_reg4]  =  data_tmp_5; 
			 A_B1R2[bn1_ma_reg5]  =  data_tmp_6; 
			 A_B0R2[bn0_ma_reg5]  =  data_tmp_7; 
			 A_B0R2[bn0_ma_reg6]  =  data_tmp_8; 
			 A_B1R2[bn1_ma_reg6]  =  data_tmp_9; 
			 A_B0R2[bn0_ma_reg7]  =  data_tmp_10;
			 A_B1R2[bn1_ma_reg7]  =  data_tmp_11;
			 A_B1R2[bn1_ma_reg8]  =  data_tmp_12;
			 A_B0R2[bn0_ma_reg8]  =  data_tmp_13;
			 //******************************************************
			 data_tmp_1  = A_B0R4[bn0_ma_reg2];
			 data_tmp_2  = A_B0R5[bn0_ma_reg2];
			 data_tmp_3  = A_B0R6[bn0_ma_reg2];
			 data_tmp_4  = A_B0R7[bn0_ma_reg2];
			 data_tmp_5  = A_B0R8[bn0_ma_reg2];
			 data_tmp_6  = A_B0R9[bn0_ma_reg2];
			 data_tmp_7  = A_B0R10[bn0_ma_reg2];
			 data_tmp_8  = A_B0R11[bn0_ma_reg2];
			 data_tmp_9  = A_B0R12[bn0_ma_reg2];
			 data_tmp_10 = A_B0R13[bn0_ma_reg2];
			 data_tmp_11 = A_B0R14[bn0_ma_reg2];
			 data_tmp_12 = A_B0R15[bn0_ma_reg2];
			 A_B0R4[bn0_ma_reg2]  =  A_B1R3[bn1_ma_reg3];
			 A_B0R5[bn0_ma_reg2]  =  A_B0R3[bn0_ma_reg3];
			 A_B0R6[bn0_ma_reg2]  =  A_B0R3[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg2]  =  A_B1R3[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg2]  =  A_B1R3[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg2]  =  A_B0R3[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg2] =  A_B0R3[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg2] =  A_B1R3[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg2] =  A_B0R3[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg2] =  A_B1R3[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg2] =  A_B1R3[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg2] =  A_B0R3[bn0_ma_reg8];
			 A_B1R3[bn1_ma_reg3]  =  data_tmp_1;  
			 A_B0R3[bn0_ma_reg3]  =  data_tmp_2;  
			 A_B0R3[bn0_ma_reg4]  =  data_tmp_3;  
			 A_B1R3[bn1_ma_reg4]  =  data_tmp_4;  
			 A_B1R3[bn1_ma_reg5]  =  data_tmp_5;  
			 A_B0R3[bn0_ma_reg5]  =  data_tmp_6;  
			 A_B0R3[bn0_ma_reg6]  =  data_tmp_7;  
			 A_B1R3[bn1_ma_reg6]  =  data_tmp_8;  
			 A_B0R3[bn0_ma_reg7]  =  data_tmp_9;  
			 A_B1R3[bn1_ma_reg7]  =  data_tmp_10; 
			 A_B1R3[bn1_ma_reg8]  =  data_tmp_11; 
			 A_B0R3[bn0_ma_reg8]  =  data_tmp_12; 
			 //----------------------------------------------------------------------
             data_tmp_1  = A_B1R5[bn1_ma_reg3];
			 data_tmp_2  = A_B1R6[bn1_ma_reg3];
			 data_tmp_3  = A_B1R7[bn1_ma_reg3];
			 data_tmp_4  = A_B1R8[bn1_ma_reg3];
			 data_tmp_5  = A_B1R9[bn1_ma_reg3];
			 data_tmp_6  = A_B1R10[bn1_ma_reg3];
			 data_tmp_7  = A_B1R11[bn1_ma_reg3];
			 data_tmp_8  = A_B1R12[bn1_ma_reg3];
			 data_tmp_9  = A_B1R13[bn1_ma_reg3];
			 data_tmp_10 = A_B1R14[bn1_ma_reg3];
			 data_tmp_11 = A_B1R15[bn1_ma_reg3];
			 A_B1R5[bn1_ma_reg3]  =  A_B0R4[bn0_ma_reg3];
			 A_B1R6[bn1_ma_reg3]  =  A_B0R4[bn0_ma_reg4];
			 A_B1R7[bn1_ma_reg3]  =  A_B1R4[bn1_ma_reg4];
			 A_B1R8[bn1_ma_reg3]  =  A_B1R4[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg3]  =  A_B0R4[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg3] =  A_B0R4[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg3] =  A_B1R4[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg3] =  A_B0R4[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg3] =  A_B1R4[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg3] =  A_B1R4[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg3] =  A_B0R4[bn0_ma_reg8];
			 A_B0R4[bn0_ma_reg3]  =  data_tmp_1 ;
			 A_B0R4[bn0_ma_reg4]  =  data_tmp_2 ;
			 A_B1R4[bn1_ma_reg4]  =  data_tmp_3 ;
			 A_B1R4[bn1_ma_reg5]  =  data_tmp_4 ;
			 A_B0R4[bn0_ma_reg5]  =  data_tmp_5 ;
			 A_B0R4[bn0_ma_reg6]  =  data_tmp_6 ;
			 A_B1R4[bn1_ma_reg6]  =  data_tmp_7 ;
			 A_B0R4[bn0_ma_reg7]  =  data_tmp_8 ;
			 A_B1R4[bn1_ma_reg7]  =  data_tmp_9 ;
			 A_B1R4[bn1_ma_reg8]  =  data_tmp_10;
			 A_B0R4[bn0_ma_reg8]  =  data_tmp_11;
			 //************************************************************************
			 data_tmp_1  = A_B0R6[bn0_ma_reg3];
			 data_tmp_2  = A_B0R7[bn0_ma_reg3];
			 data_tmp_3  = A_B0R8[bn0_ma_reg3];
			 data_tmp_4  = A_B0R9[bn0_ma_reg3];
			 data_tmp_5  = A_B0R10[bn0_ma_reg3];
			 data_tmp_6  = A_B0R11[bn0_ma_reg3];
			 data_tmp_7  = A_B0R12[bn0_ma_reg3];
			 data_tmp_8  = A_B0R13[bn0_ma_reg3];
			 data_tmp_9  = A_B0R14[bn0_ma_reg3];
			 data_tmp_10 = A_B0R15[bn0_ma_reg3];
			 A_B0R6[bn0_ma_reg3]  =  A_B0R5[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg3]  =  A_B1R5[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg3]  =  A_B1R5[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg3]  =  A_B0R5[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg3] =  A_B0R5[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg3] =  A_B1R5[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg3] =  A_B0R5[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg3] =  A_B1R5[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg3] =  A_B1R5[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg3] =  A_B0R5[bn0_ma_reg8];
			 A_B0R5[bn0_ma_reg4]  =  data_tmp_1; 
			 A_B1R5[bn1_ma_reg4]  =  data_tmp_2; 
			 A_B1R5[bn1_ma_reg5]  =  data_tmp_3; 
			 A_B0R5[bn0_ma_reg5]  =  data_tmp_4; 
			 A_B0R5[bn0_ma_reg6]  =  data_tmp_5; 
			 A_B1R5[bn1_ma_reg6]  =  data_tmp_6; 
			 A_B0R5[bn0_ma_reg7]  =  data_tmp_7; 
			 A_B1R5[bn1_ma_reg7]  =  data_tmp_8; 
			 A_B1R5[bn1_ma_reg8]  =  data_tmp_9; 
			 A_B0R5[bn0_ma_reg8]  =  data_tmp_10;
			 //-----------------------------------------------------------------------
			 data_tmp_1  = A_B0R7[bn0_ma_reg4];
			 data_tmp_2  = A_B0R8[bn0_ma_reg4];
			 data_tmp_3  = A_B0R9[bn0_ma_reg4];
			 data_tmp_4  = A_B0R10[bn0_ma_reg4];
			 data_tmp_5  = A_B0R11[bn0_ma_reg4];
			 data_tmp_6  = A_B0R12[bn0_ma_reg4];
			 data_tmp_7  = A_B0R13[bn0_ma_reg4];
			 data_tmp_8  = A_B0R14[bn0_ma_reg4];
			 data_tmp_9  = A_B0R15[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg4]  =  A_B1R6[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg4]  =  A_B1R6[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg4]  =  A_B0R6[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg4] =  A_B0R6[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg4] =  A_B1R6[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg4] =  A_B0R6[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg4] =  A_B1R6[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg4] =  A_B1R6[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg4] =  A_B0R6[bn0_ma_reg8];
			 A_B1R6[bn1_ma_reg4]  =  data_tmp_1;
			 A_B1R6[bn1_ma_reg5]  =  data_tmp_2;
			 A_B0R6[bn0_ma_reg5]  =  data_tmp_3;
			 A_B0R6[bn0_ma_reg6]  =  data_tmp_4;
			 A_B1R6[bn1_ma_reg6]  =  data_tmp_5;
			 A_B0R6[bn0_ma_reg7]  =  data_tmp_6;
			 A_B1R6[bn1_ma_reg7]  =  data_tmp_7;
			 A_B1R6[bn1_ma_reg8]  =  data_tmp_8;
			 A_B0R6[bn0_ma_reg8]  =  data_tmp_9;
			 //----------------------------------------------------------------------
			 data_tmp_1  = A_B1R8[bn1_ma_reg4];
			 data_tmp_2  = A_B1R9[bn1_ma_reg4];
			 data_tmp_3  = A_B1R10[bn1_ma_reg4];
			 data_tmp_4  = A_B1R11[bn1_ma_reg4];
			 data_tmp_5  = A_B1R12[bn1_ma_reg4];
			 data_tmp_6  = A_B1R13[bn1_ma_reg4];
			 data_tmp_7  = A_B1R14[bn1_ma_reg4];
			 data_tmp_8  = A_B1R15[bn1_ma_reg4];
			 A_B1R8[bn1_ma_reg4]  =  A_B1R7[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg4]  =  A_B0R7[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg4] =  A_B0R7[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg4] =  A_B1R7[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg4] =  A_B0R7[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg4] =  A_B1R7[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg4] =  A_B1R7[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg4] =  A_B0R7[bn0_ma_reg8];
			 A_B1R7[bn1_ma_reg5]  =  data_tmp_1;
			 A_B0R7[bn0_ma_reg5]  =  data_tmp_2;
			 A_B0R7[bn0_ma_reg6]  =  data_tmp_3;
			 A_B1R7[bn1_ma_reg6]  =  data_tmp_4;
			 A_B0R7[bn0_ma_reg7]  =  data_tmp_5;
			 A_B1R7[bn1_ma_reg7]  =  data_tmp_6;
			 A_B1R7[bn1_ma_reg8]  =  data_tmp_7;
			 A_B0R7[bn0_ma_reg8]  =  data_tmp_8;
			 //----------------------------------------------------------------------
			 data_tmp_1 = A_B1R9[bn1_ma_reg5];
			 data_tmp_2 = A_B1R10[bn1_ma_reg5];
			 data_tmp_3 = A_B1R11[bn1_ma_reg5];
			 data_tmp_4 = A_B1R12[bn1_ma_reg5];
			 data_tmp_5 = A_B1R13[bn1_ma_reg5];
			 data_tmp_6 = A_B1R14[bn1_ma_reg5];
			 data_tmp_7 = A_B1R15[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg5]  = A_B0R8[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg5] = A_B0R8[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg5] = A_B1R8[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg5] = A_B0R8[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg5] = A_B1R8[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg5] = A_B1R8[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg5] = A_B0R8[bn0_ma_reg8];
			 A_B0R8[bn0_ma_reg5]  = data_tmp_1;
			 A_B0R8[bn0_ma_reg6]  = data_tmp_2;
			 A_B1R8[bn1_ma_reg6]  = data_tmp_3;
			 A_B0R8[bn0_ma_reg7]  = data_tmp_4;
			 A_B1R8[bn1_ma_reg7]  = data_tmp_5;
			 A_B1R8[bn1_ma_reg8]  = data_tmp_6;
			 A_B0R8[bn0_ma_reg8]  = data_tmp_7;
			 //---------------------------------------------------------------------
			 data_tmp_1  = A_B0R10[bn0_ma_reg5];
			 data_tmp_2  = A_B0R11[bn0_ma_reg5];
			 data_tmp_3  = A_B0R12[bn0_ma_reg5];
			 data_tmp_4  = A_B0R13[bn0_ma_reg5];
			 data_tmp_5  = A_B0R14[bn0_ma_reg5];
			 data_tmp_6  = A_B0R15[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg5] = A_B0R9[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg5] = A_B1R9[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg5] = A_B0R9[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg5] = A_B1R9[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg5] = A_B1R9[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg5] = A_B0R9[bn0_ma_reg8];
			 A_B0R9[bn0_ma_reg6]  = data_tmp_1;
			 A_B1R9[bn1_ma_reg6]  = data_tmp_2;
			 A_B0R9[bn0_ma_reg7]  = data_tmp_3;
			 A_B1R9[bn1_ma_reg7]  = data_tmp_4;
			 A_B1R9[bn1_ma_reg8]  = data_tmp_5;
			 A_B0R9[bn0_ma_reg8]  = data_tmp_6;
			 //--------------------------------------------------------------------
			 data_tmp_1  = A_B0R11[bn0_ma_reg6];
			 data_tmp_2  = A_B0R12[bn0_ma_reg6];
			 data_tmp_3  = A_B0R13[bn0_ma_reg6];
			 data_tmp_4  = A_B0R14[bn0_ma_reg6];
			 data_tmp_5  = A_B0R15[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg6] = A_B1R10[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg6] = A_B0R10[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg6] = A_B1R10[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg6] = A_B1R10[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg6] = A_B0R10[bn0_ma_reg8];
			 A_B1R10[bn1_ma_reg6] = data_tmp_1;
			 A_B0R10[bn0_ma_reg7] = data_tmp_2;
			 A_B1R10[bn1_ma_reg7] = data_tmp_3;
			 A_B1R10[bn1_ma_reg8] = data_tmp_4;
			 A_B0R10[bn0_ma_reg8] = data_tmp_5;
			 //--------------------------------------------------------------------
			 data_tmp_1  = A_B1R12[bn1_ma_reg6];
			 data_tmp_2  = A_B1R13[bn1_ma_reg6];
			 data_tmp_3  = A_B1R14[bn1_ma_reg6];
			 data_tmp_4  = A_B1R15[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg6] = A_B0R11[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg6] = A_B1R11[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg6] = A_B1R11[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg6] = A_B0R11[bn0_ma_reg8];
			 A_B0R11[bn0_ma_reg7] = data_tmp_1;
			 A_B1R11[bn1_ma_reg7] = data_tmp_2;
			 A_B1R11[bn1_ma_reg8] = data_tmp_3;
			 A_B0R11[bn0_ma_reg8] = data_tmp_4;
			 //--------------------------------------------------------------------
			 data_tmp_1 = A_B0R13[bn0_ma_reg7];
			 data_tmp_2 = A_B0R14[bn0_ma_reg7];
			 data_tmp_3 = A_B0R15[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg7] = A_B1R12[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg7] = A_B1R12[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg7] = A_B0R12[bn0_ma_reg8];
			 A_B1R12[bn1_ma_reg7] = data_tmp_1;
			 A_B1R12[bn1_ma_reg8] = data_tmp_2;
			 A_B0R12[bn0_ma_reg8] = data_tmp_3;
			 //--------------------------------------------------------------------
			 data_tmp_1 = A_B1R14[bn1_ma_reg7];
			 data_tmp_2 = A_B1R15[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg7] = A_B1R13[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg7] = A_B0R13[bn0_ma_reg8];
			 A_B1R13[bn1_ma_reg8] = data_tmp_1;
			 A_B0R13[bn0_ma_reg8] = data_tmp_2;
			 //--------------------------------------------------------------------
			 data_tmp_1 = A_B1R15[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg8] = A_B0R14[bn0_ma_reg8];
			 A_B0R14[bn0_ma_reg8] = data_tmp_1; 
		  }else{	 
			 data_tmp_1  = A_B1R1[bn1_ma_reg1];
			 data_tmp_2  = A_B1R2[bn1_ma_reg1];
			 data_tmp_3  = A_B1R3[bn1_ma_reg1];
			 data_tmp_4  = A_B1R4[bn1_ma_reg1];
			 data_tmp_5  = A_B1R5[bn1_ma_reg1];
			 data_tmp_6  = A_B1R6[bn1_ma_reg1];
			 data_tmp_7  = A_B1R7[bn1_ma_reg1];
			 data_tmp_8  = A_B1R8[bn1_ma_reg1];
			 data_tmp_9  = A_B1R9[bn1_ma_reg1];
			 data_tmp_10 = A_B1R10[bn1_ma_reg1];
			 data_tmp_11 = A_B1R11[bn1_ma_reg1];
			 data_tmp_12 = A_B1R12[bn1_ma_reg1];
			 data_tmp_13 = A_B1R13[bn1_ma_reg1];
			 data_tmp_14 = A_B1R14[bn1_ma_reg1];
			 data_tmp_15 = A_B1R15[bn1_ma_reg1];
             A_B1R1[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg1];
             A_B1R2[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg2];
             A_B1R3[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg2];
             A_B1R4[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg3];
             A_B1R5[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg3];
             A_B1R6[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg4];
             A_B1R7[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg4];
             A_B1R8[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg5];
             A_B1R9[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg5];
             A_B1R10[bn1_ma_reg1] =  A_B1R0[bn1_ma_reg6];
             A_B1R11[bn1_ma_reg1] =  A_B0R0[bn0_ma_reg6];
             A_B1R12[bn1_ma_reg1] =  A_B1R0[bn1_ma_reg7];
             A_B1R13[bn1_ma_reg1] =  A_B0R0[bn0_ma_reg7];
             A_B1R14[bn1_ma_reg1] =  A_B0R0[bn0_ma_reg8];
             A_B1R15[bn1_ma_reg1] =  A_B1R0[bn1_ma_reg8];
             A_B0R0[bn0_ma_reg1]  =  data_tmp_1; 
             A_B0R0[bn0_ma_reg2]  =  data_tmp_2; 
             A_B1R0[bn1_ma_reg2]  =  data_tmp_3; 
             A_B0R0[bn0_ma_reg3]  =  data_tmp_4; 
             A_B1R0[bn1_ma_reg3]  =  data_tmp_5; 
             A_B1R0[bn1_ma_reg4]  =  data_tmp_6; 
             A_B0R0[bn0_ma_reg4]  =  data_tmp_7; 
             A_B0R0[bn0_ma_reg5]  =  data_tmp_8; 
             A_B1R0[bn1_ma_reg5]  =  data_tmp_9; 
             A_B1R0[bn1_ma_reg6]  =  data_tmp_10;
             A_B0R0[bn0_ma_reg6]  =  data_tmp_11;
             A_B1R0[bn1_ma_reg7]  =  data_tmp_12;
             A_B0R0[bn0_ma_reg7]  =  data_tmp_13;
             A_B0R0[bn0_ma_reg8]  =  data_tmp_14;
             A_B1R0[bn1_ma_reg8]  =  data_tmp_15;
             //-----------------------------------------------------------
             data_tmp_1  = A_B0R2[bn0_ma_reg1];
			 data_tmp_2  = A_B0R3[bn0_ma_reg1];
			 data_tmp_3  = A_B0R4[bn0_ma_reg1];
			 data_tmp_4  = A_B0R5[bn0_ma_reg1];
			 data_tmp_5  = A_B0R6[bn0_ma_reg1];
			 data_tmp_6  = A_B0R7[bn0_ma_reg1];
			 data_tmp_7  = A_B0R8[bn0_ma_reg1];
			 data_tmp_8  = A_B0R9[bn0_ma_reg1];
			 data_tmp_9  = A_B0R10[bn0_ma_reg1];
			 data_tmp_10 = A_B0R11[bn0_ma_reg1];
			 data_tmp_11 = A_B0R12[bn0_ma_reg1];
			 data_tmp_12 = A_B0R13[bn0_ma_reg1];
			 data_tmp_13 = A_B0R14[bn0_ma_reg1];
			 data_tmp_14 = A_B0R15[bn0_ma_reg1];
			 A_B0R2[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg2];
			 A_B0R3[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg2];
			 A_B0R4[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg3];
			 A_B0R5[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg3];
			 A_B0R6[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg4];
			 A_B0R7[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg1] =  A_B1R1[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg1] =  A_B0R1[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg1] =  A_B1R1[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg1] =  A_B0R1[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg1] =  A_B0R1[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg1] =  A_B1R1[bn1_ma_reg8];
             A_B0R1[bn0_ma_reg2]  =  data_tmp_1; 
             A_B1R1[bn1_ma_reg2]  =  data_tmp_2; 
             A_B0R1[bn0_ma_reg3]  =  data_tmp_3; 
             A_B1R1[bn1_ma_reg3]  =  data_tmp_4; 
             A_B1R1[bn1_ma_reg4]  =  data_tmp_5; 
             A_B0R1[bn0_ma_reg4]  =  data_tmp_6; 
             A_B0R1[bn0_ma_reg5]  =  data_tmp_7; 
             A_B1R1[bn1_ma_reg5]  =  data_tmp_8; 
             A_B1R1[bn1_ma_reg6]  =  data_tmp_9; 
             A_B0R1[bn0_ma_reg6]  =  data_tmp_10;
             A_B1R1[bn1_ma_reg7]  =  data_tmp_11;
             A_B0R1[bn0_ma_reg7]  =  data_tmp_12;
             A_B0R1[bn0_ma_reg8]  =  data_tmp_13;
             A_B1R1[bn1_ma_reg8]  =  data_tmp_14;
             //------------------------------------------------------------
             data_tmp_1  =  A_B0R3[bn0_ma_reg2];
			 data_tmp_2  =  A_B0R4[bn0_ma_reg2];
			 data_tmp_3  =  A_B0R5[bn0_ma_reg2];
			 data_tmp_4  =  A_B0R6[bn0_ma_reg2];
			 data_tmp_5  =  A_B0R7[bn0_ma_reg2];
			 data_tmp_6  =  A_B0R8[bn0_ma_reg2];
			 data_tmp_7  =  A_B0R9[bn0_ma_reg2];
			 data_tmp_8  =  A_B0R10[bn0_ma_reg2];
			 data_tmp_9  =  A_B0R11[bn0_ma_reg2];
			 data_tmp_10 =  A_B0R12[bn0_ma_reg2];
			 data_tmp_11 =  A_B0R13[bn0_ma_reg2];
			 data_tmp_12 =  A_B0R14[bn0_ma_reg2];
			 data_tmp_13 =  A_B0R15[bn0_ma_reg2];
			 A_B0R3[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg2];
			 A_B0R4[bn0_ma_reg2]  =  A_B0R2[bn0_ma_reg3];
			 A_B0R5[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg3];
			 A_B0R6[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg4];
			 A_B0R7[bn0_ma_reg2]  =  A_B0R2[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg2]  =  A_B0R2[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg2] =  A_B1R2[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg2] =  A_B0R2[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg2] =  A_B1R2[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg2] =  A_B0R2[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg2] =  A_B0R2[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg2] =  A_B1R2[bn1_ma_reg8];
			 A_B1R2[bn1_ma_reg2]  =  data_tmp_1; 
			 A_B0R2[bn0_ma_reg3]  =  data_tmp_2; 
			 A_B1R2[bn1_ma_reg3]  =  data_tmp_3; 
			 A_B1R2[bn1_ma_reg4]  =  data_tmp_4; 
			 A_B0R2[bn0_ma_reg4]  =  data_tmp_5; 
			 A_B0R2[bn0_ma_reg5]  =  data_tmp_6; 
			 A_B1R2[bn1_ma_reg5]  =  data_tmp_7; 
			 A_B1R2[bn1_ma_reg6]  =  data_tmp_8; 
			 A_B0R2[bn0_ma_reg6]  =  data_tmp_9; 
			 A_B1R2[bn1_ma_reg7]  =  data_tmp_10;
			 A_B0R2[bn0_ma_reg7]  =  data_tmp_11;
			 A_B0R2[bn0_ma_reg8]  =  data_tmp_12;
			 A_B1R2[bn1_ma_reg8]  =  data_tmp_13;
			 //-----------------------------------------------------------                        
			 data_tmp_1  = A_B1R4[bn1_ma_reg2];
			 data_tmp_2  = A_B1R5[bn1_ma_reg2];
			 data_tmp_3  = A_B1R6[bn1_ma_reg2];
			 data_tmp_4  = A_B1R7[bn1_ma_reg2];
			 data_tmp_5  = A_B1R8[bn1_ma_reg2];
			 data_tmp_6  = A_B1R9[bn1_ma_reg2];
			 data_tmp_7  = A_B1R10[bn1_ma_reg2];
			 data_tmp_8  = A_B1R11[bn1_ma_reg2];
			 data_tmp_9  = A_B1R12[bn1_ma_reg2];
			 data_tmp_10 = A_B1R13[bn1_ma_reg2];
			 data_tmp_11 = A_B1R14[bn1_ma_reg2];
			 data_tmp_12 = A_B1R15[bn1_ma_reg2];
			 A_B1R4[bn1_ma_reg2]  = A_B0R3[bn0_ma_reg3];
			 A_B1R5[bn1_ma_reg2]  = A_B1R3[bn1_ma_reg3];
			 A_B1R6[bn1_ma_reg2]  = A_B1R3[bn1_ma_reg4];
			 A_B1R7[bn1_ma_reg2]  = A_B0R3[bn0_ma_reg4];
			 A_B1R8[bn1_ma_reg2]  = A_B0R3[bn0_ma_reg5];
			 A_B1R9[bn1_ma_reg2]  = A_B1R3[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg2] = A_B1R3[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg2] = A_B0R3[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg2] = A_B1R3[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg2] = A_B0R3[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg2] = A_B0R3[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg2] = A_B1R3[bn1_ma_reg8];
			 A_B0R3[bn0_ma_reg3]  = data_tmp_1; 
			 A_B1R3[bn1_ma_reg3]  = data_tmp_2; 
			 A_B1R3[bn1_ma_reg4]  = data_tmp_3; 
			 A_B0R3[bn0_ma_reg4]  = data_tmp_4; 
			 A_B0R3[bn0_ma_reg5]  = data_tmp_5; 
			 A_B1R3[bn1_ma_reg5]  = data_tmp_6; 
			 A_B1R3[bn1_ma_reg6]  = data_tmp_7; 
			 A_B0R3[bn0_ma_reg6]  = data_tmp_8; 
			 A_B1R3[bn1_ma_reg7]  = data_tmp_9; 
			 A_B0R3[bn0_ma_reg7]  = data_tmp_10;
			 A_B0R3[bn0_ma_reg8]  = data_tmp_11;
			 A_B1R3[bn1_ma_reg8]  = data_tmp_12;
			 //------------------------------------------------------------
			 data_tmp_1  =  A_B0R5[bn0_ma_reg3];
			 data_tmp_2  =  A_B0R6[bn0_ma_reg3];
			 data_tmp_3  =  A_B0R7[bn0_ma_reg3];
			 data_tmp_4  =  A_B0R8[bn0_ma_reg3];
			 data_tmp_5  =  A_B0R9[bn0_ma_reg3];
			 data_tmp_6  =  A_B0R10[bn0_ma_reg3];
			 data_tmp_7  =  A_B0R11[bn0_ma_reg3];
			 data_tmp_8  =  A_B0R12[bn0_ma_reg3];
			 data_tmp_9  =  A_B0R13[bn0_ma_reg3];
			 data_tmp_10 =  A_B0R14[bn0_ma_reg3];
			 data_tmp_11 =  A_B0R15[bn0_ma_reg3];
			 A_B0R5[bn0_ma_reg3]  =  A_B1R4[bn1_ma_reg3];
			 A_B0R6[bn0_ma_reg3]  =  A_B1R4[bn1_ma_reg4];
			 A_B0R7[bn0_ma_reg3]  =  A_B0R4[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg3]  =  A_B0R4[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg3]  =  A_B1R4[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg3] =  A_B1R4[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg3] =  A_B0R4[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg3] =  A_B1R4[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg3] =  A_B0R4[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg3] =  A_B0R4[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg3] =  A_B1R4[bn1_ma_reg8];
			 A_B1R4[bn1_ma_reg3]  =  data_tmp_1; 
			 A_B1R4[bn1_ma_reg4]  =  data_tmp_2; 
			 A_B0R4[bn0_ma_reg4]  =  data_tmp_3; 
			 A_B0R4[bn0_ma_reg5]  =  data_tmp_4; 
			 A_B1R4[bn1_ma_reg5]  =  data_tmp_5; 
			 A_B1R4[bn1_ma_reg6]  =  data_tmp_6; 
			 A_B0R4[bn0_ma_reg6]  =  data_tmp_7; 
			 A_B1R4[bn1_ma_reg7]  =  data_tmp_8; 
			 A_B0R4[bn0_ma_reg7]  =  data_tmp_9; 
			 A_B0R4[bn0_ma_reg8]  =  data_tmp_10;
			 A_B1R4[bn1_ma_reg8]  =  data_tmp_11;
			 //-------------------------------------------------------------
			 data_tmp_1  =  A_B1R6[bn1_ma_reg3];
			 data_tmp_2  =  A_B1R7[bn1_ma_reg3];
			 data_tmp_3  =  A_B1R8[bn1_ma_reg3];
			 data_tmp_4  =  A_B1R9[bn1_ma_reg3];
			 data_tmp_5  =  A_B1R10[bn1_ma_reg3];
			 data_tmp_6  =  A_B1R11[bn1_ma_reg3];
			 data_tmp_7  =  A_B1R12[bn1_ma_reg3];
			 data_tmp_8  =  A_B1R13[bn1_ma_reg3];
			 data_tmp_9  =  A_B1R14[bn1_ma_reg3];
			 data_tmp_10 =  A_B1R15[bn1_ma_reg3];
			 A_B1R6[bn1_ma_reg3]  =  A_B1R5[bn1_ma_reg4];
			 A_B1R7[bn1_ma_reg3]  =  A_B0R5[bn0_ma_reg4];
			 A_B1R8[bn1_ma_reg3]  =  A_B0R5[bn0_ma_reg5];
			 A_B1R9[bn1_ma_reg3]  =  A_B1R5[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg3] =  A_B1R5[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg3] =  A_B0R5[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg3] =  A_B1R5[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg3] =  A_B0R5[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg3] =  A_B0R5[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg3] =  A_B1R5[bn1_ma_reg8];
			 A_B1R5[bn1_ma_reg4]  =  data_tmp_1; 
			 A_B0R5[bn0_ma_reg4]  =  data_tmp_2; 
			 A_B0R5[bn0_ma_reg5]  =  data_tmp_3; 
			 A_B1R5[bn1_ma_reg5]  =  data_tmp_4; 
			 A_B1R5[bn1_ma_reg6]  =  data_tmp_5; 
			 A_B0R5[bn0_ma_reg6]  =  data_tmp_6; 
			 A_B1R5[bn1_ma_reg7]  =  data_tmp_7; 
			 A_B0R5[bn0_ma_reg7]  =  data_tmp_8; 
			 A_B0R5[bn0_ma_reg8]  =  data_tmp_9; 
			 A_B1R5[bn1_ma_reg8]  =  data_tmp_10;
			 //----------------------------------------------------------------
			 data_tmp_1  =  A_B1R7[bn1_ma_reg4];
			 data_tmp_2  =  A_B1R8[bn1_ma_reg4];
			 data_tmp_3  =  A_B1R9[bn1_ma_reg4];
			 data_tmp_4  =  A_B1R10[bn1_ma_reg4];
			 data_tmp_5  =  A_B1R11[bn1_ma_reg4];
			 data_tmp_6  =  A_B1R12[bn1_ma_reg4];
			 data_tmp_7  =  A_B1R13[bn1_ma_reg4];
			 data_tmp_8  =  A_B1R14[bn1_ma_reg4];
			 data_tmp_9  =  A_B1R15[bn1_ma_reg4];
			 A_B1R7[bn1_ma_reg4]  =  A_B0R6[bn0_ma_reg4];
			 A_B1R8[bn1_ma_reg4]  =  A_B0R6[bn0_ma_reg5];
			 A_B1R9[bn1_ma_reg4]  =  A_B1R6[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg4] =  A_B1R6[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg4] =  A_B0R6[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg4] =  A_B1R6[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg4] =  A_B0R6[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg4] =  A_B0R6[bn0_ma_reg8];
             A_B1R15[bn1_ma_reg4] =  A_B1R6[bn1_ma_reg8];
             A_B0R6[bn0_ma_reg4]  =  data_tmp_1;
             A_B0R6[bn0_ma_reg5]  =  data_tmp_2;
             A_B1R6[bn1_ma_reg5]  =  data_tmp_3;
             A_B1R6[bn1_ma_reg6]  =  data_tmp_4;
             A_B0R6[bn0_ma_reg6]  =  data_tmp_5;
             A_B1R6[bn1_ma_reg7]  =  data_tmp_6;
             A_B0R6[bn0_ma_reg7]  =  data_tmp_7;
             A_B0R6[bn0_ma_reg8]  =  data_tmp_8;
             A_B1R6[bn1_ma_reg8]  =  data_tmp_9;
             //------------------------------------------------------------------
			 data_tmp_1  = A_B0R8[bn0_ma_reg4];
			 data_tmp_2  = A_B0R9[bn0_ma_reg4];
			 data_tmp_3  = A_B0R10[bn0_ma_reg4];
			 data_tmp_4  = A_B0R11[bn0_ma_reg4];
			 data_tmp_5  = A_B0R12[bn0_ma_reg4];
			 data_tmp_6  = A_B0R13[bn0_ma_reg4];
			 data_tmp_7  = A_B0R14[bn0_ma_reg4];
			 data_tmp_8  = A_B0R15[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg4]  =  A_B0R7[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg4]  =  A_B1R7[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg4] =  A_B1R7[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg4] =  A_B0R7[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg4] =  A_B1R7[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg4] =  A_B0R7[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg4] =  A_B0R7[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg4] =  A_B1R7[bn1_ma_reg8];
			 A_B0R7[bn0_ma_reg5]  =  data_tmp_1;
			 A_B1R7[bn1_ma_reg5]  =  data_tmp_2;
			 A_B1R7[bn1_ma_reg6]  =  data_tmp_3;
			 A_B0R7[bn0_ma_reg6]  =  data_tmp_4;
			 A_B1R7[bn1_ma_reg7]  =  data_tmp_5;
			 A_B0R7[bn0_ma_reg7]  =  data_tmp_6;
			 A_B0R7[bn0_ma_reg8]  =  data_tmp_7;
			 A_B1R7[bn1_ma_reg8]  =  data_tmp_8;
			 //----------------------------------------------------------------
			 data_tmp_1  = A_B0R9[bn0_ma_reg5];
			 data_tmp_2  = A_B0R10[bn0_ma_reg5];
			 data_tmp_3  = A_B0R11[bn0_ma_reg5];
			 data_tmp_4  = A_B0R12[bn0_ma_reg5];
			 data_tmp_5  = A_B0R13[bn0_ma_reg5];
			 data_tmp_6  = A_B0R14[bn0_ma_reg5];
			 data_tmp_7  = A_B0R15[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg5]  =  A_B1R8[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg8];
			 A_B1R8[bn1_ma_reg5]  =  data_tmp_1;
			 A_B1R8[bn1_ma_reg6]  =  data_tmp_2;
			 A_B0R8[bn0_ma_reg6]  =  data_tmp_3;
			 A_B1R8[bn1_ma_reg7]  =  data_tmp_4;
			 A_B0R8[bn0_ma_reg7]  =  data_tmp_5;
			 A_B0R8[bn0_ma_reg8]  =  data_tmp_6;
			 A_B1R8[bn1_ma_reg8]  =  data_tmp_7;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B1R10[bn1_ma_reg5];
			 data_tmp_2 = A_B1R11[bn1_ma_reg5];
			 data_tmp_3 = A_B1R12[bn1_ma_reg5];
			 data_tmp_4 = A_B1R13[bn1_ma_reg5];
			 data_tmp_5 = A_B1R14[bn1_ma_reg5];
			 data_tmp_6 = A_B1R15[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg5] = A_B1R9[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg5] = A_B0R9[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg5] = A_B1R9[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg5] = A_B0R9[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg5] = A_B0R9[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg5] = A_B1R9[bn1_ma_reg8];
			 A_B1R9[bn1_ma_reg6]  = data_tmp_1;
			 A_B0R9[bn0_ma_reg6]  = data_tmp_2;
			 A_B1R9[bn1_ma_reg7]  = data_tmp_3;
			 A_B0R9[bn0_ma_reg7]  = data_tmp_4;
			 A_B0R9[bn0_ma_reg8]  = data_tmp_5;
			 A_B1R9[bn1_ma_reg8]  = data_tmp_6;
			 //----------------------------------------------------------------
			 data_tmp_1 =  A_B1R11[bn1_ma_reg6];
			 data_tmp_2 =  A_B1R12[bn1_ma_reg6];
			 data_tmp_3 =  A_B1R13[bn1_ma_reg6];
			 data_tmp_4 =  A_B1R14[bn1_ma_reg6];
			 data_tmp_5 =  A_B1R15[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg6]  = A_B0R10[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg6]  = A_B1R10[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg6]  = A_B0R10[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg6]  = A_B0R10[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg6]  = A_B1R10[bn1_ma_reg8];
			 A_B0R10[bn0_ma_reg6]  = data_tmp_1;
			 A_B1R10[bn1_ma_reg7]  = data_tmp_2;
			 A_B0R10[bn0_ma_reg7]  = data_tmp_3;
			 A_B0R10[bn0_ma_reg8]  = data_tmp_4;
			 A_B1R10[bn1_ma_reg8]  = data_tmp_5;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B0R12[bn0_ma_reg6];
			 data_tmp_2 = A_B0R13[bn0_ma_reg6];
			 data_tmp_3 = A_B0R14[bn0_ma_reg6];
			 data_tmp_4 = A_B0R15[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg6] = A_B1R11[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg6] = A_B0R11[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg6] = A_B0R11[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg6] = A_B1R11[bn1_ma_reg8];
			 A_B1R11[bn1_ma_reg7] = data_tmp_1;
			 A_B0R11[bn0_ma_reg7] = data_tmp_2;
			 A_B0R11[bn0_ma_reg8] = data_tmp_3;
			 A_B1R11[bn1_ma_reg8] = data_tmp_4;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B1R13[bn1_ma_reg7];
			 data_tmp_2 = A_B1R14[bn1_ma_reg7];
			 data_tmp_3 = A_B1R15[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg7] = A_B0R12[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg7] = A_B0R12[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg7] = A_B1R12[bn1_ma_reg8];
			 A_B0R12[bn0_ma_reg7] = data_tmp_1;
			 A_B0R12[bn0_ma_reg8] = data_tmp_2;
			 A_B1R12[bn1_ma_reg8] = data_tmp_3;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B0R14[bn0_ma_reg7];
			 data_tmp_2 = A_B0R15[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg7] = A_B0R13[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg7] = A_B1R13[bn1_ma_reg8];
			 A_B0R13[bn0_ma_reg8] = data_tmp_1;
			 A_B1R13[bn1_ma_reg8] = data_tmp_2;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B0R15[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg8] = A_B1R14[bn1_ma_reg8];
			 A_B1R14[bn1_ma_reg8] = data_tmp_1;
			 //-----------------------------------------------------------------
		  }
		 }
		 else{
		  //radix-16 final stage data relocation
		  //use two radix-8 dc
		    if(bn1_bc_tmp > bn0_bc_tmp){
			   data_tmp_1  = A_B0R1[bn0_ma_reg1];
			   data_tmp_2  = A_B0R2[bn0_ma_reg1];
			   data_tmp_3  = A_B0R3[bn0_ma_reg1];
			   data_tmp_4  = A_B0R4[bn0_ma_reg1];
			   data_tmp_5  = A_B0R5[bn0_ma_reg1];
			   data_tmp_6  = A_B0R6[bn0_ma_reg1];
			   data_tmp_7  = A_B0R7[bn0_ma_reg1];
			   data_tmp_9  = A_B0R9[bn0_ma_reg1];
			   data_tmp_10 = A_B0R10[bn0_ma_reg1];
			   data_tmp_11 = A_B0R11[bn0_ma_reg1];
			   data_tmp_12 = A_B0R12[bn0_ma_reg1];
			   data_tmp_13 = A_B0R13[bn0_ma_reg1];
			   data_tmp_14 = A_B0R14[bn0_ma_reg1];
			   data_tmp_15 = A_B0R15[bn0_ma_reg1];
			   A_B0R1[bn0_ma_reg1]  =  A_B1R0[bn1_ma_reg1];
			   A_B0R2[bn0_ma_reg1]  =  A_B1R0[bn1_ma_reg2];
			   A_B0R3[bn0_ma_reg1]  =  A_B0R0[bn0_ma_reg2];
			   A_B0R4[bn0_ma_reg1]  =  A_B1R0[bn1_ma_reg3];
			   A_B0R5[bn0_ma_reg1]  =  A_B0R0[bn0_ma_reg3];
			   A_B0R6[bn0_ma_reg1]  =  A_B0R0[bn0_ma_reg4];
			   A_B0R7[bn0_ma_reg1]  =  A_B1R0[bn1_ma_reg4];
			   A_B0R9[bn0_ma_reg1]  =  A_B1R8[bn1_ma_reg1];
			   A_B0R10[bn0_ma_reg1] =  A_B1R8[bn1_ma_reg2];
			   A_B0R11[bn0_ma_reg1] =  A_B0R8[bn0_ma_reg2];
			   A_B0R12[bn0_ma_reg1] =  A_B1R8[bn1_ma_reg3];
			   A_B0R13[bn0_ma_reg1] =  A_B0R8[bn0_ma_reg3];
			   A_B0R14[bn0_ma_reg1] =  A_B0R8[bn0_ma_reg4];
			   A_B0R15[bn0_ma_reg1] =  A_B1R8[bn1_ma_reg4];
			   A_B1R0[bn1_ma_reg1]  =  data_tmp_1;
			   A_B1R0[bn1_ma_reg2]  =  data_tmp_2;
			   A_B0R0[bn0_ma_reg2]  =  data_tmp_3;
			   A_B1R0[bn1_ma_reg3]  =  data_tmp_4;
			   A_B0R0[bn0_ma_reg3]  =  data_tmp_5;
			   A_B0R0[bn0_ma_reg4]  =  data_tmp_6;
			   A_B1R0[bn1_ma_reg4]  =  data_tmp_7;
			   A_B1R8[bn1_ma_reg1]  =  data_tmp_9;
			   A_B1R8[bn1_ma_reg2]  =  data_tmp_10; 
			   A_B0R8[bn0_ma_reg2]  =  data_tmp_11; 
			   A_B1R8[bn1_ma_reg3]  =  data_tmp_12; 
			   A_B0R8[bn0_ma_reg3]  =  data_tmp_13; 
			   A_B0R8[bn0_ma_reg4]  =  data_tmp_14; 
			   A_B1R8[bn1_ma_reg4]  =  data_tmp_15; 
			   //**********************************
			   data_tmp_1  = A_B1R2[bn1_ma_reg1];
			   data_tmp_2  = A_B1R3[bn1_ma_reg1];
			   data_tmp_3  = A_B1R4[bn1_ma_reg1];
			   data_tmp_4  = A_B1R5[bn1_ma_reg1];
			   data_tmp_5  = A_B1R6[bn1_ma_reg1];
			   data_tmp_6  = A_B1R7[bn1_ma_reg1];
			   data_tmp_9  = A_B1R10[bn1_ma_reg1];
			   data_tmp_10 = A_B1R11[bn1_ma_reg1];
			   data_tmp_11 = A_B1R12[bn1_ma_reg1];
			   data_tmp_12 = A_B1R13[bn1_ma_reg1];
			   data_tmp_13 = A_B1R14[bn1_ma_reg1];
			   data_tmp_14 = A_B1R15[bn1_ma_reg1];
			   A_B1R2[bn1_ma_reg1]   =  A_B1R1[bn1_ma_reg2];
			   A_B1R3[bn1_ma_reg1]   =  A_B0R1[bn0_ma_reg2];
			   A_B1R4[bn1_ma_reg1]   =  A_B1R1[bn1_ma_reg3];
			   A_B1R5[bn1_ma_reg1]   =  A_B0R1[bn0_ma_reg3];
			   A_B1R6[bn1_ma_reg1]   =  A_B0R1[bn0_ma_reg4];
			   A_B1R7[bn1_ma_reg1]   =  A_B1R1[bn1_ma_reg4];
			   A_B1R10[bn1_ma_reg1]  =  A_B1R9[bn1_ma_reg2];
			   A_B1R11[bn1_ma_reg1]  =  A_B0R9[bn0_ma_reg2];
			   A_B1R12[bn1_ma_reg1]  =  A_B1R9[bn1_ma_reg3];
			   A_B1R13[bn1_ma_reg1]  =  A_B0R9[bn0_ma_reg3];
			   A_B1R14[bn1_ma_reg1]  =  A_B0R9[bn0_ma_reg4];
			   A_B1R15[bn1_ma_reg1]  =  A_B1R9[bn1_ma_reg4];
			   A_B1R1[bn1_ma_reg2] =  data_tmp_1;
			   A_B0R1[bn0_ma_reg2] =  data_tmp_2;
			   A_B1R1[bn1_ma_reg3] =  data_tmp_3;
			   A_B0R1[bn0_ma_reg3] =  data_tmp_4;
			   A_B0R1[bn0_ma_reg4] =  data_tmp_5;
			   A_B1R1[bn1_ma_reg4] =  data_tmp_6;
			   A_B1R9[bn1_ma_reg2] =  data_tmp_9;
			   A_B0R9[bn0_ma_reg2] =  data_tmp_10; 
			   A_B1R9[bn1_ma_reg3] =  data_tmp_11; 
			   A_B0R9[bn0_ma_reg3] =  data_tmp_12; 
			   A_B0R9[bn0_ma_reg4] =  data_tmp_13; 
			   A_B1R9[bn1_ma_reg4] =  data_tmp_14; 
			   //*************************************
			   data_tmp_1  =  A_B1R3[bn1_ma_reg2];
			   data_tmp_2  =  A_B1R4[bn1_ma_reg2];
			   data_tmp_3  =  A_B1R5[bn1_ma_reg2];
			   data_tmp_4  =  A_B1R6[bn1_ma_reg2];
			   data_tmp_5  =  A_B1R7[bn1_ma_reg2];
			   data_tmp_9  =  A_B1R11[bn1_ma_reg2];
			   data_tmp_10 =  A_B1R12[bn1_ma_reg2];
			   data_tmp_11 =  A_B1R13[bn1_ma_reg2];
			   data_tmp_12 =  A_B1R14[bn1_ma_reg2];
			   data_tmp_13 =  A_B1R15[bn1_ma_reg2];
			   A_B1R3[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg2];
			   A_B1R4[bn1_ma_reg2]  =  A_B1R2[bn1_ma_reg3];
			   A_B1R5[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg3];
			   A_B1R6[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg4];
			   A_B1R7[bn1_ma_reg2]  =  A_B1R2[bn1_ma_reg4];
			   A_B1R11[bn1_ma_reg2] =  A_B0R10[bn0_ma_reg2];
			   A_B1R12[bn1_ma_reg2] =  A_B1R10[bn1_ma_reg3];
			   A_B1R13[bn1_ma_reg2] =  A_B0R10[bn0_ma_reg3];
			   A_B1R14[bn1_ma_reg2] =  A_B0R10[bn0_ma_reg4];
			   A_B1R15[bn1_ma_reg2] =  A_B1R10[bn1_ma_reg4];
			   A_B0R2[bn0_ma_reg2]  =  data_tmp_1;
			   A_B1R2[bn1_ma_reg3]  =  data_tmp_2;
			   A_B0R2[bn0_ma_reg3]  =  data_tmp_3;
			   A_B0R2[bn0_ma_reg4]  =  data_tmp_4;
			   A_B1R2[bn1_ma_reg4]  =  data_tmp_5;
			   A_B0R10[bn0_ma_reg2] =  data_tmp_9;
			   A_B1R10[bn1_ma_reg3] =  data_tmp_10;
			   A_B0R10[bn0_ma_reg3] =  data_tmp_11;
			   A_B0R10[bn0_ma_reg4] =  data_tmp_12;
			   A_B1R10[bn1_ma_reg4] =  data_tmp_13;
			   //************************************
			   data_tmp_1  =  A_B0R4[bn0_ma_reg2];
			   data_tmp_2  =  A_B0R5[bn0_ma_reg2];
			   data_tmp_3  =  A_B0R6[bn0_ma_reg2];
			   data_tmp_4  =  A_B0R7[bn0_ma_reg2];
			   data_tmp_9  =  A_B0R12[bn0_ma_reg2];
			   data_tmp_10 =  A_B0R13[bn0_ma_reg2];
 			   data_tmp_11 =  A_B0R14[bn0_ma_reg2];
			   data_tmp_12 =  A_B0R15[bn0_ma_reg2];
			   A_B0R4[bn0_ma_reg2]  =  A_B1R3[bn1_ma_reg3];
			   A_B0R5[bn0_ma_reg2]  =  A_B0R3[bn0_ma_reg3];
			   A_B0R6[bn0_ma_reg2]  =  A_B0R3[bn0_ma_reg4];
			   A_B0R7[bn0_ma_reg2]  =  A_B1R3[bn1_ma_reg4];
			   A_B0R12[bn0_ma_reg2] =  A_B1R11[bn1_ma_reg3];
			   A_B0R13[bn0_ma_reg2] =  A_B0R11[bn0_ma_reg3];
			   A_B0R14[bn0_ma_reg2] =  A_B0R11[bn0_ma_reg4];
			   A_B0R15[bn0_ma_reg2] =  A_B1R11[bn1_ma_reg4];
			   A_B1R3[bn1_ma_reg3]  =  data_tmp_1;
			   A_B0R3[bn0_ma_reg3]  =  data_tmp_2;
			   A_B0R3[bn0_ma_reg4]  =  data_tmp_3;
			   A_B1R3[bn1_ma_reg4]  =  data_tmp_4;
			   A_B1R11[bn1_ma_reg3] =  data_tmp_9;
			   A_B0R11[bn0_ma_reg3] =  data_tmp_10; 
			   A_B0R11[bn0_ma_reg4] =  data_tmp_11; 
			   A_B1R11[bn1_ma_reg4] =  data_tmp_12; 
			   //***********************************
			   data_tmp_1  = A_B1R5[bn1_ma_reg3];
			   data_tmp_2  = A_B1R6[bn1_ma_reg3];
			   data_tmp_3  = A_B1R7[bn1_ma_reg3];
			   data_tmp_9  = A_B1R13[bn1_ma_reg3];
			   data_tmp_10 = A_B1R14[bn1_ma_reg3];
			   data_tmp_11 = A_B1R15[bn1_ma_reg3];
			   A_B1R5[bn1_ma_reg3]  = A_B0R4[bn0_ma_reg3];
			   A_B1R6[bn1_ma_reg3]  = A_B0R4[bn0_ma_reg4];
			   A_B1R7[bn1_ma_reg3]  = A_B1R4[bn1_ma_reg4];
			   A_B1R13[bn1_ma_reg3] = A_B0R12[bn0_ma_reg3];
			   A_B1R14[bn1_ma_reg3] = A_B0R12[bn0_ma_reg4];
			   A_B1R15[bn1_ma_reg3] = A_B1R12[bn1_ma_reg4];
			   A_B0R4[bn0_ma_reg3]  = data_tmp_1;
			   A_B0R4[bn0_ma_reg4]  = data_tmp_2;
			   A_B1R4[bn1_ma_reg4]  = data_tmp_3;
			   A_B0R12[bn0_ma_reg3] = data_tmp_9;
			   A_B0R12[bn0_ma_reg4] = data_tmp_10;
			   A_B1R12[bn1_ma_reg4] = data_tmp_11;
			   //*********************************
			   data_tmp_1  =  A_B0R6[bn0_ma_reg3];
			   data_tmp_2  =  A_B0R7[bn0_ma_reg3];
			   data_tmp_9  =  A_B0R14[bn0_ma_reg3];
			   data_tmp_10 =  A_B0R15[bn0_ma_reg3];
			   A_B0R6[bn0_ma_reg3]  = A_B0R5[bn0_ma_reg4];
			   A_B0R7[bn0_ma_reg3]  = A_B1R5[bn1_ma_reg4];
			   A_B0R14[bn0_ma_reg3] = A_B0R13[bn0_ma_reg4];
			   A_B0R15[bn0_ma_reg3] = A_B1R13[bn1_ma_reg4];
			   A_B0R5[bn0_ma_reg4]  = data_tmp_1;
			   A_B1R5[bn1_ma_reg4]  = data_tmp_2;
			   A_B0R13[bn0_ma_reg4] = data_tmp_9;
			   A_B1R13[bn1_ma_reg4] = data_tmp_10;
			   //********************************
			   data_tmp_1 = A_B0R7[bn0_ma_reg4];
			   data_tmp_9 = A_B0R15[bn0_ma_reg4];
			   A_B0R7[bn0_ma_reg4]  = A_B1R6[bn1_ma_reg4];
			   A_B0R15[bn0_ma_reg4] = A_B1R14[bn1_ma_reg4];
			   A_B1R6[bn1_ma_reg4]  = data_tmp_1;
			   A_B1R14[bn1_ma_reg4] = data_tmp_9;
			   //-----------------------------------------------------------
			   //***********************************************************
			   //-----------------------------------------------------------
			   data_tmp_1  = A_B1R1[bn1_ma_reg5];
			   data_tmp_2  = A_B1R2[bn1_ma_reg5];
			   data_tmp_3  = A_B1R3[bn1_ma_reg5];
			   data_tmp_4  = A_B1R4[bn1_ma_reg5];
			   data_tmp_5  = A_B1R5[bn1_ma_reg5];
			   data_tmp_6  = A_B1R6[bn1_ma_reg5];
			   data_tmp_7  = A_B1R7[bn1_ma_reg5];
			   data_tmp_9  = A_B1R9[bn1_ma_reg5];
			   data_tmp_10 = A_B1R10[bn1_ma_reg5];
			   data_tmp_11 = A_B1R11[bn1_ma_reg5];
			   data_tmp_12 = A_B1R12[bn1_ma_reg5];
			   data_tmp_13 = A_B1R13[bn1_ma_reg5];
			   data_tmp_14 = A_B1R14[bn1_ma_reg5];
			   data_tmp_15 = A_B1R15[bn1_ma_reg5];
			   A_B1R1[bn1_ma_reg5]  =  A_B0R0[bn0_ma_reg5];
			   A_B1R2[bn1_ma_reg5]  =  A_B0R0[bn0_ma_reg6];
			   A_B1R3[bn1_ma_reg5]  =  A_B1R0[bn1_ma_reg6];
			   A_B1R4[bn1_ma_reg5]  =  A_B0R0[bn0_ma_reg7];
			   A_B1R5[bn1_ma_reg5]  =  A_B1R0[bn1_ma_reg7];
			   A_B1R6[bn1_ma_reg5]  =  A_B1R0[bn1_ma_reg8];
			   A_B1R7[bn1_ma_reg5]  =  A_B0R0[bn0_ma_reg8];
			   A_B1R9[bn1_ma_reg5]  =  A_B0R8[bn0_ma_reg5];
			   A_B1R10[bn1_ma_reg5] =  A_B0R8[bn0_ma_reg6];
			   A_B1R11[bn1_ma_reg5] =  A_B1R8[bn1_ma_reg6];
			   A_B1R12[bn1_ma_reg5] =  A_B0R8[bn0_ma_reg7];
			   A_B1R13[bn1_ma_reg5] =  A_B1R8[bn1_ma_reg7];
			   A_B1R14[bn1_ma_reg5] =  A_B1R8[bn1_ma_reg8];
			   A_B1R15[bn1_ma_reg5] =  A_B0R8[bn0_ma_reg8];
			   A_B0R0[bn0_ma_reg5]  =  data_tmp_1;
			   A_B0R0[bn0_ma_reg6]  =  data_tmp_2;
			   A_B1R0[bn1_ma_reg6]  =  data_tmp_3;
			   A_B0R0[bn0_ma_reg7]  =  data_tmp_4;
			   A_B1R0[bn1_ma_reg7]  =  data_tmp_5;
			   A_B1R0[bn1_ma_reg8]  =  data_tmp_6;
			   A_B0R0[bn0_ma_reg8]  =  data_tmp_7;
			   A_B0R8[bn0_ma_reg5]  =  data_tmp_9;
			   A_B0R8[bn0_ma_reg6]  =  data_tmp_10; 
			   A_B1R8[bn1_ma_reg6]  =  data_tmp_11; 
			   A_B0R8[bn0_ma_reg7]  =  data_tmp_12; 
			   A_B1R8[bn1_ma_reg7]  =  data_tmp_13; 
			   A_B1R8[bn1_ma_reg8]  =  data_tmp_14; 
			   A_B0R8[bn0_ma_reg8]  =  data_tmp_15; 
			   //**********************************
			   data_tmp_1  = A_B0R2[bn0_ma_reg5];
			   data_tmp_2  = A_B0R3[bn0_ma_reg5];
			   data_tmp_3  = A_B0R4[bn0_ma_reg5];
			   data_tmp_4  = A_B0R5[bn0_ma_reg5];
			   data_tmp_5  = A_B0R6[bn0_ma_reg5];
			   data_tmp_6  = A_B0R7[bn0_ma_reg5];
			   data_tmp_9  = A_B0R10[bn0_ma_reg5];
			   data_tmp_10 = A_B0R11[bn0_ma_reg5];
			   data_tmp_11 = A_B0R12[bn0_ma_reg5];
			   data_tmp_12 = A_B0R13[bn0_ma_reg5];
			   data_tmp_13 = A_B0R14[bn0_ma_reg5];
			   data_tmp_14 = A_B0R15[bn0_ma_reg5];
			   A_B0R2[bn0_ma_reg5]  =  A_B0R1[bn0_ma_reg6];
			   A_B0R3[bn0_ma_reg5]  =  A_B1R1[bn1_ma_reg6];
			   A_B0R4[bn0_ma_reg5]  =  A_B0R1[bn0_ma_reg7];
			   A_B0R5[bn0_ma_reg5]  =  A_B1R1[bn1_ma_reg7];
			   A_B0R6[bn0_ma_reg5]  =  A_B1R1[bn1_ma_reg8];
			   A_B0R7[bn0_ma_reg5]  =  A_B0R1[bn0_ma_reg8];
			   A_B0R10[bn0_ma_reg5] =  A_B0R9[bn0_ma_reg6];
			   A_B0R11[bn0_ma_reg5] =  A_B1R9[bn1_ma_reg6];
			   A_B0R12[bn0_ma_reg5] =  A_B0R9[bn0_ma_reg7];
			   A_B0R13[bn0_ma_reg5] =  A_B1R9[bn1_ma_reg7];
			   A_B0R14[bn0_ma_reg5] =  A_B1R9[bn1_ma_reg8];
			   A_B0R15[bn0_ma_reg5] =  A_B0R9[bn0_ma_reg8];
			   A_B0R1[bn0_ma_reg6]  =  data_tmp_1;
			   A_B1R1[bn1_ma_reg6]  =  data_tmp_2;
			   A_B0R1[bn0_ma_reg7]  =  data_tmp_3;
			   A_B1R1[bn1_ma_reg7]  =  data_tmp_4;
			   A_B1R1[bn1_ma_reg8]  =  data_tmp_5;
			   A_B0R1[bn0_ma_reg8]  =  data_tmp_6;
			   A_B0R9[bn0_ma_reg6]  =  data_tmp_9;
			   A_B1R9[bn1_ma_reg6]  =  data_tmp_10;
			   A_B0R9[bn0_ma_reg7]  =  data_tmp_11;
			   A_B1R9[bn1_ma_reg7]  =  data_tmp_12;
			   A_B1R9[bn1_ma_reg8]  =  data_tmp_13;
			   A_B0R9[bn0_ma_reg8]  =  data_tmp_14;
			   //************************************
			   data_tmp_1  =  A_B0R3[bn0_ma_reg6];
			   data_tmp_2  =  A_B0R4[bn0_ma_reg6];
			   data_tmp_3  =  A_B0R5[bn0_ma_reg6];
			   data_tmp_4  =  A_B0R6[bn0_ma_reg6];
			   data_tmp_5  =  A_B0R7[bn0_ma_reg6];
			   data_tmp_9  =  A_B0R11[bn0_ma_reg6];
			   data_tmp_10 =  A_B0R12[bn0_ma_reg6];
			   data_tmp_11 =  A_B0R13[bn0_ma_reg6];
			   data_tmp_12 =  A_B0R14[bn0_ma_reg6];
			   data_tmp_13 =  A_B0R15[bn0_ma_reg6];
			   A_B0R3[bn0_ma_reg6]  =  A_B1R2[bn1_ma_reg6];
			   A_B0R4[bn0_ma_reg6]  =  A_B0R2[bn0_ma_reg7];
			   A_B0R5[bn0_ma_reg6]  =  A_B1R2[bn1_ma_reg7];
			   A_B0R6[bn0_ma_reg6]  =  A_B1R2[bn1_ma_reg8];
			   A_B0R7[bn0_ma_reg6]  =  A_B0R2[bn0_ma_reg8];
			   A_B0R11[bn0_ma_reg6] =  A_B1R10[bn1_ma_reg6];
			   A_B0R12[bn0_ma_reg6] =  A_B0R10[bn0_ma_reg7];
			   A_B0R13[bn0_ma_reg6] =  A_B1R10[bn1_ma_reg7];
			   A_B0R14[bn0_ma_reg6] =  A_B1R10[bn1_ma_reg8];
			   A_B0R15[bn0_ma_reg6] =  A_B0R10[bn0_ma_reg8];
			   A_B1R2[bn1_ma_reg6]  =  data_tmp_1;
			   A_B0R2[bn0_ma_reg7]  =  data_tmp_2;
			   A_B1R2[bn1_ma_reg7]  =  data_tmp_3;
			   A_B1R2[bn1_ma_reg8]  =  data_tmp_4;
			   A_B0R2[bn0_ma_reg8]  =  data_tmp_5;
			   A_B1R10[bn1_ma_reg6]  =  data_tmp_9;
			   A_B0R10[bn0_ma_reg7]  =  data_tmp_10;
			   A_B1R10[bn1_ma_reg7]  =  data_tmp_11;
			   A_B1R10[bn1_ma_reg8]  =  data_tmp_12;
			   A_B0R10[bn0_ma_reg8]  =  data_tmp_13;
			   //***********************************
               data_tmp_1  =  A_B1R4[bn1_ma_reg6];
               data_tmp_2  =  A_B1R5[bn1_ma_reg6];
               data_tmp_3  =  A_B1R6[bn1_ma_reg6];
               data_tmp_4  =  A_B1R7[bn1_ma_reg6];
               data_tmp_9  =  A_B1R12[bn1_ma_reg6];
               data_tmp_10 =  A_B1R13[bn1_ma_reg6];
               data_tmp_11 =  A_B1R14[bn1_ma_reg6];
               data_tmp_12 =  A_B1R15[bn1_ma_reg6];
               A_B1R4[bn1_ma_reg6]  =   A_B0R3[bn0_ma_reg7];
			   A_B1R5[bn1_ma_reg6]  =   A_B1R3[bn1_ma_reg7];
			   A_B1R6[bn1_ma_reg6]  =   A_B1R3[bn1_ma_reg8];
			   A_B1R7[bn1_ma_reg6]  =   A_B0R3[bn0_ma_reg8];
			   A_B1R12[bn1_ma_reg6] =   A_B0R11[bn0_ma_reg7];
               A_B1R13[bn1_ma_reg6] =   A_B1R11[bn1_ma_reg7];
               A_B1R14[bn1_ma_reg6] =   A_B1R11[bn1_ma_reg8];
               A_B1R15[bn1_ma_reg6] =   A_B0R11[bn0_ma_reg8];
               A_B0R3[bn0_ma_reg7]  =   data_tmp_1; 
			   A_B1R3[bn1_ma_reg7]  =   data_tmp_2; 
			   A_B1R3[bn1_ma_reg8]  =   data_tmp_3; 
			   A_B0R3[bn0_ma_reg8]  =   data_tmp_4; 
			   A_B0R11[bn0_ma_reg7] =   data_tmp_9; 
			   A_B1R11[bn1_ma_reg7] =   data_tmp_10;
               A_B1R11[bn1_ma_reg8] =   data_tmp_11;
               A_B0R11[bn0_ma_reg8] =   data_tmp_12;
               //*******************************************
			   data_tmp_1  =  A_B0R5[bn0_ma_reg7];
			   data_tmp_2  =  A_B0R6[bn0_ma_reg7];
			   data_tmp_3  =  A_B0R7[bn0_ma_reg7];
			   data_tmp_9  =  A_B0R13[bn0_ma_reg7];
			   data_tmp_10 =  A_B0R14[bn0_ma_reg7];
			   data_tmp_11 =  A_B0R15[bn0_ma_reg7];
			   A_B0R5[bn0_ma_reg7]  =  A_B1R4[bn1_ma_reg7];
			   A_B0R6[bn0_ma_reg7]  =  A_B1R4[bn1_ma_reg8];
			   A_B0R7[bn0_ma_reg7]  =  A_B0R4[bn0_ma_reg8];
			   A_B0R13[bn0_ma_reg7] =  A_B1R12[bn1_ma_reg7];
			   A_B0R14[bn0_ma_reg7] =  A_B1R12[bn1_ma_reg8];
			   A_B0R15[bn0_ma_reg7] =  A_B0R12[bn0_ma_reg8];
			   A_B1R4[bn1_ma_reg7]  =  data_tmp_1; 
			   A_B1R4[bn1_ma_reg8]  =  data_tmp_2; 
			   A_B0R4[bn0_ma_reg8]  =  data_tmp_3; 
			   A_B1R12[bn1_ma_reg7] =  data_tmp_9; 
			   A_B1R12[bn1_ma_reg8] =  data_tmp_10;
			   A_B0R12[bn0_ma_reg8] =  data_tmp_11;
			   //*******************************************
			   data_tmp_1  =  A_B1R6[bn1_ma_reg7];
			   data_tmp_2  =  A_B1R7[bn1_ma_reg7];
			   data_tmp_9  =  A_B1R14[bn1_ma_reg7];
			   data_tmp_10 =  A_B1R15[bn1_ma_reg7];
			   A_B1R6[bn1_ma_reg7]  = A_B1R5[bn1_ma_reg8];
			   A_B1R7[bn1_ma_reg7]  = A_B0R5[bn0_ma_reg8];
			   A_B1R14[bn1_ma_reg7] = A_B1R13[bn1_ma_reg8];
			   A_B1R15[bn1_ma_reg7] = A_B0R13[bn0_ma_reg8];
               A_B1R5[bn1_ma_reg8]  = data_tmp_1;  
               A_B0R5[bn0_ma_reg8]  = data_tmp_2;  
               A_B1R13[bn1_ma_reg8] = data_tmp_9;  
			   A_B0R13[bn0_ma_reg8] = data_tmp_10;
               //******************************************
               data_tmp_1  =  A_B1R7[bn1_ma_reg8];
			   data_tmp_9  =  A_B1R15[bn1_ma_reg8];
               A_B1R7[bn1_ma_reg8]   =  A_B0R6[bn0_ma_reg8];
			   A_B1R15[bn1_ma_reg8]  =  A_B0R14[bn0_ma_reg8];
               A_B0R6[bn0_ma_reg8]   =  data_tmp_1;
			   A_B0R14[bn0_ma_reg8]  =  data_tmp_9;
			}
            else{
			   data_tmp_1  = A_B1R1[bn1_ma_reg1];
			   data_tmp_2  = A_B1R2[bn1_ma_reg1];
			   data_tmp_3  = A_B1R3[bn1_ma_reg1];
			   data_tmp_4  = A_B1R4[bn1_ma_reg1];
			   data_tmp_5  = A_B1R5[bn1_ma_reg1];
			   data_tmp_6  = A_B1R6[bn1_ma_reg1];
			   data_tmp_7  = A_B1R7[bn1_ma_reg1];
			   data_tmp_9  = A_B1R9[bn1_ma_reg1];
			   data_tmp_10 = A_B1R10[bn1_ma_reg1];
			   data_tmp_11 = A_B1R11[bn1_ma_reg1];
			   data_tmp_12 = A_B1R12[bn1_ma_reg1];
			   data_tmp_13 = A_B1R13[bn1_ma_reg1];
			   data_tmp_14 = A_B1R14[bn1_ma_reg1];
			   data_tmp_15 = A_B1R15[bn1_ma_reg1];
			   A_B1R1[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg1];
			   A_B1R2[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg2];
			   A_B1R3[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg2];
			   A_B1R4[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg3];
			   A_B1R5[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg3];
			   A_B1R6[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg4];
			   A_B1R7[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg4];
			   A_B1R9[bn1_ma_reg1]  =  A_B0R8[bn0_ma_reg1];
			   A_B1R10[bn1_ma_reg1] =  A_B0R8[bn0_ma_reg2];
			   A_B1R11[bn1_ma_reg1] =  A_B1R8[bn1_ma_reg2];
			   A_B1R12[bn1_ma_reg1] =  A_B0R8[bn0_ma_reg3];
			   A_B1R13[bn1_ma_reg1] =  A_B1R8[bn1_ma_reg3];
			   A_B1R14[bn1_ma_reg1] =  A_B1R8[bn1_ma_reg4];
			   A_B1R15[bn1_ma_reg1] =  A_B0R8[bn0_ma_reg4];
			   A_B0R0[bn0_ma_reg1]  =  data_tmp_1;
			   A_B0R0[bn0_ma_reg2]  =  data_tmp_2;
			   A_B1R0[bn1_ma_reg2]  =  data_tmp_3;
			   A_B0R0[bn0_ma_reg3]  =  data_tmp_4;
			   A_B1R0[bn1_ma_reg3]  =  data_tmp_5;
			   A_B1R0[bn1_ma_reg4]  =  data_tmp_6;
			   A_B0R0[bn0_ma_reg4]  =  data_tmp_7;
			   A_B0R8[bn0_ma_reg1]  =  data_tmp_9;
			   A_B0R8[bn0_ma_reg2]  =  data_tmp_10; 
			   A_B1R8[bn1_ma_reg2]  =  data_tmp_11; 
			   A_B0R8[bn0_ma_reg3]  =  data_tmp_12; 
			   A_B1R8[bn1_ma_reg3]  =  data_tmp_13; 
			   A_B1R8[bn1_ma_reg4]  =  data_tmp_14; 
			   A_B0R8[bn0_ma_reg4]  =  data_tmp_15; 
			   //**********************************
			   data_tmp_1  = A_B0R2[bn0_ma_reg1];
			   data_tmp_2  = A_B0R3[bn0_ma_reg1];
			   data_tmp_3  = A_B0R4[bn0_ma_reg1];
			   data_tmp_4  = A_B0R5[bn0_ma_reg1];
			   data_tmp_5  = A_B0R6[bn0_ma_reg1];
			   data_tmp_6  = A_B0R7[bn0_ma_reg1];
			   data_tmp_9  = A_B0R10[bn0_ma_reg1];
			   data_tmp_10 = A_B0R11[bn0_ma_reg1];
			   data_tmp_11 = A_B0R12[bn0_ma_reg1];
			   data_tmp_12 = A_B0R13[bn0_ma_reg1];
			   data_tmp_13 = A_B0R14[bn0_ma_reg1];
			   data_tmp_14 = A_B0R15[bn0_ma_reg1];
			   A_B0R2[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg2];
			   A_B0R3[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg2];
			   A_B0R4[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg3];
			   A_B0R5[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg3];
			   A_B0R6[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg4];
			   A_B0R7[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg4];
			   A_B0R10[bn0_ma_reg1] =  A_B0R9[bn0_ma_reg2];
			   A_B0R11[bn0_ma_reg1] =  A_B1R9[bn1_ma_reg2];
			   A_B0R12[bn0_ma_reg1] =  A_B0R9[bn0_ma_reg3];
			   A_B0R13[bn0_ma_reg1] =  A_B1R9[bn1_ma_reg3];
			   A_B0R14[bn0_ma_reg1] =  A_B1R9[bn1_ma_reg4];
			   A_B0R15[bn0_ma_reg1] =  A_B0R9[bn0_ma_reg4];
			   A_B0R1[bn0_ma_reg2]  =  data_tmp_1;
			   A_B1R1[bn1_ma_reg2]  =  data_tmp_2;
			   A_B0R1[bn0_ma_reg3]  =  data_tmp_3;
			   A_B1R1[bn1_ma_reg3]  =  data_tmp_4;
			   A_B1R1[bn1_ma_reg4]  =  data_tmp_5;
			   A_B0R1[bn0_ma_reg4]  =  data_tmp_6;
			   A_B0R9[bn0_ma_reg2]  =  data_tmp_9;
			   A_B1R9[bn1_ma_reg2]  =  data_tmp_10;
			   A_B0R9[bn0_ma_reg3]  =  data_tmp_11;
			   A_B1R9[bn1_ma_reg3]  =  data_tmp_12;
			   A_B1R9[bn1_ma_reg4]  =  data_tmp_13;
			   A_B0R9[bn0_ma_reg4]  =  data_tmp_14;
			   //************************************
			   data_tmp_1  =  A_B0R3[bn0_ma_reg2];
			   data_tmp_2  =  A_B0R4[bn0_ma_reg2];
			   data_tmp_3  =  A_B0R5[bn0_ma_reg2];
			   data_tmp_4  =  A_B0R6[bn0_ma_reg2];
			   data_tmp_5  =  A_B0R7[bn0_ma_reg2];
			   data_tmp_9  =  A_B0R11[bn0_ma_reg2];
			   data_tmp_10 =  A_B0R12[bn0_ma_reg2];
			   data_tmp_11 =  A_B0R13[bn0_ma_reg2];
			   data_tmp_12 =  A_B0R14[bn0_ma_reg2];
			   data_tmp_13 =  A_B0R15[bn0_ma_reg2];
			   A_B0R3[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg2];
			   A_B0R4[bn0_ma_reg2]  =  A_B0R2[bn0_ma_reg3];
			   A_B0R5[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg3];
			   A_B0R6[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg4];
			   A_B0R7[bn0_ma_reg2]  =  A_B0R2[bn0_ma_reg4];
			   A_B0R11[bn0_ma_reg2] =  A_B1R10[bn1_ma_reg2];
			   A_B0R12[bn0_ma_reg2] =  A_B0R10[bn0_ma_reg3];
			   A_B0R13[bn0_ma_reg2] =  A_B1R10[bn1_ma_reg3];
			   A_B0R14[bn0_ma_reg2] =  A_B1R10[bn1_ma_reg4];
			   A_B0R15[bn0_ma_reg2] =  A_B0R10[bn0_ma_reg4];
			   A_B1R2[bn1_ma_reg2]  =  data_tmp_1;
			   A_B0R2[bn0_ma_reg3]  =  data_tmp_2;
			   A_B1R2[bn1_ma_reg3]  =  data_tmp_3;
			   A_B1R2[bn1_ma_reg4]  =  data_tmp_4;
			   A_B0R2[bn0_ma_reg4]  =  data_tmp_5;
			   A_B1R10[bn1_ma_reg2]  =  data_tmp_9;
			   A_B0R10[bn0_ma_reg3]  =  data_tmp_10;
			   A_B1R10[bn1_ma_reg3]  =  data_tmp_11;
			   A_B1R10[bn1_ma_reg4]  =  data_tmp_12;
			   A_B0R10[bn0_ma_reg4]  =  data_tmp_13;
			   //***********************************
               data_tmp_1  =  A_B1R4[bn1_ma_reg2];
               data_tmp_2  =  A_B1R5[bn1_ma_reg2];
               data_tmp_3  =  A_B1R6[bn1_ma_reg2];
               data_tmp_4  =  A_B1R7[bn1_ma_reg2];
               data_tmp_9  =  A_B1R12[bn1_ma_reg2];
               data_tmp_10 =  A_B1R13[bn1_ma_reg2];
               data_tmp_11 =  A_B1R14[bn1_ma_reg2];
               data_tmp_12 =  A_B1R15[bn1_ma_reg2];
               A_B1R4[bn1_ma_reg2]  =   A_B0R3[bn0_ma_reg3];
			   A_B1R5[bn1_ma_reg2]  =   A_B1R3[bn1_ma_reg3];
			   A_B1R6[bn1_ma_reg2]  =   A_B1R3[bn1_ma_reg4];
			   A_B1R7[bn1_ma_reg2]  =   A_B0R3[bn0_ma_reg4];
			   A_B1R12[bn1_ma_reg2] =   A_B0R11[bn0_ma_reg3];
               A_B1R13[bn1_ma_reg2] =   A_B1R11[bn1_ma_reg3];
               A_B1R14[bn1_ma_reg2] =   A_B1R11[bn1_ma_reg4];
               A_B1R15[bn1_ma_reg2] =   A_B0R11[bn0_ma_reg4];
               A_B0R3[bn0_ma_reg3]  =   data_tmp_1; 
			   A_B1R3[bn1_ma_reg3]  =   data_tmp_2; 
			   A_B1R3[bn1_ma_reg4]  =   data_tmp_3; 
			   A_B0R3[bn0_ma_reg4]  =   data_tmp_4; 
			   A_B0R11[bn0_ma_reg3] =   data_tmp_9; 
			   A_B1R11[bn1_ma_reg3] =   data_tmp_10;
               A_B1R11[bn1_ma_reg4] =   data_tmp_11;
               A_B0R11[bn0_ma_reg4] =   data_tmp_12;
               //*******************************************
			   data_tmp_1  =  A_B0R5[bn0_ma_reg3];
			   data_tmp_2  =  A_B0R6[bn0_ma_reg3];
			   data_tmp_3  =  A_B0R7[bn0_ma_reg3];
			   data_tmp_9  =  A_B0R13[bn0_ma_reg3];
			   data_tmp_10 =  A_B0R14[bn0_ma_reg3];
			   data_tmp_11 =  A_B0R15[bn0_ma_reg3];
			   A_B0R5[bn0_ma_reg3]  =  A_B1R4[bn1_ma_reg3];
			   A_B0R6[bn0_ma_reg3]  =  A_B1R4[bn1_ma_reg4];
			   A_B0R7[bn0_ma_reg3]  =  A_B0R4[bn0_ma_reg4];
			   A_B0R13[bn0_ma_reg3] =  A_B1R12[bn1_ma_reg3];
			   A_B0R14[bn0_ma_reg3] =  A_B1R12[bn1_ma_reg4];
			   A_B0R15[bn0_ma_reg3] =  A_B0R12[bn0_ma_reg4];
			   A_B1R4[bn1_ma_reg3]  =  data_tmp_1; 
			   A_B1R4[bn1_ma_reg4]  =  data_tmp_2; 
			   A_B0R4[bn0_ma_reg4]  =  data_tmp_3; 
			   A_B1R12[bn1_ma_reg3] =  data_tmp_9; 
			   A_B1R12[bn1_ma_reg4] =  data_tmp_10;
			   A_B0R12[bn0_ma_reg4] =  data_tmp_11;
			   //*******************************************
			   data_tmp_1  =  A_B1R6[bn1_ma_reg3];
			   data_tmp_2  =  A_B1R7[bn1_ma_reg3];
			   data_tmp_9  =  A_B1R14[bn1_ma_reg3];
			   data_tmp_10 =  A_B1R15[bn1_ma_reg3];
			   A_B1R6[bn1_ma_reg3]  = A_B1R5[bn1_ma_reg4];
			   A_B1R7[bn1_ma_reg3]  = A_B0R5[bn0_ma_reg4];
			   A_B1R14[bn1_ma_reg3] = A_B1R13[bn1_ma_reg4];
			   A_B1R15[bn1_ma_reg3] = A_B0R13[bn0_ma_reg4];
               A_B1R5[bn1_ma_reg4]  = data_tmp_1;  
               A_B0R5[bn0_ma_reg4]  = data_tmp_2;  
               A_B1R13[bn1_ma_reg4] = data_tmp_9;  
			   A_B0R13[bn0_ma_reg4] = data_tmp_10;
               //******************************************
               data_tmp_1  =  A_B1R7[bn1_ma_reg4];
			   data_tmp_9  =  A_B1R15[bn1_ma_reg4];
               A_B1R7[bn1_ma_reg4]   =  A_B0R6[bn0_ma_reg4];
			   A_B1R15[bn1_ma_reg4]  =  A_B0R14[bn0_ma_reg4];
               A_B0R6[bn0_ma_reg4]   =  data_tmp_1;
			   A_B0R14[bn0_ma_reg4]  =  data_tmp_9;	
			   //-----------------------------------------------------
			   //-----------------------------------------------------
			   data_tmp_1  = A_B0R1[bn0_ma_reg5];
			   data_tmp_2  = A_B0R2[bn0_ma_reg5];
			   data_tmp_3  = A_B0R3[bn0_ma_reg5];
			   data_tmp_4  = A_B0R4[bn0_ma_reg5];
			   data_tmp_5  = A_B0R5[bn0_ma_reg5];
			   data_tmp_6  = A_B0R6[bn0_ma_reg5];
			   data_tmp_7  = A_B0R7[bn0_ma_reg5];
			   data_tmp_9  = A_B0R9[bn0_ma_reg5];
			   data_tmp_10 = A_B0R10[bn0_ma_reg5];
			   data_tmp_11 = A_B0R11[bn0_ma_reg5];
			   data_tmp_12 = A_B0R12[bn0_ma_reg5];
			   data_tmp_13 = A_B0R13[bn0_ma_reg5];
			   data_tmp_14 = A_B0R14[bn0_ma_reg5];
			   data_tmp_15 = A_B0R15[bn0_ma_reg5];
			   A_B0R1[bn0_ma_reg5]  =  A_B1R0[bn1_ma_reg5];
			   A_B0R2[bn0_ma_reg5]  =  A_B1R0[bn1_ma_reg6];
			   A_B0R3[bn0_ma_reg5]  =  A_B0R0[bn0_ma_reg6];
			   A_B0R4[bn0_ma_reg5]  =  A_B1R0[bn1_ma_reg7];
			   A_B0R5[bn0_ma_reg5]  =  A_B0R0[bn0_ma_reg7];
			   A_B0R6[bn0_ma_reg5]  =  A_B0R0[bn0_ma_reg8];
			   A_B0R7[bn0_ma_reg5]  =  A_B1R0[bn1_ma_reg8];
			   A_B0R9[bn0_ma_reg5]  =  A_B1R8[bn1_ma_reg5];
			   A_B0R10[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg6];
			   A_B0R11[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg6];
			   A_B0R12[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg7];
			   A_B0R13[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg7];
			   A_B0R14[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg8];
			   A_B0R15[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg8];
			   A_B1R0[bn1_ma_reg5]  =  data_tmp_1;
			   A_B1R0[bn1_ma_reg6]  =  data_tmp_2;
			   A_B0R0[bn0_ma_reg6]  =  data_tmp_3;
			   A_B1R0[bn1_ma_reg7]  =  data_tmp_4;
			   A_B0R0[bn0_ma_reg7]  =  data_tmp_5;
			   A_B0R0[bn0_ma_reg8]  =  data_tmp_6;
			   A_B1R0[bn1_ma_reg8]  =  data_tmp_7;
			   A_B1R8[bn1_ma_reg5]  =  data_tmp_9;
			   A_B1R8[bn1_ma_reg6]  =  data_tmp_10; 
			   A_B0R8[bn0_ma_reg6]  =  data_tmp_11; 
			   A_B1R8[bn1_ma_reg7]  =  data_tmp_12; 
			   A_B0R8[bn0_ma_reg7]  =  data_tmp_13; 
			   A_B0R8[bn0_ma_reg8]  =  data_tmp_14; 
			   A_B1R8[bn1_ma_reg8]  =  data_tmp_15; 
			   //**********************************
			   data_tmp_1  = A_B1R2[bn1_ma_reg5];
			   data_tmp_2  = A_B1R3[bn1_ma_reg5];
			   data_tmp_3  = A_B1R4[bn1_ma_reg5];
			   data_tmp_4  = A_B1R5[bn1_ma_reg5];
			   data_tmp_5  = A_B1R6[bn1_ma_reg5];
			   data_tmp_6  = A_B1R7[bn1_ma_reg5];
			   data_tmp_9  = A_B1R10[bn1_ma_reg5];
			   data_tmp_10 = A_B1R11[bn1_ma_reg5];
			   data_tmp_11 = A_B1R12[bn1_ma_reg5];
			   data_tmp_12 = A_B1R13[bn1_ma_reg5];
			   data_tmp_13 = A_B1R14[bn1_ma_reg5];
			   data_tmp_14 = A_B1R15[bn1_ma_reg5];
			   A_B1R2[bn1_ma_reg5]   =  A_B1R1[bn1_ma_reg6];
			   A_B1R3[bn1_ma_reg5]   =  A_B0R1[bn0_ma_reg6];
			   A_B1R4[bn1_ma_reg5]   =  A_B1R1[bn1_ma_reg7];
			   A_B1R5[bn1_ma_reg5]   =  A_B0R1[bn0_ma_reg7];
			   A_B1R6[bn1_ma_reg5]   =  A_B0R1[bn0_ma_reg8];
			   A_B1R7[bn1_ma_reg5]   =  A_B1R1[bn1_ma_reg8];
			   A_B1R10[bn1_ma_reg5]  =  A_B1R9[bn1_ma_reg6];
			   A_B1R11[bn1_ma_reg5]  =  A_B0R9[bn0_ma_reg6];
			   A_B1R12[bn1_ma_reg5]  =  A_B1R9[bn1_ma_reg7];
			   A_B1R13[bn1_ma_reg5]  =  A_B0R9[bn0_ma_reg7];
			   A_B1R14[bn1_ma_reg5]  =  A_B0R9[bn0_ma_reg8];
			   A_B1R15[bn1_ma_reg5]  =  A_B1R9[bn1_ma_reg8];
			   A_B1R1[bn1_ma_reg6] =  data_tmp_1;
			   A_B0R1[bn0_ma_reg6] =  data_tmp_2;
			   A_B1R1[bn1_ma_reg7] =  data_tmp_3;
			   A_B0R1[bn0_ma_reg7] =  data_tmp_4;
			   A_B0R1[bn0_ma_reg8] =  data_tmp_5;
			   A_B1R1[bn1_ma_reg8] =  data_tmp_6;
			   A_B1R9[bn1_ma_reg6] =  data_tmp_9;
			   A_B0R9[bn0_ma_reg6] =  data_tmp_10; 
			   A_B1R9[bn1_ma_reg7] =  data_tmp_11; 
			   A_B0R9[bn0_ma_reg7] =  data_tmp_12; 
			   A_B0R9[bn0_ma_reg8] =  data_tmp_13; 
			   A_B1R9[bn1_ma_reg8] =  data_tmp_14; 
			   //*************************************
			   data_tmp_1  =  A_B1R3[bn1_ma_reg6];
			   data_tmp_2  =  A_B1R4[bn1_ma_reg6];
			   data_tmp_3  =  A_B1R5[bn1_ma_reg6];
			   data_tmp_4  =  A_B1R6[bn1_ma_reg6];
			   data_tmp_5  =  A_B1R7[bn1_ma_reg6];
			   data_tmp_9  =  A_B1R11[bn1_ma_reg6];
			   data_tmp_10 =  A_B1R12[bn1_ma_reg6];
			   data_tmp_11 =  A_B1R13[bn1_ma_reg6];
			   data_tmp_12 =  A_B1R14[bn1_ma_reg6];
			   data_tmp_13 =  A_B1R15[bn1_ma_reg6];
			   A_B1R3[bn1_ma_reg6]  =  A_B0R2[bn0_ma_reg6];
			   A_B1R4[bn1_ma_reg6]  =  A_B1R2[bn1_ma_reg7];
			   A_B1R5[bn1_ma_reg6]  =  A_B0R2[bn0_ma_reg7];
			   A_B1R6[bn1_ma_reg6]  =  A_B0R2[bn0_ma_reg8];
			   A_B1R7[bn1_ma_reg6]  =  A_B1R2[bn1_ma_reg8];
			   A_B1R11[bn1_ma_reg6] =  A_B0R10[bn0_ma_reg6];
			   A_B1R12[bn1_ma_reg6] =  A_B1R10[bn1_ma_reg7];
			   A_B1R13[bn1_ma_reg6] =  A_B0R10[bn0_ma_reg7];
			   A_B1R14[bn1_ma_reg6] =  A_B0R10[bn0_ma_reg8];
			   A_B1R15[bn1_ma_reg6] =  A_B1R10[bn1_ma_reg8];
			   A_B0R2[bn0_ma_reg6]  =  data_tmp_1;
			   A_B1R2[bn1_ma_reg7]  =  data_tmp_2;
			   A_B0R2[bn0_ma_reg7]  =  data_tmp_3;
			   A_B0R2[bn0_ma_reg8]  =  data_tmp_4;
			   A_B1R2[bn1_ma_reg8]  =  data_tmp_5;
			   A_B0R10[bn0_ma_reg6] =  data_tmp_9;
			   A_B1R10[bn1_ma_reg7] =  data_tmp_10;
			   A_B0R10[bn0_ma_reg7] =  data_tmp_11;
			   A_B0R10[bn0_ma_reg8] =  data_tmp_12;
			   A_B1R10[bn1_ma_reg8] =  data_tmp_13;
			   //************************************
			   data_tmp_1  =  A_B0R4[bn0_ma_reg6];
			   data_tmp_2  =  A_B0R5[bn0_ma_reg6];
			   data_tmp_3  =  A_B0R6[bn0_ma_reg6];
			   data_tmp_4  =  A_B0R7[bn0_ma_reg6];
			   data_tmp_9  =  A_B0R12[bn0_ma_reg6];
			   data_tmp_10 =  A_B0R13[bn0_ma_reg6];
 			   data_tmp_11 =  A_B0R14[bn0_ma_reg6];
			   data_tmp_12 =  A_B0R15[bn0_ma_reg6];
			   A_B0R4[bn0_ma_reg6]  =  A_B1R3[bn1_ma_reg7];
			   A_B0R5[bn0_ma_reg6]  =  A_B0R3[bn0_ma_reg7];
			   A_B0R6[bn0_ma_reg6]  =  A_B0R3[bn0_ma_reg8];
			   A_B0R7[bn0_ma_reg6]  =  A_B1R3[bn1_ma_reg8];
			   A_B0R12[bn0_ma_reg6] =  A_B1R11[bn1_ma_reg7];
			   A_B0R13[bn0_ma_reg6] =  A_B0R11[bn0_ma_reg7];
			   A_B0R14[bn0_ma_reg6] =  A_B0R11[bn0_ma_reg8];
			   A_B0R15[bn0_ma_reg6] =  A_B1R11[bn1_ma_reg8];
			   A_B1R3[bn1_ma_reg7]  =  data_tmp_1;
			   A_B0R3[bn0_ma_reg7]  =  data_tmp_2;
			   A_B0R3[bn0_ma_reg8]  =  data_tmp_3;
			   A_B1R3[bn1_ma_reg8]  =  data_tmp_4;
			   A_B1R11[bn1_ma_reg7] =  data_tmp_9;
			   A_B0R11[bn0_ma_reg7] =  data_tmp_10; 
			   A_B0R11[bn0_ma_reg8] =  data_tmp_11; 
			   A_B1R11[bn1_ma_reg8] =  data_tmp_12; 
			   //***********************************
			   data_tmp_1  = A_B1R5[bn1_ma_reg7];
			   data_tmp_2  = A_B1R6[bn1_ma_reg7];
			   data_tmp_3  = A_B1R7[bn1_ma_reg7];
			   data_tmp_9  = A_B1R13[bn1_ma_reg7];
			   data_tmp_10 = A_B1R14[bn1_ma_reg7];
			   data_tmp_11 = A_B1R15[bn1_ma_reg7];
			   A_B1R5[bn1_ma_reg7]  = A_B0R4[bn0_ma_reg7];
			   A_B1R6[bn1_ma_reg7]  = A_B0R4[bn0_ma_reg8];
			   A_B1R7[bn1_ma_reg7]  = A_B1R4[bn1_ma_reg8];
			   A_B1R13[bn1_ma_reg7] = A_B0R12[bn0_ma_reg7];
			   A_B1R14[bn1_ma_reg7] = A_B0R12[bn0_ma_reg8];
			   A_B1R15[bn1_ma_reg7] = A_B1R12[bn1_ma_reg8];
			   A_B0R4[bn0_ma_reg7]  = data_tmp_1;
			   A_B0R4[bn0_ma_reg8]  = data_tmp_2;
			   A_B1R4[bn1_ma_reg8]  = data_tmp_3;
			   A_B0R12[bn0_ma_reg7] = data_tmp_9;
			   A_B0R12[bn0_ma_reg8] = data_tmp_10;
			   A_B1R12[bn1_ma_reg8] = data_tmp_11;
			   //*********************************
			   data_tmp_1  =  A_B0R6[bn0_ma_reg7];
			   data_tmp_2  =  A_B0R7[bn0_ma_reg7];
			   data_tmp_9  =  A_B0R14[bn0_ma_reg7];
			   data_tmp_10 =  A_B0R15[bn0_ma_reg7];
			   A_B0R6[bn0_ma_reg7]  = A_B0R5[bn0_ma_reg8];
			   A_B0R7[bn0_ma_reg7]  = A_B1R5[bn1_ma_reg8];
			   A_B0R14[bn0_ma_reg7] = A_B0R13[bn0_ma_reg8];
			   A_B0R15[bn0_ma_reg7] = A_B1R13[bn1_ma_reg8];
			   A_B0R5[bn0_ma_reg8]  = data_tmp_1;
			   A_B1R5[bn1_ma_reg8]  = data_tmp_2;
			   A_B0R13[bn0_ma_reg8] = data_tmp_9;
			   A_B1R13[bn1_ma_reg8] = data_tmp_10;
			   //********************************
			   data_tmp_1 = A_B0R7[bn0_ma_reg8];
			   data_tmp_9 = A_B0R15[bn0_ma_reg8];
			   A_B0R7[bn0_ma_reg8]  = A_B1R6[bn1_ma_reg8];
			   A_B0R15[bn0_ma_reg8] = A_B1R14[bn1_ma_reg8];
			   A_B1R6[bn1_ma_reg8]  = data_tmp_1;
			   A_B1R14[bn1_ma_reg8] = data_tmp_9;			   
            }
		 } 
		}
	}
	
	std::cout << "radix-16 FFT computing over!!\n";
		
	//if(display == 1)DATARECORD <<" ----------------------------------------------- \n";
	//if(display == 1)DATARECORD <<" Radix-8 FFT computing start!!! \n";
	
	int    right_shift_bits;
	double relocation_difference;
	
	relocation_difference = (double)group / 4;
	right_shift_bits = (int)log2(relocation_difference);
	
	for(int i = 0;i < group; i++){
        for(int j=0; j < radix; j++){
        	gray_i  = Gray(i,group);
        	BC_tmp  = j * group + gray_i;
			RR_R16_R8(BC_tmp,right_shift_bits,BC);
        	AGU_R16(BC,bn_tmp,ma_tmp);
			//if(display == 1)DATARECORD <<" ----------------------------------------------- \n";
			//if(display == 1)DATARECORD <<"BC: "<< BC <<" \n";
			//if(display == 1)DATARECORD <<"bn_tmp: "<< bn_tmp <<" \n";
			//if(display == 1)DATARECORD <<"ma_tmp: "<< ma_tmp <<" \n";
        	if(bn_tmp == 0){
                //if(display == 1)DATARECORD <<" Before radix-8 BU operation!! \n";					
                //if(display == 1)DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";
                Radix8_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp],A_B0R2[ma_tmp],A_B0R3[ma_tmp],A_B0R4[ma_tmp],A_B0R5[ma_tmp],A_B0R6[ma_tmp],A_B0R7[ma_tmp]);
        		Radix8_BU(A_B0R8[ma_tmp],A_B0R9[ma_tmp],A_B0R10[ma_tmp],A_B0R11[ma_tmp],A_B0R12[ma_tmp],A_B0R13[ma_tmp],A_B0R14[ma_tmp],A_B0R15[ma_tmp]);
                //if(display == 1)DATARECORD <<"After radix-8 BU operation!! \n";					
                //if(display == 1)DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";						 
        	}else {
                //if(display == 1)DATARECORD <<"Before radix-8 BU operation!! \n";					
                //if(display == 1)DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";				
        	    Radix8_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp],A_B1R2[ma_tmp],A_B1R3[ma_tmp],A_B1R4[ma_tmp],A_B1R5[ma_tmp],A_B1R6[ma_tmp],A_B1R7[ma_tmp]);
        	    Radix8_BU(A_B1R8[ma_tmp],A_B1R9[ma_tmp],A_B1R10[ma_tmp],A_B1R11[ma_tmp],A_B1R12[ma_tmp],A_B1R13[ma_tmp],A_B1R14[ma_tmp],A_B1R15[ma_tmp]);
                //if(display == 1)DATARECORD <<"After radix-8 BU operation!! \n";					
                //if(display == 1)DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";					
                //if(display == 1)DATARECORD <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";							  
        	}			
        }		
	}
	
	int index0_i;
	int index1_i;
	int index2_i;
	int index3_i;
	int index4_i;
	int index5_i;
	int index6_i;
	int index7_i;
	int index8_i;
	int index9_i;
	int index10_i;
	int index11_i;
	int index12_i;
	int index13_i;
	int index14_i;
	int index15_i;
	
	int index0;
    int index1;
    int index2;
    int index3;
    int index4;
    int index5;
    int index6;
    int index7;
    int index8;
    int index9;
    int index10;
    int index11;
    int index12;
    int index13;
    int index14;
    int index15;
	
	int index_128_flag;
	int index_3_BITS_FLAG;
	
	DATARECORD <<"***************************************************************\n";
	DATARECORD <<"***** DATA OUTPUT!!                                         ***\n";
	DATARECORD <<"***************************************************************\n";
	//data output
	// SPMB data output , output function is "NTT_REORDERINDEX_R16_R8_OUT".
	for(int i = 0;i < group;i++){
		for(int j = 0;j < radix;j++){
			//gray_i  = Gray(i,group);
			BC_tmp  = j * group + i;
			if(display == 1)DATARECORD <<"---------------------------------------------------\n";
			if(display == 1)DATARECORD <<" BC_tmp: "<< BC_tmp <<"\n";
			RR_R16_R8(BC_tmp,((4 * (Stage-1)) - 1),BC);
			if(display == 1)DATARECORD <<" BC: "<< BC <<"\n";
			index_128_flag   = BC >> 3;
			index_3_BITS_FLAG = BC % 8;
			AGU_R16(BC,bn_tmp,ma_tmp);
			//if(display == 1)DATARECORD <<" bn_tmp:  "<< bn_tmp <<"\n";
			//if(display == 1)DATARECORD <<" ma_tmp:  "<< ma_tmp <<"\n";
			index0_i  = index_128_flag * 128 + index_3_BITS_FLAG * 8 + 0;
			index1_i  = index_128_flag * 128 + index_3_BITS_FLAG * 8 + 1;
			index2_i  = index_128_flag * 128 + index_3_BITS_FLAG * 8 + 2;
			index3_i  = index_128_flag * 128 + index_3_BITS_FLAG * 8 + 3;
			index4_i  = index_128_flag * 128 + index_3_BITS_FLAG * 8 + 4;
			index5_i  = index_128_flag * 128 + index_3_BITS_FLAG * 8 + 5;
			index6_i  = index_128_flag * 128 + index_3_BITS_FLAG * 8 + 6;
			index7_i  = index_128_flag * 128 + index_3_BITS_FLAG * 8 + 7;
			index8_i  = index_128_flag * 128 + index_3_BITS_FLAG * 8 + 64;
			index9_i  = index_128_flag * 128 + index_3_BITS_FLAG * 8 + 65;
			index10_i = index_128_flag * 128 + index_3_BITS_FLAG * 8 + 66;
			index11_i = index_128_flag * 128 + index_3_BITS_FLAG * 8 + 67;
			index12_i = index_128_flag * 128 + index_3_BITS_FLAG * 8 + 68;
			index13_i = index_128_flag * 128 + index_3_BITS_FLAG * 8 + 69;
			index14_i = index_128_flag * 128 + index_3_BITS_FLAG * 8 + 70;
			index15_i = index_128_flag * 128 + index_3_BITS_FLAG * 8 + 71;

			//if(display == 1)DATARECORD <<" index0_i:  "<< index0_i <<"\n";
			//if(display == 1)DATARECORD <<" index1_i:  "<< index1_i <<"\n";
			//if(display == 1)DATARECORD <<" index2_i:  "<< index2_i <<"\n";
			//if(display == 1)DATARECORD <<" index3_i:  "<< index3_i <<"\n";
			//if(display == 1)DATARECORD <<" index4_i:  "<< index4_i <<"\n";
			//if(display == 1)DATARECORD <<" index5_i:  "<< index5_i <<"\n";
			//if(display == 1)DATARECORD <<" index6_i:  "<< index6_i <<"\n";
			//if(display == 1)DATARECORD <<" index7_i:  "<< index7_i <<"\n";
			//if(display == 1)DATARECORD <<" index8_i:  "<< index8_i <<"\n";
			//if(display == 1)DATARECORD <<" index9_i:  "<< index9_i <<"\n";
			//if(display == 1)DATARECORD <<" index10_i: "<< index10_i <<"\n";
			//if(display == 1)DATARECORD <<" index11_i: "<< index11_i <<"\n";
			//if(display == 1)DATARECORD <<" index12_i: "<< index12_i <<"\n";
			//if(display == 1)DATARECORD <<" index13_i: "<< index13_i <<"\n";
			//if(display == 1)DATARECORD <<" index14_i: "<< index14_i <<"\n";
			//if(display == 1)DATARECORD <<" index15_i: "<< index15_i <<"\n";
            NTT_REORDERINDEX_R16_R8_OUT(index0_i,index0);			
            NTT_REORDERINDEX_R16_R8_OUT(index1_i,index1);			
            NTT_REORDERINDEX_R16_R8_OUT(index2_i,index2);			
            NTT_REORDERINDEX_R16_R8_OUT(index3_i,index3);			
            NTT_REORDERINDEX_R16_R8_OUT(index4_i,index4);			
            NTT_REORDERINDEX_R16_R8_OUT(index5_i,index5);			
            NTT_REORDERINDEX_R16_R8_OUT(index6_i,index6);			
            NTT_REORDERINDEX_R16_R8_OUT(index7_i,index7);			
            NTT_REORDERINDEX_R16_R8_OUT(index8_i,index8);			
            NTT_REORDERINDEX_R16_R8_OUT(index9_i,index9);			
            NTT_REORDERINDEX_R16_R8_OUT(index10_i,index10);			
            NTT_REORDERINDEX_R16_R8_OUT(index11_i,index11);			
            NTT_REORDERINDEX_R16_R8_OUT(index12_i,index12);			
            NTT_REORDERINDEX_R16_R8_OUT(index13_i,index13);			
            NTT_REORDERINDEX_R16_R8_OUT(index14_i,index14);			
            NTT_REORDERINDEX_R16_R8_OUT(index15_i,index15);		

            //if(display == 1)DATARECORD <<" index0:  "<< index0 <<"\n";
			//if(display == 1)DATARECORD <<" index1:  "<< index1 <<"\n";
			//if(display == 1)DATARECORD <<" index2:  "<< index2 <<"\n";
			//if(display == 1)DATARECORD <<" index3:  "<< index3 <<"\n";
			//if(display == 1)DATARECORD <<" index4:  "<< index4 <<"\n";
			//if(display == 1)DATARECORD <<" index5:  "<< index5 <<"\n";
			//if(display == 1)DATARECORD <<" index6:  "<< index6 <<"\n";
			//if(display == 1)DATARECORD <<" index7:  "<< index7 <<"\n";
			//if(display == 1)DATARECORD <<" index8:  "<< index8 <<"\n";
			//if(display == 1)DATARECORD <<" index9:  "<< index9 <<"\n";
			//if(display == 1)DATARECORD <<" index10: "<< index10 <<"\n";
			//if(display == 1)DATARECORD <<" index11: "<< index11 <<"\n";
			//if(display == 1)DATARECORD <<" index12: "<< index12 <<"\n";
			//if(display == 1)DATARECORD <<" index13: "<< index13 <<"\n";
			//if(display == 1)DATARECORD <<" index14: "<< index14 <<"\n";
			//if(display == 1)DATARECORD <<" index15: "<< index15 <<"\n";
			
			if(bn_tmp == 0){
               A[index0]  = A_B0R0[ma_tmp];
			   A[index1]  = A_B0R1[ma_tmp];
			   A[index2]  = A_B0R2[ma_tmp];
			   A[index3]  = A_B0R3[ma_tmp];
			   A[index4]  = A_B0R4[ma_tmp];
			   A[index5]  = A_B0R5[ma_tmp];
			   A[index6]  = A_B0R6[ma_tmp];
			   A[index7]  = A_B0R7[ma_tmp];
			   A[index8]  = A_B0R8[ma_tmp];
			   A[index9]  = A_B0R9[ma_tmp];
			   A[index10] = A_B0R10[ma_tmp];
			   A[index11] = A_B0R11[ma_tmp];
			   A[index12] = A_B0R12[ma_tmp];
			   A[index13] = A_B0R13[ma_tmp];
			   A[index14] = A_B0R14[ma_tmp];
			   A[index15] = A_B0R15[ma_tmp];
			   //if(display == 1)DATARECORD << "A_B0R0[" << ma_tmp <<"]: "<< A_B0R0[ma_tmp] << "\n";
			   //if(display == 1)DATARECORD << "A_B0R1[" << ma_tmp <<"]: "<< A_B0R1[ma_tmp] << "\n";
			   //if(display == 1)DATARECORD << "A_B0R2[" << ma_tmp <<"]: "<< A_B0R2[ma_tmp] << "\n";
			   //if(display == 1)DATARECORD << "A_B0R3[" << ma_tmp <<"]: "<< A_B0R3[ma_tmp] << "\n";
			   //if(display == 1)DATARECORD << "A_B0R4[" << ma_tmp <<"]: "<< A_B0R4[ma_tmp] << "\n";
			   //if(display == 1)DATARECORD << "A_B0R5[" << ma_tmp <<"]: "<< A_B0R5[ma_tmp] << "\n";
			   //if(display == 1)DATARECORD << "A_B0R6[" << ma_tmp <<"]: "<< A_B0R6[ma_tmp] << "\n";
			   //if(display == 1)DATARECORD << "A_B0R7[" << ma_tmp <<"]: "<< A_B0R7[ma_tmp] << "\n";
			   //if(display == 1)DATARECORD << "A_B0R8[" << ma_tmp <<"]: "<< A_B0R8[ma_tmp] << "\n";
			   //if(display == 1)DATARECORD << "A_B0R9[" << ma_tmp <<"]: "<< A_B0R9[ma_tmp] << "\n";
			   //if(display == 1)DATARECORD << "A_B0R10[" << ma_tmp <<"]: "<< A_B0R10[ma_tmp] << "\n";					
			   //if(display == 1)DATARECORD << "A_B0R11[" << ma_tmp <<"]: "<< A_B0R11[ma_tmp] << "\n";					
			   //if(display == 1)DATARECORD << "A_B0R12[" << ma_tmp <<"]: "<< A_B0R12[ma_tmp] << "\n";					
			   //if(display == 1)DATARECORD << "A_B0R13[" << ma_tmp <<"]: "<< A_B0R13[ma_tmp] << "\n";					
			   //if(display == 1)DATARECORD << "A_B0R14[" << ma_tmp <<"]: "<< A_B0R14[ma_tmp] << "\n";					
			   //if(display == 1)DATARECORD << "A_B0R15[" << ma_tmp <<"]: "<< A_B0R15[ma_tmp] <<"\n";
			}
			else {
               A[index0]     = A_B1R0[ma_tmp];
               A[index1]     = A_B1R1[ma_tmp];
               A[index2]     = A_B1R2[ma_tmp];
               A[index3]     = A_B1R3[ma_tmp];
               A[index4]     = A_B1R4[ma_tmp];
               A[index5]     = A_B1R5[ma_tmp];
               A[index6]     = A_B1R6[ma_tmp];
               A[index7]     = A_B1R7[ma_tmp];
               A[index8]     = A_B1R8[ma_tmp];
               A[index9]     = A_B1R9[ma_tmp];
               A[index10]    = A_B1R10[ma_tmp];
               A[index11]    = A_B1R11[ma_tmp];
               A[index12]    = A_B1R12[ma_tmp];
               A[index13]    = A_B1R13[ma_tmp];
               A[index14]    = A_B1R14[ma_tmp];
               A[index15]    = A_B1R15[ma_tmp];
			   //if(display == 1)DATARECORD << "A_B1R0[" << ma_tmp <<"]: "<< A_B1R0[ma_tmp] << "\n";
			   //if(display == 1)DATARECORD << "A_B1R1[" << ma_tmp <<"]: "<< A_B1R1[ma_tmp] << "\n";
			   //if(display == 1)DATARECORD << "A_B1R2[" << ma_tmp <<"]: "<< A_B1R2[ma_tmp] << "\n";
			   //if(display == 1)DATARECORD << "A_B1R3[" << ma_tmp <<"]: "<< A_B1R3[ma_tmp] << "\n";
			   //if(display == 1)DATARECORD << "A_B1R4[" << ma_tmp <<"]: "<< A_B1R4[ma_tmp] << "\n";
			   //if(display == 1)DATARECORD << "A_B1R5[" << ma_tmp <<"]: "<< A_B1R5[ma_tmp] << "\n";
			   //if(display == 1)DATARECORD << "A_B1R6[" << ma_tmp <<"]: "<< A_B1R6[ma_tmp] << "\n";
			   //if(display == 1)DATARECORD << "A_B1R7[" << ma_tmp <<"]: "<< A_B1R7[ma_tmp] << "\n";
			   //if(display == 1)DATARECORD << "A_B1R8[" << ma_tmp <<"]: "<< A_B1R8[ma_tmp] << "\n";
			   //if(display == 1)DATARECORD << "A_B1R9[" << ma_tmp <<"]: "<< A_B1R9[ma_tmp] << "\n";
			   //if(display == 1)DATARECORD << "A_B1R10[" << ma_tmp <<"]: "<< A_B1R10[ma_tmp] << "\n";					
			   //if(display == 1)DATARECORD << "A_B1R11[" << ma_tmp <<"]: "<< A_B1R11[ma_tmp] << "\n";					
			   //if(display == 1)DATARECORD << "A_B1R12[" << ma_tmp <<"]: "<< A_B1R12[ma_tmp] << "\n";					
			   //if(display == 1)DATARECORD << "A_B1R13[" << ma_tmp <<"]: "<< A_B1R13[ma_tmp] << "\n";					
			   //if(display == 1)DATARECORD << "A_B1R14[" << ma_tmp <<"]: "<< A_B1R14[ma_tmp] << "\n";					
			   //if(display == 1)DATARECORD << "A_B1R15[" << ma_tmp <<"]: "<< A_B1R15[ma_tmp] <<"\n";			   
			}			
		}		
	}
	
	std::ofstream BEFORE_DATA_RELOCATION("./NTT_R16_R8_Before_DATA_RELOCATION.txt");
	std::ofstream AFTER_DATA_RELOCATION("./NTT_R16_R8_After_DATA_RELOCATION.txt");
	
	BEFORE_DATA_RELOCATION <<"***************************************************************\n";
	BEFORE_DATA_RELOCATION <<"DATA Relocation!!\n";
	BEFORE_DATA_RELOCATION <<"***************************************************************\n";
	
	BEFORE_DATA_RELOCATION <<"***************************************************************\n";
	BEFORE_DATA_RELOCATION << "Relocation right shift bits : " << right_shift_bits << "\n";
	BEFORE_DATA_RELOCATION <<"***************************************************************\n";
    //data relocation for INTT input
	
	for(int i = 0; i < group;i++){
		for(int j = 0; j < radix; j++){
			gray_i  = Gray(i,group);
        	BC_tmp  = j * group + gray_i;
			BEFORE_DATA_RELOCATION <<"---------------------------------------------------\n";
			BEFORE_DATA_RELOCATION <<" BC_tmp: "<< BC_tmp <<"\n";
			RR_R16_R8(BC_tmp,right_shift_bits,BC);
			BEFORE_DATA_RELOCATION <<" BC: "<< BC <<"\n";
            AGU_R16(BC,bn_tmp,ma_tmp);
            BEFORE_DATA_RELOCATION <<" bn_tmp:  "<< bn_tmp <<"\n";
			BEFORE_DATA_RELOCATION <<" ma_tmp:  "<< ma_tmp <<"\n";
			if(bn_tmp == 0){
				BEFORE_DATA_RELOCATION <<"**********************************************************\n";
		        BEFORE_DATA_RELOCATION <<"A_B0R0["<<ma_tmp<<"]:  " << A_B0R0[ma_tmp]<<"\n";
		        BEFORE_DATA_RELOCATION <<"A_B0R1["<<ma_tmp<<"]:  " << A_B0R1[ma_tmp]<<"\n";
		        BEFORE_DATA_RELOCATION <<"A_B0R2["<<ma_tmp<<"]:  " << A_B0R2[ma_tmp]<<"\n";
		        BEFORE_DATA_RELOCATION <<"A_B0R3["<<ma_tmp<<"]:  " << A_B0R3[ma_tmp]<<"\n";
		        BEFORE_DATA_RELOCATION <<"A_B0R4["<<ma_tmp<<"]:  " << A_B0R4[ma_tmp]<<"\n";
		        BEFORE_DATA_RELOCATION <<"A_B0R5["<<ma_tmp<<"]:  " << A_B0R5[ma_tmp]<<"\n";
		        BEFORE_DATA_RELOCATION <<"A_B0R6["<<ma_tmp<<"]:  " << A_B0R6[ma_tmp]<<"\n";
		        BEFORE_DATA_RELOCATION <<"A_B0R7["<<ma_tmp<<"]:  " << A_B0R7[ma_tmp]<<"\n";
		        BEFORE_DATA_RELOCATION <<"A_B0R8["<<ma_tmp<<"]:  " << A_B0R8[ma_tmp]<<"\n";
		        BEFORE_DATA_RELOCATION <<"A_B0R9["<<ma_tmp<<"]:  " << A_B0R9[ma_tmp]<<"\n";
		        BEFORE_DATA_RELOCATION <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";
		        BEFORE_DATA_RELOCATION <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";
		        BEFORE_DATA_RELOCATION <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";
		        BEFORE_DATA_RELOCATION <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";
		        BEFORE_DATA_RELOCATION <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";
		        BEFORE_DATA_RELOCATION <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";
				//-----------------------------------------------------------------------------
                data_tmp_1   = A_B0R1[ma_tmp];
                data_tmp_2   = A_B0R2[ma_tmp];
                data_tmp_3   = A_B0R3[ma_tmp];
                data_tmp_4   = A_B0R4[ma_tmp];
                data_tmp_5   = A_B0R5[ma_tmp];
                data_tmp_6   = A_B0R6[ma_tmp];
                data_tmp_7   = A_B0R7[ma_tmp];
				data_tmp_8   = A_B0R8[ma_tmp];
				data_tmp_9   = A_B0R9[ma_tmp];
				data_tmp_10  = A_B0R10[ma_tmp];
				data_tmp_11  = A_B0R11[ma_tmp];
				data_tmp_12  = A_B0R12[ma_tmp];
				data_tmp_13  = A_B0R13[ma_tmp];
				data_tmp_14  = A_B0R14[ma_tmp];
				A_B0R1[ma_tmp]  = data_tmp_8;
				A_B0R2[ma_tmp]  = data_tmp_1;
				A_B0R3[ma_tmp]  = data_tmp_9;
				A_B0R4[ma_tmp]  = data_tmp_2;
				A_B0R5[ma_tmp]  = data_tmp_10;
				A_B0R6[ma_tmp]  = data_tmp_3;
				A_B0R7[ma_tmp]  = data_tmp_11;
				A_B0R8[ma_tmp]  = data_tmp_4;
				A_B0R9[ma_tmp]  = data_tmp_12;
				A_B0R10[ma_tmp] = data_tmp_5;
				A_B0R11[ma_tmp] = data_tmp_13;
				A_B0R12[ma_tmp] = data_tmp_6;
                A_B0R13[ma_tmp] = data_tmp_14;
				A_B0R14[ma_tmp] = data_tmp_7;
				
                //------------------------------------------------------------------------------
				AFTER_DATA_RELOCATION <<"**********************************************************\n";
		        AFTER_DATA_RELOCATION <<"A_B0R0["<<ma_tmp<<"]:  " << A_B0R0[ma_tmp]<<"\n";
		        AFTER_DATA_RELOCATION <<"A_B0R1["<<ma_tmp<<"]:  " << A_B0R1[ma_tmp]<<"\n";
		        AFTER_DATA_RELOCATION <<"A_B0R2["<<ma_tmp<<"]:  " << A_B0R2[ma_tmp]<<"\n";
		        AFTER_DATA_RELOCATION <<"A_B0R3["<<ma_tmp<<"]:  " << A_B0R3[ma_tmp]<<"\n";
		        AFTER_DATA_RELOCATION <<"A_B0R4["<<ma_tmp<<"]:  " << A_B0R4[ma_tmp]<<"\n";
		        AFTER_DATA_RELOCATION <<"A_B0R5["<<ma_tmp<<"]:  " << A_B0R5[ma_tmp]<<"\n";
		        AFTER_DATA_RELOCATION <<"A_B0R6["<<ma_tmp<<"]:  " << A_B0R6[ma_tmp]<<"\n";
		        AFTER_DATA_RELOCATION <<"A_B0R7["<<ma_tmp<<"]:  " << A_B0R7[ma_tmp]<<"\n";
		        AFTER_DATA_RELOCATION <<"A_B0R8["<<ma_tmp<<"]:  " << A_B0R8[ma_tmp]<<"\n";
		        AFTER_DATA_RELOCATION <<"A_B0R9["<<ma_tmp<<"]:  " << A_B0R9[ma_tmp]<<"\n";
		        AFTER_DATA_RELOCATION <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";
		        AFTER_DATA_RELOCATION <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";
		        AFTER_DATA_RELOCATION <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";
		        AFTER_DATA_RELOCATION <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";
		        AFTER_DATA_RELOCATION <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";
		        AFTER_DATA_RELOCATION <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";	
			}else{
		        BEFORE_DATA_RELOCATION <<"*********************************************************\n";
		        BEFORE_DATA_RELOCATION <<"A_B1R0["<<ma_tmp<<"]:  " << A_B1R0[ma_tmp]<<"\n";
		        BEFORE_DATA_RELOCATION <<"A_B1R1["<<ma_tmp<<"]:  " << A_B1R1[ma_tmp]<<"\n";
		        BEFORE_DATA_RELOCATION <<"A_B1R2["<<ma_tmp<<"]:  " << A_B1R2[ma_tmp]<<"\n";
		        BEFORE_DATA_RELOCATION <<"A_B1R3["<<ma_tmp<<"]:  " << A_B1R3[ma_tmp]<<"\n";
		        BEFORE_DATA_RELOCATION <<"A_B1R4["<<ma_tmp<<"]:  " << A_B1R4[ma_tmp]<<"\n";
		        BEFORE_DATA_RELOCATION <<"A_B1R5["<<ma_tmp<<"]:  " << A_B1R5[ma_tmp]<<"\n";
		        BEFORE_DATA_RELOCATION <<"A_B1R6["<<ma_tmp<<"]:  " << A_B1R6[ma_tmp]<<"\n";
		        BEFORE_DATA_RELOCATION <<"A_B1R7["<<ma_tmp<<"]:  " << A_B1R7[ma_tmp]<<"\n";
		        BEFORE_DATA_RELOCATION <<"A_B1R8["<<ma_tmp<<"]:  " << A_B1R8[ma_tmp]<<"\n";
		        BEFORE_DATA_RELOCATION <<"A_B1R9["<<ma_tmp<<"]:  " << A_B1R9[ma_tmp]<<"\n";
		        BEFORE_DATA_RELOCATION <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";
		        BEFORE_DATA_RELOCATION <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";
		        BEFORE_DATA_RELOCATION <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";
		        BEFORE_DATA_RELOCATION <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";
		        BEFORE_DATA_RELOCATION <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";
		        BEFORE_DATA_RELOCATION <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";	
				//---------------------------------------------------------------------------
                data_tmp_1   = A_B1R1[ma_tmp];
                data_tmp_2   = A_B1R2[ma_tmp];
                data_tmp_3   = A_B1R3[ma_tmp];
                data_tmp_4   = A_B1R4[ma_tmp];
                data_tmp_5   = A_B1R5[ma_tmp];
                data_tmp_6   = A_B1R6[ma_tmp];
                data_tmp_7   = A_B1R7[ma_tmp];
				data_tmp_8   = A_B1R8[ma_tmp];
				data_tmp_9   = A_B1R9[ma_tmp];
				data_tmp_10  = A_B1R10[ma_tmp];
				data_tmp_11  = A_B1R11[ma_tmp];
				data_tmp_12  = A_B1R12[ma_tmp];
				data_tmp_13  = A_B1R13[ma_tmp];
				data_tmp_14  = A_B1R14[ma_tmp];
				A_B1R1[ma_tmp]  = data_tmp_8;
				A_B1R2[ma_tmp]  = data_tmp_1;
				A_B1R3[ma_tmp]  = data_tmp_9;
				A_B1R4[ma_tmp]  = data_tmp_2;
				A_B1R5[ma_tmp]  = data_tmp_10;
				A_B1R6[ma_tmp]  = data_tmp_3;
				A_B1R7[ma_tmp]  = data_tmp_11;
				A_B1R8[ma_tmp]  = data_tmp_4;
				A_B1R9[ma_tmp]  = data_tmp_12;
				A_B1R10[ma_tmp] = data_tmp_5;
				A_B1R11[ma_tmp] = data_tmp_13;
				A_B1R12[ma_tmp] = data_tmp_6;
                A_B1R13[ma_tmp] = data_tmp_14;
				A_B1R14[ma_tmp] = data_tmp_7;
								
				//--------------------------------------------------------------------------
		        AFTER_DATA_RELOCATION <<"*********************************************************\n";
		        AFTER_DATA_RELOCATION <<"A_B1R0["<<ma_tmp<<"]:  " << A_B1R0[ma_tmp]<<"\n";
		        AFTER_DATA_RELOCATION <<"A_B1R1["<<ma_tmp<<"]:  " << A_B1R1[ma_tmp]<<"\n";
		        AFTER_DATA_RELOCATION <<"A_B1R2["<<ma_tmp<<"]:  " << A_B1R2[ma_tmp]<<"\n";
		        AFTER_DATA_RELOCATION <<"A_B1R3["<<ma_tmp<<"]:  " << A_B1R3[ma_tmp]<<"\n";
		        AFTER_DATA_RELOCATION <<"A_B1R4["<<ma_tmp<<"]:  " << A_B1R4[ma_tmp]<<"\n";
		        AFTER_DATA_RELOCATION <<"A_B1R5["<<ma_tmp<<"]:  " << A_B1R5[ma_tmp]<<"\n";
		        AFTER_DATA_RELOCATION <<"A_B1R6["<<ma_tmp<<"]:  " << A_B1R6[ma_tmp]<<"\n";
		        AFTER_DATA_RELOCATION <<"A_B1R7["<<ma_tmp<<"]:  " << A_B1R7[ma_tmp]<<"\n";
		        AFTER_DATA_RELOCATION <<"A_B1R8["<<ma_tmp<<"]:  " << A_B1R8[ma_tmp]<<"\n";
		        AFTER_DATA_RELOCATION <<"A_B1R9["<<ma_tmp<<"]:  " << A_B1R9[ma_tmp]<<"\n";
		        AFTER_DATA_RELOCATION <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";
		        AFTER_DATA_RELOCATION <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";
		        AFTER_DATA_RELOCATION <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";
		        AFTER_DATA_RELOCATION <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";
		        AFTER_DATA_RELOCATION <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";
		        AFTER_DATA_RELOCATION <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";					
			}
		
		}
	}
    
	for(int ss = 0; ss < word_size; ss++){
		//bn0
		B0R0[ss]  = A_B0R0[ss];
		B0R1[ss]  = A_B0R1[ss];
		B0R2[ss]  = A_B0R2[ss];
		B0R3[ss]  = A_B0R3[ss];
		B0R4[ss]  = A_B0R4[ss];
		B0R5[ss]  = A_B0R5[ss];
		B0R6[ss]  = A_B0R6[ss];
		B0R7[ss]  = A_B0R7[ss];
		B0R8[ss]  = A_B0R8[ss];
		B0R9[ss]  = A_B0R9[ss];
		B0R10[ss] = A_B0R10[ss];
		B0R11[ss] = A_B0R11[ss];
		B0R12[ss] = A_B0R12[ss];
		B0R13[ss] = A_B0R13[ss];
		B0R14[ss] = A_B0R14[ss];
		B0R15[ss] = A_B0R15[ss];
		//----------------------------
		//bn1
		B1R0[ss]  = A_B1R0[ss];
		B1R1[ss]  = A_B1R1[ss];
		B1R2[ss]  = A_B1R2[ss];
		B1R3[ss]  = A_B1R3[ss];
		B1R4[ss]  = A_B1R4[ss];
		B1R5[ss]  = A_B1R5[ss];
		B1R6[ss]  = A_B1R6[ss];
		B1R7[ss]  = A_B1R7[ss];
		B1R8[ss]  = A_B1R8[ss];
		B1R9[ss]  = A_B1R9[ss];
		B1R10[ss] = A_B1R10[ss];
		B1R11[ss] = A_B1R11[ss];
		B1R12[ss] = A_B1R12[ss];
		B1R13[ss] = A_B1R13[ss];
		B1R14[ss] = A_B1R14[ss];
		B1R15[ss] = A_B1R15[ss];		
	}	
}
void NTTSPMB::INTT_r16_r8(std::vector<ZZ> &A,std::vector<ZZ> &B0R0,std::vector<ZZ> &B0R1,std::vector<ZZ> &B0R2,std::vector<ZZ> &B0R3,
std::vector<ZZ> &B0R4,std::vector<ZZ> &B0R5,std::vector<ZZ> &B0R6,std::vector<ZZ> &B0R7,std::vector<ZZ> &B0R8,std::vector<ZZ> &B0R9,
std::vector<ZZ> &B0R10,std::vector<ZZ> &B0R11,std::vector<ZZ> &B0R12,std::vector<ZZ> &B0R13,std::vector<ZZ> &B0R14,std::vector<ZZ> &B0R15,
std::vector<ZZ> &B1R0,std::vector<ZZ> &B1R1,std::vector<ZZ> &B1R2,std::vector<ZZ> &B1R3,std::vector<ZZ> &B1R4,std::vector<ZZ> &B1R5,
std::vector<ZZ> &B1R6,std::vector<ZZ> &B1R7,std::vector<ZZ> &B1R8,std::vector<ZZ> &B1R9,std::vector<ZZ> &B1R10,std::vector<ZZ> &B1R11,
std::vector<ZZ> &B1R12,std::vector<ZZ> &B1R13,std::vector<ZZ> &B1R14,std::vector<ZZ> &B1R15){
	
	int            Stage;//stage	
	int            word_size;
	int            offset;
	int            BC_RR;
	int            BC_tmp;
	int            BC_Reorder;
    int            group;
    int            ma_tmp;
	int            bn_tmp;
    int            bit_tmp;
	int            bn0_bc_tmp;//frist in bc data
	int            bn1_bc_tmp;//frist in bc data
	int            bn0_ma_reg1;
	int            bn0_ma_reg2;
	int            bn0_ma_reg3;
	int            bn0_ma_reg4;
	int            bn0_ma_reg5;
	int            bn0_ma_reg6;
	int            bn0_ma_reg7;
	int            bn0_ma_reg8;
	int            bn1_ma_reg1;
	int            bn1_ma_reg2;
	int            bn1_ma_reg3;
	int            bn1_ma_reg4;
	int            bn1_ma_reg5;
	int            bn1_ma_reg6;
	int            bn1_ma_reg7;
	int            bn1_ma_reg8;
	int            gray_i;
    int            BC_WIDTH;
	int            tw_modulus;
    int            tw_modulus_tmp;
	//display parameter
	int            display;
    double         Stage_double;	
    std::vector<int> bit_array_tmp;
	
	std::ofstream DATARECORD("./INTT_R16_R8_SPMB.txt");
	
	//-------------------------------------
	display = 1;
	//-------------------------------------

	//radix-16 Stage
    Stage_double  = log2(N);
    Stage       =  (int)floor(Stage_double/4);
	BC_WIDTH    =  (int)ceil(log2(N/16));
	offset      =  (int)N /  16;
	word_size   =  (int)N / (2 * 16);
	group       =  (int)N / (256);
	tw_modulus  =  (int)N /  16;
	bit_array_tmp.resize(BC_WIDTH);

    DATARECORD << "group: "    << group << "\n";
    DATARECORD << "BC_WIDTH: " << BC_WIDTH << "\n";
	
	ZZ               data_tmp_1;
	ZZ               data_tmp_2;
	ZZ               data_tmp_3;
	ZZ               data_tmp_4;
	ZZ               data_tmp_5;
	ZZ               data_tmp_6;
	ZZ               data_tmp_7;
	ZZ               data_tmp_8;
	ZZ               data_tmp_9;
	ZZ               data_tmp_10;
	ZZ               data_tmp_11;
	ZZ               data_tmp_12;
	ZZ               data_tmp_13;
	ZZ               data_tmp_14;
	ZZ               data_tmp_15;
	//data tmp for index 8 of each array
	//16 * 8 matrix
    std::vector<ZZ>  data_tmp_array_1;  
    std::vector<ZZ>  data_tmp_array_2;  
    std::vector<ZZ>  data_tmp_array_3;  
    std::vector<ZZ>  data_tmp_array_4;  
    std::vector<ZZ>  data_tmp_array_5;  
    std::vector<ZZ>  data_tmp_array_6;  
    std::vector<ZZ>  data_tmp_array_7;  
    std::vector<ZZ>  data_tmp_array_8;

    data_tmp_array_1.resize(16);
    data_tmp_array_2.resize(16);
    data_tmp_array_3.resize(16);
    data_tmp_array_4.resize(16);
    data_tmp_array_5.resize(16);
    data_tmp_array_6.resize(16);
	data_tmp_array_7.resize(16);
	data_tmp_array_8.resize(16);
	//---------------------------------
	std::vector<ZZ>  A_B0R0;
	std::vector<ZZ>  A_B0R1;
	std::vector<ZZ>  A_B0R2;
	std::vector<ZZ>  A_B0R3;
	std::vector<ZZ>  A_B0R4;
	std::vector<ZZ>  A_B0R5;
	std::vector<ZZ>  A_B0R6;
	std::vector<ZZ>  A_B0R7;
	std::vector<ZZ>  A_B0R8;
	std::vector<ZZ>  A_B0R9;
	std::vector<ZZ>  A_B0R10;
	std::vector<ZZ>  A_B0R11;
	std::vector<ZZ>  A_B0R12;
	std::vector<ZZ>  A_B0R13;
	std::vector<ZZ>  A_B0R14;
	std::vector<ZZ>  A_B0R15;
	std::vector<ZZ>  A_B1R0;
	std::vector<ZZ>  A_B1R1;
	std::vector<ZZ>  A_B1R2;
	std::vector<ZZ>  A_B1R3;
	std::vector<ZZ>  A_B1R4;
	std::vector<ZZ>  A_B1R5;
	std::vector<ZZ>  A_B1R6;
	std::vector<ZZ>  A_B1R7;
	std::vector<ZZ>  A_B1R8;
	std::vector<ZZ>  A_B1R9;
	std::vector<ZZ>  A_B1R10;
	std::vector<ZZ>  A_B1R11;
	std::vector<ZZ>  A_B1R12;
	std::vector<ZZ>  A_B1R13;
	std::vector<ZZ>  A_B1R14;
	std::vector<ZZ>  A_B1R15;
	
	A_B0R0.resize(word_size);
	A_B0R1.resize(word_size);
	A_B0R2.resize(word_size);
	A_B0R3.resize(word_size);
	A_B0R4.resize(word_size);
	A_B0R5.resize(word_size);
	A_B0R6.resize(word_size);
	A_B0R7.resize(word_size);
	A_B0R8.resize(word_size);
	A_B0R9.resize(word_size);
	A_B0R10.resize(word_size);
	A_B0R11.resize(word_size);
	A_B0R12.resize(word_size);
	A_B0R13.resize(word_size);
	A_B0R14.resize(word_size);
	A_B0R15.resize(word_size);
	A_B1R0.resize(word_size);
	A_B1R1.resize(word_size);
	A_B1R2.resize(word_size);
	A_B1R3.resize(word_size);
	A_B1R4.resize(word_size);
	A_B1R5.resize(word_size);
	A_B1R6.resize(word_size);
	A_B1R7.resize(word_size);
	A_B1R8.resize(word_size);
	A_B1R9.resize(word_size);
	A_B1R10.resize(word_size);
	A_B1R11.resize(word_size);
	A_B1R12.resize(word_size);
	A_B1R13.resize(word_size);
	A_B1R14.resize(word_size);
	A_B1R15.resize(word_size);
	//----------------------------------------------------
	int length;
	ZZ  factor;   //base factor
	ZZ  factor_t; //acctually mul factor
	ZZ  factor_2t;
	ZZ  factor_3t;
	ZZ  factor_4t;
	ZZ  factor_5t;
	ZZ  factor_6t;
	ZZ  factor_7t;
	ZZ  factor_8t;
	ZZ  factor_9t;
	ZZ  factor_10t;
	ZZ  factor_11t;
	ZZ  factor_12t;
	ZZ  factor_13t;
	ZZ  factor_14t;
	ZZ  factor_15t;
    
	//init load data
    for(int i = 0; i < word_size; i++){ 
		A_B0R0[i]  =  B0R0[i];
        A_B0R1[i]  =  B0R1[i];
        A_B0R2[i]  =  B0R2[i];
        A_B0R3[i]  =  B0R3[i];
        A_B0R4[i]  =  B0R4[i];
        A_B0R5[i]  =  B0R5[i];
        A_B0R6[i]  =  B0R6[i];
        A_B0R7[i]  =  B0R7[i];
        A_B0R8[i]  =  B0R8[i];
        A_B0R9[i]  =  B0R9[i];
        A_B0R10[i] =  B0R10[i];
        A_B0R11[i] =  B0R11[i];
        A_B0R12[i] =  B0R12[i];
        A_B0R13[i] =  B0R13[i];
        A_B0R14[i] =  B0R14[i];
        A_B0R15[i] =  B0R15[i];
        A_B1R0[i]  =  B1R0[i];
        A_B1R1[i]  =  B1R1[i];
        A_B1R2[i]  =  B1R2[i];
        A_B1R3[i]  =  B1R3[i];
        A_B1R4[i]  =  B1R4[i];
        A_B1R5[i]  =  B1R5[i];
        A_B1R6[i]  =  B1R6[i];
        A_B1R7[i]  =  B1R7[i];
        A_B1R8[i]  =  B1R8[i];
        A_B1R9[i]  =  B1R9[i];
        A_B1R10[i] =  B1R10[i];
        A_B1R11[i] =  B1R11[i];
        A_B1R12[i] =  B1R12[i];
        A_B1R13[i] =  B1R13[i];
        A_B1R14[i] =  B1R14[i];
        A_B1R15[i] =  B1R15[i];
	}
	
	ma_tmp      = 0;
	bn_tmp      = 0;
	BC_RR       = 0;
	BC_Reorder  = 0;
	
	std::cout << "init load over! \n";
	DATARECORD <<"Stage: " << Stage << "\n";		
	//need modify to mult by twiddle factor
	for(int s = 0; s < Stage; s++){
		if(s == 0)factor = IW;
		else {
			SqrMod(factor,factor,p);
			SqrMod(factor,factor,p);
			SqrMod(factor,factor,p);
			SqrMod(factor,factor,p);
		}
		if(display == 1)DATARECORD <<"***************************************************\n";					
		if(display == 1)DATARECORD <<"NOW stage: " << s << "!!!!\n";					
		tw_modulus_tmp  = tw_modulus >> ( 4 * s);
		for(int i = 0 ;i < group;i++){
			bn0_bc_tmp  = 0;
			bn1_bc_tmp  = 0;
			for(int j = 0;j < radix;j++){
				gray_i  = Gray(i,group);
			    BC_tmp  = j * group + gray_i;
				if(display == 1)DATARECORD <<" -------------------------------------\n";					
				if(display == 1)DATARECORD <<"Now stage: " << s << "!!!!\n";
				if(display == 1)DATARECORD <<"i: " << i << "!!!!\n";
				if(display == 1)DATARECORD <<"j : "<< j <<"\n";
				if(display == 1)DATARECORD <<"gray_i : "<< gray_i <<"\n";
				if(display == 1)DATARECORD <<"BC_tmp: "<< BC_tmp <<"\n";
				if(s == Stage - 1) RR_R16_R8(BC_tmp,(4 * s - 1),BC_RR);
				else RR_R16(BC_tmp,s,BC_RR);
				if(display == 1)DATARECORD <<"BC_RR: "<< BC_RR <<"\n";
				length = BC_RR % tw_modulus_tmp;
				if(display == 1)DATARECORD <<"length: "<< length <<"\n";
				PowerMod(factor_t,factor,length,p);
				BC_IFFT_Reorder_R16_R8(BC_RR,BC_Reorder);
				if(display == 1)DATARECORD <<"BC_Reorder: "<< BC_Reorder <<"\n";
				AGU_R16(BC_Reorder,bn_tmp,ma_tmp);
				if(display == 1)DATARECORD <<"bn_tmp: "<< bn_tmp <<"\n";
				if(display == 1)DATARECORD <<"ma_tmp: "<< ma_tmp <<"\n";
				if(bn_tmp == 0){
					if(j < 2)bn0_bc_tmp = BC_tmp;
                    if(display == 1)DATARECORD <<" Before radix-16 BU operation!! \n";					
                    if(display == 1)DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";						
					PowerMod(factor_2t,factor_t,2,p);
					PowerMod(factor_3t,factor_t,3,p);
					PowerMod(factor_4t,factor_t,4,p);
					PowerMod(factor_5t,factor_t,5,p);
					PowerMod(factor_6t,factor_t,6,p);
					PowerMod(factor_7t,factor_t,7,p);
					PowerMod(factor_8t,factor_t,8,p);
					PowerMod(factor_9t,factor_t,9,p);
					PowerMod(factor_10t,factor_t,10,p);
					PowerMod(factor_11t,factor_t,11,p);
					PowerMod(factor_12t,factor_t,12,p);
					PowerMod(factor_13t,factor_t,13,p);
					PowerMod(factor_14t,factor_t,14,p);
					PowerMod(factor_15t,factor_t,15,p);
					Radix16_BU_INTT(A_B0R0[ma_tmp],A_B0R1[ma_tmp], A_B0R2[ma_tmp], A_B0R3[ma_tmp],A_B0R4[ma_tmp],
							   A_B0R5[ma_tmp],A_B0R6[ma_tmp], A_B0R7[ma_tmp], A_B0R8[ma_tmp],A_B0R9[ma_tmp],
							   A_B0R10[ma_tmp],A_B0R11[ma_tmp],A_B0R12[ma_tmp],A_B0R13[ma_tmp],A_B0R14[ma_tmp],
							   A_B0R15[ma_tmp]);
                    if(display == 1)DATARECORD <<" After radix-16 BU Operation!!!\n";												   
                    if(display == 1)DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";													   
					MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],factor_t,p);
					MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],factor_2t,p);
					MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],factor_3t,p);
					MulMod(A_B0R4[ma_tmp],A_B0R4[ma_tmp],factor_4t,p);
					MulMod(A_B0R5[ma_tmp],A_B0R5[ma_tmp],factor_5t,p);
					MulMod(A_B0R6[ma_tmp],A_B0R6[ma_tmp],factor_6t,p);
					MulMod(A_B0R7[ma_tmp],A_B0R7[ma_tmp],factor_7t,p);
					MulMod(A_B0R8[ma_tmp],A_B0R8[ma_tmp],factor_8t,p);
					MulMod(A_B0R9[ma_tmp],A_B0R9[ma_tmp],factor_9t,p);
					MulMod(A_B0R10[ma_tmp],A_B0R10[ma_tmp],factor_10t,p);
					MulMod(A_B0R11[ma_tmp],A_B0R11[ma_tmp],factor_11t,p);
					MulMod(A_B0R12[ma_tmp],A_B0R12[ma_tmp],factor_12t,p);
					MulMod(A_B0R13[ma_tmp],A_B0R13[ma_tmp],factor_13t,p);
					MulMod(A_B0R14[ma_tmp],A_B0R14[ma_tmp],factor_14t,p);
					MulMod(A_B0R15[ma_tmp],A_B0R15[ma_tmp],factor_15t,p);
                    if(display == 1)DATARECORD <<" After multiplication!!!\n";					
                    if(display == 1)DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";											
					if(j <  2)bn0_ma_reg1 = ma_tmp;
					if((j >= 2)  && (j < 4))bn0_ma_reg2 = ma_tmp;
					if((j >= 4)  && (j < 6))bn0_ma_reg3 = ma_tmp;
					if((j >= 6)  && (j < 8))bn0_ma_reg4 = ma_tmp;
					if((j >= 8)  && (j < 10))bn0_ma_reg5 = ma_tmp;
					if((j >= 10) && (j < 12))bn0_ma_reg6 = ma_tmp;
					if((j >= 12) && (j < 14))bn0_ma_reg7 = ma_tmp;
					if((j >= 14) && (j < 16))bn0_ma_reg8 = ma_tmp;
				}
			    else {
					if(j < 2)bn1_bc_tmp = BC_tmp;
                    if(display == 1)DATARECORD <<" Before radix-16 BU operation!! \n";					
                    if(display == 1)DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";						
					PowerMod(factor_2t,factor_t,2,p);
					PowerMod(factor_3t,factor_t,3,p);
					PowerMod(factor_4t,factor_t,4,p);
					PowerMod(factor_5t,factor_t,5,p);
					PowerMod(factor_6t,factor_t,6,p);
					PowerMod(factor_7t,factor_t,7,p);
					PowerMod(factor_8t,factor_t,8,p);
					PowerMod(factor_9t,factor_t,9,p);
					PowerMod(factor_10t,factor_t,10,p);
					PowerMod(factor_11t,factor_t,11,p);
					PowerMod(factor_12t,factor_t,12,p);
					PowerMod(factor_13t,factor_t,13,p);
					PowerMod(factor_14t,factor_t,14,p);
					PowerMod(factor_15t,factor_t,15,p);					
					Radix16_BU_INTT(A_B1R0[ma_tmp],A_B1R1[ma_tmp], A_B1R2[ma_tmp], A_B1R3[ma_tmp],A_B1R4[ma_tmp],
							   A_B1R5[ma_tmp],A_B1R6[ma_tmp], A_B1R7[ma_tmp], A_B1R8[ma_tmp],A_B1R9[ma_tmp],
							   A_B1R10[ma_tmp],A_B1R11[ma_tmp],A_B1R12[ma_tmp],A_B1R13[ma_tmp],A_B1R14[ma_tmp],
							   A_B1R15[ma_tmp]);
                    if(display == 1)DATARECORD <<" After radix-16 BU Operation!!!\n";												   
                    if(display == 1)DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";								   
					MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],factor_t,p);
					MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],factor_2t,p);
					MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],factor_3t,p);
					MulMod(A_B1R4[ma_tmp],A_B1R4[ma_tmp],factor_4t,p);
					MulMod(A_B1R5[ma_tmp],A_B1R5[ma_tmp],factor_5t,p);
					MulMod(A_B1R6[ma_tmp],A_B1R6[ma_tmp],factor_6t,p);
					MulMod(A_B1R7[ma_tmp],A_B1R7[ma_tmp],factor_7t,p);
					MulMod(A_B1R8[ma_tmp],A_B1R8[ma_tmp],factor_8t,p);
					MulMod(A_B1R9[ma_tmp],A_B1R9[ma_tmp],factor_9t,p);
					MulMod(A_B1R10[ma_tmp],A_B1R10[ma_tmp],factor_10t,p);
					MulMod(A_B1R11[ma_tmp],A_B1R11[ma_tmp],factor_11t,p);
					MulMod(A_B1R12[ma_tmp],A_B1R12[ma_tmp],factor_12t,p);
					MulMod(A_B1R13[ma_tmp],A_B1R13[ma_tmp],factor_13t,p);
					MulMod(A_B1R14[ma_tmp],A_B1R14[ma_tmp],factor_14t,p);
					MulMod(A_B1R15[ma_tmp],A_B1R15[ma_tmp],factor_15t,p);
                    if(display == 1)DATARECORD <<" After multiplication!!!\n";					
                    if(display == 1)DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";					
                    if(display == 1)DATARECORD <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";						
                    if(j <  2)bn1_ma_reg1 = ma_tmp;					
                    if((j >= 2)  && (j < 4))bn1_ma_reg2 = ma_tmp;
                    if((j >= 4)  && (j < 6))bn1_ma_reg3 = ma_tmp;
                    if((j >= 6)  && (j < 8))bn1_ma_reg4 = ma_tmp;
                    if((j >= 8)  && (j < 10))bn1_ma_reg5 = ma_tmp;
                    if((j >= 10) && (j < 12))bn1_ma_reg6 = ma_tmp;
                    if((j >= 12) && (j < 14))bn1_ma_reg7 = ma_tmp;
                    if((j >= 14) && (j < 16))bn1_ma_reg8 = ma_tmp;					
				}
			}
		//data relocation
		 if(s < Stage-1){
		  if(bn1_bc_tmp > bn0_bc_tmp){
			 data_tmp_1  = A_B0R1[bn0_ma_reg1];
			 data_tmp_2  = A_B0R2[bn0_ma_reg1];
			 data_tmp_3  = A_B0R3[bn0_ma_reg1];
			 data_tmp_4  = A_B0R4[bn0_ma_reg1];
			 data_tmp_5  = A_B0R5[bn0_ma_reg1];
			 data_tmp_6  = A_B0R6[bn0_ma_reg1];
			 data_tmp_7  = A_B0R7[bn0_ma_reg1];
			 data_tmp_8  = A_B0R8[bn0_ma_reg1];
			 data_tmp_9  = A_B0R9[bn0_ma_reg1];
			 data_tmp_10 = A_B0R10[bn0_ma_reg1];
			 data_tmp_11 = A_B0R11[bn0_ma_reg1];
			 data_tmp_12 = A_B0R12[bn0_ma_reg1];
			 data_tmp_13 = A_B0R13[bn0_ma_reg1];
			 data_tmp_14 = A_B0R14[bn0_ma_reg1];
			 data_tmp_15 = A_B0R15[bn0_ma_reg1];
			 A_B0R1[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg1];
			 A_B0R2[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg2];
			 A_B0R3[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg2];
			 A_B0R4[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg3];
			 A_B0R5[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg3];
			 A_B0R6[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg1] = A_B0R0[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg1] = A_B1R0[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg1] = A_B0R0[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg1] = A_B1R0[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg1] = A_B1R0[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg1] = A_B0R0[bn0_ma_reg8];
			 A_B1R0[bn1_ma_reg1]  = data_tmp_1; 
			 A_B1R0[bn1_ma_reg2]  = data_tmp_2; 
			 A_B0R0[bn0_ma_reg2]  = data_tmp_3; 
			 A_B1R0[bn1_ma_reg3]  = data_tmp_4; 
			 A_B0R0[bn0_ma_reg3]  = data_tmp_5; 
			 A_B0R0[bn0_ma_reg4]  = data_tmp_6; 
			 A_B1R0[bn1_ma_reg4]  = data_tmp_7; 
			 A_B1R0[bn1_ma_reg5]  = data_tmp_8; 
			 A_B0R0[bn0_ma_reg5]  = data_tmp_9; 
			 A_B0R0[bn0_ma_reg6]  = data_tmp_10;
			 A_B1R0[bn1_ma_reg6]  = data_tmp_11;
			 A_B0R0[bn0_ma_reg7]  = data_tmp_12;
			 A_B1R0[bn1_ma_reg7]  = data_tmp_13;
			 A_B1R0[bn1_ma_reg8]  = data_tmp_14;
			 A_B0R0[bn0_ma_reg8]  = data_tmp_15;
			 /*********************************************/
			 data_tmp_1  = A_B1R2[bn1_ma_reg1];
			 data_tmp_2  = A_B1R3[bn1_ma_reg1];
			 data_tmp_3  = A_B1R4[bn1_ma_reg1];
			 data_tmp_4  = A_B1R5[bn1_ma_reg1];
			 data_tmp_5  = A_B1R6[bn1_ma_reg1];
			 data_tmp_6  = A_B1R7[bn1_ma_reg1];
			 data_tmp_7  = A_B1R8[bn1_ma_reg1];
			 data_tmp_8  = A_B1R9[bn1_ma_reg1];
			 data_tmp_9  = A_B1R10[bn1_ma_reg1];
			 data_tmp_10 = A_B1R11[bn1_ma_reg1];
			 data_tmp_11 = A_B1R12[bn1_ma_reg1];
			 data_tmp_12 = A_B1R13[bn1_ma_reg1];
			 data_tmp_13 = A_B1R14[bn1_ma_reg1];
			 data_tmp_14 = A_B1R15[bn1_ma_reg1];
             A_B1R2[bn1_ma_reg1]  =	 A_B1R1[bn1_ma_reg2];
             A_B1R3[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg2];
             A_B1R4[bn1_ma_reg1]  =  A_B1R1[bn1_ma_reg3];
             A_B1R5[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg3];
             A_B1R6[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg4];
             A_B1R7[bn1_ma_reg1]  =  A_B1R1[bn1_ma_reg4];
             A_B1R8[bn1_ma_reg1]  =  A_B1R1[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg1] =  A_B0R1[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg1] =  A_B1R1[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg1] =  A_B0R1[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg1] =  A_B1R1[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg1] =  A_B1R1[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg1] =  A_B0R1[bn0_ma_reg8];
			 A_B1R1[bn1_ma_reg2]  =  data_tmp_1; 
			 A_B0R1[bn0_ma_reg2]  =  data_tmp_2; 
			 A_B1R1[bn1_ma_reg3]  =  data_tmp_3; 
			 A_B0R1[bn0_ma_reg3]  =  data_tmp_4; 
			 A_B0R1[bn0_ma_reg4]  =  data_tmp_5; 
			 A_B1R1[bn1_ma_reg4]  =  data_tmp_6; 
			 A_B1R1[bn1_ma_reg5]  =  data_tmp_7; 
			 A_B0R1[bn0_ma_reg5]  =  data_tmp_8; 
			 A_B0R1[bn0_ma_reg6]  =  data_tmp_9; 
			 A_B1R1[bn1_ma_reg6]  =  data_tmp_10;
			 A_B0R1[bn0_ma_reg7]  =  data_tmp_11;
			 A_B1R1[bn1_ma_reg7]  =  data_tmp_12;
			 A_B1R1[bn1_ma_reg8]  =  data_tmp_13;
			 A_B0R1[bn0_ma_reg8]  =  data_tmp_14;
			/************************************************************/ 
			 data_tmp_1  =   A_B1R3[bn1_ma_reg2];
			 data_tmp_2  =   A_B1R4[bn1_ma_reg2];
			 data_tmp_3  =   A_B1R5[bn1_ma_reg2];
			 data_tmp_4  =   A_B1R6[bn1_ma_reg2];
			 data_tmp_5  =   A_B1R7[bn1_ma_reg2];
			 data_tmp_6  =   A_B1R8[bn1_ma_reg2];
			 data_tmp_7  =   A_B1R9[bn1_ma_reg2];
			 data_tmp_8  =   A_B1R10[bn1_ma_reg2];
			 data_tmp_9  =   A_B1R11[bn1_ma_reg2];
			 data_tmp_10 =   A_B1R12[bn1_ma_reg2];
			 data_tmp_11 =   A_B1R13[bn1_ma_reg2];
			 data_tmp_12 =   A_B1R14[bn1_ma_reg2];
			 data_tmp_13 =   A_B1R15[bn1_ma_reg2];
			 A_B1R3[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg2];
			 A_B1R4[bn1_ma_reg2]  =  A_B1R2[bn1_ma_reg3];
			 A_B1R5[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg3];
			 A_B1R6[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg4];
			 A_B1R7[bn1_ma_reg2]  =  A_B1R2[bn1_ma_reg4];
			 A_B1R8[bn1_ma_reg2]  =  A_B1R2[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg2] =  A_B0R2[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg2] =  A_B1R2[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg2] =  A_B0R2[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg2] =  A_B1R2[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg2] =  A_B1R2[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg2] =  A_B0R2[bn0_ma_reg8];
			 A_B0R2[bn0_ma_reg2]  =  data_tmp_1; 
			 A_B1R2[bn1_ma_reg3]  =  data_tmp_2; 
			 A_B0R2[bn0_ma_reg3]  =  data_tmp_3; 
			 A_B0R2[bn0_ma_reg4]  =  data_tmp_4; 
			 A_B1R2[bn1_ma_reg4]  =  data_tmp_5; 
			 A_B1R2[bn1_ma_reg5]  =  data_tmp_6; 
			 A_B0R2[bn0_ma_reg5]  =  data_tmp_7; 
			 A_B0R2[bn0_ma_reg6]  =  data_tmp_8; 
			 A_B1R2[bn1_ma_reg6]  =  data_tmp_9; 
			 A_B0R2[bn0_ma_reg7]  =  data_tmp_10;
			 A_B1R2[bn1_ma_reg7]  =  data_tmp_11;
			 A_B1R2[bn1_ma_reg8]  =  data_tmp_12;
			 A_B0R2[bn0_ma_reg8]  =  data_tmp_13;
			 //******************************************************
			 data_tmp_1  = A_B0R4[bn0_ma_reg2];
			 data_tmp_2  = A_B0R5[bn0_ma_reg2];
			 data_tmp_3  = A_B0R6[bn0_ma_reg2];
			 data_tmp_4  = A_B0R7[bn0_ma_reg2];
			 data_tmp_5  = A_B0R8[bn0_ma_reg2];
			 data_tmp_6  = A_B0R9[bn0_ma_reg2];
			 data_tmp_7  = A_B0R10[bn0_ma_reg2];
			 data_tmp_8  = A_B0R11[bn0_ma_reg2];
			 data_tmp_9  = A_B0R12[bn0_ma_reg2];
			 data_tmp_10 = A_B0R13[bn0_ma_reg2];
			 data_tmp_11 = A_B0R14[bn0_ma_reg2];
			 data_tmp_12 = A_B0R15[bn0_ma_reg2];
			 A_B0R4[bn0_ma_reg2]  =  A_B1R3[bn1_ma_reg3];
			 A_B0R5[bn0_ma_reg2]  =  A_B0R3[bn0_ma_reg3];
			 A_B0R6[bn0_ma_reg2]  =  A_B0R3[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg2]  =  A_B1R3[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg2]  =  A_B1R3[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg2]  =  A_B0R3[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg2] =  A_B0R3[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg2] =  A_B1R3[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg2] =  A_B0R3[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg2] =  A_B1R3[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg2] =  A_B1R3[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg2] =  A_B0R3[bn0_ma_reg8];
			 A_B1R3[bn1_ma_reg3]  =  data_tmp_1;  
			 A_B0R3[bn0_ma_reg3]  =  data_tmp_2;  
			 A_B0R3[bn0_ma_reg4]  =  data_tmp_3;  
			 A_B1R3[bn1_ma_reg4]  =  data_tmp_4;  
			 A_B1R3[bn1_ma_reg5]  =  data_tmp_5;  
			 A_B0R3[bn0_ma_reg5]  =  data_tmp_6;  
			 A_B0R3[bn0_ma_reg6]  =  data_tmp_7;  
			 A_B1R3[bn1_ma_reg6]  =  data_tmp_8;  
			 A_B0R3[bn0_ma_reg7]  =  data_tmp_9;  
			 A_B1R3[bn1_ma_reg7]  =  data_tmp_10; 
			 A_B1R3[bn1_ma_reg8]  =  data_tmp_11; 
			 A_B0R3[bn0_ma_reg8]  =  data_tmp_12; 
			 //----------------------------------------------------------------------
             data_tmp_1  = A_B1R5[bn1_ma_reg3];
			 data_tmp_2  = A_B1R6[bn1_ma_reg3];
			 data_tmp_3  = A_B1R7[bn1_ma_reg3];
			 data_tmp_4  = A_B1R8[bn1_ma_reg3];
			 data_tmp_5  = A_B1R9[bn1_ma_reg3];
			 data_tmp_6  = A_B1R10[bn1_ma_reg3];
			 data_tmp_7  = A_B1R11[bn1_ma_reg3];
			 data_tmp_8  = A_B1R12[bn1_ma_reg3];
			 data_tmp_9  = A_B1R13[bn1_ma_reg3];
			 data_tmp_10 = A_B1R14[bn1_ma_reg3];
			 data_tmp_11 = A_B1R15[bn1_ma_reg3];
			 A_B1R5[bn1_ma_reg3]  =  A_B0R4[bn0_ma_reg3];
			 A_B1R6[bn1_ma_reg3]  =  A_B0R4[bn0_ma_reg4];
			 A_B1R7[bn1_ma_reg3]  =  A_B1R4[bn1_ma_reg4];
			 A_B1R8[bn1_ma_reg3]  =  A_B1R4[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg3]  =  A_B0R4[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg3] =  A_B0R4[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg3] =  A_B1R4[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg3] =  A_B0R4[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg3] =  A_B1R4[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg3] =  A_B1R4[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg3] =  A_B0R4[bn0_ma_reg8];
			 A_B0R4[bn0_ma_reg3]  =  data_tmp_1 ;
			 A_B0R4[bn0_ma_reg4]  =  data_tmp_2 ;
			 A_B1R4[bn1_ma_reg4]  =  data_tmp_3 ;
			 A_B1R4[bn1_ma_reg5]  =  data_tmp_4 ;
			 A_B0R4[bn0_ma_reg5]  =  data_tmp_5 ;
			 A_B0R4[bn0_ma_reg6]  =  data_tmp_6 ;
			 A_B1R4[bn1_ma_reg6]  =  data_tmp_7 ;
			 A_B0R4[bn0_ma_reg7]  =  data_tmp_8 ;
			 A_B1R4[bn1_ma_reg7]  =  data_tmp_9 ;
			 A_B1R4[bn1_ma_reg8]  =  data_tmp_10;
			 A_B0R4[bn0_ma_reg8]  =  data_tmp_11;
			 //************************************************************************
			 data_tmp_1  = A_B0R6[bn0_ma_reg3];
			 data_tmp_2  = A_B0R7[bn0_ma_reg3];
			 data_tmp_3  = A_B0R8[bn0_ma_reg3];
			 data_tmp_4  = A_B0R9[bn0_ma_reg3];
			 data_tmp_5  = A_B0R10[bn0_ma_reg3];
			 data_tmp_6  = A_B0R11[bn0_ma_reg3];
			 data_tmp_7  = A_B0R12[bn0_ma_reg3];
			 data_tmp_8  = A_B0R13[bn0_ma_reg3];
			 data_tmp_9  = A_B0R14[bn0_ma_reg3];
			 data_tmp_10 = A_B0R15[bn0_ma_reg3];
			 A_B0R6[bn0_ma_reg3]  =  A_B0R5[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg3]  =  A_B1R5[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg3]  =  A_B1R5[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg3]  =  A_B0R5[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg3] =  A_B0R5[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg3] =  A_B1R5[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg3] =  A_B0R5[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg3] =  A_B1R5[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg3] =  A_B1R5[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg3] =  A_B0R5[bn0_ma_reg8];
			 A_B0R5[bn0_ma_reg4]  =  data_tmp_1; 
			 A_B1R5[bn1_ma_reg4]  =  data_tmp_2; 
			 A_B1R5[bn1_ma_reg5]  =  data_tmp_3; 
			 A_B0R5[bn0_ma_reg5]  =  data_tmp_4; 
			 A_B0R5[bn0_ma_reg6]  =  data_tmp_5; 
			 A_B1R5[bn1_ma_reg6]  =  data_tmp_6; 
			 A_B0R5[bn0_ma_reg7]  =  data_tmp_7; 
			 A_B1R5[bn1_ma_reg7]  =  data_tmp_8; 
			 A_B1R5[bn1_ma_reg8]  =  data_tmp_9; 
			 A_B0R5[bn0_ma_reg8]  =  data_tmp_10;
			 //-----------------------------------------------------------------------
			 data_tmp_1  = A_B0R7[bn0_ma_reg4];
			 data_tmp_2  = A_B0R8[bn0_ma_reg4];
			 data_tmp_3  = A_B0R9[bn0_ma_reg4];
			 data_tmp_4  = A_B0R10[bn0_ma_reg4];
			 data_tmp_5  = A_B0R11[bn0_ma_reg4];
			 data_tmp_6  = A_B0R12[bn0_ma_reg4];
			 data_tmp_7  = A_B0R13[bn0_ma_reg4];
			 data_tmp_8  = A_B0R14[bn0_ma_reg4];
			 data_tmp_9  = A_B0R15[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg4]  =  A_B1R6[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg4]  =  A_B1R6[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg4]  =  A_B0R6[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg4] =  A_B0R6[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg4] =  A_B1R6[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg4] =  A_B0R6[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg4] =  A_B1R6[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg4] =  A_B1R6[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg4] =  A_B0R6[bn0_ma_reg8];
			 A_B1R6[bn1_ma_reg4]  =  data_tmp_1;
			 A_B1R6[bn1_ma_reg5]  =  data_tmp_2;
			 A_B0R6[bn0_ma_reg5]  =  data_tmp_3;
			 A_B0R6[bn0_ma_reg6]  =  data_tmp_4;
			 A_B1R6[bn1_ma_reg6]  =  data_tmp_5;
			 A_B0R6[bn0_ma_reg7]  =  data_tmp_6;
			 A_B1R6[bn1_ma_reg7]  =  data_tmp_7;
			 A_B1R6[bn1_ma_reg8]  =  data_tmp_8;
			 A_B0R6[bn0_ma_reg8]  =  data_tmp_9;
			 //----------------------------------------------------------------------
			 data_tmp_1  = A_B1R8[bn1_ma_reg4];
			 data_tmp_2  = A_B1R9[bn1_ma_reg4];
			 data_tmp_3  = A_B1R10[bn1_ma_reg4];
			 data_tmp_4  = A_B1R11[bn1_ma_reg4];
			 data_tmp_5  = A_B1R12[bn1_ma_reg4];
			 data_tmp_6  = A_B1R13[bn1_ma_reg4];
			 data_tmp_7  = A_B1R14[bn1_ma_reg4];
			 data_tmp_8  = A_B1R15[bn1_ma_reg4];
			 A_B1R8[bn1_ma_reg4]  =  A_B1R7[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg4]  =  A_B0R7[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg4] =  A_B0R7[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg4] =  A_B1R7[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg4] =  A_B0R7[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg4] =  A_B1R7[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg4] =  A_B1R7[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg4] =  A_B0R7[bn0_ma_reg8];
			 A_B1R7[bn1_ma_reg5]  =  data_tmp_1;
			 A_B0R7[bn0_ma_reg5]  =  data_tmp_2;
			 A_B0R7[bn0_ma_reg6]  =  data_tmp_3;
			 A_B1R7[bn1_ma_reg6]  =  data_tmp_4;
			 A_B0R7[bn0_ma_reg7]  =  data_tmp_5;
			 A_B1R7[bn1_ma_reg7]  =  data_tmp_6;
			 A_B1R7[bn1_ma_reg8]  =  data_tmp_7;
			 A_B0R7[bn0_ma_reg8]  =  data_tmp_8;
			 //----------------------------------------------------------------------
			 data_tmp_1 = A_B1R9[bn1_ma_reg5];
			 data_tmp_2 = A_B1R10[bn1_ma_reg5];
			 data_tmp_3 = A_B1R11[bn1_ma_reg5];
			 data_tmp_4 = A_B1R12[bn1_ma_reg5];
			 data_tmp_5 = A_B1R13[bn1_ma_reg5];
			 data_tmp_6 = A_B1R14[bn1_ma_reg5];
			 data_tmp_7 = A_B1R15[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg5]  = A_B0R8[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg5] = A_B0R8[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg5] = A_B1R8[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg5] = A_B0R8[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg5] = A_B1R8[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg5] = A_B1R8[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg5] = A_B0R8[bn0_ma_reg8];
			 A_B0R8[bn0_ma_reg5]  = data_tmp_1;
			 A_B0R8[bn0_ma_reg6]  = data_tmp_2;
			 A_B1R8[bn1_ma_reg6]  = data_tmp_3;
			 A_B0R8[bn0_ma_reg7]  = data_tmp_4;
			 A_B1R8[bn1_ma_reg7]  = data_tmp_5;
			 A_B1R8[bn1_ma_reg8]  = data_tmp_6;
			 A_B0R8[bn0_ma_reg8]  = data_tmp_7;
			 //---------------------------------------------------------------------
			 data_tmp_1  = A_B0R10[bn0_ma_reg5];
			 data_tmp_2  = A_B0R11[bn0_ma_reg5];
			 data_tmp_3  = A_B0R12[bn0_ma_reg5];
			 data_tmp_4  = A_B0R13[bn0_ma_reg5];
			 data_tmp_5  = A_B0R14[bn0_ma_reg5];
			 data_tmp_6  = A_B0R15[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg5] = A_B0R9[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg5] = A_B1R9[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg5] = A_B0R9[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg5] = A_B1R9[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg5] = A_B1R9[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg5] = A_B0R9[bn0_ma_reg8];
			 A_B0R9[bn0_ma_reg6]  = data_tmp_1;
			 A_B1R9[bn1_ma_reg6]  = data_tmp_2;
			 A_B0R9[bn0_ma_reg7]  = data_tmp_3;
			 A_B1R9[bn1_ma_reg7]  = data_tmp_4;
			 A_B1R9[bn1_ma_reg8]  = data_tmp_5;
			 A_B0R9[bn0_ma_reg8]  = data_tmp_6;
			 //--------------------------------------------------------------------
			 data_tmp_1  = A_B0R11[bn0_ma_reg6];
			 data_tmp_2  = A_B0R12[bn0_ma_reg6];
			 data_tmp_3  = A_B0R13[bn0_ma_reg6];
			 data_tmp_4  = A_B0R14[bn0_ma_reg6];
			 data_tmp_5  = A_B0R15[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg6] = A_B1R10[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg6] = A_B0R10[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg6] = A_B1R10[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg6] = A_B1R10[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg6] = A_B0R10[bn0_ma_reg8];
			 A_B1R10[bn1_ma_reg6] = data_tmp_1;
			 A_B0R10[bn0_ma_reg7] = data_tmp_2;
			 A_B1R10[bn1_ma_reg7] = data_tmp_3;
			 A_B1R10[bn1_ma_reg8] = data_tmp_4;
			 A_B0R10[bn0_ma_reg8] = data_tmp_5;
			 //--------------------------------------------------------------------
			 data_tmp_1  = A_B1R12[bn1_ma_reg6];
			 data_tmp_2  = A_B1R13[bn1_ma_reg6];
			 data_tmp_3  = A_B1R14[bn1_ma_reg6];
			 data_tmp_4  = A_B1R15[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg6] = A_B0R11[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg6] = A_B1R11[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg6] = A_B1R11[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg6] = A_B0R11[bn0_ma_reg8];
			 A_B0R11[bn0_ma_reg7] = data_tmp_1;
			 A_B1R11[bn1_ma_reg7] = data_tmp_2;
			 A_B1R11[bn1_ma_reg8] = data_tmp_3;
			 A_B0R11[bn0_ma_reg8] = data_tmp_4;
			 //--------------------------------------------------------------------
			 data_tmp_1 = A_B0R13[bn0_ma_reg7];
			 data_tmp_2 = A_B0R14[bn0_ma_reg7];
			 data_tmp_3 = A_B0R15[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg7] = A_B1R12[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg7] = A_B1R12[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg7] = A_B0R12[bn0_ma_reg8];
			 A_B1R12[bn1_ma_reg7] = data_tmp_1;
			 A_B1R12[bn1_ma_reg8] = data_tmp_2;
			 A_B0R12[bn0_ma_reg8] = data_tmp_3;
			 //--------------------------------------------------------------------
			 data_tmp_1 = A_B1R14[bn1_ma_reg7];
			 data_tmp_2 = A_B1R15[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg7] = A_B1R13[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg7] = A_B0R13[bn0_ma_reg8];
			 A_B1R13[bn1_ma_reg8] = data_tmp_1;
			 A_B0R13[bn0_ma_reg8] = data_tmp_2;
			 //--------------------------------------------------------------------
			 data_tmp_1 = A_B1R15[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg8] = A_B0R14[bn0_ma_reg8];
			 A_B0R14[bn0_ma_reg8] = data_tmp_1; 
		  }else{	 
			 data_tmp_1  = A_B1R1[bn1_ma_reg1];
			 data_tmp_2  = A_B1R2[bn1_ma_reg1];
			 data_tmp_3  = A_B1R3[bn1_ma_reg1];
			 data_tmp_4  = A_B1R4[bn1_ma_reg1];
			 data_tmp_5  = A_B1R5[bn1_ma_reg1];
			 data_tmp_6  = A_B1R6[bn1_ma_reg1];
			 data_tmp_7  = A_B1R7[bn1_ma_reg1];
			 data_tmp_8  = A_B1R8[bn1_ma_reg1];
			 data_tmp_9  = A_B1R9[bn1_ma_reg1];
			 data_tmp_10 = A_B1R10[bn1_ma_reg1];
			 data_tmp_11 = A_B1R11[bn1_ma_reg1];
			 data_tmp_12 = A_B1R12[bn1_ma_reg1];
			 data_tmp_13 = A_B1R13[bn1_ma_reg1];
			 data_tmp_14 = A_B1R14[bn1_ma_reg1];
			 data_tmp_15 = A_B1R15[bn1_ma_reg1];
             A_B1R1[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg1];
             A_B1R2[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg2];
             A_B1R3[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg2];
             A_B1R4[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg3];
             A_B1R5[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg3];
             A_B1R6[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg4];
             A_B1R7[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg4];
             A_B1R8[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg5];
             A_B1R9[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg5];
             A_B1R10[bn1_ma_reg1] =  A_B1R0[bn1_ma_reg6];
             A_B1R11[bn1_ma_reg1] =  A_B0R0[bn0_ma_reg6];
             A_B1R12[bn1_ma_reg1] =  A_B1R0[bn1_ma_reg7];
             A_B1R13[bn1_ma_reg1] =  A_B0R0[bn0_ma_reg7];
             A_B1R14[bn1_ma_reg1] =  A_B0R0[bn0_ma_reg8];
             A_B1R15[bn1_ma_reg1] =  A_B1R0[bn1_ma_reg8];
             A_B0R0[bn0_ma_reg1]  =  data_tmp_1; 
             A_B0R0[bn0_ma_reg2]  =  data_tmp_2; 
             A_B1R0[bn1_ma_reg2]  =  data_tmp_3; 
             A_B0R0[bn0_ma_reg3]  =  data_tmp_4; 
             A_B1R0[bn1_ma_reg3]  =  data_tmp_5; 
             A_B1R0[bn1_ma_reg4]  =  data_tmp_6; 
             A_B0R0[bn0_ma_reg4]  =  data_tmp_7; 
             A_B0R0[bn0_ma_reg5]  =  data_tmp_8; 
             A_B1R0[bn1_ma_reg5]  =  data_tmp_9; 
             A_B1R0[bn1_ma_reg6]  =  data_tmp_10;
             A_B0R0[bn0_ma_reg6]  =  data_tmp_11;
             A_B1R0[bn1_ma_reg7]  =  data_tmp_12;
             A_B0R0[bn0_ma_reg7]  =  data_tmp_13;
             A_B0R0[bn0_ma_reg8]  =  data_tmp_14;
             A_B1R0[bn1_ma_reg8]  =  data_tmp_15;
             //-----------------------------------------------------------
             data_tmp_1  = A_B0R2[bn0_ma_reg1];
			 data_tmp_2  = A_B0R3[bn0_ma_reg1];
			 data_tmp_3  = A_B0R4[bn0_ma_reg1];
			 data_tmp_4  = A_B0R5[bn0_ma_reg1];
			 data_tmp_5  = A_B0R6[bn0_ma_reg1];
			 data_tmp_6  = A_B0R7[bn0_ma_reg1];
			 data_tmp_7  = A_B0R8[bn0_ma_reg1];
			 data_tmp_8  = A_B0R9[bn0_ma_reg1];
			 data_tmp_9  = A_B0R10[bn0_ma_reg1];
			 data_tmp_10 = A_B0R11[bn0_ma_reg1];
			 data_tmp_11 = A_B0R12[bn0_ma_reg1];
			 data_tmp_12 = A_B0R13[bn0_ma_reg1];
			 data_tmp_13 = A_B0R14[bn0_ma_reg1];
			 data_tmp_14 = A_B0R15[bn0_ma_reg1];
			 A_B0R2[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg2];
			 A_B0R3[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg2];
			 A_B0R4[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg3];
			 A_B0R5[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg3];
			 A_B0R6[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg4];
			 A_B0R7[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg1] =  A_B1R1[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg1] =  A_B0R1[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg1] =  A_B1R1[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg1] =  A_B0R1[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg1] =  A_B0R1[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg1] =  A_B1R1[bn1_ma_reg8];
             A_B0R1[bn0_ma_reg2]  =  data_tmp_1; 
             A_B1R1[bn1_ma_reg2]  =  data_tmp_2; 
             A_B0R1[bn0_ma_reg3]  =  data_tmp_3; 
             A_B1R1[bn1_ma_reg3]  =  data_tmp_4; 
             A_B1R1[bn1_ma_reg4]  =  data_tmp_5; 
             A_B0R1[bn0_ma_reg4]  =  data_tmp_6; 
             A_B0R1[bn0_ma_reg5]  =  data_tmp_7; 
             A_B1R1[bn1_ma_reg5]  =  data_tmp_8; 
             A_B1R1[bn1_ma_reg6]  =  data_tmp_9; 
             A_B0R1[bn0_ma_reg6]  =  data_tmp_10;
             A_B1R1[bn1_ma_reg7]  =  data_tmp_11;
             A_B0R1[bn0_ma_reg7]  =  data_tmp_12;
             A_B0R1[bn0_ma_reg8]  =  data_tmp_13;
             A_B1R1[bn1_ma_reg8]  =  data_tmp_14;
             //------------------------------------------------------------
             data_tmp_1  =  A_B0R3[bn0_ma_reg2];
			 data_tmp_2  =  A_B0R4[bn0_ma_reg2];
			 data_tmp_3  =  A_B0R5[bn0_ma_reg2];
			 data_tmp_4  =  A_B0R6[bn0_ma_reg2];
			 data_tmp_5  =  A_B0R7[bn0_ma_reg2];
			 data_tmp_6  =  A_B0R8[bn0_ma_reg2];
			 data_tmp_7  =  A_B0R9[bn0_ma_reg2];
			 data_tmp_8  =  A_B0R10[bn0_ma_reg2];
			 data_tmp_9  =  A_B0R11[bn0_ma_reg2];
			 data_tmp_10 =  A_B0R12[bn0_ma_reg2];
			 data_tmp_11 =  A_B0R13[bn0_ma_reg2];
			 data_tmp_12 =  A_B0R14[bn0_ma_reg2];
			 data_tmp_13 =  A_B0R15[bn0_ma_reg2];
			 A_B0R3[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg2];
			 A_B0R4[bn0_ma_reg2]  =  A_B0R2[bn0_ma_reg3];
			 A_B0R5[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg3];
			 A_B0R6[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg4];
			 A_B0R7[bn0_ma_reg2]  =  A_B0R2[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg2]  =  A_B0R2[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg2] =  A_B1R2[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg2] =  A_B0R2[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg2] =  A_B1R2[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg2] =  A_B0R2[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg2] =  A_B0R2[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg2] =  A_B1R2[bn1_ma_reg8];
			 A_B1R2[bn1_ma_reg2]  =  data_tmp_1; 
			 A_B0R2[bn0_ma_reg3]  =  data_tmp_2; 
			 A_B1R2[bn1_ma_reg3]  =  data_tmp_3; 
			 A_B1R2[bn1_ma_reg4]  =  data_tmp_4; 
			 A_B0R2[bn0_ma_reg4]  =  data_tmp_5; 
			 A_B0R2[bn0_ma_reg5]  =  data_tmp_6; 
			 A_B1R2[bn1_ma_reg5]  =  data_tmp_7; 
			 A_B1R2[bn1_ma_reg6]  =  data_tmp_8; 
			 A_B0R2[bn0_ma_reg6]  =  data_tmp_9; 
			 A_B1R2[bn1_ma_reg7]  =  data_tmp_10;
			 A_B0R2[bn0_ma_reg7]  =  data_tmp_11;
			 A_B0R2[bn0_ma_reg8]  =  data_tmp_12;
			 A_B1R2[bn1_ma_reg8]  =  data_tmp_13;
			 //-----------------------------------------------------------                        
			 data_tmp_1  = A_B1R4[bn1_ma_reg2];
			 data_tmp_2  = A_B1R5[bn1_ma_reg2];
			 data_tmp_3  = A_B1R6[bn1_ma_reg2];
			 data_tmp_4  = A_B1R7[bn1_ma_reg2];
			 data_tmp_5  = A_B1R8[bn1_ma_reg2];
			 data_tmp_6  = A_B1R9[bn1_ma_reg2];
			 data_tmp_7  = A_B1R10[bn1_ma_reg2];
			 data_tmp_8  = A_B1R11[bn1_ma_reg2];
			 data_tmp_9  = A_B1R12[bn1_ma_reg2];
			 data_tmp_10 = A_B1R13[bn1_ma_reg2];
			 data_tmp_11 = A_B1R14[bn1_ma_reg2];
			 data_tmp_12 = A_B1R15[bn1_ma_reg2];
			 A_B1R4[bn1_ma_reg2]  = A_B0R3[bn0_ma_reg3];
			 A_B1R5[bn1_ma_reg2]  = A_B1R3[bn1_ma_reg3];
			 A_B1R6[bn1_ma_reg2]  = A_B1R3[bn1_ma_reg4];
			 A_B1R7[bn1_ma_reg2]  = A_B0R3[bn0_ma_reg4];
			 A_B1R8[bn1_ma_reg2]  = A_B0R3[bn0_ma_reg5];
			 A_B1R9[bn1_ma_reg2]  = A_B1R3[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg2] = A_B1R3[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg2] = A_B0R3[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg2] = A_B1R3[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg2] = A_B0R3[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg2] = A_B0R3[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg2] = A_B1R3[bn1_ma_reg8];
			 A_B0R3[bn0_ma_reg3]  = data_tmp_1; 
			 A_B1R3[bn1_ma_reg3]  = data_tmp_2; 
			 A_B1R3[bn1_ma_reg4]  = data_tmp_3; 
			 A_B0R3[bn0_ma_reg4]  = data_tmp_4; 
			 A_B0R3[bn0_ma_reg5]  = data_tmp_5; 
			 A_B1R3[bn1_ma_reg5]  = data_tmp_6; 
			 A_B1R3[bn1_ma_reg6]  = data_tmp_7; 
			 A_B0R3[bn0_ma_reg6]  = data_tmp_8; 
			 A_B1R3[bn1_ma_reg7]  = data_tmp_9; 
			 A_B0R3[bn0_ma_reg7]  = data_tmp_10;
			 A_B0R3[bn0_ma_reg8]  = data_tmp_11;
			 A_B1R3[bn1_ma_reg8]  = data_tmp_12;
			 //------------------------------------------------------------
			 data_tmp_1  =  A_B0R5[bn0_ma_reg3];
			 data_tmp_2  =  A_B0R6[bn0_ma_reg3];
			 data_tmp_3  =  A_B0R7[bn0_ma_reg3];
			 data_tmp_4  =  A_B0R8[bn0_ma_reg3];
			 data_tmp_5  =  A_B0R9[bn0_ma_reg3];
			 data_tmp_6  =  A_B0R10[bn0_ma_reg3];
			 data_tmp_7  =  A_B0R11[bn0_ma_reg3];
			 data_tmp_8  =  A_B0R12[bn0_ma_reg3];
			 data_tmp_9  =  A_B0R13[bn0_ma_reg3];
			 data_tmp_10 =  A_B0R14[bn0_ma_reg3];
			 data_tmp_11 =  A_B0R15[bn0_ma_reg3];
			 A_B0R5[bn0_ma_reg3]  =  A_B1R4[bn1_ma_reg3];
			 A_B0R6[bn0_ma_reg3]  =  A_B1R4[bn1_ma_reg4];
			 A_B0R7[bn0_ma_reg3]  =  A_B0R4[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg3]  =  A_B0R4[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg3]  =  A_B1R4[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg3] =  A_B1R4[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg3] =  A_B0R4[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg3] =  A_B1R4[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg3] =  A_B0R4[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg3] =  A_B0R4[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg3] =  A_B1R4[bn1_ma_reg8];
			 A_B1R4[bn1_ma_reg3]  =  data_tmp_1; 
			 A_B1R4[bn1_ma_reg4]  =  data_tmp_2; 
			 A_B0R4[bn0_ma_reg4]  =  data_tmp_3; 
			 A_B0R4[bn0_ma_reg5]  =  data_tmp_4; 
			 A_B1R4[bn1_ma_reg5]  =  data_tmp_5; 
			 A_B1R4[bn1_ma_reg6]  =  data_tmp_6; 
			 A_B0R4[bn0_ma_reg6]  =  data_tmp_7; 
			 A_B1R4[bn1_ma_reg7]  =  data_tmp_8; 
			 A_B0R4[bn0_ma_reg7]  =  data_tmp_9; 
			 A_B0R4[bn0_ma_reg8]  =  data_tmp_10;
			 A_B1R4[bn1_ma_reg8]  =  data_tmp_11;
			 //-------------------------------------------------------------
			 data_tmp_1  =  A_B1R6[bn1_ma_reg3];
			 data_tmp_2  =  A_B1R7[bn1_ma_reg3];
			 data_tmp_3  =  A_B1R8[bn1_ma_reg3];
			 data_tmp_4  =  A_B1R9[bn1_ma_reg3];
			 data_tmp_5  =  A_B1R10[bn1_ma_reg3];
			 data_tmp_6  =  A_B1R11[bn1_ma_reg3];
			 data_tmp_7  =  A_B1R12[bn1_ma_reg3];
			 data_tmp_8  =  A_B1R13[bn1_ma_reg3];
			 data_tmp_9  =  A_B1R14[bn1_ma_reg3];
			 data_tmp_10 =  A_B1R15[bn1_ma_reg3];
			 A_B1R6[bn1_ma_reg3]  =  A_B1R5[bn1_ma_reg4];
			 A_B1R7[bn1_ma_reg3]  =  A_B0R5[bn0_ma_reg4];
			 A_B1R8[bn1_ma_reg3]  =  A_B0R5[bn0_ma_reg5];
			 A_B1R9[bn1_ma_reg3]  =  A_B1R5[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg3] =  A_B1R5[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg3] =  A_B0R5[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg3] =  A_B1R5[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg3] =  A_B0R5[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg3] =  A_B0R5[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg3] =  A_B1R5[bn1_ma_reg8];
			 A_B1R5[bn1_ma_reg4]  =  data_tmp_1; 
			 A_B0R5[bn0_ma_reg4]  =  data_tmp_2; 
			 A_B0R5[bn0_ma_reg5]  =  data_tmp_3; 
			 A_B1R5[bn1_ma_reg5]  =  data_tmp_4; 
			 A_B1R5[bn1_ma_reg6]  =  data_tmp_5; 
			 A_B0R5[bn0_ma_reg6]  =  data_tmp_6; 
			 A_B1R5[bn1_ma_reg7]  =  data_tmp_7; 
			 A_B0R5[bn0_ma_reg7]  =  data_tmp_8; 
			 A_B0R5[bn0_ma_reg8]  =  data_tmp_9; 
			 A_B1R5[bn1_ma_reg8]  =  data_tmp_10;
			 //----------------------------------------------------------------
			 data_tmp_1  =  A_B1R7[bn1_ma_reg4];
			 data_tmp_2  =  A_B1R8[bn1_ma_reg4];
			 data_tmp_3  =  A_B1R9[bn1_ma_reg4];
			 data_tmp_4  =  A_B1R10[bn1_ma_reg4];
			 data_tmp_5  =  A_B1R11[bn1_ma_reg4];
			 data_tmp_6  =  A_B1R12[bn1_ma_reg4];
			 data_tmp_7  =  A_B1R13[bn1_ma_reg4];
			 data_tmp_8  =  A_B1R14[bn1_ma_reg4];
			 data_tmp_9  =  A_B1R15[bn1_ma_reg4];
			 A_B1R7[bn1_ma_reg4]  =  A_B0R6[bn0_ma_reg4];
			 A_B1R8[bn1_ma_reg4]  =  A_B0R6[bn0_ma_reg5];
			 A_B1R9[bn1_ma_reg4]  =  A_B1R6[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg4] =  A_B1R6[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg4] =  A_B0R6[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg4] =  A_B1R6[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg4] =  A_B0R6[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg4] =  A_B0R6[bn0_ma_reg8];
             A_B1R15[bn1_ma_reg4] =  A_B1R6[bn1_ma_reg8];
             A_B0R6[bn0_ma_reg4]  =  data_tmp_1;
             A_B0R6[bn0_ma_reg5]  =  data_tmp_2;
             A_B1R6[bn1_ma_reg5]  =  data_tmp_3;
             A_B1R6[bn1_ma_reg6]  =  data_tmp_4;
             A_B0R6[bn0_ma_reg6]  =  data_tmp_5;
             A_B1R6[bn1_ma_reg7]  =  data_tmp_6;
             A_B0R6[bn0_ma_reg7]  =  data_tmp_7;
             A_B0R6[bn0_ma_reg8]  =  data_tmp_8;
             A_B1R6[bn1_ma_reg8]  =  data_tmp_9;
             //------------------------------------------------------------------
			 data_tmp_1  = A_B0R8[bn0_ma_reg4];
			 data_tmp_2  = A_B0R9[bn0_ma_reg4];
			 data_tmp_3  = A_B0R10[bn0_ma_reg4];
			 data_tmp_4  = A_B0R11[bn0_ma_reg4];
			 data_tmp_5  = A_B0R12[bn0_ma_reg4];
			 data_tmp_6  = A_B0R13[bn0_ma_reg4];
			 data_tmp_7  = A_B0R14[bn0_ma_reg4];
			 data_tmp_8  = A_B0R15[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg4]  =  A_B0R7[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg4]  =  A_B1R7[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg4] =  A_B1R7[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg4] =  A_B0R7[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg4] =  A_B1R7[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg4] =  A_B0R7[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg4] =  A_B0R7[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg4] =  A_B1R7[bn1_ma_reg8];
			 A_B0R7[bn0_ma_reg5]  =  data_tmp_1;
			 A_B1R7[bn1_ma_reg5]  =  data_tmp_2;
			 A_B1R7[bn1_ma_reg6]  =  data_tmp_3;
			 A_B0R7[bn0_ma_reg6]  =  data_tmp_4;
			 A_B1R7[bn1_ma_reg7]  =  data_tmp_5;
			 A_B0R7[bn0_ma_reg7]  =  data_tmp_6;
			 A_B0R7[bn0_ma_reg8]  =  data_tmp_7;
			 A_B1R7[bn1_ma_reg8]  =  data_tmp_8;
			 //----------------------------------------------------------------
			 data_tmp_1  = A_B0R9[bn0_ma_reg5];
			 data_tmp_2  = A_B0R10[bn0_ma_reg5];
			 data_tmp_3  = A_B0R11[bn0_ma_reg5];
			 data_tmp_4  = A_B0R12[bn0_ma_reg5];
			 data_tmp_5  = A_B0R13[bn0_ma_reg5];
			 data_tmp_6  = A_B0R14[bn0_ma_reg5];
			 data_tmp_7  = A_B0R15[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg5]  =  A_B1R8[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg8];
			 A_B1R8[bn1_ma_reg5]  =  data_tmp_1;
			 A_B1R8[bn1_ma_reg6]  =  data_tmp_2;
			 A_B0R8[bn0_ma_reg6]  =  data_tmp_3;
			 A_B1R8[bn1_ma_reg7]  =  data_tmp_4;
			 A_B0R8[bn0_ma_reg7]  =  data_tmp_5;
			 A_B0R8[bn0_ma_reg8]  =  data_tmp_6;
			 A_B1R8[bn1_ma_reg8]  =  data_tmp_7;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B1R10[bn1_ma_reg5];
			 data_tmp_2 = A_B1R11[bn1_ma_reg5];
			 data_tmp_3 = A_B1R12[bn1_ma_reg5];
			 data_tmp_4 = A_B1R13[bn1_ma_reg5];
			 data_tmp_5 = A_B1R14[bn1_ma_reg5];
			 data_tmp_6 = A_B1R15[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg5] = A_B1R9[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg5] = A_B0R9[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg5] = A_B1R9[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg5] = A_B0R9[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg5] = A_B0R9[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg5] = A_B1R9[bn1_ma_reg8];
			 A_B1R9[bn1_ma_reg6]  = data_tmp_1;
			 A_B0R9[bn0_ma_reg6]  = data_tmp_2;
			 A_B1R9[bn1_ma_reg7]  = data_tmp_3;
			 A_B0R9[bn0_ma_reg7]  = data_tmp_4;
			 A_B0R9[bn0_ma_reg8]  = data_tmp_5;
			 A_B1R9[bn1_ma_reg8]  = data_tmp_6;
			 //----------------------------------------------------------------
			 data_tmp_1 =  A_B1R11[bn1_ma_reg6];
			 data_tmp_2 =  A_B1R12[bn1_ma_reg6];
			 data_tmp_3 =  A_B1R13[bn1_ma_reg6];
			 data_tmp_4 =  A_B1R14[bn1_ma_reg6];
			 data_tmp_5 =  A_B1R15[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg6]  = A_B0R10[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg6]  = A_B1R10[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg6]  = A_B0R10[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg6]  = A_B0R10[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg6]  = A_B1R10[bn1_ma_reg8];
			 A_B0R10[bn0_ma_reg6]  = data_tmp_1;
			 A_B1R10[bn1_ma_reg7]  = data_tmp_2;
			 A_B0R10[bn0_ma_reg7]  = data_tmp_3;
			 A_B0R10[bn0_ma_reg8]  = data_tmp_4;
			 A_B1R10[bn1_ma_reg8]  = data_tmp_5;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B0R12[bn0_ma_reg6];
			 data_tmp_2 = A_B0R13[bn0_ma_reg6];
			 data_tmp_3 = A_B0R14[bn0_ma_reg6];
			 data_tmp_4 = A_B0R15[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg6] = A_B1R11[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg6] = A_B0R11[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg6] = A_B0R11[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg6] = A_B1R11[bn1_ma_reg8];
			 A_B1R11[bn1_ma_reg7] = data_tmp_1;
			 A_B0R11[bn0_ma_reg7] = data_tmp_2;
			 A_B0R11[bn0_ma_reg8] = data_tmp_3;
			 A_B1R11[bn1_ma_reg8] = data_tmp_4;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B1R13[bn1_ma_reg7];
			 data_tmp_2 = A_B1R14[bn1_ma_reg7];
			 data_tmp_3 = A_B1R15[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg7] = A_B0R12[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg7] = A_B0R12[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg7] = A_B1R12[bn1_ma_reg8];
			 A_B0R12[bn0_ma_reg7] = data_tmp_1;
			 A_B0R12[bn0_ma_reg8] = data_tmp_2;
			 A_B1R12[bn1_ma_reg8] = data_tmp_3;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B0R14[bn0_ma_reg7];
			 data_tmp_2 = A_B0R15[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg7] = A_B0R13[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg7] = A_B1R13[bn1_ma_reg8];
			 A_B0R13[bn0_ma_reg8] = data_tmp_1;
			 A_B1R13[bn1_ma_reg8] = data_tmp_2;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B0R15[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg8] = A_B1R14[bn1_ma_reg8];
			 A_B1R14[bn1_ma_reg8] = data_tmp_1;
			 //-----------------------------------------------------------------
		  }
		 }
		 else{
		  //radix-16 final stage data relocation
		    if(bn1_bc_tmp > bn0_bc_tmp){
			   data_tmp_1  = A_B0R1[bn0_ma_reg1];
			   data_tmp_2  = A_B0R2[bn0_ma_reg1];
			   data_tmp_3  = A_B0R3[bn0_ma_reg1];
			   data_tmp_4  = A_B0R4[bn0_ma_reg1];
			   data_tmp_5  = A_B0R5[bn0_ma_reg1];
			   data_tmp_6  = A_B0R6[bn0_ma_reg1];
			   data_tmp_7  = A_B0R7[bn0_ma_reg1];
			   data_tmp_9  = A_B0R9[bn0_ma_reg1];
			   data_tmp_10 = A_B0R10[bn0_ma_reg1];
			   data_tmp_11 = A_B0R11[bn0_ma_reg1];
			   data_tmp_12 = A_B0R12[bn0_ma_reg1];
			   data_tmp_13 = A_B0R13[bn0_ma_reg1];
			   data_tmp_14 = A_B0R14[bn0_ma_reg1];
			   data_tmp_15 = A_B0R15[bn0_ma_reg1];
			   A_B0R1[bn0_ma_reg1]  =  A_B1R0[bn1_ma_reg1];
			   A_B0R2[bn0_ma_reg1]  =  A_B1R0[bn1_ma_reg2];
			   A_B0R3[bn0_ma_reg1]  =  A_B0R0[bn0_ma_reg2];
			   A_B0R4[bn0_ma_reg1]  =  A_B1R0[bn1_ma_reg3];
			   A_B0R5[bn0_ma_reg1]  =  A_B0R0[bn0_ma_reg3];
			   A_B0R6[bn0_ma_reg1]  =  A_B0R0[bn0_ma_reg4];
			   A_B0R7[bn0_ma_reg1]  =  A_B1R0[bn1_ma_reg4];
			   A_B0R9[bn0_ma_reg1]  =  A_B1R8[bn1_ma_reg1];
			   A_B0R10[bn0_ma_reg1] =  A_B1R8[bn1_ma_reg2];
			   A_B0R11[bn0_ma_reg1] =  A_B0R8[bn0_ma_reg2];
			   A_B0R12[bn0_ma_reg1] =  A_B1R8[bn1_ma_reg3];
			   A_B0R13[bn0_ma_reg1] =  A_B0R8[bn0_ma_reg3];
			   A_B0R14[bn0_ma_reg1] =  A_B0R8[bn0_ma_reg4];
			   A_B0R15[bn0_ma_reg1] =  A_B1R8[bn1_ma_reg4];
			   A_B1R0[bn1_ma_reg1]  =  data_tmp_1;
			   A_B1R0[bn1_ma_reg2]  =  data_tmp_2;
			   A_B0R0[bn0_ma_reg2]  =  data_tmp_3;
			   A_B1R0[bn1_ma_reg3]  =  data_tmp_4;
			   A_B0R0[bn0_ma_reg3]  =  data_tmp_5;
			   A_B0R0[bn0_ma_reg4]  =  data_tmp_6;
			   A_B1R0[bn1_ma_reg4]  =  data_tmp_7;
			   A_B1R8[bn1_ma_reg1]  =  data_tmp_9;
			   A_B1R8[bn1_ma_reg2]  =  data_tmp_10; 
			   A_B0R8[bn0_ma_reg2]  =  data_tmp_11; 
			   A_B1R8[bn1_ma_reg3]  =  data_tmp_12; 
			   A_B0R8[bn0_ma_reg3]  =  data_tmp_13; 
			   A_B0R8[bn0_ma_reg4]  =  data_tmp_14; 
			   A_B1R8[bn1_ma_reg4]  =  data_tmp_15; 
			   //**********************************
			   data_tmp_1  = A_B1R2[bn1_ma_reg1];
			   data_tmp_2  = A_B1R3[bn1_ma_reg1];
			   data_tmp_3  = A_B1R4[bn1_ma_reg1];
			   data_tmp_4  = A_B1R5[bn1_ma_reg1];
			   data_tmp_5  = A_B1R6[bn1_ma_reg1];
			   data_tmp_6  = A_B1R7[bn1_ma_reg1];
			   data_tmp_9  = A_B1R10[bn1_ma_reg1];
			   data_tmp_10 = A_B1R11[bn1_ma_reg1];
			   data_tmp_11 = A_B1R12[bn1_ma_reg1];
			   data_tmp_12 = A_B1R13[bn1_ma_reg1];
			   data_tmp_13 = A_B1R14[bn1_ma_reg1];
			   data_tmp_14 = A_B1R15[bn1_ma_reg1];
			   A_B1R2[bn1_ma_reg1]   =  A_B1R1[bn1_ma_reg2];
			   A_B1R3[bn1_ma_reg1]   =  A_B0R1[bn0_ma_reg2];
			   A_B1R4[bn1_ma_reg1]   =  A_B1R1[bn1_ma_reg3];
			   A_B1R5[bn1_ma_reg1]   =  A_B0R1[bn0_ma_reg3];
			   A_B1R6[bn1_ma_reg1]   =  A_B0R1[bn0_ma_reg4];
			   A_B1R7[bn1_ma_reg1]   =  A_B1R1[bn1_ma_reg4];
			   A_B1R10[bn1_ma_reg1]  =  A_B1R9[bn1_ma_reg2];
			   A_B1R11[bn1_ma_reg1]  =  A_B0R9[bn0_ma_reg2];
			   A_B1R12[bn1_ma_reg1]  =  A_B1R9[bn1_ma_reg3];
			   A_B1R13[bn1_ma_reg1]  =  A_B0R9[bn0_ma_reg3];
			   A_B1R14[bn1_ma_reg1]  =  A_B0R9[bn0_ma_reg4];
			   A_B1R15[bn1_ma_reg1]  =  A_B1R9[bn1_ma_reg4];
			   A_B1R1[bn1_ma_reg2] =  data_tmp_1;
			   A_B0R1[bn0_ma_reg2] =  data_tmp_2;
			   A_B1R1[bn1_ma_reg3] =  data_tmp_3;
			   A_B0R1[bn0_ma_reg3] =  data_tmp_4;
			   A_B0R1[bn0_ma_reg4] =  data_tmp_5;
			   A_B1R1[bn1_ma_reg4] =  data_tmp_6;
			   A_B1R9[bn1_ma_reg2] =  data_tmp_9;
			   A_B0R9[bn0_ma_reg2] =  data_tmp_10; 
			   A_B1R9[bn1_ma_reg3] =  data_tmp_11; 
			   A_B0R9[bn0_ma_reg3] =  data_tmp_12; 
			   A_B0R9[bn0_ma_reg4] =  data_tmp_13; 
			   A_B1R9[bn1_ma_reg4] =  data_tmp_14; 
			   //*************************************
			   data_tmp_1  =  A_B1R3[bn1_ma_reg2];
			   data_tmp_2  =  A_B1R4[bn1_ma_reg2];
			   data_tmp_3  =  A_B1R5[bn1_ma_reg2];
			   data_tmp_4  =  A_B1R6[bn1_ma_reg2];
			   data_tmp_5  =  A_B1R7[bn1_ma_reg2];
			   data_tmp_9  =  A_B1R11[bn1_ma_reg2];
			   data_tmp_10 =  A_B1R12[bn1_ma_reg2];
			   data_tmp_11 =  A_B1R13[bn1_ma_reg2];
			   data_tmp_12 =  A_B1R14[bn1_ma_reg2];
			   data_tmp_13 =  A_B1R15[bn1_ma_reg2];
			   A_B1R3[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg2];
			   A_B1R4[bn1_ma_reg2]  =  A_B1R2[bn1_ma_reg3];
			   A_B1R5[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg3];
			   A_B1R6[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg4];
			   A_B1R7[bn1_ma_reg2]  =  A_B1R2[bn1_ma_reg4];
			   A_B1R11[bn1_ma_reg2] =  A_B0R10[bn0_ma_reg2];
			   A_B1R12[bn1_ma_reg2] =  A_B1R10[bn1_ma_reg3];
			   A_B1R13[bn1_ma_reg2] =  A_B0R10[bn0_ma_reg3];
			   A_B1R14[bn1_ma_reg2] =  A_B0R10[bn0_ma_reg4];
			   A_B1R15[bn1_ma_reg2] =  A_B1R10[bn1_ma_reg4];
			   A_B0R2[bn0_ma_reg2]  =  data_tmp_1;
			   A_B1R2[bn1_ma_reg3]  =  data_tmp_2;
			   A_B0R2[bn0_ma_reg3]  =  data_tmp_3;
			   A_B0R2[bn0_ma_reg4]  =  data_tmp_4;
			   A_B1R2[bn1_ma_reg4]  =  data_tmp_5;
			   A_B0R10[bn0_ma_reg2] =  data_tmp_9;
			   A_B1R10[bn1_ma_reg3] =  data_tmp_10;
			   A_B0R10[bn0_ma_reg3] =  data_tmp_11;
			   A_B0R10[bn0_ma_reg4] =  data_tmp_12;
			   A_B1R10[bn1_ma_reg4] =  data_tmp_13;
			   //************************************
			   data_tmp_1  =  A_B0R4[bn0_ma_reg2];
			   data_tmp_2  =  A_B0R5[bn0_ma_reg2];
			   data_tmp_3  =  A_B0R6[bn0_ma_reg2];
			   data_tmp_4  =  A_B0R7[bn0_ma_reg2];
			   data_tmp_9  =  A_B0R12[bn0_ma_reg2];
			   data_tmp_10 =  A_B0R13[bn0_ma_reg2];
 			   data_tmp_11 =  A_B0R14[bn0_ma_reg2];
			   data_tmp_12 =  A_B0R15[bn0_ma_reg2];
			   A_B0R4[bn0_ma_reg2]  =  A_B1R3[bn1_ma_reg3];
			   A_B0R5[bn0_ma_reg2]  =  A_B0R3[bn0_ma_reg3];
			   A_B0R6[bn0_ma_reg2]  =  A_B0R3[bn0_ma_reg4];
			   A_B0R7[bn0_ma_reg2]  =  A_B1R3[bn1_ma_reg4];
			   A_B0R12[bn0_ma_reg2] =  A_B1R11[bn1_ma_reg3];
			   A_B0R13[bn0_ma_reg2] =  A_B0R11[bn0_ma_reg3];
			   A_B0R14[bn0_ma_reg2] =  A_B0R11[bn0_ma_reg4];
			   A_B0R15[bn0_ma_reg2] =  A_B1R11[bn1_ma_reg4];
			   A_B1R3[bn1_ma_reg3]  =  data_tmp_1;
			   A_B0R3[bn0_ma_reg3]  =  data_tmp_2;
			   A_B0R3[bn0_ma_reg4]  =  data_tmp_3;
			   A_B1R3[bn1_ma_reg4]  =  data_tmp_4;
			   A_B1R11[bn1_ma_reg3] =  data_tmp_9;
			   A_B0R11[bn0_ma_reg3] =  data_tmp_10; 
			   A_B0R11[bn0_ma_reg4] =  data_tmp_11; 
			   A_B1R11[bn1_ma_reg4] =  data_tmp_12; 
			   //***********************************
			   data_tmp_1  = A_B1R5[bn1_ma_reg3];
			   data_tmp_2  = A_B1R6[bn1_ma_reg3];
			   data_tmp_3  = A_B1R7[bn1_ma_reg3];
			   data_tmp_9  = A_B1R13[bn1_ma_reg3];
			   data_tmp_10 = A_B1R14[bn1_ma_reg3];
			   data_tmp_11 = A_B1R15[bn1_ma_reg3];
			   A_B1R5[bn1_ma_reg3]  = A_B0R4[bn0_ma_reg3];
			   A_B1R6[bn1_ma_reg3]  = A_B0R4[bn0_ma_reg4];
			   A_B1R7[bn1_ma_reg3]  = A_B1R4[bn1_ma_reg4];
			   A_B1R13[bn1_ma_reg3] = A_B0R12[bn0_ma_reg3];
			   A_B1R14[bn1_ma_reg3] = A_B0R12[bn0_ma_reg4];
			   A_B1R15[bn1_ma_reg3] = A_B1R12[bn1_ma_reg4];
			   A_B0R4[bn0_ma_reg3]  = data_tmp_1;
			   A_B0R4[bn0_ma_reg4]  = data_tmp_2;
			   A_B1R4[bn1_ma_reg4]  = data_tmp_3;
			   A_B0R12[bn0_ma_reg3] = data_tmp_9;
			   A_B0R12[bn0_ma_reg4] = data_tmp_10;
			   A_B1R12[bn1_ma_reg4] = data_tmp_11;
			   //*********************************
			   data_tmp_1  =  A_B0R6[bn0_ma_reg3];
			   data_tmp_2  =  A_B0R7[bn0_ma_reg3];
			   data_tmp_9  =  A_B0R14[bn0_ma_reg3];
			   data_tmp_10 =  A_B0R15[bn0_ma_reg3];
			   A_B0R6[bn0_ma_reg3]  = A_B0R5[bn0_ma_reg4];
			   A_B0R7[bn0_ma_reg3]  = A_B1R5[bn1_ma_reg4];
			   A_B0R14[bn0_ma_reg3] = A_B0R13[bn0_ma_reg4];
			   A_B0R15[bn0_ma_reg3] = A_B1R13[bn1_ma_reg4];
			   A_B0R5[bn0_ma_reg4]  = data_tmp_1;
			   A_B1R5[bn1_ma_reg4]  = data_tmp_2;
			   A_B0R13[bn0_ma_reg4] = data_tmp_9;
			   A_B1R13[bn1_ma_reg4] = data_tmp_10;
			   //********************************
			   data_tmp_1 = A_B0R7[bn0_ma_reg4];
			   data_tmp_9 = A_B0R15[bn0_ma_reg4];
			   A_B0R7[bn0_ma_reg4]  = A_B1R6[bn1_ma_reg4];
			   A_B0R15[bn0_ma_reg4] = A_B1R14[bn1_ma_reg4];
			   A_B1R6[bn1_ma_reg4]  = data_tmp_1;
			   A_B1R14[bn1_ma_reg4] = data_tmp_9;
			   //-----------------------------------------------------------
			   //***********************************************************
			   //-----------------------------------------------------------
			   data_tmp_1  = A_B1R1[bn1_ma_reg5];
			   data_tmp_2  = A_B1R2[bn1_ma_reg5];
			   data_tmp_3  = A_B1R3[bn1_ma_reg5];
			   data_tmp_4  = A_B1R4[bn1_ma_reg5];
			   data_tmp_5  = A_B1R5[bn1_ma_reg5];
			   data_tmp_6  = A_B1R6[bn1_ma_reg5];
			   data_tmp_7  = A_B1R7[bn1_ma_reg5];
			   data_tmp_9  = A_B1R9[bn1_ma_reg5];
			   data_tmp_10 = A_B1R10[bn1_ma_reg5];
			   data_tmp_11 = A_B1R11[bn1_ma_reg5];
			   data_tmp_12 = A_B1R12[bn1_ma_reg5];
			   data_tmp_13 = A_B1R13[bn1_ma_reg5];
			   data_tmp_14 = A_B1R14[bn1_ma_reg5];
			   data_tmp_15 = A_B1R15[bn1_ma_reg5];
			   A_B1R1[bn1_ma_reg5]  =  A_B0R0[bn0_ma_reg5];
			   A_B1R2[bn1_ma_reg5]  =  A_B0R0[bn0_ma_reg6];
			   A_B1R3[bn1_ma_reg5]  =  A_B1R0[bn1_ma_reg6];
			   A_B1R4[bn1_ma_reg5]  =  A_B0R0[bn0_ma_reg7];
			   A_B1R5[bn1_ma_reg5]  =  A_B1R0[bn1_ma_reg7];
			   A_B1R6[bn1_ma_reg5]  =  A_B1R0[bn1_ma_reg8];
			   A_B1R7[bn1_ma_reg5]  =  A_B0R0[bn0_ma_reg8];
			   A_B1R9[bn1_ma_reg5]  =  A_B0R8[bn0_ma_reg5];
			   A_B1R10[bn1_ma_reg5] =  A_B0R8[bn0_ma_reg6];
			   A_B1R11[bn1_ma_reg5] =  A_B1R8[bn1_ma_reg6];
			   A_B1R12[bn1_ma_reg5] =  A_B0R8[bn0_ma_reg7];
			   A_B1R13[bn1_ma_reg5] =  A_B1R8[bn1_ma_reg7];
			   A_B1R14[bn1_ma_reg5] =  A_B1R8[bn1_ma_reg8];
			   A_B1R15[bn1_ma_reg5] =  A_B0R8[bn0_ma_reg8];
			   A_B0R0[bn0_ma_reg5]  =  data_tmp_1;
			   A_B0R0[bn0_ma_reg6]  =  data_tmp_2;
			   A_B1R0[bn1_ma_reg6]  =  data_tmp_3;
			   A_B0R0[bn0_ma_reg7]  =  data_tmp_4;
			   A_B1R0[bn1_ma_reg7]  =  data_tmp_5;
			   A_B1R0[bn1_ma_reg8]  =  data_tmp_6;
			   A_B0R0[bn0_ma_reg8]  =  data_tmp_7;
			   A_B0R8[bn0_ma_reg5]  =  data_tmp_9;
			   A_B0R8[bn0_ma_reg6]  =  data_tmp_10; 
			   A_B1R8[bn1_ma_reg6]  =  data_tmp_11; 
			   A_B0R8[bn0_ma_reg7]  =  data_tmp_12; 
			   A_B1R8[bn1_ma_reg7]  =  data_tmp_13; 
			   A_B1R8[bn1_ma_reg8]  =  data_tmp_14; 
			   A_B0R8[bn0_ma_reg8]  =  data_tmp_15; 
			   //**********************************
			   data_tmp_1  = A_B0R2[bn0_ma_reg5];
			   data_tmp_2  = A_B0R3[bn0_ma_reg5];
			   data_tmp_3  = A_B0R4[bn0_ma_reg5];
			   data_tmp_4  = A_B0R5[bn0_ma_reg5];
			   data_tmp_5  = A_B0R6[bn0_ma_reg5];
			   data_tmp_6  = A_B0R7[bn0_ma_reg5];
			   data_tmp_9  = A_B0R10[bn0_ma_reg5];
			   data_tmp_10 = A_B0R11[bn0_ma_reg5];
			   data_tmp_11 = A_B0R12[bn0_ma_reg5];
			   data_tmp_12 = A_B0R13[bn0_ma_reg5];
			   data_tmp_13 = A_B0R14[bn0_ma_reg5];
			   data_tmp_14 = A_B0R15[bn0_ma_reg5];
			   A_B0R2[bn0_ma_reg5]  =  A_B0R1[bn0_ma_reg6];
			   A_B0R3[bn0_ma_reg5]  =  A_B1R1[bn1_ma_reg6];
			   A_B0R4[bn0_ma_reg5]  =  A_B0R1[bn0_ma_reg7];
			   A_B0R5[bn0_ma_reg5]  =  A_B1R1[bn1_ma_reg7];
			   A_B0R6[bn0_ma_reg5]  =  A_B1R1[bn1_ma_reg8];
			   A_B0R7[bn0_ma_reg5]  =  A_B0R1[bn0_ma_reg8];
			   A_B0R10[bn0_ma_reg5] =  A_B0R9[bn0_ma_reg6];
			   A_B0R11[bn0_ma_reg5] =  A_B1R9[bn1_ma_reg6];
			   A_B0R12[bn0_ma_reg5] =  A_B0R9[bn0_ma_reg7];
			   A_B0R13[bn0_ma_reg5] =  A_B1R9[bn1_ma_reg7];
			   A_B0R14[bn0_ma_reg5] =  A_B1R9[bn1_ma_reg8];
			   A_B0R15[bn0_ma_reg5] =  A_B0R9[bn0_ma_reg8];
			   A_B0R1[bn0_ma_reg6]  =  data_tmp_1;
			   A_B1R1[bn1_ma_reg6]  =  data_tmp_2;
			   A_B0R1[bn0_ma_reg7]  =  data_tmp_3;
			   A_B1R1[bn1_ma_reg7]  =  data_tmp_4;
			   A_B1R1[bn1_ma_reg8]  =  data_tmp_5;
			   A_B0R1[bn0_ma_reg8]  =  data_tmp_6;
			   A_B0R9[bn0_ma_reg6]  =  data_tmp_9;
			   A_B1R9[bn1_ma_reg6]  =  data_tmp_10;
			   A_B0R9[bn0_ma_reg7]  =  data_tmp_11;
			   A_B1R9[bn1_ma_reg7]  =  data_tmp_12;
			   A_B1R9[bn1_ma_reg8]  =  data_tmp_13;
			   A_B0R9[bn0_ma_reg8]  =  data_tmp_14;
			   //************************************
			   data_tmp_1  =  A_B0R3[bn0_ma_reg6];
			   data_tmp_2  =  A_B0R4[bn0_ma_reg6];
			   data_tmp_3  =  A_B0R5[bn0_ma_reg6];
			   data_tmp_4  =  A_B0R6[bn0_ma_reg6];
			   data_tmp_5  =  A_B0R7[bn0_ma_reg6];
			   data_tmp_9  =  A_B0R11[bn0_ma_reg6];
			   data_tmp_10 =  A_B0R12[bn0_ma_reg6];
			   data_tmp_11 =  A_B0R13[bn0_ma_reg6];
			   data_tmp_12 =  A_B0R14[bn0_ma_reg6];
			   data_tmp_13 =  A_B0R15[bn0_ma_reg6];
			   A_B0R3[bn0_ma_reg6]  =  A_B1R2[bn1_ma_reg6];
			   A_B0R4[bn0_ma_reg6]  =  A_B0R2[bn0_ma_reg7];
			   A_B0R5[bn0_ma_reg6]  =  A_B1R2[bn1_ma_reg7];
			   A_B0R6[bn0_ma_reg6]  =  A_B1R2[bn1_ma_reg8];
			   A_B0R7[bn0_ma_reg6]  =  A_B0R2[bn0_ma_reg8];
			   A_B0R11[bn0_ma_reg6] =  A_B1R10[bn1_ma_reg6];
			   A_B0R12[bn0_ma_reg6] =  A_B0R10[bn0_ma_reg7];
			   A_B0R13[bn0_ma_reg6] =  A_B1R10[bn1_ma_reg7];
			   A_B0R14[bn0_ma_reg6] =  A_B1R10[bn1_ma_reg8];
			   A_B0R15[bn0_ma_reg6] =  A_B0R10[bn0_ma_reg8];
			   A_B1R2[bn1_ma_reg6]  =  data_tmp_1;
			   A_B0R2[bn0_ma_reg7]  =  data_tmp_2;
			   A_B1R2[bn1_ma_reg7]  =  data_tmp_3;
			   A_B1R2[bn1_ma_reg8]  =  data_tmp_4;
			   A_B0R2[bn0_ma_reg8]  =  data_tmp_5;
			   A_B1R10[bn1_ma_reg6]  =  data_tmp_9;
			   A_B0R10[bn0_ma_reg7]  =  data_tmp_10;
			   A_B1R10[bn1_ma_reg7]  =  data_tmp_11;
			   A_B1R10[bn1_ma_reg8]  =  data_tmp_12;
			   A_B0R10[bn0_ma_reg8]  =  data_tmp_13;
			   //***********************************
               data_tmp_1  =  A_B1R4[bn1_ma_reg6];
               data_tmp_2  =  A_B1R5[bn1_ma_reg6];
               data_tmp_3  =  A_B1R6[bn1_ma_reg6];
               data_tmp_4  =  A_B1R7[bn1_ma_reg6];
               data_tmp_9  =  A_B1R12[bn1_ma_reg6];
               data_tmp_10 =  A_B1R13[bn1_ma_reg6];
               data_tmp_11 =  A_B1R14[bn1_ma_reg6];
               data_tmp_12 =  A_B1R15[bn1_ma_reg6];
               A_B1R4[bn1_ma_reg6]  =   A_B0R3[bn0_ma_reg7];
			   A_B1R5[bn1_ma_reg6]  =   A_B1R3[bn1_ma_reg7];
			   A_B1R6[bn1_ma_reg6]  =   A_B1R3[bn1_ma_reg8];
			   A_B1R7[bn1_ma_reg6]  =   A_B0R3[bn0_ma_reg8];
			   A_B1R12[bn1_ma_reg6] =   A_B0R11[bn0_ma_reg7];
               A_B1R13[bn1_ma_reg6] =   A_B1R11[bn1_ma_reg7];
               A_B1R14[bn1_ma_reg6] =   A_B1R11[bn1_ma_reg8];
               A_B1R15[bn1_ma_reg6] =   A_B0R11[bn0_ma_reg8];
               A_B0R3[bn0_ma_reg7]  =   data_tmp_1; 
			   A_B1R3[bn1_ma_reg7]  =   data_tmp_2; 
			   A_B1R3[bn1_ma_reg8]  =   data_tmp_3; 
			   A_B0R3[bn0_ma_reg8]  =   data_tmp_4; 
			   A_B0R11[bn0_ma_reg7] =   data_tmp_9; 
			   A_B1R11[bn1_ma_reg7] =   data_tmp_10;
               A_B1R11[bn1_ma_reg8] =   data_tmp_11;
               A_B0R11[bn0_ma_reg8] =   data_tmp_12;
               //*******************************************
			   data_tmp_1  =  A_B0R5[bn0_ma_reg7];
			   data_tmp_2  =  A_B0R6[bn0_ma_reg7];
			   data_tmp_3  =  A_B0R7[bn0_ma_reg7];
			   data_tmp_9  =  A_B0R13[bn0_ma_reg7];
			   data_tmp_10 =  A_B0R14[bn0_ma_reg7];
			   data_tmp_11 =  A_B0R15[bn0_ma_reg7];
			   A_B0R5[bn0_ma_reg7]  =  A_B1R4[bn1_ma_reg7];
			   A_B0R6[bn0_ma_reg7]  =  A_B1R4[bn1_ma_reg8];
			   A_B0R7[bn0_ma_reg7]  =  A_B0R4[bn0_ma_reg8];
			   A_B0R13[bn0_ma_reg7] =  A_B1R12[bn1_ma_reg7];
			   A_B0R14[bn0_ma_reg7] =  A_B1R12[bn1_ma_reg8];
			   A_B0R15[bn0_ma_reg7] =  A_B0R12[bn0_ma_reg8];
			   A_B1R4[bn1_ma_reg7]  =  data_tmp_1; 
			   A_B1R4[bn1_ma_reg8]  =  data_tmp_2; 
			   A_B0R4[bn0_ma_reg8]  =  data_tmp_3; 
			   A_B1R12[bn1_ma_reg7] =  data_tmp_9; 
			   A_B1R12[bn1_ma_reg8] =  data_tmp_10;
			   A_B0R12[bn0_ma_reg8] =  data_tmp_11;
			   //*******************************************
			   data_tmp_1  =  A_B1R6[bn1_ma_reg7];
			   data_tmp_2  =  A_B1R7[bn1_ma_reg7];
			   data_tmp_9  =  A_B1R14[bn1_ma_reg7];
			   data_tmp_10 =  A_B1R15[bn1_ma_reg7];
			   A_B1R6[bn1_ma_reg7]  = A_B1R5[bn1_ma_reg8];
			   A_B1R7[bn1_ma_reg7]  = A_B0R5[bn0_ma_reg8];
			   A_B1R14[bn1_ma_reg7] = A_B1R13[bn1_ma_reg8];
			   A_B1R15[bn1_ma_reg7] = A_B0R13[bn0_ma_reg8];
               A_B1R5[bn1_ma_reg8]  = data_tmp_1;  
               A_B0R5[bn0_ma_reg8]  = data_tmp_2;  
               A_B1R13[bn1_ma_reg8] = data_tmp_9;  
			   A_B0R13[bn0_ma_reg8] = data_tmp_10;
               //******************************************
               data_tmp_1  =  A_B1R7[bn1_ma_reg8];
			   data_tmp_9  =  A_B1R15[bn1_ma_reg8];
               A_B1R7[bn1_ma_reg8]   =  A_B0R6[bn0_ma_reg8];
			   A_B1R15[bn1_ma_reg8]  =  A_B0R14[bn0_ma_reg8];
               A_B0R6[bn0_ma_reg8]   =  data_tmp_1;
			   A_B0R14[bn0_ma_reg8]  =  data_tmp_9;
			}
            else{
			   data_tmp_1  = A_B1R1[bn1_ma_reg1];
			   data_tmp_2  = A_B1R2[bn1_ma_reg1];
			   data_tmp_3  = A_B1R3[bn1_ma_reg1];
			   data_tmp_4  = A_B1R4[bn1_ma_reg1];
			   data_tmp_5  = A_B1R5[bn1_ma_reg1];
			   data_tmp_6  = A_B1R6[bn1_ma_reg1];
			   data_tmp_7  = A_B1R7[bn1_ma_reg1];
			   data_tmp_9  = A_B1R9[bn1_ma_reg1];
			   data_tmp_10 = A_B1R10[bn1_ma_reg1];
			   data_tmp_11 = A_B1R11[bn1_ma_reg1];
			   data_tmp_12 = A_B1R12[bn1_ma_reg1];
			   data_tmp_13 = A_B1R13[bn1_ma_reg1];
			   data_tmp_14 = A_B1R14[bn1_ma_reg1];
			   data_tmp_15 = A_B1R15[bn1_ma_reg1];
			   A_B1R1[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg1];
			   A_B1R2[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg2];
			   A_B1R3[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg2];
			   A_B1R4[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg3];
			   A_B1R5[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg3];
			   A_B1R6[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg4];
			   A_B1R7[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg4];
			   A_B1R9[bn1_ma_reg1]  =  A_B0R8[bn0_ma_reg1];
			   A_B1R10[bn1_ma_reg1] =  A_B0R8[bn0_ma_reg2];
			   A_B1R11[bn1_ma_reg1] =  A_B1R8[bn1_ma_reg2];
			   A_B1R12[bn1_ma_reg1] =  A_B0R8[bn0_ma_reg3];
			   A_B1R13[bn1_ma_reg1] =  A_B1R8[bn1_ma_reg3];
			   A_B1R14[bn1_ma_reg1] =  A_B1R8[bn1_ma_reg4];
			   A_B1R15[bn1_ma_reg1] =  A_B0R8[bn0_ma_reg4];
			   A_B0R0[bn0_ma_reg1]  =  data_tmp_1;
			   A_B0R0[bn0_ma_reg2]  =  data_tmp_2;
			   A_B1R0[bn1_ma_reg2]  =  data_tmp_3;
			   A_B0R0[bn0_ma_reg3]  =  data_tmp_4;
			   A_B1R0[bn1_ma_reg3]  =  data_tmp_5;
			   A_B1R0[bn1_ma_reg4]  =  data_tmp_6;
			   A_B0R0[bn0_ma_reg4]  =  data_tmp_7;
			   A_B0R8[bn0_ma_reg1]  =  data_tmp_9;
			   A_B0R8[bn0_ma_reg2]  =  data_tmp_10; 
			   A_B1R8[bn1_ma_reg2]  =  data_tmp_11; 
			   A_B0R8[bn0_ma_reg3]  =  data_tmp_12; 
			   A_B1R8[bn1_ma_reg3]  =  data_tmp_13; 
			   A_B1R8[bn1_ma_reg4]  =  data_tmp_14; 
			   A_B0R8[bn0_ma_reg4]  =  data_tmp_15; 
			   //**********************************
			   data_tmp_1  = A_B0R2[bn0_ma_reg1];
			   data_tmp_2  = A_B0R3[bn0_ma_reg1];
			   data_tmp_3  = A_B0R4[bn0_ma_reg1];
			   data_tmp_4  = A_B0R5[bn0_ma_reg1];
			   data_tmp_5  = A_B0R6[bn0_ma_reg1];
			   data_tmp_6  = A_B0R7[bn0_ma_reg1];
			   data_tmp_9  = A_B0R10[bn0_ma_reg1];
			   data_tmp_10 = A_B0R11[bn0_ma_reg1];
			   data_tmp_11 = A_B0R12[bn0_ma_reg1];
			   data_tmp_12 = A_B0R13[bn0_ma_reg1];
			   data_tmp_13 = A_B0R14[bn0_ma_reg1];
			   data_tmp_14 = A_B0R15[bn0_ma_reg1];
			   A_B0R2[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg2];
			   A_B0R3[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg2];
			   A_B0R4[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg3];
			   A_B0R5[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg3];
			   A_B0R6[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg4];
			   A_B0R7[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg4];
			   A_B0R10[bn0_ma_reg1] =  A_B0R9[bn0_ma_reg2];
			   A_B0R11[bn0_ma_reg1] =  A_B1R9[bn1_ma_reg2];
			   A_B0R12[bn0_ma_reg1] =  A_B0R9[bn0_ma_reg3];
			   A_B0R13[bn0_ma_reg1] =  A_B1R9[bn1_ma_reg3];
			   A_B0R14[bn0_ma_reg1] =  A_B1R9[bn1_ma_reg4];
			   A_B0R15[bn0_ma_reg1] =  A_B0R9[bn0_ma_reg4];
			   A_B0R1[bn0_ma_reg2]  =  data_tmp_1;
			   A_B1R1[bn1_ma_reg2]  =  data_tmp_2;
			   A_B0R1[bn0_ma_reg3]  =  data_tmp_3;
			   A_B1R1[bn1_ma_reg3]  =  data_tmp_4;
			   A_B1R1[bn1_ma_reg4]  =  data_tmp_5;
			   A_B0R1[bn0_ma_reg4]  =  data_tmp_6;
			   A_B0R9[bn0_ma_reg2]  =  data_tmp_9;
			   A_B1R9[bn1_ma_reg2]  =  data_tmp_10;
			   A_B0R9[bn0_ma_reg3]  =  data_tmp_11;
			   A_B1R9[bn1_ma_reg3]  =  data_tmp_12;
			   A_B1R9[bn1_ma_reg4]  =  data_tmp_13;
			   A_B0R9[bn0_ma_reg4]  =  data_tmp_14;
			   //************************************
			   data_tmp_1  =  A_B0R3[bn0_ma_reg2];
			   data_tmp_2  =  A_B0R4[bn0_ma_reg2];
			   data_tmp_3  =  A_B0R5[bn0_ma_reg2];
			   data_tmp_4  =  A_B0R6[bn0_ma_reg2];
			   data_tmp_5  =  A_B0R7[bn0_ma_reg2];
			   data_tmp_9  =  A_B0R11[bn0_ma_reg2];
			   data_tmp_10 =  A_B0R12[bn0_ma_reg2];
			   data_tmp_11 =  A_B0R13[bn0_ma_reg2];
			   data_tmp_12 =  A_B0R14[bn0_ma_reg2];
			   data_tmp_13 =  A_B0R15[bn0_ma_reg2];
			   A_B0R3[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg2];
			   A_B0R4[bn0_ma_reg2]  =  A_B0R2[bn0_ma_reg3];
			   A_B0R5[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg3];
			   A_B0R6[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg4];
			   A_B0R7[bn0_ma_reg2]  =  A_B0R2[bn0_ma_reg4];
			   A_B0R11[bn0_ma_reg2] =  A_B1R10[bn1_ma_reg2];
			   A_B0R12[bn0_ma_reg2] =  A_B0R10[bn0_ma_reg3];
			   A_B0R13[bn0_ma_reg2] =  A_B1R10[bn1_ma_reg3];
			   A_B0R14[bn0_ma_reg2] =  A_B1R10[bn1_ma_reg4];
			   A_B0R15[bn0_ma_reg2] =  A_B0R10[bn0_ma_reg4];
			   A_B1R2[bn1_ma_reg2]  =  data_tmp_1;
			   A_B0R2[bn0_ma_reg3]  =  data_tmp_2;
			   A_B1R2[bn1_ma_reg3]  =  data_tmp_3;
			   A_B1R2[bn1_ma_reg4]  =  data_tmp_4;
			   A_B0R2[bn0_ma_reg4]  =  data_tmp_5;
			   A_B1R10[bn1_ma_reg2]  =  data_tmp_9;
			   A_B0R10[bn0_ma_reg3]  =  data_tmp_10;
			   A_B1R10[bn1_ma_reg3]  =  data_tmp_11;
			   A_B1R10[bn1_ma_reg4]  =  data_tmp_12;
			   A_B0R10[bn0_ma_reg4]  =  data_tmp_13;
			   //***********************************
               data_tmp_1  =  A_B1R4[bn1_ma_reg2];
               data_tmp_2  =  A_B1R5[bn1_ma_reg2];
               data_tmp_3  =  A_B1R6[bn1_ma_reg2];
               data_tmp_4  =  A_B1R7[bn1_ma_reg2];
               data_tmp_9  =  A_B1R12[bn1_ma_reg2];
               data_tmp_10 =  A_B1R13[bn1_ma_reg2];
               data_tmp_11 =  A_B1R14[bn1_ma_reg2];
               data_tmp_12 =  A_B1R15[bn1_ma_reg2];
               A_B1R4[bn1_ma_reg2]  =   A_B0R3[bn0_ma_reg3];
			   A_B1R5[bn1_ma_reg2]  =   A_B1R3[bn1_ma_reg3];
			   A_B1R6[bn1_ma_reg2]  =   A_B1R3[bn1_ma_reg4];
			   A_B1R7[bn1_ma_reg2]  =   A_B0R3[bn0_ma_reg4];
			   A_B1R12[bn1_ma_reg2] =   A_B0R11[bn0_ma_reg3];
               A_B1R13[bn1_ma_reg2] =   A_B1R11[bn1_ma_reg3];
               A_B1R14[bn1_ma_reg2] =   A_B1R11[bn1_ma_reg4];
               A_B1R15[bn1_ma_reg2] =   A_B0R11[bn0_ma_reg4];
               A_B0R3[bn0_ma_reg3]  =   data_tmp_1; 
			   A_B1R3[bn1_ma_reg3]  =   data_tmp_2; 
			   A_B1R3[bn1_ma_reg4]  =   data_tmp_3; 
			   A_B0R3[bn0_ma_reg4]  =   data_tmp_4; 
			   A_B0R11[bn0_ma_reg3] =   data_tmp_9; 
			   A_B1R11[bn1_ma_reg3] =   data_tmp_10;
               A_B1R11[bn1_ma_reg4] =   data_tmp_11;
               A_B0R11[bn0_ma_reg4] =   data_tmp_12;
               //*******************************************
			   data_tmp_1  =  A_B0R5[bn0_ma_reg3];
			   data_tmp_2  =  A_B0R6[bn0_ma_reg3];
			   data_tmp_3  =  A_B0R7[bn0_ma_reg3];
			   data_tmp_9  =  A_B0R13[bn0_ma_reg3];
			   data_tmp_10 =  A_B0R14[bn0_ma_reg3];
			   data_tmp_11 =  A_B0R15[bn0_ma_reg3];
			   A_B0R5[bn0_ma_reg3]  =  A_B1R4[bn1_ma_reg3];
			   A_B0R6[bn0_ma_reg3]  =  A_B1R4[bn1_ma_reg4];
			   A_B0R7[bn0_ma_reg3]  =  A_B0R4[bn0_ma_reg4];
			   A_B0R13[bn0_ma_reg3] =  A_B1R12[bn1_ma_reg3];
			   A_B0R14[bn0_ma_reg3] =  A_B1R12[bn1_ma_reg4];
			   A_B0R15[bn0_ma_reg3] =  A_B0R12[bn0_ma_reg4];
			   A_B1R4[bn1_ma_reg3]  =  data_tmp_1; 
			   A_B1R4[bn1_ma_reg4]  =  data_tmp_2; 
			   A_B0R4[bn0_ma_reg4]  =  data_tmp_3; 
			   A_B1R12[bn1_ma_reg3] =  data_tmp_9; 
			   A_B1R12[bn1_ma_reg4] =  data_tmp_10;
			   A_B0R12[bn0_ma_reg4] =  data_tmp_11;
			   //*******************************************
			   data_tmp_1  =  A_B1R6[bn1_ma_reg3];
			   data_tmp_2  =  A_B1R7[bn1_ma_reg3];
			   data_tmp_9  =  A_B1R14[bn1_ma_reg3];
			   data_tmp_10 =  A_B1R15[bn1_ma_reg3];
			   A_B1R6[bn1_ma_reg3]  = A_B1R5[bn1_ma_reg4];
			   A_B1R7[bn1_ma_reg3]  = A_B0R5[bn0_ma_reg4];
			   A_B1R14[bn1_ma_reg3] = A_B1R13[bn1_ma_reg4];
			   A_B1R15[bn1_ma_reg3] = A_B0R13[bn0_ma_reg4];
               A_B1R5[bn1_ma_reg4]  = data_tmp_1;  
               A_B0R5[bn0_ma_reg4]  = data_tmp_2;  
               A_B1R13[bn1_ma_reg4] = data_tmp_9;  
			   A_B0R13[bn0_ma_reg4] = data_tmp_10;
               //******************************************
               data_tmp_1  =  A_B1R7[bn1_ma_reg4];
			   data_tmp_9  =  A_B1R15[bn1_ma_reg4];
               A_B1R7[bn1_ma_reg4]   =  A_B0R6[bn0_ma_reg4];
			   A_B1R15[bn1_ma_reg4]  =  A_B0R14[bn0_ma_reg4];
               A_B0R6[bn0_ma_reg4]   =  data_tmp_1;
			   A_B0R14[bn0_ma_reg4]  =  data_tmp_9;	
			   //-----------------------------------------------------
			   //-----------------------------------------------------
			   data_tmp_1  = A_B0R1[bn0_ma_reg5];
			   data_tmp_2  = A_B0R2[bn0_ma_reg5];
			   data_tmp_3  = A_B0R3[bn0_ma_reg5];
			   data_tmp_4  = A_B0R4[bn0_ma_reg5];
			   data_tmp_5  = A_B0R5[bn0_ma_reg5];
			   data_tmp_6  = A_B0R6[bn0_ma_reg5];
			   data_tmp_7  = A_B0R7[bn0_ma_reg5];
			   data_tmp_9  = A_B0R9[bn0_ma_reg5];
			   data_tmp_10 = A_B0R10[bn0_ma_reg5];
			   data_tmp_11 = A_B0R11[bn0_ma_reg5];
			   data_tmp_12 = A_B0R12[bn0_ma_reg5];
			   data_tmp_13 = A_B0R13[bn0_ma_reg5];
			   data_tmp_14 = A_B0R14[bn0_ma_reg5];
			   data_tmp_15 = A_B0R15[bn0_ma_reg5];
			   A_B0R1[bn0_ma_reg5]  =  A_B1R0[bn1_ma_reg5];
			   A_B0R2[bn0_ma_reg5]  =  A_B1R0[bn1_ma_reg6];
			   A_B0R3[bn0_ma_reg5]  =  A_B0R0[bn0_ma_reg6];
			   A_B0R4[bn0_ma_reg5]  =  A_B1R0[bn1_ma_reg7];
			   A_B0R5[bn0_ma_reg5]  =  A_B0R0[bn0_ma_reg7];
			   A_B0R6[bn0_ma_reg5]  =  A_B0R0[bn0_ma_reg8];
			   A_B0R7[bn0_ma_reg5]  =  A_B1R0[bn1_ma_reg8];
			   A_B0R9[bn0_ma_reg5]  =  A_B1R8[bn1_ma_reg5];
			   A_B0R10[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg6];
			   A_B0R11[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg6];
			   A_B0R12[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg7];
			   A_B0R13[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg7];
			   A_B0R14[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg8];
			   A_B0R15[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg8];
			   A_B1R0[bn1_ma_reg5]  =  data_tmp_1;
			   A_B1R0[bn1_ma_reg6]  =  data_tmp_2;
			   A_B0R0[bn0_ma_reg6]  =  data_tmp_3;
			   A_B1R0[bn1_ma_reg7]  =  data_tmp_4;
			   A_B0R0[bn0_ma_reg7]  =  data_tmp_5;
			   A_B0R0[bn0_ma_reg8]  =  data_tmp_6;
			   A_B1R0[bn1_ma_reg8]  =  data_tmp_7;
			   A_B1R8[bn1_ma_reg5]  =  data_tmp_9;
			   A_B1R8[bn1_ma_reg6]  =  data_tmp_10; 
			   A_B0R8[bn0_ma_reg6]  =  data_tmp_11; 
			   A_B1R8[bn1_ma_reg7]  =  data_tmp_12; 
			   A_B0R8[bn0_ma_reg7]  =  data_tmp_13; 
			   A_B0R8[bn0_ma_reg8]  =  data_tmp_14; 
			   A_B1R8[bn1_ma_reg8]  =  data_tmp_15; 
			   //**********************************
			   data_tmp_1  = A_B1R2[bn1_ma_reg5];
			   data_tmp_2  = A_B1R3[bn1_ma_reg5];
			   data_tmp_3  = A_B1R4[bn1_ma_reg5];
			   data_tmp_4  = A_B1R5[bn1_ma_reg5];
			   data_tmp_5  = A_B1R6[bn1_ma_reg5];
			   data_tmp_6  = A_B1R7[bn1_ma_reg5];
			   data_tmp_9  = A_B1R10[bn1_ma_reg5];
			   data_tmp_10 = A_B1R11[bn1_ma_reg5];
			   data_tmp_11 = A_B1R12[bn1_ma_reg5];
			   data_tmp_12 = A_B1R13[bn1_ma_reg5];
			   data_tmp_13 = A_B1R14[bn1_ma_reg5];
			   data_tmp_14 = A_B1R15[bn1_ma_reg5];
			   A_B1R2[bn1_ma_reg5]   =  A_B1R1[bn1_ma_reg6];
			   A_B1R3[bn1_ma_reg5]   =  A_B0R1[bn0_ma_reg6];
			   A_B1R4[bn1_ma_reg5]   =  A_B1R1[bn1_ma_reg7];
			   A_B1R5[bn1_ma_reg5]   =  A_B0R1[bn0_ma_reg7];
			   A_B1R6[bn1_ma_reg5]   =  A_B0R1[bn0_ma_reg8];
			   A_B1R7[bn1_ma_reg5]   =  A_B1R1[bn1_ma_reg8];
			   A_B1R10[bn1_ma_reg5]  =  A_B1R9[bn1_ma_reg6];
			   A_B1R11[bn1_ma_reg5]  =  A_B0R9[bn0_ma_reg6];
			   A_B1R12[bn1_ma_reg5]  =  A_B1R9[bn1_ma_reg7];
			   A_B1R13[bn1_ma_reg5]  =  A_B0R9[bn0_ma_reg7];
			   A_B1R14[bn1_ma_reg5]  =  A_B0R9[bn0_ma_reg8];
			   A_B1R15[bn1_ma_reg5]  =  A_B1R9[bn1_ma_reg8];
			   A_B1R1[bn1_ma_reg6] =  data_tmp_1;
			   A_B0R1[bn0_ma_reg6] =  data_tmp_2;
			   A_B1R1[bn1_ma_reg7] =  data_tmp_3;
			   A_B0R1[bn0_ma_reg7] =  data_tmp_4;
			   A_B0R1[bn0_ma_reg8] =  data_tmp_5;
			   A_B1R1[bn1_ma_reg8] =  data_tmp_6;
			   A_B1R9[bn1_ma_reg6] =  data_tmp_9;
			   A_B0R9[bn0_ma_reg6] =  data_tmp_10; 
			   A_B1R9[bn1_ma_reg7] =  data_tmp_11; 
			   A_B0R9[bn0_ma_reg7] =  data_tmp_12; 
			   A_B0R9[bn0_ma_reg8] =  data_tmp_13; 
			   A_B1R9[bn1_ma_reg8] =  data_tmp_14; 
			   //*************************************
			   data_tmp_1  =  A_B1R3[bn1_ma_reg6];
			   data_tmp_2  =  A_B1R4[bn1_ma_reg6];
			   data_tmp_3  =  A_B1R5[bn1_ma_reg6];
			   data_tmp_4  =  A_B1R6[bn1_ma_reg6];
			   data_tmp_5  =  A_B1R7[bn1_ma_reg6];
			   data_tmp_9  =  A_B1R11[bn1_ma_reg6];
			   data_tmp_10 =  A_B1R12[bn1_ma_reg6];
			   data_tmp_11 =  A_B1R13[bn1_ma_reg6];
			   data_tmp_12 =  A_B1R14[bn1_ma_reg6];
			   data_tmp_13 =  A_B1R15[bn1_ma_reg6];
			   A_B1R3[bn1_ma_reg6]  =  A_B0R2[bn0_ma_reg6];
			   A_B1R4[bn1_ma_reg6]  =  A_B1R2[bn1_ma_reg7];
			   A_B1R5[bn1_ma_reg6]  =  A_B0R2[bn0_ma_reg7];
			   A_B1R6[bn1_ma_reg6]  =  A_B0R2[bn0_ma_reg8];
			   A_B1R7[bn1_ma_reg6]  =  A_B1R2[bn1_ma_reg8];
			   A_B1R11[bn1_ma_reg6] =  A_B0R10[bn0_ma_reg6];
			   A_B1R12[bn1_ma_reg6] =  A_B1R10[bn1_ma_reg7];
			   A_B1R13[bn1_ma_reg6] =  A_B0R10[bn0_ma_reg7];
			   A_B1R14[bn1_ma_reg6] =  A_B0R10[bn0_ma_reg8];
			   A_B1R15[bn1_ma_reg6] =  A_B1R10[bn1_ma_reg8];
			   A_B0R2[bn0_ma_reg6]  =  data_tmp_1;
			   A_B1R2[bn1_ma_reg7]  =  data_tmp_2;
			   A_B0R2[bn0_ma_reg7]  =  data_tmp_3;
			   A_B0R2[bn0_ma_reg8]  =  data_tmp_4;
			   A_B1R2[bn1_ma_reg8]  =  data_tmp_5;
			   A_B0R10[bn0_ma_reg6] =  data_tmp_9;
			   A_B1R10[bn1_ma_reg7] =  data_tmp_10;
			   A_B0R10[bn0_ma_reg7] =  data_tmp_11;
			   A_B0R10[bn0_ma_reg8] =  data_tmp_12;
			   A_B1R10[bn1_ma_reg8] =  data_tmp_13;
			   //************************************
			   data_tmp_1  =  A_B0R4[bn0_ma_reg6];
			   data_tmp_2  =  A_B0R5[bn0_ma_reg6];
			   data_tmp_3  =  A_B0R6[bn0_ma_reg6];
			   data_tmp_4  =  A_B0R7[bn0_ma_reg6];
			   data_tmp_9  =  A_B0R12[bn0_ma_reg6];
			   data_tmp_10 =  A_B0R13[bn0_ma_reg6];
 			   data_tmp_11 =  A_B0R14[bn0_ma_reg6];
			   data_tmp_12 =  A_B0R15[bn0_ma_reg6];
			   A_B0R4[bn0_ma_reg6]  =  A_B1R3[bn1_ma_reg7];
			   A_B0R5[bn0_ma_reg6]  =  A_B0R3[bn0_ma_reg7];
			   A_B0R6[bn0_ma_reg6]  =  A_B0R3[bn0_ma_reg8];
			   A_B0R7[bn0_ma_reg6]  =  A_B1R3[bn1_ma_reg8];
			   A_B0R12[bn0_ma_reg6] =  A_B1R11[bn1_ma_reg7];
			   A_B0R13[bn0_ma_reg6] =  A_B0R11[bn0_ma_reg7];
			   A_B0R14[bn0_ma_reg6] =  A_B0R11[bn0_ma_reg8];
			   A_B0R15[bn0_ma_reg6] =  A_B1R11[bn1_ma_reg8];
			   A_B1R3[bn1_ma_reg7]  =  data_tmp_1;
			   A_B0R3[bn0_ma_reg7]  =  data_tmp_2;
			   A_B0R3[bn0_ma_reg8]  =  data_tmp_3;
			   A_B1R3[bn1_ma_reg8]  =  data_tmp_4;
			   A_B1R11[bn1_ma_reg7] =  data_tmp_9;
			   A_B0R11[bn0_ma_reg7] =  data_tmp_10; 
			   A_B0R11[bn0_ma_reg8] =  data_tmp_11; 
			   A_B1R11[bn1_ma_reg8] =  data_tmp_12; 
			   //***********************************
			   data_tmp_1  = A_B1R5[bn1_ma_reg7];
			   data_tmp_2  = A_B1R6[bn1_ma_reg7];
			   data_tmp_3  = A_B1R7[bn1_ma_reg7];
			   data_tmp_9  = A_B1R13[bn1_ma_reg7];
			   data_tmp_10 = A_B1R14[bn1_ma_reg7];
			   data_tmp_11 = A_B1R15[bn1_ma_reg7];
			   A_B1R5[bn1_ma_reg7]  = A_B0R4[bn0_ma_reg7];
			   A_B1R6[bn1_ma_reg7]  = A_B0R4[bn0_ma_reg8];
			   A_B1R7[bn1_ma_reg7]  = A_B1R4[bn1_ma_reg8];
			   A_B1R13[bn1_ma_reg7] = A_B0R12[bn0_ma_reg7];
			   A_B1R14[bn1_ma_reg7] = A_B0R12[bn0_ma_reg8];
			   A_B1R15[bn1_ma_reg7] = A_B1R12[bn1_ma_reg8];
			   A_B0R4[bn0_ma_reg7]  = data_tmp_1;
			   A_B0R4[bn0_ma_reg8]  = data_tmp_2;
			   A_B1R4[bn1_ma_reg8]  = data_tmp_3;
			   A_B0R12[bn0_ma_reg7] = data_tmp_9;
			   A_B0R12[bn0_ma_reg8] = data_tmp_10;
			   A_B1R12[bn1_ma_reg8] = data_tmp_11;
			   //*********************************
			   data_tmp_1  =  A_B0R6[bn0_ma_reg7];
			   data_tmp_2  =  A_B0R7[bn0_ma_reg7];
			   data_tmp_9  =  A_B0R14[bn0_ma_reg7];
			   data_tmp_10 =  A_B0R15[bn0_ma_reg7];
			   A_B0R6[bn0_ma_reg7]  = A_B0R5[bn0_ma_reg8];
			   A_B0R7[bn0_ma_reg7]  = A_B1R5[bn1_ma_reg8];
			   A_B0R14[bn0_ma_reg7] = A_B0R13[bn0_ma_reg8];
			   A_B0R15[bn0_ma_reg7] = A_B1R13[bn1_ma_reg8];
			   A_B0R5[bn0_ma_reg8]  = data_tmp_1;
			   A_B1R5[bn1_ma_reg8]  = data_tmp_2;
			   A_B0R13[bn0_ma_reg8] = data_tmp_9;
			   A_B1R13[bn1_ma_reg8] = data_tmp_10;
			   //********************************
			   data_tmp_1 = A_B0R7[bn0_ma_reg8];
			   data_tmp_9 = A_B0R15[bn0_ma_reg8];
			   A_B0R7[bn0_ma_reg8]  = A_B1R6[bn1_ma_reg8];
			   A_B0R15[bn0_ma_reg8] = A_B1R14[bn1_ma_reg8];
			   A_B1R6[bn1_ma_reg8]  = data_tmp_1;
			   A_B1R14[bn1_ma_reg8] = data_tmp_9;				
            }
		 } 
		}
	}
	
	std::cout << "radix-16 FFT computing over!!\n";
	
	if(display == 1)DATARECORD <<" ----------------------------------------------- \n";
	if(display == 1)DATARECORD <<" Radix-8 FFT computing start!!! \n";
	
	for(int i = 0;i < group; i++){
        for(int j=0; j < radix; j++){
        	gray_i  = Gray(i,group);
        	BC_tmp  = j * group + gray_i;
			if(display == 1)DATARECORD <<" ----------------------------------------------- \n";
			if(display == 1)DATARECORD <<"i: "<< i <<" \n";
			if(display == 1)DATARECORD <<"j: "<< j <<" \n";
			if(display == 1)DATARECORD <<"gray_i: "<< gray_i <<" \n";
			if(display == 1)DATARECORD <<"BC_tmp: "<< BC_tmp <<" \n";
			BC_RR   = BC_tmp;
			if(display == 1)DATARECORD <<"BC_RR: "<< BC_RR <<" \n";
			BC_IFFT_Reorder_R16_R8(BC_RR,BC_Reorder);
        	AGU_R16(BC_Reorder,bn_tmp,ma_tmp);
			if(display == 1)DATARECORD <<"BC_Reorder: "<< BC_Reorder <<" \n";
			if(display == 1)DATARECORD <<"bn_tmp: "<< bn_tmp <<" \n";
			if(display == 1)DATARECORD <<"ma_tmp: "<< ma_tmp <<" \n";
        	if(bn_tmp == 0){
                if(display == 1)DATARECORD <<" Before radix-8 BU operation!! \n";					
                if(display == 1)DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";
                Radix8_BU_INTT(A_B0R0[ma_tmp],A_B0R1[ma_tmp],A_B0R2[ma_tmp],A_B0R3[ma_tmp],A_B0R4[ma_tmp],A_B0R5[ma_tmp],A_B0R6[ma_tmp],A_B0R7[ma_tmp]);
        		Radix8_BU_INTT(A_B0R8[ma_tmp],A_B0R9[ma_tmp],A_B0R10[ma_tmp],A_B0R11[ma_tmp],A_B0R12[ma_tmp],A_B0R13[ma_tmp],A_B0R14[ma_tmp],A_B0R15[ma_tmp]);
                if(display == 1)DATARECORD <<"After radix-8 BU operation!! \n";					
                if(display == 1)DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";						 
        	}else {
                if(display == 1)DATARECORD <<"Before radix-8 BU operation!! \n";					
                if(display == 1)DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";				
        	    Radix8_BU_INTT(A_B1R0[ma_tmp],A_B1R1[ma_tmp],A_B1R2[ma_tmp],A_B1R3[ma_tmp],A_B1R4[ma_tmp],A_B1R5[ma_tmp],A_B1R6[ma_tmp],A_B1R7[ma_tmp]);
        	    Radix8_BU_INTT(A_B1R8[ma_tmp],A_B1R9[ma_tmp],A_B1R10[ma_tmp],A_B1R11[ma_tmp],A_B1R12[ma_tmp],A_B1R13[ma_tmp],A_B1R14[ma_tmp],A_B1R15[ma_tmp]);
                if(display == 1)DATARECORD <<"After radix-8 BU operation!! \n";					
                if(display == 1)DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";					
                if(display == 1)DATARECORD <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";							  
        	}			
        }		
	}
	
	int index0_i;
	int index1_i;
	int index2_i;
	int index3_i;
	int index4_i;
	int index5_i;
	int index6_i;
	int index7_i;
	int index8_i;
	int index9_i;
	int index10_i;
	int index11_i;
	int index12_i;
	int index13_i;
	int index14_i;
	int index15_i;
	
	int index0;
    int index1;
    int index2;
    int index3;
    int index4;
    int index5;
    int index6;
    int index7;
    int index8;
    int index9;
    int index10;
    int index11;
    int index12;
    int index13;
    int index14;
    int index15;
	
	//Mult by 1/N 
	for(int i = 0; i < word_size; i++){
		MulMod(A_B0R0[i],A_B0R0[i],IN,p);
		MulMod(A_B0R1[i],A_B0R1[i],IN,p);
		MulMod(A_B0R2[i],A_B0R2[i],IN,p);
		MulMod(A_B0R3[i],A_B0R3[i],IN,p);
		MulMod(A_B0R4[i],A_B0R4[i],IN,p);
		MulMod(A_B0R5[i],A_B0R5[i],IN,p);
		MulMod(A_B0R6[i],A_B0R6[i],IN,p);
		MulMod(A_B0R7[i],A_B0R7[i],IN,p);
		MulMod(A_B0R8[i],A_B0R8[i],IN,p);
		MulMod(A_B0R9[i],A_B0R9[i],IN,p);
		MulMod(A_B0R10[i],A_B0R10[i],IN,p);
		MulMod(A_B0R11[i],A_B0R11[i],IN,p);
		MulMod(A_B0R12[i],A_B0R12[i],IN,p);
		MulMod(A_B0R13[i],A_B0R13[i],IN,p);
		MulMod(A_B0R14[i],A_B0R14[i],IN,p);
		MulMod(A_B0R15[i],A_B0R15[i],IN,p);
		//bank1 
		MulMod(A_B1R0[i],A_B1R0[i],IN,p);
		MulMod(A_B1R1[i],A_B1R1[i],IN,p);
		MulMod(A_B1R2[i],A_B1R2[i],IN,p);
		MulMod(A_B1R3[i],A_B1R3[i],IN,p);
		MulMod(A_B1R4[i],A_B1R4[i],IN,p);
		MulMod(A_B1R5[i],A_B1R5[i],IN,p);
		MulMod(A_B1R6[i],A_B1R6[i],IN,p);
		MulMod(A_B1R7[i],A_B1R7[i],IN,p);
		MulMod(A_B1R8[i],A_B1R8[i],IN,p);
		MulMod(A_B1R9[i],A_B1R9[i],IN,p);
		MulMod(A_B1R10[i],A_B1R10[i],IN,p);
		MulMod(A_B1R11[i],A_B1R11[i],IN,p);
		MulMod(A_B1R12[i],A_B1R12[i],IN,p);
		MulMod(A_B1R13[i],A_B1R13[i],IN,p);
		MulMod(A_B1R14[i],A_B1R14[i],IN,p);
		MulMod(A_B1R15[i],A_B1R15[i],IN,p);		
	}	
	
	std::ofstream INTT_Before_Relocation("./INTT_R16_R8_Before_Relocation.txt");
	std::ofstream INTT_After_Relocation("./INTT_R16_R8_After_Relocation.txt");
	
	//data relocation for continuous index
	INTT_Before_Relocation <<" ----------------------------------------------- \n";
	INTT_Before_Relocation <<" Data Relocation!!!!! \n";
	INTT_Before_Relocation <<" ----------------------------------------------- \n";
	for(int i = 0;i < group; i++){
        for(int j=0; j < radix; j++){
        	gray_i  = Gray(i,group);
        	BC_tmp  = j * group + gray_i;
			INTT_Before_Relocation <<" ----------------------------------------------- \n";
			INTT_Before_Relocation <<"BC_tmp: "<< BC_tmp <<" \n";
			BC_RR = BC_tmp;
			INTT_Before_Relocation <<"BC_RR: "<< BC_RR <<" \n";
			BC_IFFT_Reorder_R16_R8(BC_RR,BC_Reorder);
        	AGU_R16(BC_Reorder,bn_tmp,ma_tmp);
			INTT_Before_Relocation <<"BC_Reorder: "<< BC_Reorder <<" \n";
			INTT_Before_Relocation <<"bn_tmp: "<< bn_tmp <<" \n";
			INTT_Before_Relocation <<"ma_tmp: "<< ma_tmp <<" \n";
        	if(bn_tmp == 0){
				if(j < 2)bn0_bc_tmp = BC_tmp;
			    if(j <  2)bn0_ma_reg1 = ma_tmp;
                INTT_Before_Relocation <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]<<"\n";	
                INTT_Before_Relocation <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]<<"\n";	
                INTT_Before_Relocation <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]<<"\n";	
                INTT_Before_Relocation <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]<<"\n";	
                INTT_Before_Relocation <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]<<"\n";	
                INTT_Before_Relocation <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]<<"\n";	
                INTT_Before_Relocation <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]<<"\n";	
                INTT_Before_Relocation <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]<<"\n";	
                INTT_Before_Relocation <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]<<"\n";	
                INTT_Before_Relocation <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]<<"\n";	
                INTT_Before_Relocation <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";	
                INTT_Before_Relocation <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";	
                INTT_Before_Relocation <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";	
                INTT_Before_Relocation <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";	
                INTT_Before_Relocation <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";	
                INTT_Before_Relocation <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";
			    if((j >= 2)  && (j < 4))bn0_ma_reg2 = ma_tmp;
			    if((j >= 4)  && (j < 6))bn0_ma_reg3 = ma_tmp;
			    if((j >= 6)  && (j < 8))bn0_ma_reg4 = ma_tmp;
			    if((j >= 8)  && (j < 10))bn0_ma_reg5 = ma_tmp;
			    if((j >= 10) && (j < 12))bn0_ma_reg6 = ma_tmp;
			    if((j >= 12) && (j < 14))bn0_ma_reg7 = ma_tmp;
			    if((j >= 14) && (j < 16))bn0_ma_reg8 = ma_tmp;                						 
        	}else {
				if(j < 2)bn1_bc_tmp = BC_tmp;
                if(j <  2)bn1_ma_reg1 = ma_tmp;					
                INTT_Before_Relocation <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]<<"\n";	
                INTT_Before_Relocation <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]<<"\n";	
                INTT_Before_Relocation <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]<<"\n";	
                INTT_Before_Relocation <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]<<"\n";	
                INTT_Before_Relocation <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]<<"\n";	
                INTT_Before_Relocation <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]<<"\n";	
                INTT_Before_Relocation <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]<<"\n";	
                INTT_Before_Relocation <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]<<"\n";	
                INTT_Before_Relocation <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]<<"\n";	
                INTT_Before_Relocation <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]<<"\n";	
                INTT_Before_Relocation <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";	
                INTT_Before_Relocation <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";	
                INTT_Before_Relocation <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";	
                INTT_Before_Relocation <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";	
                INTT_Before_Relocation <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";	
                INTT_Before_Relocation <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";				
                if((j >= 2)  && (j < 4))bn1_ma_reg2 = ma_tmp;
                if((j >= 4)  && (j < 6))bn1_ma_reg3 = ma_tmp;
                if((j >= 6)  && (j < 8))bn1_ma_reg4 = ma_tmp;
                if((j >= 8)  && (j < 10))bn1_ma_reg5 = ma_tmp;
                if((j >= 10) && (j < 12))bn1_ma_reg6 = ma_tmp;
                if((j >= 12) && (j < 14))bn1_ma_reg7 = ma_tmp;
                if((j >= 14) && (j < 16))bn1_ma_reg8 = ma_tmp;								  
        	}			
        }
		// 16*16 Matrix transpose
		if(bn1_bc_tmp > bn0_bc_tmp){
			 data_tmp_1  = A_B0R1[bn0_ma_reg1];
			 data_tmp_2  = A_B0R2[bn0_ma_reg1];
			 data_tmp_3  = A_B0R3[bn0_ma_reg1];
			 data_tmp_4  = A_B0R4[bn0_ma_reg1];
			 data_tmp_5  = A_B0R5[bn0_ma_reg1];
			 data_tmp_6  = A_B0R6[bn0_ma_reg1];
			 data_tmp_7  = A_B0R7[bn0_ma_reg1];
			 data_tmp_8  = A_B0R8[bn0_ma_reg1];
			 data_tmp_9  = A_B0R9[bn0_ma_reg1];
			 data_tmp_10 = A_B0R10[bn0_ma_reg1];
			 data_tmp_11 = A_B0R11[bn0_ma_reg1];
			 data_tmp_12 = A_B0R12[bn0_ma_reg1];
			 data_tmp_13 = A_B0R13[bn0_ma_reg1];
			 data_tmp_14 = A_B0R14[bn0_ma_reg1];
			 data_tmp_15 = A_B0R15[bn0_ma_reg1];
			 A_B0R1[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg1];
			 A_B0R2[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg2];
			 A_B0R3[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg2];
			 A_B0R4[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg3];
			 A_B0R5[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg3];
			 A_B0R6[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg1]  = A_B1R0[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg1]  = A_B0R0[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg1] = A_B0R0[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg1] = A_B1R0[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg1] = A_B0R0[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg1] = A_B1R0[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg1] = A_B1R0[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg1] = A_B0R0[bn0_ma_reg8];
			 A_B1R0[bn1_ma_reg1]  = data_tmp_1; 
			 A_B1R0[bn1_ma_reg2]  = data_tmp_2; 
			 A_B0R0[bn0_ma_reg2]  = data_tmp_3; 
			 A_B1R0[bn1_ma_reg3]  = data_tmp_4; 
			 A_B0R0[bn0_ma_reg3]  = data_tmp_5; 
			 A_B0R0[bn0_ma_reg4]  = data_tmp_6; 
			 A_B1R0[bn1_ma_reg4]  = data_tmp_7; 
			 A_B1R0[bn1_ma_reg5]  = data_tmp_8; 
			 A_B0R0[bn0_ma_reg5]  = data_tmp_9; 
			 A_B0R0[bn0_ma_reg6]  = data_tmp_10;
			 A_B1R0[bn1_ma_reg6]  = data_tmp_11;
			 A_B0R0[bn0_ma_reg7]  = data_tmp_12;
			 A_B1R0[bn1_ma_reg7]  = data_tmp_13;
			 A_B1R0[bn1_ma_reg8]  = data_tmp_14;
			 A_B0R0[bn0_ma_reg8]  = data_tmp_15;
			 /*********************************************/
			 data_tmp_1  = A_B1R2[bn1_ma_reg1];
			 data_tmp_2  = A_B1R3[bn1_ma_reg1];
			 data_tmp_3  = A_B1R4[bn1_ma_reg1];
			 data_tmp_4  = A_B1R5[bn1_ma_reg1];
			 data_tmp_5  = A_B1R6[bn1_ma_reg1];
			 data_tmp_6  = A_B1R7[bn1_ma_reg1];
			 data_tmp_7  = A_B1R8[bn1_ma_reg1];
			 data_tmp_8  = A_B1R9[bn1_ma_reg1];
			 data_tmp_9  = A_B1R10[bn1_ma_reg1];
			 data_tmp_10 = A_B1R11[bn1_ma_reg1];
			 data_tmp_11 = A_B1R12[bn1_ma_reg1];
			 data_tmp_12 = A_B1R13[bn1_ma_reg1];
			 data_tmp_13 = A_B1R14[bn1_ma_reg1];
			 data_tmp_14 = A_B1R15[bn1_ma_reg1];
             A_B1R2[bn1_ma_reg1]  =	 A_B1R1[bn1_ma_reg2];
             A_B1R3[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg2];
             A_B1R4[bn1_ma_reg1]  =  A_B1R1[bn1_ma_reg3];
             A_B1R5[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg3];
             A_B1R6[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg4];
             A_B1R7[bn1_ma_reg1]  =  A_B1R1[bn1_ma_reg4];
             A_B1R8[bn1_ma_reg1]  =  A_B1R1[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg1]  =  A_B0R1[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg1] =  A_B0R1[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg1] =  A_B1R1[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg1] =  A_B0R1[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg1] =  A_B1R1[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg1] =  A_B1R1[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg1] =  A_B0R1[bn0_ma_reg8];
			 A_B1R1[bn1_ma_reg2]  =  data_tmp_1; 
			 A_B0R1[bn0_ma_reg2]  =  data_tmp_2; 
			 A_B1R1[bn1_ma_reg3]  =  data_tmp_3; 
			 A_B0R1[bn0_ma_reg3]  =  data_tmp_4; 
			 A_B0R1[bn0_ma_reg4]  =  data_tmp_5; 
			 A_B1R1[bn1_ma_reg4]  =  data_tmp_6; 
			 A_B1R1[bn1_ma_reg5]  =  data_tmp_7; 
			 A_B0R1[bn0_ma_reg5]  =  data_tmp_8; 
			 A_B0R1[bn0_ma_reg6]  =  data_tmp_9; 
			 A_B1R1[bn1_ma_reg6]  =  data_tmp_10;
			 A_B0R1[bn0_ma_reg7]  =  data_tmp_11;
			 A_B1R1[bn1_ma_reg7]  =  data_tmp_12;
			 A_B1R1[bn1_ma_reg8]  =  data_tmp_13;
			 A_B0R1[bn0_ma_reg8]  =  data_tmp_14;
			/************************************************************/ 
			 data_tmp_1  =   A_B1R3[bn1_ma_reg2];
			 data_tmp_2  =   A_B1R4[bn1_ma_reg2];
			 data_tmp_3  =   A_B1R5[bn1_ma_reg2];
			 data_tmp_4  =   A_B1R6[bn1_ma_reg2];
			 data_tmp_5  =   A_B1R7[bn1_ma_reg2];
			 data_tmp_6  =   A_B1R8[bn1_ma_reg2];
			 data_tmp_7  =   A_B1R9[bn1_ma_reg2];
			 data_tmp_8  =   A_B1R10[bn1_ma_reg2];
			 data_tmp_9  =   A_B1R11[bn1_ma_reg2];
			 data_tmp_10 =   A_B1R12[bn1_ma_reg2];
			 data_tmp_11 =   A_B1R13[bn1_ma_reg2];
			 data_tmp_12 =   A_B1R14[bn1_ma_reg2];
			 data_tmp_13 =   A_B1R15[bn1_ma_reg2];
			 A_B1R3[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg2];
			 A_B1R4[bn1_ma_reg2]  =  A_B1R2[bn1_ma_reg3];
			 A_B1R5[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg3];
			 A_B1R6[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg4];
			 A_B1R7[bn1_ma_reg2]  =  A_B1R2[bn1_ma_reg4];
			 A_B1R8[bn1_ma_reg2]  =  A_B1R2[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg2]  =  A_B0R2[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg2] =  A_B0R2[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg2] =  A_B1R2[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg2] =  A_B0R2[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg2] =  A_B1R2[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg2] =  A_B1R2[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg2] =  A_B0R2[bn0_ma_reg8];
			 A_B0R2[bn0_ma_reg2]  =  data_tmp_1; 
			 A_B1R2[bn1_ma_reg3]  =  data_tmp_2; 
			 A_B0R2[bn0_ma_reg3]  =  data_tmp_3; 
			 A_B0R2[bn0_ma_reg4]  =  data_tmp_4; 
			 A_B1R2[bn1_ma_reg4]  =  data_tmp_5; 
			 A_B1R2[bn1_ma_reg5]  =  data_tmp_6; 
			 A_B0R2[bn0_ma_reg5]  =  data_tmp_7; 
			 A_B0R2[bn0_ma_reg6]  =  data_tmp_8; 
			 A_B1R2[bn1_ma_reg6]  =  data_tmp_9; 
			 A_B0R2[bn0_ma_reg7]  =  data_tmp_10;
			 A_B1R2[bn1_ma_reg7]  =  data_tmp_11;
			 A_B1R2[bn1_ma_reg8]  =  data_tmp_12;
			 A_B0R2[bn0_ma_reg8]  =  data_tmp_13;
			 //******************************************************
			 data_tmp_1  = A_B0R4[bn0_ma_reg2];
			 data_tmp_2  = A_B0R5[bn0_ma_reg2];
			 data_tmp_3  = A_B0R6[bn0_ma_reg2];
			 data_tmp_4  = A_B0R7[bn0_ma_reg2];
			 data_tmp_5  = A_B0R8[bn0_ma_reg2];
			 data_tmp_6  = A_B0R9[bn0_ma_reg2];
			 data_tmp_7  = A_B0R10[bn0_ma_reg2];
			 data_tmp_8  = A_B0R11[bn0_ma_reg2];
			 data_tmp_9  = A_B0R12[bn0_ma_reg2];
			 data_tmp_10 = A_B0R13[bn0_ma_reg2];
			 data_tmp_11 = A_B0R14[bn0_ma_reg2];
			 data_tmp_12 = A_B0R15[bn0_ma_reg2];
			 A_B0R4[bn0_ma_reg2]  =  A_B1R3[bn1_ma_reg3];
			 A_B0R5[bn0_ma_reg2]  =  A_B0R3[bn0_ma_reg3];
			 A_B0R6[bn0_ma_reg2]  =  A_B0R3[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg2]  =  A_B1R3[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg2]  =  A_B1R3[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg2]  =  A_B0R3[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg2] =  A_B0R3[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg2] =  A_B1R3[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg2] =  A_B0R3[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg2] =  A_B1R3[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg2] =  A_B1R3[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg2] =  A_B0R3[bn0_ma_reg8];
			 A_B1R3[bn1_ma_reg3]  =  data_tmp_1;  
			 A_B0R3[bn0_ma_reg3]  =  data_tmp_2;  
			 A_B0R3[bn0_ma_reg4]  =  data_tmp_3;  
			 A_B1R3[bn1_ma_reg4]  =  data_tmp_4;  
			 A_B1R3[bn1_ma_reg5]  =  data_tmp_5;  
			 A_B0R3[bn0_ma_reg5]  =  data_tmp_6;  
			 A_B0R3[bn0_ma_reg6]  =  data_tmp_7;  
			 A_B1R3[bn1_ma_reg6]  =  data_tmp_8;  
			 A_B0R3[bn0_ma_reg7]  =  data_tmp_9;  
			 A_B1R3[bn1_ma_reg7]  =  data_tmp_10; 
			 A_B1R3[bn1_ma_reg8]  =  data_tmp_11; 
			 A_B0R3[bn0_ma_reg8]  =  data_tmp_12; 
			 //----------------------------------------------------------------------
             data_tmp_1  = A_B1R5[bn1_ma_reg3];
			 data_tmp_2  = A_B1R6[bn1_ma_reg3];
			 data_tmp_3  = A_B1R7[bn1_ma_reg3];
			 data_tmp_4  = A_B1R8[bn1_ma_reg3];
			 data_tmp_5  = A_B1R9[bn1_ma_reg3];
			 data_tmp_6  = A_B1R10[bn1_ma_reg3];
			 data_tmp_7  = A_B1R11[bn1_ma_reg3];
			 data_tmp_8  = A_B1R12[bn1_ma_reg3];
			 data_tmp_9  = A_B1R13[bn1_ma_reg3];
			 data_tmp_10 = A_B1R14[bn1_ma_reg3];
			 data_tmp_11 = A_B1R15[bn1_ma_reg3];
			 A_B1R5[bn1_ma_reg3]  =  A_B0R4[bn0_ma_reg3];
			 A_B1R6[bn1_ma_reg3]  =  A_B0R4[bn0_ma_reg4];
			 A_B1R7[bn1_ma_reg3]  =  A_B1R4[bn1_ma_reg4];
			 A_B1R8[bn1_ma_reg3]  =  A_B1R4[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg3]  =  A_B0R4[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg3] =  A_B0R4[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg3] =  A_B1R4[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg3] =  A_B0R4[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg3] =  A_B1R4[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg3] =  A_B1R4[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg3] =  A_B0R4[bn0_ma_reg8];
			 A_B0R4[bn0_ma_reg3]  =  data_tmp_1 ;
			 A_B0R4[bn0_ma_reg4]  =  data_tmp_2 ;
			 A_B1R4[bn1_ma_reg4]  =  data_tmp_3 ;
			 A_B1R4[bn1_ma_reg5]  =  data_tmp_4 ;
			 A_B0R4[bn0_ma_reg5]  =  data_tmp_5 ;
			 A_B0R4[bn0_ma_reg6]  =  data_tmp_6 ;
			 A_B1R4[bn1_ma_reg6]  =  data_tmp_7 ;
			 A_B0R4[bn0_ma_reg7]  =  data_tmp_8 ;
			 A_B1R4[bn1_ma_reg7]  =  data_tmp_9 ;
			 A_B1R4[bn1_ma_reg8]  =  data_tmp_10;
			 A_B0R4[bn0_ma_reg8]  =  data_tmp_11;
			 //************************************************************************
			 data_tmp_1  = A_B0R6[bn0_ma_reg3];
			 data_tmp_2  = A_B0R7[bn0_ma_reg3];
			 data_tmp_3  = A_B0R8[bn0_ma_reg3];
			 data_tmp_4  = A_B0R9[bn0_ma_reg3];
			 data_tmp_5  = A_B0R10[bn0_ma_reg3];
			 data_tmp_6  = A_B0R11[bn0_ma_reg3];
			 data_tmp_7  = A_B0R12[bn0_ma_reg3];
			 data_tmp_8  = A_B0R13[bn0_ma_reg3];
			 data_tmp_9  = A_B0R14[bn0_ma_reg3];
			 data_tmp_10 = A_B0R15[bn0_ma_reg3];
			 A_B0R6[bn0_ma_reg3]  =  A_B0R5[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg3]  =  A_B1R5[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg3]  =  A_B1R5[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg3]  =  A_B0R5[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg3] =  A_B0R5[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg3] =  A_B1R5[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg3] =  A_B0R5[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg3] =  A_B1R5[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg3] =  A_B1R5[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg3] =  A_B0R5[bn0_ma_reg8];
			 A_B0R5[bn0_ma_reg4]  =  data_tmp_1; 
			 A_B1R5[bn1_ma_reg4]  =  data_tmp_2; 
			 A_B1R5[bn1_ma_reg5]  =  data_tmp_3; 
			 A_B0R5[bn0_ma_reg5]  =  data_tmp_4; 
			 A_B0R5[bn0_ma_reg6]  =  data_tmp_5; 
			 A_B1R5[bn1_ma_reg6]  =  data_tmp_6; 
			 A_B0R5[bn0_ma_reg7]  =  data_tmp_7; 
			 A_B1R5[bn1_ma_reg7]  =  data_tmp_8; 
			 A_B1R5[bn1_ma_reg8]  =  data_tmp_9; 
			 A_B0R5[bn0_ma_reg8]  =  data_tmp_10;
			 //-----------------------------------------------------------------------
			 data_tmp_1  = A_B0R7[bn0_ma_reg4];
			 data_tmp_2  = A_B0R8[bn0_ma_reg4];
			 data_tmp_3  = A_B0R9[bn0_ma_reg4];
			 data_tmp_4  = A_B0R10[bn0_ma_reg4];
			 data_tmp_5  = A_B0R11[bn0_ma_reg4];
			 data_tmp_6  = A_B0R12[bn0_ma_reg4];
			 data_tmp_7  = A_B0R13[bn0_ma_reg4];
			 data_tmp_8  = A_B0R14[bn0_ma_reg4];
			 data_tmp_9  = A_B0R15[bn0_ma_reg4];
			 A_B0R7[bn0_ma_reg4]  =  A_B1R6[bn1_ma_reg4];
			 A_B0R8[bn0_ma_reg4]  =  A_B1R6[bn1_ma_reg5];
			 A_B0R9[bn0_ma_reg4]  =  A_B0R6[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg4] =  A_B0R6[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg4] =  A_B1R6[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg4] =  A_B0R6[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg4] =  A_B1R6[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg4] =  A_B1R6[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg4] =  A_B0R6[bn0_ma_reg8];
			 A_B1R6[bn1_ma_reg4]  =  data_tmp_1;
			 A_B1R6[bn1_ma_reg5]  =  data_tmp_2;
			 A_B0R6[bn0_ma_reg5]  =  data_tmp_3;
			 A_B0R6[bn0_ma_reg6]  =  data_tmp_4;
			 A_B1R6[bn1_ma_reg6]  =  data_tmp_5;
			 A_B0R6[bn0_ma_reg7]  =  data_tmp_6;
			 A_B1R6[bn1_ma_reg7]  =  data_tmp_7;
			 A_B1R6[bn1_ma_reg8]  =  data_tmp_8;
			 A_B0R6[bn0_ma_reg8]  =  data_tmp_9;
			 //----------------------------------------------------------------------
			 data_tmp_1  = A_B1R8[bn1_ma_reg4];
			 data_tmp_2  = A_B1R9[bn1_ma_reg4];
			 data_tmp_3  = A_B1R10[bn1_ma_reg4];
			 data_tmp_4  = A_B1R11[bn1_ma_reg4];
			 data_tmp_5  = A_B1R12[bn1_ma_reg4];
			 data_tmp_6  = A_B1R13[bn1_ma_reg4];
			 data_tmp_7  = A_B1R14[bn1_ma_reg4];
			 data_tmp_8  = A_B1R15[bn1_ma_reg4];
			 A_B1R8[bn1_ma_reg4]  =  A_B1R7[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg4]  =  A_B0R7[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg4] =  A_B0R7[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg4] =  A_B1R7[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg4] =  A_B0R7[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg4] =  A_B1R7[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg4] =  A_B1R7[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg4] =  A_B0R7[bn0_ma_reg8];
			 A_B1R7[bn1_ma_reg5]  =  data_tmp_1;
			 A_B0R7[bn0_ma_reg5]  =  data_tmp_2;
			 A_B0R7[bn0_ma_reg6]  =  data_tmp_3;
			 A_B1R7[bn1_ma_reg6]  =  data_tmp_4;
			 A_B0R7[bn0_ma_reg7]  =  data_tmp_5;
			 A_B1R7[bn1_ma_reg7]  =  data_tmp_6;
			 A_B1R7[bn1_ma_reg8]  =  data_tmp_7;
			 A_B0R7[bn0_ma_reg8]  =  data_tmp_8;
			 //----------------------------------------------------------------------
			 data_tmp_1 = A_B1R9[bn1_ma_reg5];
			 data_tmp_2 = A_B1R10[bn1_ma_reg5];
			 data_tmp_3 = A_B1R11[bn1_ma_reg5];
			 data_tmp_4 = A_B1R12[bn1_ma_reg5];
			 data_tmp_5 = A_B1R13[bn1_ma_reg5];
			 data_tmp_6 = A_B1R14[bn1_ma_reg5];
			 data_tmp_7 = A_B1R15[bn1_ma_reg5];
			 A_B1R9[bn1_ma_reg5]  = A_B0R8[bn0_ma_reg5];
			 A_B1R10[bn1_ma_reg5] = A_B0R8[bn0_ma_reg6];
			 A_B1R11[bn1_ma_reg5] = A_B1R8[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg5] = A_B0R8[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg5] = A_B1R8[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg5] = A_B1R8[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg5] = A_B0R8[bn0_ma_reg8];
			 A_B0R8[bn0_ma_reg5]  = data_tmp_1;
			 A_B0R8[bn0_ma_reg6]  = data_tmp_2;
			 A_B1R8[bn1_ma_reg6]  = data_tmp_3;
			 A_B0R8[bn0_ma_reg7]  = data_tmp_4;
			 A_B1R8[bn1_ma_reg7]  = data_tmp_5;
			 A_B1R8[bn1_ma_reg8]  = data_tmp_6;
			 A_B0R8[bn0_ma_reg8]  = data_tmp_7;
			 //---------------------------------------------------------------------
			 data_tmp_1  = A_B0R10[bn0_ma_reg5];
			 data_tmp_2  = A_B0R11[bn0_ma_reg5];
			 data_tmp_3  = A_B0R12[bn0_ma_reg5];
			 data_tmp_4  = A_B0R13[bn0_ma_reg5];
			 data_tmp_5  = A_B0R14[bn0_ma_reg5];
			 data_tmp_6  = A_B0R15[bn0_ma_reg5];
			 A_B0R10[bn0_ma_reg5] = A_B0R9[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg5] = A_B1R9[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg5] = A_B0R9[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg5] = A_B1R9[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg5] = A_B1R9[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg5] = A_B0R9[bn0_ma_reg8];
			 A_B0R9[bn0_ma_reg6]  = data_tmp_1;
			 A_B1R9[bn1_ma_reg6]  = data_tmp_2;
			 A_B0R9[bn0_ma_reg7]  = data_tmp_3;
			 A_B1R9[bn1_ma_reg7]  = data_tmp_4;
			 A_B1R9[bn1_ma_reg8]  = data_tmp_5;
			 A_B0R9[bn0_ma_reg8]  = data_tmp_6;
			 //--------------------------------------------------------------------
			 data_tmp_1  = A_B0R11[bn0_ma_reg6];
			 data_tmp_2  = A_B0R12[bn0_ma_reg6];
			 data_tmp_3  = A_B0R13[bn0_ma_reg6];
			 data_tmp_4  = A_B0R14[bn0_ma_reg6];
			 data_tmp_5  = A_B0R15[bn0_ma_reg6];
			 A_B0R11[bn0_ma_reg6] = A_B1R10[bn1_ma_reg6];
			 A_B0R12[bn0_ma_reg6] = A_B0R10[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg6] = A_B1R10[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg6] = A_B1R10[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg6] = A_B0R10[bn0_ma_reg8];
			 A_B1R10[bn1_ma_reg6] = data_tmp_1;
			 A_B0R10[bn0_ma_reg7] = data_tmp_2;
			 A_B1R10[bn1_ma_reg7] = data_tmp_3;
			 A_B1R10[bn1_ma_reg8] = data_tmp_4;
			 A_B0R10[bn0_ma_reg8] = data_tmp_5;
			 //--------------------------------------------------------------------
			 data_tmp_1  = A_B1R12[bn1_ma_reg6];
			 data_tmp_2  = A_B1R13[bn1_ma_reg6];
			 data_tmp_3  = A_B1R14[bn1_ma_reg6];
			 data_tmp_4  = A_B1R15[bn1_ma_reg6];
			 A_B1R12[bn1_ma_reg6] = A_B0R11[bn0_ma_reg7];
			 A_B1R13[bn1_ma_reg6] = A_B1R11[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg6] = A_B1R11[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg6] = A_B0R11[bn0_ma_reg8];
			 A_B0R11[bn0_ma_reg7] = data_tmp_1;
			 A_B1R11[bn1_ma_reg7] = data_tmp_2;
			 A_B1R11[bn1_ma_reg8] = data_tmp_3;
			 A_B0R11[bn0_ma_reg8] = data_tmp_4;
			 //--------------------------------------------------------------------
			 data_tmp_1 = A_B0R13[bn0_ma_reg7];
			 data_tmp_2 = A_B0R14[bn0_ma_reg7];
			 data_tmp_3 = A_B0R15[bn0_ma_reg7];
			 A_B0R13[bn0_ma_reg7] = A_B1R12[bn1_ma_reg7];
			 A_B0R14[bn0_ma_reg7] = A_B1R12[bn1_ma_reg8];
			 A_B0R15[bn0_ma_reg7] = A_B0R12[bn0_ma_reg8];
			 A_B1R12[bn1_ma_reg7] = data_tmp_1;
			 A_B1R12[bn1_ma_reg8] = data_tmp_2;
			 A_B0R12[bn0_ma_reg8] = data_tmp_3;
			 //--------------------------------------------------------------------
			 data_tmp_1 = A_B1R14[bn1_ma_reg7];
			 data_tmp_2 = A_B1R15[bn1_ma_reg7];
			 A_B1R14[bn1_ma_reg7] = A_B1R13[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg7] = A_B0R13[bn0_ma_reg8];
			 A_B1R13[bn1_ma_reg8] = data_tmp_1;
			 A_B0R13[bn0_ma_reg8] = data_tmp_2;
			 //--------------------------------------------------------------------
			 data_tmp_1 = A_B1R15[bn1_ma_reg8];
			 A_B1R15[bn1_ma_reg8] = A_B0R14[bn0_ma_reg8];
			 A_B0R14[bn0_ma_reg8] = data_tmp_1; 
		}else{	 
			 data_tmp_1  = A_B1R1[bn1_ma_reg1];
			 data_tmp_2  = A_B1R2[bn1_ma_reg1];
			 data_tmp_3  = A_B1R3[bn1_ma_reg1];
			 data_tmp_4  = A_B1R4[bn1_ma_reg1];
			 data_tmp_5  = A_B1R5[bn1_ma_reg1];
			 data_tmp_6  = A_B1R6[bn1_ma_reg1];
			 data_tmp_7  = A_B1R7[bn1_ma_reg1];
			 data_tmp_8  = A_B1R8[bn1_ma_reg1];
			 data_tmp_9  = A_B1R9[bn1_ma_reg1];
			 data_tmp_10 = A_B1R10[bn1_ma_reg1];
			 data_tmp_11 = A_B1R11[bn1_ma_reg1];
			 data_tmp_12 = A_B1R12[bn1_ma_reg1];
			 data_tmp_13 = A_B1R13[bn1_ma_reg1];
			 data_tmp_14 = A_B1R14[bn1_ma_reg1];
			 data_tmp_15 = A_B1R15[bn1_ma_reg1];
             A_B1R1[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg1];
             A_B1R2[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg2];
             A_B1R3[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg2];
             A_B1R4[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg3];
             A_B1R5[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg3];
             A_B1R6[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg4];
             A_B1R7[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg4];
             A_B1R8[bn1_ma_reg1]  =  A_B0R0[bn0_ma_reg5];
             A_B1R9[bn1_ma_reg1]  =  A_B1R0[bn1_ma_reg5];
             A_B1R10[bn1_ma_reg1] =  A_B1R0[bn1_ma_reg6];
             A_B1R11[bn1_ma_reg1] =  A_B0R0[bn0_ma_reg6];
             A_B1R12[bn1_ma_reg1] =  A_B1R0[bn1_ma_reg7];
             A_B1R13[bn1_ma_reg1] =  A_B0R0[bn0_ma_reg7];
             A_B1R14[bn1_ma_reg1] =  A_B0R0[bn0_ma_reg8];
             A_B1R15[bn1_ma_reg1] =  A_B1R0[bn1_ma_reg8];
             A_B0R0[bn0_ma_reg1]  =  data_tmp_1; 
             A_B0R0[bn0_ma_reg2]  =  data_tmp_2; 
             A_B1R0[bn1_ma_reg2]  =  data_tmp_3; 
             A_B0R0[bn0_ma_reg3]  =  data_tmp_4; 
             A_B1R0[bn1_ma_reg3]  =  data_tmp_5; 
             A_B1R0[bn1_ma_reg4]  =  data_tmp_6; 
             A_B0R0[bn0_ma_reg4]  =  data_tmp_7; 
             A_B0R0[bn0_ma_reg5]  =  data_tmp_8; 
             A_B1R0[bn1_ma_reg5]  =  data_tmp_9; 
             A_B1R0[bn1_ma_reg6]  =  data_tmp_10;
             A_B0R0[bn0_ma_reg6]  =  data_tmp_11;
             A_B1R0[bn1_ma_reg7]  =  data_tmp_12;
             A_B0R0[bn0_ma_reg7]  =  data_tmp_13;
             A_B0R0[bn0_ma_reg8]  =  data_tmp_14;
             A_B1R0[bn1_ma_reg8]  =  data_tmp_15;
             //-----------------------------------------------------------
             data_tmp_1  = A_B0R2[bn0_ma_reg1];
			 data_tmp_2  = A_B0R3[bn0_ma_reg1];
			 data_tmp_3  = A_B0R4[bn0_ma_reg1];
			 data_tmp_4  = A_B0R5[bn0_ma_reg1];
			 data_tmp_5  = A_B0R6[bn0_ma_reg1];
			 data_tmp_6  = A_B0R7[bn0_ma_reg1];
			 data_tmp_7  = A_B0R8[bn0_ma_reg1];
			 data_tmp_8  = A_B0R9[bn0_ma_reg1];
			 data_tmp_9  = A_B0R10[bn0_ma_reg1];
			 data_tmp_10 = A_B0R11[bn0_ma_reg1];
			 data_tmp_11 = A_B0R12[bn0_ma_reg1];
			 data_tmp_12 = A_B0R13[bn0_ma_reg1];
			 data_tmp_13 = A_B0R14[bn0_ma_reg1];
			 data_tmp_14 = A_B0R15[bn0_ma_reg1];
			 A_B0R2[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg2];
			 A_B0R3[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg2];
			 A_B0R4[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg3];
			 A_B0R5[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg3];
			 A_B0R6[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg4];
			 A_B0R7[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg1]  =  A_B0R1[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg1]  =  A_B1R1[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg1] =  A_B1R1[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg1] =  A_B0R1[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg1] =  A_B1R1[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg1] =  A_B0R1[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg1] =  A_B0R1[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg1] =  A_B1R1[bn1_ma_reg8];
             A_B0R1[bn0_ma_reg2]  =  data_tmp_1; 
             A_B1R1[bn1_ma_reg2]  =  data_tmp_2; 
             A_B0R1[bn0_ma_reg3]  =  data_tmp_3; 
             A_B1R1[bn1_ma_reg3]  =  data_tmp_4; 
             A_B1R1[bn1_ma_reg4]  =  data_tmp_5; 
             A_B0R1[bn0_ma_reg4]  =  data_tmp_6; 
             A_B0R1[bn0_ma_reg5]  =  data_tmp_7; 
             A_B1R1[bn1_ma_reg5]  =  data_tmp_8; 
             A_B1R1[bn1_ma_reg6]  =  data_tmp_9; 
             A_B0R1[bn0_ma_reg6]  =  data_tmp_10;
             A_B1R1[bn1_ma_reg7]  =  data_tmp_11;
             A_B0R1[bn0_ma_reg7]  =  data_tmp_12;
             A_B0R1[bn0_ma_reg8]  =  data_tmp_13;
             A_B1R1[bn1_ma_reg8]  =  data_tmp_14;
             //------------------------------------------------------------
             data_tmp_1  =  A_B0R3[bn0_ma_reg2];
			 data_tmp_2  =  A_B0R4[bn0_ma_reg2];
			 data_tmp_3  =  A_B0R5[bn0_ma_reg2];
			 data_tmp_4  =  A_B0R6[bn0_ma_reg2];
			 data_tmp_5  =  A_B0R7[bn0_ma_reg2];
			 data_tmp_6  =  A_B0R8[bn0_ma_reg2];
			 data_tmp_7  =  A_B0R9[bn0_ma_reg2];
			 data_tmp_8  =  A_B0R10[bn0_ma_reg2];
			 data_tmp_9  =  A_B0R11[bn0_ma_reg2];
			 data_tmp_10 =  A_B0R12[bn0_ma_reg2];
			 data_tmp_11 =  A_B0R13[bn0_ma_reg2];
			 data_tmp_12 =  A_B0R14[bn0_ma_reg2];
			 data_tmp_13 =  A_B0R15[bn0_ma_reg2];
			 A_B0R3[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg2];
			 A_B0R4[bn0_ma_reg2]  =  A_B0R2[bn0_ma_reg3];
			 A_B0R5[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg3];
			 A_B0R6[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg4];
			 A_B0R7[bn0_ma_reg2]  =  A_B0R2[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg2]  =  A_B0R2[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg2]  =  A_B1R2[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg2] =  A_B1R2[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg2] =  A_B0R2[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg2] =  A_B1R2[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg2] =  A_B0R2[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg2] =  A_B0R2[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg2] =  A_B1R2[bn1_ma_reg8];
			 A_B1R2[bn1_ma_reg2]  =  data_tmp_1; 
			 A_B0R2[bn0_ma_reg3]  =  data_tmp_2; 
			 A_B1R2[bn1_ma_reg3]  =  data_tmp_3; 
			 A_B1R2[bn1_ma_reg4]  =  data_tmp_4; 
			 A_B0R2[bn0_ma_reg4]  =  data_tmp_5; 
			 A_B0R2[bn0_ma_reg5]  =  data_tmp_6; 
			 A_B1R2[bn1_ma_reg5]  =  data_tmp_7; 
			 A_B1R2[bn1_ma_reg6]  =  data_tmp_8; 
			 A_B0R2[bn0_ma_reg6]  =  data_tmp_9; 
			 A_B1R2[bn1_ma_reg7]  =  data_tmp_10;
			 A_B0R2[bn0_ma_reg7]  =  data_tmp_11;
			 A_B0R2[bn0_ma_reg8]  =  data_tmp_12;
			 A_B1R2[bn1_ma_reg8]  =  data_tmp_13;
			 //-----------------------------------------------------------                        
			 data_tmp_1  = A_B1R4[bn1_ma_reg2];
			 data_tmp_2  = A_B1R5[bn1_ma_reg2];
			 data_tmp_3  = A_B1R6[bn1_ma_reg2];
			 data_tmp_4  = A_B1R7[bn1_ma_reg2];
			 data_tmp_5  = A_B1R8[bn1_ma_reg2];
			 data_tmp_6  = A_B1R9[bn1_ma_reg2];
			 data_tmp_7  = A_B1R10[bn1_ma_reg2];
			 data_tmp_8  = A_B1R11[bn1_ma_reg2];
			 data_tmp_9  = A_B1R12[bn1_ma_reg2];
			 data_tmp_10 = A_B1R13[bn1_ma_reg2];
			 data_tmp_11 = A_B1R14[bn1_ma_reg2];
			 data_tmp_12 = A_B1R15[bn1_ma_reg2];
			 A_B1R4[bn1_ma_reg2]  = A_B0R3[bn0_ma_reg3];
			 A_B1R5[bn1_ma_reg2]  = A_B1R3[bn1_ma_reg3];
			 A_B1R6[bn1_ma_reg2]  = A_B1R3[bn1_ma_reg4];
			 A_B1R7[bn1_ma_reg2]  = A_B0R3[bn0_ma_reg4];
			 A_B1R8[bn1_ma_reg2]  = A_B0R3[bn0_ma_reg5];
			 A_B1R9[bn1_ma_reg2]  = A_B1R3[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg2] = A_B1R3[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg2] = A_B0R3[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg2] = A_B1R3[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg2] = A_B0R3[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg2] = A_B0R3[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg2] = A_B1R3[bn1_ma_reg8];
			 A_B0R3[bn0_ma_reg3]  = data_tmp_1; 
			 A_B1R3[bn1_ma_reg3]  = data_tmp_2; 
			 A_B1R3[bn1_ma_reg4]  = data_tmp_3; 
			 A_B0R3[bn0_ma_reg4]  = data_tmp_4; 
			 A_B0R3[bn0_ma_reg5]  = data_tmp_5; 
			 A_B1R3[bn1_ma_reg5]  = data_tmp_6; 
			 A_B1R3[bn1_ma_reg6]  = data_tmp_7; 
			 A_B0R3[bn0_ma_reg6]  = data_tmp_8; 
			 A_B1R3[bn1_ma_reg7]  = data_tmp_9; 
			 A_B0R3[bn0_ma_reg7]  = data_tmp_10;
			 A_B0R3[bn0_ma_reg8]  = data_tmp_11;
			 A_B1R3[bn1_ma_reg8]  = data_tmp_12;
			 //------------------------------------------------------------
			 data_tmp_1  =  A_B0R5[bn0_ma_reg3];
			 data_tmp_2  =  A_B0R6[bn0_ma_reg3];
			 data_tmp_3  =  A_B0R7[bn0_ma_reg3];
			 data_tmp_4  =  A_B0R8[bn0_ma_reg3];
			 data_tmp_5  =  A_B0R9[bn0_ma_reg3];
			 data_tmp_6  =  A_B0R10[bn0_ma_reg3];
			 data_tmp_7  =  A_B0R11[bn0_ma_reg3];
			 data_tmp_8  =  A_B0R12[bn0_ma_reg3];
			 data_tmp_9  =  A_B0R13[bn0_ma_reg3];
			 data_tmp_10 =  A_B0R14[bn0_ma_reg3];
			 data_tmp_11 =  A_B0R15[bn0_ma_reg3];
			 A_B0R5[bn0_ma_reg3]  =  A_B1R4[bn1_ma_reg3];
			 A_B0R6[bn0_ma_reg3]  =  A_B1R4[bn1_ma_reg4];
			 A_B0R7[bn0_ma_reg3]  =  A_B0R4[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg3]  =  A_B0R4[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg3]  =  A_B1R4[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg3] =  A_B1R4[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg3] =  A_B0R4[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg3] =  A_B1R4[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg3] =  A_B0R4[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg3] =  A_B0R4[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg3] =  A_B1R4[bn1_ma_reg8];
			 A_B1R4[bn1_ma_reg3]  =  data_tmp_1; 
			 A_B1R4[bn1_ma_reg4]  =  data_tmp_2; 
			 A_B0R4[bn0_ma_reg4]  =  data_tmp_3; 
			 A_B0R4[bn0_ma_reg5]  =  data_tmp_4; 
			 A_B1R4[bn1_ma_reg5]  =  data_tmp_5; 
			 A_B1R4[bn1_ma_reg6]  =  data_tmp_6; 
			 A_B0R4[bn0_ma_reg6]  =  data_tmp_7; 
			 A_B1R4[bn1_ma_reg7]  =  data_tmp_8; 
			 A_B0R4[bn0_ma_reg7]  =  data_tmp_9; 
			 A_B0R4[bn0_ma_reg8]  =  data_tmp_10;
			 A_B1R4[bn1_ma_reg8]  =  data_tmp_11;
			 //-------------------------------------------------------------
			 data_tmp_1  =  A_B1R6[bn1_ma_reg3];
			 data_tmp_2  =  A_B1R7[bn1_ma_reg3];
			 data_tmp_3  =  A_B1R8[bn1_ma_reg3];
			 data_tmp_4  =  A_B1R9[bn1_ma_reg3];
			 data_tmp_5  =  A_B1R10[bn1_ma_reg3];
			 data_tmp_6  =  A_B1R11[bn1_ma_reg3];
			 data_tmp_7  =  A_B1R12[bn1_ma_reg3];
			 data_tmp_8  =  A_B1R13[bn1_ma_reg3];
			 data_tmp_9  =  A_B1R14[bn1_ma_reg3];
			 data_tmp_10 =  A_B1R15[bn1_ma_reg3];
			 A_B1R6[bn1_ma_reg3]  =  A_B1R5[bn1_ma_reg4];
			 A_B1R7[bn1_ma_reg3]  =  A_B0R5[bn0_ma_reg4];
			 A_B1R8[bn1_ma_reg3]  =  A_B0R5[bn0_ma_reg5];
			 A_B1R9[bn1_ma_reg3]  =  A_B1R5[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg3] =  A_B1R5[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg3] =  A_B0R5[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg3] =  A_B1R5[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg3] =  A_B0R5[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg3] =  A_B0R5[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg3] =  A_B1R5[bn1_ma_reg8];
			 A_B1R5[bn1_ma_reg4]  =  data_tmp_1; 
			 A_B0R5[bn0_ma_reg4]  =  data_tmp_2; 
			 A_B0R5[bn0_ma_reg5]  =  data_tmp_3; 
			 A_B1R5[bn1_ma_reg5]  =  data_tmp_4; 
			 A_B1R5[bn1_ma_reg6]  =  data_tmp_5; 
			 A_B0R5[bn0_ma_reg6]  =  data_tmp_6; 
			 A_B1R5[bn1_ma_reg7]  =  data_tmp_7; 
			 A_B0R5[bn0_ma_reg7]  =  data_tmp_8; 
			 A_B0R5[bn0_ma_reg8]  =  data_tmp_9; 
			 A_B1R5[bn1_ma_reg8]  =  data_tmp_10;
			 //----------------------------------------------------------------
			 data_tmp_1  =  A_B1R7[bn1_ma_reg4];
			 data_tmp_2  =  A_B1R8[bn1_ma_reg4];
			 data_tmp_3  =  A_B1R9[bn1_ma_reg4];
			 data_tmp_4  =  A_B1R10[bn1_ma_reg4];
			 data_tmp_5  =  A_B1R11[bn1_ma_reg4];
			 data_tmp_6  =  A_B1R12[bn1_ma_reg4];
			 data_tmp_7  =  A_B1R13[bn1_ma_reg4];
			 data_tmp_8  =  A_B1R14[bn1_ma_reg4];
			 data_tmp_9  =  A_B1R15[bn1_ma_reg4];
			 A_B1R7[bn1_ma_reg4]  =  A_B0R6[bn0_ma_reg4];
			 A_B1R8[bn1_ma_reg4]  =  A_B0R6[bn0_ma_reg5];
			 A_B1R9[bn1_ma_reg4]  =  A_B1R6[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg4] =  A_B1R6[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg4] =  A_B0R6[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg4] =  A_B1R6[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg4] =  A_B0R6[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg4] =  A_B0R6[bn0_ma_reg8];
             A_B1R15[bn1_ma_reg4] =  A_B1R6[bn1_ma_reg8];
             A_B0R6[bn0_ma_reg4]  =  data_tmp_1;
             A_B0R6[bn0_ma_reg5]  =  data_tmp_2;
             A_B1R6[bn1_ma_reg5]  =  data_tmp_3;
             A_B1R6[bn1_ma_reg6]  =  data_tmp_4;
             A_B0R6[bn0_ma_reg6]  =  data_tmp_5;
             A_B1R6[bn1_ma_reg7]  =  data_tmp_6;
             A_B0R6[bn0_ma_reg7]  =  data_tmp_7;
             A_B0R6[bn0_ma_reg8]  =  data_tmp_8;
             A_B1R6[bn1_ma_reg8]  =  data_tmp_9;
             //------------------------------------------------------------------
			 data_tmp_1  = A_B0R8[bn0_ma_reg4];
			 data_tmp_2  = A_B0R9[bn0_ma_reg4];
			 data_tmp_3  = A_B0R10[bn0_ma_reg4];
			 data_tmp_4  = A_B0R11[bn0_ma_reg4];
			 data_tmp_5  = A_B0R12[bn0_ma_reg4];
			 data_tmp_6  = A_B0R13[bn0_ma_reg4];
			 data_tmp_7  = A_B0R14[bn0_ma_reg4];
			 data_tmp_8  = A_B0R15[bn0_ma_reg4];
			 A_B0R8[bn0_ma_reg4]  =  A_B0R7[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg4]  =  A_B1R7[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg4] =  A_B1R7[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg4] =  A_B0R7[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg4] =  A_B1R7[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg4] =  A_B0R7[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg4] =  A_B0R7[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg4] =  A_B1R7[bn1_ma_reg8];
			 A_B0R7[bn0_ma_reg5]  =  data_tmp_1;
			 A_B1R7[bn1_ma_reg5]  =  data_tmp_2;
			 A_B1R7[bn1_ma_reg6]  =  data_tmp_3;
			 A_B0R7[bn0_ma_reg6]  =  data_tmp_4;
			 A_B1R7[bn1_ma_reg7]  =  data_tmp_5;
			 A_B0R7[bn0_ma_reg7]  =  data_tmp_6;
			 A_B0R7[bn0_ma_reg8]  =  data_tmp_7;
			 A_B1R7[bn1_ma_reg8]  =  data_tmp_8;
			 //----------------------------------------------------------------
			 data_tmp_1  = A_B0R9[bn0_ma_reg5];
			 data_tmp_2  = A_B0R10[bn0_ma_reg5];
			 data_tmp_3  = A_B0R11[bn0_ma_reg5];
			 data_tmp_4  = A_B0R12[bn0_ma_reg5];
			 data_tmp_5  = A_B0R13[bn0_ma_reg5];
			 data_tmp_6  = A_B0R14[bn0_ma_reg5];
			 data_tmp_7  = A_B0R15[bn0_ma_reg5];
			 A_B0R9[bn0_ma_reg5]  =  A_B1R8[bn1_ma_reg5];
			 A_B0R10[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg6];
			 A_B0R11[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg5] =  A_B0R8[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg5] =  A_B1R8[bn1_ma_reg8];
			 A_B1R8[bn1_ma_reg5]  =  data_tmp_1;
			 A_B1R8[bn1_ma_reg6]  =  data_tmp_2;
			 A_B0R8[bn0_ma_reg6]  =  data_tmp_3;
			 A_B1R8[bn1_ma_reg7]  =  data_tmp_4;
			 A_B0R8[bn0_ma_reg7]  =  data_tmp_5;
			 A_B0R8[bn0_ma_reg8]  =  data_tmp_6;
			 A_B1R8[bn1_ma_reg8]  =  data_tmp_7;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B1R10[bn1_ma_reg5];
			 data_tmp_2 = A_B1R11[bn1_ma_reg5];
			 data_tmp_3 = A_B1R12[bn1_ma_reg5];
			 data_tmp_4 = A_B1R13[bn1_ma_reg5];
			 data_tmp_5 = A_B1R14[bn1_ma_reg5];
			 data_tmp_6 = A_B1R15[bn1_ma_reg5];
			 A_B1R10[bn1_ma_reg5] = A_B1R9[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg5] = A_B0R9[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg5] = A_B1R9[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg5] = A_B0R9[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg5] = A_B0R9[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg5] = A_B1R9[bn1_ma_reg8];
			 A_B1R9[bn1_ma_reg6]  = data_tmp_1;
			 A_B0R9[bn0_ma_reg6]  = data_tmp_2;
			 A_B1R9[bn1_ma_reg7]  = data_tmp_3;
			 A_B0R9[bn0_ma_reg7]  = data_tmp_4;
			 A_B0R9[bn0_ma_reg8]  = data_tmp_5;
			 A_B1R9[bn1_ma_reg8]  = data_tmp_6;
			 //----------------------------------------------------------------
			 data_tmp_1 =  A_B1R11[bn1_ma_reg6];
			 data_tmp_2 =  A_B1R12[bn1_ma_reg6];
			 data_tmp_3 =  A_B1R13[bn1_ma_reg6];
			 data_tmp_4 =  A_B1R14[bn1_ma_reg6];
			 data_tmp_5 =  A_B1R15[bn1_ma_reg6];
			 A_B1R11[bn1_ma_reg6]  = A_B0R10[bn0_ma_reg6];
			 A_B1R12[bn1_ma_reg6]  = A_B1R10[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg6]  = A_B0R10[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg6]  = A_B0R10[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg6]  = A_B1R10[bn1_ma_reg8];
			 A_B0R10[bn0_ma_reg6]  = data_tmp_1;
			 A_B1R10[bn1_ma_reg7]  = data_tmp_2;
			 A_B0R10[bn0_ma_reg7]  = data_tmp_3;
			 A_B0R10[bn0_ma_reg8]  = data_tmp_4;
			 A_B1R10[bn1_ma_reg8]  = data_tmp_5;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B0R12[bn0_ma_reg6];
			 data_tmp_2 = A_B0R13[bn0_ma_reg6];
			 data_tmp_3 = A_B0R14[bn0_ma_reg6];
			 data_tmp_4 = A_B0R15[bn0_ma_reg6];
			 A_B0R12[bn0_ma_reg6] = A_B1R11[bn1_ma_reg7];
			 A_B0R13[bn0_ma_reg6] = A_B0R11[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg6] = A_B0R11[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg6] = A_B1R11[bn1_ma_reg8];
			 A_B1R11[bn1_ma_reg7] = data_tmp_1;
			 A_B0R11[bn0_ma_reg7] = data_tmp_2;
			 A_B0R11[bn0_ma_reg8] = data_tmp_3;
			 A_B1R11[bn1_ma_reg8] = data_tmp_4;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B1R13[bn1_ma_reg7];
			 data_tmp_2 = A_B1R14[bn1_ma_reg7];
			 data_tmp_3 = A_B1R15[bn1_ma_reg7];
			 A_B1R13[bn1_ma_reg7] = A_B0R12[bn0_ma_reg7];
			 A_B1R14[bn1_ma_reg7] = A_B0R12[bn0_ma_reg8];
			 A_B1R15[bn1_ma_reg7] = A_B1R12[bn1_ma_reg8];
			 A_B0R12[bn0_ma_reg7] = data_tmp_1;
			 A_B0R12[bn0_ma_reg8] = data_tmp_2;
			 A_B1R12[bn1_ma_reg8] = data_tmp_3;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B0R14[bn0_ma_reg7];
			 data_tmp_2 = A_B0R15[bn0_ma_reg7];
			 A_B0R14[bn0_ma_reg7] = A_B0R13[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg7] = A_B1R13[bn1_ma_reg8];
			 A_B0R13[bn0_ma_reg8] = data_tmp_1;
			 A_B1R13[bn1_ma_reg8] = data_tmp_2;
			 //-----------------------------------------------------------------
			 data_tmp_1 = A_B0R15[bn0_ma_reg8];
			 A_B0R15[bn0_ma_reg8] = A_B1R14[bn1_ma_reg8];
			 A_B1R14[bn1_ma_reg8] = data_tmp_1;
			 //-----------------------------------------------------------------
		}		
	}	

	for(int i = 0;i < group; i++){
        for(int j=0; j < radix; j++){
        	gray_i  = Gray(i,group);
        	BC_tmp  = j * group + gray_i;
			INTT_After_Relocation <<" ----------------------------------------------- \n";
			INTT_After_Relocation <<"BC_tmp: "<< BC_tmp <<" \n";
			BC_RR = BC_tmp;
			INTT_After_Relocation <<"BC_RR: "<< BC_RR <<" \n";
			BC_IFFT_Reorder_R16_R8(BC_RR,BC_Reorder);
        	AGU_R16(BC_Reorder,bn_tmp,ma_tmp);
			INTT_After_Relocation <<"BC_Reorder: "<< BC_Reorder <<" \n";
			INTT_After_Relocation <<"bn_tmp: "<< bn_tmp <<" \n";
			INTT_After_Relocation <<"ma_tmp: "<< ma_tmp <<" \n";
        	if(bn_tmp == 0){
                INTT_After_Relocation <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]<<"\n";	
                INTT_After_Relocation <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]<<"\n";	
                INTT_After_Relocation <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]<<"\n";	
                INTT_After_Relocation <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]<<"\n";	
                INTT_After_Relocation <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]<<"\n";	
                INTT_After_Relocation <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]<<"\n";	
                INTT_After_Relocation <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]<<"\n";	
                INTT_After_Relocation <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]<<"\n";	
                INTT_After_Relocation <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]<<"\n";	
                INTT_After_Relocation <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]<<"\n";	
                INTT_After_Relocation <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";	
                INTT_After_Relocation <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";	
                INTT_After_Relocation <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";	
                INTT_After_Relocation <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";	
                INTT_After_Relocation <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";	
                INTT_After_Relocation <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";							 
        	}else {
                INTT_After_Relocation <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]<<"\n";	
                INTT_After_Relocation <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]<<"\n";	
                INTT_After_Relocation <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]<<"\n";	
                INTT_After_Relocation <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]<<"\n";	
                INTT_After_Relocation <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]<<"\n";	
                INTT_After_Relocation <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]<<"\n";	
                INTT_After_Relocation <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]<<"\n";	
                INTT_After_Relocation <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]<<"\n";	
                INTT_After_Relocation <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]<<"\n";	
                INTT_After_Relocation <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]<<"\n";	
                INTT_After_Relocation <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";	
                INTT_After_Relocation <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";	
                INTT_After_Relocation <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";	
                INTT_After_Relocation <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";	
                INTT_After_Relocation <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";	
                INTT_After_Relocation <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";
        	}			
        }
    }


	DATARECORD <<"***************************************************************\n";
	DATARECORD <<"***** DATA OUTPUT!!                                         ***\n";
	DATARECORD <<"***************************************************************\n";
	//data output
	// SPMB data output , output function is "NTT_REORDERINDEX_R16_R8_OUT".
	for(int i = 0;i < group;i++){
		for(int j = 0;j < radix;j++){
			//gray_i  = Gray(i,group);
			BC_tmp  = j * group + i;
			DATARECORD <<"---------------------------------------------------\n";
			DATARECORD <<" BC_tmp: "<< BC_tmp <<"\n";
			BC_RR = BC_tmp;
			DATARECORD <<" BC_RR: "<< BC_RR <<"\n";
			BC_IFFT_Reorder_R16_R8_OUT(BC_RR,BC_Reorder);
			AGU_R16(BC_Reorder,bn_tmp,ma_tmp);
			DATARECORD <<" BC_Reorder: "<< BC_Reorder <<"\n";
			DATARECORD <<" bn_tmp:  "<< bn_tmp <<"\n";
			DATARECORD <<" ma_tmp:  "<< ma_tmp <<"\n";
			index0_i  = ( 256 * i ) + ( 16 * j );
			index1_i  = ( 256 * i ) + ( 16 * j )  + 1;
			index2_i  = ( 256 * i ) + ( 16 * j )  + 2;
			index3_i  = ( 256 * i ) + ( 16 * j )  + 3;
			index4_i  = ( 256 * i ) + ( 16 * j )  + 4;
			index5_i  = ( 256 * i ) + ( 16 * j )  + 5;
			index6_i  = ( 256 * i ) + ( 16 * j )  + 6;
			index7_i  = ( 256 * i ) + ( 16 * j )  + 7;
			index8_i  = ( 256 * i ) + ( 16 * j )  + 8;
			index9_i  = ( 256 * i ) + ( 16 * j )  + 9;
			index10_i = ( 256 * i ) + ( 16 * j ) + 10;
			index11_i = ( 256 * i ) + ( 16 * j ) + 11;
			index12_i = ( 256 * i ) + ( 16 * j ) + 12;
			index13_i = ( 256 * i ) + ( 16 * j ) + 13;
			index14_i = ( 256 * i ) + ( 16 * j ) + 14;
			index15_i = ( 256 * i ) + ( 16 * j ) + 15;

			DATARECORD <<" index0_i:  "<< index0_i <<"\n";
			DATARECORD <<" index1_i:  "<< index1_i <<"\n";
			DATARECORD <<" index2_i:  "<< index2_i <<"\n";
			DATARECORD <<" index3_i:  "<< index3_i <<"\n";
			DATARECORD <<" index4_i:  "<< index4_i <<"\n";
			DATARECORD <<" index5_i:  "<< index5_i <<"\n";
			DATARECORD <<" index6_i:  "<< index6_i <<"\n";
			DATARECORD <<" index7_i:  "<< index7_i <<"\n";
			DATARECORD <<" index8_i:  "<< index8_i <<"\n";
			DATARECORD <<" index9_i:  "<< index9_i <<"\n";
			DATARECORD <<" index10_i: "<< index10_i <<"\n";
			DATARECORD <<" index11_i: "<< index11_i <<"\n";
			DATARECORD <<" index12_i: "<< index12_i <<"\n";
			DATARECORD <<" index13_i: "<< index13_i <<"\n";
			DATARECORD <<" index14_i: "<< index14_i <<"\n";
			DATARECORD <<" index15_i: "<< index15_i <<"\n";
			index0  = index0_i;
			index1  = index1_i;
			index2  = index2_i;
			index3  = index3_i;
			index4  = index4_i;
			index5  = index5_i;
			index6  = index6_i;
			index7  = index7_i;
			index8  = index8_i;
			index9  = index9_i;
			index10 = index10_i;
			index11 = index11_i;
			index12 = index12_i;
			index13 = index13_i;
			index14 = index14_i;
			index15 = index15_i;
			

            DATARECORD <<" index0:  "<< index0 <<"\n";
			DATARECORD <<" index1:  "<< index1 <<"\n";
			DATARECORD <<" index2:  "<< index2 <<"\n";
			DATARECORD <<" index3:  "<< index3 <<"\n";
			DATARECORD <<" index4:  "<< index4 <<"\n";
			DATARECORD <<" index5:  "<< index5 <<"\n";
			DATARECORD <<" index6:  "<< index6 <<"\n";
			DATARECORD <<" index7:  "<< index7 <<"\n";
			DATARECORD <<" index8:  "<< index8 <<"\n";
			DATARECORD <<" index9:  "<< index9 <<"\n";
			DATARECORD <<" index10: "<< index10 <<"\n";
			DATARECORD <<" index11: "<< index11 <<"\n";
			DATARECORD <<" index12: "<< index12 <<"\n";
			DATARECORD <<" index13: "<< index13 <<"\n";
			DATARECORD <<" index14: "<< index14 <<"\n";
			DATARECORD <<" index15: "<< index15 <<"\n";
			
			if(bn_tmp == 0){
               A[index0]  = A_B0R0[ma_tmp];
			   A[index1]  = A_B0R1[ma_tmp];
			   A[index2]  = A_B0R2[ma_tmp];
			   A[index3]  = A_B0R3[ma_tmp];
			   A[index4]  = A_B0R4[ma_tmp];
			   A[index5]  = A_B0R5[ma_tmp];
			   A[index6]  = A_B0R6[ma_tmp];
			   A[index7]  = A_B0R7[ma_tmp];
			   A[index8]  = A_B0R8[ma_tmp];
			   A[index9]  = A_B0R9[ma_tmp];
			   A[index10] = A_B0R10[ma_tmp];
			   A[index11] = A_B0R11[ma_tmp];
			   A[index12] = A_B0R12[ma_tmp];
			   A[index13] = A_B0R13[ma_tmp];
			   A[index14] = A_B0R14[ma_tmp];
			   A[index15] = A_B0R15[ma_tmp];
               DATARECORD <<"radix-8 BU operation!! \n";					
               DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]<<"\n";					
               DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]<<"\n";					
               DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]<<"\n";					
               DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]<<"\n";					
               DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]<<"\n";					
               DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]<<"\n";					
               DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]<<"\n";					
               DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]<<"\n";					
               DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]<<"\n";					
               DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]<<"\n";					
               DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";					
               DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";					
               DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";					
               DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";					
               DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";					
               DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";			   
			}
			else {
               A[index0]     = A_B1R0[ma_tmp];
               A[index1]     = A_B1R1[ma_tmp];
               A[index2]     = A_B1R2[ma_tmp];
               A[index3]     = A_B1R3[ma_tmp];
               A[index4]     = A_B1R4[ma_tmp];
               A[index5]     = A_B1R5[ma_tmp];
               A[index6]     = A_B1R6[ma_tmp];
               A[index7]     = A_B1R7[ma_tmp];
               A[index8]     = A_B1R8[ma_tmp];
               A[index9]     = A_B1R9[ma_tmp];
               A[index10]    = A_B1R10[ma_tmp];
               A[index11]    = A_B1R11[ma_tmp];
               A[index12]    = A_B1R12[ma_tmp];
               A[index13]    = A_B1R13[ma_tmp];
               A[index14]    = A_B1R14[ma_tmp];
               A[index15]    = A_B1R15[ma_tmp];
               DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]<<"\n";					
               DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]<<"\n";					
               DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]<<"\n";					
               DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]<<"\n";					
               DATARECORD <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]<<"\n";					
               DATARECORD <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]<<"\n";					
               DATARECORD <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]<<"\n";					
               DATARECORD <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]<<"\n";					
               DATARECORD <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]<<"\n";					
               DATARECORD <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]<<"\n";					
               DATARECORD <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";					
               DATARECORD <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";					
               DATARECORD <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";					
               DATARECORD <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";					
               DATARECORD <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";					
               DATARECORD <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";				   
			}			
		}		
	}
}
//***********************************************
//************ Right shift **********************
//***********************************************
void NTTSPMB::RR_R2(int BC,int shift_bit,int &result){
    int    BC_WIDTH;
    int    BC_tmp;
	int    bit_tmp;
	int    weight_tmp;
	std::vector<int> bit_array;
	std::vector<int> bit_array_tmp;
	
    BC_WIDTH = (int)ceil(log2(N/2));
	bit_array.resize(BC_WIDTH);
	bit_array_tmp.resize(BC_WIDTH);
	BC_tmp = BC;
	//bit calculate
	for(int j=0; j < BC_WIDTH;j++){
        bit_tmp = BC_tmp % 2;
        BC_tmp = BC_tmp >> 1;
        bit_array[j]     = bit_tmp;
        bit_array_tmp[j] = bit_tmp;
    } 

	int index;
	result  = 0;
	for(int j=0; j < BC_WIDTH;j++){
        index =  j + shift_bit;
        if(index >= BC_WIDTH) index = index - BC_WIDTH;
		bit_array[j] =  bit_array_tmp[index];
		if(bit_array[j] == 1) weight_tmp = 1 << j;
		else weight_tmp = 0;
		result = result + weight_tmp;
	}		
}
void NTTSPMB::RR_R4_R2(int BC,int shift_bit,int &result){
    int    BC_WIDTH;
    int    BC_tmp;
	int    bit_tmp;
	int    weight_tmp;
	std::vector<int> bit_array;
	std::vector<int> bit_array_tmp;
	
    BC_WIDTH = (int)ceil(log2(N/4));
	bit_array.resize(BC_WIDTH);
	bit_array_tmp.resize(BC_WIDTH);
	BC_tmp = BC;
	//bit calculate
	for(int j=0; j < BC_WIDTH;j++){
        bit_tmp = BC_tmp % 2;
        BC_tmp = BC_tmp >> 1;
        bit_array[j]     = bit_tmp;
        bit_array_tmp[j] = bit_tmp;
    } 

	int index;
	result  = 0;
	for(int j=0; j < BC_WIDTH;j++){
        index =  j + shift_bit;
        if(index >= BC_WIDTH) index = index - BC_WIDTH;
		bit_array[j] =  bit_array_tmp[index];
		if(bit_array[j] == 1) weight_tmp = 1 << j;
		else weight_tmp = 0;
		result = result + weight_tmp;
	}		
}
void NTTSPMB::RR_R4(int BC,int shift_bit,int &result){
    int    BC_WIDTH;
    int    BC_tmp;
	int    bit_tmp;
	int    weight_tmp;
	std::vector<int> bit_array;
	std::vector<int> bit_array_tmp;
	
    BC_WIDTH = (int)ceil(log2(N/4));
	bit_array.resize(BC_WIDTH);
	bit_array_tmp.resize(BC_WIDTH);
	BC_tmp = BC;
	//bit calculate
	for(int j=0; j < BC_WIDTH;j++){
        bit_tmp = BC_tmp % 2;
        BC_tmp = BC_tmp >> 1;
        bit_array[j]     = bit_tmp;
        bit_array_tmp[j] = bit_tmp;
    } 

	int index;
	result  = 0;
	for(int j=0; j < BC_WIDTH;j++){
        index =  j + (shift_bit * 2);
        if(index >= BC_WIDTH) index = index - BC_WIDTH;
		bit_array[j] =  bit_array_tmp[index];
		if(bit_array[j] == 1) weight_tmp = 1 << j;
		else weight_tmp = 0;
		result = result + weight_tmp;
	}		
}
void NTTSPMB::RR_R16(int BC,int shift_bit,int &result){
    int    BC_WIDTH;
    int    BC_tmp;
	int    bit_tmp;
	int    weight_tmp;
	std::vector<int> bit_array;
	std::vector<int> bit_array_tmp;
	
    BC_WIDTH = (int)ceil(log2(N/16));
	bit_array.resize(BC_WIDTH);
	bit_array_tmp.resize(BC_WIDTH);
	BC_tmp = BC;
	//bit calculate
	for(int j=0; j < BC_WIDTH;j++){
        bit_tmp = BC_tmp % 2;
        BC_tmp = BC_tmp >> 1;
        bit_array[j]     = bit_tmp;
        bit_array_tmp[j] = bit_tmp;
    } 

	int index;
	result  = 0;
	for(int j=0; j < BC_WIDTH;j++){
        index =  j + (shift_bit * 4);
        if(index >= BC_WIDTH) index = index - BC_WIDTH;
		bit_array[j] =  bit_array_tmp[index];
		if(bit_array[j] == 1) weight_tmp = 1 << j;
		else weight_tmp = 0;
		result = result + weight_tmp;
	}		
}
void NTTSPMB::RR_R16_R2(int BC,int shift_bit,int &result){
    int     BC_WIDTH;
    int     BC_tmp;
	int     bit_tmp;
	int     weight_tmp;
	
	std::vector<int> bit_array;
	std::vector<int> bit_array_tmp;
	
    BC_WIDTH = (int)ceil(log2(N/16));
	bit_array.resize(BC_WIDTH);
	bit_array_tmp.resize(BC_WIDTH);
	BC_tmp = BC;
	//bit calculate
	for(int j=0; j < BC_WIDTH;j++){
        bit_tmp = BC_tmp % 2;
        BC_tmp = BC_tmp >> 1;
        bit_array[j]     = bit_tmp;
        bit_array_tmp[j] = bit_tmp;
    } 

	int index;
	result  = 0;
	for(int j=0; j < BC_WIDTH;j++){
        index =  j + shift_bit;
        if(index >= BC_WIDTH) index = index - BC_WIDTH;
		bit_array[j] =  bit_array_tmp[index];
		if(bit_array[j] == 1) weight_tmp = 1 << j;
		else weight_tmp = 0;
		result = result + weight_tmp;
	}		
}
void NTTSPMB::RR_R16_R4(int BC,int shift_bit,int &result){
    int     BC_WIDTH;
    int     BC_tmp;
	int     bit_tmp;
	int     weight_tmp;
	
	std::vector<int> bit_array;
	std::vector<int> bit_array_tmp;
	
    BC_WIDTH = (int)ceil(log2(N/16));
	bit_array.resize(BC_WIDTH);
	bit_array_tmp.resize(BC_WIDTH);
	BC_tmp = BC;
	//bit calculate
	for(int j=0; j < BC_WIDTH;j++){
        bit_tmp = BC_tmp % 2;
        BC_tmp = BC_tmp >> 1;
        bit_array[j]     = bit_tmp;
        bit_array_tmp[j] = bit_tmp;
    } 

	int index;
	result  = 0;
	for(int j=0; j < BC_WIDTH;j++){
        index =  j + shift_bit;
        if(index >= BC_WIDTH) index = index - BC_WIDTH;
		bit_array[j] =  bit_array_tmp[index];
		if(bit_array[j] == 1) weight_tmp = 1 << j;
		else weight_tmp = 0;
		result = result + weight_tmp;
	}		
}
void NTTSPMB::RR_R16_R8(int BC,int shift_bit,int &result){
    int     BC_WIDTH;
	int     shift_bit_remainder;
    int     BC_tmp;
	int     bit_tmp;
	int     weight_tmp;
	
	std::vector<int> bit_array;
	std::vector<int> bit_array_tmp;
	
    BC_WIDTH = (int)ceil(log2(N/16));
	bit_array.resize(BC_WIDTH);
	bit_array_tmp.resize(BC_WIDTH);
	BC_tmp = BC;
	//bit calculate
	for(int j=0; j < BC_WIDTH;j++){
        bit_tmp = BC_tmp % 2;
        BC_tmp = BC_tmp >> 1;
        bit_array[j]     = bit_tmp;
        bit_array_tmp[j] = bit_tmp;
    } 

	int index;
	result  = 0;
	
	shift_bit_remainder  = shift_bit  % BC_WIDTH;
	
	for(int j=0; j < BC_WIDTH;j++){
        index =  j + shift_bit_remainder;
        if(index >= BC_WIDTH) index = index - BC_WIDTH;
		bit_array[j] =  bit_array_tmp[index];
		if(bit_array[j] == 1) weight_tmp = 1 << j;
		else weight_tmp = 0;
		result = result + weight_tmp;
	}		
}
void NTTSPMB::BR_R2(int BC,int &result){
    int    fft_point_bit_size;
    int    BC_tmp;
	int    bit_tmp;
	int    weight_tmp;
	std::vector<int> bit_array;
	std::vector<int> bit_array_tmp;
	
    fft_point_bit_size = (int)ceil(log2(N));
	bit_array.resize(fft_point_bit_size);
	bit_array_tmp.resize(fft_point_bit_size);
	BC_tmp = BC;
	//bit calculate
	for(int j=0; j < fft_point_bit_size;j++){
        bit_tmp = BC_tmp % 2;
        BC_tmp = BC_tmp >> 1;
        bit_array[j]     = bit_tmp;
        bit_array_tmp[j] = bit_tmp;
    } 

	result  = 0;
	for(int j=0; j < fft_point_bit_size;j++){
		bit_array[j] =  bit_array_tmp[fft_point_bit_size-j-1];
		if(bit_array[j] == 1) weight_tmp = 1 << j;
		else weight_tmp = 0;
		result = result + weight_tmp;
	}		
}
void NTTSPMB::BR_R4(int BC,int &result){
    int    fft_point_bit_size;
    int    BC_tmp;
	int    bit_tmp;
	int    weight_tmp;
	std::vector<int> bit_array;
	std::vector<int> bit_array_tmp;
	
    fft_point_bit_size = (int)ceil(log2(N));
	bit_array.resize(fft_point_bit_size);
	bit_array_tmp.resize(fft_point_bit_size);
	BC_tmp = BC;
	//bit calculate
	for(int j=0; j < fft_point_bit_size;j++){
        bit_tmp = BC_tmp % 2;
        BC_tmp = BC_tmp >> 1;
        bit_array[j]     = bit_tmp;
        bit_array_tmp[j] = bit_tmp;
		//std::cout <<"bit_array_tmp["<<j<<"]:"<< bit_array_tmp[j] <<"\n";
    } 

	result  = 0;
	for(int j=0; j < fft_point_bit_size/2;j++){
		bit_array[2 * j] =  bit_array_tmp[(fft_point_bit_size-1)- 2*j- 1];
		bit_array[2 * j + 1] =  bit_array_tmp[(fft_point_bit_size-1)- 2*j];
	}	
    

    for(int j=0; j < fft_point_bit_size;j++){
		//std::cout <<"bit_array["<<j<<"]:"<< bit_array[j] <<"\n";
		if(bit_array[j] == 1) weight_tmp = 1 << j;
		else weight_tmp = 0;
		result = result + weight_tmp;
    }
	
}
//------------------------
void NTTSPMB::BR_R16(int BC,int &result){
    int    fft_point_bit_size;
    int    BC_tmp;
	int    bit_tmp;
	int    weight_tmp;
	std::vector<int> bit_array;
	std::vector<int> bit_array_tmp;
	
    fft_point_bit_size = (int)ceil(log2(N));
	bit_array.resize(fft_point_bit_size);
	bit_array_tmp.resize(fft_point_bit_size);
	BC_tmp = BC;
	//bit calculate
	for(int j=0; j < fft_point_bit_size;j++){
        bit_tmp = BC_tmp % 2;
        BC_tmp = BC_tmp >> 1;
        bit_array[j]     = bit_tmp;
        bit_array_tmp[j] = bit_tmp;
		//std::cout <<"bit_array_tmp["<<j<<"]:"<< bit_array_tmp[j] <<"\n";
    } 

	result  = 0;
	for(int j=0; j < fft_point_bit_size/4;j++){
		bit_array[4 * j    ] =  bit_array_tmp[(fft_point_bit_size-1)- 4 * j - 3];
		bit_array[4 * j + 1] =  bit_array_tmp[(fft_point_bit_size-1)- 4 * j - 2];
		bit_array[4 * j + 2] =  bit_array_tmp[(fft_point_bit_size-1)- 4 * j - 1];
		bit_array[4 * j + 3] =  bit_array_tmp[(fft_point_bit_size-1)- 4 * j];
	}	
    

    for(int j=0; j < fft_point_bit_size;j++){
		//std::cout <<"bit_array["<<j<<"]:"<< bit_array[j] <<"\n";
		if(bit_array[j] == 1) weight_tmp = 1 << j;
		else weight_tmp = 0;
		result = result + weight_tmp;
    }
	
}
void NTTSPMB::BR_R4_R2(int BC,int &result){
    int    BC_bit_size;
    int    BC_tmp;
	int    bit_tmp;
	int    weight_tmp;
	std::vector<int> bit_array;
	std::vector<int> bit_array_tmp;
	
    BC_bit_size = (int)ceil(log2(N/4));
	bit_array.resize(BC_bit_size);
	bit_array_tmp.resize(BC_bit_size);
	BC_tmp = BC;
	//bit calculate
	for(int j=0; j < BC_bit_size;j++){
        bit_tmp = BC_tmp % 2;
        BC_tmp = BC_tmp >> 1;
        bit_array[j]     = bit_tmp;
        bit_array_tmp[j] = bit_tmp;
		//std::cout <<"bit_array_tmp["<<j<<"]:"<< bit_array_tmp[j] <<"\n";
    } 

	result  = 0;
	//Orinal lsb ------> MSB
	//Others bit ,two bits as one digit ,then perform reverse operation
	bit_array[BC_bit_size-1] = bit_array_tmp[0];
	for(int j=0; j < (BC_bit_size-1)/2;j++){
		bit_array[2 * j]     =  bit_array_tmp[(BC_bit_size-1)- 2*j- 1];
		bit_array[2 * j + 1] =  bit_array_tmp[(BC_bit_size-1)- 2*j];
	}	
    

    for(int j=0; j < BC_bit_size;j++){
		//std::cout <<"bit_array["<<j<<"]:"<< bit_array[j] <<"\n";
		if(bit_array[j] == 1) weight_tmp = 1 << j;
		else weight_tmp = 0;
		result = result + weight_tmp;
    }
	
}
//***********************************************************************
//*************************** Radix-4 ***********************************
//***********************************************************************
void NTTSPMB::BC_IFFT_Reorder_R4_R2(int BC,int &result){
    int    BC_bit_size;
    int    BC_tmp;
	int    bit_tmp;
	int    weight_tmp;
	int    Isodd;
	int    odd_index_tmp;
	int    even_index_tmp;
	std::vector<int> bit_array;
	std::vector<int> bit_array_tmp;
	
    BC_bit_size = (int)ceil(log2(N/4));
	bit_array.resize(BC_bit_size);
	bit_array_tmp.resize(BC_bit_size);
	BC_tmp = BC;
	//bit calculate
	//std::cout << "---------------------------\n";
	//std::cout << "BC: " << BC <<"\n";
	for(int j=0; j < BC_bit_size;j++){
        bit_tmp = BC_tmp % 2;
        BC_tmp = BC_tmp >> 1;
        bit_array[j]     = bit_tmp;
        bit_array_tmp[j] = bit_tmp;
		//std::cout <<"bit_array_tmp["<<j<<"]:"<< bit_array_tmp[j] <<"\n";
    } 
    
	result  = 0;
	//Orinal lsb ------> MSB
	//Others bit ,two bits as one digit ,then perform reverse operation
	if(BC_bit_size > 3)Isodd = 1;
	else  Isodd = 0;
	even_index_tmp = 2;
	odd_index_tmp  = 3;
    for(int j=BC_bit_size-1; j >= 0;j--){
       if(j == (BC_bit_size-1)) bit_array[j] = bit_array_tmp[1];
       else if (j == (BC_bit_size-2)) bit_array[j] = bit_array_tmp[0];
	   else if (j == 0) bit_array[j] = bit_array_tmp[even_index_tmp];
       else {
		   if(Isodd == 1){
			 bit_array[j] = bit_array_tmp[odd_index_tmp];   
			 odd_index_tmp = odd_index_tmp + 2;
			 Isodd = 0;
		   }else {
			 bit_array[j] = bit_array_tmp[even_index_tmp]; 
             //std::cout << "bit_array["<<j<<"]: " << bit_array[j] <<"\n";
             //std::cout << "bit_array_tmp["<<even_index_tmp<<"]: " << bit_array_tmp[even_index_tmp] <<"\n";
			 even_index_tmp = even_index_tmp + 2;  
			 Isodd = 1;  
		   }
	   }
	}
    

    for(int j=0; j < BC_bit_size;j++){
		//std::cout <<"bit_array["<<j<<"]:"<< bit_array[j] <<"\n";
		if(bit_array[j] == 1) weight_tmp = 1 << j;
		else weight_tmp = 0;
		result = result + weight_tmp;
    }
	
}
void NTTSPMB::REORDERBC_R4_R2_OUT(int BC,int &result){
    int    BC_bit_size;
    int    BC_tmp;
	int    bit_tmp;
	int    weight_tmp;
	std::vector<int> bit_array;
	std::vector<int> bit_array_tmp;
	
	//std::cout << "---------------------------------\n";
	//std::cout << "BC: " << BC <<"\n";
    BC_bit_size = (int)ceil(log2(N/4));
	bit_array.resize(BC_bit_size);
	bit_array_tmp.resize(BC_bit_size);

    //std::cout << "BC_bit_size: " << BC_bit_size <<"\n";	
	//Frist , BC circular right shift 2 bit
	BC_tmp = BC;
	//bit calculate
	for(int j=0; j < BC_bit_size;j++){
        bit_tmp = BC_tmp % 2;
        BC_tmp = BC_tmp >> 1;
        bit_array[j]     = bit_tmp;
        bit_array_tmp[j] = bit_tmp;
		//std::cout <<"bit_array_tmp["<<j<<"]:"<< bit_array_tmp[j] <<"\n";
		//std::cout << "bit_array[" << j <<"]: "<< bit_array[j] <<"\n";	
    } 

	result  = 0;
	//Orinal msb index m
	// X[m-1] left side and right side relocation 
	// ORIGINAL : X[m] X[m-1] X[m-2] X[m-3] ...........  X[0]
	//RELOCATION: X[m-1] X[0] X[m].............. X[m-2] X[m-3]
	int right_shift_lsb_index;
	int isodd = 1;
	int odd_tmp;
	int even_tmp;
	right_shift_lsb_index = BC_bit_size - 2;
	//std::cout << "right_shift_lsb_index: " << right_shift_lsb_index << "\n";
	odd_tmp  = 3;
	even_tmp = 2;
	if(BC_bit_size > 3){
		//std::cout << "---------------------------------------------------------------------------\n";
		//std::cout << "******************   bit_array relocation!                  ***************\n";
		for(int i = BC_bit_size-1; i >= 0;i--){
         	//std::cout << "iteration: " << i << "\n";
			if(i == BC_bit_size-1)bit_array[i] = bit_array_tmp[right_shift_lsb_index-3];
			else if(i == BC_bit_size-2)bit_array[i] = bit_array_tmp[right_shift_lsb_index-1];
			else if(i == BC_bit_size-3)bit_array[i] = bit_array_tmp[right_shift_lsb_index+1];
			else if(i == BC_bit_size-4)bit_array[i] = bit_array_tmp[right_shift_lsb_index];
			else if(i == BC_bit_size-5)bit_array[i] = bit_array_tmp[1];
			else if(i == BC_bit_size-6)bit_array[i] = bit_array_tmp[0];
			else {
				  if(isodd == 1) {
					 bit_array[i] = bit_array_tmp[odd_tmp];
					 odd_tmp = odd_tmp + 2;
					 isodd = 0;
				  }else {
					 bit_array[i] = bit_array_tmp[even_tmp];
					 even_tmp = even_tmp + 2;
					 isodd = 1;
				  }
			}
			//std::cout << "bit_array[" << i <<"]: "<< bit_array[i] <<"\n";
			//std::cout << "*************************************************************************\n";
		}
    }else{
		bit_array[0] = bit_array_tmp[2];
		bit_array[1] = bit_array_tmp[0];
		bit_array[2] = bit_array_tmp[1];
	}

    for(int j=0; j < BC_bit_size;j++){
		//std::cout <<"bit_array["<<j<<"]:"<< bit_array[j] <<"\n";
		if(bit_array[j] == 1) weight_tmp = 1 << j;
		else weight_tmp = 0;
		result = result + weight_tmp;
		//std::cout << "bit_array[" << j <<"]: "<< bit_array[j] <<"\n";
    }
	//std::cout << "result: " << result <<"\n";
}
//INTT , Reorder BC index for output
//Then output data index is inorder.  0,1,2,3,4,5,.....
void NTTSPMB::INTT_REORDERBC_R4_R2_OUT(int BC,int &result){
    int    BC_bit_size;
    int    BC_tmp;
	int    bit_tmp;
	int    weight_tmp;
	std::vector<int> bit_array;
	std::vector<int> bit_array_tmp;
	
	//std::cout << "---------------------------------\n";
	//std::cout << "BC: " << BC <<"\n";
    BC_bit_size = (int)ceil(log2(N/4));
	bit_array.resize(BC_bit_size);
	bit_array_tmp.resize(BC_bit_size);

    //std::cout << "BC_bit_size: " << BC_bit_size <<"\n";	
	//Frist , BC circular right shift 2 bit
	BC_tmp = BC;
	//bit calculate
	for(int j=0; j < BC_bit_size;j++){
        bit_tmp = BC_tmp % 2;
        BC_tmp = BC_tmp >> 1;
        bit_array[j]     = bit_tmp;
        bit_array_tmp[j] = bit_tmp;
		//std::cout <<"bit_array_tmp["<<j<<"]:"<< bit_array_tmp[j] <<"\n";
		//std::cout << "bit_array[" << j <<"]: "<< bit_array[j] <<"\n";	
    } 

	result  = 0;
	//Orinal msb index m
    // For example, 32 point fft ,BC have 3 bits
	// original BC : [2] [1] [0] --------->  [0] [2] [1]
	// 128 point fft , BC  have 5 bits
	// original BC : [4] [3] [2] [1] [0] --------> [3] [1] [2] [4] [0]
	// 512 point , BC have 7 bits 
	// original BC : [6] [5] [4] [3] [2] [1] [0] ------------>  [0] [3] [5] [1] [4] [6] [2]
    // 2048 point , BC have 9 bits 
    // original BC : [8] [7] [6] [5] [4] [3] [2] [1] [0] ------------> [2] [5] [0] [3] [7] [1] [6] [8] [4]
    int reorder_msb_index;	
    int reorder_lsb_index;
	int reorder_p1_index; // reorder BC[1](poistion 1 bit) 
	int reorder_p4_index; // reorder BC[1](poistion 1 bit) 
	int Isodd;
	int odd_index_tmp;
	int even_index_tmp;
	reorder_lsb_index = BC_bit_size - 5;
	reorder_p1_index  = BC_bit_size - 1;
	reorder_p4_index  = BC_bit_size - 2;
	if(reorder_lsb_index == 0) reorder_msb_index = BC_bit_size - 2;
    else reorder_msb_index = reorder_lsb_index - 2;
	//std::cout << "reorder_msb_index: " << reorder_msb_index <<"\n";
	//std::cout << "reorder_lsb_index: " << reorder_lsb_index <<"\n";
	//---------------------------------------------
	//over index 4 of reorder sequence 
	Isodd = 1;
	odd_index_tmp  = BC_bit_size - 4;
	even_index_tmp = BC_bit_size - 9;
	//std::cout << "right_shift_lsb_index: " << right_shift_lsb_index << "\n";
	if(BC_bit_size > 3){
		//std::cout << "---------------------------------------------------------------------------\n";
		//std::cout << "******************   bit_array relocation!                  ***************\n";
		for(int i = BC_bit_size-1; i >= 0;i--){
            if(i == 0) bit_array[i] = bit_array_tmp[reorder_lsb_index];
            else if (i == 1) bit_array[i] = bit_array_tmp[reorder_p1_index];
			else if (i == 2) bit_array[i] = bit_array_tmp[reorder_lsb_index + 2];
			else if (i == 3) bit_array[i] = bit_array_tmp[1];
			else if (i == 4) bit_array[i] = bit_array_tmp[reorder_p4_index];
			else if (i == (BC_bit_size - 1)) bit_array[i] = bit_array_tmp[reorder_msb_index];
			else {
				if(Isodd == 1){
				    bit_array[i] = bit_array_tmp[odd_index_tmp];
                    odd_index_tmp = odd_index_tmp - 2 ;
                    Isodd = 0;					
				}else{
					bit_array[i] = bit_array_tmp[even_index_tmp];
					even_index_tmp = even_index_tmp - 2;
					Isodd = 1;
				}
			}
		}
    }else{
		bit_array[0] = bit_array_tmp[1];
		bit_array[1] = bit_array_tmp[2];
		bit_array[2] = bit_array_tmp[0];
	}

    for(int j=0; j < BC_bit_size;j++){
		//std::cout <<"bit_array["<<j<<"]:"<< bit_array[j] <<"\n";
		if(bit_array[j] == 1) weight_tmp = 1 << j;
		else weight_tmp = 0;
		result = result + weight_tmp;
		//std::cout << "bit_array[" << j <<"]: "<< bit_array[j] <<"\n";
    }
	//std::cout << "result: " << result <<"\n";
}
//***********************************************************************
//*************************** Radix-16 **********************************
//***********************************************************************
//*******************************************************
//               Reindex for NTT output
//*******************************************************
//NTT radix-16 and radix-2 Reindex for output
void NTTSPMB::NTT_REORDERINDEX_R16_R2_OUT(int index,int &result){
    int     fft_bit_size;
    int     index_tmp;
	int     bit_tmp;
	int     weight_tmp;
	int     Number_of_digit;

	std::vector<int> bit_array;
	std::vector<int> bit_array_tmp;
	
	//-----------------------------------------------------
    fft_bit_size = (int)ceil(log2(N));
	bit_array.resize(fft_bit_size);
	bit_array_tmp.resize(fft_bit_size);
	
    //Number of digit , four bits as one digit.
	Number_of_digit = (fft_bit_size - 1) / 4;
	//-----------------------------------------------------
	index_tmp = index;
	//bit array calculate
	for(int i = 0; i < fft_bit_size; i++){
        bit_tmp    = index_tmp % 2;
        index_tmp  = index_tmp >> 1;
        bit_array[i]     = bit_tmp;
        bit_array_tmp[i] = bit_tmp;	
	}
    
	result = 0;
	//MSB relocation , BC[MSB] = BC[0]	
	bit_array[fft_bit_size - 1] = bit_array_tmp[0];
    for(int i = 0; i < Number_of_digit; i++){
       bit_array[fft_bit_size - (4 * i) - 2] = bit_array_tmp[(4 * i) + 4];
       bit_array[fft_bit_size - (4 * i) - 3] = bit_array_tmp[(4 * i) + 3];
	   bit_array[fft_bit_size - (4 * i) - 4] = bit_array_tmp[(4 * i) + 2];
	   bit_array[fft_bit_size - (4 * i) - 5] = bit_array_tmp[(4 * i) + 1];
    }
    
    for(int j = 0; j < fft_bit_size; j++){
		if(bit_array[j] == 1) weight_tmp = 1 << j;
		else weight_tmp = 0;
		result = result + weight_tmp;
    }
    
}
//NTT radix-16 and radix-2 Reindex for output
void NTTSPMB::NTT_REORDERINDEX_R16_R4_OUT(int index,int &result){
    int     fft_bit_size;
    int     index_tmp;
	int     bit_tmp;
	int     weight_tmp;
	int     Number_of_digit;

	std::vector<int> bit_array;
	std::vector<int> bit_array_tmp;
	
	//-----------------------------------------------------
    fft_bit_size = (int)ceil(log2(N));
	bit_array.resize(fft_bit_size);
	bit_array_tmp.resize(fft_bit_size);
	
    //Number of digit , four bits as one digit.
	Number_of_digit = (fft_bit_size - 2) / 4;
	//std::cout << "fft_bit_size: " << fft_bit_size <<"\n";
	//std::cout << "Number_of_digit: " << Number_of_digit <<"\n";
	//-----------------------------------------------------
	index_tmp = index;
	//bit array calculate
	for(int i = 0; i < fft_bit_size; i++){
        bit_tmp    = index_tmp % 2;
        index_tmp  = index_tmp >> 1;
        bit_array[i]     = bit_tmp;
        bit_array_tmp[i] = bit_tmp;	
	}
    
	result = 0;
	//MSB relocation , BC[MSB] = BC[0]	
	bit_array[fft_bit_size - 1] = bit_array_tmp[1];
	bit_array[fft_bit_size - 2] = bit_array_tmp[0];
    for(int i = 0; i < Number_of_digit; i++){
	  if(i==0){
        bit_array[fft_bit_size - (4 * i) - 3] = bit_array_tmp[(4 * i) + 3];		
		bit_array[fft_bit_size - (4 * i) - 4] = bit_array_tmp[(4 * i) + 2];
		bit_array[fft_bit_size - (4 * i) - 5] = bit_array_tmp[(4 * i) + 5];  
		bit_array[fft_bit_size - (4 * i) - 6] = bit_array_tmp[(4 * i) + 4];  
	  }else{
       bit_array[fft_bit_size - (4 * i) - 3] = bit_array_tmp[(4 * i) + 5];
       bit_array[fft_bit_size - (4 * i) - 4] = bit_array_tmp[(4 * i) + 4];
	   bit_array[fft_bit_size - (4 * i) - 5] = bit_array_tmp[(4 * i) + 3];
	   bit_array[fft_bit_size - (4 * i) - 6] = bit_array_tmp[(4 * i) + 2];
	  }
    }
    
    for(int j = 0; j < fft_bit_size; j++){
		if(bit_array[j] == 1) weight_tmp = 1 << j;
		else weight_tmp = 0;
		result = result + weight_tmp;
    }
    
}
//-------
//NTT radix-16 and radix-2 Reindex for output
void NTTSPMB::NTT_REORDERINDEX_R16_R8_OUT(int index,int &result){
    int     fft_bit_size;
    int     index_tmp;
	int     bit_tmp;
	int     weight_tmp;
	int     Number_of_digit;

	std::vector<int> bit_array;
	std::vector<int> bit_array_tmp;
	
	//-----------------------------------------------------
    fft_bit_size = (int)ceil(log2(N));
	bit_array.resize(fft_bit_size);
	bit_array_tmp.resize(fft_bit_size);
	
    //Number of digit , four bits as one digit.
	Number_of_digit = (fft_bit_size - 3) / 4;
	//std::cout << "fft_bit_size: " << fft_bit_size <<"\n";
	//std::cout << "Number_of_digit: " << Number_of_digit <<"\n";
	//-----------------------------------------------------
	index_tmp = index;
	//bit array calculate
	for(int i = 0; i < fft_bit_size; i++){
        bit_tmp    = index_tmp % 2;
        index_tmp  = index_tmp >> 1;
        bit_array[i]     = bit_tmp;
        bit_array_tmp[i] = bit_tmp;	
	}
    
	result = 0;
	//MSB relocation , BC[MSB] = BC[0]	
	bit_array[fft_bit_size - 1] = bit_array_tmp[2];
	bit_array[fft_bit_size - 2] = bit_array_tmp[1];
	bit_array[fft_bit_size - 3] = bit_array_tmp[0];
    for(int i = 0; i < Number_of_digit; i++){
       bit_array[fft_bit_size - (4 * i) - 4] = bit_array_tmp[(4 * i) + 6];
       bit_array[fft_bit_size - (4 * i) - 5] = bit_array_tmp[(4 * i) + 5];
	   bit_array[fft_bit_size - (4 * i) - 6] = bit_array_tmp[(4 * i) + 4];
	   bit_array[fft_bit_size - (4 * i) - 7] = bit_array_tmp[(4 * i) + 3];
    }
    
    for(int j = 0; j < fft_bit_size; j++){
		if(bit_array[j] == 1) weight_tmp = 1 << j;
		else weight_tmp = 0;
		result = result + weight_tmp;
    }
    
}
//**************************************************************************
//IFFT                   Reorder BC  in IFFT computing 
//**************************************************************************
//radix-16  and radix-2 
void NTTSPMB::BC_IFFT_Reorder_R16_R2(int BC,int &result){
    int    BC_bit_size;
    int    BC_tmp;
	int    Number_of_digit;
	int    bit_tmp;
	int    weight_tmp;
	std::vector<int> bit_array;
	std::vector<int> bit_array_tmp;
	//parameter initial
    BC_bit_size = (int)ceil(log2(N/16));
	bit_array.resize(BC_bit_size);
	bit_array_tmp.resize(BC_bit_size);
	BC_tmp = BC;
    //bit calculate
	//std::cout << "---------------------------\n";
	//std::cout << "BC: " << BC <<"\n";
	for(int j=0; j < BC_bit_size;j++){
        bit_tmp = BC_tmp % 2;
        BC_tmp = BC_tmp >> 1;
        bit_array[j]     = bit_tmp;
        bit_array_tmp[j] = bit_tmp;
		//std::cout <<"bit_array_tmp["<<j<<"]:"<< bit_array_tmp[j] <<"\n";
    } 
    
	Number_of_digit  = (BC_bit_size - 1) / 4;
	//Orinal lsb ------> MSB
	//Others bit ,two bits as one digit ,then perform reverse operation
	bit_array[0] = bit_array_tmp[BC_bit_size-1];
	for(int j=0;j < Number_of_digit;j++){
		bit_array[4 * j + 1] = bit_array_tmp[BC_bit_size - 5 - (4 * j)];
		bit_array[4 * j + 2] = bit_array_tmp[BC_bit_size - 4 - (4 * j)];
		bit_array[4 * j + 3] = bit_array_tmp[BC_bit_size - 3 - (4 * j)];
		bit_array[4 * j + 4] = bit_array_tmp[BC_bit_size - 2 - (4 * j)];
	}
	
	result  = 0;
    for(int j=0; j < BC_bit_size;j++){
		//std::cout <<"bit_array["<<j<<"]:"<< bit_array[j] <<"\n";
		if(bit_array[j] == 1) weight_tmp = 1 << j;
		else weight_tmp = 0;
		result = result + weight_tmp;
    }
	
}
void NTTSPMB::BC_IFFT_Reorder_R16_R2_OUT(int BC,int &result){
    int    BC_bit_size;
    int    BC_tmp;
	int    bit_tmp;
	int    weight_tmp;
	std::vector<int> bit_array;
	std::vector<int> bit_array_tmp;
	
    BC_bit_size = (int)ceil(log2(N/16));
	bit_array.resize(BC_bit_size);
	bit_array_tmp.resize(BC_bit_size);

	//Frist , BC circular right shift 2 bit
	BC_tmp = BC;
	//bit calculate
	for(int j=0; j < BC_bit_size;j++){
        bit_tmp = BC_tmp % 2;
        BC_tmp = BC_tmp >> 1;
        bit_array[j]     = bit_tmp;
        bit_array_tmp[j] = bit_tmp;
    } 

	result  = 0;
	//for radix-16 , radix-2 mixed radix
	//for 8192 points, BC have 9 bits [8][7][6][5][4][3][2][1][0] -----> [7][6][5][3][2][1][0][4][8][2]
	//for 512 points, BC have 5 bits [4][3][2][1][0] -----> [2][1][0][4][3]
    // other points parameter, is not defined.
    if(BC_bit_size == 5){
		bit_array[4] = bit_array_tmp[3];
		bit_array[3] = bit_array_tmp[2];
		bit_array[2] = bit_array_tmp[0];
		bit_array[1] = bit_array_tmp[1];
		bit_array[0] = bit_array_tmp[4];
	}
	
	if(BC_bit_size == 9){
		bit_array[8] = bit_array_tmp[7];
		bit_array[7] = bit_array_tmp[6];
		bit_array[6] = bit_array_tmp[5];
		bit_array[5] = bit_array_tmp[0];
		bit_array[4] = bit_array_tmp[2];		
		bit_array[3] = bit_array_tmp[1];
		bit_array[2] = bit_array_tmp[4];
		bit_array[1] = bit_array_tmp[8];
		bit_array[0] = bit_array_tmp[3];				
	}

    for(int j=0; j < BC_bit_size;j++){
		//std::cout <<"bit_array["<<j<<"]:"<< bit_array[j] <<"\n";
		if(bit_array[j] == 1) weight_tmp = 1 << j;
		else weight_tmp = 0;
		result = result + weight_tmp;
		//std::cout << "bit_array[" << j <<"]: "<< bit_array[j] <<"\n";
    }
	
}
//radix-16 and radix-4
void NTTSPMB::BC_IFFT_Reorder_R16_R4(int BC,int &result){
    int    BC_bit_size;
    int    BC_tmp;
	int    Number_of_digit;
	int    bit_tmp;
	int    weight_tmp;
	std::vector<int> bit_array;
	std::vector<int> bit_array_tmp;
	//parameter initial
    BC_bit_size = (int)ceil(log2(N/16));
	bit_array.resize(BC_bit_size);
	bit_array_tmp.resize(BC_bit_size);
	BC_tmp = BC;
    //bit calculate
	//std::cout << "---------------------------\n";
	//std::cout << "BC: " << BC <<"\n";
	for(int j=0; j < BC_bit_size;j++){
        bit_tmp = BC_tmp % 2;
        BC_tmp = BC_tmp >> 1;
        bit_array[j]     = bit_tmp;
        bit_array_tmp[j] = bit_tmp;
		//std::cout <<"bit_array_tmp["<<j<<"]:"<< bit_array_tmp[j] <<"\n";
    } 
    
	Number_of_digit  = (BC_bit_size - 2) / 4;
	//Orinal lsb ------> MSB
	//Others bit ,two bits as one digit ,then perform reverse operation
	bit_array[1] = bit_array_tmp[BC_bit_size-1];
	bit_array[0] = bit_array_tmp[BC_bit_size-2];
	for(int j=0;j < Number_of_digit;j++){
		bit_array[4 * j + 2] = bit_array_tmp[BC_bit_size - 6 - (4 * j)];
		bit_array[4 * j + 3] = bit_array_tmp[BC_bit_size - 5 - (4 * j)];
		bit_array[4 * j + 4] = bit_array_tmp[BC_bit_size - 4 - (4 * j)];
		bit_array[4 * j + 5] = bit_array_tmp[BC_bit_size - 3 - (4 * j)];
	}
	
	result  = 0;
    for(int j=0; j < BC_bit_size;j++){
		//std::cout <<"bit_array["<<j<<"]:"<< bit_array[j] <<"\n";
		if(bit_array[j] == 1) weight_tmp = 1 << j;
		else weight_tmp = 0;
		result = result + weight_tmp;
    }
	
}
void NTTSPMB::BC_IFFT_Reorder_R16_R4_OUT(int BC,int &result){
    int    BC_bit_size;
    int    BC_tmp;
	int    bit_tmp;
	int    weight_tmp;
	std::vector<int> bit_array;
	std::vector<int> bit_array_tmp;
	
    BC_bit_size = (int)ceil(log2(N/16));
	bit_array.resize(BC_bit_size);
	bit_array_tmp.resize(BC_bit_size);

	//Frist , BC circular right shift 2 bit
	BC_tmp = BC;
	//bit calculate
	for(int j=0; j < BC_bit_size;j++){
        bit_tmp = BC_tmp % 2;
        BC_tmp = BC_tmp >> 1;
        bit_array[j]     = bit_tmp;
        bit_array_tmp[j] = bit_tmp;
    } 

	result  = 0;
	//for radix-16 , radix-2 mixed radix
	//for 8192 points, BC have 9 bits [8][7][6][5][4][3][2][1][0] -----> [7][6][5][3][2][1][0][4][8][2]
	//for 512 points, BC have 5 bits [4][3][2][1][0] -----> [2][1][0][4][3]
    // other points parameter, is not defined.
    if(BC_bit_size == 6){
		bit_array[5] = bit_array_tmp[1];
		bit_array[4] = bit_array_tmp[0];
		bit_array[3] = bit_array_tmp[3];
		bit_array[2] = bit_array_tmp[2];
		bit_array[1] = bit_array_tmp[5];
		bit_array[0] = bit_array_tmp[4];
	}


    if(BC_bit_size == 10){
		bit_array[9] = bit_array_tmp[7];
		bit_array[8] = bit_array_tmp[6];
		bit_array[7] = bit_array_tmp[1];
		bit_array[6] = bit_array_tmp[0];
		bit_array[5] = bit_array_tmp[5];
		bit_array[4] = bit_array_tmp[4];
		bit_array[3] = bit_array_tmp[9];
		bit_array[2] = bit_array_tmp[8];
		bit_array[1] = bit_array_tmp[3];
		bit_array[0] = bit_array_tmp[2];
	}	

    for(int j=0; j < BC_bit_size;j++){
		//std::cout <<"bit_array["<<j<<"]:"<< bit_array[j] <<"\n";
		if(bit_array[j] == 1) weight_tmp = 1 << j;
		else weight_tmp = 0;
		result = result + weight_tmp;
		//std::cout << "bit_array[" << j <<"]: "<< bit_array[j] <<"\n";
    }
	
}
//raidx-16 and radix-8
void NTTSPMB::BC_IFFT_Reorder_R16_R8(int BC,int &result){
    int    BC_bit_size;
    int    BC_tmp;
	int    Number_of_digit;
	int    bit_tmp;
	int    weight_tmp;
	std::vector<int> bit_array;
	std::vector<int> bit_array_tmp;
	//parameter initial
    BC_bit_size = (int)ceil(log2(N/16));
	bit_array.resize(BC_bit_size);
	bit_array_tmp.resize(BC_bit_size);
	BC_tmp = BC;
    //bit calculate
	//std::cout << "---------------------------\n";
	//std::cout << "BC: " << BC <<"\n";
	for(int j=0; j < BC_bit_size;j++){
        bit_tmp = BC_tmp % 2;
        BC_tmp = BC_tmp >> 1;
        bit_array[j]     = bit_tmp;
        bit_array_tmp[j] = bit_tmp;
		//std::cout <<"bit_array_tmp["<<j<<"]:"<< bit_array_tmp[j] <<"\n";
    } 
    
	Number_of_digit  = (BC_bit_size - 3) / 4;
	//Orinal lsb ------> MSB
	//Others bit ,two bits as one digit ,then perform reverse operation
	bit_array[2] = bit_array_tmp[BC_bit_size-1];
	bit_array[1] = bit_array_tmp[BC_bit_size-2];
	bit_array[0] = bit_array_tmp[BC_bit_size-3];
	for(int j=0;j < Number_of_digit;j++){
		bit_array[4 * j + 3] = bit_array_tmp[BC_bit_size - 7 - (4 * j)];
		bit_array[4 * j + 4] = bit_array_tmp[BC_bit_size - 6 - (4 * j)];
		bit_array[4 * j + 5] = bit_array_tmp[BC_bit_size - 5 - (4 * j)];
		bit_array[4 * j + 6] = bit_array_tmp[BC_bit_size - 4 - (4 * j)];
	}
	
	result  = 0;
    for(int j=0; j < BC_bit_size;j++){
		//std::cout <<"bit_array["<<j<<"]:"<< bit_array[j] <<"\n";
		if(bit_array[j] == 1) weight_tmp = 1 << j;
		else weight_tmp = 0;
		result = result + weight_tmp;
    }
	
}
void NTTSPMB::BC_IFFT_Reorder_R16_R8_OUT(int BC,int &result){
    int    BC_bit_size;
    int    BC_tmp;
	int    bit_tmp;
	int    weight_tmp;
	
	//------------
	std::vector<int> bit_array;
	std::vector<int> bit_array_tmp;
	//parameter initial
    BC_bit_size = (int)ceil(log2(N/16));
	bit_array.resize(BC_bit_size);
	bit_array_tmp.resize(BC_bit_size);
	BC_tmp = BC;
    //bit calculate
	//std::cout << "---------------------------\n";
	//std::cout << "BC: " << BC <<"\n";
	for(int j=0; j < BC_bit_size;j++){
        bit_tmp = BC_tmp % 2;
        BC_tmp = BC_tmp >> 1;
        bit_array[j]     = bit_tmp;
        bit_array_tmp[j] = bit_tmp;
    } 
    if(BC_bit_size == 7){
		bit_array[6]  = bit_array_tmp[0];
		bit_array[5]  = bit_array_tmp[6];
		bit_array[4]  = bit_array_tmp[5];
		bit_array[3]  = bit_array_tmp[4];
		bit_array[2]  = bit_array_tmp[1];
		bit_array[1]  = bit_array_tmp[3];
		bit_array[0]  = bit_array_tmp[2];
		
	}    

    if(BC_bit_size == 11){
		bit_array[10] = bit_array_tmp[7];
		bit_array[9]  = bit_array_tmp[2];
		bit_array[8]  = bit_array_tmp[1];
		bit_array[7]  = bit_array_tmp[0];
		bit_array[6]  = bit_array_tmp[4];
		bit_array[5]  = bit_array_tmp[10];
		bit_array[4]  = bit_array_tmp[9];
		bit_array[3]  = bit_array_tmp[8];
		bit_array[2]  = bit_array_tmp[3];
		bit_array[1]  = bit_array_tmp[6];
		bit_array[0]  = bit_array_tmp[5];
	}
	
	
	result  = 0;
    for(int j=0; j < BC_bit_size;j++){
		//std::cout <<"bit_array["<<j<<"]:"<< bit_array[j] <<"\n";
		if(bit_array[j] == 1) weight_tmp = 1 << j;
		else weight_tmp = 0;
		result = result + weight_tmp;
    }
	
}
//Gray code transform
int NTTSPMB::Gray(int index,int group){
	int result_tmp;
	int weight_tmp;
	int bit_tmp;
	int index_tmp;
	int group_bit_size;
	std::vector<int> bit_array;
	std::vector<int> bit_array_tmp;
	
	group_bit_size = (int) ceil(log2(group));
	bit_array.resize(group_bit_size);
	bit_array_tmp.resize(group_bit_size);
	
	index_tmp = index;
	for(int i = 0;i < group_bit_size ; i++){
        bit_tmp = index_tmp % 2;
        index_tmp = index_tmp >> 1;
        bit_array_tmp[i] = bit_tmp;		
	}
	
	result_tmp = 0;
	for(int i = 0;i < group_bit_size ;i++){
		if(i == (group_bit_size-1)) bit_array[i] = bit_array_tmp[i];
		else bit_array[i] = bit_array_tmp[i] ^ bit_array_tmp[i+1];
		
		if(bit_array[i] == 1) weight_tmp = 1 << i;
		else weight_tmp = 0;
		result_tmp = result_tmp + weight_tmp; 
	}
	
	return result_tmp;
}
void NTTSPMB::AGU_R2(int BC,int &BN,int &MA){
     int BC_WIDTH;
	 int BC_tmp;
	 int bit_tmp;
     
     BC_WIDTH = (unsigned long)ceil(log2(N/2));	 
	 MA     = BC / 2;
	 BC_tmp = BC;
	 
	 BN = 0;
	 for(int i = 0 ; i < BC_WIDTH;i++){
		bit_tmp = BC_tmp % 2;
		BC_tmp  = BC_tmp >> 1;
		BN = BN ^ bit_tmp; 
	 }	 
}
void NTTSPMB::AGU_R4(int BC,int &BN,int &MA){
     int BC_WIDTH;
	 int BC_tmp;
	 int bit_tmp;
     
     BC_WIDTH = (unsigned long)ceil(log2(N/4));	 
	 MA     = BC / 2;
	 BC_tmp = BC;
	 
	 BN = 0;
	 for(int i = 0 ; i < BC_WIDTH;i++){
		bit_tmp = BC_tmp % 2;
		BC_tmp  = BC_tmp >> 1;
		BN = BN ^ bit_tmp; 
	 }	 
}
void NTTSPMB::AGU_R16(int BC,int &BN,int &MA){
     int BC_WIDTH;
	 int BC_tmp;
	 int bit_tmp;
     
     BC_WIDTH = (unsigned long)ceil(log2(N/16));	 
	 MA     = BC / 2;
	 BC_tmp = BC;
	 
	 BN = 0;
	 for(int i = 0 ; i < BC_WIDTH;i++){
		bit_tmp = BC_tmp % 2;
		BC_tmp  = BC_tmp >> 1;
		BN = BN ^ bit_tmp; 
	 }	 
}

