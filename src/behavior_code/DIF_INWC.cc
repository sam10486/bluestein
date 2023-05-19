#include <iostream>
#include <fstream>
#include <NTL/ZZ.h>
#include <time.h>

#include "DIF_INWC.h"
#include "SPMB.h"

#include <vector>
#include <algorithm>
#include "BitOperate.h"
#include "DTFAG.h"

using namespace std;

void DIF_INWC::INWC_Radix2_BU(ZZ &a,ZZ &b){
	ZZ tmp_a;
	ZZ tmp_b;
	AddMod(tmp_a, a, b, p);
	if (b < 0)b = b + p;
	SubMod(tmp_b, a, b, p);
	a = tmp_a;
	b = tmp_b;
}

//radix 2^(2) DIF
void DIF_INWC::INWC_Radix4_BU(ZZ &a,ZZ &b,ZZ &c,ZZ &d){
	//a: x[n] b:[n+N/4] c:[n+N/2] d: 
	ZZ IW_1_4; // W^(N/4)
	ZZ twiddle_3_4; // W^(3N/4)
	unsigned long len_1_4; // N/4
	unsigned long len_3_4; // N/4
	len_1_4 = N / 4;
	PowerMod(IW_1_4, IW, len_1_4, p);
	
	ZZ tmp_a;
	ZZ tmp_b;
    ZZ tmp_c;
    ZZ tmp_d;
	
	//stage 0
	AddMod(tmp_a,a,c,p);
	SubMod(tmp_c,a,c,p);
	AddMod(tmp_b,b,d,p);
	SubMod(tmp_d,b,d,p);
	MulMod(tmp_d,tmp_d,IW_1_4,p);
	
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

void DIF_INWC::INWC_Radix16_BU(ZZ &a_r0,ZZ &a_r1,ZZ &a_r2,ZZ &a_r3,ZZ &a_r4,ZZ &a_r5,ZZ &a_r6,ZZ &a_r7,
ZZ &a_r8,ZZ &a_r9,ZZ &a_r10,ZZ &a_r11,ZZ &a_r12,ZZ &a_r13,ZZ &a_r14,ZZ &a_r15){
	ZZ IW_1_16; // W^(N/8)
	ZZ twiddle_r16_nk;
	ZZ data_tmp;
	unsigned long len_1_16; // N/16
	
	len_1_16 = N/16;
	PowerMod(IW_1_16,IW,len_1_16,p);
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
			PowerMod(twiddle_r16_nk,IW_1_16,nk_exp, p);
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


void DIF_INWC::DIF_INWC_radix2(std::vector<ZZ> &A){
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

	std::ofstream INWC_DATARECORD("./NWC_PrintData/INWC_R2_SPMB.txt");
	std::ofstream INWC_radix2("./NWC_PrintData/INWC_r2.txt");
	//-------------Bit Reverse--------------
	BitOperate RR, IntToVec, VecToInt;
	long long bit_width = log2(N);
	long long bit_width_s = log(radix)/log(2);
	vector<long long> bit_array;
	//----------------------------------------

    //------------DTFAG generator-------------
	DTFAG DTFAG;
	vector<ZZ > st0_Tw, st1_Tw, st2_Tw;
	st0_Tw.resize(radix);
	st1_Tw.resize(radix);
	st2_Tw.resize(radix);
	int DTFAG_t = 0;
	int DTFAG_i = 0;
	int DTFAG_j = 0;

	int fft_point = N;
	int radix_r1 = radix;
	int radix_r2 = radix;
	ZZ fft_twiddle = IW;	//***due to INWC, this ways use Inverse fft_twiddle***
	ZZ fft_prime = p;
	int debug = 0;
	vector<vector<ZZ > > ROM0;
    vector<vector<ZZ > > ROM1, ROM2;
    
    
	int arr_size = radix_r1 * radix_r1;
    ROM0.resize(radix_r1);
    for(int i=0; i<radix_r1; i++){
        ROM0[i].resize(radix_r1);
    }
    ROM1.resize(arr_size);
	for(int i=0; i<radix_r1; i++){
        ROM1[i].resize(radix_r1);
    }
    ROM2.resize(arr_size);
	for(int i=0; i<radix_r1; i++){
        ROM2[i].resize(radix_r1);
    }
	DTFAG.DTFAG_ROM_init(
        radix_r1, radix_r2, fft_twiddle, fft_prime, debug,
        ROM0, ROM1, ROM2);
	////----------------------------------------

	//-----------NWC PART-----------------------
	ZZ InvTwo;
	ZZ InvPhi_dot_IW, InvPhi_dot_IW_dot_InvTwo;
	InvMod(InvTwo, (ZZ)2, p);
	cout << "W = " << W << ", IW = " << IW <<  ", Phi = " << Phi << ", InvPhi = " << InvPhi << ", Inv_2 = " << InvTwo << endl;
	cout << "p = " << p << endl;
	//------------------------------------------

	
    Stage = (unsigned long)ceil(log2(N));
	BC_WIDTH  = (int)ceil(log2(N/2));
	offset    =  (int)N /radix;
	word_size =  (int)N / (2 * radix);
	group     =  (int)N / (radix * radix);
	bit_array_tmp.resize(BC_WIDTH);
	INWC_DATARECORD <<"Stage: "<< Stage<<"\n";

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
			if(bn_tmp == 0){
				A_B0R0[ma_tmp] = A[BC];
				A_B0R1[ma_tmp] = A[BC + offset];
			}else {
				A_B1R0[ma_tmp] = A[BC];
				A_B1R1[ma_tmp] = A[BC + offset];
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
		INWC_DATARECORD <<"Stage: "<< s<<"\n";
		INWC_DATARECORD << "factor "<< factor <<"\n";
		INWC_DATARECORD << "****************************\n";

		INWC_radix2 <<"Now Stage: "<< s <<"\n";
		INWC_radix2 <<"twiddle factor : "<< factor <<"\n";
		for(int i = 0 ;i < group;i++){
			INWC_radix2 << "----------------i =" << i << " ----------------" << std::endl;
			INWC_DATARECORD << "----------------i =" << i << " ----------------" << std::endl;
			for(int j = 0;j < radix;j++){
				INWC_DATARECORD << "---j =" << j << " ---" << std::endl;
				INWC_DATARECORD <<"twiddle factor : "<< factor <<"\n";
				INWC_DATARECORD <<"p : "<< p <<"\n";
	    		INWC_DATARECORD << "********\n";
				gray_i  = Gray(i,group);
				BC_tmp  = j * group + gray_i;
				INWC_DATARECORD << "BC_tmp = " << BC_tmp << ", BC = " << BC << ", s = " << s << endl;
				RR_R2(BC_tmp,s,BC);
				length = BC_tmp >> s;
				INWC_DATARECORD << "length: " <<  length <<"\n";
				PowerMod(factor_t,factor,length,p);
				INWC_DATARECORD << "factor_t: "<<factor_t<<"\n";

				AGU_R2(BC,bn_tmp,ma_tmp);

				//-----------compute data idx-------------
				INWC_radix2 << "BC_tmp = " << BC_tmp << std::endl;
				
				IntToVec.IntToVec(BC_tmp, N, bit_array);
				rotate(bit_array.begin(), bit_array.begin()+ bit_width_s*s , bit_array.end());
				long long Data = VecToInt.VecToInt(bit_array, N);
				INWC_radix2 << "Data_index = ";
                INWC_radix2 << "( " ;					
				for(int k = 0; k < radix ; k++ ){
					INWC_radix2 << Data + k*(1<<(bit_width-bit_width_s-bit_width_s*s)) <<" ";	
				}
				INWC_radix2 << ") " ;
				INWC_radix2 << ", (w^" << 0 << ", w^" << tw_degree * length << ")" <<std::endl;
				//-----------------------------------------
				
                //-----------DTFAG generator-------------
                DTFAG.DTFAG_SPMB_DIF_MR(
                    s, fft_point, radix_r1, radix_r2, debug,
                    ROM0, ROM1, ROM2,
					st0_Tw, st1_Tw, st2_Tw,
                    DTFAG_i, DTFAG_t, DTFAG_j);
                /*switch(s){
					case 0:
						for(int i=0; i<radix; i++){
							cout << "st0_Tw[" << i << "] = w^" << st0_Tw[i] << endl;
						}
						break;
					case 1:
						for(int i=0; i<radix; i++){
							cout << "st1_Tw[" << i << "] = w^" << st1_Tw[i] << endl;
						}
						break;
					case 2:
						for(int i=0; i<radix; i++){
							cout << "st2_Tw[" << i << "] = w^" << st2_Tw[i] << endl;
						}
						break;
				}*/
                //cout << "stage = " << s << ", DTFAG_i = " << DTFAG_i << ", DTFAG_t = " << DTFAG_t << ", DTFAG_j = " << DTFAG_j << endl;
                if(DTFAG_i == radix_r2-1 && DTFAG_j == radix_r1-1 && DTFAG_t == radix_r1-1){
					DTFAG_i = 0;
				}else if(DTFAG_t == radix_r1-1 && DTFAG_j == radix_r1-1){
					DTFAG_i++;
				}
				if(DTFAG_j == radix_r1-1 && DTFAG_t == radix_r1-1){
					DTFAG_t = 0;
				}else if(DTFAG_j == radix_r1-1){
					DTFAG_t++;
				}
				if(DTFAG_j == radix_r1-1){
					DTFAG_j = 0;
				}else{
					DTFAG_j++;
				}
                //---------------------------------------
				//---------NWC PART-------------
				ZZ InvPhi_deg = PowerMod((ZZ)2, s, p);
				ZZ InvPhi_Order = PowerMod(InvPhi, InvPhi_deg, p);
				//------------------------------
				if(bn_tmp == 0){
					INWC_DATARECORD << "bn_tmp = " << bn_tmp << ", ma_tmp = " << ma_tmp << std::endl;
					bn0_bc_tmp = BC_tmp;
                    switch(s){
                        case 0:
                            INWC_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp] << ", st0_Tw[0] = " << st0_Tw[0] << endl;
					        INWC_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp] << ", st0_Tw[1] = " << st0_Tw[1] << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_Order = " << InvPhi_Order << endl;
							if(debug) INWC_DATARECORD << "InvPhi_Order = " << InvPhi_deg << endl;
					        if(!debug) INWC_Radix2_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp]);
							// upper
                            if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],st0_Tw[0],p);
					        if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],InvTwo,p);
							// down
							//if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],st0_Tw[1],p);
							if(!debug) MulMod(InvPhi_dot_IW, InvPhi_Order, st0_Tw[1], p);
							if(!debug) MulMod(InvPhi_dot_IW_dot_InvTwo, InvPhi_dot_IW, InvTwo, p);
							if(!debug) MulMod(A_B0R1[ma_tmp], A_B0R1[ma_tmp], InvPhi_dot_IW_dot_InvTwo, p);
                            break;
                        case 1:
                            INWC_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp] << ", st1_Tw[0] = " << st1_Tw[0] << endl;
					        INWC_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp] << ", st1_Tw[1] = " << st1_Tw[1] << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_Order = " << InvPhi_Order << endl;
							if(debug) INWC_DATARECORD << "InvPhi_Order = " << InvPhi_deg << endl;
					        if(!debug) INWC_Radix2_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp]);
							// upper
                            if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],st1_Tw[0],p);
							if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],InvTwo,p);
							// down
					        //if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],st1_Tw[1],p);
							if(!debug) MulMod(InvPhi_dot_IW, InvPhi_Order, st1_Tw[1], p);
							if(!debug) MulMod(InvPhi_dot_IW_dot_InvTwo, InvPhi_dot_IW, InvTwo, p);
							if(!debug) MulMod(A_B0R1[ma_tmp], A_B0R1[ma_tmp], InvPhi_dot_IW_dot_InvTwo, p);
                            break;
                        case 2:
                            INWC_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp] << ", st2_Tw[0] = " << st2_Tw[0] << endl;
					        INWC_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp] << ", st2_Tw[1] = " << st2_Tw[1] << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_Order = " << InvPhi_Order << endl;
							if(debug) INWC_DATARECORD << "InvPhi_Order = " << InvPhi_deg << endl;
					        if(!debug) INWC_Radix2_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp]);
							// upper
                            if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],st2_Tw[0],p);
							if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],InvTwo,p);
							// down
					        //if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],st2_Tw[1],p);
							if(!debug) MulMod(InvPhi_dot_IW, InvPhi_Order, st2_Tw[1], p);
							if(!debug) MulMod(InvPhi_dot_IW_dot_InvTwo, InvPhi_dot_IW, InvTwo, p);
							if(!debug) MulMod(A_B0R1[ma_tmp], A_B0R1[ma_tmp], InvPhi_dot_IW_dot_InvTwo, p);
                            break;
                        case 3:
                            INWC_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp] << ", st3_Tw[0] = " << 1 << endl;
					        INWC_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp] << ", st3_Tw[1] = " << 1 << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_Order = " << InvPhi_Order << endl;
							if(debug) INWC_DATARECORD << "InvPhi_Order = " << InvPhi_deg << endl;
					        if(!debug) INWC_Radix2_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp]);
							// upper
                            if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],1,p);
							if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],InvTwo,p);
							// down
					        //if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],1,p);	
							if(!debug) MulMod(InvPhi_dot_IW, InvPhi_Order, 1, p);
							if(!debug) MulMod(InvPhi_dot_IW_dot_InvTwo, InvPhi_dot_IW, InvTwo, p);
							if(!debug) MulMod(A_B0R1[ma_tmp], A_B0R1[ma_tmp], InvPhi_dot_IW_dot_InvTwo, p);	
                            break;
                    }
                    
					INWC_DATARECORD << "---after BU compute---" << std::endl;
				    INWC_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp]<<"\n";
					INWC_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp]<<"\n";
					if(!debug) INWC_DATARECORD << "InvPhi_dot_IW = " << InvPhi_dot_IW << endl;
					bn0_ma_reg = ma_tmp;

				}
			    else {
					INWC_DATARECORD << "bn_tmp = " << bn_tmp << ", ma_tmp = " << ma_tmp << std::endl;
					bn1_bc_tmp = BC_tmp;
                    switch(s){
                        case 0:
                            INWC_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<< ", st0_Tw[0] = " << st0_Tw[0] << endl;
					        INWC_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<< ", st0_Tw[1] = " << st0_Tw[1] << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_Order = " << InvPhi_Order << endl;
							if(debug) INWC_DATARECORD << "InvPhi_Order = " << InvPhi_deg << endl;
					        if(!debug) INWC_Radix2_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp]);
							// upper
                            if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],st0_Tw[0],p);
							if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],InvTwo,p);
							// down
					        //if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],st0_Tw[1],p);  
							if(!debug) MulMod(InvPhi_dot_IW, InvPhi_Order, st0_Tw[1], p);
							if(!debug) MulMod(InvPhi_dot_IW_dot_InvTwo, InvPhi_dot_IW, InvTwo, p);
							if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],InvPhi_dot_IW_dot_InvTwo,p);
                            break;
                        case 1:
                            INWC_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<< ", st1_Tw[0] = " << st1_Tw[0] << endl;
					        INWC_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<< ", st1_Tw[1] = " << st1_Tw[1] << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_Order = " << InvPhi_Order << endl;
							if(debug) INWC_DATARECORD << "InvPhi_Order = " << InvPhi_deg << endl;
					        if(!debug) INWC_Radix2_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp]);
							// upper
                            if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],st1_Tw[0],p);
							if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],InvTwo,p);
							// down
					        //if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],st1_Tw[1],p);
							if(!debug) MulMod(InvPhi_dot_IW, InvPhi_Order, st1_Tw[1], p);
							if(!debug) MulMod(InvPhi_dot_IW_dot_InvTwo, InvPhi_dot_IW, InvTwo, p);
							if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],InvPhi_dot_IW_dot_InvTwo,p);
                            break;
                        case 2:
                            INWC_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<< ", st2_Tw[0] = " << st2_Tw[0] << endl;
					        INWC_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<< ", st2_Tw[1] = " << st2_Tw[1] << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_Order = " << InvPhi_Order << endl;
							if(debug) INWC_DATARECORD << "InvPhi_Order = " << InvPhi_deg << endl;
					        if(!debug) INWC_Radix2_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp]);
							// upper
                            if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],st2_Tw[0],p);
							if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],InvTwo,p);
							// down
					        //if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],st2_Tw[1],p);
							if(!debug) MulMod(InvPhi_dot_IW, InvPhi_Order, st2_Tw[1], p);
							if(!debug) MulMod(InvPhi_dot_IW_dot_InvTwo, InvPhi_dot_IW, InvTwo, p);
							if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],InvPhi_dot_IW_dot_InvTwo,p);
                            break;
                        case 3:
                            INWC_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<< ", st3_Tw[0] = " << 1 << endl;
					        INWC_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<< ", st3_Tw[1] = " << 1 << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_Order = " << InvPhi_Order << endl;
							if(debug) INWC_DATARECORD << "InvPhi_Order = " << InvPhi_deg << endl;
					        if(!debug) INWC_Radix2_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp]);
							// upper
                            if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],1,p);
							if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],InvTwo,p);
							// down
					        if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],1,p);
							if(!debug) MulMod(InvPhi_dot_IW, InvPhi_Order, 1, p);
							if(!debug) MulMod(InvPhi_dot_IW_dot_InvTwo, InvPhi_dot_IW, InvTwo, p);
							if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],InvPhi_dot_IW_dot_InvTwo,p);
                            break;
                    }
					INWC_DATARECORD << "---after BU compute---" << std::endl;
				    INWC_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<<"\n";
                    INWC_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<<"\n";
					if(!debug) INWC_DATARECORD << "InvPhi_dot_IW = " << InvPhi_dot_IW << endl;
					bn1_ma_reg = ma_tmp;
				}
				INWC_DATARECORD <<"--------------------------------------------------------------------\n";
			}
		    //data relocation
		    if(s < Stage-1){
		    	if(bn1_bc_tmp > bn0_bc_tmp){
		    		INWC_DATARECORD << "bn0_bc_tmp = " << bn0_bc_tmp << ", bn1_bc_tmp = " << bn1_bc_tmp << std::endl;
		    		data_tmp = A_B0R1[bn0_ma_reg];
		    		A_B0R1[bn0_ma_reg] = A_B1R0[bn1_ma_reg];
		    		A_B1R0[bn1_ma_reg] = data_tmp;
		    	}else {
		    		INWC_DATARECORD << "bn0_bc_tmp = " << bn0_bc_tmp << ", bn1_bc_tmp = " << bn1_bc_tmp << std::endl;
		    		data_tmp = A_B1R1[bn1_ma_reg];
		    		A_B1R1[bn1_ma_reg] = A_B0R0[bn0_ma_reg];
		    		A_B0R0[bn0_ma_reg] = data_tmp;
		    	}
		    }
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
			   BR_R2(2 * BC_tmp,index0);
			   BR_R2(2 * BC_tmp + 1,index1);
               A[index0]     = A_B0R0[ma_tmp];
			   A[index1] = A_B0R1[ma_tmp];
			}
			else {
			   BR_R2(2 * BC_tmp,index0);
			   BR_R2(2 * BC_tmp + 1,index1);
               A[index0]     = A_B1R0[ma_tmp];
               A[index1] = A_B1R1[ma_tmp];
			}			
		}		
	}
}

void DIF_INWC::DIF_INWC_radix4(std::vector<ZZ> &A){
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
	
	std::ofstream INWC_DATARECORD("./NWC_PrintData/INWC_R4_SPMB.txt");
	std::ofstream INWC_radix4("./NWC_PrintData/INWC_radix4.txt");
	//-------------Bit Reverse--------------
	BitOperate RR, IntToVec, VecToInt, BR;
	long long bit_width = log2(N);
	long long bit_width_s = log(radix)/log(2);
	vector<long long> bit_array;
	//----------------------------------------

    //------------DTFAG generator-------------
	DTFAG DTFAG;
	vector<ZZ > st0_Tw, st1_Tw, st2_Tw;
	st0_Tw.resize(radix);
	st1_Tw.resize(radix);
	st2_Tw.resize(radix);
	int DTFAG_t = 0;
	int DTFAG_i = 0;
	int DTFAG_j = 0;

	int fft_point = N;
	int radix_r1 = radix;
	int radix_r2 = radix;
	ZZ fft_twiddle = IW;  //***due to INWC, this ways use Inverse fft_twiddle***
	ZZ fft_prime = p;
	int debug = 0;
	vector<vector<ZZ > > ROM0;
    vector<vector<ZZ > > ROM1, ROM2;
    
    
	int arr_size = radix_r1 * radix_r1;
    ROM0.resize(radix_r1);
    for(int i=0; i<radix_r1; i++){
        ROM0[i].resize(radix_r1);
    }
    ROM1.resize(arr_size);
	for(int i=0; i<radix_r1; i++){
        ROM1[i].resize(radix_r1);
    }
    ROM2.resize(arr_size);
	for(int i=0; i<radix_r1; i++){
        ROM2[i].resize(radix_r1);
    }
	DTFAG.DTFAG_ROM_init(
        radix_r1, radix_r2, fft_twiddle, fft_prime, debug,
        ROM0, ROM1, ROM2);
	////----------------------------------------

	//-----------NWC PART-----------------------
	ZZ InvFour;
	ZZ InvPhi_0t_dot_IW, InvPhi_0t_dot_IW_dot_InvTwo;
	ZZ InvPhi_1t_dot_IW, InvPhi_1t_dot_IW_dot_InvTwo;
	ZZ InvPhi_2t_dot_IW, InvPhi_2t_dot_IW_dot_InvTwo;
	ZZ InvPhi_3t_dot_IW, InvPhi_3t_dot_IW_dot_InvTwo;
	InvMod(InvFour, (ZZ)4, p);
	cout << "W = " << W << ", IW = " << IW <<  ", Phi = " << Phi << ", InvPhi = " << InvPhi << ", Inv_4 = " << InvFour << endl;
	cout << "p = " << p << endl;
	//------------------------------------------
	
    Stage = (unsigned long)ceil(log2(N));
	Stage = (unsigned long)Stage / 2;
	BC_WIDTH  = (int)ceil(log2(N/4));
	offset    =  (int)N /radix;
	word_size =  (int)N / (2 * radix);
	group     =  (int)N / (radix * radix);
	bit_array_tmp.resize(BC_WIDTH);
	INWC_DATARECORD <<"Stage: "<< Stage<<"\n";
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
		INWC_DATARECORD <<"---------------------------------\n";
		INWC_DATARECORD <<"Now Stage: "<< s <<"\n";
		

		INWC_radix4 <<"Now Stage: "<< s <<"\n";
		INWC_radix4 <<"twiddle factor : "<< factor <<"\n";
	    INWC_radix4 << "********\n";
		for(int i = 0 ;i < group;i++){
			INWC_DATARECORD <<"twiddle factor : "<< factor <<"\n";
			INWC_DATARECORD <<"p : "<< p <<"\n";
			INWC_DATARECORD << "W = " << W << ", IW = " << IW <<  ", Phi = " << Phi << ", InvPhi = " << InvPhi << ", Inv_2 = " << InvFour << endl;
			INWC_DATARECORD << "p = " << p << endl;
	    	INWC_DATARECORD << "********\n";
			bn0_bc_tmp  = 0;
			bn1_bc_tmp  = 0;
			INWC_radix4 <<"--------------i = " << i << "----------------\n";
			for(int j = 0;j < radix;j++){	
				gray_i  = Gray(i,group);
				BC_tmp  = j * group + gray_i;
				INWC_DATARECORD << "i: " << i <<"\n";
				INWC_DATARECORD << "gray_i: " << gray_i <<"\n";
				INWC_DATARECORD << "BC_tmp: " << BC_tmp <<"\n";
				RR_R4(BC_tmp,s,BC);
				INWC_DATARECORD << "BC: " << BC <<"\n";
				length = BC_tmp >> (2*s);
				INWC_DATARECORD << "length: " <<  length <<"\n";

				//-----------compute data idx-------------
				INWC_radix4 << "BC_tmp = " << BC_tmp << std::endl;
				
				IntToVec.IntToVec(BC_tmp, N, bit_array);
				rotate(bit_array.begin(), bit_array.begin()+ bit_width_s*s , bit_array.end());
				long long Data = VecToInt.VecToInt(bit_array, N);
				INWC_radix4 << "Data_index = ";
                INWC_radix4 << "( " ;		
				int spmb_radix4_arr[radix];		
				for(long long k = 0; k < radix ; k++ ){
					long long BR_tmp = BR.BitReserve(k, log2(radix));
					spmb_radix4_arr[BR_tmp] = Data + k*(1<<(bit_width-bit_width_s-bit_width_s*s));	
				}
				for(int k=0; k<radix; k++) {
					INWC_radix4 << spmb_radix4_arr[k] << ", ";
				}
				INWC_radix4 << ") " ;
				INWC_radix4 << ", (w^" << 0 << ", w^" << tw_degree * length 
				<< ", w^"  << tw_degree * length * 2 << ", w^" << tw_degree * length * 3 << ")" << std::endl;
				//-----------------------------------------
				
				PowerMod(factor_t,factor,length,p);
				INWC_DATARECORD << "factor_t: "<<factor_t<<"\n";
				AGU_R4(BC,bn_tmp,ma_tmp);
				INWC_DATARECORD << "bn_tmp: "<<bn_tmp<<"\n";
				//-----------DTFAG generator-------------
                DTFAG.DTFAG_SPMB_DIF_MR(
                    s, fft_point, radix_r1, radix_r2, debug,
					ROM0, ROM1, ROM2,
                    st0_Tw, st1_Tw, st2_Tw,
                    DTFAG_i, DTFAG_t, DTFAG_j);
                /*switch(s){
					case 0:
						for(int i=0; i<radix; i++){
							cout << "st0_Tw[" << i << "] = w^" << st0_Tw[i] << endl;
						}
						break;
					case 1:
						for(int i=0; i<radix; i++){
							cout << "st1_Tw[" << i << "] = w^" << st1_Tw[i] << endl;
						}
						break;
					case 2:
						for(int i=0; i<radix; i++){
							cout << "st2_Tw[" << i << "] = w^" << st2_Tw[i] << endl;
						}
						break;
				}*/
                //cout << "stage = " << s << ", DTFAG_i = " << DTFAG_i << ", DTFAG_t = " << DTFAG_t << ", DTFAG_j = " << DTFAG_j << endl;
                if(DTFAG_i == radix_r2-1 && DTFAG_j == radix_r1-1 && DTFAG_t == radix_r1-1){
					DTFAG_i = 0;
				}else if(DTFAG_t == radix_r1-1 && DTFAG_j == radix_r1-1){
					DTFAG_i++;
				}
				if(DTFAG_j == radix_r1-1 && DTFAG_t == radix_r1-1){
					DTFAG_t = 0;
				}else if(DTFAG_j == radix_r1-1){
					DTFAG_t++;
				}
				if(DTFAG_j == radix_r1-1){
					DTFAG_j = 0;
				}else{
					DTFAG_j++;
				}
                //---------------------------------------
                //---------NWC PART-------------
				ZZ InvPhi_0t, InvPhi_1t, InvPhi_2t, InvPhi_3t;
				ZZ InvPhi_0t_Order, InvPhi_1t_Order, InvPhi_2t_Order, InvPhi_3t_Order;
				ZZ InvPhi_deg = PowerMod((ZZ)4, s, p);
				InvPhi_0t = PowerMod(InvPhi, 0, p);
				InvPhi_1t = PowerMod(InvPhi, 1, p);
				InvPhi_2t = PowerMod(InvPhi, 2, p);
				InvPhi_3t = PowerMod(InvPhi, 3, p);
				InvPhi_0t_Order = PowerMod(InvPhi_0t, InvPhi_deg, p);
				InvPhi_1t_Order = PowerMod(InvPhi_1t, InvPhi_deg, p);
				InvPhi_2t_Order = PowerMod(InvPhi_2t, InvPhi_deg, p);
				InvPhi_3t_Order = PowerMod(InvPhi_3t, InvPhi_deg, p);
				//------------------------------
                if(bn_tmp == 0){
					if(j < 2)bn0_bc_tmp = BC_tmp;
					PowerMod(factor_2t,factor_t,2,p);
					PowerMod(factor_3t,factor_t,3,p);
                    switch(s){
                        case 0:
                            INWC_DATARECORD <<"Before butterfly unit operation! \n";
                            INWC_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: " <<A_B0R0[ma_tmp] << ", st0_Tw[0] = " << st0_Tw[0] << endl;
					        INWC_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: " <<A_B0R1[ma_tmp] << ", st0_Tw[1] = " << st0_Tw[1] << endl;
					        INWC_DATARECORD <<"A_B0R2["<<ma_tmp<<"]: " <<A_B0R2[ma_tmp] << ", st0_Tw[2] = " << st0_Tw[2] << endl;
					        INWC_DATARECORD <<"A_B0R3["<<ma_tmp<<"]: " <<A_B0R3[ma_tmp] << ", st0_Tw[3] = " << st0_Tw[3] << endl;   
							if(!debug) INWC_DATARECORD << "InvPhi_0t_Order = " << InvPhi_0t_Order << endl; else INWC_DATARECORD << "InvPhi_0t_Order = " << InvPhi_deg*0 << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_1t_Order = " << InvPhi_1t_Order << endl; else INWC_DATARECORD << "InvPhi_1t_Order = " << InvPhi_deg*1 << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_2t_Order = " << InvPhi_2t_Order << endl; else INWC_DATARECORD << "InvPhi_2t_Order = " << InvPhi_deg*2 << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_3t_Order = " << InvPhi_3t_Order << endl; else INWC_DATARECORD << "InvPhi_3t_Order = " << InvPhi_deg*3 << endl;   
							if(!debug) INWC_Radix4_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp],A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
					        INWC_DATARECORD <<"-------------------" << std::endl;
							//------------compute for INWC---------------
							if(!debug) MulMod(InvPhi_0t_dot_IW, InvPhi_0t_Order, st0_Tw[0], p);
							if(!debug) MulMod(InvPhi_1t_dot_IW, InvPhi_1t_Order, st0_Tw[1], p);
							if(!debug) MulMod(InvPhi_2t_dot_IW, InvPhi_2t_Order, st0_Tw[2], p);
							if(!debug) MulMod(InvPhi_3t_dot_IW, InvPhi_3t_Order, st0_Tw[3], p);

							if(!debug) MulMod(InvPhi_0t_dot_IW_dot_InvTwo, InvPhi_0t_dot_IW, InvFour, p);
							if(!debug) MulMod(InvPhi_1t_dot_IW_dot_InvTwo, InvPhi_1t_dot_IW, InvFour, p);
							if(!debug) MulMod(InvPhi_2t_dot_IW_dot_InvTwo, InvPhi_2t_dot_IW, InvFour, p);
							if(!debug) MulMod(InvPhi_3t_dot_IW_dot_InvTwo, InvPhi_3t_dot_IW, InvFour, p);
							//-------------------------------------------
                            if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],InvPhi_0t_dot_IW_dot_InvTwo,p);
                            if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],InvPhi_1t_dot_IW_dot_InvTwo,p);
					        if(!debug) MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],InvPhi_2t_dot_IW_dot_InvTwo,p);
					        if(!debug) MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],InvPhi_3t_dot_IW_dot_InvTwo,p);			
                            break;
                        case 1:
                            INWC_DATARECORD <<"Before butterfly unit operation! \n";
					        INWC_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: " <<A_B0R0[ma_tmp] << ", st1_Tw[0] = " << st1_Tw[0] << endl;
					        INWC_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: " <<A_B0R1[ma_tmp] << ", st1_Tw[1] = " << st1_Tw[1] << endl;
					        INWC_DATARECORD <<"A_B0R2["<<ma_tmp<<"]: " <<A_B0R2[ma_tmp] << ", st1_Tw[2] = " << st1_Tw[2] << endl;
					        INWC_DATARECORD <<"A_B0R3["<<ma_tmp<<"]: " <<A_B0R3[ma_tmp] << ", st1_Tw[3] = " << st1_Tw[3] << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_0t_Order = " << InvPhi_0t_Order << endl; else INWC_DATARECORD << "InvPhi_0t_Order = " << InvPhi_deg*0 << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_1t_Order = " << InvPhi_1t_Order << endl; else INWC_DATARECORD << "InvPhi_1t_Order = " << InvPhi_deg*1 << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_2t_Order = " << InvPhi_2t_Order << endl; else INWC_DATARECORD << "InvPhi_2t_Order = " << InvPhi_deg*2 << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_3t_Order = " << InvPhi_3t_Order << endl; else INWC_DATARECORD << "InvPhi_3t_Order = " << InvPhi_deg*3 << endl;   
					        if(!debug) INWC_Radix4_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp],A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
					        INWC_DATARECORD <<"-------------------" << std::endl;
							//------------compute for INWC---------------
							if(!debug) MulMod(InvPhi_0t_dot_IW, InvPhi_0t_Order, st1_Tw[0], p);
							if(!debug) MulMod(InvPhi_1t_dot_IW, InvPhi_1t_Order, st1_Tw[1], p);
							if(!debug) MulMod(InvPhi_2t_dot_IW, InvPhi_2t_Order, st1_Tw[2], p);
							if(!debug) MulMod(InvPhi_3t_dot_IW, InvPhi_3t_Order, st1_Tw[3], p);

							if(!debug) MulMod(InvPhi_0t_dot_IW_dot_InvTwo, InvPhi_0t_dot_IW, InvFour, p);
							if(!debug) MulMod(InvPhi_1t_dot_IW_dot_InvTwo, InvPhi_1t_dot_IW, InvFour, p);
							if(!debug) MulMod(InvPhi_2t_dot_IW_dot_InvTwo, InvPhi_2t_dot_IW, InvFour, p);
							if(!debug) MulMod(InvPhi_3t_dot_IW_dot_InvTwo, InvPhi_3t_dot_IW, InvFour, p);
							//-------------------------------------------
                            if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],InvPhi_0t_dot_IW_dot_InvTwo,p);
                            if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],InvPhi_1t_dot_IW_dot_InvTwo,p);
					        if(!debug) MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],InvPhi_2t_dot_IW_dot_InvTwo,p);
					        if(!debug) MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],InvPhi_3t_dot_IW_dot_InvTwo,p);
                            break;
                        case 2:
                            INWC_DATARECORD <<"Before butterfly unit operation! \n";
					        INWC_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: " <<A_B0R0[ma_tmp] << ", st2_Tw[0] = " << st2_Tw[0] << endl;
					        INWC_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: " <<A_B0R1[ma_tmp] << ", st2_Tw[1] = " << st2_Tw[1] << endl;
					        INWC_DATARECORD <<"A_B0R2["<<ma_tmp<<"]: " <<A_B0R2[ma_tmp] << ", st2_Tw[2] = " << st2_Tw[2] << endl;
					        INWC_DATARECORD <<"A_B0R3["<<ma_tmp<<"]: " <<A_B0R3[ma_tmp] << ", st2_Tw[3] = " << st2_Tw[3] << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_0t_Order = " << InvPhi_0t_Order << endl; else INWC_DATARECORD << "InvPhi_0t_Order = " << InvPhi_deg*0 << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_1t_Order = " << InvPhi_1t_Order << endl; else INWC_DATARECORD << "InvPhi_1t_Order = " << InvPhi_deg*1 << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_2t_Order = " << InvPhi_2t_Order << endl; else INWC_DATARECORD << "InvPhi_2t_Order = " << InvPhi_deg*2 << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_3t_Order = " << InvPhi_3t_Order << endl; else INWC_DATARECORD << "InvPhi_3t_Order = " << InvPhi_deg*3 << endl;
					        if(!debug) INWC_Radix4_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp],A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
					        INWC_DATARECORD <<"-------------------" << std::endl;
							//------------compute for INWC---------------
							if(!debug) MulMod(InvPhi_0t_dot_IW, InvPhi_0t_Order, st2_Tw[0], p);
							if(!debug) MulMod(InvPhi_1t_dot_IW, InvPhi_1t_Order, st2_Tw[1], p);
							if(!debug) MulMod(InvPhi_2t_dot_IW, InvPhi_2t_Order, st2_Tw[2], p);
							if(!debug) MulMod(InvPhi_3t_dot_IW, InvPhi_3t_Order, st2_Tw[3], p);

							if(!debug) MulMod(InvPhi_0t_dot_IW_dot_InvTwo, InvPhi_0t_dot_IW, InvFour, p);
							if(!debug) MulMod(InvPhi_1t_dot_IW_dot_InvTwo, InvPhi_1t_dot_IW, InvFour, p);
							if(!debug) MulMod(InvPhi_2t_dot_IW_dot_InvTwo, InvPhi_2t_dot_IW, InvFour, p);
							if(!debug) MulMod(InvPhi_3t_dot_IW_dot_InvTwo, InvPhi_3t_dot_IW, InvFour, p);
							//-------------------------------------------
                            if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],InvPhi_0t_dot_IW_dot_InvTwo,p);
                            if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],InvPhi_1t_dot_IW_dot_InvTwo,p);
					        if(!debug) MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],InvPhi_2t_dot_IW_dot_InvTwo,p);
					        if(!debug) MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],InvPhi_3t_dot_IW_dot_InvTwo,p);
                            break;
                        case 3:
                            INWC_DATARECORD <<"Before butterfly unit operation! \n";
					        INWC_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: " <<A_B0R0[ma_tmp] << ", st3_Tw[0] = " << 1 << endl;
					        INWC_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: " <<A_B0R1[ma_tmp] << ", st3_Tw[1] = " << 1 << endl;
					        INWC_DATARECORD <<"A_B0R2["<<ma_tmp<<"]: " <<A_B0R2[ma_tmp] << ", st3_Tw[2] = " << 1 << endl;
					        INWC_DATARECORD <<"A_B0R3["<<ma_tmp<<"]: " <<A_B0R3[ma_tmp] << ", st3_Tw[3] = " << 1 << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_0t_Order = " << InvPhi_0t_Order << endl; else INWC_DATARECORD << "InvPhi_0t_Order = " << InvPhi_deg*0 << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_1t_Order = " << InvPhi_1t_Order << endl; else INWC_DATARECORD << "InvPhi_1t_Order = " << InvPhi_deg*1 << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_2t_Order = " << InvPhi_2t_Order << endl; else INWC_DATARECORD << "InvPhi_2t_Order = " << InvPhi_deg*2 << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_3t_Order = " << InvPhi_3t_Order << endl; else INWC_DATARECORD << "InvPhi_3t_Order = " << InvPhi_deg*3 << endl;
					        if(!debug) INWC_Radix4_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp],A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
					        INWC_DATARECORD <<"-------------------" << std::endl;
							//------------compute for INWC---------------
							if(!debug) MulMod(InvPhi_0t_dot_IW, InvPhi_0t_Order, 1, p);
							if(!debug) MulMod(InvPhi_1t_dot_IW, InvPhi_1t_Order, 1, p);
							if(!debug) MulMod(InvPhi_2t_dot_IW, InvPhi_2t_Order, 1, p);
							if(!debug) MulMod(InvPhi_3t_dot_IW, InvPhi_3t_Order, 1, p);

							if(!debug) MulMod(InvPhi_0t_dot_IW_dot_InvTwo, InvPhi_0t_dot_IW, InvFour, p);
							if(!debug) MulMod(InvPhi_1t_dot_IW_dot_InvTwo, InvPhi_1t_dot_IW, InvFour, p);
							if(!debug) MulMod(InvPhi_2t_dot_IW_dot_InvTwo, InvPhi_2t_dot_IW, InvFour, p);
							if(!debug) MulMod(InvPhi_3t_dot_IW_dot_InvTwo, InvPhi_3t_dot_IW, InvFour, p);
							//-------------------------------------------
                            if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],InvPhi_0t_dot_IW_dot_InvTwo,p);
                            if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],InvPhi_1t_dot_IW_dot_InvTwo,p);
					        if(!debug) MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],InvPhi_2t_dot_IW_dot_InvTwo,p);
					        if(!debug) MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],InvPhi_3t_dot_IW_dot_InvTwo,p);
                            break;
                    }
					INWC_DATARECORD <<"------after BU compute and Mul-------" << std::endl;
					INWC_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp] << endl;
					INWC_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp] << endl;
					INWC_DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<<A_B0R2[ma_tmp] << endl;
					INWC_DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<<A_B0R3[ma_tmp] << endl;
				
					
					INWC_DATARECORD <<"--------------------------------------------------------------------\n";
					if(j <  2)bn0_ma_reg1 = ma_tmp;
					if(j >= 2)bn0_ma_reg2 = ma_tmp;
				}
			    else {
					if(j < 2)bn1_bc_tmp = BC_tmp;
					PowerMod(factor_2t,factor_t,2,p);
					PowerMod(factor_3t,factor_t,3,p);
                    switch(s){
                        case 0:
							INWC_DATARECORD <<"Before butterfly unit operation! \n";
                            INWC_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<< ", st0_Tw[0] = " << st0_Tw[0] << endl;
					        INWC_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<< ", st0_Tw[1] = " << st0_Tw[1] << endl;
					        INWC_DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<<A_B1R2[ma_tmp]<< ", st0_Tw[2] = " << st0_Tw[2] << endl;
					        INWC_DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<<A_B1R3[ma_tmp]<< ", st0_Tw[3] = " << st0_Tw[3] << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_0t_Order = " << InvPhi_0t_Order << endl; else INWC_DATARECORD << "InvPhi_0t_Order = " << InvPhi_deg*0 << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_1t_Order = " << InvPhi_1t_Order << endl; else INWC_DATARECORD << "InvPhi_1t_Order = " << InvPhi_deg*1 << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_2t_Order = " << InvPhi_2t_Order << endl; else INWC_DATARECORD << "InvPhi_2t_Order = " << InvPhi_deg*2 << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_3t_Order = " << InvPhi_3t_Order << endl; else INWC_DATARECORD << "InvPhi_3t_Order = " << InvPhi_deg*3 << endl;
                            if(!debug) INWC_Radix4_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp],A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
					        INWC_DATARECORD <<"-------------------" << std::endl;
							//------------compute for INWC---------------
							if(!debug) MulMod(InvPhi_0t_dot_IW, InvPhi_0t_Order, st0_Tw[0], p);
							if(!debug) MulMod(InvPhi_1t_dot_IW, InvPhi_1t_Order, st0_Tw[1], p);
							if(!debug) MulMod(InvPhi_2t_dot_IW, InvPhi_2t_Order, st0_Tw[2], p);
							if(!debug) MulMod(InvPhi_3t_dot_IW, InvPhi_3t_Order, st0_Tw[3], p);

							if(!debug) MulMod(InvPhi_0t_dot_IW_dot_InvTwo, InvPhi_0t_dot_IW, InvFour, p);
							if(!debug) MulMod(InvPhi_1t_dot_IW_dot_InvTwo, InvPhi_1t_dot_IW, InvFour, p);
							if(!debug) MulMod(InvPhi_2t_dot_IW_dot_InvTwo, InvPhi_2t_dot_IW, InvFour, p);
							if(!debug) MulMod(InvPhi_3t_dot_IW_dot_InvTwo, InvPhi_3t_dot_IW, InvFour, p);
							//-------------------------------------------
                            if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],InvPhi_0t_dot_IW_dot_InvTwo,p);
					        if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],InvPhi_1t_dot_IW_dot_InvTwo,p);
					        if(!debug) MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],InvPhi_2t_dot_IW_dot_InvTwo,p);
					        if(!debug) MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],InvPhi_3t_dot_IW_dot_InvTwo,p);
                            break;
                        case 1:
							INWC_DATARECORD <<"Before butterfly unit operation! \n";
                            INWC_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<< ", st1_Tw[0] = " << st1_Tw[0] << endl;
					        INWC_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<< ", st1_Tw[1] = " << st1_Tw[1] << endl;
					        INWC_DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<<A_B1R2[ma_tmp]<< ", st1_Tw[2] = " << st1_Tw[2] << endl;
					        INWC_DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<<A_B1R3[ma_tmp]<< ", st1_Tw[3] = " << st1_Tw[3] << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_0t_Order = " << InvPhi_0t_Order << endl; else INWC_DATARECORD << "InvPhi_0t_Order = " << InvPhi_deg*0 << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_1t_Order = " << InvPhi_1t_Order << endl; else INWC_DATARECORD << "InvPhi_1t_Order = " << InvPhi_deg*1 << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_2t_Order = " << InvPhi_2t_Order << endl; else INWC_DATARECORD << "InvPhi_2t_Order = " << InvPhi_deg*2 << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_3t_Order = " << InvPhi_3t_Order << endl; else INWC_DATARECORD << "InvPhi_3t_Order = " << InvPhi_deg*3 << endl;
                            if(!debug) INWC_Radix4_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp],A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
					        INWC_DATARECORD <<"-------------------" << std::endl;
							//------------compute for INWC---------------
							if(!debug) MulMod(InvPhi_0t_dot_IW, InvPhi_0t_Order, st1_Tw[0], p);
							if(!debug) MulMod(InvPhi_1t_dot_IW, InvPhi_1t_Order, st1_Tw[1], p);
							if(!debug) MulMod(InvPhi_2t_dot_IW, InvPhi_2t_Order, st1_Tw[2], p);
							if(!debug) MulMod(InvPhi_3t_dot_IW, InvPhi_3t_Order, st1_Tw[3], p);

							if(!debug) MulMod(InvPhi_0t_dot_IW_dot_InvTwo, InvPhi_0t_dot_IW, InvFour, p);
							if(!debug) MulMod(InvPhi_1t_dot_IW_dot_InvTwo, InvPhi_1t_dot_IW, InvFour, p);
							if(!debug) MulMod(InvPhi_2t_dot_IW_dot_InvTwo, InvPhi_2t_dot_IW, InvFour, p);
							if(!debug) MulMod(InvPhi_3t_dot_IW_dot_InvTwo, InvPhi_3t_dot_IW, InvFour, p);
							//-------------------------------------------
                            if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],InvPhi_0t_dot_IW_dot_InvTwo,p);
					        if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],InvPhi_1t_dot_IW_dot_InvTwo,p);
					        if(!debug) MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],InvPhi_2t_dot_IW_dot_InvTwo,p);
					        if(!debug) MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],InvPhi_3t_dot_IW_dot_InvTwo,p);
                            break;
                        case 2:
							INWC_DATARECORD <<"Before butterfly unit operation! \n";
                            INWC_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<< ", st2_Tw[0] = " << st2_Tw[0] << endl;
					        INWC_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<< ", st2_Tw[1] = " << st2_Tw[1] << endl;
					        INWC_DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<<A_B1R2[ma_tmp]<< ", st2_Tw[2] = " << st2_Tw[2] << endl;
					        INWC_DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<<A_B1R3[ma_tmp]<< ", st2_Tw[3] = " << st2_Tw[3] << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_0t_Order = " << InvPhi_0t_Order << endl; else INWC_DATARECORD << "InvPhi_0t_Order = " << InvPhi_deg*0 << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_1t_Order = " << InvPhi_1t_Order << endl; else INWC_DATARECORD << "InvPhi_1t_Order = " << InvPhi_deg*1 << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_2t_Order = " << InvPhi_2t_Order << endl; else INWC_DATARECORD << "InvPhi_2t_Order = " << InvPhi_deg*2 << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_3t_Order = " << InvPhi_3t_Order << endl; else INWC_DATARECORD << "InvPhi_3t_Order = " << InvPhi_deg*3 << endl;
                            if(!debug) INWC_Radix4_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp],A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
					        INWC_DATARECORD <<"-------------------" << std::endl;
							//------------compute for INWC---------------
							if(!debug) MulMod(InvPhi_0t_dot_IW, InvPhi_0t_Order, st2_Tw[0], p);
							if(!debug) MulMod(InvPhi_1t_dot_IW, InvPhi_1t_Order, st2_Tw[1], p);
							if(!debug) MulMod(InvPhi_2t_dot_IW, InvPhi_2t_Order, st2_Tw[2], p);
							if(!debug) MulMod(InvPhi_3t_dot_IW, InvPhi_3t_Order, st2_Tw[3], p);

							if(!debug) MulMod(InvPhi_0t_dot_IW_dot_InvTwo, InvPhi_0t_dot_IW, InvFour, p);
							if(!debug) MulMod(InvPhi_1t_dot_IW_dot_InvTwo, InvPhi_1t_dot_IW, InvFour, p);
							if(!debug) MulMod(InvPhi_2t_dot_IW_dot_InvTwo, InvPhi_2t_dot_IW, InvFour, p);
							if(!debug) MulMod(InvPhi_3t_dot_IW_dot_InvTwo, InvPhi_3t_dot_IW, InvFour, p);
							//-------------------------------------------
                            if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],InvPhi_0t_dot_IW_dot_InvTwo,p);
					        if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],InvPhi_1t_dot_IW_dot_InvTwo,p);
					        if(!debug) MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],InvPhi_2t_dot_IW_dot_InvTwo,p);
					        if(!debug) MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],InvPhi_3t_dot_IW_dot_InvTwo,p);
                            break;
                        case 3:
							INWC_DATARECORD <<"Before butterfly unit operation! \n";
                            INWC_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<< ", st3_Tw[0] = " << 1 << endl;
					        INWC_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<< ", st3_Tw[1] = " << 1 << endl;
					        INWC_DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<<A_B1R2[ma_tmp]<< ", st3_Tw[2] = " << 1 << endl;
					        INWC_DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<<A_B1R3[ma_tmp]<< ", st3_Tw[3] = " << 1 << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_0t_Order = " << InvPhi_0t_Order << endl; else INWC_DATARECORD << "InvPhi_0t_Order = " << InvPhi_deg*0 << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_1t_Order = " << InvPhi_1t_Order << endl; else INWC_DATARECORD << "InvPhi_1t_Order = " << InvPhi_deg*1 << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_2t_Order = " << InvPhi_2t_Order << endl; else INWC_DATARECORD << "InvPhi_2t_Order = " << InvPhi_deg*2 << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_3t_Order = " << InvPhi_3t_Order << endl; else INWC_DATARECORD << "InvPhi_3t_Order = " << InvPhi_deg*3 << endl;
                            if(!debug) INWC_Radix4_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp],A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
					        INWC_DATARECORD <<"-------------------" << std::endl;
							//------------compute for INWC---------------
							if(!debug) MulMod(InvPhi_0t_dot_IW, InvPhi_0t_Order, 1, p);
							if(!debug) MulMod(InvPhi_1t_dot_IW, InvPhi_1t_Order, 1, p);
							if(!debug) MulMod(InvPhi_2t_dot_IW, InvPhi_2t_Order, 1, p);
							if(!debug) MulMod(InvPhi_3t_dot_IW, InvPhi_3t_Order, 1, p);

							if(!debug) MulMod(InvPhi_0t_dot_IW_dot_InvTwo, InvPhi_0t_dot_IW, InvFour, p);
							if(!debug) MulMod(InvPhi_1t_dot_IW_dot_InvTwo, InvPhi_1t_dot_IW, InvFour, p);
							if(!debug) MulMod(InvPhi_2t_dot_IW_dot_InvTwo, InvPhi_2t_dot_IW, InvFour, p);
							if(!debug) MulMod(InvPhi_3t_dot_IW_dot_InvTwo, InvPhi_3t_dot_IW, InvFour, p);
							//-------------------------------------------
                            if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],InvPhi_0t_dot_IW_dot_InvTwo,p);
					        if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],InvPhi_1t_dot_IW_dot_InvTwo,p);
					        if(!debug) MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],InvPhi_2t_dot_IW_dot_InvTwo,p);
					        if(!debug) MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],InvPhi_3t_dot_IW_dot_InvTwo,p);
                            break;
                    }
                    INWC_DATARECORD <<"------after BU compute and Mul-------" << std::endl;
					INWC_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: " << A_B1R0[ma_tmp] << endl;
					INWC_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: " << A_B1R1[ma_tmp] << endl;
					INWC_DATARECORD <<"A_B1R2["<<ma_tmp<<"]: " << A_B1R2[ma_tmp] << endl;
					INWC_DATARECORD <<"A_B1R3["<<ma_tmp<<"]: " << A_B1R3[ma_tmp] << endl;
					
					INWC_DATARECORD <<"--------------------------------------------------------------------\n";
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

void DIF_INWC::DIF_INWC_radix16(std::vector<ZZ> &A){
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
	
	std::ofstream INWC_DATARECORD("./NWC_PrintData/INWC_R16_SPMB.txt");

    //------------DTFAG generator-------------
	DTFAG DTFAG;
	vector<ZZ > st0_Tw, st1_Tw, st2_Tw;
	st0_Tw.resize(radix);
	st1_Tw.resize(radix);
	st2_Tw.resize(radix);
	int DTFAG_t = 0;
	int DTFAG_i = 0;
	int DTFAG_j = 0;

	int fft_point = N;
	int radix_r1 = radix;
	int radix_r2 = radix;
	ZZ fft_twiddle = IW; //***due to INWC, this ways use Inverse fft_twiddle***
	ZZ fft_prime = p;
	int debug = 0;
	vector<vector<ZZ > > ROM0;
    vector<vector<ZZ > > ROM1, ROM2;
    
    
	int arr_size = radix_r1 * radix_r1;
    ROM0.resize(radix_r1);
    for(int i=0; i<radix_r1; i++){
        ROM0[i].resize(radix_r1);
    }
    ROM1.resize(arr_size);
	for(int i=0; i<radix_r1; i++){
        ROM1[i].resize(radix_r1);
    }
    ROM2.resize(arr_size);
	for(int i=0; i<radix_r1; i++){
        ROM2[i].resize(radix_r1);
    }
	DTFAG.DTFAG_ROM_init(
        radix_r1, radix_r2, fft_twiddle, fft_prime, debug,
        ROM0, ROM1, ROM2);
	////----------------------------------------

	//-----------NWC PART-----------------------
	ZZ InvSixTeen;
	ZZ InvPhi_0t_dot_IW	, InvPhi_0t_dot_IW_dot_InvSixTeen	;
	ZZ InvPhi_1t_dot_IW	, InvPhi_1t_dot_IW_dot_InvSixTeen	;
	ZZ InvPhi_2t_dot_IW	, InvPhi_2t_dot_IW_dot_InvSixTeen	;
	ZZ InvPhi_3t_dot_IW	, InvPhi_3t_dot_IW_dot_InvSixTeen	;
	ZZ InvPhi_4t_dot_IW	, InvPhi_4t_dot_IW_dot_InvSixTeen	;
	ZZ InvPhi_5t_dot_IW	, InvPhi_5t_dot_IW_dot_InvSixTeen	;
	ZZ InvPhi_6t_dot_IW	, InvPhi_6t_dot_IW_dot_InvSixTeen	;
	ZZ InvPhi_7t_dot_IW	, InvPhi_7t_dot_IW_dot_InvSixTeen	;
	ZZ InvPhi_8t_dot_IW	, InvPhi_8t_dot_IW_dot_InvSixTeen	;
	ZZ InvPhi_9t_dot_IW	, InvPhi_9t_dot_IW_dot_InvSixTeen	;
	ZZ InvPhi_10t_dot_IW, InvPhi_10t_dot_IW_dot_InvSixTeen	;
	ZZ InvPhi_11t_dot_IW, InvPhi_11t_dot_IW_dot_InvSixTeen	;
	ZZ InvPhi_12t_dot_IW, InvPhi_12t_dot_IW_dot_InvSixTeen	;
	ZZ InvPhi_13t_dot_IW, InvPhi_13t_dot_IW_dot_InvSixTeen	;
	ZZ InvPhi_14t_dot_IW, InvPhi_14t_dot_IW_dot_InvSixTeen	;
	ZZ InvPhi_15t_dot_IW, InvPhi_15t_dot_IW_dot_InvSixTeen	;
	InvMod(InvSixTeen, (ZZ)16, p);
	cout << "W = " << W << ", IW = " << IW <<  ", Phi = " << Phi << ", InvPhi = " << InvPhi << ", Inv_16 = " << InvSixTeen << endl;
	cout << "p = " << p << endl;
	//------------------------------------------
	
    Stage = (unsigned long)ceil(log2(N));
	Stage = (unsigned long)Stage / 4;
	BC_WIDTH  = (int)ceil(log2(N/16));
	offset    =  (int)N /radix;
	word_size =  (int)N / (2 * radix);
	group     =  (int)N / (radix * radix);
	bit_array_tmp.resize(BC_WIDTH);
	INWC_DATARECORD <<"Stage: "<< Stage<<"\n";

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
		INWC_DATARECORD <<"---------------------------------\n";
		INWC_DATARECORD <<"Now Stage: "<< s <<"\n";
		INWC_DATARECORD <<"twiddle factor : "<< factor <<"\n";

		std::cout << "twiddle factor : "<< factor <<"\n";
	    INWC_DATARECORD << "********\n";
		for(int i = 0 ;i < group;i++){
			bn0_bc_tmp  = 0;
			bn1_bc_tmp  = 0;	
			for(int j = 0;j < radix;j++){
				gray_i  = Gray(i,group);
				BC_tmp  = j * group + gray_i;
				INWC_DATARECORD << "i: " << i <<"\n";	
				INWC_DATARECORD << "gray_i: " << gray_i <<"\n";
				INWC_DATARECORD << "BC_tmp: " << BC_tmp <<"\n";
				RR_R16(BC_tmp,s,BC);
				INWC_DATARECORD << "BC: " << BC <<"\n";
				length = BC_tmp >> (4*s);
				INWC_DATARECORD << "length: " <<  length <<"\n";
				
				PowerMod(factor_t,factor,length,p);
				INWC_DATARECORD << "factor_t: "<<factor_t<<"\n";
				AGU_R16(BC,bn_tmp,ma_tmp);
				INWC_DATARECORD << "bn_tmp: "<<bn_tmp<<"\n";

                //-----------DTFAG generator-------------
                DTFAG.DTFAG_SPMB_DIF_MR(
                    s, fft_point, radix_r1, radix_r2, debug,
					ROM0, ROM1, ROM2,
                    st0_Tw, st1_Tw, st2_Tw,
                    DTFAG_i, DTFAG_t, DTFAG_j);
                /*switch(s){
					case 0:
						for(int i=0; i<radix; i++){
							cout << "st0_Tw[" << i << "] = w^" << st0_Tw[i] << endl;
						}
						break;
					case 1:
						for(int i=0; i<radix; i++){
							cout << "st1_Tw[" << i << "] = w^" << st1_Tw[i] << endl;
						}
						break;
					case 2:
						for(int i=0; i<radix; i++){
							cout << "st2_Tw[" << i << "] = w^" << st2_Tw[i] << endl;
						}
						break;
				}*/
                INWC_DATARECORD << "stage = " << s << ", DTFAG_i = " << DTFAG_i << ", DTFAG_t = " << DTFAG_t << ", DTFAG_j = " << DTFAG_j << endl;
                if(DTFAG_i == radix_r2-1 && DTFAG_j == radix_r1-1 && DTFAG_t == radix_r1-1){
					DTFAG_i = 0;
				}else if(DTFAG_t == radix_r1-1 && DTFAG_j == radix_r1-1){
					DTFAG_i++;
				}
				if(DTFAG_j == radix_r1-1 && DTFAG_t == radix_r1-1){
					DTFAG_t = 0;
				}else if(DTFAG_j == radix_r1-1){
					DTFAG_t++;
				}
				if(DTFAG_j == radix_r1-1){
					DTFAG_j = 0;
				}else{
					DTFAG_j++;
				}
                //---------------------------------------
				//---------NWC PART-------------
				ZZ 	InvPhi_0t, InvPhi_1t, InvPhi_2t, InvPhi_3t,
					InvPhi_4t, InvPhi_9t, InvPhi_6t, InvPhi_7t,
					InvPhi_8t, InvPhi_5t, InvPhi_10t, InvPhi_11t,
					InvPhi_12t, InvPhi_13t, InvPhi_14t, InvPhi_15t;
				ZZ 	InvPhi_0t_Order, InvPhi_1t_Order, InvPhi_2t_Order, InvPhi_3t_Order,
					InvPhi_4t_Order, InvPhi_5t_Order, InvPhi_6t_Order, InvPhi_7t_Order,
					InvPhi_8t_Order, InvPhi_9t_Order, InvPhi_10t_Order, InvPhi_11t_Order,
					InvPhi_12t_Order, InvPhi_13t_Order, InvPhi_14t_Order, InvPhi_15t_Order;
				ZZ InvPhi_deg = PowerMod((ZZ)16, s, p);
				InvPhi_0t  = PowerMod(InvPhi, 0, p);
				InvPhi_1t  = PowerMod(InvPhi, 1, p);
				InvPhi_2t  = PowerMod(InvPhi, 2, p);
				InvPhi_3t  = PowerMod(InvPhi, 3, p);
				InvPhi_4t  = PowerMod(InvPhi, 4, p);
				InvPhi_5t  = PowerMod(InvPhi, 5, p);
				InvPhi_6t  = PowerMod(InvPhi, 6, p);
				InvPhi_7t  = PowerMod(InvPhi, 7, p);
				InvPhi_8t  = PowerMod(InvPhi, 8, p);
				InvPhi_9t  = PowerMod(InvPhi, 9, p);
				InvPhi_10t = PowerMod(InvPhi, 10, p);
				InvPhi_11t = PowerMod(InvPhi, 11, p);
				InvPhi_12t = PowerMod(InvPhi, 12, p);
				InvPhi_13t = PowerMod(InvPhi, 13, p);
				InvPhi_14t = PowerMod(InvPhi, 14, p);
				InvPhi_15t = PowerMod(InvPhi, 15, p);
				InvPhi_0t_Order  = PowerMod(InvPhi_0t, InvPhi_deg, p);
				InvPhi_1t_Order  = PowerMod(InvPhi_1t, InvPhi_deg, p);
				InvPhi_2t_Order  = PowerMod(InvPhi_2t, InvPhi_deg, p);
				InvPhi_3t_Order  = PowerMod(InvPhi_3t, InvPhi_deg, p);
				InvPhi_4t_Order  = PowerMod(InvPhi_4t, InvPhi_deg, p);
				InvPhi_5t_Order  = PowerMod(InvPhi_5t, InvPhi_deg, p);
				InvPhi_6t_Order  = PowerMod(InvPhi_6t, InvPhi_deg, p);
				InvPhi_7t_Order  = PowerMod(InvPhi_7t, InvPhi_deg, p);
				InvPhi_8t_Order  = PowerMod(InvPhi_8t, InvPhi_deg, p);
				InvPhi_9t_Order  = PowerMod(InvPhi_9t, InvPhi_deg, p);
				InvPhi_10t_Order = PowerMod(InvPhi_10t, InvPhi_deg, p);
				InvPhi_11t_Order = PowerMod(InvPhi_11t, InvPhi_deg, p);
				InvPhi_12t_Order = PowerMod(InvPhi_12t, InvPhi_deg, p);
				InvPhi_13t_Order = PowerMod(InvPhi_13t, InvPhi_deg, p);
				InvPhi_14t_Order = PowerMod(InvPhi_14t, InvPhi_deg, p);
				InvPhi_15t_Order = PowerMod(InvPhi_15t, InvPhi_deg, p);
				//------------------------------
				if(bn_tmp == 0){
					if(j < 2)bn0_bc_tmp = BC_tmp;
					switch(s){
                        case 0:
                            INWC_DATARECORD <<" A_B0R0["<<ma_tmp<<"]:  "<<A_B0R0[ma_tmp]  << ", st0_Tw[0] = " << st0_Tw[0] << endl;
				            INWC_DATARECORD <<" A_B0R1["<<ma_tmp<<"]:  "<<A_B0R1[ma_tmp]  << ", st0_Tw[1] = " << st0_Tw[1] << endl;
					        INWC_DATARECORD <<" A_B0R2["<<ma_tmp<<"]:  "<<A_B0R2[ma_tmp]  << ", st0_Tw[2] = " << st0_Tw[2] << endl;
					        INWC_DATARECORD <<" A_B0R3["<<ma_tmp<<"]:  "<<A_B0R3[ma_tmp]  << ", st0_Tw[3] = " << st0_Tw[3] << endl;
					        INWC_DATARECORD <<" A_B0R4["<<ma_tmp<<"]:  "<<A_B0R4[ma_tmp]  << ", st0_Tw[4] = " << st0_Tw[4] << endl;
					        INWC_DATARECORD <<" A_B0R5["<<ma_tmp<<"]:  "<<A_B0R5[ma_tmp]  << ", st0_Tw[5] = " << st0_Tw[5] << endl;
					        INWC_DATARECORD <<" A_B0R6["<<ma_tmp<<"]:  "<<A_B0R6[ma_tmp]  << ", st0_Tw[6] = " << st0_Tw[6] << endl;
					        INWC_DATARECORD <<" A_B0R7["<<ma_tmp<<"]:  "<<A_B0R7[ma_tmp]  << ", st0_Tw[7] = " << st0_Tw[7] << endl;
					        INWC_DATARECORD <<" A_B0R8["<<ma_tmp<<"]:  "<<A_B0R8[ma_tmp]  << ", st0_Tw[8] = " << st0_Tw[8] << endl;
					        INWC_DATARECORD <<" A_B0R9["<<ma_tmp<<"]:  "<<A_B0R9[ma_tmp]  << ", st0_Tw[9] = " << st0_Tw[9] << endl;
					        INWC_DATARECORD <<" A_B0R10["<<ma_tmp<<"]: "<<A_B0R10[ma_tmp] << ", st0_Tw[10] = " << st0_Tw[10] << endl;   
					        INWC_DATARECORD <<" A_B0R11["<<ma_tmp<<"]: "<<A_B0R11[ma_tmp] << ", st0_Tw[11] = " << st0_Tw[11] << endl;   
					        INWC_DATARECORD <<" A_B0R12["<<ma_tmp<<"]: "<<A_B0R12[ma_tmp] << ", st0_Tw[12] = " << st0_Tw[12] << endl;   
					        INWC_DATARECORD <<" A_B0R13["<<ma_tmp<<"]: "<<A_B0R13[ma_tmp] << ", st0_Tw[13] = " << st0_Tw[13] << endl;   
					        INWC_DATARECORD <<" A_B0R14["<<ma_tmp<<"]: "<<A_B0R14[ma_tmp] << ", st0_Tw[14] = " << st0_Tw[14] << endl;   
					        INWC_DATARECORD <<" A_B0R15["<<ma_tmp<<"]: "<<A_B0R15[ma_tmp] << ", st0_Tw[15] = " << st0_Tw[15] << endl;   
                            if(!debug) INWC_DATARECORD << "InvPhi_0t_Order = " 	<< InvPhi_0t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_0t_Order = " 	<< InvPhi_deg*0 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_1t_Order = " 	<< InvPhi_1t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_1t_Order = " 	<< InvPhi_deg*1 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_2t_Order = " 	<< InvPhi_2t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_2t_Order = " 	<< InvPhi_deg*2 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_3t_Order = " 	<< InvPhi_3t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_3t_Order = " 	<< InvPhi_deg*3 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_4t_Order = " 	<< InvPhi_4t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_4t_Order = " 	<< InvPhi_deg*4 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_5t_Order = " 	<< InvPhi_5t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_5t_Order = " 	<< InvPhi_deg*5 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_6t_Order = " 	<< InvPhi_6t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_6t_Order = " 	<< InvPhi_deg*6 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_7t_Order = " 	<< InvPhi_7t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_7t_Order = " 	<< InvPhi_deg*7 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_8t_Order = " 	<< InvPhi_8t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_8t_Order = " 	<< InvPhi_deg*8 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_9t_Order = " 	<< InvPhi_9t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_9t_Order = " 	<< InvPhi_deg*9 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_10t_Order = " << InvPhi_10t_Order << endl; else INWC_DATARECORD << "InvPhi_10t_Order = " 	<< InvPhi_deg*10 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_11t_Order = " << InvPhi_11t_Order << endl; else INWC_DATARECORD << "InvPhi_11t_Order = " 	<< InvPhi_deg*11 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_12t_Order = " << InvPhi_12t_Order << endl; else INWC_DATARECORD << "InvPhi_12t_Order = " 	<< InvPhi_deg*12 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_13t_Order = " << InvPhi_13t_Order << endl; else INWC_DATARECORD << "InvPhi_13t_Order = " 	<< InvPhi_deg*13 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_14t_Order = " << InvPhi_14t_Order << endl; else INWC_DATARECORD << "InvPhi_14t_Order = " 	<< InvPhi_deg*14 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_15t_Order = " << InvPhi_15t_Order << endl; else INWC_DATARECORD << "InvPhi_15t_Order = " 	<< InvPhi_deg*15 	<< endl;
							INWC_Radix16_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp], A_B0R2[ma_tmp], A_B0R3[ma_tmp],A_B0R4[ma_tmp],
							   A_B0R5[ma_tmp],A_B0R6[ma_tmp], A_B0R7[ma_tmp], A_B0R8[ma_tmp],A_B0R9[ma_tmp],
							   A_B0R10[ma_tmp],A_B0R11[ma_tmp],A_B0R12[ma_tmp],A_B0R13[ma_tmp],A_B0R14[ma_tmp],
							   A_B0R15[ma_tmp]);
							//----------------------compute for INWC--------------------------
                            if(!debug) MulMod(InvPhi_0t_dot_IW, InvPhi_0t_Order, st0_Tw[0], p);
							if(!debug) MulMod(InvPhi_1t_dot_IW, InvPhi_1t_Order, st0_Tw[1], p);
							if(!debug) MulMod(InvPhi_2t_dot_IW, InvPhi_2t_Order, st0_Tw[2], p);
							if(!debug) MulMod(InvPhi_3t_dot_IW, InvPhi_3t_Order, st0_Tw[3], p);
							if(!debug) MulMod(InvPhi_4t_dot_IW, InvPhi_4t_Order, st0_Tw[4], p);
							if(!debug) MulMod(InvPhi_5t_dot_IW, InvPhi_5t_Order, st0_Tw[5], p);
							if(!debug) MulMod(InvPhi_6t_dot_IW, InvPhi_6t_Order, st0_Tw[6], p);
							if(!debug) MulMod(InvPhi_7t_dot_IW, InvPhi_7t_Order, st0_Tw[7], p);
							if(!debug) MulMod(InvPhi_8t_dot_IW, InvPhi_8t_Order, st0_Tw[8], p);
							if(!debug) MulMod(InvPhi_9t_dot_IW, InvPhi_9t_Order, st0_Tw[9], p);
							if(!debug) MulMod(InvPhi_10t_dot_IW, InvPhi_10t_Order, st0_Tw[10], p);
							if(!debug) MulMod(InvPhi_11t_dot_IW, InvPhi_11t_Order, st0_Tw[11], p);
							if(!debug) MulMod(InvPhi_12t_dot_IW, InvPhi_12t_Order, st0_Tw[12], p);
							if(!debug) MulMod(InvPhi_13t_dot_IW, InvPhi_13t_Order, st0_Tw[13], p);
							if(!debug) MulMod(InvPhi_14t_dot_IW, InvPhi_14t_Order, st0_Tw[14], p);
							if(!debug) MulMod(InvPhi_15t_dot_IW, InvPhi_15t_Order, st0_Tw[15], p);

							if(!debug) MulMod(InvPhi_0t_dot_IW_dot_InvSixTeen, InvPhi_0t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_1t_dot_IW_dot_InvSixTeen, InvPhi_1t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_2t_dot_IW_dot_InvSixTeen, InvPhi_2t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_3t_dot_IW_dot_InvSixTeen, InvPhi_3t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_4t_dot_IW_dot_InvSixTeen, InvPhi_4t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_5t_dot_IW_dot_InvSixTeen, InvPhi_5t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_6t_dot_IW_dot_InvSixTeen, InvPhi_6t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_7t_dot_IW_dot_InvSixTeen, InvPhi_7t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_8t_dot_IW_dot_InvSixTeen, InvPhi_8t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_9t_dot_IW_dot_InvSixTeen, InvPhi_9t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_10t_dot_IW_dot_InvSixTeen, InvPhi_10t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_11t_dot_IW_dot_InvSixTeen, InvPhi_11t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_12t_dot_IW_dot_InvSixTeen, InvPhi_12t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_13t_dot_IW_dot_InvSixTeen, InvPhi_13t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_14t_dot_IW_dot_InvSixTeen, InvPhi_14t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_15t_dot_IW_dot_InvSixTeen, InvPhi_15t_dot_IW, InvSixTeen, p);
							//----------------------------------------------------------------
							MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],InvPhi_0t_dot_IW_dot_InvSixTeen,p);
                            MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],InvPhi_1t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],InvPhi_2t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],InvPhi_3t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R4[ma_tmp],A_B0R4[ma_tmp],InvPhi_4t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R5[ma_tmp],A_B0R5[ma_tmp],InvPhi_5t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R6[ma_tmp],A_B0R6[ma_tmp],InvPhi_6t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R7[ma_tmp],A_B0R7[ma_tmp],InvPhi_7t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R8[ma_tmp],A_B0R8[ma_tmp],InvPhi_8t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R9[ma_tmp],A_B0R9[ma_tmp],InvPhi_9t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R10[ma_tmp],A_B0R10[ma_tmp],InvPhi_10t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R11[ma_tmp],A_B0R11[ma_tmp],InvPhi_11t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R12[ma_tmp],A_B0R12[ma_tmp],InvPhi_12t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R13[ma_tmp],A_B0R13[ma_tmp],InvPhi_13t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R14[ma_tmp],A_B0R14[ma_tmp],InvPhi_14t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R15[ma_tmp],A_B0R15[ma_tmp],InvPhi_15t_dot_IW_dot_InvSixTeen,p);
                            break;
                        case 1:
                            INWC_DATARECORD <<" A_B0R0["<<ma_tmp<<"]:  "<<A_B0R0[ma_tmp]  << ", st1_Tw[0] = " << st1_Tw[0] << endl;
				            INWC_DATARECORD <<" A_B0R1["<<ma_tmp<<"]:  "<<A_B0R1[ma_tmp]  << ", st1_Tw[1] = " << st1_Tw[1] << endl;
					        INWC_DATARECORD <<" A_B0R2["<<ma_tmp<<"]:  "<<A_B0R2[ma_tmp]  << ", st1_Tw[2] = " << st1_Tw[2] << endl;
					        INWC_DATARECORD <<" A_B0R3["<<ma_tmp<<"]:  "<<A_B0R3[ma_tmp]  << ", st1_Tw[3] = " << st1_Tw[3] << endl;
					        INWC_DATARECORD <<" A_B0R4["<<ma_tmp<<"]:  "<<A_B0R4[ma_tmp]  << ", st1_Tw[4] = " << st1_Tw[4] << endl;
					        INWC_DATARECORD <<" A_B0R5["<<ma_tmp<<"]:  "<<A_B0R5[ma_tmp]  << ", st1_Tw[5] = " << st1_Tw[5] << endl;
					        INWC_DATARECORD <<" A_B0R6["<<ma_tmp<<"]:  "<<A_B0R6[ma_tmp]  << ", st1_Tw[6] = " << st1_Tw[6] << endl;
					        INWC_DATARECORD <<" A_B0R7["<<ma_tmp<<"]:  "<<A_B0R7[ma_tmp]  << ", st1_Tw[7] = " << st1_Tw[7] << endl;
					        INWC_DATARECORD <<" A_B0R8["<<ma_tmp<<"]:  "<<A_B0R8[ma_tmp]  << ", st1_Tw[8] = " << st1_Tw[8] << endl;
					        INWC_DATARECORD <<" A_B0R9["<<ma_tmp<<"]:  "<<A_B0R9[ma_tmp]  << ", st1_Tw[9] = " << st1_Tw[9] << endl;
					        INWC_DATARECORD <<" A_B0R10["<<ma_tmp<<"]: "<<A_B0R10[ma_tmp] << ", st1_Tw[10] = " << st1_Tw[10] << endl;   
					        INWC_DATARECORD <<" A_B0R11["<<ma_tmp<<"]: "<<A_B0R11[ma_tmp] << ", st1_Tw[11] = " << st1_Tw[11] << endl;   
					        INWC_DATARECORD <<" A_B0R12["<<ma_tmp<<"]: "<<A_B0R12[ma_tmp] << ", st1_Tw[12] = " << st1_Tw[12] << endl;   
					        INWC_DATARECORD <<" A_B0R13["<<ma_tmp<<"]: "<<A_B0R13[ma_tmp] << ", st1_Tw[13] = " << st1_Tw[13] << endl;   
					        INWC_DATARECORD <<" A_B0R14["<<ma_tmp<<"]: "<<A_B0R14[ma_tmp] << ", st1_Tw[14] = " << st1_Tw[14] << endl;   
					        INWC_DATARECORD <<" A_B0R15["<<ma_tmp<<"]: "<<A_B0R15[ma_tmp] << ", st1_Tw[15] = " << st1_Tw[15] << endl;   
                            if(!debug) INWC_DATARECORD << "InvPhi_0t_Order = " 	<< InvPhi_0t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_0t_Order = " 	<< InvPhi_deg*0 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_1t_Order = " 	<< InvPhi_1t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_1t_Order = " 	<< InvPhi_deg*1 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_2t_Order = " 	<< InvPhi_2t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_2t_Order = " 	<< InvPhi_deg*2 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_3t_Order = " 	<< InvPhi_3t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_3t_Order = " 	<< InvPhi_deg*3 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_4t_Order = " 	<< InvPhi_4t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_4t_Order = " 	<< InvPhi_deg*4 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_5t_Order = " 	<< InvPhi_5t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_5t_Order = " 	<< InvPhi_deg*5 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_6t_Order = " 	<< InvPhi_6t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_6t_Order = " 	<< InvPhi_deg*6 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_7t_Order = " 	<< InvPhi_7t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_7t_Order = " 	<< InvPhi_deg*7 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_8t_Order = " 	<< InvPhi_8t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_8t_Order = " 	<< InvPhi_deg*8 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_9t_Order = " 	<< InvPhi_9t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_9t_Order = " 	<< InvPhi_deg*9 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_10t_Order = " << InvPhi_10t_Order << endl; else INWC_DATARECORD << "InvPhi_10t_Order = " 	<< InvPhi_deg*10 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_11t_Order = " << InvPhi_11t_Order << endl; else INWC_DATARECORD << "InvPhi_11t_Order = " 	<< InvPhi_deg*11 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_12t_Order = " << InvPhi_12t_Order << endl; else INWC_DATARECORD << "InvPhi_12t_Order = " 	<< InvPhi_deg*12 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_13t_Order = " << InvPhi_13t_Order << endl; else INWC_DATARECORD << "InvPhi_13t_Order = " 	<< InvPhi_deg*13 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_14t_Order = " << InvPhi_14t_Order << endl; else INWC_DATARECORD << "InvPhi_14t_Order = " 	<< InvPhi_deg*14 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_15t_Order = " << InvPhi_15t_Order << endl; else INWC_DATARECORD << "InvPhi_15t_Order = " 	<< InvPhi_deg*15 	<< endl;
							INWC_Radix16_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp], A_B0R2[ma_tmp], A_B0R3[ma_tmp],A_B0R4[ma_tmp],
							   A_B0R5[ma_tmp],A_B0R6[ma_tmp], A_B0R7[ma_tmp], A_B0R8[ma_tmp],A_B0R9[ma_tmp],
							   A_B0R10[ma_tmp],A_B0R11[ma_tmp],A_B0R12[ma_tmp],A_B0R13[ma_tmp],A_B0R14[ma_tmp],
							   A_B0R15[ma_tmp]);
							//----------------------compute for INWC--------------------------
                            if(!debug) MulMod(InvPhi_0t_dot_IW, InvPhi_0t_Order, st1_Tw[0], p);
							if(!debug) MulMod(InvPhi_1t_dot_IW, InvPhi_1t_Order, st1_Tw[1], p);
							if(!debug) MulMod(InvPhi_2t_dot_IW, InvPhi_2t_Order, st1_Tw[2], p);
							if(!debug) MulMod(InvPhi_3t_dot_IW, InvPhi_3t_Order, st1_Tw[3], p);
							if(!debug) MulMod(InvPhi_4t_dot_IW, InvPhi_4t_Order, st1_Tw[4], p);
							if(!debug) MulMod(InvPhi_5t_dot_IW, InvPhi_5t_Order, st1_Tw[5], p);
							if(!debug) MulMod(InvPhi_6t_dot_IW, InvPhi_6t_Order, st1_Tw[6], p);
							if(!debug) MulMod(InvPhi_7t_dot_IW, InvPhi_7t_Order, st1_Tw[7], p);
							if(!debug) MulMod(InvPhi_8t_dot_IW, InvPhi_8t_Order, st1_Tw[8], p);
							if(!debug) MulMod(InvPhi_9t_dot_IW, InvPhi_9t_Order, st1_Tw[9], p);
							if(!debug) MulMod(InvPhi_10t_dot_IW, InvPhi_10t_Order, st1_Tw[10], p);
							if(!debug) MulMod(InvPhi_11t_dot_IW, InvPhi_11t_Order, st1_Tw[11], p);
							if(!debug) MulMod(InvPhi_12t_dot_IW, InvPhi_12t_Order, st1_Tw[12], p);
							if(!debug) MulMod(InvPhi_13t_dot_IW, InvPhi_13t_Order, st1_Tw[13], p);
							if(!debug) MulMod(InvPhi_14t_dot_IW, InvPhi_14t_Order, st1_Tw[14], p);
							if(!debug) MulMod(InvPhi_15t_dot_IW, InvPhi_15t_Order, st1_Tw[15], p);

							if(!debug) MulMod(InvPhi_0t_dot_IW_dot_InvSixTeen, InvPhi_0t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_1t_dot_IW_dot_InvSixTeen, InvPhi_1t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_2t_dot_IW_dot_InvSixTeen, InvPhi_2t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_3t_dot_IW_dot_InvSixTeen, InvPhi_3t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_4t_dot_IW_dot_InvSixTeen, InvPhi_4t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_5t_dot_IW_dot_InvSixTeen, InvPhi_5t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_6t_dot_IW_dot_InvSixTeen, InvPhi_6t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_7t_dot_IW_dot_InvSixTeen, InvPhi_7t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_8t_dot_IW_dot_InvSixTeen, InvPhi_8t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_9t_dot_IW_dot_InvSixTeen, InvPhi_9t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_10t_dot_IW_dot_InvSixTeen, InvPhi_10t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_11t_dot_IW_dot_InvSixTeen, InvPhi_11t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_12t_dot_IW_dot_InvSixTeen, InvPhi_12t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_13t_dot_IW_dot_InvSixTeen, InvPhi_13t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_14t_dot_IW_dot_InvSixTeen, InvPhi_14t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_15t_dot_IW_dot_InvSixTeen, InvPhi_15t_dot_IW, InvSixTeen, p);
							//----------------------------------------------------------------
                            MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],InvPhi_0t_dot_IW_dot_InvSixTeen,p);
                            MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],InvPhi_1t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],InvPhi_2t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],InvPhi_3t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R4[ma_tmp],A_B0R4[ma_tmp],InvPhi_4t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R5[ma_tmp],A_B0R5[ma_tmp],InvPhi_5t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R6[ma_tmp],A_B0R6[ma_tmp],InvPhi_6t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R7[ma_tmp],A_B0R7[ma_tmp],InvPhi_7t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R8[ma_tmp],A_B0R8[ma_tmp],InvPhi_8t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R9[ma_tmp],A_B0R9[ma_tmp],InvPhi_9t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R10[ma_tmp],A_B0R10[ma_tmp],InvPhi_10t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R11[ma_tmp],A_B0R11[ma_tmp],InvPhi_11t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R12[ma_tmp],A_B0R12[ma_tmp],InvPhi_12t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R13[ma_tmp],A_B0R13[ma_tmp],InvPhi_13t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R14[ma_tmp],A_B0R14[ma_tmp],InvPhi_14t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R15[ma_tmp],A_B0R15[ma_tmp],InvPhi_15t_dot_IW_dot_InvSixTeen,p);
                            break;
                        case 2:
                            INWC_DATARECORD <<" A_B0R0["<<ma_tmp<<"]:  "<<A_B0R0[ma_tmp]  << ", st2_Tw[0] = " << st2_Tw[0] << endl;
				            INWC_DATARECORD <<" A_B0R1["<<ma_tmp<<"]:  "<<A_B0R1[ma_tmp]  << ", st2_Tw[1] = " << st2_Tw[1] << endl;
					        INWC_DATARECORD <<" A_B0R2["<<ma_tmp<<"]:  "<<A_B0R2[ma_tmp]  << ", st2_Tw[2] = " << st2_Tw[2] << endl;
					        INWC_DATARECORD <<" A_B0R3["<<ma_tmp<<"]:  "<<A_B0R3[ma_tmp]  << ", st2_Tw[3] = " << st2_Tw[3] << endl;
					        INWC_DATARECORD <<" A_B0R4["<<ma_tmp<<"]:  "<<A_B0R4[ma_tmp]  << ", st2_Tw[4] = " << st2_Tw[4] << endl;
					        INWC_DATARECORD <<" A_B0R5["<<ma_tmp<<"]:  "<<A_B0R5[ma_tmp]  << ", st2_Tw[5] = " << st2_Tw[5] << endl;
					        INWC_DATARECORD <<" A_B0R6["<<ma_tmp<<"]:  "<<A_B0R6[ma_tmp]  << ", st2_Tw[6] = " << st2_Tw[6] << endl;
					        INWC_DATARECORD <<" A_B0R7["<<ma_tmp<<"]:  "<<A_B0R7[ma_tmp]  << ", st2_Tw[7] = " << st2_Tw[7] << endl;
					        INWC_DATARECORD <<" A_B0R8["<<ma_tmp<<"]:  "<<A_B0R8[ma_tmp]  << ", st2_Tw[8] = " << st2_Tw[8] << endl;
					        INWC_DATARECORD <<" A_B0R9["<<ma_tmp<<"]:  "<<A_B0R9[ma_tmp]  << ", st2_Tw[9] = " << st2_Tw[9] << endl;
					        INWC_DATARECORD <<" A_B0R10["<<ma_tmp<<"]: "<<A_B0R10[ma_tmp] << ", st2_Tw[10] = " << st2_Tw[10] << endl;   
					        INWC_DATARECORD <<" A_B0R11["<<ma_tmp<<"]: "<<A_B0R11[ma_tmp] << ", st2_Tw[11] = " << st2_Tw[11] << endl;   
					        INWC_DATARECORD <<" A_B0R12["<<ma_tmp<<"]: "<<A_B0R12[ma_tmp] << ", st2_Tw[12] = " << st2_Tw[12] << endl;   
					        INWC_DATARECORD <<" A_B0R13["<<ma_tmp<<"]: "<<A_B0R13[ma_tmp] << ", st2_Tw[13] = " << st2_Tw[13] << endl;   
					        INWC_DATARECORD <<" A_B0R14["<<ma_tmp<<"]: "<<A_B0R14[ma_tmp] << ", st2_Tw[14] = " << st2_Tw[14] << endl;   
					        INWC_DATARECORD <<" A_B0R15["<<ma_tmp<<"]: "<<A_B0R15[ma_tmp] << ", st2_Tw[15] = " << st2_Tw[15] << endl;   
                            if(!debug) INWC_DATARECORD << "InvPhi_0t_Order = " 	<< InvPhi_0t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_0t_Order = " 	<< InvPhi_deg*0 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_1t_Order = " 	<< InvPhi_1t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_1t_Order = " 	<< InvPhi_deg*1 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_2t_Order = " 	<< InvPhi_2t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_2t_Order = " 	<< InvPhi_deg*2 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_3t_Order = " 	<< InvPhi_3t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_3t_Order = " 	<< InvPhi_deg*3 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_4t_Order = " 	<< InvPhi_4t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_4t_Order = " 	<< InvPhi_deg*4 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_5t_Order = " 	<< InvPhi_5t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_5t_Order = " 	<< InvPhi_deg*5 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_6t_Order = " 	<< InvPhi_6t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_6t_Order = " 	<< InvPhi_deg*6 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_7t_Order = " 	<< InvPhi_7t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_7t_Order = " 	<< InvPhi_deg*7 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_8t_Order = " 	<< InvPhi_8t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_8t_Order = " 	<< InvPhi_deg*8 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_9t_Order = " 	<< InvPhi_9t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_9t_Order = " 	<< InvPhi_deg*9 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_10t_Order = " << InvPhi_10t_Order << endl; else INWC_DATARECORD << "InvPhi_10t_Order = " 	<< InvPhi_deg*10 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_11t_Order = " << InvPhi_11t_Order << endl; else INWC_DATARECORD << "InvPhi_11t_Order = " 	<< InvPhi_deg*11 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_12t_Order = " << InvPhi_12t_Order << endl; else INWC_DATARECORD << "InvPhi_12t_Order = " 	<< InvPhi_deg*12 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_13t_Order = " << InvPhi_13t_Order << endl; else INWC_DATARECORD << "InvPhi_13t_Order = " 	<< InvPhi_deg*13 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_14t_Order = " << InvPhi_14t_Order << endl; else INWC_DATARECORD << "InvPhi_14t_Order = " 	<< InvPhi_deg*14 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_15t_Order = " << InvPhi_15t_Order << endl; else INWC_DATARECORD << "InvPhi_15t_Order = " 	<< InvPhi_deg*15 	<< endl;
							INWC_Radix16_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp], A_B0R2[ma_tmp], A_B0R3[ma_tmp],A_B0R4[ma_tmp],
							   A_B0R5[ma_tmp],A_B0R6[ma_tmp], A_B0R7[ma_tmp], A_B0R8[ma_tmp],A_B0R9[ma_tmp],
							   A_B0R10[ma_tmp],A_B0R11[ma_tmp],A_B0R12[ma_tmp],A_B0R13[ma_tmp],A_B0R14[ma_tmp],
							   A_B0R15[ma_tmp]);
							//----------------------compute for INWC--------------------------
                            if(!debug) MulMod(InvPhi_0t_dot_IW, InvPhi_0t_Order, st2_Tw[0], p);
							if(!debug) MulMod(InvPhi_1t_dot_IW, InvPhi_1t_Order, st2_Tw[1], p);
							if(!debug) MulMod(InvPhi_2t_dot_IW, InvPhi_2t_Order, st2_Tw[2], p);
							if(!debug) MulMod(InvPhi_3t_dot_IW, InvPhi_3t_Order, st2_Tw[3], p);
							if(!debug) MulMod(InvPhi_4t_dot_IW, InvPhi_4t_Order, st2_Tw[4], p);
							if(!debug) MulMod(InvPhi_5t_dot_IW, InvPhi_5t_Order, st2_Tw[5], p);
							if(!debug) MulMod(InvPhi_6t_dot_IW, InvPhi_6t_Order, st2_Tw[6], p);
							if(!debug) MulMod(InvPhi_7t_dot_IW, InvPhi_7t_Order, st2_Tw[7], p);
							if(!debug) MulMod(InvPhi_8t_dot_IW, InvPhi_8t_Order, st2_Tw[8], p);
							if(!debug) MulMod(InvPhi_9t_dot_IW, InvPhi_9t_Order, st2_Tw[9], p);
							if(!debug) MulMod(InvPhi_10t_dot_IW, InvPhi_10t_Order, st2_Tw[10], p);
							if(!debug) MulMod(InvPhi_11t_dot_IW, InvPhi_11t_Order, st2_Tw[11], p);
							if(!debug) MulMod(InvPhi_12t_dot_IW, InvPhi_12t_Order, st2_Tw[12], p);
							if(!debug) MulMod(InvPhi_13t_dot_IW, InvPhi_13t_Order, st2_Tw[13], p);
							if(!debug) MulMod(InvPhi_14t_dot_IW, InvPhi_14t_Order, st2_Tw[14], p);
							if(!debug) MulMod(InvPhi_15t_dot_IW, InvPhi_15t_Order, st2_Tw[15], p);

							if(!debug) MulMod(InvPhi_0t_dot_IW_dot_InvSixTeen, InvPhi_0t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_1t_dot_IW_dot_InvSixTeen, InvPhi_1t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_2t_dot_IW_dot_InvSixTeen, InvPhi_2t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_3t_dot_IW_dot_InvSixTeen, InvPhi_3t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_4t_dot_IW_dot_InvSixTeen, InvPhi_4t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_5t_dot_IW_dot_InvSixTeen, InvPhi_5t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_6t_dot_IW_dot_InvSixTeen, InvPhi_6t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_7t_dot_IW_dot_InvSixTeen, InvPhi_7t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_8t_dot_IW_dot_InvSixTeen, InvPhi_8t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_9t_dot_IW_dot_InvSixTeen, InvPhi_9t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_10t_dot_IW_dot_InvSixTeen, InvPhi_10t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_11t_dot_IW_dot_InvSixTeen, InvPhi_11t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_12t_dot_IW_dot_InvSixTeen, InvPhi_12t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_13t_dot_IW_dot_InvSixTeen, InvPhi_13t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_14t_dot_IW_dot_InvSixTeen, InvPhi_14t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_15t_dot_IW_dot_InvSixTeen, InvPhi_15t_dot_IW, InvSixTeen, p);
							//----------------------------------------------------------------
                            MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],InvPhi_0t_dot_IW_dot_InvSixTeen,p);
                            MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],InvPhi_1t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],InvPhi_2t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],InvPhi_3t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R4[ma_tmp],A_B0R4[ma_tmp],InvPhi_4t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R5[ma_tmp],A_B0R5[ma_tmp],InvPhi_5t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R6[ma_tmp],A_B0R6[ma_tmp],InvPhi_6t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R7[ma_tmp],A_B0R7[ma_tmp],InvPhi_7t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R8[ma_tmp],A_B0R8[ma_tmp],InvPhi_8t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R9[ma_tmp],A_B0R9[ma_tmp],InvPhi_9t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R10[ma_tmp],A_B0R10[ma_tmp],InvPhi_10t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R11[ma_tmp],A_B0R11[ma_tmp],InvPhi_11t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R12[ma_tmp],A_B0R12[ma_tmp],InvPhi_12t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R13[ma_tmp],A_B0R13[ma_tmp],InvPhi_13t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R14[ma_tmp],A_B0R14[ma_tmp],InvPhi_14t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R15[ma_tmp],A_B0R15[ma_tmp],InvPhi_15t_dot_IW_dot_InvSixTeen,p);
                            break;
                        case 3:
                            INWC_DATARECORD <<" A_B0R0["<<ma_tmp<<"]:  "<<A_B0R0[ma_tmp]  << ", st3_Tw[0] = " << 1 << endl;
				            INWC_DATARECORD <<" A_B0R1["<<ma_tmp<<"]:  "<<A_B0R1[ma_tmp]  << ", st3_Tw[1] = " << 1 << endl;
					        INWC_DATARECORD <<" A_B0R2["<<ma_tmp<<"]:  "<<A_B0R2[ma_tmp]  << ", st3_Tw[2] = " << 1 << endl;
					        INWC_DATARECORD <<" A_B0R3["<<ma_tmp<<"]:  "<<A_B0R3[ma_tmp]  << ", st3_Tw[3] = " << 1 << endl;
					        INWC_DATARECORD <<" A_B0R4["<<ma_tmp<<"]:  "<<A_B0R4[ma_tmp]  << ", st3_Tw[4] = " << 1 << endl;
					        INWC_DATARECORD <<" A_B0R5["<<ma_tmp<<"]:  "<<A_B0R5[ma_tmp]  << ", st3_Tw[5] = " << 1 << endl;
					        INWC_DATARECORD <<" A_B0R6["<<ma_tmp<<"]:  "<<A_B0R6[ma_tmp]  << ", st3_Tw[6] = " << 1 << endl;
					        INWC_DATARECORD <<" A_B0R7["<<ma_tmp<<"]:  "<<A_B0R7[ma_tmp]  << ", st3_Tw[7] = " << 1 << endl;
					        INWC_DATARECORD <<" A_B0R8["<<ma_tmp<<"]:  "<<A_B0R8[ma_tmp]  << ", st3_Tw[8] = " << 1 << endl;
					        INWC_DATARECORD <<" A_B0R9["<<ma_tmp<<"]:  "<<A_B0R9[ma_tmp]  << ", st3_Tw[9] = " << 1 << endl;
					        INWC_DATARECORD <<" A_B0R10["<<ma_tmp<<"]: "<<A_B0R10[ma_tmp] << ", st3_Tw[10] = " << 1 << endl;   
					        INWC_DATARECORD <<" A_B0R11["<<ma_tmp<<"]: "<<A_B0R11[ma_tmp] << ", st3_Tw[11] = " << 1 << endl;   
					        INWC_DATARECORD <<" A_B0R12["<<ma_tmp<<"]: "<<A_B0R12[ma_tmp] << ", st3_Tw[12] = " << 1 << endl;   
					        INWC_DATARECORD <<" A_B0R13["<<ma_tmp<<"]: "<<A_B0R13[ma_tmp] << ", st3_Tw[13] = " << 1 << endl;   
					        INWC_DATARECORD <<" A_B0R14["<<ma_tmp<<"]: "<<A_B0R14[ma_tmp] << ", st3_Tw[14] = " << 1 << endl;   
					        INWC_DATARECORD <<" A_B0R15["<<ma_tmp<<"]: "<<A_B0R15[ma_tmp] << ", st3_Tw[15] = " << 1 << endl;   
                            if(!debug) INWC_DATARECORD << "InvPhi_0t_Order = " 	<< InvPhi_0t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_0t_Order = " 	<< InvPhi_deg*0 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_1t_Order = " 	<< InvPhi_1t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_1t_Order = " 	<< InvPhi_deg*1 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_2t_Order = " 	<< InvPhi_2t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_2t_Order = " 	<< InvPhi_deg*2 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_3t_Order = " 	<< InvPhi_3t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_3t_Order = " 	<< InvPhi_deg*3 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_4t_Order = " 	<< InvPhi_4t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_4t_Order = " 	<< InvPhi_deg*4 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_5t_Order = " 	<< InvPhi_5t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_5t_Order = " 	<< InvPhi_deg*5 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_6t_Order = " 	<< InvPhi_6t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_6t_Order = " 	<< InvPhi_deg*6 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_7t_Order = " 	<< InvPhi_7t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_7t_Order = " 	<< InvPhi_deg*7 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_8t_Order = " 	<< InvPhi_8t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_8t_Order = " 	<< InvPhi_deg*8 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_9t_Order = " 	<< InvPhi_9t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_9t_Order = " 	<< InvPhi_deg*9 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_10t_Order = " << InvPhi_10t_Order << endl; else INWC_DATARECORD << "InvPhi_10t_Order = " 	<< InvPhi_deg*10 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_11t_Order = " << InvPhi_11t_Order << endl; else INWC_DATARECORD << "InvPhi_11t_Order = " 	<< InvPhi_deg*11 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_12t_Order = " << InvPhi_12t_Order << endl; else INWC_DATARECORD << "InvPhi_12t_Order = " 	<< InvPhi_deg*12 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_13t_Order = " << InvPhi_13t_Order << endl; else INWC_DATARECORD << "InvPhi_13t_Order = " 	<< InvPhi_deg*13 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_14t_Order = " << InvPhi_14t_Order << endl; else INWC_DATARECORD << "InvPhi_14t_Order = " 	<< InvPhi_deg*14 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_15t_Order = " << InvPhi_15t_Order << endl; else INWC_DATARECORD << "InvPhi_15t_Order = " 	<< InvPhi_deg*15 	<< endl;
							INWC_Radix16_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp], A_B0R2[ma_tmp], A_B0R3[ma_tmp],A_B0R4[ma_tmp],
							   A_B0R5[ma_tmp],A_B0R6[ma_tmp], A_B0R7[ma_tmp], A_B0R8[ma_tmp],A_B0R9[ma_tmp],
							   A_B0R10[ma_tmp],A_B0R11[ma_tmp],A_B0R12[ma_tmp],A_B0R13[ma_tmp],A_B0R14[ma_tmp],
							   A_B0R15[ma_tmp]);
							//----------------------compute for INWC--------------------------
                            if(!debug) MulMod(InvPhi_0t_dot_IW, InvPhi_0t_Order, 1, p);
							if(!debug) MulMod(InvPhi_1t_dot_IW, InvPhi_1t_Order, 1, p);
							if(!debug) MulMod(InvPhi_2t_dot_IW, InvPhi_2t_Order, 1, p);
							if(!debug) MulMod(InvPhi_3t_dot_IW, InvPhi_3t_Order, 1, p);
							if(!debug) MulMod(InvPhi_4t_dot_IW, InvPhi_4t_Order, 1, p);
							if(!debug) MulMod(InvPhi_5t_dot_IW, InvPhi_5t_Order, 1, p);
							if(!debug) MulMod(InvPhi_6t_dot_IW, InvPhi_6t_Order, 1, p);
							if(!debug) MulMod(InvPhi_7t_dot_IW, InvPhi_7t_Order, 1, p);
							if(!debug) MulMod(InvPhi_8t_dot_IW, InvPhi_8t_Order, 1, p);
							if(!debug) MulMod(InvPhi_9t_dot_IW, InvPhi_9t_Order, 1, p);
							if(!debug) MulMod(InvPhi_10t_dot_IW, InvPhi_10t_Order, 1, p);
							if(!debug) MulMod(InvPhi_11t_dot_IW, InvPhi_11t_Order, 1, p);
							if(!debug) MulMod(InvPhi_12t_dot_IW, InvPhi_12t_Order, 1, p);
							if(!debug) MulMod(InvPhi_13t_dot_IW, InvPhi_13t_Order, 1, p);
							if(!debug) MulMod(InvPhi_14t_dot_IW, InvPhi_14t_Order, 1, p);
							if(!debug) MulMod(InvPhi_15t_dot_IW, InvPhi_15t_Order, 1, p);

							if(!debug) MulMod(InvPhi_0t_dot_IW_dot_InvSixTeen, InvPhi_0t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_1t_dot_IW_dot_InvSixTeen, InvPhi_1t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_2t_dot_IW_dot_InvSixTeen, InvPhi_2t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_3t_dot_IW_dot_InvSixTeen, InvPhi_3t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_4t_dot_IW_dot_InvSixTeen, InvPhi_4t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_5t_dot_IW_dot_InvSixTeen, InvPhi_5t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_6t_dot_IW_dot_InvSixTeen, InvPhi_6t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_7t_dot_IW_dot_InvSixTeen, InvPhi_7t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_8t_dot_IW_dot_InvSixTeen, InvPhi_8t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_9t_dot_IW_dot_InvSixTeen, InvPhi_9t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_10t_dot_IW_dot_InvSixTeen, InvPhi_10t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_11t_dot_IW_dot_InvSixTeen, InvPhi_11t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_12t_dot_IW_dot_InvSixTeen, InvPhi_12t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_13t_dot_IW_dot_InvSixTeen, InvPhi_13t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_14t_dot_IW_dot_InvSixTeen, InvPhi_14t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_15t_dot_IW_dot_InvSixTeen, InvPhi_15t_dot_IW, InvSixTeen, p);
							//----------------------------------------------------------------
                            MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],InvPhi_0t_dot_IW_dot_InvSixTeen,p);
                            MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],InvPhi_1t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],InvPhi_2t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],InvPhi_3t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R4[ma_tmp],A_B0R4[ma_tmp],InvPhi_4t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R5[ma_tmp],A_B0R5[ma_tmp],InvPhi_5t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R6[ma_tmp],A_B0R6[ma_tmp],InvPhi_6t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R7[ma_tmp],A_B0R7[ma_tmp],InvPhi_7t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R8[ma_tmp],A_B0R8[ma_tmp],InvPhi_8t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R9[ma_tmp],A_B0R9[ma_tmp],InvPhi_9t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R10[ma_tmp],A_B0R10[ma_tmp],InvPhi_10t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R11[ma_tmp],A_B0R11[ma_tmp],InvPhi_11t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R12[ma_tmp],A_B0R12[ma_tmp],InvPhi_12t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R13[ma_tmp],A_B0R13[ma_tmp],InvPhi_13t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R14[ma_tmp],A_B0R14[ma_tmp],InvPhi_14t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B0R15[ma_tmp],A_B0R15[ma_tmp],InvPhi_15t_dot_IW_dot_InvSixTeen,p);
                            break;
                    }
                    INWC_DATARECORD <<" ------after BU compute and Mul-------" << std::endl;
				    INWC_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp]<<"\n";
				    INWC_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp]<<"\n";
					INWC_DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<<A_B0R2[ma_tmp]<<"\n";
					INWC_DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<<A_B0R3[ma_tmp]<<"\n";
					INWC_DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<<A_B0R4[ma_tmp]<<"\n";
					INWC_DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<<A_B0R5[ma_tmp]<<"\n";
					INWC_DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<<A_B0R6[ma_tmp]<<"\n";
					INWC_DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<<A_B0R7[ma_tmp]<<"\n";
					INWC_DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<<A_B0R8[ma_tmp]<<"\n";
					INWC_DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<<A_B0R9[ma_tmp]<<"\n";
					INWC_DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<<A_B0R10[ma_tmp]<<"\n";
					INWC_DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<<A_B0R11[ma_tmp]<<"\n";
					INWC_DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<<A_B0R12[ma_tmp]<<"\n";
					INWC_DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<<A_B0R13[ma_tmp]<<"\n";
					INWC_DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<<A_B0R14[ma_tmp]<<"\n";
					INWC_DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<<A_B0R15[ma_tmp]<<"\n";
					INWC_DATARECORD <<"--------------------------------------------------------------------\n";
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
                    switch(s){
                        case 0:
                            INWC_DATARECORD <<" A_B1R0["<<ma_tmp<<"]:  "<<A_B1R0[ma_tmp]  << ", st0_Tw[0] = " << st0_Tw[0] << endl;     
				            INWC_DATARECORD <<" A_B1R1["<<ma_tmp<<"]:  "<<A_B1R1[ma_tmp]  << ", st0_Tw[1] = " << st0_Tw[1] << endl;     
					        INWC_DATARECORD <<" A_B1R2["<<ma_tmp<<"]:  "<<A_B1R2[ma_tmp]  << ", st0_Tw[2] = " << st0_Tw[2] << endl;     
					        INWC_DATARECORD <<" A_B1R3["<<ma_tmp<<"]:  "<<A_B1R3[ma_tmp]  << ", st0_Tw[3] = " << st0_Tw[3] << endl;     
					        INWC_DATARECORD <<" A_B1R4["<<ma_tmp<<"]:  "<<A_B1R4[ma_tmp]  << ", st0_Tw[4] = " << st0_Tw[4] << endl;     
					        INWC_DATARECORD <<" A_B1R5["<<ma_tmp<<"]:  "<<A_B1R5[ma_tmp]  << ", st0_Tw[5] = " << st0_Tw[5] << endl;     
					        INWC_DATARECORD <<" A_B1R6["<<ma_tmp<<"]:  "<<A_B1R6[ma_tmp]  << ", st0_Tw[6] = " << st0_Tw[6] << endl;     
					        INWC_DATARECORD <<" A_B1R7["<<ma_tmp<<"]:  "<<A_B1R7[ma_tmp]  << ", st0_Tw[7] = " << st0_Tw[7] << endl;     
					        INWC_DATARECORD <<" A_B1R8["<<ma_tmp<<"]:  "<<A_B1R8[ma_tmp]  << ", st0_Tw[8] = " << st0_Tw[8] << endl;     
					        INWC_DATARECORD <<" A_B1R9["<<ma_tmp<<"]:  "<<A_B1R9[ma_tmp]  << ", st0_Tw[9] = " << st0_Tw[9] << endl;     
					        INWC_DATARECORD <<" A_B1R10["<<ma_tmp<<"]: "<<A_B1R10[ma_tmp] << ", st0_Tw[10] = " << st0_Tw[10] << endl;    
					        INWC_DATARECORD <<" A_B1R11["<<ma_tmp<<"]: "<<A_B1R11[ma_tmp] << ", st0_Tw[11] = " << st0_Tw[11] << endl;    
					        INWC_DATARECORD <<" A_B1R12["<<ma_tmp<<"]: "<<A_B1R12[ma_tmp] << ", st0_Tw[12] = " << st0_Tw[12] << endl;    
					        INWC_DATARECORD <<" A_B1R13["<<ma_tmp<<"]: "<<A_B1R13[ma_tmp] << ", st0_Tw[13] = " << st0_Tw[13] << endl;    
					        INWC_DATARECORD <<" A_B1R14["<<ma_tmp<<"]: "<<A_B1R14[ma_tmp] << ", st0_Tw[14] = " << st0_Tw[14] << endl;    
					        INWC_DATARECORD <<" A_B1R15["<<ma_tmp<<"]: "<<A_B1R15[ma_tmp] << ", st0_Tw[15] = " << st0_Tw[15] << endl;    
                            if(!debug) INWC_DATARECORD << "InvPhi_0t_Order = " 	<< InvPhi_0t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_0t_Order = " 	<< InvPhi_deg*0 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_1t_Order = " 	<< InvPhi_1t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_1t_Order = " 	<< InvPhi_deg*1 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_2t_Order = " 	<< InvPhi_2t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_2t_Order = " 	<< InvPhi_deg*2 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_3t_Order = " 	<< InvPhi_3t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_3t_Order = " 	<< InvPhi_deg*3 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_4t_Order = " 	<< InvPhi_4t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_4t_Order = " 	<< InvPhi_deg*4 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_5t_Order = " 	<< InvPhi_5t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_5t_Order = " 	<< InvPhi_deg*5 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_6t_Order = " 	<< InvPhi_6t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_6t_Order = " 	<< InvPhi_deg*6 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_7t_Order = " 	<< InvPhi_7t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_7t_Order = " 	<< InvPhi_deg*7 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_8t_Order = " 	<< InvPhi_8t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_8t_Order = " 	<< InvPhi_deg*8 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_9t_Order = " 	<< InvPhi_9t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_9t_Order = " 	<< InvPhi_deg*9 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_10t_Order = " << InvPhi_10t_Order << endl; else INWC_DATARECORD << "InvPhi_10t_Order = " 	<< InvPhi_deg*10 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_11t_Order = " << InvPhi_11t_Order << endl; else INWC_DATARECORD << "InvPhi_11t_Order = " 	<< InvPhi_deg*11 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_12t_Order = " << InvPhi_12t_Order << endl; else INWC_DATARECORD << "InvPhi_12t_Order = " 	<< InvPhi_deg*12 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_13t_Order = " << InvPhi_13t_Order << endl; else INWC_DATARECORD << "InvPhi_13t_Order = " 	<< InvPhi_deg*13 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_14t_Order = " << InvPhi_14t_Order << endl; else INWC_DATARECORD << "InvPhi_14t_Order = " 	<< InvPhi_deg*14 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_15t_Order = " << InvPhi_15t_Order << endl; else INWC_DATARECORD << "InvPhi_15t_Order = " 	<< InvPhi_deg*15 	<< endl;
							INWC_Radix16_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp], A_B1R2[ma_tmp], A_B1R3[ma_tmp],A_B1R4[ma_tmp],
					        		   A_B1R5[ma_tmp],A_B1R6[ma_tmp], A_B1R7[ma_tmp], A_B1R8[ma_tmp],A_B1R9[ma_tmp],
					        		   A_B1R10[ma_tmp],A_B1R11[ma_tmp],A_B1R12[ma_tmp],A_B1R13[ma_tmp],A_B1R14[ma_tmp],
					        		   A_B1R15[ma_tmp]);
							//----------------------compute for INWC--------------------------
                            if(!debug) MulMod(InvPhi_0t_dot_IW, InvPhi_0t_Order, st0_Tw[0], p);
							if(!debug) MulMod(InvPhi_1t_dot_IW, InvPhi_1t_Order, st0_Tw[1], p);
							if(!debug) MulMod(InvPhi_2t_dot_IW, InvPhi_2t_Order, st0_Tw[2], p);
							if(!debug) MulMod(InvPhi_3t_dot_IW, InvPhi_3t_Order, st0_Tw[3], p);
							if(!debug) MulMod(InvPhi_4t_dot_IW, InvPhi_4t_Order, st0_Tw[4], p);
							if(!debug) MulMod(InvPhi_5t_dot_IW, InvPhi_5t_Order, st0_Tw[5], p);
							if(!debug) MulMod(InvPhi_6t_dot_IW, InvPhi_6t_Order, st0_Tw[6], p);
							if(!debug) MulMod(InvPhi_7t_dot_IW, InvPhi_7t_Order, st0_Tw[7], p);
							if(!debug) MulMod(InvPhi_8t_dot_IW, InvPhi_8t_Order, st0_Tw[8], p);
							if(!debug) MulMod(InvPhi_9t_dot_IW, InvPhi_9t_Order, st0_Tw[9], p);
							if(!debug) MulMod(InvPhi_10t_dot_IW, InvPhi_10t_Order, st0_Tw[10], p);
							if(!debug) MulMod(InvPhi_11t_dot_IW, InvPhi_11t_Order, st0_Tw[11], p);
							if(!debug) MulMod(InvPhi_12t_dot_IW, InvPhi_12t_Order, st0_Tw[12], p);
							if(!debug) MulMod(InvPhi_13t_dot_IW, InvPhi_13t_Order, st0_Tw[13], p);
							if(!debug) MulMod(InvPhi_14t_dot_IW, InvPhi_14t_Order, st0_Tw[14], p);
							if(!debug) MulMod(InvPhi_15t_dot_IW, InvPhi_15t_Order, st0_Tw[15], p);

							if(!debug) MulMod(InvPhi_0t_dot_IW_dot_InvSixTeen, InvPhi_0t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_1t_dot_IW_dot_InvSixTeen, InvPhi_1t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_2t_dot_IW_dot_InvSixTeen, InvPhi_2t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_3t_dot_IW_dot_InvSixTeen, InvPhi_3t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_4t_dot_IW_dot_InvSixTeen, InvPhi_4t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_5t_dot_IW_dot_InvSixTeen, InvPhi_5t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_6t_dot_IW_dot_InvSixTeen, InvPhi_6t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_7t_dot_IW_dot_InvSixTeen, InvPhi_7t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_8t_dot_IW_dot_InvSixTeen, InvPhi_8t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_9t_dot_IW_dot_InvSixTeen, InvPhi_9t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_10t_dot_IW_dot_InvSixTeen, InvPhi_10t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_11t_dot_IW_dot_InvSixTeen, InvPhi_11t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_12t_dot_IW_dot_InvSixTeen, InvPhi_12t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_13t_dot_IW_dot_InvSixTeen, InvPhi_13t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_14t_dot_IW_dot_InvSixTeen, InvPhi_14t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_15t_dot_IW_dot_InvSixTeen, InvPhi_15t_dot_IW, InvSixTeen, p);
							//----------------------------------------------------------------
                            MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],InvPhi_0t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],InvPhi_1t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],InvPhi_2t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],InvPhi_3t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R4[ma_tmp],A_B1R4[ma_tmp],InvPhi_4t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R5[ma_tmp],A_B1R5[ma_tmp],InvPhi_5t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R6[ma_tmp],A_B1R6[ma_tmp],InvPhi_6t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R7[ma_tmp],A_B1R7[ma_tmp],InvPhi_7t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R8[ma_tmp],A_B1R8[ma_tmp],InvPhi_8t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R9[ma_tmp],A_B1R9[ma_tmp],InvPhi_9t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R10[ma_tmp],A_B1R10[ma_tmp],InvPhi_10t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R11[ma_tmp],A_B1R11[ma_tmp],InvPhi_11t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R12[ma_tmp],A_B1R12[ma_tmp],InvPhi_12t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R13[ma_tmp],A_B1R13[ma_tmp],InvPhi_13t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R14[ma_tmp],A_B1R14[ma_tmp],InvPhi_14t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R15[ma_tmp],A_B1R15[ma_tmp],InvPhi_15t_dot_IW_dot_InvSixTeen,p);
                            break;
                        case 1:
                            INWC_DATARECORD <<" A_B1R0["<<ma_tmp<<"]:  "<<A_B1R0[ma_tmp]  << ", st1_Tw[0] = " << st1_Tw[0] << endl;     
				            INWC_DATARECORD <<" A_B1R1["<<ma_tmp<<"]:  "<<A_B1R1[ma_tmp]  << ", st1_Tw[1] = " << st1_Tw[1] << endl;     
					        INWC_DATARECORD <<" A_B1R2["<<ma_tmp<<"]:  "<<A_B1R2[ma_tmp]  << ", st1_Tw[2] = " << st1_Tw[2] << endl;     
					        INWC_DATARECORD <<" A_B1R3["<<ma_tmp<<"]:  "<<A_B1R3[ma_tmp]  << ", st1_Tw[3] = " << st1_Tw[3] << endl;     
					        INWC_DATARECORD <<" A_B1R4["<<ma_tmp<<"]:  "<<A_B1R4[ma_tmp]  << ", st1_Tw[4] = " << st1_Tw[4] << endl;     
					        INWC_DATARECORD <<" A_B1R5["<<ma_tmp<<"]:  "<<A_B1R5[ma_tmp]  << ", st1_Tw[5] = " << st1_Tw[5] << endl;     
					        INWC_DATARECORD <<" A_B1R6["<<ma_tmp<<"]:  "<<A_B1R6[ma_tmp]  << ", st1_Tw[6] = " << st1_Tw[6] << endl;     
					        INWC_DATARECORD <<" A_B1R7["<<ma_tmp<<"]:  "<<A_B1R7[ma_tmp]  << ", st1_Tw[7] = " << st1_Tw[7] << endl;     
					        INWC_DATARECORD <<" A_B1R8["<<ma_tmp<<"]:  "<<A_B1R8[ma_tmp]  << ", st1_Tw[8] = " << st1_Tw[8] << endl;     
					        INWC_DATARECORD <<" A_B1R9["<<ma_tmp<<"]:  "<<A_B1R9[ma_tmp]  << ", st1_Tw[9] = " << st1_Tw[9] << endl;     
					        INWC_DATARECORD <<" A_B1R10["<<ma_tmp<<"]: "<<A_B1R10[ma_tmp] << ", st1_Tw[10] = " << st1_Tw[10] << endl;    
					        INWC_DATARECORD <<" A_B1R11["<<ma_tmp<<"]: "<<A_B1R11[ma_tmp] << ", st1_Tw[11] = " << st1_Tw[11] << endl;    
					        INWC_DATARECORD <<" A_B1R12["<<ma_tmp<<"]: "<<A_B1R12[ma_tmp] << ", st1_Tw[12] = " << st1_Tw[12] << endl;    
					        INWC_DATARECORD <<" A_B1R13["<<ma_tmp<<"]: "<<A_B1R13[ma_tmp] << ", st1_Tw[13] = " << st1_Tw[13] << endl;    
					        INWC_DATARECORD <<" A_B1R14["<<ma_tmp<<"]: "<<A_B1R14[ma_tmp] << ", st1_Tw[14] = " << st1_Tw[14] << endl;    
					        INWC_DATARECORD <<" A_B1R15["<<ma_tmp<<"]: "<<A_B1R15[ma_tmp] << ", st1_Tw[15] = " << st1_Tw[15] << endl;    
                            if(!debug) INWC_DATARECORD << "InvPhi_0t_Order = " 	<< InvPhi_0t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_0t_Order = " 	<< InvPhi_deg*0 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_1t_Order = " 	<< InvPhi_1t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_1t_Order = " 	<< InvPhi_deg*1 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_2t_Order = " 	<< InvPhi_2t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_2t_Order = " 	<< InvPhi_deg*2 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_3t_Order = " 	<< InvPhi_3t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_3t_Order = " 	<< InvPhi_deg*3 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_4t_Order = " 	<< InvPhi_4t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_4t_Order = " 	<< InvPhi_deg*4 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_5t_Order = " 	<< InvPhi_5t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_5t_Order = " 	<< InvPhi_deg*5 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_6t_Order = " 	<< InvPhi_6t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_6t_Order = " 	<< InvPhi_deg*6 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_7t_Order = " 	<< InvPhi_7t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_7t_Order = " 	<< InvPhi_deg*7 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_8t_Order = " 	<< InvPhi_8t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_8t_Order = " 	<< InvPhi_deg*8 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_9t_Order = " 	<< InvPhi_9t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_9t_Order = " 	<< InvPhi_deg*9 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_10t_Order = " << InvPhi_10t_Order << endl; else INWC_DATARECORD << "InvPhi_10t_Order = " 	<< InvPhi_deg*10 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_11t_Order = " << InvPhi_11t_Order << endl; else INWC_DATARECORD << "InvPhi_11t_Order = " 	<< InvPhi_deg*11 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_12t_Order = " << InvPhi_12t_Order << endl; else INWC_DATARECORD << "InvPhi_12t_Order = " 	<< InvPhi_deg*12 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_13t_Order = " << InvPhi_13t_Order << endl; else INWC_DATARECORD << "InvPhi_13t_Order = " 	<< InvPhi_deg*13 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_14t_Order = " << InvPhi_14t_Order << endl; else INWC_DATARECORD << "InvPhi_14t_Order = " 	<< InvPhi_deg*14 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_15t_Order = " << InvPhi_15t_Order << endl; else INWC_DATARECORD << "InvPhi_15t_Order = " 	<< InvPhi_deg*15 	<< endl;
							INWC_Radix16_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp], A_B1R2[ma_tmp], A_B1R3[ma_tmp],A_B1R4[ma_tmp],
					        		   A_B1R5[ma_tmp],A_B1R6[ma_tmp], A_B1R7[ma_tmp], A_B1R8[ma_tmp],A_B1R9[ma_tmp],
					        		   A_B1R10[ma_tmp],A_B1R11[ma_tmp],A_B1R12[ma_tmp],A_B1R13[ma_tmp],A_B1R14[ma_tmp],
					        		   A_B1R15[ma_tmp]);
							//----------------------compute for INWC--------------------------
                            if(!debug) MulMod(InvPhi_0t_dot_IW, InvPhi_0t_Order, st1_Tw[0], p);
							if(!debug) MulMod(InvPhi_1t_dot_IW, InvPhi_1t_Order, st1_Tw[1], p);
							if(!debug) MulMod(InvPhi_2t_dot_IW, InvPhi_2t_Order, st1_Tw[2], p);
							if(!debug) MulMod(InvPhi_3t_dot_IW, InvPhi_3t_Order, st1_Tw[3], p);
							if(!debug) MulMod(InvPhi_4t_dot_IW, InvPhi_4t_Order, st1_Tw[4], p);
							if(!debug) MulMod(InvPhi_5t_dot_IW, InvPhi_5t_Order, st1_Tw[5], p);
							if(!debug) MulMod(InvPhi_6t_dot_IW, InvPhi_6t_Order, st1_Tw[6], p);
							if(!debug) MulMod(InvPhi_7t_dot_IW, InvPhi_7t_Order, st1_Tw[7], p);
							if(!debug) MulMod(InvPhi_8t_dot_IW, InvPhi_8t_Order, st1_Tw[8], p);
							if(!debug) MulMod(InvPhi_9t_dot_IW, InvPhi_9t_Order, st1_Tw[9], p);
							if(!debug) MulMod(InvPhi_10t_dot_IW, InvPhi_10t_Order, st1_Tw[10], p);
							if(!debug) MulMod(InvPhi_11t_dot_IW, InvPhi_11t_Order, st1_Tw[11], p);
							if(!debug) MulMod(InvPhi_12t_dot_IW, InvPhi_12t_Order, st1_Tw[12], p);
							if(!debug) MulMod(InvPhi_13t_dot_IW, InvPhi_13t_Order, st1_Tw[13], p);
							if(!debug) MulMod(InvPhi_14t_dot_IW, InvPhi_14t_Order, st1_Tw[14], p);
							if(!debug) MulMod(InvPhi_15t_dot_IW, InvPhi_15t_Order, st1_Tw[15], p);

							if(!debug) MulMod(InvPhi_0t_dot_IW_dot_InvSixTeen, InvPhi_0t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_1t_dot_IW_dot_InvSixTeen, InvPhi_1t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_2t_dot_IW_dot_InvSixTeen, InvPhi_2t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_3t_dot_IW_dot_InvSixTeen, InvPhi_3t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_4t_dot_IW_dot_InvSixTeen, InvPhi_4t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_5t_dot_IW_dot_InvSixTeen, InvPhi_5t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_6t_dot_IW_dot_InvSixTeen, InvPhi_6t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_7t_dot_IW_dot_InvSixTeen, InvPhi_7t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_8t_dot_IW_dot_InvSixTeen, InvPhi_8t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_9t_dot_IW_dot_InvSixTeen, InvPhi_9t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_10t_dot_IW_dot_InvSixTeen, InvPhi_10t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_11t_dot_IW_dot_InvSixTeen, InvPhi_11t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_12t_dot_IW_dot_InvSixTeen, InvPhi_12t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_13t_dot_IW_dot_InvSixTeen, InvPhi_13t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_14t_dot_IW_dot_InvSixTeen, InvPhi_14t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_15t_dot_IW_dot_InvSixTeen, InvPhi_15t_dot_IW, InvSixTeen, p);
							//----------------------------------------------------------------
                            MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],InvPhi_0t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],InvPhi_1t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],InvPhi_2t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],InvPhi_3t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R4[ma_tmp],A_B1R4[ma_tmp],InvPhi_4t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R5[ma_tmp],A_B1R5[ma_tmp],InvPhi_5t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R6[ma_tmp],A_B1R6[ma_tmp],InvPhi_6t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R7[ma_tmp],A_B1R7[ma_tmp],InvPhi_7t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R8[ma_tmp],A_B1R8[ma_tmp],InvPhi_8t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R9[ma_tmp],A_B1R9[ma_tmp],InvPhi_9t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R10[ma_tmp],A_B1R10[ma_tmp],InvPhi_10t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R11[ma_tmp],A_B1R11[ma_tmp],InvPhi_11t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R12[ma_tmp],A_B1R12[ma_tmp],InvPhi_12t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R13[ma_tmp],A_B1R13[ma_tmp],InvPhi_13t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R14[ma_tmp],A_B1R14[ma_tmp],InvPhi_14t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R15[ma_tmp],A_B1R15[ma_tmp],InvPhi_15t_dot_IW_dot_InvSixTeen,p);
                            break;
                        case 2:
                            INWC_DATARECORD <<" A_B1R0["<<ma_tmp<<"]:  "<<A_B1R0[ma_tmp]  << ", st2_Tw[0] = " << st2_Tw[0] << endl;     
				            INWC_DATARECORD <<" A_B1R1["<<ma_tmp<<"]:  "<<A_B1R1[ma_tmp]  << ", st2_Tw[1] = " << st2_Tw[1] << endl;     
					        INWC_DATARECORD <<" A_B1R2["<<ma_tmp<<"]:  "<<A_B1R2[ma_tmp]  << ", st2_Tw[2] = " << st2_Tw[2] << endl;     
					        INWC_DATARECORD <<" A_B1R3["<<ma_tmp<<"]:  "<<A_B1R3[ma_tmp]  << ", st2_Tw[3] = " << st2_Tw[3] << endl;     
					        INWC_DATARECORD <<" A_B1R4["<<ma_tmp<<"]:  "<<A_B1R4[ma_tmp]  << ", st2_Tw[4] = " << st2_Tw[4] << endl;     
					        INWC_DATARECORD <<" A_B1R5["<<ma_tmp<<"]:  "<<A_B1R5[ma_tmp]  << ", st2_Tw[5] = " << st2_Tw[5] << endl;     
					        INWC_DATARECORD <<" A_B1R6["<<ma_tmp<<"]:  "<<A_B1R6[ma_tmp]  << ", st2_Tw[6] = " << st2_Tw[6] << endl;     
					        INWC_DATARECORD <<" A_B1R7["<<ma_tmp<<"]:  "<<A_B1R7[ma_tmp]  << ", st2_Tw[7] = " << st2_Tw[7] << endl;     
					        INWC_DATARECORD <<" A_B1R8["<<ma_tmp<<"]:  "<<A_B1R8[ma_tmp]  << ", st2_Tw[8] = " << st2_Tw[8] << endl;     
					        INWC_DATARECORD <<" A_B1R9["<<ma_tmp<<"]:  "<<A_B1R9[ma_tmp]  << ", st2_Tw[9] = " << st2_Tw[9] << endl;     
					        INWC_DATARECORD <<" A_B1R10["<<ma_tmp<<"]: "<<A_B1R10[ma_tmp] << ", st2_Tw[10] = " << st2_Tw[10] << endl;    
					        INWC_DATARECORD <<" A_B1R11["<<ma_tmp<<"]: "<<A_B1R11[ma_tmp] << ", st2_Tw[11] = " << st2_Tw[11] << endl;    
					        INWC_DATARECORD <<" A_B1R12["<<ma_tmp<<"]: "<<A_B1R12[ma_tmp] << ", st2_Tw[12] = " << st2_Tw[12] << endl;    
					        INWC_DATARECORD <<" A_B1R13["<<ma_tmp<<"]: "<<A_B1R13[ma_tmp] << ", st2_Tw[13] = " << st2_Tw[13] << endl;    
					        INWC_DATARECORD <<" A_B1R14["<<ma_tmp<<"]: "<<A_B1R14[ma_tmp] << ", st2_Tw[14] = " << st2_Tw[14] << endl;    
					        INWC_DATARECORD <<" A_B1R15["<<ma_tmp<<"]: "<<A_B1R15[ma_tmp] << ", st2_Tw[15] = " << st2_Tw[15] << endl;    
                            if(!debug) INWC_DATARECORD << "InvPhi_0t_Order = " 	<< InvPhi_0t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_0t_Order = " 	<< InvPhi_deg*0 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_1t_Order = " 	<< InvPhi_1t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_1t_Order = " 	<< InvPhi_deg*1 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_2t_Order = " 	<< InvPhi_2t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_2t_Order = " 	<< InvPhi_deg*2 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_3t_Order = " 	<< InvPhi_3t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_3t_Order = " 	<< InvPhi_deg*3 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_4t_Order = " 	<< InvPhi_4t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_4t_Order = " 	<< InvPhi_deg*4 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_5t_Order = " 	<< InvPhi_5t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_5t_Order = " 	<< InvPhi_deg*5 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_6t_Order = " 	<< InvPhi_6t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_6t_Order = " 	<< InvPhi_deg*6 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_7t_Order = " 	<< InvPhi_7t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_7t_Order = " 	<< InvPhi_deg*7 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_8t_Order = " 	<< InvPhi_8t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_8t_Order = " 	<< InvPhi_deg*8 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_9t_Order = " 	<< InvPhi_9t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_9t_Order = " 	<< InvPhi_deg*9 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_10t_Order = " << InvPhi_10t_Order << endl; else INWC_DATARECORD << "InvPhi_10t_Order = " 	<< InvPhi_deg*10 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_11t_Order = " << InvPhi_11t_Order << endl; else INWC_DATARECORD << "InvPhi_11t_Order = " 	<< InvPhi_deg*11 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_12t_Order = " << InvPhi_12t_Order << endl; else INWC_DATARECORD << "InvPhi_12t_Order = " 	<< InvPhi_deg*12 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_13t_Order = " << InvPhi_13t_Order << endl; else INWC_DATARECORD << "InvPhi_13t_Order = " 	<< InvPhi_deg*13 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_14t_Order = " << InvPhi_14t_Order << endl; else INWC_DATARECORD << "InvPhi_14t_Order = " 	<< InvPhi_deg*14 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_15t_Order = " << InvPhi_15t_Order << endl; else INWC_DATARECORD << "InvPhi_15t_Order = " 	<< InvPhi_deg*15 	<< endl;
							INWC_Radix16_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp], A_B1R2[ma_tmp], A_B1R3[ma_tmp],A_B1R4[ma_tmp],
					        		   A_B1R5[ma_tmp],A_B1R6[ma_tmp], A_B1R7[ma_tmp], A_B1R8[ma_tmp],A_B1R9[ma_tmp],
					        		   A_B1R10[ma_tmp],A_B1R11[ma_tmp],A_B1R12[ma_tmp],A_B1R13[ma_tmp],A_B1R14[ma_tmp],
					        		   A_B1R15[ma_tmp]);
							//----------------------compute for INWC--------------------------
                            if(!debug) MulMod(InvPhi_0t_dot_IW, InvPhi_0t_Order, st2_Tw[0], p);
							if(!debug) MulMod(InvPhi_1t_dot_IW, InvPhi_1t_Order, st2_Tw[1], p);
							if(!debug) MulMod(InvPhi_2t_dot_IW, InvPhi_2t_Order, st2_Tw[2], p);
							if(!debug) MulMod(InvPhi_3t_dot_IW, InvPhi_3t_Order, st2_Tw[3], p);
							if(!debug) MulMod(InvPhi_4t_dot_IW, InvPhi_4t_Order, st2_Tw[4], p);
							if(!debug) MulMod(InvPhi_5t_dot_IW, InvPhi_5t_Order, st2_Tw[5], p);
							if(!debug) MulMod(InvPhi_6t_dot_IW, InvPhi_6t_Order, st2_Tw[6], p);
							if(!debug) MulMod(InvPhi_7t_dot_IW, InvPhi_7t_Order, st2_Tw[7], p);
							if(!debug) MulMod(InvPhi_8t_dot_IW, InvPhi_8t_Order, st2_Tw[8], p);
							if(!debug) MulMod(InvPhi_9t_dot_IW, InvPhi_9t_Order, st2_Tw[9], p);
							if(!debug) MulMod(InvPhi_10t_dot_IW, InvPhi_10t_Order, st2_Tw[10], p);
							if(!debug) MulMod(InvPhi_11t_dot_IW, InvPhi_11t_Order, st2_Tw[11], p);
							if(!debug) MulMod(InvPhi_12t_dot_IW, InvPhi_12t_Order, st2_Tw[12], p);
							if(!debug) MulMod(InvPhi_13t_dot_IW, InvPhi_13t_Order, st2_Tw[13], p);
							if(!debug) MulMod(InvPhi_14t_dot_IW, InvPhi_14t_Order, st2_Tw[14], p);
							if(!debug) MulMod(InvPhi_15t_dot_IW, InvPhi_15t_Order, st2_Tw[15], p);

							if(!debug) MulMod(InvPhi_0t_dot_IW_dot_InvSixTeen, InvPhi_0t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_1t_dot_IW_dot_InvSixTeen, InvPhi_1t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_2t_dot_IW_dot_InvSixTeen, InvPhi_2t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_3t_dot_IW_dot_InvSixTeen, InvPhi_3t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_4t_dot_IW_dot_InvSixTeen, InvPhi_4t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_5t_dot_IW_dot_InvSixTeen, InvPhi_5t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_6t_dot_IW_dot_InvSixTeen, InvPhi_6t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_7t_dot_IW_dot_InvSixTeen, InvPhi_7t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_8t_dot_IW_dot_InvSixTeen, InvPhi_8t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_9t_dot_IW_dot_InvSixTeen, InvPhi_9t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_10t_dot_IW_dot_InvSixTeen, InvPhi_10t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_11t_dot_IW_dot_InvSixTeen, InvPhi_11t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_12t_dot_IW_dot_InvSixTeen, InvPhi_12t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_13t_dot_IW_dot_InvSixTeen, InvPhi_13t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_14t_dot_IW_dot_InvSixTeen, InvPhi_14t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_15t_dot_IW_dot_InvSixTeen, InvPhi_15t_dot_IW, InvSixTeen, p);
							//----------------------------------------------------------------
                            MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],InvPhi_0t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],InvPhi_1t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],InvPhi_2t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],InvPhi_3t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R4[ma_tmp],A_B1R4[ma_tmp],InvPhi_4t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R5[ma_tmp],A_B1R5[ma_tmp],InvPhi_5t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R6[ma_tmp],A_B1R6[ma_tmp],InvPhi_6t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R7[ma_tmp],A_B1R7[ma_tmp],InvPhi_7t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R8[ma_tmp],A_B1R8[ma_tmp],InvPhi_8t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R9[ma_tmp],A_B1R9[ma_tmp],InvPhi_9t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R10[ma_tmp],A_B1R10[ma_tmp],InvPhi_10t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R11[ma_tmp],A_B1R11[ma_tmp],InvPhi_11t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R12[ma_tmp],A_B1R12[ma_tmp],InvPhi_12t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R13[ma_tmp],A_B1R13[ma_tmp],InvPhi_13t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R14[ma_tmp],A_B1R14[ma_tmp],InvPhi_14t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R15[ma_tmp],A_B1R15[ma_tmp],InvPhi_15t_dot_IW_dot_InvSixTeen,p);
                            break;
                        case 3:
                            INWC_DATARECORD <<" A_B1R0["<<ma_tmp<<"]:  "<<A_B1R0[ma_tmp]  << ", st3_Tw[0] = " << 1 << endl;     
				            INWC_DATARECORD <<" A_B1R1["<<ma_tmp<<"]:  "<<A_B1R1[ma_tmp]  << ", st3_Tw[1] = " << 1 << endl;     
					        INWC_DATARECORD <<" A_B1R2["<<ma_tmp<<"]:  "<<A_B1R2[ma_tmp]  << ", st3_Tw[2] = " << 1 << endl;     
					        INWC_DATARECORD <<" A_B1R3["<<ma_tmp<<"]:  "<<A_B1R3[ma_tmp]  << ", st3_Tw[3] = " << 1 << endl;     
					        INWC_DATARECORD <<" A_B1R4["<<ma_tmp<<"]:  "<<A_B1R4[ma_tmp]  << ", st3_Tw[4] = " << 1 << endl;     
					        INWC_DATARECORD <<" A_B1R5["<<ma_tmp<<"]:  "<<A_B1R5[ma_tmp]  << ", st3_Tw[5] = " << 1 << endl;     
					        INWC_DATARECORD <<" A_B1R6["<<ma_tmp<<"]:  "<<A_B1R6[ma_tmp]  << ", st3_Tw[6] = " << 1 << endl;     
					        INWC_DATARECORD <<" A_B1R7["<<ma_tmp<<"]:  "<<A_B1R7[ma_tmp]  << ", st3_Tw[7] = " << 1 << endl;     
					        INWC_DATARECORD <<" A_B1R8["<<ma_tmp<<"]:  "<<A_B1R8[ma_tmp]  << ", st3_Tw[8] = " << 1 << endl;     
					        INWC_DATARECORD <<" A_B1R9["<<ma_tmp<<"]:  "<<A_B1R9[ma_tmp]  << ", st3_Tw[9] = " << 1 << endl;     
					        INWC_DATARECORD <<" A_B1R10["<<ma_tmp<<"]: "<<A_B1R10[ma_tmp] << ", st3_Tw[10] = " << 1 << endl;    
					        INWC_DATARECORD <<" A_B1R11["<<ma_tmp<<"]: "<<A_B1R11[ma_tmp] << ", st3_Tw[11] = " << 1 << endl;    
					        INWC_DATARECORD <<" A_B1R12["<<ma_tmp<<"]: "<<A_B1R12[ma_tmp] << ", st3_Tw[12] = " << 1 << endl;    
					        INWC_DATARECORD <<" A_B1R13["<<ma_tmp<<"]: "<<A_B1R13[ma_tmp] << ", st3_Tw[13] = " << 1 << endl;    
					        INWC_DATARECORD <<" A_B1R14["<<ma_tmp<<"]: "<<A_B1R14[ma_tmp] << ", st3_Tw[14] = " << 1 << endl;    
					        INWC_DATARECORD <<" A_B1R15["<<ma_tmp<<"]: "<<A_B1R15[ma_tmp] << ", st3_Tw[15] = " << 1 << endl;    
                            if(!debug) INWC_DATARECORD << "InvPhi_0t_Order = " 	<< InvPhi_0t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_0t_Order = " 	<< InvPhi_deg*0 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_1t_Order = " 	<< InvPhi_1t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_1t_Order = " 	<< InvPhi_deg*1 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_2t_Order = " 	<< InvPhi_2t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_2t_Order = " 	<< InvPhi_deg*2 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_3t_Order = " 	<< InvPhi_3t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_3t_Order = " 	<< InvPhi_deg*3 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_4t_Order = " 	<< InvPhi_4t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_4t_Order = " 	<< InvPhi_deg*4 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_5t_Order = " 	<< InvPhi_5t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_5t_Order = " 	<< InvPhi_deg*5 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_6t_Order = " 	<< InvPhi_6t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_6t_Order = " 	<< InvPhi_deg*6 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_7t_Order = " 	<< InvPhi_7t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_7t_Order = " 	<< InvPhi_deg*7 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_8t_Order = " 	<< InvPhi_8t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_8t_Order = " 	<< InvPhi_deg*8 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_9t_Order = " 	<< InvPhi_9t_Order 	<< endl; else INWC_DATARECORD << "InvPhi_9t_Order = " 	<< InvPhi_deg*9 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_10t_Order = " << InvPhi_10t_Order << endl; else INWC_DATARECORD << "InvPhi_10t_Order = " 	<< InvPhi_deg*10 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_11t_Order = " << InvPhi_11t_Order << endl; else INWC_DATARECORD << "InvPhi_11t_Order = " 	<< InvPhi_deg*11 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_12t_Order = " << InvPhi_12t_Order << endl; else INWC_DATARECORD << "InvPhi_12t_Order = " 	<< InvPhi_deg*12 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_13t_Order = " << InvPhi_13t_Order << endl; else INWC_DATARECORD << "InvPhi_13t_Order = " 	<< InvPhi_deg*13 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_14t_Order = " << InvPhi_14t_Order << endl; else INWC_DATARECORD << "InvPhi_14t_Order = " 	<< InvPhi_deg*14 	<< endl;
							if(!debug) INWC_DATARECORD << "InvPhi_15t_Order = " << InvPhi_15t_Order << endl; else INWC_DATARECORD << "InvPhi_15t_Order = " 	<< InvPhi_deg*15 	<< endl;
							INWC_Radix16_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp], A_B1R2[ma_tmp], A_B1R3[ma_tmp],A_B1R4[ma_tmp],
					        		   A_B1R5[ma_tmp],A_B1R6[ma_tmp], A_B1R7[ma_tmp], A_B1R8[ma_tmp],A_B1R9[ma_tmp],
					        		   A_B1R10[ma_tmp],A_B1R11[ma_tmp],A_B1R12[ma_tmp],A_B1R13[ma_tmp],A_B1R14[ma_tmp],
					        		   A_B1R15[ma_tmp]);
							//----------------------compute for INWC--------------------------
                            if(!debug) MulMod(InvPhi_0t_dot_IW, InvPhi_0t_Order, 1, p);
							if(!debug) MulMod(InvPhi_1t_dot_IW, InvPhi_1t_Order, 1, p);
							if(!debug) MulMod(InvPhi_2t_dot_IW, InvPhi_2t_Order, 1, p);
							if(!debug) MulMod(InvPhi_3t_dot_IW, InvPhi_3t_Order, 1, p);
							if(!debug) MulMod(InvPhi_4t_dot_IW, InvPhi_4t_Order, 1, p);
							if(!debug) MulMod(InvPhi_5t_dot_IW, InvPhi_5t_Order, 1, p);
							if(!debug) MulMod(InvPhi_6t_dot_IW, InvPhi_6t_Order, 1, p);
							if(!debug) MulMod(InvPhi_7t_dot_IW, InvPhi_7t_Order, 1, p);
							if(!debug) MulMod(InvPhi_8t_dot_IW, InvPhi_8t_Order, 1, p);
							if(!debug) MulMod(InvPhi_9t_dot_IW, InvPhi_9t_Order, 1, p);
							if(!debug) MulMod(InvPhi_10t_dot_IW, InvPhi_10t_Order, 1, p);
							if(!debug) MulMod(InvPhi_11t_dot_IW, InvPhi_11t_Order, 1, p);
							if(!debug) MulMod(InvPhi_12t_dot_IW, InvPhi_12t_Order, 1, p);
							if(!debug) MulMod(InvPhi_13t_dot_IW, InvPhi_13t_Order, 1, p);
							if(!debug) MulMod(InvPhi_14t_dot_IW, InvPhi_14t_Order, 1, p);
							if(!debug) MulMod(InvPhi_15t_dot_IW, InvPhi_15t_Order, 1, p);

							if(!debug) MulMod(InvPhi_0t_dot_IW_dot_InvSixTeen, InvPhi_0t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_1t_dot_IW_dot_InvSixTeen, InvPhi_1t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_2t_dot_IW_dot_InvSixTeen, InvPhi_2t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_3t_dot_IW_dot_InvSixTeen, InvPhi_3t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_4t_dot_IW_dot_InvSixTeen, InvPhi_4t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_5t_dot_IW_dot_InvSixTeen, InvPhi_5t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_6t_dot_IW_dot_InvSixTeen, InvPhi_6t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_7t_dot_IW_dot_InvSixTeen, InvPhi_7t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_8t_dot_IW_dot_InvSixTeen, InvPhi_8t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_9t_dot_IW_dot_InvSixTeen, InvPhi_9t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_10t_dot_IW_dot_InvSixTeen, InvPhi_10t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_11t_dot_IW_dot_InvSixTeen, InvPhi_11t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_12t_dot_IW_dot_InvSixTeen, InvPhi_12t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_13t_dot_IW_dot_InvSixTeen, InvPhi_13t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_14t_dot_IW_dot_InvSixTeen, InvPhi_14t_dot_IW, InvSixTeen, p);
							if(!debug) MulMod(InvPhi_15t_dot_IW_dot_InvSixTeen, InvPhi_15t_dot_IW, InvSixTeen, p);
							//----------------------------------------------------------------
                            MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],InvPhi_0t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],InvPhi_1t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],InvPhi_2t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],InvPhi_3t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R4[ma_tmp],A_B1R4[ma_tmp],InvPhi_4t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R5[ma_tmp],A_B1R5[ma_tmp],InvPhi_5t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R6[ma_tmp],A_B1R6[ma_tmp],InvPhi_6t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R7[ma_tmp],A_B1R7[ma_tmp],InvPhi_7t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R8[ma_tmp],A_B1R8[ma_tmp],InvPhi_8t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R9[ma_tmp],A_B1R9[ma_tmp],InvPhi_9t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R10[ma_tmp],A_B1R10[ma_tmp],InvPhi_10t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R11[ma_tmp],A_B1R11[ma_tmp],InvPhi_11t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R12[ma_tmp],A_B1R12[ma_tmp],InvPhi_12t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R13[ma_tmp],A_B1R13[ma_tmp],InvPhi_13t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R14[ma_tmp],A_B1R14[ma_tmp],InvPhi_14t_dot_IW_dot_InvSixTeen,p);
					        MulMod(A_B1R15[ma_tmp],A_B1R15[ma_tmp],InvPhi_15t_dot_IW_dot_InvSixTeen,p);
                            break;
                    }
					
                    INWC_DATARECORD <<" ------after BU compute and Mul-------" << std::endl;
				    INWC_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<<"\n";
				    INWC_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<<"\n";
					INWC_DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<<A_B1R2[ma_tmp]<<"\n";
					INWC_DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<<A_B1R3[ma_tmp]<<"\n";
					INWC_DATARECORD <<"A_B1R4["<<ma_tmp<<"]: "<<A_B1R4[ma_tmp]<<"\n";
					INWC_DATARECORD <<"A_B1R5["<<ma_tmp<<"]: "<<A_B1R5[ma_tmp]<<"\n";
					INWC_DATARECORD <<"A_B1R6["<<ma_tmp<<"]: "<<A_B1R6[ma_tmp]<<"\n";
					INWC_DATARECORD <<"A_B1R7["<<ma_tmp<<"]: "<<A_B1R7[ma_tmp]<<"\n";
					INWC_DATARECORD <<"A_B1R8["<<ma_tmp<<"]: "<<A_B1R8[ma_tmp]<<"\n";
					INWC_DATARECORD <<"A_B1R9["<<ma_tmp<<"]: "<<A_B1R9[ma_tmp]<<"\n";
					INWC_DATARECORD <<"A_B1R10["<<ma_tmp<<"]: "<<A_B1R10[ma_tmp]<<"\n";
					INWC_DATARECORD <<"A_B1R11["<<ma_tmp<<"]: "<<A_B1R11[ma_tmp]<<"\n";
					INWC_DATARECORD <<"A_B1R12["<<ma_tmp<<"]: "<<A_B1R12[ma_tmp]<<"\n";
					INWC_DATARECORD <<"A_B1R13["<<ma_tmp<<"]: "<<A_B1R13[ma_tmp]<<"\n";
					INWC_DATARECORD <<"A_B1R14["<<ma_tmp<<"]: "<<A_B1R14[ma_tmp]<<"\n";
					INWC_DATARECORD <<"A_B1R15["<<ma_tmp<<"]: "<<A_B1R15[ma_tmp]<<"\n";
					INWC_DATARECORD <<"--------------------------------------------------------------------\n";
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
