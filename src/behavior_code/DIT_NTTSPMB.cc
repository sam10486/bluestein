#include <iostream>
#include <fstream>
#include <NTL/ZZ.h>
#include <time.h>

#include "DIT_NTTSPMB.h"
#include "SPMB.h"

#include <vector>
#include <algorithm>
#include "BitOperate.h"
#include "DTFAG.h"

using namespace std;


void DIT_NTTSPMB::DIT_NTT_radix2(std::vector<ZZ> &A){
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

	std::ofstream DIT_DATARECORD("./DIT_NTT_R2_SPMB.txt");
	std::ofstream DIT_spmb_radix2("./SPMB_tw/DIT_spmb_radix2.txt");
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
	ZZ fft_twiddle = W;
	ZZ fft_prime = p;
	int debug = 0;
	vector<vector<ZZ > > ROM0;
    vector<ZZ > ROM1, ROM2;


	int arr_size = radix_r1 * radix_r1;
    ROM0.resize(radix_r1);
    for(int i=0; i<radix_r1; i++){
        ROM0[i].resize(radix_r1);
    }
    ROM1.resize(arr_size);
    ROM2.resize(arr_size);
	DTFAG.DTFAG_ROM_init(
                        radix_r1, radix_r2, fft_twiddle, fft_prime, debug,
                        ROM0, ROM1, ROM2);
	//----------------------------------------

	
    Stage = (unsigned long)ceil(log2(N));
	BC_WIDTH  = (int)ceil(log2(N/2));
	offset    =  (int)N /radix;
	word_size =  (int)N / (2 * radix);
	group     =  (int)N / (radix * radix);
	bit_array_tmp.resize(BC_WIDTH);

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
				std::cout <<"A_B0R0["<<ma_tmp<<"]"<<A[BC]<<"\n";
				std::cout <<"A_B0R1["<<ma_tmp<<"]"<<A[BC+ offset]<<"\n";
			}else {
				A_B1R0[ma_tmp] = A[BC];
				A_B1R1[ma_tmp] = A[BC + offset];
				std::cout <<"A_B1R0["<<ma_tmp<<"]"<<A[BC]<<"\n";
				std::cout <<"A_B1R1["<<ma_tmp<<"]"<<A[BC+ offset]<<"\n";
			}
		}
	}
	
	ma_tmp = 0;
	bn_tmp = 0;
	BC     = 0;
	int tw_degree = pow(2, Stage-1); // siang
	std::cout << "init load over! \n";
	//need modify to mult by twiddle factor
	for(int s = 0; s < Stage; s++){
		//----------DTFAG generate---------------
		int DTFAG_t = 0;
		int DTFAG_i = 0;
		int DTFAG_j = 0;
		//----------DTFAG fin--------------------
        int s_complement = number_complement(s, Stage);
		if(s == 0){
        }
		else {
			tw_degree = tw_degree / radix;
		}
        
		DIT_DATARECORD <<"Stage: "<< s<<"\n";
		DIT_DATARECORD << "factor "<< factor <<"\n";
		DIT_DATARECORD << "****************************\n";

		DIT_spmb_radix2 <<"Now Stage: "<< s <<"\n";
		DIT_spmb_radix2 <<"twiddle factor : "<< factor <<"\n";
		for(int i = 0 ;i < group;i++){
			DIT_spmb_radix2 << "----------------i =" << i << " ----------------" << std::endl;
            DIT_DATARECORD << "----------------i =" << i << " ----------------" << std::endl;
			for(int j = 0;j < radix;j++){
                DIT_DATARECORD << "---j =" << j << " ---" << std::endl;
				DIT_DATARECORD <<"twiddle factor : "<< factor <<"\n";
				DIT_DATARECORD <<"p : "<< p <<"\n";
	    		DIT_DATARECORD << "********\n";
				gray_i  = Gray(i,group);
				BC_tmp  = j * group + gray_i;
				RR_R2(BC_tmp,s,BC);
                DIT_DATARECORD << "BC_tmp = " << BC_tmp << ", BC = " << BC << ", s_complement = " << s_complement << endl;
				length = BC_tmp >> s_complement;
				DIT_DATARECORD << "length: " <<  length <<"\n";
				PowerMod(factor_t,factor,length,p);
				DIT_DATARECORD << "factor_t: "<<factor_t<<"\n";

				AGU_R2(BC,bn_tmp,ma_tmp);

				//-----------compute data idx-------------
				DIT_spmb_radix2 << "BC_tmp = " << BC_tmp << std::endl;
				
				IntToVec.IntToVec(BC_tmp, N, bit_array);
				rotate(bit_array.begin(), bit_array.begin()+ bit_width_s*s , bit_array.end());
				long long Data = VecToInt.VecToInt(bit_array, N);
				DIT_spmb_radix2 << "Data_index = ";
                DIT_spmb_radix2 << "( " ;					
				for(int k = 0; k < radix ; k++ ){
					DIT_spmb_radix2 << Data + k*(1<<(bit_width-bit_width_s-bit_width_s*s)) <<" ";	
				}
				DIT_spmb_radix2 << ") " ;
				DIT_spmb_radix2 << ", (w^" << 0 << ", w^" << tw_degree * length << ")" <<std::endl;
				//-----------------------------------------

				//-----------DFTAG generator--------------
				DTFAG.DTFAG_SPMB_DIT(
								s, fft_point, radix_r1, radix_r2, debug,
								ROM0, ROM1, ROM2,
								st0_Tw, st1_Tw, st2_Tw, 
								DTFAG_t, DTFAG_i, DTFAG_j);
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
				if(DTFAG_i == radix-1 && DTFAG_j == radix-1 && DTFAG_t == radix-1){
					DTFAG_t = 0;
				}else if(DTFAG_i == radix-1 && DTFAG_j == radix-1){
					DTFAG_t++;
				}
				if(DTFAG_j == radix-1 && DTFAG_i == radix-1){
					DTFAG_i = 0;
				}else if(DTFAG_j == radix-1){
					DTFAG_i++;
				}
				if(DTFAG_j == radix-1){
					DTFAG_j = 0;
				}else{
					DTFAG_j++;
				}
				//-----------DTFAG fin--------------------
				if(bn_tmp == 0){
                    DIT_DATARECORD << "bn_tmp = " << bn_tmp << ", ma_tmp = " << ma_tmp << std::endl;
					bn0_bc_tmp = BC_tmp;
					switch(s){
						case 0:
							DIT_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp]<< ", st0_Tw[0] = " << st0_Tw[0] << "\n";
							DIT_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp]<< ", st0_Tw[1] = " << st0_Tw[1] << "\n";
							if(!debug) Radix2_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp]);
							if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],st0_Tw[0],p);
							if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],st0_Tw[1],p);
							break;
						case 1:
							DIT_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp]<< ", st1_Tw[0] = " << st1_Tw[0] << "\n";
							DIT_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp]<< ", st1_Tw[1] = " << st1_Tw[1] << "\n";
							if(!debug) Radix2_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp]);
							if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],st1_Tw[0],p);
							if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],st1_Tw[1],p);
							break;
						case 2:
							DIT_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp]<< ", st2_Tw[0] = " << st2_Tw[0] << "\n";
							DIT_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp]<< ", st2_Tw[1] = " << st2_Tw[1] << "\n";
							if(!debug) Radix2_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp]);
							if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],st2_Tw[0],p);
							if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],st2_Tw[1],p);
							break;
						case 3:
							DIT_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp]<< ", st3_Tw[0] = " << 1 << "\n";
							DIT_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp]<< ", st3_Tw[1] = " << 1 << "\n";
							if(!debug) Radix2_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp]);
							if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],1,p);
							if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],1,p);
							break;
					}
					
					DIT_DATARECORD << "---after BU compute---" << std::endl;
                    DIT_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp]<<"\n";
					DIT_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp]<<"\n";
				    
					bn0_ma_reg = ma_tmp;

				}
			    else {
                    DIT_DATARECORD << "bn_tmp = " << bn_tmp << ", ma_tmp = " << ma_tmp << std::endl;
					bn1_bc_tmp = BC_tmp;
					switch(s){
						case 0:
							DIT_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<< ", st0_Tw[0] = " << st0_Tw[0] << "\n";
							DIT_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<< ", st0_Tw[1] = " << st0_Tw[1] << "\n";
							if(!debug) Radix2_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp]);
							if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],st0_Tw[0],p);
							if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],st0_Tw[1],p);
							break;
						case 1:
							DIT_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<< ", st1_Tw[0] = " << st1_Tw[0] << "\n";
							DIT_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<< ", st1_Tw[1] = " << st1_Tw[1] << "\n";
							if(!debug) Radix2_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp]);
							if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],st1_Tw[0],p);
							if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],st1_Tw[1],p);
							break;
						case 2:
							DIT_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<< ", st2_Tw[0] = " << st2_Tw[0] << "\n";
							DIT_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<< ", st2_Tw[1] = " << st2_Tw[1] << "\n";
							if(!debug) Radix2_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp]);
							if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],st2_Tw[0],p);
							if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],st2_Tw[1],p);
							break;
						case 3:
							DIT_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<< ", st3_Tw[0] = " << 1 << "\n";
							DIT_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<< ", st3_Tw[1] = " << 1 << "\n";
							if(!debug) Radix2_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp]);
							if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],1,p);
							if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],1,p);
							break;
					}
					
                    DIT_DATARECORD << "---after BU compute---" << std::endl;
                    DIT_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<<"\n";
					DIT_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<<"\n";
				    
					bn1_ma_reg = ma_tmp;
				}
				DIT_DATARECORD <<"--------------------------------------------------------------------\n";
			}
		//data relocation
		if(s < Stage-1){
		    if(bn1_bc_tmp > bn0_bc_tmp){
                DIT_DATARECORD << "bn0_bc_tmp = " << bn0_bc_tmp << ", bn1_bc_tmp = " << bn1_bc_tmp << std::endl;
		        data_tmp = A_B0R1[bn0_ma_reg];
		        A_B0R1[bn0_ma_reg] = A_B1R0[bn1_ma_reg];
		        A_B1R0[bn1_ma_reg] = data_tmp;
		    }else {
                DIT_DATARECORD << "bn0_bc_tmp = " << bn0_bc_tmp << ", bn1_bc_tmp = " << bn1_bc_tmp << std::endl;
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

void DIT_NTTSPMB::DIT_NTT_radix4(std::vector<ZZ> &A){
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
	
	std::ofstream DATARECORD("./DIT_NTT_R4_SPMB.txt");
	std::ofstream DIT_spmb_radix4("./SPMB_tw/DIT_spmb_radix4.txt");
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
	ZZ fft_twiddle = W;
	ZZ fft_prime = p;
	int debug = 0;
	vector<vector<ZZ > > ROM0;
    vector<ZZ > ROM1, ROM2;


	int arr_size = radix_r1 * radix_r1;
    ROM0.resize(radix_r1);
    for(int i=0; i<radix_r1; i++){
        ROM0[i].resize(radix_r1);
    }
    ROM1.resize(arr_size);
    ROM2.resize(arr_size);
	DTFAG.DTFAG_ROM_init(
                        radix_r1, radix_r2, fft_twiddle, fft_prime, debug,
                        ROM0, ROM1, ROM2);
	//----------------------------------------
	
    Stage = (unsigned long)ceil(log2(N));
	Stage = (unsigned long)Stage / 2;
	BC_WIDTH  = (int)ceil(log2(N/4));
	offset    =  (int)N /radix;
	word_size =  (int)N / (2 * radix);
	group     =  (int)N / (radix * radix);
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
		

		DIT_spmb_radix4 <<"Now Stage: "<< s <<"\n";
		DIT_spmb_radix4 <<"twiddle factor : "<< factor <<"\n";
	    DIT_spmb_radix4 << "********\n";
		for(int i = 0 ;i < group;i++){
			DATARECORD <<"twiddle factor : "<< factor <<"\n";
			DATARECORD <<"p : "<< p <<"\n";
	    	DATARECORD << "********\n";
			bn0_bc_tmp  = 0;
			bn1_bc_tmp  = 0;
			DIT_spmb_radix4 <<"--------------i = " << i << "----------------\n";
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
				DIT_spmb_radix4 << "BC_tmp = " << BC_tmp << std::endl;
				
				IntToVec.IntToVec(BC_tmp, N, bit_array);
				rotate(bit_array.begin(), bit_array.begin()+ bit_width_s*s , bit_array.end());
				long long Data = VecToInt.VecToInt(bit_array, N);
				DIT_spmb_radix4 << "Data_index = ";
                DIT_spmb_radix4 << "( " ;		
				int spmb_radix4_arr[radix];		
				for(long long k = 0; k < radix ; k++ ){
					long long BR_tmp = BR.BitReserve(k, log2(radix));
					spmb_radix4_arr[BR_tmp] = Data + k*(1<<(bit_width-bit_width_s-bit_width_s*s));	
				}
				for(int k=0; k<radix; k++) {
					DIT_spmb_radix4 << spmb_radix4_arr[k] << ", ";
				}
				DIT_spmb_radix4 << ") " ;
				DIT_spmb_radix4 << ", (w^" << 0 << ", w^" << tw_degree * length 
				<< ", w^"  << tw_degree * length * 2 << ", w^" << tw_degree * length * 3 << ")" << std::endl;
				//-----------------------------------------
				
				//-----------DFTAG generator--------------
				DTFAG.DTFAG_SPMB_DIT(
								s, fft_point, radix_r1, radix_r2, debug,
								ROM0, ROM1, ROM2,
								st0_Tw, st1_Tw, st2_Tw, 
								DTFAG_t, DTFAG_i, DTFAG_j);
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
				if(DTFAG_i == radix-1 && DTFAG_j == radix-1 && DTFAG_t == radix-1){
					DTFAG_t = 0;
				}else if(DTFAG_i == radix-1 && DTFAG_j == radix-1){
					DTFAG_t++;
				}
				if(DTFAG_j == radix-1 && DTFAG_i == radix-1){
					DTFAG_i = 0;
				}else if(DTFAG_j == radix-1){
					DTFAG_i++;
				}
				if(DTFAG_j == radix-1){
					DTFAG_j = 0;
				}else{
					DTFAG_j++;
				}
				//-----------DTFAG fin--------------------
				PowerMod(factor_t,factor,length,p);
				DATARECORD << "factor_t: "<<factor_t<<"\n";
				AGU_R4(BC,bn_tmp,ma_tmp);
				DATARECORD << "bn_tmp: "<<bn_tmp<<"\n";
				if(bn_tmp == 0){
					if(j < 2)bn0_bc_tmp = BC_tmp;
					switch(s){
						case 0:
							DATARECORD <<" A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp] << ", st0_Tw[0] = " << st0_Tw[0] << endl;
				    		DATARECORD <<" A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp] << ", st0_Tw[1] = " << st0_Tw[1] << endl;
							DATARECORD <<" A_B0R2["<<ma_tmp<<"]: "<<A_B0R2[ma_tmp] << ", st0_Tw[2] = " << st0_Tw[2] << endl;
							DATARECORD <<" A_B0R3["<<ma_tmp<<"]: "<<A_B0R3[ma_tmp] << ", st0_Tw[3] = " << st0_Tw[3] << endl;
							if(!debug) Radix4_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp],A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
							if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],st0_Tw[0],p);
							if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],st0_Tw[1],p);
							if(!debug) MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],st0_Tw[2],p);
							if(!debug) MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],st0_Tw[3],p);
							break;
						case 1:
							DATARECORD <<" A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp] << ", st1_Tw[0] = " << st1_Tw[0] << endl;
							DATARECORD <<" A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp] << ", st1_Tw[1] = " << st1_Tw[1] << endl;
							DATARECORD <<" A_B0R2["<<ma_tmp<<"]: "<<A_B0R2[ma_tmp] << ", st1_Tw[2] = " << st1_Tw[2] << endl;
							DATARECORD <<" A_B0R3["<<ma_tmp<<"]: "<<A_B0R3[ma_tmp] << ", st1_Tw[3] = " << st1_Tw[3] << endl;
							if(!debug) Radix4_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp],A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
							if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],st1_Tw[0],p);
							if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],st1_Tw[1],p);
							if(!debug) MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],st1_Tw[2],p);
							if(!debug) MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],st1_Tw[3],p);
							break;
						case 2:
							DATARECORD <<" A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp] << ", st2_Tw[0] = " << st2_Tw[0] << endl;
							DATARECORD <<" A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp] << ", st2_Tw[1] = " << st2_Tw[1] << endl;
							DATARECORD <<" A_B0R2["<<ma_tmp<<"]: "<<A_B0R2[ma_tmp] << ", st2_Tw[2] = " << st2_Tw[2] << endl;
							DATARECORD <<" A_B0R3["<<ma_tmp<<"]: "<<A_B0R3[ma_tmp] << ", st2_Tw[3] = " << st2_Tw[3] << endl;
							if(!debug) Radix4_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp],A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
							if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],st2_Tw[0],p);
							if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],st2_Tw[1],p);
							if(!debug) MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],st2_Tw[2],p);
							if(!debug) MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],st2_Tw[3],p);
							break;
						case 3:
							DATARECORD <<" A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp] << ", st3_Tw[0] = " << 1 << endl;
							DATARECORD <<" A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp] << ", st3_Tw[1] = " << 1 << endl;
							DATARECORD <<" A_B0R2["<<ma_tmp<<"]: "<<A_B0R2[ma_tmp] << ", st3_Tw[2] = " << 1 << endl;
							DATARECORD <<" A_B0R3["<<ma_tmp<<"]: "<<A_B0R3[ma_tmp] << ", st3_Tw[3] = " << 1 << endl;
							if(!debug) Radix4_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp],A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
							if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],1,p);
							if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],1,p);
							if(!debug) MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],1,p);
							if(!debug) MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],1,p);
							break;
					}
					DATARECORD <<" ------after BU compute and Mul-------" << std::endl;
					DATARECORD <<" A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp] << endl;
					DATARECORD <<" A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp] << endl;
					DATARECORD <<" A_B0R2["<<ma_tmp<<"]: "<<A_B0R2[ma_tmp] << endl;
					DATARECORD <<" A_B0R3["<<ma_tmp<<"]: "<<A_B0R3[ma_tmp] << endl;
					
					if(j <  2)bn0_ma_reg1 = ma_tmp;
					if(j >= 2)bn0_ma_reg2 = ma_tmp;
				}
			    else {
					if(j < 2)bn1_bc_tmp = BC_tmp;
					switch(s){
						case 0:
							DATARECORD <<" A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp] << ", st0_Tw[0] = " << st0_Tw[0] << endl;
							DATARECORD <<" A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp] << ", st0_Tw[1] = " << st0_Tw[1] << endl;
							DATARECORD <<" A_B1R2["<<ma_tmp<<"]: "<<A_B1R2[ma_tmp] << ", st0_Tw[2] = " << st0_Tw[2] << endl;
							DATARECORD <<" A_B1R3["<<ma_tmp<<"]: "<<A_B1R3[ma_tmp] << ", st0_Tw[3] = " << st0_Tw[3] << endl;
							if(!debug) Radix4_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp],A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
							if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],st0_Tw[0],p);
							if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],st0_Tw[1],p);
							if(!debug) MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],st0_Tw[2],p);
							if(!debug) MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],st0_Tw[3],p);
							break;
						case 1:
							DATARECORD <<" A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp] << ", st1_Tw[0] = " << st1_Tw[0] << endl;
							DATARECORD <<" A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp] << ", st1_Tw[1] = " << st1_Tw[1] << endl;
							DATARECORD <<" A_B1R2["<<ma_tmp<<"]: "<<A_B1R2[ma_tmp] << ", st1_Tw[2] = " << st1_Tw[2] << endl;
							DATARECORD <<" A_B1R3["<<ma_tmp<<"]: "<<A_B1R3[ma_tmp] << ", st1_Tw[3] = " << st1_Tw[3] << endl;
							if(!debug) Radix4_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp],A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
							if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],st1_Tw[0],p);
							if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],st1_Tw[1],p);
							if(!debug) MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],st1_Tw[2],p);
							if(!debug) MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],st1_Tw[3],p);
							break;
						case 2:
							DATARECORD <<" A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp] << ", st2_Tw[0] = " << st2_Tw[0] << endl;
							DATARECORD <<" A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp] << ", st2_Tw[1] = " << st2_Tw[1] << endl;
							DATARECORD <<" A_B1R2["<<ma_tmp<<"]: "<<A_B1R2[ma_tmp] << ", st2_Tw[2] = " << st2_Tw[2] << endl;
							DATARECORD <<" A_B1R3["<<ma_tmp<<"]: "<<A_B1R3[ma_tmp] << ", st2_Tw[3] = " << st2_Tw[3] << endl;
							if(!debug) Radix4_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp],A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
							if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],st2_Tw[0],p);
							if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],st2_Tw[1],p);
							if(!debug) MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],st2_Tw[2],p);
							if(!debug) MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],st2_Tw[3],p);
							break;
						case 3:
							DATARECORD <<" A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp] << ", st3_Tw[0] = " << 1 << endl;
							DATARECORD <<" A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp] << ", st3_Tw[1] = " << 1 << endl;
							DATARECORD <<" A_B1R2["<<ma_tmp<<"]: "<<A_B1R2[ma_tmp] << ", st3_Tw[2] = " << 1 << endl;
							DATARECORD <<" A_B1R3["<<ma_tmp<<"]: "<<A_B1R3[ma_tmp] << ", st3_Tw[3] = " << 1 << endl;
							if(!debug) Radix4_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp],A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
							if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],1,p);
							if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],1,p);
							if(!debug) MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],1,p);
							if(!debug) MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],1,p);
							break;
					}
			
					DATARECORD <<" ------after BU compute and Mul-------" << std::endl;
					DATARECORD <<" A_B1R0["<<ma_tmp<<"]: " << A_B1R0[ma_tmp] << endl;
					DATARECORD <<" A_B1R1["<<ma_tmp<<"]: " << A_B1R1[ma_tmp] << endl;
					DATARECORD <<" A_B1R2["<<ma_tmp<<"]: " << A_B1R2[ma_tmp] << endl;
					DATARECORD <<" A_B1R3["<<ma_tmp<<"]: " << A_B1R3[ma_tmp] << endl;
					
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

void DIT_NTTSPMB::DIT_NTT_radix16(std::vector<ZZ> &A){
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

    Stage = (unsigned long)ceil(log2(N));
	Stage = (unsigned long)Stage / 4;
	BC_WIDTH  = (int)ceil(log2(N/16));
	offset    =  (int)N /radix;
	word_size =  (int)N / (2 * radix);
	group     =  (int)N / (radix * radix);
	bit_array_tmp.resize(BC_WIDTH);
	DATARECORD <<"Stage: "<< Stage<<"\n";

	ZZ  factor;   //base factor
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
	ZZ fft_twiddle = W;
	ZZ fft_prime = p;
	int debug = 0;
	vector<vector<ZZ > > ROM0;
    vector<ZZ > ROM1, ROM2;


	int arr_size = radix_r1 * radix_r1;
    ROM0.resize(radix_r1);
    for(int i=0; i<radix_r1; i++){
        ROM0[i].resize(radix_r1);
    }
    ROM1.resize(arr_size);
    ROM2.resize(arr_size);
	DTFAG.DTFAG_ROM_init(
                        radix_r1, radix_r2, fft_twiddle, fft_prime, debug,
                        ROM0, ROM1, ROM2);
	//----------------------------------------
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
		std::cout << "twiddle factor : "<< factor <<"\n";
	    DATARECORD << "********\n";
		for(int i = 0 ;i < group;i++){
			bn0_bc_tmp  = 0;
			bn1_bc_tmp  = 0;
			for(int j = 0;j < radix;j++){
				gray_i  = Gray(i,group);
				BC_tmp  = j * group + gray_i;
				DATARECORD << "i: " << i <<"\n";
				DATARECORD << "gray_i: " << gray_i <<"\n";
				DATARECORD << "BC_tmp: " << BC_tmp <<"\n";
				RR_R16(BC_tmp,s,BC);
				DATARECORD << "BC: " << BC <<"\n";
				length = BC_tmp >> (4*s);
				DATARECORD << "length: " <<  length <<"\n";
				
				AGU_R16(BC,bn_tmp,ma_tmp);
				DATARECORD << "bn_tmp: "<<bn_tmp<<"\n";

				//-----------DFTAG generator--------------
				DTFAG.DTFAG_SPMB_DIT(
								s, fft_point, radix_r1, radix_r2, debug,
								ROM0, ROM1, ROM2,
								st0_Tw, st1_Tw, st2_Tw, 
								DTFAG_t, DTFAG_i, DTFAG_j);
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
				if(DTFAG_i == radix-1 && DTFAG_j == radix-1 && DTFAG_t == radix-1){
					DTFAG_t = 0;
				}else if(DTFAG_i == radix-1 && DTFAG_j == radix-1){
					DTFAG_t++;
				}
				if(DTFAG_j == radix-1 && DTFAG_i == radix-1){
					DTFAG_i = 0;
				}else if(DTFAG_j == radix-1){
					DTFAG_i++;
				}
				if(DTFAG_j == radix-1){
					DTFAG_j = 0;
				}else{
					DTFAG_j++;
				}
				//-----------DTFAG fin--------------------
				if(bn_tmp == 0){
					if(j < 2)bn0_bc_tmp = BC_tmp;
					switch(s){
						case 0:
							DATARECORD <<" A_B0R0["<<ma_tmp<<"]:  "<<A_B0R0[ma_tmp]	<< ", st0_Tw[0] = "  << st0_Tw[0] << endl;
				    		DATARECORD <<" A_B0R1["<<ma_tmp<<"]:  "<<A_B0R1[ma_tmp]	<< ", st0_Tw[1] = "  << st0_Tw[1] << endl;
							DATARECORD <<" A_B0R2["<<ma_tmp<<"]:  "<<A_B0R2[ma_tmp]	<< ", st0_Tw[2] = "  << st0_Tw[2] << endl;
							DATARECORD <<" A_B0R3["<<ma_tmp<<"]:  "<<A_B0R3[ma_tmp]	<< ", st0_Tw[3] = "  << st0_Tw[3] << endl;
							DATARECORD <<" A_B0R4["<<ma_tmp<<"]:  "<<A_B0R4[ma_tmp]	<< ", st0_Tw[4] = "  << st0_Tw[4] << endl;
							DATARECORD <<" A_B0R5["<<ma_tmp<<"]:  "<<A_B0R5[ma_tmp]	<< ", st0_Tw[5] = "  << st0_Tw[5] << endl;
							DATARECORD <<" A_B0R6["<<ma_tmp<<"]:  "<<A_B0R6[ma_tmp]	<< ", st0_Tw[6] = "  << st0_Tw[6] << endl;
							DATARECORD <<" A_B0R7["<<ma_tmp<<"]:  "<<A_B0R7[ma_tmp]	<< ", st0_Tw[7] = "  << st0_Tw[7] << endl;
							DATARECORD <<" A_B0R8["<<ma_tmp<<"]:  "<<A_B0R8[ma_tmp]	<< ", st0_Tw[8] = "  << st0_Tw[8] << endl;
							DATARECORD <<" A_B0R9["<<ma_tmp<<"]:  "<<A_B0R9[ma_tmp]	<< ", st0_Tw[9] = "  << st0_Tw[9] << endl;
							DATARECORD <<" A_B0R10["<<ma_tmp<<"]: "<<A_B0R10[ma_tmp]<< ", st0_Tw[10] = " << st0_Tw[10] << endl;
							DATARECORD <<" A_B0R11["<<ma_tmp<<"]: "<<A_B0R11[ma_tmp]<< ", st0_Tw[11] = " << st0_Tw[11] << endl;
							DATARECORD <<" A_B0R12["<<ma_tmp<<"]: "<<A_B0R12[ma_tmp]<< ", st0_Tw[12] = " << st0_Tw[12] << endl;
							DATARECORD <<" A_B0R13["<<ma_tmp<<"]: "<<A_B0R13[ma_tmp]<< ", st0_Tw[13] = " << st0_Tw[13] << endl;
							DATARECORD <<" A_B0R14["<<ma_tmp<<"]: "<<A_B0R14[ma_tmp]<< ", st0_Tw[14] = " << st0_Tw[14] << endl;
							DATARECORD <<" A_B0R15["<<ma_tmp<<"]: "<<A_B0R15[ma_tmp]<< ", st0_Tw[15] = " << st0_Tw[15] << endl;
							Radix16_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp], A_B0R2[ma_tmp], A_B0R3[ma_tmp],A_B0R4[ma_tmp],
							   A_B0R5[ma_tmp],A_B0R6[ma_tmp], A_B0R7[ma_tmp], A_B0R8[ma_tmp],A_B0R9[ma_tmp],
							   A_B0R10[ma_tmp],A_B0R11[ma_tmp],A_B0R12[ma_tmp],A_B0R13[ma_tmp],A_B0R14[ma_tmp],
							   A_B0R15[ma_tmp]);
							MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],st0_Tw[0],p);   
							MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],st0_Tw[1],p);
							MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],st0_Tw[2],p);
							MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],st0_Tw[3],p);
							MulMod(A_B0R4[ma_tmp],A_B0R4[ma_tmp],st0_Tw[4],p);
							MulMod(A_B0R5[ma_tmp],A_B0R5[ma_tmp],st0_Tw[5],p);
							MulMod(A_B0R6[ma_tmp],A_B0R6[ma_tmp],st0_Tw[6],p);
							MulMod(A_B0R7[ma_tmp],A_B0R7[ma_tmp],st0_Tw[7],p);
							MulMod(A_B0R8[ma_tmp],A_B0R8[ma_tmp],st0_Tw[8],p);
							MulMod(A_B0R9[ma_tmp],A_B0R9[ma_tmp],st0_Tw[9],p);
							MulMod(A_B0R10[ma_tmp],A_B0R10[ma_tmp],st0_Tw[10],p);
							MulMod(A_B0R11[ma_tmp],A_B0R11[ma_tmp],st0_Tw[11],p);
							MulMod(A_B0R12[ma_tmp],A_B0R12[ma_tmp],st0_Tw[12],p);
							MulMod(A_B0R13[ma_tmp],A_B0R13[ma_tmp],st0_Tw[13],p);
							MulMod(A_B0R14[ma_tmp],A_B0R14[ma_tmp],st0_Tw[14],p);
							MulMod(A_B0R15[ma_tmp],A_B0R15[ma_tmp],st0_Tw[15],p);
							break;
						case 1:
							DATARECORD <<" A_B0R0["<<ma_tmp<<"]:  "<<A_B0R0[ma_tmp]	 << ", st1_Tw[0] = "  << st1_Tw[0] << endl;
				    		DATARECORD <<" A_B0R1["<<ma_tmp<<"]:  "<<A_B0R1[ma_tmp]	 << ", st1_Tw[1] = "  << st1_Tw[1] << endl;
							DATARECORD <<" A_B0R2["<<ma_tmp<<"]:  "<<A_B0R2[ma_tmp]	 << ", st1_Tw[2] = "  << st1_Tw[2] << endl;
							DATARECORD <<" A_B0R3["<<ma_tmp<<"]:  "<<A_B0R3[ma_tmp]	 << ", st1_Tw[3] = "  << st1_Tw[3] << endl;
							DATARECORD <<" A_B0R4["<<ma_tmp<<"]:  "<<A_B0R4[ma_tmp]	 << ", st1_Tw[4] = "  << st1_Tw[4] << endl;
							DATARECORD <<" A_B0R5["<<ma_tmp<<"]:  "<<A_B0R5[ma_tmp]	 << ", st1_Tw[5] = "  << st1_Tw[5] << endl;
							DATARECORD <<" A_B0R6["<<ma_tmp<<"]:  "<<A_B0R6[ma_tmp]	 << ", st1_Tw[6] = "  << st1_Tw[6] << endl;
							DATARECORD <<" A_B0R7["<<ma_tmp<<"]:  "<<A_B0R7[ma_tmp]	 << ", st1_Tw[7] = "  << st1_Tw[7] << endl;
							DATARECORD <<" A_B0R8["<<ma_tmp<<"]:  "<<A_B0R8[ma_tmp]	 << ", st1_Tw[8] = "  << st1_Tw[8] << endl;
							DATARECORD <<" A_B0R9["<<ma_tmp<<"]:  "<<A_B0R9[ma_tmp]	 << ", st1_Tw[9] = "  << st1_Tw[9] << endl;
							DATARECORD <<" A_B0R10["<<ma_tmp<<"]: "<<A_B0R10[ma_tmp] << ", st1_Tw[10] = " << st1_Tw[10] << endl;
							DATARECORD <<" A_B0R11["<<ma_tmp<<"]: "<<A_B0R11[ma_tmp] << ", st1_Tw[11] = " << st1_Tw[11] << endl;
							DATARECORD <<" A_B0R12["<<ma_tmp<<"]: "<<A_B0R12[ma_tmp] << ", st1_Tw[12] = " << st1_Tw[12] << endl;
							DATARECORD <<" A_B0R13["<<ma_tmp<<"]: "<<A_B0R13[ma_tmp] << ", st1_Tw[13] = " << st1_Tw[13] << endl;
							DATARECORD <<" A_B0R14["<<ma_tmp<<"]: "<<A_B0R14[ma_tmp] << ", st1_Tw[14] = " << st1_Tw[14] << endl;
							DATARECORD <<" A_B0R15["<<ma_tmp<<"]: "<<A_B0R15[ma_tmp] << ", st1_Tw[15] = " << st1_Tw[15] << endl;
							Radix16_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp], A_B0R2[ma_tmp], A_B0R3[ma_tmp],A_B0R4[ma_tmp],
							   A_B0R5[ma_tmp],A_B0R6[ma_tmp], A_B0R7[ma_tmp], A_B0R8[ma_tmp],A_B0R9[ma_tmp],
							   A_B0R10[ma_tmp],A_B0R11[ma_tmp],A_B0R12[ma_tmp],A_B0R13[ma_tmp],A_B0R14[ma_tmp],
							   A_B0R15[ma_tmp]);
							MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],st1_Tw[0],p);   
							MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],st1_Tw[1],p);
							MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],st1_Tw[2],p);
							MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],st1_Tw[3],p);
							MulMod(A_B0R4[ma_tmp],A_B0R4[ma_tmp],st1_Tw[4],p);
							MulMod(A_B0R5[ma_tmp],A_B0R5[ma_tmp],st1_Tw[5],p);
							MulMod(A_B0R6[ma_tmp],A_B0R6[ma_tmp],st1_Tw[6],p);
							MulMod(A_B0R7[ma_tmp],A_B0R7[ma_tmp],st1_Tw[7],p);
							MulMod(A_B0R8[ma_tmp],A_B0R8[ma_tmp],st1_Tw[8],p);
							MulMod(A_B0R9[ma_tmp],A_B0R9[ma_tmp],st1_Tw[9],p);
							MulMod(A_B0R10[ma_tmp],A_B0R10[ma_tmp],st1_Tw[10],p);
							MulMod(A_B0R11[ma_tmp],A_B0R11[ma_tmp],st1_Tw[11],p);
							MulMod(A_B0R12[ma_tmp],A_B0R12[ma_tmp],st1_Tw[12],p);
							MulMod(A_B0R13[ma_tmp],A_B0R13[ma_tmp],st1_Tw[13],p);
							MulMod(A_B0R14[ma_tmp],A_B0R14[ma_tmp],st1_Tw[14],p);
							MulMod(A_B0R15[ma_tmp],A_B0R15[ma_tmp],st1_Tw[15],p);
							break;
						case 2:
							DATARECORD <<" A_B0R0["<<ma_tmp<<"]:  "<<A_B0R0[ma_tmp]  << ", st2_Tw[0] = "  << st2_Tw[0] << endl;
				    		DATARECORD <<" A_B0R1["<<ma_tmp<<"]:  "<<A_B0R1[ma_tmp]  << ", st2_Tw[1] = "  << st2_Tw[1] << endl;
							DATARECORD <<" A_B0R2["<<ma_tmp<<"]:  "<<A_B0R2[ma_tmp]  << ", st2_Tw[2] = "  << st2_Tw[2] << endl;
							DATARECORD <<" A_B0R3["<<ma_tmp<<"]:  "<<A_B0R3[ma_tmp]  << ", st2_Tw[3] = "  << st2_Tw[3] << endl;
							DATARECORD <<" A_B0R4["<<ma_tmp<<"]:  "<<A_B0R4[ma_tmp]  << ", st2_Tw[4] = "  << st2_Tw[4] << endl;
							DATARECORD <<" A_B0R5["<<ma_tmp<<"]:  "<<A_B0R5[ma_tmp]  << ", st2_Tw[5] = "  << st2_Tw[5] << endl;
							DATARECORD <<" A_B0R6["<<ma_tmp<<"]:  "<<A_B0R6[ma_tmp]  << ", st2_Tw[6] = "  << st2_Tw[6] << endl;
							DATARECORD <<" A_B0R7["<<ma_tmp<<"]:  "<<A_B0R7[ma_tmp]  << ", st2_Tw[7] = "  << st2_Tw[7] << endl;
							DATARECORD <<" A_B0R8["<<ma_tmp<<"]:  "<<A_B0R8[ma_tmp]  << ", st2_Tw[8] = "  << st2_Tw[8] << endl;
							DATARECORD <<" A_B0R9["<<ma_tmp<<"]:  "<<A_B0R9[ma_tmp]  << ", st2_Tw[9] = "  << st2_Tw[9] << endl;
							DATARECORD <<" A_B0R10["<<ma_tmp<<"]: "<<A_B0R10[ma_tmp] << ", st2_Tw[10] = " << st2_Tw[10] << endl;
							DATARECORD <<" A_B0R11["<<ma_tmp<<"]: "<<A_B0R11[ma_tmp] << ", st2_Tw[11] = " << st2_Tw[11] << endl;
							DATARECORD <<" A_B0R12["<<ma_tmp<<"]: "<<A_B0R12[ma_tmp] << ", st2_Tw[12] = " << st2_Tw[12] << endl;
							DATARECORD <<" A_B0R13["<<ma_tmp<<"]: "<<A_B0R13[ma_tmp] << ", st2_Tw[13] = " << st2_Tw[13] << endl;
							DATARECORD <<" A_B0R14["<<ma_tmp<<"]: "<<A_B0R14[ma_tmp] << ", st2_Tw[14] = " << st2_Tw[14] << endl;
							DATARECORD <<" A_B0R15["<<ma_tmp<<"]: "<<A_B0R15[ma_tmp] << ", st2_Tw[15] = " << st2_Tw[15] << endl;
							Radix16_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp], A_B0R2[ma_tmp], A_B0R3[ma_tmp],A_B0R4[ma_tmp],
							   A_B0R5[ma_tmp],A_B0R6[ma_tmp], A_B0R7[ma_tmp], A_B0R8[ma_tmp],A_B0R9[ma_tmp],
							   A_B0R10[ma_tmp],A_B0R11[ma_tmp],A_B0R12[ma_tmp],A_B0R13[ma_tmp],A_B0R14[ma_tmp],
							   A_B0R15[ma_tmp]);
							MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],st2_Tw[0],p);   
							MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],st2_Tw[1],p);
							MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],st2_Tw[2],p);
							MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],st2_Tw[3],p);
							MulMod(A_B0R4[ma_tmp],A_B0R4[ma_tmp],st2_Tw[4],p);
							MulMod(A_B0R5[ma_tmp],A_B0R5[ma_tmp],st2_Tw[5],p);
							MulMod(A_B0R6[ma_tmp],A_B0R6[ma_tmp],st2_Tw[6],p);
							MulMod(A_B0R7[ma_tmp],A_B0R7[ma_tmp],st2_Tw[7],p);
							MulMod(A_B0R8[ma_tmp],A_B0R8[ma_tmp],st2_Tw[8],p);
							MulMod(A_B0R9[ma_tmp],A_B0R9[ma_tmp],st2_Tw[9],p);
							MulMod(A_B0R10[ma_tmp],A_B0R10[ma_tmp],st2_Tw[10],p);
							MulMod(A_B0R11[ma_tmp],A_B0R11[ma_tmp],st2_Tw[11],p);
							MulMod(A_B0R12[ma_tmp],A_B0R12[ma_tmp],st2_Tw[12],p);
							MulMod(A_B0R13[ma_tmp],A_B0R13[ma_tmp],st2_Tw[13],p);
							MulMod(A_B0R14[ma_tmp],A_B0R14[ma_tmp],st2_Tw[14],p);
							MulMod(A_B0R15[ma_tmp],A_B0R15[ma_tmp],st2_Tw[15],p);
							break;
						case 3:
							DATARECORD <<" A_B0R0["<<ma_tmp<<"]:  "<<A_B0R0[ma_tmp]  << ", st3_Tw[0] = "  << 1 << endl;
				    		DATARECORD <<" A_B0R1["<<ma_tmp<<"]:  "<<A_B0R1[ma_tmp]  << ", st3_Tw[1] = "  << 1 << endl;
							DATARECORD <<" A_B0R2["<<ma_tmp<<"]:  "<<A_B0R2[ma_tmp]  << ", st3_Tw[2] = "  << 1 << endl;
							DATARECORD <<" A_B0R3["<<ma_tmp<<"]:  "<<A_B0R3[ma_tmp]  << ", st3_Tw[3] = "  << 1 << endl;
							DATARECORD <<" A_B0R4["<<ma_tmp<<"]:  "<<A_B0R4[ma_tmp]  << ", st3_Tw[4] = "  << 1 << endl;
							DATARECORD <<" A_B0R5["<<ma_tmp<<"]:  "<<A_B0R5[ma_tmp]  << ", st3_Tw[5] = "  << 1 << endl;
							DATARECORD <<" A_B0R6["<<ma_tmp<<"]:  "<<A_B0R6[ma_tmp]  << ", st3_Tw[6] = "  << 1 << endl;
							DATARECORD <<" A_B0R7["<<ma_tmp<<"]:  "<<A_B0R7[ma_tmp]  << ", st3_Tw[7] = "  << 1 << endl;
							DATARECORD <<" A_B0R8["<<ma_tmp<<"]:  "<<A_B0R8[ma_tmp]  << ", st3_Tw[8] = "  << 1 << endl;
							DATARECORD <<" A_B0R9["<<ma_tmp<<"]:  "<<A_B0R9[ma_tmp]  << ", st3_Tw[9] = "  << 1 << endl;
							DATARECORD <<" A_B0R10["<<ma_tmp<<"]: "<<A_B0R10[ma_tmp] << ", st3_Tw[10] = " << 1 << endl;
							DATARECORD <<" A_B0R11["<<ma_tmp<<"]: "<<A_B0R11[ma_tmp] << ", st3_Tw[11] = " << 1 << endl;
							DATARECORD <<" A_B0R12["<<ma_tmp<<"]: "<<A_B0R12[ma_tmp] << ", st3_Tw[12] = " << 1 << endl;
							DATARECORD <<" A_B0R13["<<ma_tmp<<"]: "<<A_B0R13[ma_tmp] << ", st3_Tw[13] = " << 1 << endl;
							DATARECORD <<" A_B0R14["<<ma_tmp<<"]: "<<A_B0R14[ma_tmp] << ", st3_Tw[14] = " << 1 << endl;
							DATARECORD <<" A_B0R15["<<ma_tmp<<"]: "<<A_B0R15[ma_tmp] << ", st3_Tw[15] = " << 1 << endl;
							Radix16_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp], A_B0R2[ma_tmp], A_B0R3[ma_tmp],A_B0R4[ma_tmp],
							   A_B0R5[ma_tmp],A_B0R6[ma_tmp], A_B0R7[ma_tmp], A_B0R8[ma_tmp],A_B0R9[ma_tmp],
							   A_B0R10[ma_tmp],A_B0R11[ma_tmp],A_B0R12[ma_tmp],A_B0R13[ma_tmp],A_B0R14[ma_tmp],
							   A_B0R15[ma_tmp]);
							MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],1,p);   
							MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],1,p);
							MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],1,p);
							MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],1,p);
							MulMod(A_B0R4[ma_tmp],A_B0R4[ma_tmp],1,p);
							MulMod(A_B0R5[ma_tmp],A_B0R5[ma_tmp],1,p);
							MulMod(A_B0R6[ma_tmp],A_B0R6[ma_tmp],1,p);
							MulMod(A_B0R7[ma_tmp],A_B0R7[ma_tmp],1,p);
							MulMod(A_B0R8[ma_tmp],A_B0R8[ma_tmp],1,p);
							MulMod(A_B0R9[ma_tmp],A_B0R9[ma_tmp],1,p);
							MulMod(A_B0R10[ma_tmp],A_B0R10[ma_tmp],1,p);
							MulMod(A_B0R11[ma_tmp],A_B0R11[ma_tmp],1,p);
							MulMod(A_B0R12[ma_tmp],A_B0R12[ma_tmp],1,p);
							MulMod(A_B0R13[ma_tmp],A_B0R13[ma_tmp],1,p);
							MulMod(A_B0R14[ma_tmp],A_B0R14[ma_tmp],1,p);
							MulMod(A_B0R15[ma_tmp],A_B0R15[ma_tmp],1,p);
							break;
					}
					DATARECORD <<" ------after BU compute and Mul-------" << std::endl;
				    DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp]<<"\n";
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
					DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<<A_B0R15[ma_tmp]<<"\n";
					DATARECORD <<"--------------------------------------------------------------------\n";
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
							DATARECORD <<" A_B1R0["<<ma_tmp<<"]:  "<<A_B1R0[ma_tmp]  << ", st0_Tw[0] = "  << st0_Tw[0]  << endl;
				    		DATARECORD <<" A_B1R1["<<ma_tmp<<"]:  "<<A_B1R1[ma_tmp]  << ", st0_Tw[1] = "  << st0_Tw[1]  << endl;
							DATARECORD <<" A_B1R2["<<ma_tmp<<"]:  "<<A_B1R2[ma_tmp]  << ", st0_Tw[2] = "  << st0_Tw[2]  << endl;
							DATARECORD <<" A_B1R3["<<ma_tmp<<"]:  "<<A_B1R3[ma_tmp]  << ", st0_Tw[3] = "  << st0_Tw[3]  << endl;
							DATARECORD <<" A_B1R4["<<ma_tmp<<"]:  "<<A_B1R4[ma_tmp]  << ", st0_Tw[4] = "  << st0_Tw[4]  << endl;
							DATARECORD <<" A_B1R5["<<ma_tmp<<"]:  "<<A_B1R5[ma_tmp]  << ", st0_Tw[5] = "  << st0_Tw[5]  << endl;
							DATARECORD <<" A_B1R6["<<ma_tmp<<"]:  "<<A_B1R6[ma_tmp]  << ", st0_Tw[6] = "  << st0_Tw[6]  << endl;
							DATARECORD <<" A_B1R7["<<ma_tmp<<"]:  "<<A_B1R7[ma_tmp]  << ", st0_Tw[7] = "  << st0_Tw[7]  << endl;
							DATARECORD <<" A_B1R8["<<ma_tmp<<"]:  "<<A_B1R8[ma_tmp]  << ", st0_Tw[8] = "  << st0_Tw[8]  << endl;
							DATARECORD <<" A_B1R9["<<ma_tmp<<"]:  "<<A_B1R9[ma_tmp]  << ", st0_Tw[9] = "  << st0_Tw[9]  << endl;
							DATARECORD <<" A_B1R10["<<ma_tmp<<"]: "<<A_B1R10[ma_tmp] << ", st0_Tw[10] = " << st0_Tw[10] << endl;
							DATARECORD <<" A_B1R11["<<ma_tmp<<"]: "<<A_B1R11[ma_tmp] << ", st0_Tw[11] = " << st0_Tw[11] << endl;
							DATARECORD <<" A_B1R12["<<ma_tmp<<"]: "<<A_B1R12[ma_tmp] << ", st0_Tw[12] = " << st0_Tw[12] << endl;
							DATARECORD <<" A_B1R13["<<ma_tmp<<"]: "<<A_B1R13[ma_tmp] << ", st0_Tw[13] = " << st0_Tw[13] << endl;
							DATARECORD <<" A_B1R14["<<ma_tmp<<"]: "<<A_B1R14[ma_tmp] << ", st0_Tw[14] = " << st0_Tw[14] << endl;
							DATARECORD <<" A_B1R15["<<ma_tmp<<"]: "<<A_B1R15[ma_tmp] << ", st0_Tw[15] = " << st0_Tw[15] << endl;
							Radix16_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp], A_B1R2[ma_tmp], A_B1R3[ma_tmp],A_B1R4[ma_tmp],
							   A_B1R5[ma_tmp],A_B1R6[ma_tmp], A_B1R7[ma_tmp], A_B1R8[ma_tmp],A_B1R9[ma_tmp],
							   A_B1R10[ma_tmp],A_B1R11[ma_tmp],A_B1R12[ma_tmp],A_B1R13[ma_tmp],A_B1R14[ma_tmp],
							   A_B1R15[ma_tmp]);
							MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],st0_Tw[0],p);
							MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],st0_Tw[1],p);
							MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],st0_Tw[2],p);
							MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],st0_Tw[3],p);
							MulMod(A_B1R4[ma_tmp],A_B1R4[ma_tmp],st0_Tw[4],p);
							MulMod(A_B1R5[ma_tmp],A_B1R5[ma_tmp],st0_Tw[5],p);
							MulMod(A_B1R6[ma_tmp],A_B1R6[ma_tmp],st0_Tw[6],p);
							MulMod(A_B1R7[ma_tmp],A_B1R7[ma_tmp],st0_Tw[7],p);
							MulMod(A_B1R8[ma_tmp],A_B1R8[ma_tmp],st0_Tw[8],p);
							MulMod(A_B1R9[ma_tmp],A_B1R9[ma_tmp],st0_Tw[9],p);
							MulMod(A_B1R10[ma_tmp],A_B1R10[ma_tmp],st0_Tw[10],p);
							MulMod(A_B1R11[ma_tmp],A_B1R11[ma_tmp],st0_Tw[11],p);
							MulMod(A_B1R12[ma_tmp],A_B1R12[ma_tmp],st0_Tw[12],p);
							MulMod(A_B1R13[ma_tmp],A_B1R13[ma_tmp],st0_Tw[13],p);
							MulMod(A_B1R14[ma_tmp],A_B1R14[ma_tmp],st0_Tw[14],p);
							MulMod(A_B1R15[ma_tmp],A_B1R15[ma_tmp],st0_Tw[15],p);
							break;
						case 1:
							DATARECORD <<" A_B1R0["<<ma_tmp<<"]:  "<<A_B1R0[ma_tmp]  << ", st1_Tw[0] = "  << st1_Tw[0]  << endl;
				    		DATARECORD <<" A_B1R1["<<ma_tmp<<"]:  "<<A_B1R1[ma_tmp]  << ", st1_Tw[1] = "  << st1_Tw[1]  << endl;
							DATARECORD <<" A_B1R2["<<ma_tmp<<"]:  "<<A_B1R2[ma_tmp]  << ", st1_Tw[2] = "  << st1_Tw[2]  << endl;
							DATARECORD <<" A_B1R3["<<ma_tmp<<"]:  "<<A_B1R3[ma_tmp]  << ", st1_Tw[3] = "  << st1_Tw[3]  << endl;
							DATARECORD <<" A_B1R4["<<ma_tmp<<"]:  "<<A_B1R4[ma_tmp]  << ", st1_Tw[4] = "  << st1_Tw[4]  << endl;
							DATARECORD <<" A_B1R5["<<ma_tmp<<"]:  "<<A_B1R5[ma_tmp]  << ", st1_Tw[5] = "  << st1_Tw[5]  << endl;
							DATARECORD <<" A_B1R6["<<ma_tmp<<"]:  "<<A_B1R6[ma_tmp]  << ", st1_Tw[6] = "  << st1_Tw[6]  << endl;
							DATARECORD <<" A_B1R7["<<ma_tmp<<"]:  "<<A_B1R7[ma_tmp]  << ", st1_Tw[7] = "  << st1_Tw[7]  << endl;
							DATARECORD <<" A_B1R8["<<ma_tmp<<"]:  "<<A_B1R8[ma_tmp]  << ", st1_Tw[8] = "  << st1_Tw[8]  << endl;
							DATARECORD <<" A_B1R9["<<ma_tmp<<"]:  "<<A_B1R9[ma_tmp]  << ", st1_Tw[9] = "  << st1_Tw[9]  << endl;
							DATARECORD <<" A_B1R10["<<ma_tmp<<"]: "<<A_B1R10[ma_tmp] << ", st1_Tw[10] = " << st1_Tw[10] << endl;
							DATARECORD <<" A_B1R11["<<ma_tmp<<"]: "<<A_B1R11[ma_tmp] << ", st1_Tw[11] = " << st1_Tw[11] << endl;
							DATARECORD <<" A_B1R12["<<ma_tmp<<"]: "<<A_B1R12[ma_tmp] << ", st1_Tw[12] = " << st1_Tw[12] << endl;
							DATARECORD <<" A_B1R13["<<ma_tmp<<"]: "<<A_B1R13[ma_tmp] << ", st1_Tw[13] = " << st1_Tw[13] << endl;
							DATARECORD <<" A_B1R14["<<ma_tmp<<"]: "<<A_B1R14[ma_tmp] << ", st1_Tw[14] = " << st1_Tw[14] << endl;
							DATARECORD <<" A_B1R15["<<ma_tmp<<"]: "<<A_B1R15[ma_tmp] << ", st1_Tw[15] = " << st1_Tw[15] << endl;
							Radix16_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp], A_B1R2[ma_tmp], A_B1R3[ma_tmp],A_B1R4[ma_tmp],
							   A_B1R5[ma_tmp],A_B1R6[ma_tmp], A_B1R7[ma_tmp], A_B1R8[ma_tmp],A_B1R9[ma_tmp],
							   A_B1R10[ma_tmp],A_B1R11[ma_tmp],A_B1R12[ma_tmp],A_B1R13[ma_tmp],A_B1R14[ma_tmp],
							   A_B1R15[ma_tmp]);
							MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],st1_Tw[0],p);
							MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],st1_Tw[1],p);
							MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],st1_Tw[2],p);
							MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],st1_Tw[3],p);
							MulMod(A_B1R4[ma_tmp],A_B1R4[ma_tmp],st1_Tw[4],p);
							MulMod(A_B1R5[ma_tmp],A_B1R5[ma_tmp],st1_Tw[5],p);
							MulMod(A_B1R6[ma_tmp],A_B1R6[ma_tmp],st1_Tw[6],p);
							MulMod(A_B1R7[ma_tmp],A_B1R7[ma_tmp],st1_Tw[7],p);
							MulMod(A_B1R8[ma_tmp],A_B1R8[ma_tmp],st1_Tw[8],p);
							MulMod(A_B1R9[ma_tmp],A_B1R9[ma_tmp],st1_Tw[9],p);
							MulMod(A_B1R10[ma_tmp],A_B1R10[ma_tmp],st1_Tw[10],p);
							MulMod(A_B1R11[ma_tmp],A_B1R11[ma_tmp],st1_Tw[11],p);
							MulMod(A_B1R12[ma_tmp],A_B1R12[ma_tmp],st1_Tw[12],p);
							MulMod(A_B1R13[ma_tmp],A_B1R13[ma_tmp],st1_Tw[13],p);
							MulMod(A_B1R14[ma_tmp],A_B1R14[ma_tmp],st1_Tw[14],p);
							MulMod(A_B1R15[ma_tmp],A_B1R15[ma_tmp],st1_Tw[15],p);
							break;
						case 2:
							DATARECORD <<" A_B1R0["<<ma_tmp<<"]:  "<<A_B1R0[ma_tmp]  << ", st2_Tw[0] = "  << st2_Tw[0]  << endl;
				    		DATARECORD <<" A_B1R1["<<ma_tmp<<"]:  "<<A_B1R1[ma_tmp]  << ", st2_Tw[1] = "  << st2_Tw[1]  << endl;
							DATARECORD <<" A_B1R2["<<ma_tmp<<"]:  "<<A_B1R2[ma_tmp]  << ", st2_Tw[2] = "  << st2_Tw[2]  << endl;
							DATARECORD <<" A_B1R3["<<ma_tmp<<"]:  "<<A_B1R3[ma_tmp]  << ", st2_Tw[3] = "  << st2_Tw[3]  << endl;
							DATARECORD <<" A_B1R4["<<ma_tmp<<"]:  "<<A_B1R4[ma_tmp]  << ", st2_Tw[4] = "  << st2_Tw[4]  << endl;
							DATARECORD <<" A_B1R5["<<ma_tmp<<"]:  "<<A_B1R5[ma_tmp]  << ", st2_Tw[5] = "  << st2_Tw[5]  << endl;
							DATARECORD <<" A_B1R6["<<ma_tmp<<"]:  "<<A_B1R6[ma_tmp]  << ", st2_Tw[6] = "  << st2_Tw[6]  << endl;
							DATARECORD <<" A_B1R7["<<ma_tmp<<"]:  "<<A_B1R7[ma_tmp]  << ", st2_Tw[7] = "  << st2_Tw[7]  << endl;
							DATARECORD <<" A_B1R8["<<ma_tmp<<"]:  "<<A_B1R8[ma_tmp]  << ", st2_Tw[8] = "  << st2_Tw[8]  << endl;
							DATARECORD <<" A_B1R9["<<ma_tmp<<"]:  "<<A_B1R9[ma_tmp]  << ", st2_Tw[9] = "  << st2_Tw[9]  << endl;
							DATARECORD <<" A_B1R10["<<ma_tmp<<"]: "<<A_B1R10[ma_tmp] << ", st2_Tw[10] = " << st2_Tw[10] << endl;
							DATARECORD <<" A_B1R11["<<ma_tmp<<"]: "<<A_B1R11[ma_tmp] << ", st2_Tw[11] = " << st2_Tw[11] << endl;
							DATARECORD <<" A_B1R12["<<ma_tmp<<"]: "<<A_B1R12[ma_tmp] << ", st2_Tw[12] = " << st2_Tw[12] << endl;
							DATARECORD <<" A_B1R13["<<ma_tmp<<"]: "<<A_B1R13[ma_tmp] << ", st2_Tw[13] = " << st2_Tw[13] << endl;
							DATARECORD <<" A_B1R14["<<ma_tmp<<"]: "<<A_B1R14[ma_tmp] << ", st2_Tw[14] = " << st2_Tw[14] << endl;
							DATARECORD <<" A_B1R15["<<ma_tmp<<"]: "<<A_B1R15[ma_tmp] << ", st2_Tw[15] = " << st2_Tw[15] << endl;
							Radix16_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp], A_B1R2[ma_tmp], A_B1R3[ma_tmp],A_B1R4[ma_tmp],
							   A_B1R5[ma_tmp],A_B1R6[ma_tmp], A_B1R7[ma_tmp], A_B1R8[ma_tmp],A_B1R9[ma_tmp],
							   A_B1R10[ma_tmp],A_B1R11[ma_tmp],A_B1R12[ma_tmp],A_B1R13[ma_tmp],A_B1R14[ma_tmp],
							   A_B1R15[ma_tmp]);
							MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],st2_Tw[0],p);
							MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],st2_Tw[1],p);
							MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],st2_Tw[2],p);
							MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],st2_Tw[3],p);
							MulMod(A_B1R4[ma_tmp],A_B1R4[ma_tmp],st2_Tw[4],p);
							MulMod(A_B1R5[ma_tmp],A_B1R5[ma_tmp],st2_Tw[5],p);
							MulMod(A_B1R6[ma_tmp],A_B1R6[ma_tmp],st2_Tw[6],p);
							MulMod(A_B1R7[ma_tmp],A_B1R7[ma_tmp],st2_Tw[7],p);
							MulMod(A_B1R8[ma_tmp],A_B1R8[ma_tmp],st2_Tw[8],p);
							MulMod(A_B1R9[ma_tmp],A_B1R9[ma_tmp],st2_Tw[9],p);
							MulMod(A_B1R10[ma_tmp],A_B1R10[ma_tmp],st2_Tw[10],p);
							MulMod(A_B1R11[ma_tmp],A_B1R11[ma_tmp],st2_Tw[11],p);
							MulMod(A_B1R12[ma_tmp],A_B1R12[ma_tmp],st2_Tw[12],p);
							MulMod(A_B1R13[ma_tmp],A_B1R13[ma_tmp],st2_Tw[13],p);
							MulMod(A_B1R14[ma_tmp],A_B1R14[ma_tmp],st2_Tw[14],p);
							MulMod(A_B1R15[ma_tmp],A_B1R15[ma_tmp],st2_Tw[15],p);
							break;
						case 3:
							DATARECORD <<" A_B1R0["<<ma_tmp<<"]:  "<<A_B1R0[ma_tmp]  << ", st3_Tw[0] = "  << 1 << endl;
				    		DATARECORD <<" A_B1R1["<<ma_tmp<<"]:  "<<A_B1R1[ma_tmp]  << ", st3_Tw[1] = "  << 1 << endl;
							DATARECORD <<" A_B1R2["<<ma_tmp<<"]:  "<<A_B1R2[ma_tmp]  << ", st3_Tw[2] = "  << 1 << endl;
							DATARECORD <<" A_B1R3["<<ma_tmp<<"]:  "<<A_B1R3[ma_tmp]  << ", st3_Tw[3] = "  << 1 << endl;
							DATARECORD <<" A_B1R4["<<ma_tmp<<"]:  "<<A_B1R4[ma_tmp]  << ", st3_Tw[4] = "  << 1 << endl;
							DATARECORD <<" A_B1R5["<<ma_tmp<<"]:  "<<A_B1R5[ma_tmp]  << ", st3_Tw[5] = "  << 1 << endl;
							DATARECORD <<" A_B1R6["<<ma_tmp<<"]:  "<<A_B1R6[ma_tmp]  << ", st3_Tw[6] = "  << 1 << endl;
							DATARECORD <<" A_B1R7["<<ma_tmp<<"]:  "<<A_B1R7[ma_tmp]  << ", st3_Tw[7] = "  << 1 << endl;
							DATARECORD <<" A_B1R8["<<ma_tmp<<"]:  "<<A_B1R8[ma_tmp]  << ", st3_Tw[8] = "  << 1 << endl;
							DATARECORD <<" A_B1R9["<<ma_tmp<<"]:  "<<A_B1R9[ma_tmp]  << ", st3_Tw[9] = "  << 1 << endl;
							DATARECORD <<" A_B1R10["<<ma_tmp<<"]: "<<A_B1R10[ma_tmp] << ", st3_Tw[10] = " << 1 << endl;
							DATARECORD <<" A_B1R11["<<ma_tmp<<"]: "<<A_B1R11[ma_tmp] << ", st3_Tw[11] = " << 1 << endl;
							DATARECORD <<" A_B1R12["<<ma_tmp<<"]: "<<A_B1R12[ma_tmp] << ", st3_Tw[12] = " << 1 << endl;
							DATARECORD <<" A_B1R13["<<ma_tmp<<"]: "<<A_B1R13[ma_tmp] << ", st3_Tw[13] = " << 1 << endl;
							DATARECORD <<" A_B1R14["<<ma_tmp<<"]: "<<A_B1R14[ma_tmp] << ", st3_Tw[14] = " << 1 << endl;
							DATARECORD <<" A_B1R15["<<ma_tmp<<"]: "<<A_B1R15[ma_tmp] << ", st3_Tw[15] = " << 1 << endl;
							Radix16_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp], A_B1R2[ma_tmp], A_B1R3[ma_tmp],A_B1R4[ma_tmp],
							   A_B1R5[ma_tmp],A_B1R6[ma_tmp], A_B1R7[ma_tmp], A_B1R8[ma_tmp],A_B1R9[ma_tmp],
							   A_B1R10[ma_tmp],A_B1R11[ma_tmp],A_B1R12[ma_tmp],A_B1R13[ma_tmp],A_B1R14[ma_tmp],
							   A_B1R15[ma_tmp]);
							MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],1,p);
							MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],1,p);
							MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],1,p);
							MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],1,p);
							MulMod(A_B1R4[ma_tmp],A_B1R4[ma_tmp],1,p);
							MulMod(A_B1R5[ma_tmp],A_B1R5[ma_tmp],1,p);
							MulMod(A_B1R6[ma_tmp],A_B1R6[ma_tmp],1,p);
							MulMod(A_B1R7[ma_tmp],A_B1R7[ma_tmp],1,p);
							MulMod(A_B1R8[ma_tmp],A_B1R8[ma_tmp],1,p);
							MulMod(A_B1R9[ma_tmp],A_B1R9[ma_tmp],1,p);
							MulMod(A_B1R10[ma_tmp],A_B1R10[ma_tmp],1,p);
							MulMod(A_B1R11[ma_tmp],A_B1R11[ma_tmp],1,p);
							MulMod(A_B1R12[ma_tmp],A_B1R12[ma_tmp],1,p);
							MulMod(A_B1R13[ma_tmp],A_B1R13[ma_tmp],1,p);
							MulMod(A_B1R14[ma_tmp],A_B1R14[ma_tmp],1,p);
							MulMod(A_B1R15[ma_tmp],A_B1R15[ma_tmp],1,p);
							break;
					}
					DATARECORD <<" ------after BU compute and Mul-------" << std::endl;
					DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp]<<"\n";
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
					DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<<A_B0R15[ma_tmp]<<"\n";	

					DATARECORD <<"--------------------------------------------------------------------\n";
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