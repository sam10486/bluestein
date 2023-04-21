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
				int BU_index_0;
				int BU_index_1;
				BU_index_0 = RR.RR(BC_tmp, s, N, log2(radix));
				BU_index_1 = BU_index_0 + ((N/radix)>>s);
				if(bn_tmp == 0){
                    DIT_DATARECORD << "bn_tmp = " << bn_tmp << ", ma_tmp = " << ma_tmp << std::endl;
					bn0_bc_tmp = BC_tmp;
					switch(s){
						case 0:
							DIT_DATARECORD << "idx:" << BU_index_0 << ", A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp]<< ", st0_Tw[0] = " << st0_Tw[0] << "\n";
							DIT_DATARECORD << "idx:" << BU_index_1 << ", A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp]<< ", st0_Tw[1] = " << st0_Tw[1] << "\n";
							if(!debug) Radix2_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp]);
							DIT_DATARECORD << "---after BU compute---" << std::endl;
                    		DIT_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp]<<"\n";
							DIT_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp]<<"\n";
							if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],st0_Tw[0],p);
							if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],st0_Tw[1],p);
							break;
						case 1:
							DIT_DATARECORD << "idx:" << BU_index_0 << ", A_B0R0["<<ma_tmp<<"]: " << A_B0R0[ma_tmp]<< ", st1_Tw[0] = " << st1_Tw[0] << "\n";
							DIT_DATARECORD << "idx:" << BU_index_1 << ", A_B0R1["<<ma_tmp<<"]: " << A_B0R1[ma_tmp]<< ", st1_Tw[1] = " << st1_Tw[1] << "\n";
							if(!debug) Radix2_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp]);
							DIT_DATARECORD << "---after BU compute---" << std::endl;
                    		DIT_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp]<<"\n";
							DIT_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp]<<"\n";
							if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],st1_Tw[0],p);
							if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],st1_Tw[1],p);
							break;
						case 2:
							DIT_DATARECORD << "idx:" << BU_index_0 << ", A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp]<< ", st2_Tw[0] = " << st2_Tw[0] << "\n";
							DIT_DATARECORD << "idx:" << BU_index_1 << ", A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp]<< ", st2_Tw[1] = " << st2_Tw[1] << "\n";
							if(!debug) Radix2_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp]);
							DIT_DATARECORD << "---after BU compute---" << std::endl;
                    		DIT_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp]<<"\n";
							DIT_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp]<<"\n";
							if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],st2_Tw[0],p);
							if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],st2_Tw[1],p);
							break;
						case 3:
							DIT_DATARECORD << "idx:" << BU_index_0 << ", A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp]<< ", st3_Tw[0] = " << 1 << "\n";
							DIT_DATARECORD << "idx:" << BU_index_1 << ", A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp]<< ", st3_Tw[1] = " << 1 << "\n";
							if(!debug) Radix2_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp]);
							DIT_DATARECORD << "---after BU compute---" << std::endl;
                    		DIT_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp]<<"\n";
							DIT_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp]<<"\n";
							if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],1,p);
							if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],1,p);
							break;
					}
					
				    
					bn0_ma_reg = ma_tmp;

				}
			    else {
                    DIT_DATARECORD << "bn_tmp = " << bn_tmp << ", ma_tmp = " << ma_tmp << std::endl;
					bn1_bc_tmp = BC_tmp;
					switch(s){
						case 0:
							DIT_DATARECORD << "idx:" << BU_index_0 << ", A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<< ", st0_Tw[0] = " << st0_Tw[0] << "\n";
							DIT_DATARECORD << "idx:" << BU_index_1 << ", A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<< ", st0_Tw[1] = " << st0_Tw[1] << "\n";
							if(!debug) Radix2_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp]);
							DIT_DATARECORD << "---after BU compute---" << std::endl;
                    		DIT_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<<"\n";
							DIT_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<<"\n";
							if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],st0_Tw[0],p);
							if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],st0_Tw[1],p);
							break;
						case 1:
							DIT_DATARECORD << "idx:" << BU_index_0 << ", A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<< ", st1_Tw[0] = " << st1_Tw[0] << "\n";
							DIT_DATARECORD << "idx:" << BU_index_1 << ", A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<< ", st1_Tw[1] = " << st1_Tw[1] << "\n";
							if(!debug) Radix2_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp]);
							DIT_DATARECORD << "---after BU compute---" << std::endl;
                    		DIT_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<<"\n";
							DIT_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<<"\n";
							if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],st1_Tw[0],p);
							if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],st1_Tw[1],p);
							break;
						case 2:
							DIT_DATARECORD << "idx:" << BU_index_0 << ", A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<< ", st2_Tw[0] = " << st2_Tw[0] << "\n";
							DIT_DATARECORD << "idx:" << BU_index_1 << ", A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<< ", st2_Tw[1] = " << st2_Tw[1] << "\n";
							if(!debug) Radix2_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp]);
							DIT_DATARECORD << "---after BU compute---" << std::endl;
                    		DIT_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<<"\n";
							DIT_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<<"\n";
							if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],st2_Tw[0],p);
							if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],st2_Tw[1],p);
							break;
						case 3:
							DIT_DATARECORD << "idx:" << BU_index_0 << ", A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<< ", st3_Tw[0] = " << 1 << "\n";
							DIT_DATARECORD << "idx:" << BU_index_1 << ", A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<< ", st3_Tw[1] = " << 1 << "\n";
							if(!debug) Radix2_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp]);
							DIT_DATARECORD << "---after BU compute---" << std::endl;
                    		DIT_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<<"\n";
							DIT_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<<"\n";
							if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],1,p);
							if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],1,p);
							break;
					}
				    
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
	
	std::ofstream DIT_DATARECORD("./DIT_NTT_R4_SPMB.txt");
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
	//----------------------------------------
	
    Stage = (unsigned long)ceil(log2(N));
	Stage = (unsigned long)Stage / 2;
	BC_WIDTH  = (int)ceil(log2(N/4));
	offset    =  (int)N /radix;
	word_size =  (int)N / (2 * radix);
	group     =  (int)N / (radix * radix);
	bit_array_tmp.resize(BC_WIDTH);
	DIT_DATARECORD <<"Stage: "<< Stage<<"\n";
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
		DIT_DATARECORD <<"---------------------------------\n";
		DIT_DATARECORD <<"Now Stage: "<< s <<"\n";
		

		DIT_spmb_radix4 <<"Now Stage: "<< s <<"\n";
		DIT_spmb_radix4 <<"twiddle factor : "<< factor <<"\n";
	    DIT_spmb_radix4 << "********\n";
		for(int i = 0 ;i < group;i++){
			DIT_DATARECORD <<"twiddle factor : "<< factor <<"\n";
			DIT_DATARECORD <<"p : "<< p <<"\n";
	    	DIT_DATARECORD << "********\n";
			bn0_bc_tmp  = 0;
			bn1_bc_tmp  = 0;
			DIT_spmb_radix4 <<"--------------i = " << i << "----------------\n";
			for(int j = 0;j < radix;j++){
				gray_i  = Gray(i,group);
				BC_tmp  = j * group + gray_i;
				DIT_DATARECORD << "i: " << i <<"\n";
				DIT_DATARECORD << "gray_i: " << gray_i <<"\n";
				DIT_DATARECORD << "BC_tmp: " << BC_tmp <<"\n";
				RR_R4(BC_tmp,s,BC);
				DIT_DATARECORD << "BC: " << BC <<"\n";
				length = BC_tmp >> (2*s);
				DIT_DATARECORD << "length: " <<  length <<"\n";

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
				DIT_DATARECORD << "factor_t: "<<factor_t<<"\n";
				AGU_R4(BC,bn_tmp,ma_tmp);
				DIT_DATARECORD << "bn_tmp: "<<bn_tmp<<"\n";
				if(bn_tmp == 0){
					if(j < 2)bn0_bc_tmp = BC_tmp;
					switch(s){
						case 0:
							DIT_DATARECORD <<" A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp] << ", st0_Tw[0] = " << st0_Tw[0] << endl;
				    		DIT_DATARECORD <<" A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp] << ", st0_Tw[1] = " << st0_Tw[1] << endl;
							DIT_DATARECORD <<" A_B0R2["<<ma_tmp<<"]: "<<A_B0R2[ma_tmp] << ", st0_Tw[2] = " << st0_Tw[2] << endl;
							DIT_DATARECORD <<" A_B0R3["<<ma_tmp<<"]: "<<A_B0R3[ma_tmp] << ", st0_Tw[3] = " << st0_Tw[3] << endl;
							if(!debug) Radix4_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp],A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
							DIT_DATARECORD <<" ------after BU compute and Mul-------" << std::endl;
							DIT_DATARECORD <<" A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp] << endl;
							DIT_DATARECORD <<" A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp] << endl;
							DIT_DATARECORD <<" A_B0R2["<<ma_tmp<<"]: "<<A_B0R2[ma_tmp] << endl;
							DIT_DATARECORD <<" A_B0R3["<<ma_tmp<<"]: "<<A_B0R3[ma_tmp] << endl;
							if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],st0_Tw[0],p);
							if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],st0_Tw[1],p);
							if(!debug) MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],st0_Tw[2],p);
							if(!debug) MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],st0_Tw[3],p);
							break;
						case 1:
							DIT_DATARECORD <<" A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp] << ", st1_Tw[0] = " << st1_Tw[0] << endl;
							DIT_DATARECORD <<" A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp] << ", st1_Tw[1] = " << st1_Tw[1] << endl;
							DIT_DATARECORD <<" A_B0R2["<<ma_tmp<<"]: "<<A_B0R2[ma_tmp] << ", st1_Tw[2] = " << st1_Tw[2] << endl;
							DIT_DATARECORD <<" A_B0R3["<<ma_tmp<<"]: "<<A_B0R3[ma_tmp] << ", st1_Tw[3] = " << st1_Tw[3] << endl;
							if(!debug) Radix4_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp],A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
							DIT_DATARECORD <<" ------after BU compute and Mul-------" << std::endl;
							DIT_DATARECORD <<" A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp] << endl;
							DIT_DATARECORD <<" A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp] << endl;
							DIT_DATARECORD <<" A_B0R2["<<ma_tmp<<"]: "<<A_B0R2[ma_tmp] << endl;
							DIT_DATARECORD <<" A_B0R3["<<ma_tmp<<"]: "<<A_B0R3[ma_tmp] << endl;
							if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],st1_Tw[0],p);
							if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],st1_Tw[1],p);
							if(!debug) MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],st1_Tw[2],p);
							if(!debug) MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],st1_Tw[3],p);
							break;
						case 2:
							DIT_DATARECORD <<" A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp] << ", st2_Tw[0] = " << st2_Tw[0] << endl;
							DIT_DATARECORD <<" A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp] << ", st2_Tw[1] = " << st2_Tw[1] << endl;
							DIT_DATARECORD <<" A_B0R2["<<ma_tmp<<"]: "<<A_B0R2[ma_tmp] << ", st2_Tw[2] = " << st2_Tw[2] << endl;
							DIT_DATARECORD <<" A_B0R3["<<ma_tmp<<"]: "<<A_B0R3[ma_tmp] << ", st2_Tw[3] = " << st2_Tw[3] << endl;
							if(!debug) Radix4_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp],A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
							DIT_DATARECORD <<" ------after BU compute and Mul-------" << std::endl;
							DIT_DATARECORD <<" A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp] << endl;
							DIT_DATARECORD <<" A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp] << endl;
							DIT_DATARECORD <<" A_B0R2["<<ma_tmp<<"]: "<<A_B0R2[ma_tmp] << endl;
							DIT_DATARECORD <<" A_B0R3["<<ma_tmp<<"]: "<<A_B0R3[ma_tmp] << endl;
							if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],st2_Tw[0],p);
							if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],st2_Tw[1],p);
							if(!debug) MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],st2_Tw[2],p);
							if(!debug) MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],st2_Tw[3],p);
							break;
						case 3:
							DIT_DATARECORD <<" A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp] << ", st3_Tw[0] = " << 1 << endl;
							DIT_DATARECORD <<" A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp] << ", st3_Tw[1] = " << 1 << endl;
							DIT_DATARECORD <<" A_B0R2["<<ma_tmp<<"]: "<<A_B0R2[ma_tmp] << ", st3_Tw[2] = " << 1 << endl;
							DIT_DATARECORD <<" A_B0R3["<<ma_tmp<<"]: "<<A_B0R3[ma_tmp] << ", st3_Tw[3] = " << 1 << endl;
							if(!debug) Radix4_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp],A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
							DIT_DATARECORD <<" ------after BU compute and Mul-------" << std::endl;
							DIT_DATARECORD <<" A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp] << endl;
							DIT_DATARECORD <<" A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp] << endl;
							DIT_DATARECORD <<" A_B0R2["<<ma_tmp<<"]: "<<A_B0R2[ma_tmp] << endl;
							DIT_DATARECORD <<" A_B0R3["<<ma_tmp<<"]: "<<A_B0R3[ma_tmp] << endl;
							if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],1,p);
							if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],1,p);
							if(!debug) MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],1,p);
							if(!debug) MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],1,p);
							break;
					}
					//DIT_DATARECORD <<" ------after BU compute and Mul-------" << std::endl;
					//DIT_DATARECORD <<" A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp] << endl;
					//DIT_DATARECORD <<" A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp] << endl;
					//DIT_DATARECORD <<" A_B0R2["<<ma_tmp<<"]: "<<A_B0R2[ma_tmp] << endl;
					//DIT_DATARECORD <<" A_B0R3["<<ma_tmp<<"]: "<<A_B0R3[ma_tmp] << endl;
					DIT_DATARECORD <<"--------------------------------------------------------------------\n";
					if(j <  2)bn0_ma_reg1 = ma_tmp;
					if(j >= 2)bn0_ma_reg2 = ma_tmp;
				}
			    else {
					if(j < 2)bn1_bc_tmp = BC_tmp;
					switch(s){
						case 0:
							DIT_DATARECORD <<" A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp] << ", st0_Tw[0] = " << st0_Tw[0] << endl;
							DIT_DATARECORD <<" A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp] << ", st0_Tw[1] = " << st0_Tw[1] << endl;
							DIT_DATARECORD <<" A_B1R2["<<ma_tmp<<"]: "<<A_B1R2[ma_tmp] << ", st0_Tw[2] = " << st0_Tw[2] << endl;
							DIT_DATARECORD <<" A_B1R3["<<ma_tmp<<"]: "<<A_B1R3[ma_tmp] << ", st0_Tw[3] = " << st0_Tw[3] << endl;
							if(!debug) Radix4_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp],A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
							DIT_DATARECORD <<" ------after BU compute and Mul-------" << std::endl;
							DIT_DATARECORD <<" A_B1R0["<<ma_tmp<<"]: " << A_B1R0[ma_tmp] << endl;
							DIT_DATARECORD <<" A_B1R1["<<ma_tmp<<"]: " << A_B1R1[ma_tmp] << endl;
							DIT_DATARECORD <<" A_B1R2["<<ma_tmp<<"]: " << A_B1R2[ma_tmp] << endl;
							DIT_DATARECORD <<" A_B1R3["<<ma_tmp<<"]: " << A_B1R3[ma_tmp] << endl;
							if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],st0_Tw[0],p);
							if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],st0_Tw[1],p);
							if(!debug) MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],st0_Tw[2],p);
							if(!debug) MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],st0_Tw[3],p);
							break;
						case 1:
							DIT_DATARECORD <<" A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp] << ", st1_Tw[0] = " << st1_Tw[0] << endl;
							DIT_DATARECORD <<" A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp] << ", st1_Tw[1] = " << st1_Tw[1] << endl;
							DIT_DATARECORD <<" A_B1R2["<<ma_tmp<<"]: "<<A_B1R2[ma_tmp] << ", st1_Tw[2] = " << st1_Tw[2] << endl;
							DIT_DATARECORD <<" A_B1R3["<<ma_tmp<<"]: "<<A_B1R3[ma_tmp] << ", st1_Tw[3] = " << st1_Tw[3] << endl;
							if(!debug) Radix4_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp],A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
							DIT_DATARECORD <<" ------after BU compute and Mul-------" << std::endl;
							DIT_DATARECORD <<" A_B1R0["<<ma_tmp<<"]: " << A_B1R0[ma_tmp] << endl;
							DIT_DATARECORD <<" A_B1R1["<<ma_tmp<<"]: " << A_B1R1[ma_tmp] << endl;
							DIT_DATARECORD <<" A_B1R2["<<ma_tmp<<"]: " << A_B1R2[ma_tmp] << endl;
							DIT_DATARECORD <<" A_B1R3["<<ma_tmp<<"]: " << A_B1R3[ma_tmp] << endl;
							if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],st1_Tw[0],p);
							if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],st1_Tw[1],p);
							if(!debug) MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],st1_Tw[2],p);
							if(!debug) MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],st1_Tw[3],p);
							break;
						case 2:
							DIT_DATARECORD <<" A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp] << ", st2_Tw[0] = " << st2_Tw[0] << endl;
							DIT_DATARECORD <<" A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp] << ", st2_Tw[1] = " << st2_Tw[1] << endl;
							DIT_DATARECORD <<" A_B1R2["<<ma_tmp<<"]: "<<A_B1R2[ma_tmp] << ", st2_Tw[2] = " << st2_Tw[2] << endl;
							DIT_DATARECORD <<" A_B1R3["<<ma_tmp<<"]: "<<A_B1R3[ma_tmp] << ", st2_Tw[3] = " << st2_Tw[3] << endl;
							if(!debug) Radix4_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp],A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
							DIT_DATARECORD <<" ------after BU compute and Mul-------" << std::endl;
							DIT_DATARECORD <<" A_B1R0["<<ma_tmp<<"]: " << A_B1R0[ma_tmp] << endl;
							DIT_DATARECORD <<" A_B1R1["<<ma_tmp<<"]: " << A_B1R1[ma_tmp] << endl;
							DIT_DATARECORD <<" A_B1R2["<<ma_tmp<<"]: " << A_B1R2[ma_tmp] << endl;
							DIT_DATARECORD <<" A_B1R3["<<ma_tmp<<"]: " << A_B1R3[ma_tmp] << endl;
							if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],st2_Tw[0],p);
							if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],st2_Tw[1],p);
							if(!debug) MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],st2_Tw[2],p);
							if(!debug) MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],st2_Tw[3],p);
							break;
						case 3:
							DIT_DATARECORD <<" A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp] << ", st3_Tw[0] = " << 1 << endl;
							DIT_DATARECORD <<" A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp] << ", st3_Tw[1] = " << 1 << endl;
							DIT_DATARECORD <<" A_B1R2["<<ma_tmp<<"]: "<<A_B1R2[ma_tmp] << ", st3_Tw[2] = " << 1 << endl;
							DIT_DATARECORD <<" A_B1R3["<<ma_tmp<<"]: "<<A_B1R3[ma_tmp] << ", st3_Tw[3] = " << 1 << endl;
							if(!debug) Radix4_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp],A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
							DIT_DATARECORD <<" ------after BU compute and Mul-------" << std::endl;
							DIT_DATARECORD <<" A_B1R0["<<ma_tmp<<"]: " << A_B1R0[ma_tmp] << endl;
							DIT_DATARECORD <<" A_B1R1["<<ma_tmp<<"]: " << A_B1R1[ma_tmp] << endl;
							DIT_DATARECORD <<" A_B1R2["<<ma_tmp<<"]: " << A_B1R2[ma_tmp] << endl;
							DIT_DATARECORD <<" A_B1R3["<<ma_tmp<<"]: " << A_B1R3[ma_tmp] << endl;
							if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],1,p);
							if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],1,p);
							if(!debug) MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],1,p);
							if(!debug) MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],1,p);
							break;
					}
			
					//DIT_DATARECORD <<" ------after BU compute and Mul-------" << std::endl;
					//DIT_DATARECORD <<" A_B1R0["<<ma_tmp<<"]: " << A_B1R0[ma_tmp] << endl;
					//DIT_DATARECORD <<" A_B1R1["<<ma_tmp<<"]: " << A_B1R1[ma_tmp] << endl;
					//DIT_DATARECORD <<" A_B1R2["<<ma_tmp<<"]: " << A_B1R2[ma_tmp] << endl;
					//DIT_DATARECORD <<" A_B1R3["<<ma_tmp<<"]: " << A_B1R3[ma_tmp] << endl;
					DIT_DATARECORD <<"--------------------------------------------------------------------\n";
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
	
	std::ofstream DIT_DATARECORD("./DIT_NTT_R16_SPMB.txt");

    Stage = (unsigned long)ceil(log2(N));
	Stage = (unsigned long)Stage / 4;
	BC_WIDTH  = (int)ceil(log2(N/16));
	offset    =  (int)N /radix;
	word_size =  (int)N / (2 * radix);
	group     =  (int)N / (radix * radix);
	bit_array_tmp.resize(BC_WIDTH);
	DIT_DATARECORD <<"Stage: "<< Stage<<"\n";

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
		DIT_DATARECORD <<"---------------------------------\n";
		DIT_DATARECORD <<"Now Stage: "<< s <<"\n";
		DIT_DATARECORD <<"twiddle factor : "<< factor <<"\n";
		std::cout << "twiddle factor : "<< factor <<"\n";
	    DIT_DATARECORD << "********\n";
		for(int i = 0 ;i < group;i++){
			bn0_bc_tmp  = 0;
			bn1_bc_tmp  = 0;
			for(int j = 0;j < radix;j++){
				gray_i  = Gray(i,group);
				BC_tmp  = j * group + gray_i;
				DIT_DATARECORD << "i: " << i <<"\n";
				DIT_DATARECORD << "gray_i: " << gray_i <<"\n";
				DIT_DATARECORD << "BC_tmp: " << BC_tmp <<"\n";
				RR_R16(BC_tmp,s,BC);
				DIT_DATARECORD << "BC: " << BC <<"\n";
				length = BC_tmp >> (4*s);
				DIT_DATARECORD << "length: " <<  length <<"\n";
				
				AGU_R16(BC,bn_tmp,ma_tmp);
				DIT_DATARECORD << "bn_tmp: "<<bn_tmp<<"\n";

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
							DIT_DATARECORD <<" A_B0R0["<<ma_tmp<<"]:  "<<A_B0R0[ma_tmp]	<< ", st0_Tw[0] = "  << st0_Tw[0] << endl;
				    		DIT_DATARECORD <<" A_B0R1["<<ma_tmp<<"]:  "<<A_B0R1[ma_tmp]	<< ", st0_Tw[1] = "  << st0_Tw[1] << endl;
							DIT_DATARECORD <<" A_B0R2["<<ma_tmp<<"]:  "<<A_B0R2[ma_tmp]	<< ", st0_Tw[2] = "  << st0_Tw[2] << endl;
							DIT_DATARECORD <<" A_B0R3["<<ma_tmp<<"]:  "<<A_B0R3[ma_tmp]	<< ", st0_Tw[3] = "  << st0_Tw[3] << endl;
							DIT_DATARECORD <<" A_B0R4["<<ma_tmp<<"]:  "<<A_B0R4[ma_tmp]	<< ", st0_Tw[4] = "  << st0_Tw[4] << endl;
							DIT_DATARECORD <<" A_B0R5["<<ma_tmp<<"]:  "<<A_B0R5[ma_tmp]	<< ", st0_Tw[5] = "  << st0_Tw[5] << endl;
							DIT_DATARECORD <<" A_B0R6["<<ma_tmp<<"]:  "<<A_B0R6[ma_tmp]	<< ", st0_Tw[6] = "  << st0_Tw[6] << endl;
							DIT_DATARECORD <<" A_B0R7["<<ma_tmp<<"]:  "<<A_B0R7[ma_tmp]	<< ", st0_Tw[7] = "  << st0_Tw[7] << endl;
							DIT_DATARECORD <<" A_B0R8["<<ma_tmp<<"]:  "<<A_B0R8[ma_tmp]	<< ", st0_Tw[8] = "  << st0_Tw[8] << endl;
							DIT_DATARECORD <<" A_B0R9["<<ma_tmp<<"]:  "<<A_B0R9[ma_tmp]	<< ", st0_Tw[9] = "  << st0_Tw[9] << endl;
							DIT_DATARECORD <<" A_B0R10["<<ma_tmp<<"]: "<<A_B0R10[ma_tmp]<< ", st0_Tw[10] = " << st0_Tw[10] << endl;
							DIT_DATARECORD <<" A_B0R11["<<ma_tmp<<"]: "<<A_B0R11[ma_tmp]<< ", st0_Tw[11] = " << st0_Tw[11] << endl;
							DIT_DATARECORD <<" A_B0R12["<<ma_tmp<<"]: "<<A_B0R12[ma_tmp]<< ", st0_Tw[12] = " << st0_Tw[12] << endl;
							DIT_DATARECORD <<" A_B0R13["<<ma_tmp<<"]: "<<A_B0R13[ma_tmp]<< ", st0_Tw[13] = " << st0_Tw[13] << endl;
							DIT_DATARECORD <<" A_B0R14["<<ma_tmp<<"]: "<<A_B0R14[ma_tmp]<< ", st0_Tw[14] = " << st0_Tw[14] << endl;
							DIT_DATARECORD <<" A_B0R15["<<ma_tmp<<"]: "<<A_B0R15[ma_tmp]<< ", st0_Tw[15] = " << st0_Tw[15] << endl;
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
							DIT_DATARECORD <<" A_B0R0["<<ma_tmp<<"]:  "<<A_B0R0[ma_tmp]	 << ", st1_Tw[0] = "  << st1_Tw[0] << endl;
				    		DIT_DATARECORD <<" A_B0R1["<<ma_tmp<<"]:  "<<A_B0R1[ma_tmp]	 << ", st1_Tw[1] = "  << st1_Tw[1] << endl;
							DIT_DATARECORD <<" A_B0R2["<<ma_tmp<<"]:  "<<A_B0R2[ma_tmp]	 << ", st1_Tw[2] = "  << st1_Tw[2] << endl;
							DIT_DATARECORD <<" A_B0R3["<<ma_tmp<<"]:  "<<A_B0R3[ma_tmp]	 << ", st1_Tw[3] = "  << st1_Tw[3] << endl;
							DIT_DATARECORD <<" A_B0R4["<<ma_tmp<<"]:  "<<A_B0R4[ma_tmp]	 << ", st1_Tw[4] = "  << st1_Tw[4] << endl;
							DIT_DATARECORD <<" A_B0R5["<<ma_tmp<<"]:  "<<A_B0R5[ma_tmp]	 << ", st1_Tw[5] = "  << st1_Tw[5] << endl;
							DIT_DATARECORD <<" A_B0R6["<<ma_tmp<<"]:  "<<A_B0R6[ma_tmp]	 << ", st1_Tw[6] = "  << st1_Tw[6] << endl;
							DIT_DATARECORD <<" A_B0R7["<<ma_tmp<<"]:  "<<A_B0R7[ma_tmp]	 << ", st1_Tw[7] = "  << st1_Tw[7] << endl;
							DIT_DATARECORD <<" A_B0R8["<<ma_tmp<<"]:  "<<A_B0R8[ma_tmp]	 << ", st1_Tw[8] = "  << st1_Tw[8] << endl;
							DIT_DATARECORD <<" A_B0R9["<<ma_tmp<<"]:  "<<A_B0R9[ma_tmp]	 << ", st1_Tw[9] = "  << st1_Tw[9] << endl;
							DIT_DATARECORD <<" A_B0R10["<<ma_tmp<<"]: "<<A_B0R10[ma_tmp] << ", st1_Tw[10] = " << st1_Tw[10] << endl;
							DIT_DATARECORD <<" A_B0R11["<<ma_tmp<<"]: "<<A_B0R11[ma_tmp] << ", st1_Tw[11] = " << st1_Tw[11] << endl;
							DIT_DATARECORD <<" A_B0R12["<<ma_tmp<<"]: "<<A_B0R12[ma_tmp] << ", st1_Tw[12] = " << st1_Tw[12] << endl;
							DIT_DATARECORD <<" A_B0R13["<<ma_tmp<<"]: "<<A_B0R13[ma_tmp] << ", st1_Tw[13] = " << st1_Tw[13] << endl;
							DIT_DATARECORD <<" A_B0R14["<<ma_tmp<<"]: "<<A_B0R14[ma_tmp] << ", st1_Tw[14] = " << st1_Tw[14] << endl;
							DIT_DATARECORD <<" A_B0R15["<<ma_tmp<<"]: "<<A_B0R15[ma_tmp] << ", st1_Tw[15] = " << st1_Tw[15] << endl;
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
							DIT_DATARECORD <<" A_B0R0["<<ma_tmp<<"]:  "<<A_B0R0[ma_tmp]  << ", st2_Tw[0] = "  << st2_Tw[0] << endl;
				    		DIT_DATARECORD <<" A_B0R1["<<ma_tmp<<"]:  "<<A_B0R1[ma_tmp]  << ", st2_Tw[1] = "  << st2_Tw[1] << endl;
							DIT_DATARECORD <<" A_B0R2["<<ma_tmp<<"]:  "<<A_B0R2[ma_tmp]  << ", st2_Tw[2] = "  << st2_Tw[2] << endl;
							DIT_DATARECORD <<" A_B0R3["<<ma_tmp<<"]:  "<<A_B0R3[ma_tmp]  << ", st2_Tw[3] = "  << st2_Tw[3] << endl;
							DIT_DATARECORD <<" A_B0R4["<<ma_tmp<<"]:  "<<A_B0R4[ma_tmp]  << ", st2_Tw[4] = "  << st2_Tw[4] << endl;
							DIT_DATARECORD <<" A_B0R5["<<ma_tmp<<"]:  "<<A_B0R5[ma_tmp]  << ", st2_Tw[5] = "  << st2_Tw[5] << endl;
							DIT_DATARECORD <<" A_B0R6["<<ma_tmp<<"]:  "<<A_B0R6[ma_tmp]  << ", st2_Tw[6] = "  << st2_Tw[6] << endl;
							DIT_DATARECORD <<" A_B0R7["<<ma_tmp<<"]:  "<<A_B0R7[ma_tmp]  << ", st2_Tw[7] = "  << st2_Tw[7] << endl;
							DIT_DATARECORD <<" A_B0R8["<<ma_tmp<<"]:  "<<A_B0R8[ma_tmp]  << ", st2_Tw[8] = "  << st2_Tw[8] << endl;
							DIT_DATARECORD <<" A_B0R9["<<ma_tmp<<"]:  "<<A_B0R9[ma_tmp]  << ", st2_Tw[9] = "  << st2_Tw[9] << endl;
							DIT_DATARECORD <<" A_B0R10["<<ma_tmp<<"]: "<<A_B0R10[ma_tmp] << ", st2_Tw[10] = " << st2_Tw[10] << endl;
							DIT_DATARECORD <<" A_B0R11["<<ma_tmp<<"]: "<<A_B0R11[ma_tmp] << ", st2_Tw[11] = " << st2_Tw[11] << endl;
							DIT_DATARECORD <<" A_B0R12["<<ma_tmp<<"]: "<<A_B0R12[ma_tmp] << ", st2_Tw[12] = " << st2_Tw[12] << endl;
							DIT_DATARECORD <<" A_B0R13["<<ma_tmp<<"]: "<<A_B0R13[ma_tmp] << ", st2_Tw[13] = " << st2_Tw[13] << endl;
							DIT_DATARECORD <<" A_B0R14["<<ma_tmp<<"]: "<<A_B0R14[ma_tmp] << ", st2_Tw[14] = " << st2_Tw[14] << endl;
							DIT_DATARECORD <<" A_B0R15["<<ma_tmp<<"]: "<<A_B0R15[ma_tmp] << ", st2_Tw[15] = " << st2_Tw[15] << endl;
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
							DIT_DATARECORD <<" A_B0R0["<<ma_tmp<<"]:  "<<A_B0R0[ma_tmp]  << ", st3_Tw[0] = "  << 1 << endl;
				    		DIT_DATARECORD <<" A_B0R1["<<ma_tmp<<"]:  "<<A_B0R1[ma_tmp]  << ", st3_Tw[1] = "  << 1 << endl;
							DIT_DATARECORD <<" A_B0R2["<<ma_tmp<<"]:  "<<A_B0R2[ma_tmp]  << ", st3_Tw[2] = "  << 1 << endl;
							DIT_DATARECORD <<" A_B0R3["<<ma_tmp<<"]:  "<<A_B0R3[ma_tmp]  << ", st3_Tw[3] = "  << 1 << endl;
							DIT_DATARECORD <<" A_B0R4["<<ma_tmp<<"]:  "<<A_B0R4[ma_tmp]  << ", st3_Tw[4] = "  << 1 << endl;
							DIT_DATARECORD <<" A_B0R5["<<ma_tmp<<"]:  "<<A_B0R5[ma_tmp]  << ", st3_Tw[5] = "  << 1 << endl;
							DIT_DATARECORD <<" A_B0R6["<<ma_tmp<<"]:  "<<A_B0R6[ma_tmp]  << ", st3_Tw[6] = "  << 1 << endl;
							DIT_DATARECORD <<" A_B0R7["<<ma_tmp<<"]:  "<<A_B0R7[ma_tmp]  << ", st3_Tw[7] = "  << 1 << endl;
							DIT_DATARECORD <<" A_B0R8["<<ma_tmp<<"]:  "<<A_B0R8[ma_tmp]  << ", st3_Tw[8] = "  << 1 << endl;
							DIT_DATARECORD <<" A_B0R9["<<ma_tmp<<"]:  "<<A_B0R9[ma_tmp]  << ", st3_Tw[9] = "  << 1 << endl;
							DIT_DATARECORD <<" A_B0R10["<<ma_tmp<<"]: "<<A_B0R10[ma_tmp] << ", st3_Tw[10] = " << 1 << endl;
							DIT_DATARECORD <<" A_B0R11["<<ma_tmp<<"]: "<<A_B0R11[ma_tmp] << ", st3_Tw[11] = " << 1 << endl;
							DIT_DATARECORD <<" A_B0R12["<<ma_tmp<<"]: "<<A_B0R12[ma_tmp] << ", st3_Tw[12] = " << 1 << endl;
							DIT_DATARECORD <<" A_B0R13["<<ma_tmp<<"]: "<<A_B0R13[ma_tmp] << ", st3_Tw[13] = " << 1 << endl;
							DIT_DATARECORD <<" A_B0R14["<<ma_tmp<<"]: "<<A_B0R14[ma_tmp] << ", st3_Tw[14] = " << 1 << endl;
							DIT_DATARECORD <<" A_B0R15["<<ma_tmp<<"]: "<<A_B0R15[ma_tmp] << ", st3_Tw[15] = " << 1 << endl;
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
					DIT_DATARECORD <<" ------after BU compute and Mul-------" << std::endl;
				    DIT_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp]<<"\n";
				    DIT_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp]<<"\n";
					DIT_DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<<A_B0R2[ma_tmp]<<"\n";
					DIT_DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<<A_B0R3[ma_tmp]<<"\n";
					DIT_DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<<A_B0R4[ma_tmp]<<"\n";
					DIT_DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<<A_B0R5[ma_tmp]<<"\n";
					DIT_DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<<A_B0R6[ma_tmp]<<"\n";
					DIT_DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<<A_B0R7[ma_tmp]<<"\n";
					DIT_DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<<A_B0R8[ma_tmp]<<"\n";
					DIT_DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<<A_B0R9[ma_tmp]<<"\n";
					DIT_DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<<A_B0R10[ma_tmp]<<"\n";
					DIT_DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<<A_B0R11[ma_tmp]<<"\n";
					DIT_DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<<A_B0R12[ma_tmp]<<"\n";
					DIT_DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<<A_B0R13[ma_tmp]<<"\n";
					DIT_DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<<A_B0R14[ma_tmp]<<"\n";
					DIT_DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<<A_B0R15[ma_tmp]<<"\n";
					DIT_DATARECORD <<"--------------------------------------------------------------------\n";
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
							DIT_DATARECORD <<" A_B1R0["<<ma_tmp<<"]:  "<<A_B1R0[ma_tmp]  << ", st0_Tw[0] = "  << st0_Tw[0]  << endl;
				    		DIT_DATARECORD <<" A_B1R1["<<ma_tmp<<"]:  "<<A_B1R1[ma_tmp]  << ", st0_Tw[1] = "  << st0_Tw[1]  << endl;
							DIT_DATARECORD <<" A_B1R2["<<ma_tmp<<"]:  "<<A_B1R2[ma_tmp]  << ", st0_Tw[2] = "  << st0_Tw[2]  << endl;
							DIT_DATARECORD <<" A_B1R3["<<ma_tmp<<"]:  "<<A_B1R3[ma_tmp]  << ", st0_Tw[3] = "  << st0_Tw[3]  << endl;
							DIT_DATARECORD <<" A_B1R4["<<ma_tmp<<"]:  "<<A_B1R4[ma_tmp]  << ", st0_Tw[4] = "  << st0_Tw[4]  << endl;
							DIT_DATARECORD <<" A_B1R5["<<ma_tmp<<"]:  "<<A_B1R5[ma_tmp]  << ", st0_Tw[5] = "  << st0_Tw[5]  << endl;
							DIT_DATARECORD <<" A_B1R6["<<ma_tmp<<"]:  "<<A_B1R6[ma_tmp]  << ", st0_Tw[6] = "  << st0_Tw[6]  << endl;
							DIT_DATARECORD <<" A_B1R7["<<ma_tmp<<"]:  "<<A_B1R7[ma_tmp]  << ", st0_Tw[7] = "  << st0_Tw[7]  << endl;
							DIT_DATARECORD <<" A_B1R8["<<ma_tmp<<"]:  "<<A_B1R8[ma_tmp]  << ", st0_Tw[8] = "  << st0_Tw[8]  << endl;
							DIT_DATARECORD <<" A_B1R9["<<ma_tmp<<"]:  "<<A_B1R9[ma_tmp]  << ", st0_Tw[9] = "  << st0_Tw[9]  << endl;
							DIT_DATARECORD <<" A_B1R10["<<ma_tmp<<"]: "<<A_B1R10[ma_tmp] << ", st0_Tw[10] = " << st0_Tw[10] << endl;
							DIT_DATARECORD <<" A_B1R11["<<ma_tmp<<"]: "<<A_B1R11[ma_tmp] << ", st0_Tw[11] = " << st0_Tw[11] << endl;
							DIT_DATARECORD <<" A_B1R12["<<ma_tmp<<"]: "<<A_B1R12[ma_tmp] << ", st0_Tw[12] = " << st0_Tw[12] << endl;
							DIT_DATARECORD <<" A_B1R13["<<ma_tmp<<"]: "<<A_B1R13[ma_tmp] << ", st0_Tw[13] = " << st0_Tw[13] << endl;
							DIT_DATARECORD <<" A_B1R14["<<ma_tmp<<"]: "<<A_B1R14[ma_tmp] << ", st0_Tw[14] = " << st0_Tw[14] << endl;
							DIT_DATARECORD <<" A_B1R15["<<ma_tmp<<"]: "<<A_B1R15[ma_tmp] << ", st0_Tw[15] = " << st0_Tw[15] << endl;
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
							DIT_DATARECORD <<" A_B1R0["<<ma_tmp<<"]:  "<<A_B1R0[ma_tmp]  << ", st1_Tw[0] = "  << st1_Tw[0]  << endl;
				    		DIT_DATARECORD <<" A_B1R1["<<ma_tmp<<"]:  "<<A_B1R1[ma_tmp]  << ", st1_Tw[1] = "  << st1_Tw[1]  << endl;
							DIT_DATARECORD <<" A_B1R2["<<ma_tmp<<"]:  "<<A_B1R2[ma_tmp]  << ", st1_Tw[2] = "  << st1_Tw[2]  << endl;
							DIT_DATARECORD <<" A_B1R3["<<ma_tmp<<"]:  "<<A_B1R3[ma_tmp]  << ", st1_Tw[3] = "  << st1_Tw[3]  << endl;
							DIT_DATARECORD <<" A_B1R4["<<ma_tmp<<"]:  "<<A_B1R4[ma_tmp]  << ", st1_Tw[4] = "  << st1_Tw[4]  << endl;
							DIT_DATARECORD <<" A_B1R5["<<ma_tmp<<"]:  "<<A_B1R5[ma_tmp]  << ", st1_Tw[5] = "  << st1_Tw[5]  << endl;
							DIT_DATARECORD <<" A_B1R6["<<ma_tmp<<"]:  "<<A_B1R6[ma_tmp]  << ", st1_Tw[6] = "  << st1_Tw[6]  << endl;
							DIT_DATARECORD <<" A_B1R7["<<ma_tmp<<"]:  "<<A_B1R7[ma_tmp]  << ", st1_Tw[7] = "  << st1_Tw[7]  << endl;
							DIT_DATARECORD <<" A_B1R8["<<ma_tmp<<"]:  "<<A_B1R8[ma_tmp]  << ", st1_Tw[8] = "  << st1_Tw[8]  << endl;
							DIT_DATARECORD <<" A_B1R9["<<ma_tmp<<"]:  "<<A_B1R9[ma_tmp]  << ", st1_Tw[9] = "  << st1_Tw[9]  << endl;
							DIT_DATARECORD <<" A_B1R10["<<ma_tmp<<"]: "<<A_B1R10[ma_tmp] << ", st1_Tw[10] = " << st1_Tw[10] << endl;
							DIT_DATARECORD <<" A_B1R11["<<ma_tmp<<"]: "<<A_B1R11[ma_tmp] << ", st1_Tw[11] = " << st1_Tw[11] << endl;
							DIT_DATARECORD <<" A_B1R12["<<ma_tmp<<"]: "<<A_B1R12[ma_tmp] << ", st1_Tw[12] = " << st1_Tw[12] << endl;
							DIT_DATARECORD <<" A_B1R13["<<ma_tmp<<"]: "<<A_B1R13[ma_tmp] << ", st1_Tw[13] = " << st1_Tw[13] << endl;
							DIT_DATARECORD <<" A_B1R14["<<ma_tmp<<"]: "<<A_B1R14[ma_tmp] << ", st1_Tw[14] = " << st1_Tw[14] << endl;
							DIT_DATARECORD <<" A_B1R15["<<ma_tmp<<"]: "<<A_B1R15[ma_tmp] << ", st1_Tw[15] = " << st1_Tw[15] << endl;
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
							DIT_DATARECORD <<" A_B1R0["<<ma_tmp<<"]:  "<<A_B1R0[ma_tmp]  << ", st2_Tw[0] = "  << st2_Tw[0]  << endl;
				    		DIT_DATARECORD <<" A_B1R1["<<ma_tmp<<"]:  "<<A_B1R1[ma_tmp]  << ", st2_Tw[1] = "  << st2_Tw[1]  << endl;
							DIT_DATARECORD <<" A_B1R2["<<ma_tmp<<"]:  "<<A_B1R2[ma_tmp]  << ", st2_Tw[2] = "  << st2_Tw[2]  << endl;
							DIT_DATARECORD <<" A_B1R3["<<ma_tmp<<"]:  "<<A_B1R3[ma_tmp]  << ", st2_Tw[3] = "  << st2_Tw[3]  << endl;
							DIT_DATARECORD <<" A_B1R4["<<ma_tmp<<"]:  "<<A_B1R4[ma_tmp]  << ", st2_Tw[4] = "  << st2_Tw[4]  << endl;
							DIT_DATARECORD <<" A_B1R5["<<ma_tmp<<"]:  "<<A_B1R5[ma_tmp]  << ", st2_Tw[5] = "  << st2_Tw[5]  << endl;
							DIT_DATARECORD <<" A_B1R6["<<ma_tmp<<"]:  "<<A_B1R6[ma_tmp]  << ", st2_Tw[6] = "  << st2_Tw[6]  << endl;
							DIT_DATARECORD <<" A_B1R7["<<ma_tmp<<"]:  "<<A_B1R7[ma_tmp]  << ", st2_Tw[7] = "  << st2_Tw[7]  << endl;
							DIT_DATARECORD <<" A_B1R8["<<ma_tmp<<"]:  "<<A_B1R8[ma_tmp]  << ", st2_Tw[8] = "  << st2_Tw[8]  << endl;
							DIT_DATARECORD <<" A_B1R9["<<ma_tmp<<"]:  "<<A_B1R9[ma_tmp]  << ", st2_Tw[9] = "  << st2_Tw[9]  << endl;
							DIT_DATARECORD <<" A_B1R10["<<ma_tmp<<"]: "<<A_B1R10[ma_tmp] << ", st2_Tw[10] = " << st2_Tw[10] << endl;
							DIT_DATARECORD <<" A_B1R11["<<ma_tmp<<"]: "<<A_B1R11[ma_tmp] << ", st2_Tw[11] = " << st2_Tw[11] << endl;
							DIT_DATARECORD <<" A_B1R12["<<ma_tmp<<"]: "<<A_B1R12[ma_tmp] << ", st2_Tw[12] = " << st2_Tw[12] << endl;
							DIT_DATARECORD <<" A_B1R13["<<ma_tmp<<"]: "<<A_B1R13[ma_tmp] << ", st2_Tw[13] = " << st2_Tw[13] << endl;
							DIT_DATARECORD <<" A_B1R14["<<ma_tmp<<"]: "<<A_B1R14[ma_tmp] << ", st2_Tw[14] = " << st2_Tw[14] << endl;
							DIT_DATARECORD <<" A_B1R15["<<ma_tmp<<"]: "<<A_B1R15[ma_tmp] << ", st2_Tw[15] = " << st2_Tw[15] << endl;
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
							DIT_DATARECORD <<" A_B1R0["<<ma_tmp<<"]:  "<<A_B1R0[ma_tmp]  << ", st3_Tw[0] = "  << 1 << endl;
				    		DIT_DATARECORD <<" A_B1R1["<<ma_tmp<<"]:  "<<A_B1R1[ma_tmp]  << ", st3_Tw[1] = "  << 1 << endl;
							DIT_DATARECORD <<" A_B1R2["<<ma_tmp<<"]:  "<<A_B1R2[ma_tmp]  << ", st3_Tw[2] = "  << 1 << endl;
							DIT_DATARECORD <<" A_B1R3["<<ma_tmp<<"]:  "<<A_B1R3[ma_tmp]  << ", st3_Tw[3] = "  << 1 << endl;
							DIT_DATARECORD <<" A_B1R4["<<ma_tmp<<"]:  "<<A_B1R4[ma_tmp]  << ", st3_Tw[4] = "  << 1 << endl;
							DIT_DATARECORD <<" A_B1R5["<<ma_tmp<<"]:  "<<A_B1R5[ma_tmp]  << ", st3_Tw[5] = "  << 1 << endl;
							DIT_DATARECORD <<" A_B1R6["<<ma_tmp<<"]:  "<<A_B1R6[ma_tmp]  << ", st3_Tw[6] = "  << 1 << endl;
							DIT_DATARECORD <<" A_B1R7["<<ma_tmp<<"]:  "<<A_B1R7[ma_tmp]  << ", st3_Tw[7] = "  << 1 << endl;
							DIT_DATARECORD <<" A_B1R8["<<ma_tmp<<"]:  "<<A_B1R8[ma_tmp]  << ", st3_Tw[8] = "  << 1 << endl;
							DIT_DATARECORD <<" A_B1R9["<<ma_tmp<<"]:  "<<A_B1R9[ma_tmp]  << ", st3_Tw[9] = "  << 1 << endl;
							DIT_DATARECORD <<" A_B1R10["<<ma_tmp<<"]: "<<A_B1R10[ma_tmp] << ", st3_Tw[10] = " << 1 << endl;
							DIT_DATARECORD <<" A_B1R11["<<ma_tmp<<"]: "<<A_B1R11[ma_tmp] << ", st3_Tw[11] = " << 1 << endl;
							DIT_DATARECORD <<" A_B1R12["<<ma_tmp<<"]: "<<A_B1R12[ma_tmp] << ", st3_Tw[12] = " << 1 << endl;
							DIT_DATARECORD <<" A_B1R13["<<ma_tmp<<"]: "<<A_B1R13[ma_tmp] << ", st3_Tw[13] = " << 1 << endl;
							DIT_DATARECORD <<" A_B1R14["<<ma_tmp<<"]: "<<A_B1R14[ma_tmp] << ", st3_Tw[14] = " << 1 << endl;
							DIT_DATARECORD <<" A_B1R15["<<ma_tmp<<"]: "<<A_B1R15[ma_tmp] << ", st3_Tw[15] = " << 1 << endl;
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
					DIT_DATARECORD <<" ------after BU compute and Mul-------" << std::endl;
					DIT_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp]<<"\n";
					DIT_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp]<<"\n";
					DIT_DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<<A_B0R2[ma_tmp]<<"\n";
					DIT_DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<<A_B0R3[ma_tmp]<<"\n";
					DIT_DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<<A_B0R4[ma_tmp]<<"\n";
					DIT_DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<<A_B0R5[ma_tmp]<<"\n";
					DIT_DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<<A_B0R6[ma_tmp]<<"\n";
					DIT_DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<<A_B0R7[ma_tmp]<<"\n";
					DIT_DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<<A_B0R8[ma_tmp]<<"\n";
					DIT_DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<<A_B0R9[ma_tmp]<<"\n";
					DIT_DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<<A_B0R10[ma_tmp]<<"\n";
					DIT_DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<<A_B0R11[ma_tmp]<<"\n";
					DIT_DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<<A_B0R12[ma_tmp]<<"\n";
					DIT_DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<<A_B0R13[ma_tmp]<<"\n";
					DIT_DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<<A_B0R14[ma_tmp]<<"\n";
					DIT_DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<<A_B0R15[ma_tmp]<<"\n";	

					DIT_DATARECORD <<"--------------------------------------------------------------------\n";
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

void DIT_NTTSPMB::DIT_NTT_r4_r2(std::vector<ZZ> &A,
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
	
	std::ofstream DIT_DATARECORD("./DIT_NTT_R4_R2_SPMB.txt");
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
	int radix_r2 = 2;
	ZZ fft_twiddle = W;
	ZZ fft_prime = p;
	int debug = 1;
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
	DIT_DATARECORD <<"Stage: "<< Stage<<"\n";
	
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
	DIT_DATARECORD << "init load over! \n";
	int tw_degree = 1; // siang
	//need modify to mult by twiddle factor
	for(int s = 0; s < Stage; s++){
		DIT_DATARECORD <<"----------------------------"<<"\n";
		DIT_DATARECORD <<"Stage: "<< s<<"\n";
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
			for(int j = 0;j < radix;j++){
				DIT_DATARECORD <<"----------------------------"<<"\n";
				DIT_DATARECORD <<"i: "<< i << ", j = " << j <<"\n";
				DIT_DATARECORD <<"twiddle factor =  " << factor << std::endl;
				gray_i  = Gray(i,group);
				BC_tmp  = j * group + gray_i;
				if(s == Stage-1)RR_R4_R2(BC_tmp,(2*s - 1),BC);
				else RR_R4(BC_tmp,s,BC);
				length = BC % tw_modulus_tmp;
				PowerMod(factor_t,factor,length,p);
				AGU_R4(BC,bn_tmp,ma_tmp);

				DIT_DATARECORD << "BC_tmp = " << BC_tmp << ", length: " << length <<
				", tw_dg: " <<  tw_degree * length << " = " << tw_degree << " * " << length <<"\n";	// TF_record
				DIT_DATARECORD <<"factor_t =  " << factor_t << ", p = " <<  p <<"\n";
				//-----------DFTAG generator--------------
				DTFAG.DTFAG_SPMB_DIT(
								s, fft_point, radix_r1, radix_r2, debug,
								ROM0, ROM1, ROM2,
								st0_Tw, st1_Tw, st2_Tw, 
								DTFAG_t, DTFAG_i, DTFAG_j);
				switch(s){
					case 0:
						for(int i=0; i<radix; i++){
							//cout << "st0_Tw[" << i << "] = w^" << st0_Tw[i] << endl;
						}
						break;
					case 1:
						for(int i=0; i<radix; i++){
							//cout << "st1_Tw[" << i << "] = w^" << st1_Tw[i] << endl;
						}
						break;
					case 2:
						for(int i=0; i<radix; i++){
							//cout << "st2_Tw[" << i << "] = w^" << st2_Tw[i] << endl;
						}
						break;
				}
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
							DIT_DATARECORD <<"A_B0R0[" << ma_tmp << "] = " << A_B0R0[ma_tmp] << ", st0_Tw[0] = " << st0_Tw[0] << endl;
							DIT_DATARECORD <<"A_B0R1[" << ma_tmp << "] = " << A_B0R1[ma_tmp] << ", st0_Tw[1] = " << st0_Tw[1] << endl;
							DIT_DATARECORD <<"A_B0R2[" << ma_tmp << "] = " << A_B0R2[ma_tmp] << ", st0_Tw[2] = " << st0_Tw[2] << endl;
							DIT_DATARECORD <<"A_B0R3[" << ma_tmp << "] = " << A_B0R3[ma_tmp] << ", st0_Tw[3] = " << st0_Tw[3] << endl;
							if(!debug) Radix4_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp],A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
							if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],st0_Tw[0],p);
							if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],st0_Tw[1],p);
							if(!debug) MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],st0_Tw[2],p);
							if(!debug) MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],st0_Tw[3],p);
							break;
						case 1:
							DIT_DATARECORD <<"A_B0R0[" << ma_tmp << "] = " << A_B0R0[ma_tmp] << ", st1_Tw[0] = " << st1_Tw[0] << endl;
							DIT_DATARECORD <<"A_B0R1[" << ma_tmp << "] = " << A_B0R1[ma_tmp] << ", st1_Tw[1] = " << st1_Tw[1] << endl;
							DIT_DATARECORD <<"A_B0R2[" << ma_tmp << "] = " << A_B0R2[ma_tmp] << ", st1_Tw[2] = " << st1_Tw[2] << endl;
							DIT_DATARECORD <<"A_B0R3[" << ma_tmp << "] = " << A_B0R3[ma_tmp] << ", st1_Tw[3] = " << st1_Tw[3] << endl;
							if(!debug) Radix4_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp],A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
							if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],st1_Tw[0],p);
							if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],st1_Tw[1],p);
							if(!debug) MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],st1_Tw[2],p);
							if(!debug) MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],st1_Tw[3],p);
							break;
						case 2:
							DIT_DATARECORD <<"A_B0R0[" << ma_tmp << "] = " << A_B0R0[ma_tmp] << ", st2_Tw[0] = " << st2_Tw[0] << endl;
							DIT_DATARECORD <<"A_B0R1[" << ma_tmp << "] = " << A_B0R1[ma_tmp] << ", st2_Tw[1] = " << st2_Tw[1] << endl;
							DIT_DATARECORD <<"A_B0R2[" << ma_tmp << "] = " << A_B0R2[ma_tmp] << ", st2_Tw[2] = " << st2_Tw[2] << endl;
							DIT_DATARECORD <<"A_B0R3[" << ma_tmp << "] = " << A_B0R3[ma_tmp] << ", st2_Tw[3] = " << st2_Tw[3] << endl;
							if(!debug) Radix4_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp],A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
							if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],st2_Tw[0],p);
							if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],st2_Tw[1],p);
							if(!debug) MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],st2_Tw[2],p);
							if(!debug) MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],st2_Tw[3],p);
							break;
						case 3:
							DIT_DATARECORD <<"A_B0R0[" << ma_tmp << "] = " << A_B0R0[ma_tmp] << ", st3_Tw[0] = " << 1 << endl;
							DIT_DATARECORD <<"A_B0R1[" << ma_tmp << "] = " << A_B0R1[ma_tmp] << ", st3_Tw[1] = " << 1 << endl;
							DIT_DATARECORD <<"A_B0R2[" << ma_tmp << "] = " << A_B0R2[ma_tmp] << ", st3_Tw[2] = " << 1 << endl;
							DIT_DATARECORD <<"A_B0R3[" << ma_tmp << "] = " << A_B0R3[ma_tmp] << ", st3_Tw[3] = " << 1 << endl;
							if(!debug) Radix4_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp],A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
							if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],1,p);
							if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],1,p);
							if(!debug) MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],1,p);
							if(!debug) MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],1,p);
							break;
					}
					
					DIT_DATARECORD << "*****after Mul*****" << std::endl;
					DIT_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp] << endl;
					DIT_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp] << endl;
					DIT_DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp] << endl;
					DIT_DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp] << endl;

					if(j <  2)bn0_ma_reg1 = ma_tmp;
					if(j >= 2)bn0_ma_reg2 = ma_tmp;
				}
			    else {
					if(j < 2)bn1_bc_tmp = BC_tmp;
					switch(s){
						case 0:
							DIT_DATARECORD <<"A_B1R0[" << ma_tmp << "] = " << A_B1R0[ma_tmp] << ", st0_Tw[0] = " << st0_Tw[0] << endl;
							DIT_DATARECORD <<"A_B1R1[" << ma_tmp << "] = " << A_B1R1[ma_tmp] << ", st0_Tw[1] = " << st0_Tw[1] << endl;
							DIT_DATARECORD <<"A_B1R2[" << ma_tmp << "] = " << A_B1R2[ma_tmp] << ", st0_Tw[2] = " << st0_Tw[2] << endl;
							DIT_DATARECORD <<"A_B1R3[" << ma_tmp << "] = " << A_B1R3[ma_tmp] << ", st0_Tw[3] = " << st0_Tw[3] << endl;
							if(!debug) Radix4_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp],A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
							if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],st0_Tw[0],p);
							if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],st0_Tw[1],p);
							if(!debug) MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],st0_Tw[2],p);
							if(!debug) MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],st0_Tw[3],p);
							break;
						case 1:
							DIT_DATARECORD <<"A_B1R0[" << ma_tmp << "] = " << A_B1R0[ma_tmp] << ", st1_Tw[0] = " << st1_Tw[0] << endl;
							DIT_DATARECORD <<"A_B1R1[" << ma_tmp << "] = " << A_B1R1[ma_tmp] << ", st1_Tw[1] = " << st1_Tw[1] << endl;
							DIT_DATARECORD <<"A_B1R2[" << ma_tmp << "] = " << A_B1R2[ma_tmp] << ", st1_Tw[2] = " << st1_Tw[2] << endl;
							DIT_DATARECORD <<"A_B1R3[" << ma_tmp << "] = " << A_B1R3[ma_tmp] << ", st1_Tw[3] = " << st1_Tw[3] << endl;
							if(!debug) Radix4_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp],A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
							if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],st1_Tw[0],p);
							if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],st1_Tw[1],p);
							if(!debug) MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],st1_Tw[2],p);
							if(!debug) MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],st1_Tw[3],p);
							break;
						case 2:
							DIT_DATARECORD <<"A_B1R0[" << ma_tmp << "] = " << A_B1R0[ma_tmp] << ", st2_Tw[0] = " << st2_Tw[0] << endl;
							DIT_DATARECORD <<"A_B1R1[" << ma_tmp << "] = " << A_B1R1[ma_tmp] << ", st2_Tw[1] = " << st2_Tw[1] << endl;
							DIT_DATARECORD <<"A_B1R2[" << ma_tmp << "] = " << A_B1R2[ma_tmp] << ", st2_Tw[2] = " << st2_Tw[2] << endl;
							DIT_DATARECORD <<"A_B1R3[" << ma_tmp << "] = " << A_B1R3[ma_tmp] << ", st2_Tw[3] = " << st2_Tw[3] << endl;
							if(!debug) Radix4_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp],A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
							if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],st2_Tw[0],p);
							if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],st2_Tw[1],p);
							if(!debug) MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],st2_Tw[2],p);
							if(!debug) MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],st2_Tw[3],p);
							break;
						case 3:
							DIT_DATARECORD <<"A_B1R0[" << ma_tmp << "] = " << A_B1R0[ma_tmp] << ", st3_Tw[0] = " << 1 << endl;
							DIT_DATARECORD <<"A_B1R1[" << ma_tmp << "] = " << A_B1R1[ma_tmp] << ", st3_Tw[1] = " << 1 << endl;
							DIT_DATARECORD <<"A_B1R2[" << ma_tmp << "] = " << A_B1R2[ma_tmp] << ", st3_Tw[2] = " << 1 << endl;
							DIT_DATARECORD <<"A_B1R3[" << ma_tmp << "] = " << A_B1R3[ma_tmp] << ", st3_Tw[3] = " << 1 << endl;
							if(!debug) Radix4_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp],A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
							if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],1,p);
							if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],1,p);
							if(!debug) MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],1,p);
							if(!debug) MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],1,p);
							break;
					}
					
					DIT_DATARECORD << "*****after mul*****" << endl;
					DIT_DATARECORD << "A_B1R0[" <<ma_tmp  << "] = " << A_B1R0[ma_tmp] << endl;
					DIT_DATARECORD << "A_B1R1[" << ma_tmp << "] = " << A_B1R1[ma_tmp] << endl;
					DIT_DATARECORD << "A_B1R2[" << ma_tmp << "] = " << A_B1R2[ma_tmp] << endl;
					DIT_DATARECORD << "A_B1R3[" << ma_tmp << "] = " << A_B1R3[ma_tmp] << endl;
					if(j < 2) bn1_ma_reg1 = ma_tmp;
					if(j >= 2)bn1_ma_reg2 = ma_tmp;
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
	DIT_DATARECORD <<"------------radix_2 part----------------"<< std::endl;	
	// radix-2 FFT compute
	for(int i = 0; i < group ; i++){
		for(int j=0; j < radix; j++){
			DIT_DATARECORD <<"------------------"<< std::endl;	
			gray_i  = Gray(i,group);
			BC_tmp  = j * group + gray_i;
			BC = BC_tmp;
			AGU_R4(BC,bn_tmp,ma_tmp);
			if(bn_tmp == 0){
				if(j < 2)bn0_bc_tmp = BC_tmp;
				DIT_DATARECORD <<"A_B0R0[" << ma_tmp << "] = " << A_B0R0[ma_tmp] <<", A_B0R1[" << ma_tmp << "] = " << A_B0R1[ma_tmp] << std::endl;
				DIT_DATARECORD <<"A_B0R2[" << ma_tmp << "] = " << A_B0R2[ma_tmp] <<", A_B0R3[" << ma_tmp << "] = " << A_B0R3[ma_tmp] << std::endl;
				Radix2_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp]);
				Radix2_BU(A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
				if(j <  2)bn0_ma_reg1 = ma_tmp;
				if(j >= 2)bn0_ma_reg2 = ma_tmp;
			}else {
			    if(j < 2)bn1_bc_tmp = BC_tmp;
				DIT_DATARECORD <<"A_B1R0[" << ma_tmp << "] = " << A_B1R0[ma_tmp] <<", A_B1R1[" << ma_tmp << "] = " << A_B1R1[ma_tmp] << std::endl;
				DIT_DATARECORD <<"A_B1R2[" << ma_tmp << "] = " << A_B1R2[ma_tmp] <<", A_B1R3[" << ma_tmp << "] = " << A_B1R3[ma_tmp] << std::endl;
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

    DIT_DATARECORD << "---------------------------------------------------------\n";
	DIT_DATARECORD << "*********************************************************\n";
	DIT_DATARECORD << "                  FINAL!!!!                              \n";
	DIT_DATARECORD << "*********************************************************\n";
	DIT_DATARECORD << "---------------------------------------------------------\n";
	
    int BC_bit_size;
	BC_bit_size = (int)ceil(log2(N/4));
	DIT_DATARECORD << "BC_bit_size: "<< BC_bit_size <<"\n";
	int index_out;
	for(int i = 0; i < group; i++){
    	for(int j = 0;j < radix;j++){
			DIT_DATARECORD << "-------------------------" << std::endl;
			gray_i  = Gray(i,group);
    		BC_tmp  = j * group + gray_i;
			REORDERBC_R4_R2_OUT(BC_tmp,BC);
			BC = BC_tmp;
			DIT_DATARECORD << "BC_tmp: " << BC_tmp <<"\n";
			DIT_DATARECORD << "BC: " << BC <<"\n";
    		AGU_R4(BC,bn_tmp,ma_tmp);
			DIT_DATARECORD << "bn_tmp: " << bn_tmp <<"\n";
			DIT_DATARECORD << "ma_tmp: " << ma_tmp <<"\n";
			index_out = (int)i * ( N / group) + 4 * j ;
			DIT_DATARECORD << "index_out: " << index_out << std::endl;
    		if(bn_tmp == 0){
     		   A[index_out+0] = A_B0R0[ma_tmp];
    		   A[index_out+1] = A_B0R1[ma_tmp];
    		   A[index_out+2] = A_B0R2[ma_tmp];
    		   A[index_out+3] = A_B0R3[ma_tmp];
			   DIT_DATARECORD << "BN: 0 \n";
			   DIT_DATARECORD << "A["<< index_out      << "]: " << A[index_out+0] << " \n";
			   DIT_DATARECORD << "A["<< index_out + 1  << "]: " << A[index_out+1] << " \n";
			   DIT_DATARECORD << "A["<< index_out + 2  << "]: " << A[index_out+2] << " \n";
			   DIT_DATARECORD << "A["<< index_out + 3  << "]: " << A[index_out+3] << " \n";
    		}
    	    else {
    		   A[index_out+0] = A_B1R0[ma_tmp];
    		   A[index_out+1] = A_B1R1[ma_tmp];
    		   A[index_out+2] = A_B1R2[ma_tmp];
    		   A[index_out+3] = A_B1R3[ma_tmp];
			   DIT_DATARECORD << "BN: 1 \n";
			   DIT_DATARECORD << "A["<< index_out      << "]: " << A[index_out+0] << " \n";
			   DIT_DATARECORD << "A["<< index_out + 1  << "]: " << A[index_out+1] << " \n";
			   DIT_DATARECORD << "A["<< index_out + 2  << "]: " << A[index_out+2] << " \n";
			   DIT_DATARECORD << "A["<< index_out + 3  << "]: " << A[index_out+3] << " \n";
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
	DIT_DATARECORD << "---------------------------------------------------------\n";
	DIT_DATARECORD << "----                  OVER!!!!!!!                   -----\n";
	DIT_DATARECORD << "---------------------------------------------------------\n";
	DIT_DATARECORD << "---------------------------------------------------------\n";
	DIT_DATARECORD.close();
}

void DIT_NTTSPMB::test_radix2(std::vector<ZZ> &A){
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

	std::ofstream DATARECORD("./test_R2_SPMB.txt");
	std::ofstream test_radix2("./SPMB_tw/test_radix2.txt");
	//-------------Bit Reverse--------------
	BitOperate RR, IntToVec, VecToInt, BR, num_complement;
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
	int tw_degree = pow(2, Stage-1); // siang
	ZZ divisor;
	std::cout << "init load over! \n";
	//need modify to mult by twiddle factor
	for(int s = 0; s < Stage; s++){
		int s_complement = num_complement.number_complement(s, Stage);
		if(s == 0){
			//factor = PowerMod(W, (pow(radix, s_complement)), p);
		}
		else {
			//factor = PowerMod(W, (pow(radix, s_complement)), p);
			tw_degree = tw_degree / 2;
		}
		DATARECORD <<"	Stage: "<< s<<"\n";
		DATARECORD << "factor "<< factor <<"\n";
		DATARECORD << "****************************\n";

		test_radix2 <<"Now Stage: "<< s <<"\n";
		test_radix2 <<"twiddle factor : "<< factor <<"\n";
		for(int i = 0 ;i < group;i++){
			test_radix2 << "----------------i =" << i << " ----------------" << std::endl;
			DATARECORD << "----------------i =" << i << " ----------------" << std::endl;
			for(int j = 0;j < radix;j++){
				DATARECORD << "---j =" << j << " ---" << std::endl;
				DATARECORD <<"twiddle factor : "<< factor <<"\n";
				DATARECORD <<"p : "<< p <<"\n";
	    		DATARECORD << "********\n";
				gray_i  = Gray(i,group);
				BC_tmp  = j * group + gray_i;
				DATARECORD << "BC_tmp = " << BC_tmp << endl;
				RR_R2(BC_tmp,s,BC);
				int BC_tmp_mod = (BC_tmp % (1<<s));
				DATARECORD << "BC_tmp_mod = " << BC_tmp_mod << ", D_width = " << log2((1<<s)) << ", BR_value = " << BR.BitReserve(BC_tmp_mod, log2((1<<s))) << endl;
				length =  BR.BitReserve(BC_tmp_mod, log2((1<<s))) * pow(radix, s_complement);
				DATARECORD << "length: " <<  length <<"\n";
				PowerMod(factor_t,W,length,p);
				DATARECORD << "factor_t: "<<factor_t<<"\n";

				AGU_R2(BC,bn_tmp,ma_tmp);

				//-----------compute data idx-------------
				test_radix2 << "BC_tmp = " << BC_tmp << std::endl;
				
				IntToVec.IntToVec(BC_tmp, N, bit_array);
				rotate(bit_array.begin(), bit_array.begin()+ bit_width_s*s , bit_array.end());
				long long Data = VecToInt.VecToInt(bit_array, N);
				test_radix2 << "Data_index = ";
                test_radix2 << "( " ;					
				for(int k = 0; k < radix ; k++ ){
					test_radix2 << Data + k*(1<<(bit_width-bit_width_s-bit_width_s*s)) <<" ";	
				}
				test_radix2 << ") " ;
				test_radix2 << ", (w^" << 0 << ", w^" << length << ")" <<std::endl;
				//-----------------------------------------
				
				if(bn_tmp == 0){
					bn0_bc_tmp = BC_tmp;
					DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp]<<"\n";
					DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp]<< ", factor_t = " << factor_t << ", w^" << length << "\n";
					MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],factor_t,p);
					Radix2_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp]);
					DATARECORD << "---after BU compute---" << std::endl;
				    DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp]<<"\n";
					DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp]<<"\n";
					bn0_ma_reg = ma_tmp;

				}
			    else {
					DATARECORD << "bn_tmp = " << bn_tmp << ", ma_tmp = " << ma_tmp << std::endl;
					bn1_bc_tmp = BC_tmp;
					DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<<"\n";
					DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<< ", factor_t = " << factor_t << ", w^" << length << "\n";
					MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],factor_t,p);
					Radix2_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp]);
					DATARECORD << "---after BU compute---" << std::endl;
				    DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<<"\n";
					DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<<"\n";
					bn1_ma_reg = ma_tmp;
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


void DIT_NTTSPMB::test_radix4(std::vector<ZZ> &A){
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
	
	std::ofstream DATARECORD("./test_R4_SPMB.txt");
	std::ofstream spmb_radix4("./SPMB_tw/test_radix4.txt");
	//-------------Bit Reverse--------------
	BitOperate RR, IntToVec, VecToInt, BR, num_complement;
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
	int tw_degree = pow(radix, Stage-1); // siang
	std::cout << "init load over! \n";
	//need modify to mult by twiddle factor
	for(int s = 0; s < Stage; s++){
		int s_complement = num_complement.number_complement(s, Stage);
		if(s == 0){
			//factor = PowerMod(W, (pow(radix, s_complement)), p);
		}
		else {
			//factor = PowerMod(W, (pow(radix, s_complement)), p);
			tw_degree = tw_degree / radix;
		}
		DATARECORD <<"---------------------------------\n";
		DATARECORD <<"Now Stage: "<< s <<"\n";
		

		spmb_radix4 <<"Now Stage: "<< s <<"\n";
		spmb_radix4 <<"twiddle factor : "<< W <<"\n";
	    spmb_radix4 << "********\n";
		for(int i = 0 ;i < group;i++){
	    	DATARECORD << "********\n";
			bn0_bc_tmp  = 0;
			bn1_bc_tmp  = 0;
			spmb_radix4 <<"--------------i = " << i << "----------------\n";
			for(int j = 0;j < radix;j++){
				DATARECORD <<"twiddle factor : "<< W <<"\n";
				DATARECORD <<"p : "<< p <<"\n";
				gray_i  = Gray(i,group);
				BC_tmp  = j * group + gray_i;
				DATARECORD << "i: " << i <<"\n";
				DATARECORD << "BC_tmp: " << BC_tmp <<"\n";
				RR_R4(BC_tmp,s,BC);
				int shift_val = log2(radix) * s;
				int mod_group = pow(radix, s);
				int BC_tmp_mod = (BC_tmp % mod_group);
				DATARECORD << "BC_tmp_mod = " << BC_tmp_mod << ", A_B0R1[ma_tmp] = " << ma_tmp
							<< ", rotate = " << BR.left_rotate(BC_tmp_mod, 2, 64)
							<< ", base = " << (pow(radix, s_complement)) << endl;
				if(s == 1){
					length =  BC_tmp_mod * (pow(radix, s_complement));
				}else if (s == 2){
					length =  BR.left_rotate(BC_tmp_mod, s, pow(radix, s)) * (pow(radix, s_complement));
				}else if(s == 3){
					length =  BR.left_rotate(BC_tmp_mod, 4, 64) * (pow(radix, s_complement));
				}
				else if (s == 0){
					length = 0;
				}
				
				DATARECORD << "length: " <<  length <<"\n";
				PowerMod(factor_t,W,length,p);
				DATARECORD << "factor_t: "<<factor_t<<"\n";
				AGU_R4(BC,bn_tmp,ma_tmp);
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
				
				
				int debug = 0;
				if(bn_tmp == 0){
					if(j < 2)bn0_bc_tmp = BC_tmp;
					PowerMod(factor_2t,factor_t,2,p);
					PowerMod(factor_3t,factor_t,3,p);
					
					if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],factor_t,p);
					if(!debug) MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],factor_2t,p);
					if(!debug) MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],factor_3t,p);
					DATARECORD <<"Before butterfly unit operation! \n";
					//DATARECORD << "factor_0t = " << 1 << endl;
					//DATARECORD << "factor_1t = " << factor_t << endl;
					//DATARECORD << "factor_2t = " << factor_2t << endl;
					//DATARECORD << "factor_3t = " << factor_3t << endl;
					DATARECORD <<" A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp] << ", factor_0 : w^" << 0<<"\n";
				    DATARECORD <<" A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp] << ", facotr_1 : w^" << length  <<"\n";
					DATARECORD <<" A_B0R2["<<ma_tmp<<"]: "<<A_B0R2[ma_tmp] << ", facotr_2 : w^" << length*2 <<"\n";
					DATARECORD <<" A_B0R3["<<ma_tmp<<"]: "<<A_B0R3[ma_tmp] << ", facotr_3 : w^" << length*3 <<"\n";
					if(!debug) Radix4_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp],A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
					DATARECORD <<" -------------------" << std::endl;
					DATARECORD <<" A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp] << endl;
					DATARECORD <<" A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp] << endl;
					DATARECORD <<" A_B0R2["<<ma_tmp<<"]: "<<A_B0R2[ma_tmp] << endl;
					DATARECORD <<" A_B0R3["<<ma_tmp<<"]: "<<A_B0R3[ma_tmp] << endl;
				
			
					DATARECORD <<"--------------------------------------------------------------------\n";
					if(j <  2)bn0_ma_reg1 = ma_tmp;
					if(j >= 2)bn0_ma_reg2 = ma_tmp;
				}
			    else {
					if(j < 2)bn1_bc_tmp = BC_tmp;
					PowerMod(factor_2t,factor_t,2,p);
					PowerMod(factor_3t,factor_t,3,p);
					
					if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],factor_t,p);
					if(!debug) MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],factor_2t,p);
					if(!debug) MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],factor_3t,p);
					DATARECORD <<"Before butterfly unit operation! \n";
					//DATARECORD << "factor_0t = " << 1 << endl;
					//DATARECORD << "factor_1t = " << factor_t << endl;
					//DATARECORD << "factor_2t = " << factor_2t << endl;
					//DATARECORD << "factor_3t = " << factor_3t << endl;
					DATARECORD <<" A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp] << ", factor_0 : w^" << 0<<"\n";
					DATARECORD <<" A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp] << ", facotr_1 : w^" << length  <<"\n";
					DATARECORD <<" A_B1R2["<<ma_tmp<<"]: "<<A_B1R2[ma_tmp] << ", facotr_2 : w^" << length*2 <<"\n";
					DATARECORD <<" A_B1R3["<<ma_tmp<<"]: "<<A_B1R3[ma_tmp] << ", facotr_3 : w^" << length*3 <<"\n";
					if(!debug) Radix4_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp],A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
					DATARECORD <<" -------------------" << std::endl;
					DATARECORD <<" A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp] << endl;
					DATARECORD <<" A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp] << endl;
					DATARECORD <<" A_B1R2["<<ma_tmp<<"]: "<<A_B1R2[ma_tmp] << endl;
					DATARECORD <<" A_B1R3["<<ma_tmp<<"]: "<<A_B1R3[ma_tmp] << endl;

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