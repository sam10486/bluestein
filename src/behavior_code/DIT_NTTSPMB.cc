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
	std::ofstream DIT_DTFAG_golden_st0("./SPMB_tw/DIT_DTFAG_golden_st0_16.txt");
	std::ofstream DIT_DTFAG_golden_st1("./SPMB_tw/DIT_DTFAG_golden_st1_16.txt");
	std::ofstream DIT_DTFAG_golden_st2("./SPMB_tw/DIT_DTFAG_golden_st2_16.txt");
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
	//----------------------------------------

	
    Stage = (unsigned long)ceil(log2(N));
	BC_WIDTH  = (int)ceil(log2(N/2));
	offset    =  (int)N /radix;
	word_size =  (int)N / (2 * radix);
	group     =  (int)N / (radix * radix);
	bit_array_tmp.resize(BC_WIDTH);
	//DIT_DATARECORD <<"Stage: "<< Stage<<"\n";

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
    ZZ  factor_tmp;   
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
                //cout << "ma_tmp = " << ma_tmp << ", BC = " << BC << endl;
				A_B0R0[ma_tmp] = A[BC];
				A_B0R1[ma_tmp] = A[BC + offset];
				std::cout <<"A_B0R0["<<ma_tmp<<"]"<<A[BC]<<"\n";
				std::cout <<"A_B0R1["<<ma_tmp<<"]"<<A[BC+ offset]<<"\n";
			}else {
                //cout << "ma_tmp = " << ma_tmp << ", BC = " << BC << endl;
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
    factor_tmp = W;
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
            PowerMod(factor, factor_tmp, pow(2,s_complement), p);
        }
		else {
			PowerMod(factor, factor_tmp, pow(2,s_complement), p);
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
				//cout << "stage = " << s << ", DTFAG_t: " << DTFAG_t << ", DTFAG_i: " << DTFAG_i << ", DTFAG_j: " << DTFAG_j << endl;
				DTFAG.DTFAG_SPMB_DIT(s, 
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
				int debug = 0;
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

					//----------------DTFAG golden pattern------------------
					switch(s){
						case 0:
							DIT_DTFAG_golden_st0 << 1 		  << " \n ";
							DIT_DTFAG_golden_st0 << factor_t  << " \n ";
							break;
						case 1:
							DIT_DTFAG_golden_st1 << 1 		  << " \n ";
							DIT_DTFAG_golden_st1 << factor_t  << " \n ";
							break;
						case 2:
							DIT_DTFAG_golden_st2 << 1 		  << " \n ";
							DIT_DTFAG_golden_st2 << factor_t  << " \n ";
							break;
						default:
							break;
					}
					//------------------------------------------------------
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
					//----------------DTFAG golden pattern------------------
					switch(s){
						case 0:
							DIT_DTFAG_golden_st0 << 1 		  << " \n ";
							DIT_DTFAG_golden_st0 << factor_t  << " \n ";
							break;
						case 1:
							DIT_DTFAG_golden_st1 << 1 		  << " \n ";
							DIT_DTFAG_golden_st1 << factor_t  << " \n ";
							break;
						case 2:
							DIT_DTFAG_golden_st2 << 1 		  << " \n ";
							DIT_DTFAG_golden_st2 << factor_t  << " \n ";
							break;
						default:
							break;
					}
					//------------------------------------------------------
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