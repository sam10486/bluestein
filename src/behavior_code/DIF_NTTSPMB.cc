#include <iostream>
#include <fstream>
#include <NTL/ZZ.h>
#include <time.h>

#include "DIF_NTTSPMB.h"
#include "SPMB.h"

#include <vector>
#include <algorithm>
#include "BitOperate.h"
#include "DTFAG.h"

using namespace std;

void DIF_NTTSPMB::DIF_NTT_radix2(std::vector<ZZ> &A){
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

	std::ofstream DIF_DATARECORD("./DIF_NTT_R2_SPMB.txt");
	std::ofstream DIF_spmb_radix2("./SPMB_tw/DIF_spmb_radix2.txt");
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
	////----------------------------------------

	
    Stage = (unsigned long)ceil(log2(N));
	BC_WIDTH  = (int)ceil(log2(N/2));
	offset    =  (int)N /radix;
	word_size =  (int)N / (2 * radix);
	group     =  (int)N / (radix * radix);
	bit_array_tmp.resize(BC_WIDTH);
	DIF_DATARECORD <<"Stage: "<< Stage<<"\n";

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
		DIF_DATARECORD <<"Stage: "<< s<<"\n";
		DIF_DATARECORD << "factor "<< factor <<"\n";
		DIF_DATARECORD << "****************************\n";

		DIF_spmb_radix2 <<"Now Stage: "<< s <<"\n";
		DIF_spmb_radix2 <<"twiddle factor : "<< factor <<"\n";
		for(int i = 0 ;i < group;i++){
			DIF_spmb_radix2 << "----------------i =" << i << " ----------------" << std::endl;
			DIF_DATARECORD << "----------------i =" << i << " ----------------" << std::endl;
			for(int j = 0;j < radix;j++){
				DIF_DATARECORD << "---j =" << j << " ---" << std::endl;
				DIF_DATARECORD <<"twiddle factor : "<< factor <<"\n";
				DIF_DATARECORD <<"p : "<< p <<"\n";
	    		DIF_DATARECORD << "********\n";
				gray_i  = Gray(i,group);
				BC_tmp  = j * group + gray_i;
				DIF_DATARECORD << "BC_tmp = " << BC_tmp << ", BC = " << BC << ", s = " << s << endl;
				RR_R2(BC_tmp,s,BC);
				//std::cout << "BC: " << BC <<"\n";
				length = BC_tmp >> s;
				DIF_DATARECORD << "length: " <<  length <<"\n";
				PowerMod(factor_t,factor,length,p);
				DIF_DATARECORD << "factor_t: "<<factor_t<<"\n";

				AGU_R2(BC,bn_tmp,ma_tmp);

				//-----------compute data idx-------------
				DIF_spmb_radix2 << "BC_tmp = " << BC_tmp << std::endl;
				
				IntToVec.IntToVec(BC_tmp, N, bit_array);
				rotate(bit_array.begin(), bit_array.begin()+ bit_width_s*s , bit_array.end());
				long long Data = VecToInt.VecToInt(bit_array, N);
				DIF_spmb_radix2 << "Data_index = ";
                DIF_spmb_radix2 << "( " ;					
				for(int k = 0; k < radix ; k++ ){
					DIF_spmb_radix2 << Data + k*(1<<(bit_width-bit_width_s-bit_width_s*s)) <<" ";	
				}
				DIF_spmb_radix2 << ") " ;
				DIF_spmb_radix2 << ", (w^" << 0 << ", w^" << tw_degree * length << ")" <<std::endl;
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

				if(bn_tmp == 0){
					DIF_DATARECORD << "bn_tmp = " << bn_tmp << ", ma_tmp = " << ma_tmp << std::endl;
					bn0_bc_tmp = BC_tmp;
					
                    switch(s){
                        case 0:
                            DIF_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp] << ", st0_Tw[0] = " << st0_Tw[0] << endl;
					        DIF_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp] << ", st0_Tw[1] = " << st0_Tw[1] << endl;
					        if(!debug) Radix2_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp]);
                            if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],st0_Tw[0],p);
					        if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],st0_Tw[1],p);
                            break;
                        case 1:
                            DIF_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp] << ", st1_Tw[0] = " << st1_Tw[0] << endl;
					        DIF_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp] << ", st1_Tw[1] = " << st1_Tw[1] << endl;
					        if(!debug) Radix2_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp]);
                            if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],st1_Tw[0],p);
					        if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],st1_Tw[1],p);
                            break;
                        case 2:
                            DIF_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp] << ", st2_Tw[0] = " << st2_Tw[0] << endl;
					        DIF_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp] << ", st2_Tw[1] = " << st2_Tw[1] << endl;
					        if(!debug) Radix2_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp]);
                            if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],st2_Tw[0],p);
					        if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],st2_Tw[1],p);
                            break;
                        case 3:
                            DIF_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp] << ", st3_Tw[0] = " << 1 << endl;
					        DIF_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp] << ", st3_Tw[1] = " << 1 << endl;
					        if(!debug) Radix2_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp]);
                            if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],1,p);
					        if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],1,p);			
                            break;
                    }
                    
					DIF_DATARECORD << "---after BU compute---" << std::endl;
				    DIF_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp]<<"\n";
					DIF_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp]<<"\n";
					bn0_ma_reg = ma_tmp;

				}
			    else {
					DIF_DATARECORD << "bn_tmp = " << bn_tmp << ", ma_tmp = " << ma_tmp << std::endl;
					bn1_bc_tmp = BC_tmp;
                    switch(s){
                        case 0:
                            DIF_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<< ", st0_Tw[0] = " << st0_Tw[0] << endl;
					        DIF_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<< ", st0_Tw[1] = " << st0_Tw[1] << endl;
					        if(!debug) Radix2_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp]);
                            if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],st0_Tw[0],p);
					        if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],st0_Tw[1],p);    
                            break;
                        case 1:
                            DIF_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<< ", st1_Tw[0] = " << st1_Tw[0] << endl;
					        DIF_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<< ", st1_Tw[1] = " << st1_Tw[1] << endl;
					        if(!debug) Radix2_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp]);
                            if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],st1_Tw[0],p);
					        if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],st1_Tw[1],p);
                            break;
                        case 2:
                            DIF_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<< ", st2_Tw[0] = " << st2_Tw[0] << endl;
					        DIF_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<< ", st2_Tw[1] = " << st2_Tw[1] << endl;
					        if(!debug) Radix2_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp]);
                            if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],st2_Tw[0],p);
					        if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],st2_Tw[1],p);
                            break;
                        case 3:
                            DIF_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<< ", st3_Tw[0] = " << 1 << endl;
					        DIF_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<< ", st3_Tw[1] = " << 1 << endl;
					        if(!debug) Radix2_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp]);
                            if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],1,p);
					        if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],1,p);
                            break;
                    }
					DIF_DATARECORD << "---after BU compute---" << std::endl;
				    DIF_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<<"\n";
                    DIF_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<<"\n";
					bn1_ma_reg = ma_tmp;
				}
				DIF_DATARECORD <<"--------------------------------------------------------------------\n";
			}
		    //data relocation
		    if(s < Stage-1){
		    	if(bn1_bc_tmp > bn0_bc_tmp){
		    		DIF_DATARECORD << "bn0_bc_tmp = " << bn0_bc_tmp << ", bn1_bc_tmp = " << bn1_bc_tmp << std::endl;
		    		data_tmp = A_B0R1[bn0_ma_reg];
		    		A_B0R1[bn0_ma_reg] = A_B1R0[bn1_ma_reg];
		    		A_B1R0[bn1_ma_reg] = data_tmp;
		    	}else {
		    		DIF_DATARECORD << "bn0_bc_tmp = " << bn0_bc_tmp << ", bn1_bc_tmp = " << bn1_bc_tmp << std::endl;
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

void DIF_NTTSPMB::DIF_NTT_radix4(std::vector<ZZ> &A){
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
	
	std::ofstream DIF_DATARECORD("./DIF_NTT_R4_SPMB.txt");
	std::ofstream DIF_spmb_radix4("./SPMB_tw/DIF_spmb_radix4.txt");
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
	////----------------------------------------
	
    Stage = (unsigned long)ceil(log2(N));
	Stage = (unsigned long)Stage / 2;
	BC_WIDTH  = (int)ceil(log2(N/4));
	offset    =  (int)N /radix;
	word_size =  (int)N / (2 * radix);
	group     =  (int)N / (radix * radix);
	bit_array_tmp.resize(BC_WIDTH);
	DIF_DATARECORD <<"Stage: "<< Stage<<"\n";
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
		DIF_DATARECORD <<"---------------------------------\n";
		DIF_DATARECORD <<"Now Stage: "<< s <<"\n";
		

		DIF_spmb_radix4 <<"Now Stage: "<< s <<"\n";
		DIF_spmb_radix4 <<"twiddle factor : "<< factor <<"\n";
	    DIF_spmb_radix4 << "********\n";
		for(int i = 0 ;i < group;i++){
			DIF_DATARECORD <<"twiddle factor : "<< factor <<"\n";
			DIF_DATARECORD <<"p : "<< p <<"\n";
	    	DIF_DATARECORD << "********\n";
			bn0_bc_tmp  = 0;
			bn1_bc_tmp  = 0;
			DIF_spmb_radix4 <<"--------------i = " << i << "----------------\n";
			for(int j = 0;j < radix;j++){	
				gray_i  = Gray(i,group);
				BC_tmp  = j * group + gray_i;
				DIF_DATARECORD << "i: " << i <<"\n";
				DIF_DATARECORD << "gray_i: " << gray_i <<"\n";
				DIF_DATARECORD << "BC_tmp: " << BC_tmp <<"\n";
				RR_R4(BC_tmp,s,BC);
				DIF_DATARECORD << "BC: " << BC <<"\n";
				length = BC_tmp >> (2*s);
				DIF_DATARECORD << "length: " <<  length <<"\n";

				//-----------compute data idx-------------
				DIF_spmb_radix4 << "BC_tmp = " << BC_tmp << std::endl;
				
				IntToVec.IntToVec(BC_tmp, N, bit_array);
				rotate(bit_array.begin(), bit_array.begin()+ bit_width_s*s , bit_array.end());
				long long Data = VecToInt.VecToInt(bit_array, N);
				DIF_spmb_radix4 << "Data_index = ";
                DIF_spmb_radix4 << "( " ;		
				int spmb_radix4_arr[radix];		
				for(long long k = 0; k < radix ; k++ ){
					long long BR_tmp = BR.BitReserve(k, log2(radix));
					spmb_radix4_arr[BR_tmp] = Data + k*(1<<(bit_width-bit_width_s-bit_width_s*s));	
				}
				for(int k=0; k<radix; k++) {
					DIF_spmb_radix4 << spmb_radix4_arr[k] << ", ";
				}
				DIF_spmb_radix4 << ") " ;
				DIF_spmb_radix4 << ", (w^" << 0 << ", w^" << tw_degree * length 
				<< ", w^"  << tw_degree * length * 2 << ", w^" << tw_degree * length * 3 << ")" << std::endl;
				//-----------------------------------------
				
				PowerMod(factor_t,factor,length,p);
				DIF_DATARECORD << "factor_t: "<<factor_t<<"\n";
				AGU_R4(BC,bn_tmp,ma_tmp);
				DIF_DATARECORD << "bn_tmp: "<<bn_tmp<<"\n";
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
                
                if(bn_tmp == 0){
					if(j < 2)bn0_bc_tmp = BC_tmp;
					PowerMod(factor_2t,factor_t,2,p);
					PowerMod(factor_3t,factor_t,3,p);
                    switch(s){
                        case 0:
                            if(debug) DIF_DATARECORD <<"Before butterfly unit operation! \n";
                            if(debug) DIF_DATARECORD <<" A_B0R0["<<ma_tmp<<"]: " <<A_B0R0[ma_tmp] << ", st0_Tw[0] = " << st0_Tw[0] << endl;
					        if(debug) DIF_DATARECORD <<" A_B0R1["<<ma_tmp<<"]: " <<A_B0R1[ma_tmp] << ", st0_Tw[1] = " << st0_Tw[1] << endl;
					        if(debug) DIF_DATARECORD <<" A_B0R2["<<ma_tmp<<"]: " <<A_B0R2[ma_tmp] << ", st0_Tw[2] = " << st0_Tw[2] << endl;
					        if(debug) DIF_DATARECORD <<" A_B0R3["<<ma_tmp<<"]: " <<A_B0R3[ma_tmp] << ", st0_Tw[3] = " << st0_Tw[3] << endl;   
					        if(!debug) Radix4_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp],A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
					        DIF_DATARECORD <<" -------------------" << std::endl;
                            if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],st0_Tw[0],p);
                            if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],st0_Tw[1],p);
					        if(!debug) MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],st0_Tw[2],p);
					        if(!debug) MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],st0_Tw[3],p);
                            break;
                        case 1:
                            if(debug) DIF_DATARECORD <<"Before butterfly unit operation! \n";
					        if(debug) DIF_DATARECORD <<" A_B0R0["<<ma_tmp<<"]: " <<A_B0R0[ma_tmp] << ", st1_Tw[0] = " << st1_Tw[0] << endl;
					        if(debug) DIF_DATARECORD <<" A_B0R1["<<ma_tmp<<"]: " <<A_B0R1[ma_tmp] << ", st1_Tw[1] = " << st1_Tw[1] << endl;
					        if(debug) DIF_DATARECORD <<" A_B0R2["<<ma_tmp<<"]: " <<A_B0R2[ma_tmp] << ", st1_Tw[2] = " << st1_Tw[2] << endl;
					        if(debug) DIF_DATARECORD <<" A_B0R3["<<ma_tmp<<"]: " <<A_B0R3[ma_tmp] << ", st1_Tw[3] = " << st1_Tw[3] << endl;
					        if(!debug) Radix4_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp],A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
					        DIF_DATARECORD <<" -------------------" << std::endl;
                            if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],st1_Tw[0],p);
                            if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],st1_Tw[1],p);
					        if(!debug) MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],st1_Tw[2],p);
					        if(!debug) MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],st1_Tw[3],p);
                            break;
                        case 2:
                            if(debug) DIF_DATARECORD <<"Before butterfly unit operation! \n";
					        if(debug) DIF_DATARECORD <<" A_B0R0["<<ma_tmp<<"]: " <<A_B0R0[ma_tmp] << ", st2_Tw[0] = " << st2_Tw[0] << endl;
					        if(debug) DIF_DATARECORD <<" A_B0R1["<<ma_tmp<<"]: " <<A_B0R1[ma_tmp] << ", st2_Tw[1] = " << st2_Tw[1] << endl;
					        if(debug) DIF_DATARECORD <<" A_B0R2["<<ma_tmp<<"]: " <<A_B0R2[ma_tmp] << ", st2_Tw[2] = " << st2_Tw[2] << endl;
					        if(debug) DIF_DATARECORD <<" A_B0R3["<<ma_tmp<<"]: " <<A_B0R3[ma_tmp] << ", st2_Tw[3] = " << st2_Tw[3] << endl;
					        if(!debug) Radix4_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp],A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
					        DIF_DATARECORD <<" -------------------" << std::endl;
                            if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],st2_Tw[0],p);
                            if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],st2_Tw[1],p);
					        if(!debug) MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],st2_Tw[2],p);
					        if(!debug) MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],st2_Tw[3],p);
                            break;
                        case 3:
                            if(debug) DIF_DATARECORD <<"Before butterfly unit operation! \n";
					        if(debug) DIF_DATARECORD <<" A_B0R0["<<ma_tmp<<"]: " <<A_B0R0[ma_tmp] << ", st3_Tw[0] = " << 1 << endl;
					        if(debug) DIF_DATARECORD <<" A_B0R1["<<ma_tmp<<"]: " <<A_B0R1[ma_tmp] << ", st3_Tw[1] = " << 1 << endl;
					        if(debug) DIF_DATARECORD <<" A_B0R2["<<ma_tmp<<"]: " <<A_B0R2[ma_tmp] << ", st3_Tw[2] = " << 1 << endl;
					        if(debug) DIF_DATARECORD <<" A_B0R3["<<ma_tmp<<"]: " <<A_B0R3[ma_tmp] << ", st3_Tw[3] = " << 1 << endl;
					        if(!debug) Radix4_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp],A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
					        DIF_DATARECORD <<" -------------------" << std::endl;
                            if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],1,p);
                            if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],1,p);
					        if(!debug) MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],1,p);
					        if(!debug) MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],1,p);
                            break;
                    }
					DIF_DATARECORD <<" ------after BU compute and Mul-------" << std::endl;
					DIF_DATARECORD <<" A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp] << endl;
					DIF_DATARECORD <<" A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp] << endl;
					DIF_DATARECORD <<" A_B0R2["<<ma_tmp<<"]: "<<A_B0R2[ma_tmp] << endl;
					DIF_DATARECORD <<" A_B0R3["<<ma_tmp<<"]: "<<A_B0R3[ma_tmp] << endl;
				
					
					DIF_DATARECORD <<"--------------------------------------------------------------------\n";
					if(j <  2)bn0_ma_reg1 = ma_tmp;
					if(j >= 2)bn0_ma_reg2 = ma_tmp;
				}
			    else {
					if(j < 2)bn1_bc_tmp = BC_tmp;
					PowerMod(factor_2t,factor_t,2,p);
					PowerMod(factor_3t,factor_t,3,p);
                    switch(s){
                        case 0:
                            if(debug) DIF_DATARECORD <<" A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<< ", st0_Tw[0] = " << st0_Tw[0] << endl;
					        if(debug) DIF_DATARECORD <<" A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<< ", st0_Tw[1] = " << st0_Tw[1] << endl;
					        if(debug) DIF_DATARECORD <<" A_B1R2["<<ma_tmp<<"]: "<<A_B1R2[ma_tmp]<< ", st0_Tw[2] = " << st0_Tw[2] << endl;
					        if(debug) DIF_DATARECORD <<" A_B1R3["<<ma_tmp<<"]: "<<A_B1R3[ma_tmp]<< ", st0_Tw[3] = " << st0_Tw[3] << endl;
                            if(!debug) Radix4_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp],A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
					        DIF_DATARECORD <<" -------------------" << std::endl;
                            if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],st0_Tw[0],p);
					        if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],st0_Tw[1],p);
					        if(!debug) MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],st0_Tw[2],p);
					        if(!debug) MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],st0_Tw[3],p);
                            break;
                        case 1:
                            if(debug) DIF_DATARECORD <<" A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<< ", st1_Tw[0] = " << st1_Tw[0] << endl;
					        if(debug) DIF_DATARECORD <<" A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<< ", st1_Tw[1] = " << st1_Tw[1] << endl;
					        if(debug) DIF_DATARECORD <<" A_B1R2["<<ma_tmp<<"]: "<<A_B1R2[ma_tmp]<< ", st1_Tw[2] = " << st1_Tw[2] << endl;
					        if(debug) DIF_DATARECORD <<" A_B1R3["<<ma_tmp<<"]: "<<A_B1R3[ma_tmp]<< ", st1_Tw[3] = " << st1_Tw[3] << endl;
                            if(!debug) Radix4_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp],A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
					        DIF_DATARECORD <<" -------------------" << std::endl;
                            if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],st1_Tw[0],p);
					        if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],st1_Tw[1],p);
					        if(!debug) MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],st1_Tw[2],p);
					        if(!debug) MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],st1_Tw[3],p);
                            break;
                        case 2:
                            if(debug) DIF_DATARECORD <<" A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<< ", st2_Tw[0] = " << st2_Tw[0] << endl;
					        if(debug) DIF_DATARECORD <<" A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<< ", st2_Tw[1] = " << st2_Tw[1] << endl;
					        if(debug) DIF_DATARECORD <<" A_B1R2["<<ma_tmp<<"]: "<<A_B1R2[ma_tmp]<< ", st2_Tw[2] = " << st2_Tw[2] << endl;
					        if(debug) DIF_DATARECORD <<" A_B1R3["<<ma_tmp<<"]: "<<A_B1R3[ma_tmp]<< ", st2_Tw[3] = " << st2_Tw[3] << endl;
                            if(!debug) Radix4_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp],A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
					        DIF_DATARECORD <<" -------------------" << std::endl;
                            if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],st2_Tw[0],p);
					        if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],st2_Tw[1],p);
					        if(!debug) MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],st2_Tw[2],p);
					        if(!debug) MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],st2_Tw[3],p);
                            break;
                        case 3:
                            if(debug) DIF_DATARECORD <<" A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<< ", st3_Tw[0] = " << 1 << endl;
					        if(debug) DIF_DATARECORD <<" A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<< ", st3_Tw[1] = " << 1 << endl;
					        if(debug) DIF_DATARECORD <<" A_B1R2["<<ma_tmp<<"]: "<<A_B1R2[ma_tmp]<< ", st3_Tw[2] = " << 1 << endl;
					        if(debug) DIF_DATARECORD <<" A_B1R3["<<ma_tmp<<"]: "<<A_B1R3[ma_tmp]<< ", st3_Tw[3] = " << 1 << endl;
                            if(!debug) Radix4_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp],A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
					        DIF_DATARECORD <<" -------------------" << std::endl;
                            if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],1,p);
					        if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],1,p);
					        if(!debug) MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],1,p);
					        if(!debug) MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],1,p);
                            break;
                    }
                    DIF_DATARECORD <<" ------after BU compute and Mul-------" << std::endl;
					DIF_DATARECORD <<" A_B1R0["<<ma_tmp<<"]: " << A_B1R0[ma_tmp] << endl;
					DIF_DATARECORD <<" A_B1R1["<<ma_tmp<<"]: " << A_B1R1[ma_tmp] << endl;
					DIF_DATARECORD <<" A_B1R2["<<ma_tmp<<"]: " << A_B1R2[ma_tmp] << endl;
					DIF_DATARECORD <<" A_B1R3["<<ma_tmp<<"]: " << A_B1R3[ma_tmp] << endl;
					
					DIF_DATARECORD <<"--------------------------------------------------------------------\n";
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

void DIF_NTTSPMB::DIF_NTT_radix16(std::vector<ZZ> &A){
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
	
	std::ofstream DIF_DATARECORD("./DIF_NTT_R16_SPMB.txt");

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
	////----------------------------------------
	
    Stage = (unsigned long)ceil(log2(N));
	Stage = (unsigned long)Stage / 4;
	BC_WIDTH  = (int)ceil(log2(N/16));
	offset    =  (int)N /radix;
	word_size =  (int)N / (2 * radix);
	group     =  (int)N / (radix * radix);
	bit_array_tmp.resize(BC_WIDTH);
	DIF_DATARECORD <<"Stage: "<< Stage<<"\n";

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
		DIF_DATARECORD <<"---------------------------------\n";
		DIF_DATARECORD <<"Now Stage: "<< s <<"\n";
		DIF_DATARECORD <<"twiddle factor : "<< factor <<"\n";

		std::cout << "twiddle factor : "<< factor <<"\n";
	    DIF_DATARECORD << "********\n";
		for(int i = 0 ;i < group;i++){
			bn0_bc_tmp  = 0;
			bn1_bc_tmp  = 0;	
			for(int j = 0;j < radix;j++){
				gray_i  = Gray(i,group);
				BC_tmp  = j * group + gray_i;
				DIF_DATARECORD << "i: " << i <<"\n";	
				DIF_DATARECORD << "gray_i: " << gray_i <<"\n";
				DIF_DATARECORD << "BC_tmp: " << BC_tmp <<"\n";
				RR_R16(BC_tmp,s,BC);
				DIF_DATARECORD << "BC: " << BC <<"\n";
				length = BC_tmp >> (4*s);
				DIF_DATARECORD << "length: " <<  length <<"\n";
				
				PowerMod(factor_t,factor,length,p);
				DIF_DATARECORD << "factor_t: "<<factor_t<<"\n";
				AGU_R16(BC,bn_tmp,ma_tmp);
				DIF_DATARECORD << "bn_tmp: "<<bn_tmp<<"\n";

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
                DIF_DATARECORD << "stage = " << s << ", DTFAG_i = " << DTFAG_i << ", DTFAG_t = " << DTFAG_t << ", DTFAG_j = " << DTFAG_j << endl;
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

				if(bn_tmp == 0){
					if(j < 2)bn0_bc_tmp = BC_tmp;
					switch(s){
                        case 0:
                            DIF_DATARECORD <<" A_B0R0["<<ma_tmp<<"]:  "<<A_B0R0[ma_tmp]  << ", st0_Tw[0] = " << st0_Tw[0] << endl;
				            DIF_DATARECORD <<" A_B0R1["<<ma_tmp<<"]:  "<<A_B0R1[ma_tmp]  << ", st0_Tw[1] = " << st0_Tw[1] << endl;
					        DIF_DATARECORD <<" A_B0R2["<<ma_tmp<<"]:  "<<A_B0R2[ma_tmp]  << ", st0_Tw[2] = " << st0_Tw[2] << endl;
					        DIF_DATARECORD <<" A_B0R3["<<ma_tmp<<"]:  "<<A_B0R3[ma_tmp]  << ", st0_Tw[3] = " << st0_Tw[3] << endl;
					        DIF_DATARECORD <<" A_B0R4["<<ma_tmp<<"]:  "<<A_B0R4[ma_tmp]  << ", st0_Tw[4] = " << st0_Tw[4] << endl;
					        DIF_DATARECORD <<" A_B0R5["<<ma_tmp<<"]:  "<<A_B0R5[ma_tmp]  << ", st0_Tw[5] = " << st0_Tw[5] << endl;
					        DIF_DATARECORD <<" A_B0R6["<<ma_tmp<<"]:  "<<A_B0R6[ma_tmp]  << ", st0_Tw[6] = " << st0_Tw[6] << endl;
					        DIF_DATARECORD <<" A_B0R7["<<ma_tmp<<"]:  "<<A_B0R7[ma_tmp]  << ", st0_Tw[7] = " << st0_Tw[7] << endl;
					        DIF_DATARECORD <<" A_B0R8["<<ma_tmp<<"]:  "<<A_B0R8[ma_tmp]  << ", st0_Tw[8] = " << st0_Tw[8] << endl;
					        DIF_DATARECORD <<" A_B0R9["<<ma_tmp<<"]:  "<<A_B0R9[ma_tmp]  << ", st0_Tw[9] = " << st0_Tw[9] << endl;
					        DIF_DATARECORD <<" A_B0R10["<<ma_tmp<<"]: "<<A_B0R10[ma_tmp] << ", st0_Tw[10] = " << st0_Tw[10] << endl;   
					        DIF_DATARECORD <<" A_B0R11["<<ma_tmp<<"]: "<<A_B0R11[ma_tmp] << ", st0_Tw[11] = " << st0_Tw[11] << endl;   
					        DIF_DATARECORD <<" A_B0R12["<<ma_tmp<<"]: "<<A_B0R12[ma_tmp] << ", st0_Tw[12] = " << st0_Tw[12] << endl;   
					        DIF_DATARECORD <<" A_B0R13["<<ma_tmp<<"]: "<<A_B0R13[ma_tmp] << ", st0_Tw[13] = " << st0_Tw[13] << endl;   
					        DIF_DATARECORD <<" A_B0R14["<<ma_tmp<<"]: "<<A_B0R14[ma_tmp] << ", st0_Tw[14] = " << st0_Tw[14] << endl;   
					        DIF_DATARECORD <<" A_B0R15["<<ma_tmp<<"]: "<<A_B0R15[ma_tmp] << ", st0_Tw[15] = " << st0_Tw[15] << endl;   
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
                            DIF_DATARECORD <<" A_B0R0["<<ma_tmp<<"]:  "<<A_B0R0[ma_tmp]  << ", st1_Tw[0] = " << st1_Tw[0] << endl;
				            DIF_DATARECORD <<" A_B0R1["<<ma_tmp<<"]:  "<<A_B0R1[ma_tmp]  << ", st1_Tw[1] = " << st1_Tw[1] << endl;
					        DIF_DATARECORD <<" A_B0R2["<<ma_tmp<<"]:  "<<A_B0R2[ma_tmp]  << ", st1_Tw[2] = " << st1_Tw[2] << endl;
					        DIF_DATARECORD <<" A_B0R3["<<ma_tmp<<"]:  "<<A_B0R3[ma_tmp]  << ", st1_Tw[3] = " << st1_Tw[3] << endl;
					        DIF_DATARECORD <<" A_B0R4["<<ma_tmp<<"]:  "<<A_B0R4[ma_tmp]  << ", st1_Tw[4] = " << st1_Tw[4] << endl;
					        DIF_DATARECORD <<" A_B0R5["<<ma_tmp<<"]:  "<<A_B0R5[ma_tmp]  << ", st1_Tw[5] = " << st1_Tw[5] << endl;
					        DIF_DATARECORD <<" A_B0R6["<<ma_tmp<<"]:  "<<A_B0R6[ma_tmp]  << ", st1_Tw[6] = " << st1_Tw[6] << endl;
					        DIF_DATARECORD <<" A_B0R7["<<ma_tmp<<"]:  "<<A_B0R7[ma_tmp]  << ", st1_Tw[7] = " << st1_Tw[7] << endl;
					        DIF_DATARECORD <<" A_B0R8["<<ma_tmp<<"]:  "<<A_B0R8[ma_tmp]  << ", st1_Tw[8] = " << st1_Tw[8] << endl;
					        DIF_DATARECORD <<" A_B0R9["<<ma_tmp<<"]:  "<<A_B0R9[ma_tmp]  << ", st1_Tw[9] = " << st1_Tw[9] << endl;
					        DIF_DATARECORD <<" A_B0R10["<<ma_tmp<<"]: "<<A_B0R10[ma_tmp] << ", st1_Tw[10] = " << st1_Tw[10] << endl;   
					        DIF_DATARECORD <<" A_B0R11["<<ma_tmp<<"]: "<<A_B0R11[ma_tmp] << ", st1_Tw[11] = " << st1_Tw[11] << endl;   
					        DIF_DATARECORD <<" A_B0R12["<<ma_tmp<<"]: "<<A_B0R12[ma_tmp] << ", st1_Tw[12] = " << st1_Tw[12] << endl;   
					        DIF_DATARECORD <<" A_B0R13["<<ma_tmp<<"]: "<<A_B0R13[ma_tmp] << ", st1_Tw[13] = " << st1_Tw[13] << endl;   
					        DIF_DATARECORD <<" A_B0R14["<<ma_tmp<<"]: "<<A_B0R14[ma_tmp] << ", st1_Tw[14] = " << st1_Tw[14] << endl;   
					        DIF_DATARECORD <<" A_B0R15["<<ma_tmp<<"]: "<<A_B0R15[ma_tmp] << ", st1_Tw[15] = " << st1_Tw[15] << endl;   
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
                            DIF_DATARECORD <<" A_B0R0["<<ma_tmp<<"]:  "<<A_B0R0[ma_tmp]  << ", st2_Tw[0] = " << st2_Tw[0] << endl;
				            DIF_DATARECORD <<" A_B0R1["<<ma_tmp<<"]:  "<<A_B0R1[ma_tmp]  << ", st2_Tw[1] = " << st2_Tw[1] << endl;
					        DIF_DATARECORD <<" A_B0R2["<<ma_tmp<<"]:  "<<A_B0R2[ma_tmp]  << ", st2_Tw[2] = " << st2_Tw[2] << endl;
					        DIF_DATARECORD <<" A_B0R3["<<ma_tmp<<"]:  "<<A_B0R3[ma_tmp]  << ", st2_Tw[3] = " << st2_Tw[3] << endl;
					        DIF_DATARECORD <<" A_B0R4["<<ma_tmp<<"]:  "<<A_B0R4[ma_tmp]  << ", st2_Tw[4] = " << st2_Tw[4] << endl;
					        DIF_DATARECORD <<" A_B0R5["<<ma_tmp<<"]:  "<<A_B0R5[ma_tmp]  << ", st2_Tw[5] = " << st2_Tw[5] << endl;
					        DIF_DATARECORD <<" A_B0R6["<<ma_tmp<<"]:  "<<A_B0R6[ma_tmp]  << ", st2_Tw[6] = " << st2_Tw[6] << endl;
					        DIF_DATARECORD <<" A_B0R7["<<ma_tmp<<"]:  "<<A_B0R7[ma_tmp]  << ", st2_Tw[7] = " << st2_Tw[7] << endl;
					        DIF_DATARECORD <<" A_B0R8["<<ma_tmp<<"]:  "<<A_B0R8[ma_tmp]  << ", st2_Tw[8] = " << st2_Tw[8] << endl;
					        DIF_DATARECORD <<" A_B0R9["<<ma_tmp<<"]:  "<<A_B0R9[ma_tmp]  << ", st2_Tw[9] = " << st2_Tw[9] << endl;
					        DIF_DATARECORD <<" A_B0R10["<<ma_tmp<<"]: "<<A_B0R10[ma_tmp] << ", st2_Tw[10] = " << st2_Tw[10] << endl;   
					        DIF_DATARECORD <<" A_B0R11["<<ma_tmp<<"]: "<<A_B0R11[ma_tmp] << ", st2_Tw[11] = " << st2_Tw[11] << endl;   
					        DIF_DATARECORD <<" A_B0R12["<<ma_tmp<<"]: "<<A_B0R12[ma_tmp] << ", st2_Tw[12] = " << st2_Tw[12] << endl;   
					        DIF_DATARECORD <<" A_B0R13["<<ma_tmp<<"]: "<<A_B0R13[ma_tmp] << ", st2_Tw[13] = " << st2_Tw[13] << endl;   
					        DIF_DATARECORD <<" A_B0R14["<<ma_tmp<<"]: "<<A_B0R14[ma_tmp] << ", st2_Tw[14] = " << st2_Tw[14] << endl;   
					        DIF_DATARECORD <<" A_B0R15["<<ma_tmp<<"]: "<<A_B0R15[ma_tmp] << ", st2_Tw[15] = " << st2_Tw[15] << endl;   
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
                            DIF_DATARECORD <<" A_B0R0["<<ma_tmp<<"]:  "<<A_B0R0[ma_tmp]  << ", st3_Tw[0] = " << 1 << endl;
				            DIF_DATARECORD <<" A_B0R1["<<ma_tmp<<"]:  "<<A_B0R1[ma_tmp]  << ", st3_Tw[1] = " << 1 << endl;
					        DIF_DATARECORD <<" A_B0R2["<<ma_tmp<<"]:  "<<A_B0R2[ma_tmp]  << ", st3_Tw[2] = " << 1 << endl;
					        DIF_DATARECORD <<" A_B0R3["<<ma_tmp<<"]:  "<<A_B0R3[ma_tmp]  << ", st3_Tw[3] = " << 1 << endl;
					        DIF_DATARECORD <<" A_B0R4["<<ma_tmp<<"]:  "<<A_B0R4[ma_tmp]  << ", st3_Tw[4] = " << 1 << endl;
					        DIF_DATARECORD <<" A_B0R5["<<ma_tmp<<"]:  "<<A_B0R5[ma_tmp]  << ", st3_Tw[5] = " << 1 << endl;
					        DIF_DATARECORD <<" A_B0R6["<<ma_tmp<<"]:  "<<A_B0R6[ma_tmp]  << ", st3_Tw[6] = " << 1 << endl;
					        DIF_DATARECORD <<" A_B0R7["<<ma_tmp<<"]:  "<<A_B0R7[ma_tmp]  << ", st3_Tw[7] = " << 1 << endl;
					        DIF_DATARECORD <<" A_B0R8["<<ma_tmp<<"]:  "<<A_B0R8[ma_tmp]  << ", st3_Tw[8] = " << 1 << endl;
					        DIF_DATARECORD <<" A_B0R9["<<ma_tmp<<"]:  "<<A_B0R9[ma_tmp]  << ", st3_Tw[9] = " << 1 << endl;
					        DIF_DATARECORD <<" A_B0R10["<<ma_tmp<<"]: "<<A_B0R10[ma_tmp] << ", st3_Tw[10] = " << 1 << endl;   
					        DIF_DATARECORD <<" A_B0R11["<<ma_tmp<<"]: "<<A_B0R11[ma_tmp] << ", st3_Tw[11] = " << 1 << endl;   
					        DIF_DATARECORD <<" A_B0R12["<<ma_tmp<<"]: "<<A_B0R12[ma_tmp] << ", st3_Tw[12] = " << 1 << endl;   
					        DIF_DATARECORD <<" A_B0R13["<<ma_tmp<<"]: "<<A_B0R13[ma_tmp] << ", st3_Tw[13] = " << 1 << endl;   
					        DIF_DATARECORD <<" A_B0R14["<<ma_tmp<<"]: "<<A_B0R14[ma_tmp] << ", st3_Tw[14] = " << 1 << endl;   
					        DIF_DATARECORD <<" A_B0R15["<<ma_tmp<<"]: "<<A_B0R15[ma_tmp] << ", st3_Tw[15] = " << 1 << endl;   
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
                    DIF_DATARECORD <<" ------after BU compute and Mul-------" << std::endl;
				    DIF_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp]<<"\n";
				    DIF_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp]<<"\n";
					DIF_DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<<A_B0R2[ma_tmp]<<"\n";
					DIF_DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<<A_B0R3[ma_tmp]<<"\n";
					DIF_DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<<A_B0R4[ma_tmp]<<"\n";
					DIF_DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<<A_B0R5[ma_tmp]<<"\n";
					DIF_DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<<A_B0R6[ma_tmp]<<"\n";
					DIF_DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<<A_B0R7[ma_tmp]<<"\n";
					DIF_DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<<A_B0R8[ma_tmp]<<"\n";
					DIF_DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<<A_B0R9[ma_tmp]<<"\n";
					DIF_DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<<A_B0R10[ma_tmp]<<"\n";
					DIF_DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<<A_B0R11[ma_tmp]<<"\n";
					DIF_DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<<A_B0R12[ma_tmp]<<"\n";
					DIF_DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<<A_B0R13[ma_tmp]<<"\n";
					DIF_DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<<A_B0R14[ma_tmp]<<"\n";
					DIF_DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<<A_B0R15[ma_tmp]<<"\n";
					DIF_DATARECORD <<"--------------------------------------------------------------------\n";
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
                            DIF_DATARECORD <<" A_B1R0["<<ma_tmp<<"]:  "<<A_B1R0[ma_tmp]  << ", st0_Tw[0] = " << st0_Tw[0] << endl;     
				            DIF_DATARECORD <<" A_B1R1["<<ma_tmp<<"]:  "<<A_B1R1[ma_tmp]  << ", st0_Tw[1] = " << st0_Tw[1] << endl;     
					        DIF_DATARECORD <<" A_B1R2["<<ma_tmp<<"]:  "<<A_B1R2[ma_tmp]  << ", st0_Tw[2] = " << st0_Tw[2] << endl;     
					        DIF_DATARECORD <<" A_B1R3["<<ma_tmp<<"]:  "<<A_B1R3[ma_tmp]  << ", st0_Tw[3] = " << st0_Tw[3] << endl;     
					        DIF_DATARECORD <<" A_B1R4["<<ma_tmp<<"]:  "<<A_B1R4[ma_tmp]  << ", st0_Tw[4] = " << st0_Tw[4] << endl;     
					        DIF_DATARECORD <<" A_B1R5["<<ma_tmp<<"]:  "<<A_B1R5[ma_tmp]  << ", st0_Tw[5] = " << st0_Tw[5] << endl;     
					        DIF_DATARECORD <<" A_B1R6["<<ma_tmp<<"]:  "<<A_B1R6[ma_tmp]  << ", st0_Tw[6] = " << st0_Tw[6] << endl;     
					        DIF_DATARECORD <<" A_B1R7["<<ma_tmp<<"]:  "<<A_B1R7[ma_tmp]  << ", st0_Tw[7] = " << st0_Tw[7] << endl;     
					        DIF_DATARECORD <<" A_B1R8["<<ma_tmp<<"]:  "<<A_B1R8[ma_tmp]  << ", st0_Tw[8] = " << st0_Tw[8] << endl;     
					        DIF_DATARECORD <<" A_B1R9["<<ma_tmp<<"]:  "<<A_B1R9[ma_tmp]  << ", st0_Tw[9] = " << st0_Tw[9] << endl;     
					        DIF_DATARECORD <<" A_B1R10["<<ma_tmp<<"]: "<<A_B1R10[ma_tmp] << ", st0_Tw[10] = " << st0_Tw[10] << endl;    
					        DIF_DATARECORD <<" A_B1R11["<<ma_tmp<<"]: "<<A_B1R11[ma_tmp] << ", st0_Tw[11] = " << st0_Tw[11] << endl;    
					        DIF_DATARECORD <<" A_B1R12["<<ma_tmp<<"]: "<<A_B1R12[ma_tmp] << ", st0_Tw[12] = " << st0_Tw[12] << endl;    
					        DIF_DATARECORD <<" A_B1R13["<<ma_tmp<<"]: "<<A_B1R13[ma_tmp] << ", st0_Tw[13] = " << st0_Tw[13] << endl;    
					        DIF_DATARECORD <<" A_B1R14["<<ma_tmp<<"]: "<<A_B1R14[ma_tmp] << ", st0_Tw[14] = " << st0_Tw[14] << endl;    
					        DIF_DATARECORD <<" A_B1R15["<<ma_tmp<<"]: "<<A_B1R15[ma_tmp] << ", st0_Tw[15] = " << st0_Tw[15] << endl;    
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
                            DIF_DATARECORD <<" A_B1R0["<<ma_tmp<<"]:  "<<A_B1R0[ma_tmp]  << ", st1_Tw[0] = " << st1_Tw[0] << endl;     
				            DIF_DATARECORD <<" A_B1R1["<<ma_tmp<<"]:  "<<A_B1R1[ma_tmp]  << ", st1_Tw[1] = " << st1_Tw[1] << endl;     
					        DIF_DATARECORD <<" A_B1R2["<<ma_tmp<<"]:  "<<A_B1R2[ma_tmp]  << ", st1_Tw[2] = " << st1_Tw[2] << endl;     
					        DIF_DATARECORD <<" A_B1R3["<<ma_tmp<<"]:  "<<A_B1R3[ma_tmp]  << ", st1_Tw[3] = " << st1_Tw[3] << endl;     
					        DIF_DATARECORD <<" A_B1R4["<<ma_tmp<<"]:  "<<A_B1R4[ma_tmp]  << ", st1_Tw[4] = " << st1_Tw[4] << endl;     
					        DIF_DATARECORD <<" A_B1R5["<<ma_tmp<<"]:  "<<A_B1R5[ma_tmp]  << ", st1_Tw[5] = " << st1_Tw[5] << endl;     
					        DIF_DATARECORD <<" A_B1R6["<<ma_tmp<<"]:  "<<A_B1R6[ma_tmp]  << ", st1_Tw[6] = " << st1_Tw[6] << endl;     
					        DIF_DATARECORD <<" A_B1R7["<<ma_tmp<<"]:  "<<A_B1R7[ma_tmp]  << ", st1_Tw[7] = " << st1_Tw[7] << endl;     
					        DIF_DATARECORD <<" A_B1R8["<<ma_tmp<<"]:  "<<A_B1R8[ma_tmp]  << ", st1_Tw[8] = " << st1_Tw[8] << endl;     
					        DIF_DATARECORD <<" A_B1R9["<<ma_tmp<<"]:  "<<A_B1R9[ma_tmp]  << ", st1_Tw[9] = " << st1_Tw[9] << endl;     
					        DIF_DATARECORD <<" A_B1R10["<<ma_tmp<<"]: "<<A_B1R10[ma_tmp] << ", st1_Tw[10] = " << st1_Tw[10] << endl;    
					        DIF_DATARECORD <<" A_B1R11["<<ma_tmp<<"]: "<<A_B1R11[ma_tmp] << ", st1_Tw[11] = " << st1_Tw[11] << endl;    
					        DIF_DATARECORD <<" A_B1R12["<<ma_tmp<<"]: "<<A_B1R12[ma_tmp] << ", st1_Tw[12] = " << st1_Tw[12] << endl;    
					        DIF_DATARECORD <<" A_B1R13["<<ma_tmp<<"]: "<<A_B1R13[ma_tmp] << ", st1_Tw[13] = " << st1_Tw[13] << endl;    
					        DIF_DATARECORD <<" A_B1R14["<<ma_tmp<<"]: "<<A_B1R14[ma_tmp] << ", st1_Tw[14] = " << st1_Tw[14] << endl;    
					        DIF_DATARECORD <<" A_B1R15["<<ma_tmp<<"]: "<<A_B1R15[ma_tmp] << ", st1_Tw[15] = " << st1_Tw[15] << endl;    
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
                            DIF_DATARECORD <<" A_B1R0["<<ma_tmp<<"]:  "<<A_B1R0[ma_tmp]  << ", st2_Tw[0] = " << st2_Tw[0] << endl;     
				            DIF_DATARECORD <<" A_B1R1["<<ma_tmp<<"]:  "<<A_B1R1[ma_tmp]  << ", st2_Tw[1] = " << st2_Tw[1] << endl;     
					        DIF_DATARECORD <<" A_B1R2["<<ma_tmp<<"]:  "<<A_B1R2[ma_tmp]  << ", st2_Tw[2] = " << st2_Tw[2] << endl;     
					        DIF_DATARECORD <<" A_B1R3["<<ma_tmp<<"]:  "<<A_B1R3[ma_tmp]  << ", st2_Tw[3] = " << st2_Tw[3] << endl;     
					        DIF_DATARECORD <<" A_B1R4["<<ma_tmp<<"]:  "<<A_B1R4[ma_tmp]  << ", st2_Tw[4] = " << st2_Tw[4] << endl;     
					        DIF_DATARECORD <<" A_B1R5["<<ma_tmp<<"]:  "<<A_B1R5[ma_tmp]  << ", st2_Tw[5] = " << st2_Tw[5] << endl;     
					        DIF_DATARECORD <<" A_B1R6["<<ma_tmp<<"]:  "<<A_B1R6[ma_tmp]  << ", st2_Tw[6] = " << st2_Tw[6] << endl;     
					        DIF_DATARECORD <<" A_B1R7["<<ma_tmp<<"]:  "<<A_B1R7[ma_tmp]  << ", st2_Tw[7] = " << st2_Tw[7] << endl;     
					        DIF_DATARECORD <<" A_B1R8["<<ma_tmp<<"]:  "<<A_B1R8[ma_tmp]  << ", st2_Tw[8] = " << st2_Tw[8] << endl;     
					        DIF_DATARECORD <<" A_B1R9["<<ma_tmp<<"]:  "<<A_B1R9[ma_tmp]  << ", st2_Tw[9] = " << st2_Tw[9] << endl;     
					        DIF_DATARECORD <<" A_B1R10["<<ma_tmp<<"]: "<<A_B1R10[ma_tmp] << ", st2_Tw[10] = " << st2_Tw[10] << endl;    
					        DIF_DATARECORD <<" A_B1R11["<<ma_tmp<<"]: "<<A_B1R11[ma_tmp] << ", st2_Tw[11] = " << st2_Tw[11] << endl;    
					        DIF_DATARECORD <<" A_B1R12["<<ma_tmp<<"]: "<<A_B1R12[ma_tmp] << ", st2_Tw[12] = " << st2_Tw[12] << endl;    
					        DIF_DATARECORD <<" A_B1R13["<<ma_tmp<<"]: "<<A_B1R13[ma_tmp] << ", st2_Tw[13] = " << st2_Tw[13] << endl;    
					        DIF_DATARECORD <<" A_B1R14["<<ma_tmp<<"]: "<<A_B1R14[ma_tmp] << ", st2_Tw[14] = " << st2_Tw[14] << endl;    
					        DIF_DATARECORD <<" A_B1R15["<<ma_tmp<<"]: "<<A_B1R15[ma_tmp] << ", st2_Tw[15] = " << st2_Tw[15] << endl;    
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
                            DIF_DATARECORD <<" A_B1R0["<<ma_tmp<<"]:  "<<A_B1R0[ma_tmp]  << ", st3_Tw[0] = " << 1 << endl;     
				            DIF_DATARECORD <<" A_B1R1["<<ma_tmp<<"]:  "<<A_B1R1[ma_tmp]  << ", st3_Tw[1] = " << 1 << endl;     
					        DIF_DATARECORD <<" A_B1R2["<<ma_tmp<<"]:  "<<A_B1R2[ma_tmp]  << ", st3_Tw[2] = " << 1 << endl;     
					        DIF_DATARECORD <<" A_B1R3["<<ma_tmp<<"]:  "<<A_B1R3[ma_tmp]  << ", st3_Tw[3] = " << 1 << endl;     
					        DIF_DATARECORD <<" A_B1R4["<<ma_tmp<<"]:  "<<A_B1R4[ma_tmp]  << ", st3_Tw[4] = " << 1 << endl;     
					        DIF_DATARECORD <<" A_B1R5["<<ma_tmp<<"]:  "<<A_B1R5[ma_tmp]  << ", st3_Tw[5] = " << 1 << endl;     
					        DIF_DATARECORD <<" A_B1R6["<<ma_tmp<<"]:  "<<A_B1R6[ma_tmp]  << ", st3_Tw[6] = " << 1 << endl;     
					        DIF_DATARECORD <<" A_B1R7["<<ma_tmp<<"]:  "<<A_B1R7[ma_tmp]  << ", st3_Tw[7] = " << 1 << endl;     
					        DIF_DATARECORD <<" A_B1R8["<<ma_tmp<<"]:  "<<A_B1R8[ma_tmp]  << ", st3_Tw[8] = " << 1 << endl;     
					        DIF_DATARECORD <<" A_B1R9["<<ma_tmp<<"]:  "<<A_B1R9[ma_tmp]  << ", st3_Tw[9] = " << 1 << endl;     
					        DIF_DATARECORD <<" A_B1R10["<<ma_tmp<<"]: "<<A_B1R10[ma_tmp] << ", st3_Tw[10] = " << 1 << endl;    
					        DIF_DATARECORD <<" A_B1R11["<<ma_tmp<<"]: "<<A_B1R11[ma_tmp] << ", st3_Tw[11] = " << 1 << endl;    
					        DIF_DATARECORD <<" A_B1R12["<<ma_tmp<<"]: "<<A_B1R12[ma_tmp] << ", st3_Tw[12] = " << 1 << endl;    
					        DIF_DATARECORD <<" A_B1R13["<<ma_tmp<<"]: "<<A_B1R13[ma_tmp] << ", st3_Tw[13] = " << 1 << endl;    
					        DIF_DATARECORD <<" A_B1R14["<<ma_tmp<<"]: "<<A_B1R14[ma_tmp] << ", st3_Tw[14] = " << 1 << endl;    
					        DIF_DATARECORD <<" A_B1R15["<<ma_tmp<<"]: "<<A_B1R15[ma_tmp] << ", st3_Tw[15] = " << 1 << endl;    
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
					
                    DIF_DATARECORD <<" ------after BU compute and Mul-------" << std::endl;
				    DIF_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<<"\n";
				    DIF_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<<"\n";
					DIF_DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<<A_B1R2[ma_tmp]<<"\n";
					DIF_DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<<A_B1R3[ma_tmp]<<"\n";
					DIF_DATARECORD <<"A_B1R4["<<ma_tmp<<"]: "<<A_B1R4[ma_tmp]<<"\n";
					DIF_DATARECORD <<"A_B1R5["<<ma_tmp<<"]: "<<A_B1R5[ma_tmp]<<"\n";
					DIF_DATARECORD <<"A_B1R6["<<ma_tmp<<"]: "<<A_B1R6[ma_tmp]<<"\n";
					DIF_DATARECORD <<"A_B1R7["<<ma_tmp<<"]: "<<A_B1R7[ma_tmp]<<"\n";
					DIF_DATARECORD <<"A_B1R8["<<ma_tmp<<"]: "<<A_B1R8[ma_tmp]<<"\n";
					DIF_DATARECORD <<"A_B1R9["<<ma_tmp<<"]: "<<A_B1R9[ma_tmp]<<"\n";
					DIF_DATARECORD <<"A_B1R10["<<ma_tmp<<"]: "<<A_B1R10[ma_tmp]<<"\n";
					DIF_DATARECORD <<"A_B1R11["<<ma_tmp<<"]: "<<A_B1R11[ma_tmp]<<"\n";
					DIF_DATARECORD <<"A_B1R12["<<ma_tmp<<"]: "<<A_B1R12[ma_tmp]<<"\n";
					DIF_DATARECORD <<"A_B1R13["<<ma_tmp<<"]: "<<A_B1R13[ma_tmp]<<"\n";
					DIF_DATARECORD <<"A_B1R14["<<ma_tmp<<"]: "<<A_B1R14[ma_tmp]<<"\n";
					DIF_DATARECORD <<"A_B1R15["<<ma_tmp<<"]: "<<A_B1R15[ma_tmp]<<"\n";
					DIF_DATARECORD <<"--------------------------------------------------------------------\n";
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
			   //std::cout << "BC_tmp = " << BC_tmp << ", idx0 = " << index0 << ", MA = " << ma_tmp << endl;
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
			   //std::cout << "BC_tmp = " << BC_tmp << ", idx0 = " << index0 << ", MA = " << ma_tmp << endl;
			}			
		}		
	}
	
}

void DIF_NTTSPMB::DIF_NTT_r4_r2(std::vector<ZZ> &A,
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
	
	std::ofstream DIF_DATARECORD("./DIF_NTT_R4_R2_SPMB.txt");
	//-------------Bit Reverse--------------
	BitOperate RR, IntToVec, VecToInt;
	long long bit_width = log2(N);
	long long bit_width_s = log(radix)/log(2);
	vector<long long> bit_array;
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
	DIF_DATARECORD <<"Stage: "<< Stage<<"\n";
	
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
	int debug = 0;
	vector<vector<ZZ > > ROM0;
    vector<vector<ZZ > > ROM1, ROM2;
    
    cout << "radix_r1 = " << radix_r1 << ", radix_r2 = " << radix_r2 << endl;
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
	DIF_DATARECORD << "init load over! \n";
	int tw_degree = 1; // siang
	//need modify to mult by twiddle factor
	for(int s = 0; s < Stage; s++){
		DIF_DATARECORD <<"----------------------------"<<"\n";
		DIF_DATARECORD <<"Stage: "<< s<<"\n";
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
				DIF_DATARECORD <<"----------------------------"<<"\n";
				DIF_DATARECORD <<"i: "<< i << ", j = " << j <<"\n";
				DIF_DATARECORD <<"twiddle factor =  " << factor << std::endl;
				gray_i  = Gray(i,group);
				BC_tmp  = j * group + gray_i;
				if(s == Stage-1)RR_R4_R2(BC_tmp,(2*s - 1),BC);
				else RR_R4(BC_tmp,s,BC);
				length = BC % tw_modulus_tmp;
				PowerMod(factor_t,factor,length,p);
				AGU_R4(BC,bn_tmp,ma_tmp);

				DIF_DATARECORD << "BC_tmp = " << BC_tmp << ", length: " << length <<
				", tw_dg: " <<  tw_degree * length << " = " << tw_degree << " * " << length <<"\n";	// TF_record
				DIF_DATARECORD <<"factor_t =  " << factor_t << ", p = " <<  p <<"\n";

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

				if(bn_tmp == 0){	
					if(j < 2)bn0_bc_tmp = BC_tmp;
					switch(s){
						case 0:
							DIF_DATARECORD <<"A_B0R0[" << ma_tmp << "] = " << A_B0R0[ma_tmp] << ", st0_Tw[0] = " << st0_Tw[0] << endl;
							DIF_DATARECORD <<"A_B0R1[" << ma_tmp << "] = " << A_B0R1[ma_tmp] << ", st0_Tw[1] = " << st0_Tw[1] << endl;
							DIF_DATARECORD <<"A_B0R2[" << ma_tmp << "] = " << A_B0R2[ma_tmp] << ", st0_Tw[2] = " << st0_Tw[2] << endl;
							DIF_DATARECORD <<"A_B0R3[" << ma_tmp << "] = " << A_B0R3[ma_tmp] << ", st0_Tw[3] = " << st0_Tw[3] << endl;
							if(!debug) Radix4_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp],A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
							if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],st0_Tw[0],p);
							if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],st0_Tw[1],p);
							if(!debug) MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],st0_Tw[2],p);
							if(!debug) MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],st0_Tw[3],p);
							break;
						case 1:
							DIF_DATARECORD <<"A_B0R0[" << ma_tmp << "] = " << A_B0R0[ma_tmp] << ", st1_Tw[0] = " << st1_Tw[0] << endl;
							DIF_DATARECORD <<"A_B0R1[" << ma_tmp << "] = " << A_B0R1[ma_tmp] << ", st1_Tw[1] = " << st1_Tw[1] << endl;
							DIF_DATARECORD <<"A_B0R2[" << ma_tmp << "] = " << A_B0R2[ma_tmp] << ", st1_Tw[2] = " << st1_Tw[2] << endl;
							DIF_DATARECORD <<"A_B0R3[" << ma_tmp << "] = " << A_B0R3[ma_tmp] << ", st1_Tw[3] = " << st1_Tw[3] << endl;
							if(!debug) Radix4_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp],A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
							if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],st1_Tw[0],p);
							if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],st1_Tw[1],p);
							if(!debug) MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],st1_Tw[2],p);
							if(!debug) MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],st1_Tw[3],p);
							break;
						case 2:
							DIF_DATARECORD <<"A_B0R0[" << ma_tmp << "] = " << A_B0R0[ma_tmp] << ", st2_Tw[0] = " << st2_Tw[0] << endl;
							DIF_DATARECORD <<"A_B0R1[" << ma_tmp << "] = " << A_B0R1[ma_tmp] << ", st2_Tw[1] = " << st2_Tw[1] << endl;
							DIF_DATARECORD <<"A_B0R2[" << ma_tmp << "] = " << A_B0R2[ma_tmp] << ", st2_Tw[2] = " << st2_Tw[2] << endl;
							DIF_DATARECORD <<"A_B0R3[" << ma_tmp << "] = " << A_B0R3[ma_tmp] << ", st2_Tw[3] = " << st2_Tw[3] << endl;
							if(!debug) Radix4_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp],A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
							if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],st2_Tw[0],p);
							if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],st2_Tw[1],p);
							if(!debug) MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],st2_Tw[2],p);
							if(!debug) MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],st2_Tw[3],p);
							break;
						case 3:
							DIF_DATARECORD <<"A_B0R0[" << ma_tmp << "] = " << A_B0R0[ma_tmp] << ", st3_Tw[0] = " << 1 << endl;
							DIF_DATARECORD <<"A_B0R1[" << ma_tmp << "] = " << A_B0R1[ma_tmp] << ", st3_Tw[1] = " << 1 << endl;
							DIF_DATARECORD <<"A_B0R2[" << ma_tmp << "] = " << A_B0R2[ma_tmp] << ", st3_Tw[2] = " << 1 << endl;
							DIF_DATARECORD <<"A_B0R3[" << ma_tmp << "] = " << A_B0R3[ma_tmp] << ", st3_Tw[3] = " << 1 << endl;
							if(!debug) Radix4_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp],A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
							if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],1,p);
							if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],1,p);
							if(!debug) MulMod(A_B0R2[ma_tmp],A_B0R2[ma_tmp],1,p);
							if(!debug) MulMod(A_B0R3[ma_tmp],A_B0R3[ma_tmp],1,p);
							break;
						default:
							break;
					}			
					DIF_DATARECORD << "*****after Mul*****" << std::endl;
					DIF_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp] << endl;
					DIF_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp] << endl;
					DIF_DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp] << endl;
					DIF_DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp] << endl;
					if(j <  2)bn0_ma_reg1 = ma_tmp;
					if(j >= 2)bn0_ma_reg2 = ma_tmp;
					
				}
			    else {
					if(j < 2)bn1_bc_tmp = BC_tmp;
					switch(s){
						case 0:
							DIF_DATARECORD <<"A_B1R0[" << ma_tmp << "] = " << A_B1R0[ma_tmp] << ", st0_Tw[0] = " << st0_Tw[0] << endl;
							DIF_DATARECORD <<"A_B1R1[" << ma_tmp << "] = " << A_B1R1[ma_tmp] << ", st0_Tw[1] = " << st0_Tw[1] << endl;
							DIF_DATARECORD <<"A_B1R2[" << ma_tmp << "] = " << A_B1R2[ma_tmp] << ", st0_Tw[2] = " << st0_Tw[2] << endl;
							DIF_DATARECORD <<"A_B1R3[" << ma_tmp << "] = " << A_B1R3[ma_tmp] << ", st0_Tw[3] = " << st0_Tw[3] << endl;
							if(!debug) Radix4_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp],A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
							if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],st0_Tw[0],p);
							if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],st0_Tw[1],p);
							if(!debug) MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],st0_Tw[2],p);
							if(!debug) MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],st0_Tw[3],p);
							break;
						case 1:
							DIF_DATARECORD <<"A_B1R0[" << ma_tmp << "] = " << A_B1R0[ma_tmp] << ", st1_Tw[0] = " << st1_Tw[0] << endl;
							DIF_DATARECORD <<"A_B1R1[" << ma_tmp << "] = " << A_B1R1[ma_tmp] << ", st1_Tw[1] = " << st1_Tw[1] << endl;
							DIF_DATARECORD <<"A_B1R2[" << ma_tmp << "] = " << A_B1R2[ma_tmp] << ", st1_Tw[2] = " << st1_Tw[2] << endl;
							DIF_DATARECORD <<"A_B1R3[" << ma_tmp << "] = " << A_B1R3[ma_tmp] << ", st1_Tw[3] = " << st1_Tw[3] << endl;
							if(!debug) Radix4_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp],A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
							if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],st1_Tw[0],p);
							if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],st1_Tw[1],p);
							if(!debug) MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],st1_Tw[2],p);
							if(!debug) MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],st1_Tw[3],p);
							break;
						case 2:
							DIF_DATARECORD <<"A_B1R0[" << ma_tmp << "] = " << A_B1R0[ma_tmp] << ", st2_Tw[0] = " << st2_Tw[0] << endl;
							DIF_DATARECORD <<"A_B1R1[" << ma_tmp << "] = " << A_B1R1[ma_tmp] << ", st2_Tw[1] = " << st2_Tw[1] << endl;
							DIF_DATARECORD <<"A_B1R2[" << ma_tmp << "] = " << A_B1R2[ma_tmp] << ", st2_Tw[2] = " << st2_Tw[2] << endl;
							DIF_DATARECORD <<"A_B1R3[" << ma_tmp << "] = " << A_B1R3[ma_tmp] << ", st2_Tw[3] = " << st2_Tw[3] << endl;
							if(!debug) Radix4_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp],A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
							if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],st2_Tw[0],p);
							if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],st2_Tw[1],p);
							if(!debug) MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],st2_Tw[2],p);
							if(!debug) MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],st2_Tw[3],p);
							break;
						case 3:
							DIF_DATARECORD <<"A_B1R0[" << ma_tmp << "] = " << A_B1R0[ma_tmp] << ", st3_Tw[0] = " << 1 << endl;
							DIF_DATARECORD <<"A_B1R1[" << ma_tmp << "] = " << A_B1R1[ma_tmp] << ", st3_Tw[1] = " << 1 << endl;
							DIF_DATARECORD <<"A_B1R2[" << ma_tmp << "] = " << A_B1R2[ma_tmp] << ", st3_Tw[2] = " << 1 << endl;
							DIF_DATARECORD <<"A_B1R3[" << ma_tmp << "] = " << A_B1R3[ma_tmp] << ", st3_Tw[3] = " << 1 << endl;
							if(!debug) Radix4_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp],A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
							if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],1,p);
							if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],1,p);
							if(!debug) MulMod(A_B1R2[ma_tmp],A_B1R2[ma_tmp],1,p);
							if(!debug) MulMod(A_B1R3[ma_tmp],A_B1R3[ma_tmp],1,p);
							break;
						default:
							break;
					}
					
					DIF_DATARECORD << "*****after mul*****" << std::endl;
					DIF_DATARECORD << "A_B1R0[" << ma_tmp << "] = " << A_B1R0[ma_tmp] << endl;
					DIF_DATARECORD << "A_B1R1[" << ma_tmp << "] = " << A_B1R1[ma_tmp] << endl;
					DIF_DATARECORD << "A_B1R2[" << ma_tmp << "] = " << A_B1R2[ma_tmp] << endl;
					DIF_DATARECORD << "A_B1R3[" << ma_tmp << "] = " << A_B1R3[ma_tmp] << endl;
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
	DIF_DATARECORD <<"------------radix_2 part----------------"<< std::endl;	
	// radix-2 FFT compute
	for(int i = 0; i < group ; i++){
		for(int j=0; j < radix; j++){
			DIF_DATARECORD <<"------------------"<< std::endl;	
			gray_i  = Gray(i,group);
			BC_tmp  = j * group + gray_i;
			BC = BC_tmp;
			AGU_R4(BC,bn_tmp,ma_tmp);
			if(bn_tmp == 0){
				if(j < 2)bn0_bc_tmp = BC_tmp;
				DIF_DATARECORD <<"A_B0R0[" << ma_tmp << "] = " << A_B0R0[ma_tmp] <<", A_B0R1[" << ma_tmp << "] = " << A_B0R1[ma_tmp] << std::endl;
				DIF_DATARECORD <<"A_B0R2[" << ma_tmp << "] = " << A_B0R2[ma_tmp] <<", A_B0R3[" << ma_tmp << "] = " << A_B0R3[ma_tmp] << std::endl;
				Radix2_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp]);
				Radix2_BU(A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
				if(j <  2)bn0_ma_reg1 = ma_tmp;
				if(j >= 2)bn0_ma_reg2 = ma_tmp;
			}else {
			    if(j < 2)bn1_bc_tmp = BC_tmp;
				DIF_DATARECORD <<"A_B1R0[" << ma_tmp << "] = " << A_B1R0[ma_tmp] <<", A_B1R1[" << ma_tmp << "] = " << A_B1R1[ma_tmp] << std::endl;
				DIF_DATARECORD <<"A_B1R2[" << ma_tmp << "] = " << A_B1R2[ma_tmp] <<", A_B1R3[" << ma_tmp << "] = " << A_B1R3[ma_tmp] << std::endl;
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

    DIF_DATARECORD << "---------------------------------------------------------\n";
	DIF_DATARECORD << "*********************************************************\n";
	DIF_DATARECORD << "                  FINAL!!!!                              \n";
	DIF_DATARECORD << "*********************************************************\n";
	DIF_DATARECORD << "---------------------------------------------------------\n";
	
    int BC_bit_size;
	BC_bit_size = (int)ceil(log2(N/4));
	DIF_DATARECORD << "BC_bit_size: "<< BC_bit_size <<"\n";
	int index_out;
	for(int i = 0; i < group; i++){
    	for(int j = 0;j < radix;j++){
			DIF_DATARECORD << "-------------------------" << std::endl;
			gray_i  = Gray(i,group);
    		BC_tmp  = j * group + gray_i;
			REORDERBC_R4_R2_OUT(BC_tmp,BC);
			BC = BC_tmp;
			DIF_DATARECORD << "BC_tmp: " << BC_tmp <<"\n";
			DIF_DATARECORD << "BC: " << BC <<"\n";
    		AGU_R4(BC,bn_tmp,ma_tmp);
			DIF_DATARECORD << "bn_tmp: " << bn_tmp <<"\n";
			DIF_DATARECORD << "ma_tmp: " << ma_tmp <<"\n";
			index_out = (int)i * ( N / group) + 4 * j ;
			DIF_DATARECORD << "index_out: " << index_out << std::endl;
    		if(bn_tmp == 0){
     		   A[index_out+0] = A_B0R0[ma_tmp];
    		   A[index_out+1] = A_B0R1[ma_tmp];
    		   A[index_out+2] = A_B0R2[ma_tmp];
    		   A[index_out+3] = A_B0R3[ma_tmp];
			   DIF_DATARECORD << "BN: 0 \n";
			   DIF_DATARECORD << "A["<< index_out      << "]: " << A[index_out+0] << " \n";
			   DIF_DATARECORD << "A["<< index_out + 1  << "]: " << A[index_out+1] << " \n";
			   DIF_DATARECORD << "A["<< index_out + 2  << "]: " << A[index_out+2] << " \n";
			   DIF_DATARECORD << "A["<< index_out + 3  << "]: " << A[index_out+3] << " \n";
    		}
    	    else {
    		   A[index_out+0] = A_B1R0[ma_tmp];
    		   A[index_out+1] = A_B1R1[ma_tmp];
    		   A[index_out+2] = A_B1R2[ma_tmp];
    		   A[index_out+3] = A_B1R3[ma_tmp];
			   DIF_DATARECORD << "BN: 1 \n";
			   DIF_DATARECORD << "A["<< index_out      << "]: " << A[index_out+0] << " \n";
			   DIF_DATARECORD << "A["<< index_out + 1  << "]: " << A[index_out+1] << " \n";
			   DIF_DATARECORD << "A["<< index_out + 2  << "]: " << A[index_out+2] << " \n";
			   DIF_DATARECORD << "A["<< index_out + 3  << "]: " << A[index_out+3] << " \n";
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
	DIF_DATARECORD << "---------------------------------------------------------\n";
	DIF_DATARECORD << "----                  OVER!!!!!!!                   -----\n";
	DIF_DATARECORD << "---------------------------------------------------------\n";
	DIF_DATARECORD << "---------------------------------------------------------\n";
	DIF_DATARECORD.close();
}

void DIF_NTTSPMB::DIF_NTT_r16_r2(std::vector<ZZ> &A,std::vector<ZZ> &B0R0,std::vector<ZZ> &B0R1,std::vector<ZZ> &B0R2,std::vector<ZZ> &B0R3,
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
	
	std::ofstream DIF_DATARECORD("./DIF_NTT_R16_R2_SPMB.txt");
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
	int debug = 0;
	vector<vector<ZZ > > ROM0;
    vector<vector<ZZ > > ROM1, ROM2;
    
    cout << "radix_r1 = " << radix_r1 << ", radix_r2 = " << radix_r2 << endl;
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

	//radix-16 Stage
    Stage_double  = log2(N);
    Stage       =  (int)floor(Stage_double/4);
	BC_WIDTH    =  (int)ceil(log2(N/16));
	offset      =  (int)N /  16;
	word_size   =  (int)N / (2 * 16);
	group       =  (int)N / (256);
	tw_modulus  =  (int)N /  16;
	bit_array_tmp.resize(BC_WIDTH);

    DIF_DATARECORD << "group: "    << group << "\n";
    DIF_DATARECORD << "BC_WIDTH: " << BC_WIDTH << "\n";
	
	
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
	DIF_DATARECORD <<"radix-16 computing stage:  "<< Stage <<"\n";
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
		DIF_DATARECORD <<"---------------------------------\n";
		DIF_DATARECORD <<"Now Stage: "<< s <<"\n";
		tw_modulus_tmp  = tw_modulus >> ( 4 * s);
		for(int i = 0 ;i < group;i++){
			bn0_bc_tmp  = 0;
			bn1_bc_tmp  = 0;
			for(int j = 0;j < radix;j++){
				DIF_DATARECORD << "-------------------------------------"<< std::endl;
				DIF_DATARECORD << "i = " << i << ", j = " << j << std::endl;
				gray_i  = Gray(i,group);
				BC_tmp  = j * group + gray_i;
				DIF_DATARECORD << "BC_tmp: " << BC_tmp << "\n";
				if(s == Stage - 1) RR_R16_R2(BC_tmp,(4 * s - 3),BC);
				else RR_R16(BC_tmp,s,BC);
				DIF_DATARECORD << "After RR_R16 , BC : " << BC << "\n";
				length = BC % tw_modulus_tmp;
				PowerMod(factor_t,factor,length,p);
				AGU_R16(BC,bn_tmp,ma_tmp);			

				DIF_DATARECORD << "BN : " << bn_tmp << "\n";
				DIF_DATARECORD << "MA : " << ma_tmp << "\n";

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

				if(bn_tmp == 0){
					if(j < 2)bn0_bc_tmp = BC_tmp;
					switch(s){
						case 0:
							DIF_DATARECORD << "A_B0R0["<< ma_tmp <<"]: " << A_B0R0[ma_tmp]	<< ", st0_Tw[0] = " << st0_Tw[0] << endl;
				    		DIF_DATARECORD << "A_B0R1["<< ma_tmp <<"]: " << A_B0R1[ma_tmp]	<< ", st0_Tw[1] = " << st0_Tw[1] << endl;
				    		DIF_DATARECORD << "A_B0R2["<< ma_tmp <<"]: " << A_B0R2[ma_tmp]	<< ", st0_Tw[2] = " << st0_Tw[2] << endl;
				    		DIF_DATARECORD << "A_B0R3["<< ma_tmp <<"]: " << A_B0R3[ma_tmp]	<< ", st0_Tw[3] = " << st0_Tw[3] << endl;
				    		DIF_DATARECORD << "A_B0R4["<< ma_tmp <<"]: " << A_B0R4[ma_tmp]	<< ", st0_Tw[4] = " << st0_Tw[4] << endl;
				    		DIF_DATARECORD << "A_B0R5["<< ma_tmp <<"]: " << A_B0R5[ma_tmp]	<< ", st0_Tw[5] = " << st0_Tw[5] << endl;
				    		DIF_DATARECORD << "A_B0R6["<< ma_tmp <<"]: " << A_B0R6[ma_tmp]	<< ", st0_Tw[6] = " << st0_Tw[6] << endl;
				    		DIF_DATARECORD << "A_B0R7["<< ma_tmp <<"]: " << A_B0R7[ma_tmp]	<< ", st0_Tw[7] = " << st0_Tw[7] << endl;
				    		DIF_DATARECORD << "A_B0R8["<< ma_tmp <<"]: " << A_B0R8[ma_tmp]	<< ", st0_Tw[8] = " << st0_Tw[8] << endl;
				    		DIF_DATARECORD << "A_B0R9["<< ma_tmp <<"]: " << A_B0R9[ma_tmp]	<< ", st0_Tw[9] = " << st0_Tw[9] << endl;
				    		DIF_DATARECORD << "A_B0R10["<< ma_tmp <<"]: "<< A_B0R10[ma_tmp]	<< ", st0_Tw[10] = " << st0_Tw[10] << endl;
				    		DIF_DATARECORD << "A_B0R11["<< ma_tmp <<"]: "<< A_B0R11[ma_tmp]	<< ", st0_Tw[11] = " << st0_Tw[11] << endl;
				    		DIF_DATARECORD << "A_B0R12["<< ma_tmp <<"]: "<< A_B0R12[ma_tmp]	<< ", st0_Tw[12] = " << st0_Tw[12] << endl;
				    		DIF_DATARECORD << "A_B0R13["<< ma_tmp <<"]: "<< A_B0R13[ma_tmp]	<< ", st0_Tw[13] = " << st0_Tw[13] << endl;
				    		DIF_DATARECORD << "A_B0R14["<< ma_tmp <<"]: "<< A_B0R14[ma_tmp]	<< ", st0_Tw[14] = " << st0_Tw[14] << endl;
				    		DIF_DATARECORD << "A_B0R15["<< ma_tmp <<"]: "<< A_B0R15[ma_tmp]	<< ", st0_Tw[15] = " << st0_Tw[15] << endl;
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
							DIF_DATARECORD << "A_B0R0["<< ma_tmp <<"]: " << A_B0R0[ma_tmp]	<< ", st1_Tw[0] = " << st1_Tw[0] << endl;
				    		DIF_DATARECORD << "A_B0R1["<< ma_tmp <<"]: " << A_B0R1[ma_tmp]	<< ", st1_Tw[1] = " << st1_Tw[1] << endl;
				    		DIF_DATARECORD << "A_B0R2["<< ma_tmp <<"]: " << A_B0R2[ma_tmp]	<< ", st1_Tw[2] = " << st1_Tw[2] << endl;
				    		DIF_DATARECORD << "A_B0R3["<< ma_tmp <<"]: " << A_B0R3[ma_tmp]	<< ", st1_Tw[3] = " << st1_Tw[3] << endl;
				    		DIF_DATARECORD << "A_B0R4["<< ma_tmp <<"]: " << A_B0R4[ma_tmp]	<< ", st1_Tw[4] = " << st1_Tw[4] << endl;
				    		DIF_DATARECORD << "A_B0R5["<< ma_tmp <<"]: " << A_B0R5[ma_tmp]	<< ", st1_Tw[5] = " << st1_Tw[5] << endl;
				    		DIF_DATARECORD << "A_B0R6["<< ma_tmp <<"]: " << A_B0R6[ma_tmp]	<< ", st1_Tw[6] = " << st1_Tw[6] << endl;
				    		DIF_DATARECORD << "A_B0R7["<< ma_tmp <<"]: " << A_B0R7[ma_tmp]	<< ", st1_Tw[7] = " << st1_Tw[7] << endl;
				    		DIF_DATARECORD << "A_B0R8["<< ma_tmp <<"]: " << A_B0R8[ma_tmp]	<< ", st1_Tw[8] = " << st1_Tw[8] << endl;
				    		DIF_DATARECORD << "A_B0R9["<< ma_tmp <<"]: " << A_B0R9[ma_tmp]	<< ", st1_Tw[9] = " << st1_Tw[9] << endl;
				    		DIF_DATARECORD << "A_B0R10["<< ma_tmp <<"]: "<< A_B0R10[ma_tmp]	<< ", st1_Tw[10] = " << st1_Tw[10] << endl;
				    		DIF_DATARECORD << "A_B0R11["<< ma_tmp <<"]: "<< A_B0R11[ma_tmp]	<< ", st1_Tw[11] = " << st1_Tw[11] << endl;
				    		DIF_DATARECORD << "A_B0R12["<< ma_tmp <<"]: "<< A_B0R12[ma_tmp]	<< ", st1_Tw[12] = " << st1_Tw[12] << endl;
				    		DIF_DATARECORD << "A_B0R13["<< ma_tmp <<"]: "<< A_B0R13[ma_tmp]	<< ", st1_Tw[13] = " << st1_Tw[13] << endl;
				    		DIF_DATARECORD << "A_B0R14["<< ma_tmp <<"]: "<< A_B0R14[ma_tmp]	<< ", st1_Tw[14] = " << st1_Tw[14] << endl;
				    		DIF_DATARECORD << "A_B0R15["<< ma_tmp <<"]: "<< A_B0R15[ma_tmp]	<< ", st1_Tw[15] = " << st1_Tw[15] << endl;
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
							DIF_DATARECORD << "A_B0R0["<< ma_tmp <<"]: " << A_B0R0[ma_tmp]	<< ", st2_Tw[0] = " << st2_Tw[0] << endl;
				    		DIF_DATARECORD << "A_B0R1["<< ma_tmp <<"]: " << A_B0R1[ma_tmp]	<< ", st2_Tw[1] = " << st2_Tw[1] << endl;
				    		DIF_DATARECORD << "A_B0R2["<< ma_tmp <<"]: " << A_B0R2[ma_tmp]	<< ", st2_Tw[2] = " << st2_Tw[2] << endl;
				    		DIF_DATARECORD << "A_B0R3["<< ma_tmp <<"]: " << A_B0R3[ma_tmp]	<< ", st2_Tw[3] = " << st2_Tw[3] << endl;
				    		DIF_DATARECORD << "A_B0R4["<< ma_tmp <<"]: " << A_B0R4[ma_tmp]	<< ", st2_Tw[4] = " << st2_Tw[4] << endl;
				    		DIF_DATARECORD << "A_B0R5["<< ma_tmp <<"]: " << A_B0R5[ma_tmp]	<< ", st2_Tw[5] = " << st2_Tw[5] << endl;
				    		DIF_DATARECORD << "A_B0R6["<< ma_tmp <<"]: " << A_B0R6[ma_tmp]	<< ", st2_Tw[6] = " << st2_Tw[6] << endl;
				    		DIF_DATARECORD << "A_B0R7["<< ma_tmp <<"]: " << A_B0R7[ma_tmp]	<< ", st2_Tw[7] = " << st2_Tw[7] << endl;
				    		DIF_DATARECORD << "A_B0R8["<< ma_tmp <<"]: " << A_B0R8[ma_tmp]	<< ", st2_Tw[8] = " << st2_Tw[8] << endl;
				    		DIF_DATARECORD << "A_B0R9["<< ma_tmp <<"]: " << A_B0R9[ma_tmp]	<< ", st2_Tw[9] = " << st2_Tw[9] << endl;
				    		DIF_DATARECORD << "A_B0R10["<< ma_tmp <<"]: "<< A_B0R10[ma_tmp]	<< ", st2_Tw[10] = " << st2_Tw[10] << endl;
				    		DIF_DATARECORD << "A_B0R11["<< ma_tmp <<"]: "<< A_B0R11[ma_tmp]	<< ", st2_Tw[11] = " << st2_Tw[11] << endl;
				    		DIF_DATARECORD << "A_B0R12["<< ma_tmp <<"]: "<< A_B0R12[ma_tmp]	<< ", st2_Tw[12] = " << st2_Tw[12] << endl;
				    		DIF_DATARECORD << "A_B0R13["<< ma_tmp <<"]: "<< A_B0R13[ma_tmp]	<< ", st2_Tw[13] = " << st2_Tw[13] << endl;
				    		DIF_DATARECORD << "A_B0R14["<< ma_tmp <<"]: "<< A_B0R14[ma_tmp]	<< ", st2_Tw[14] = " << st2_Tw[14] << endl;
				    		DIF_DATARECORD << "A_B0R15["<< ma_tmp <<"]: "<< A_B0R15[ma_tmp]	<< ", st2_Tw[15] = " << st2_Tw[15] << endl;
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
							DIF_DATARECORD << "A_B0R0["<< ma_tmp <<"]: " << A_B0R0[ma_tmp]	<< ", st3_Tw[0] = " << 1 << endl;
				    		DIF_DATARECORD << "A_B0R1["<< ma_tmp <<"]: " << A_B0R1[ma_tmp]	<< ", st3_Tw[1] = " << 1 << endl;
				    		DIF_DATARECORD << "A_B0R2["<< ma_tmp <<"]: " << A_B0R2[ma_tmp]	<< ", st3_Tw[2] = " << 1 << endl;
				    		DIF_DATARECORD << "A_B0R3["<< ma_tmp <<"]: " << A_B0R3[ma_tmp]	<< ", st3_Tw[3] = " << 1 << endl;
				    		DIF_DATARECORD << "A_B0R4["<< ma_tmp <<"]: " << A_B0R4[ma_tmp]	<< ", st3_Tw[4] = " << 1 << endl;
				    		DIF_DATARECORD << "A_B0R5["<< ma_tmp <<"]: " << A_B0R5[ma_tmp]	<< ", st3_Tw[5] = " << 1 << endl;
				    		DIF_DATARECORD << "A_B0R6["<< ma_tmp <<"]: " << A_B0R6[ma_tmp]	<< ", st3_Tw[6] = " << 1 << endl;
				    		DIF_DATARECORD << "A_B0R7["<< ma_tmp <<"]: " << A_B0R7[ma_tmp]	<< ", st3_Tw[7] = " << 1 << endl;
				    		DIF_DATARECORD << "A_B0R8["<< ma_tmp <<"]: " << A_B0R8[ma_tmp]	<< ", st3_Tw[8] = " << 1 << endl;
				    		DIF_DATARECORD << "A_B0R9["<< ma_tmp <<"]: " << A_B0R9[ma_tmp]	<< ", st3_Tw[9] = " << 1 << endl;
				    		DIF_DATARECORD << "A_B0R10["<< ma_tmp <<"]: "<< A_B0R10[ma_tmp]	<< ", st3_Tw[10] = " << 1 << endl;
				    		DIF_DATARECORD << "A_B0R11["<< ma_tmp <<"]: "<< A_B0R11[ma_tmp]	<< ", st3_Tw[11] = " << 1 << endl;
				    		DIF_DATARECORD << "A_B0R12["<< ma_tmp <<"]: "<< A_B0R12[ma_tmp]	<< ", st3_Tw[12] = " << 1 << endl;
				    		DIF_DATARECORD << "A_B0R13["<< ma_tmp <<"]: "<< A_B0R13[ma_tmp]	<< ", st3_Tw[13] = " << 1 << endl;
				    		DIF_DATARECORD << "A_B0R14["<< ma_tmp <<"]: "<< A_B0R14[ma_tmp]	<< ", st3_Tw[14] = " << 1 << endl;
				    		DIF_DATARECORD << "A_B0R15["<< ma_tmp <<"]: "<< A_B0R15[ma_tmp]	<< ", st3_Tw[15] = " << 1 << endl;
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
				    DIF_DATARECORD << "A_B1R0["<< ma_tmp <<"]: " << A_B1R0[ma_tmp]<<"\n";
				    DIF_DATARECORD << "A_B1R1["<< ma_tmp <<"]: " << A_B1R1[ma_tmp]<<"\n";
				    DIF_DATARECORD << "A_B1R2["<< ma_tmp <<"]: " << A_B1R2[ma_tmp]<<"\n";
				    DIF_DATARECORD << "A_B1R3["<< ma_tmp <<"]: " << A_B1R3[ma_tmp]<<"\n";
				    DIF_DATARECORD << "A_B1R4["<< ma_tmp <<"]: " << A_B1R4[ma_tmp]<<"\n";
				    DIF_DATARECORD << "A_B1R5["<< ma_tmp <<"]: " << A_B1R5[ma_tmp]<<"\n";
				    DIF_DATARECORD << "A_B1R6["<< ma_tmp <<"]: " << A_B1R6[ma_tmp]<<"\n";
				    DIF_DATARECORD << "A_B1R7["<< ma_tmp <<"]: " << A_B1R7[ma_tmp]<<"\n";
				    DIF_DATARECORD << "A_B1R8["<< ma_tmp <<"]: " << A_B1R8[ma_tmp]<<"\n";
				    DIF_DATARECORD << "A_B1R9["<< ma_tmp <<"]: " << A_B1R9[ma_tmp]<<"\n";
				    DIF_DATARECORD << "A_B1R10["<< ma_tmp <<"]: "<< A_B1R10[ma_tmp]<<"\n";
				    DIF_DATARECORD << "A_B1R11["<< ma_tmp <<"]: "<< A_B1R11[ma_tmp]<<"\n";
				    DIF_DATARECORD << "A_B1R12["<< ma_tmp <<"]: "<< A_B1R12[ma_tmp]<<"\n";
				    DIF_DATARECORD << "A_B1R13["<< ma_tmp <<"]: "<< A_B1R13[ma_tmp]<<"\n";
				    DIF_DATARECORD << "A_B1R14["<< ma_tmp <<"]: "<< A_B1R14[ma_tmp]<<"\n";
				    DIF_DATARECORD << "A_B1R15["<< ma_tmp <<"]: "<< A_B1R15[ma_tmp]<<"\n";
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
							DIF_DATARECORD << "A_B1R0["<< ma_tmp <<"]: " << A_B1R0[ma_tmp]	<< ", st0_Tw[0] = " << st0_Tw[0] << endl;
				    		DIF_DATARECORD << "A_B1R1["<< ma_tmp <<"]: " << A_B1R1[ma_tmp]	<< ", st0_Tw[1] = " << st0_Tw[1] << endl;
				    		DIF_DATARECORD << "A_B1R2["<< ma_tmp <<"]: " << A_B1R2[ma_tmp]	<< ", st0_Tw[2] = " << st0_Tw[2] << endl;
				    		DIF_DATARECORD << "A_B1R3["<< ma_tmp <<"]: " << A_B1R3[ma_tmp]	<< ", st0_Tw[3] = " << st0_Tw[3] << endl;
				    		DIF_DATARECORD << "A_B1R4["<< ma_tmp <<"]: " << A_B1R4[ma_tmp]	<< ", st0_Tw[4] = " << st0_Tw[4] << endl;
				    		DIF_DATARECORD << "A_B1R5["<< ma_tmp <<"]: " << A_B1R5[ma_tmp]	<< ", st0_Tw[5] = " << st0_Tw[5] << endl;
				    		DIF_DATARECORD << "A_B1R6["<< ma_tmp <<"]: " << A_B1R6[ma_tmp]	<< ", st0_Tw[6] = " << st0_Tw[6] << endl;
				    		DIF_DATARECORD << "A_B1R7["<< ma_tmp <<"]: " << A_B1R7[ma_tmp]	<< ", st0_Tw[7] = " << st0_Tw[7] << endl;
				    		DIF_DATARECORD << "A_B1R8["<< ma_tmp <<"]: " << A_B1R8[ma_tmp]	<< ", st0_Tw[8] = " << st0_Tw[8] << endl;
				    		DIF_DATARECORD << "A_B1R9["<< ma_tmp <<"]: " << A_B1R9[ma_tmp]	<< ", st0_Tw[9] = " << st0_Tw[9] << endl;
				    		DIF_DATARECORD << "A_B1R10["<< ma_tmp <<"]: "<< A_B1R10[ma_tmp]	<< ", st0_Tw[10] = " << st0_Tw[10] << endl;
				    		DIF_DATARECORD << "A_B1R11["<< ma_tmp <<"]: "<< A_B1R11[ma_tmp]	<< ", st0_Tw[11] = " << st0_Tw[11] << endl;
				    		DIF_DATARECORD << "A_B1R12["<< ma_tmp <<"]: "<< A_B1R12[ma_tmp]	<< ", st0_Tw[12] = " << st0_Tw[12] << endl;
				    		DIF_DATARECORD << "A_B1R13["<< ma_tmp <<"]: "<< A_B1R13[ma_tmp]	<< ", st0_Tw[13] = " << st0_Tw[13] << endl;
				    		DIF_DATARECORD << "A_B1R14["<< ma_tmp <<"]: "<< A_B1R14[ma_tmp]	<< ", st0_Tw[14] = " << st0_Tw[14] << endl;
				    		DIF_DATARECORD << "A_B1R15["<< ma_tmp <<"]: "<< A_B1R15[ma_tmp]	<< ", st0_Tw[15] = " << st0_Tw[15] << endl;					
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
							DIF_DATARECORD << "A_B1R0["<< ma_tmp <<"]: " << A_B1R0[ma_tmp]	<< ", st1_Tw[0] = " << st1_Tw[0] << endl;
				    		DIF_DATARECORD << "A_B1R1["<< ma_tmp <<"]: " << A_B1R1[ma_tmp]	<< ", st1_Tw[1] = " << st1_Tw[1] << endl;
				    		DIF_DATARECORD << "A_B1R2["<< ma_tmp <<"]: " << A_B1R2[ma_tmp]	<< ", st1_Tw[2] = " << st1_Tw[2] << endl;
				    		DIF_DATARECORD << "A_B1R3["<< ma_tmp <<"]: " << A_B1R3[ma_tmp]	<< ", st1_Tw[3] = " << st1_Tw[3] << endl;
				    		DIF_DATARECORD << "A_B1R4["<< ma_tmp <<"]: " << A_B1R4[ma_tmp]	<< ", st1_Tw[4] = " << st1_Tw[4] << endl;
				    		DIF_DATARECORD << "A_B1R5["<< ma_tmp <<"]: " << A_B1R5[ma_tmp]	<< ", st1_Tw[5] = " << st1_Tw[5] << endl;
				    		DIF_DATARECORD << "A_B1R6["<< ma_tmp <<"]: " << A_B1R6[ma_tmp]	<< ", st1_Tw[6] = " << st1_Tw[6] << endl;
				    		DIF_DATARECORD << "A_B1R7["<< ma_tmp <<"]: " << A_B1R7[ma_tmp]	<< ", st1_Tw[7] = " << st1_Tw[7] << endl;
				    		DIF_DATARECORD << "A_B1R8["<< ma_tmp <<"]: " << A_B1R8[ma_tmp]	<< ", st1_Tw[8] = " << st1_Tw[8] << endl;
				    		DIF_DATARECORD << "A_B1R9["<< ma_tmp <<"]: " << A_B1R9[ma_tmp]	<< ", st1_Tw[9] = " << st1_Tw[9] << endl;
				    		DIF_DATARECORD << "A_B1R10["<< ma_tmp <<"]: "<< A_B1R10[ma_tmp]	<< ", st1_Tw[10] = " << st1_Tw[10] << endl;
				    		DIF_DATARECORD << "A_B1R11["<< ma_tmp <<"]: "<< A_B1R11[ma_tmp]	<< ", st1_Tw[11] = " << st1_Tw[11] << endl;
				    		DIF_DATARECORD << "A_B1R12["<< ma_tmp <<"]: "<< A_B1R12[ma_tmp]	<< ", st1_Tw[12] = " << st1_Tw[12] << endl;
				    		DIF_DATARECORD << "A_B1R13["<< ma_tmp <<"]: "<< A_B1R13[ma_tmp]	<< ", st1_Tw[13] = " << st1_Tw[13] << endl;
				    		DIF_DATARECORD << "A_B1R14["<< ma_tmp <<"]: "<< A_B1R14[ma_tmp]	<< ", st1_Tw[14] = " << st1_Tw[14] << endl;
				    		DIF_DATARECORD << "A_B1R15["<< ma_tmp <<"]: "<< A_B1R15[ma_tmp]	<< ", st1_Tw[15] = " << st1_Tw[15] << endl;					
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
							DIF_DATARECORD << "A_B1R0["<< ma_tmp <<"]: " << A_B1R0[ma_tmp]	<< ", st2_Tw[0] = " << st2_Tw[0] << endl;
				    		DIF_DATARECORD << "A_B1R1["<< ma_tmp <<"]: " << A_B1R1[ma_tmp]	<< ", st2_Tw[1] = " << st2_Tw[1] << endl;
				    		DIF_DATARECORD << "A_B1R2["<< ma_tmp <<"]: " << A_B1R2[ma_tmp]	<< ", st2_Tw[2] = " << st2_Tw[2] << endl;
				    		DIF_DATARECORD << "A_B1R3["<< ma_tmp <<"]: " << A_B1R3[ma_tmp]	<< ", st2_Tw[3] = " << st2_Tw[3] << endl;
				    		DIF_DATARECORD << "A_B1R4["<< ma_tmp <<"]: " << A_B1R4[ma_tmp]	<< ", st2_Tw[4] = " << st2_Tw[4] << endl;
				    		DIF_DATARECORD << "A_B1R5["<< ma_tmp <<"]: " << A_B1R5[ma_tmp]	<< ", st2_Tw[5] = " << st2_Tw[5] << endl;
				    		DIF_DATARECORD << "A_B1R6["<< ma_tmp <<"]: " << A_B1R6[ma_tmp]	<< ", st2_Tw[6] = " << st2_Tw[6] << endl;
				    		DIF_DATARECORD << "A_B1R7["<< ma_tmp <<"]: " << A_B1R7[ma_tmp]	<< ", st2_Tw[7] = " << st2_Tw[7] << endl;
				    		DIF_DATARECORD << "A_B1R8["<< ma_tmp <<"]: " << A_B1R8[ma_tmp]	<< ", st2_Tw[8] = " << st2_Tw[8] << endl;
				    		DIF_DATARECORD << "A_B1R9["<< ma_tmp <<"]: " << A_B1R9[ma_tmp]	<< ", st2_Tw[9] = " << st2_Tw[9] << endl;
				    		DIF_DATARECORD << "A_B1R10["<< ma_tmp <<"]: "<< A_B1R10[ma_tmp]	<< ", st2_Tw[10] = " << st2_Tw[10] << endl;
				    		DIF_DATARECORD << "A_B1R11["<< ma_tmp <<"]: "<< A_B1R11[ma_tmp]	<< ", st2_Tw[11] = " << st2_Tw[11] << endl;
				    		DIF_DATARECORD << "A_B1R12["<< ma_tmp <<"]: "<< A_B1R12[ma_tmp]	<< ", st2_Tw[12] = " << st2_Tw[12] << endl;
				    		DIF_DATARECORD << "A_B1R13["<< ma_tmp <<"]: "<< A_B1R13[ma_tmp]	<< ", st2_Tw[13] = " << st2_Tw[13] << endl;
				    		DIF_DATARECORD << "A_B1R14["<< ma_tmp <<"]: "<< A_B1R14[ma_tmp]	<< ", st2_Tw[14] = " << st2_Tw[14] << endl;
				    		DIF_DATARECORD << "A_B1R15["<< ma_tmp <<"]: "<< A_B1R15[ma_tmp]	<< ", st2_Tw[15] = " << st2_Tw[15] << endl;					
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
							DIF_DATARECORD << "A_B1R0["<< ma_tmp <<"]: " << A_B1R0[ma_tmp]	<< ", st3_Tw[0] = " << 1 << endl;
				    		DIF_DATARECORD << "A_B1R1["<< ma_tmp <<"]: " << A_B1R1[ma_tmp]	<< ", st3_Tw[1] = " << 1 << endl;
				    		DIF_DATARECORD << "A_B1R2["<< ma_tmp <<"]: " << A_B1R2[ma_tmp]	<< ", st3_Tw[2] = " << 1 << endl;
				    		DIF_DATARECORD << "A_B1R3["<< ma_tmp <<"]: " << A_B1R3[ma_tmp]	<< ", st3_Tw[3] = " << 1 << endl;
				    		DIF_DATARECORD << "A_B1R4["<< ma_tmp <<"]: " << A_B1R4[ma_tmp]	<< ", st3_Tw[4] = " << 1 << endl;
				    		DIF_DATARECORD << "A_B1R5["<< ma_tmp <<"]: " << A_B1R5[ma_tmp]	<< ", st3_Tw[5] = " << 1 << endl;
				    		DIF_DATARECORD << "A_B1R6["<< ma_tmp <<"]: " << A_B1R6[ma_tmp]	<< ", st3_Tw[6] = " << 1 << endl;
				    		DIF_DATARECORD << "A_B1R7["<< ma_tmp <<"]: " << A_B1R7[ma_tmp]	<< ", st3_Tw[7] = " << 1 << endl;
				    		DIF_DATARECORD << "A_B1R8["<< ma_tmp <<"]: " << A_B1R8[ma_tmp]	<< ", st3_Tw[8] = " << 1 << endl;
				    		DIF_DATARECORD << "A_B1R9["<< ma_tmp <<"]: " << A_B1R9[ma_tmp]	<< ", st3_Tw[9] = " << 1 << endl;
				    		DIF_DATARECORD << "A_B1R10["<< ma_tmp <<"]: "<< A_B1R10[ma_tmp]	<< ", st3_Tw[10] = " << 1 << endl;
				    		DIF_DATARECORD << "A_B1R11["<< ma_tmp <<"]: "<< A_B1R11[ma_tmp]	<< ", st3_Tw[11] = " << 1 << endl;
				    		DIF_DATARECORD << "A_B1R12["<< ma_tmp <<"]: "<< A_B1R12[ma_tmp]	<< ", st3_Tw[12] = " << 1 << endl;
				    		DIF_DATARECORD << "A_B1R13["<< ma_tmp <<"]: "<< A_B1R13[ma_tmp]	<< ", st3_Tw[13] = " << 1 << endl;
				    		DIF_DATARECORD << "A_B1R14["<< ma_tmp <<"]: "<< A_B1R14[ma_tmp]	<< ", st3_Tw[14] = " << 1 << endl;
				    		DIF_DATARECORD << "A_B1R15["<< ma_tmp <<"]: "<< A_B1R15[ma_tmp]	<< ", st3_Tw[15] = " << 1 << endl;					
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
				    
					DIF_DATARECORD << "A_B1R0["<< ma_tmp <<"]: " << A_B1R0[ma_tmp]<<"\n";
				    DIF_DATARECORD << "A_B1R1["<< ma_tmp <<"]: " << A_B1R1[ma_tmp]<<"\n";
				    DIF_DATARECORD << "A_B1R2["<< ma_tmp <<"]: " << A_B1R2[ma_tmp]<<"\n";
				    DIF_DATARECORD << "A_B1R3["<< ma_tmp <<"]: " << A_B1R3[ma_tmp]<<"\n";
				    DIF_DATARECORD << "A_B1R4["<< ma_tmp <<"]: " << A_B1R4[ma_tmp]<<"\n";
				    DIF_DATARECORD << "A_B1R5["<< ma_tmp <<"]: " << A_B1R5[ma_tmp]<<"\n";
				    DIF_DATARECORD << "A_B1R6["<< ma_tmp <<"]: " << A_B1R6[ma_tmp]<<"\n";
				    DIF_DATARECORD << "A_B1R7["<< ma_tmp <<"]: " << A_B1R7[ma_tmp]<<"\n";
				    DIF_DATARECORD << "A_B1R8["<< ma_tmp <<"]: " << A_B1R8[ma_tmp]<<"\n";
				    DIF_DATARECORD << "A_B1R9["<< ma_tmp <<"]: " << A_B1R9[ma_tmp]<<"\n";
				    DIF_DATARECORD << "A_B1R10["<< ma_tmp <<"]: "<< A_B1R10[ma_tmp]<<"\n";
				    DIF_DATARECORD << "A_B1R11["<< ma_tmp <<"]: "<< A_B1R11[ma_tmp]<<"\n";
				    DIF_DATARECORD << "A_B1R12["<< ma_tmp <<"]: "<< A_B1R12[ma_tmp]<<"\n";
				    DIF_DATARECORD << "A_B1R13["<< ma_tmp <<"]: "<< A_B1R13[ma_tmp]<<"\n";
				    DIF_DATARECORD << "A_B1R14["<< ma_tmp <<"]: "<< A_B1R14[ma_tmp]<<"\n";
				    DIF_DATARECORD << "A_B1R15["<< ma_tmp <<"]: "<< A_B1R15[ma_tmp]<<"\n";
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
	
	//DIF_DATARECORD <<"----------------------------------------------------------------------------- \n";
	//DIF_DATARECORD <<" Radix-2 FFT computing start!!! \n";
	
	
	for(int i = 0;i < group; i++){
        for(int j=0; j < radix; j++){
        	gray_i  = Gray(i,group);
        	BC_tmp  = j * group + gray_i;
			//DIF_DATARECORD << "BC_tmp: " << BC_tmp << "\n";
        	BC = BC_tmp;
			//DIF_DATARECORD << "After RR_R16 , BC : " << BC << "\n";
        	AGU_R16(BC,bn_tmp,ma_tmp);
		    //DIF_DATARECORD << "BN : " << bn_tmp << "\n";
		    //DIF_DATARECORD << "MA : " << ma_tmp << "\n";			
        	if(bn_tmp == 0){
        		if(j < 2)bn0_bc_tmp = BC_tmp;
				
        		Radix2_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp]);
        		Radix2_BU(A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
        		Radix2_BU(A_B0R4[ma_tmp],A_B0R5[ma_tmp]);
        		Radix2_BU(A_B0R6[ma_tmp],A_B0R7[ma_tmp]);
        		Radix2_BU(A_B0R8[ma_tmp],A_B0R9[ma_tmp]);
        		Radix2_BU(A_B0R10[ma_tmp],A_B0R11[ma_tmp]);
        		Radix2_BU(A_B0R12[ma_tmp],A_B0R13[ma_tmp]);
        		Radix2_BU(A_B0R14[ma_tmp],A_B0R15[ma_tmp]);
				
        	}else {
				
        	    Radix2_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp]);
        	    Radix2_BU(A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
        	    Radix2_BU(A_B1R4[ma_tmp],A_B1R5[ma_tmp]);
        	    Radix2_BU(A_B1R6[ma_tmp],A_B1R7[ma_tmp]);
        	    Radix2_BU(A_B1R8[ma_tmp],A_B1R9[ma_tmp]);
        	    Radix2_BU(A_B1R10[ma_tmp],A_B1R11[ma_tmp]);
        	    Radix2_BU(A_B1R12[ma_tmp],A_B1R13[ma_tmp]);
        	    Radix2_BU(A_B1R14[ma_tmp],A_B1R15[ma_tmp]);
				
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
	
	//DIF_DATARECORD <<"***************************************************************\n";
	//DIF_DATARECORD <<"***** DATA OUTPUT!!                                         ***\n";
	//DIF_DATARECORD <<"***************************************************************\n";
	//data output
	// SPMB data output , output function is "NTT_REORDERINDEX_R16_R2_OUT".
	for(int i = 0;i < group;i++){
		for(int j = 0;j < radix;j++){
			gray_i  = Gray(i,group);
			BC_tmp  = j * group + gray_i;
			//DIF_DATARECORD <<"---------------------------------------------------\n";
			//DIF_DATARECORD <<" BC_tmp: "<< BC_tmp <<"\n";
			RR_R16_R2(BC_tmp,((4 * (Stage-1)) - 3),BC);
            BC_INDEX_32FLAG = BC >> 1;
			BC_INDEX_2FLAG  = BC % 2;

			AGU_R16(BC,bn_tmp,ma_tmp);
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

void DIF_NTTSPMB::DIF_NTT_r16_r4(std::vector<ZZ> &A,std::vector<ZZ> &B0R0,std::vector<ZZ> &B0R1,std::vector<ZZ> &B0R2,std::vector<ZZ> &B0R3,
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
	
	std::ofstream DIF_DATARECORD("./DIF_NTT_R16_R4_SPMB.txt");
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
	int radix_r2 = 4;
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
	////----------------------------------------
	
	//radix-16 Stage
    Stage_double  = log2(N);
    Stage       =  (int)floor(Stage_double/4);
	BC_WIDTH    =  (int)ceil(log2(N/16));
	offset      =  (int)N /  16;
	word_size   =  (int)N / (2 * 16);
	group       =  (int)N / (256);
	tw_modulus  =  (int)N /  16;
	bit_array_tmp.resize(BC_WIDTH);

    DIF_DATARECORD << "group: "    << group << "\n";
    DIF_DATARECORD << "BC_WIDTH: " << BC_WIDTH << "\n";

	
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
	if(display == 1)DIF_DATARECORD <<"radix-16 computing stage:  "<< Stage <<"\n";
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
		if(display == 1)DIF_DATARECORD <<"---------------------------------\n";
		if(display == 1)DIF_DATARECORD <<"Now Stage: "<< s <<"\n";		
		tw_modulus_tmp  = tw_modulus >> ( 4 * s);
		int tw_record = 0; // siang_record
		int cnt = 0;//siang
		for(int i = 0 ;i < group;i++){
			bn0_bc_tmp  = 0;
			bn1_bc_tmp  = 0;
			for(int j = 0;j < radix;j++){
				gray_i  = Gray(i,group);
			    BC_tmp  = j * group + gray_i;
				if(display == 1)DIF_DATARECORD <<"---------------------------------\n";
				if(display == 1)DIF_DATARECORD << "BC_tmp: " << BC_tmp << "\n";	
				if(s == Stage - 1) RR_R16_R4(BC_tmp,(4 * s - 2),BC);
				else RR_R16(BC_tmp,s,BC);
				if(display == 1)DIF_DATARECORD << "After RR_R16 , BC : " << BC << "\n";		
				length = BC % tw_modulus_tmp;
				PowerMod(factor_t,factor,length,p);
				AGU_R16(BC,bn_tmp,ma_tmp);
				if(display == 1)DIF_DATARECORD << "BN : " << bn_tmp << "\n";
				if(display == 1)DIF_DATARECORD << "MA : " << ma_tmp << "\n";	

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
				if(bn_tmp == 0){
					if(j < 2)bn0_bc_tmp = BC_tmp;
					switch(s){
						case 0:
							if(display == 1)DIF_DATARECORD << "A_B0R0["<< ma_tmp <<"]: " << A_B0R0[ma_tmp]	<< ", st0_Tw[0] = " << st0_Tw[0] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R1["<< ma_tmp <<"]: " << A_B0R1[ma_tmp]	<< ", st0_Tw[1] = " << st0_Tw[1] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R2["<< ma_tmp <<"]: " << A_B0R2[ma_tmp]	<< ", st0_Tw[2] = " << st0_Tw[2] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R3["<< ma_tmp <<"]: " << A_B0R3[ma_tmp]	<< ", st0_Tw[3] = " << st0_Tw[3] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R4["<< ma_tmp <<"]: " << A_B0R4[ma_tmp]	<< ", st0_Tw[4] = " << st0_Tw[4] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R5["<< ma_tmp <<"]: " << A_B0R5[ma_tmp]	<< ", st0_Tw[5] = " << st0_Tw[5] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R6["<< ma_tmp <<"]: " << A_B0R6[ma_tmp]	<< ", st0_Tw[6] = " << st0_Tw[6] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R7["<< ma_tmp <<"]: " << A_B0R7[ma_tmp]	<< ", st0_Tw[7] = " << st0_Tw[7] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R8["<< ma_tmp <<"]: " << A_B0R8[ma_tmp]	<< ", st0_Tw[8] = " << st0_Tw[8] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R9["<< ma_tmp <<"]: " << A_B0R9[ma_tmp]	<< ", st0_Tw[9] = " << st0_Tw[9] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R10["<< ma_tmp <<"]: "<< A_B0R10[ma_tmp]	<< ", st0_Tw[10] = " << st0_Tw[10] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R11["<< ma_tmp <<"]: "<< A_B0R11[ma_tmp]	<< ", st0_Tw[11] = " << st0_Tw[11] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R12["<< ma_tmp <<"]: "<< A_B0R12[ma_tmp]	<< ", st0_Tw[12] = " << st0_Tw[12] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R13["<< ma_tmp <<"]: "<< A_B0R13[ma_tmp]	<< ", st0_Tw[13] = " << st0_Tw[13] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R14["<< ma_tmp <<"]: "<< A_B0R14[ma_tmp]	<< ", st0_Tw[14] = " << st0_Tw[14] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R15["<< ma_tmp <<"]: "<< A_B0R15[ma_tmp]	<< ", st0_Tw[15] = " << st0_Tw[15] << endl;						
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
							if(display == 1)DIF_DATARECORD << "A_B0R0["<< ma_tmp <<"]: " << A_B0R0[ma_tmp]	<< ", st1_Tw[0] = " << st1_Tw[0] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R1["<< ma_tmp <<"]: " << A_B0R1[ma_tmp]	<< ", st1_Tw[1] = " << st1_Tw[1] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R2["<< ma_tmp <<"]: " << A_B0R2[ma_tmp]	<< ", st1_Tw[2] = " << st1_Tw[2] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R3["<< ma_tmp <<"]: " << A_B0R3[ma_tmp]	<< ", st1_Tw[3] = " << st1_Tw[3] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R4["<< ma_tmp <<"]: " << A_B0R4[ma_tmp]	<< ", st1_Tw[4] = " << st1_Tw[4] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R5["<< ma_tmp <<"]: " << A_B0R5[ma_tmp]	<< ", st1_Tw[5] = " << st1_Tw[5] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R6["<< ma_tmp <<"]: " << A_B0R6[ma_tmp]	<< ", st1_Tw[6] = " << st1_Tw[6] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R7["<< ma_tmp <<"]: " << A_B0R7[ma_tmp]	<< ", st1_Tw[7] = " << st1_Tw[7] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R8["<< ma_tmp <<"]: " << A_B0R8[ma_tmp]	<< ", st1_Tw[8] = " << st1_Tw[8] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R9["<< ma_tmp <<"]: " << A_B0R9[ma_tmp]	<< ", st1_Tw[9] = " << st1_Tw[9] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R10["<< ma_tmp <<"]: "<< A_B0R10[ma_tmp]	<< ", st1_Tw[10] = " << st1_Tw[10] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R11["<< ma_tmp <<"]: "<< A_B0R11[ma_tmp]	<< ", st1_Tw[11] = " << st1_Tw[11] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R12["<< ma_tmp <<"]: "<< A_B0R12[ma_tmp]	<< ", st1_Tw[12] = " << st1_Tw[12] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R13["<< ma_tmp <<"]: "<< A_B0R13[ma_tmp]	<< ", st1_Tw[13] = " << st1_Tw[13] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R14["<< ma_tmp <<"]: "<< A_B0R14[ma_tmp]	<< ", st1_Tw[14] = " << st1_Tw[14] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R15["<< ma_tmp <<"]: "<< A_B0R15[ma_tmp]	<< ", st1_Tw[15] = " << st1_Tw[15] << endl;						
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
							if(display == 1)DIF_DATARECORD << "A_B0R0["<< ma_tmp <<"]: " << A_B0R0[ma_tmp]	<< ", st2_Tw[0] = " << st2_Tw[0] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R1["<< ma_tmp <<"]: " << A_B0R1[ma_tmp]	<< ", st2_Tw[1] = " << st2_Tw[1] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R2["<< ma_tmp <<"]: " << A_B0R2[ma_tmp]	<< ", st2_Tw[2] = " << st2_Tw[2] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R3["<< ma_tmp <<"]: " << A_B0R3[ma_tmp]	<< ", st2_Tw[3] = " << st2_Tw[3] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R4["<< ma_tmp <<"]: " << A_B0R4[ma_tmp]	<< ", st2_Tw[4] = " << st2_Tw[4] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R5["<< ma_tmp <<"]: " << A_B0R5[ma_tmp]	<< ", st2_Tw[5] = " << st2_Tw[5] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R6["<< ma_tmp <<"]: " << A_B0R6[ma_tmp]	<< ", st2_Tw[6] = " << st2_Tw[6] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R7["<< ma_tmp <<"]: " << A_B0R7[ma_tmp]	<< ", st2_Tw[7] = " << st2_Tw[7] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R8["<< ma_tmp <<"]: " << A_B0R8[ma_tmp]	<< ", st2_Tw[8] = " << st2_Tw[8] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R9["<< ma_tmp <<"]: " << A_B0R9[ma_tmp]	<< ", st2_Tw[9] = " << st2_Tw[9] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R10["<< ma_tmp <<"]: "<< A_B0R10[ma_tmp]	<< ", st2_Tw[10] = " << st2_Tw[10] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R11["<< ma_tmp <<"]: "<< A_B0R11[ma_tmp]	<< ", st2_Tw[11] = " << st2_Tw[11] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R12["<< ma_tmp <<"]: "<< A_B0R12[ma_tmp]	<< ", st2_Tw[12] = " << st2_Tw[12] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R13["<< ma_tmp <<"]: "<< A_B0R13[ma_tmp]	<< ", st2_Tw[13] = " << st2_Tw[13] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R14["<< ma_tmp <<"]: "<< A_B0R14[ma_tmp]	<< ", st2_Tw[14] = " << st2_Tw[14] << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R15["<< ma_tmp <<"]: "<< A_B0R15[ma_tmp]	<< ", st2_Tw[15] = " << st2_Tw[15] << endl;						
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
							if(display == 1)DIF_DATARECORD << "A_B0R0["<< ma_tmp <<"]: " << A_B0R0[ma_tmp]	<< ", st3_Tw[0] = " << 1 << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R1["<< ma_tmp <<"]: " << A_B0R1[ma_tmp]	<< ", st3_Tw[1] = " << 1 << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R2["<< ma_tmp <<"]: " << A_B0R2[ma_tmp]	<< ", st3_Tw[2] = " << 1 << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R3["<< ma_tmp <<"]: " << A_B0R3[ma_tmp]	<< ", st3_Tw[3] = " << 1 << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R4["<< ma_tmp <<"]: " << A_B0R4[ma_tmp]	<< ", st3_Tw[4] = " << 1 << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R5["<< ma_tmp <<"]: " << A_B0R5[ma_tmp]	<< ", st3_Tw[5] = " << 1 << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R6["<< ma_tmp <<"]: " << A_B0R6[ma_tmp]	<< ", st3_Tw[6] = " << 1 << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R7["<< ma_tmp <<"]: " << A_B0R7[ma_tmp]	<< ", st3_Tw[7] = " << 1 << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R8["<< ma_tmp <<"]: " << A_B0R8[ma_tmp]	<< ", st3_Tw[8] = " << 1 << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R9["<< ma_tmp <<"]: " << A_B0R9[ma_tmp]	<< ", st3_Tw[9] = " << 1 << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R10["<< ma_tmp <<"]: "<< A_B0R10[ma_tmp]	<< ", st3_Tw[10] = " << 1 << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R11["<< ma_tmp <<"]: "<< A_B0R11[ma_tmp]	<< ", st3_Tw[11] = " << 1 << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R12["<< ma_tmp <<"]: "<< A_B0R12[ma_tmp]	<< ", st3_Tw[12] = " << 1 << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R13["<< ma_tmp <<"]: "<< A_B0R13[ma_tmp]	<< ", st3_Tw[13] = " << 1 << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R14["<< ma_tmp <<"]: "<< A_B0R14[ma_tmp]	<< ", st3_Tw[14] = " << 1 << endl;	
				    		if(display == 1)DIF_DATARECORD << "A_B0R15["<< ma_tmp <<"]: "<< A_B0R15[ma_tmp]	<< ", st3_Tw[15] = " << 1 << endl;						
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
					if(display == 1)DIF_DATARECORD <<" After mult by twiddle factor!!! \n";
					if(display == 1)DIF_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]<<"\n";
					if(display == 1)DIF_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]<<"\n";
					if(display == 1)DIF_DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]<<"\n";
					if(display == 1)DIF_DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]<<"\n";
					if(display == 1)DIF_DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]<<"\n";
					if(display == 1)DIF_DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]<<"\n";
					if(display == 1)DIF_DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]<<"\n";
					if(display == 1)DIF_DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]<<"\n";
					if(display == 1)DIF_DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]<<"\n";
					if(display == 1)DIF_DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]<<"\n";
					if(display == 1)DIF_DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";
					if(display == 1)DIF_DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";
					if(display == 1)DIF_DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";
					if(display == 1)DIF_DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";
					if(display == 1)DIF_DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";
					if(display == 1)DIF_DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";	
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
							if(display == 1)DIF_DATARECORD << "A_B1R0["<< ma_tmp <<"]: " << A_B1R0[ma_tmp]	<<	", st0_Tw[0] = " << st0_Tw[0] << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R1["<< ma_tmp <<"]: " << A_B1R1[ma_tmp]	<<	", st0_Tw[1] = " << st0_Tw[1] << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R2["<< ma_tmp <<"]: " << A_B1R2[ma_tmp]	<<	", st0_Tw[2] = " << st0_Tw[2] << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R3["<< ma_tmp <<"]: " << A_B1R3[ma_tmp]	<<	", st0_Tw[3] = " << st0_Tw[3] << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R4["<< ma_tmp <<"]: " << A_B1R4[ma_tmp]	<<	", st0_Tw[4] = " << st0_Tw[4] << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R5["<< ma_tmp <<"]: " << A_B1R5[ma_tmp]	<<	", st0_Tw[5] = " << st0_Tw[5] << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R6["<< ma_tmp <<"]: " << A_B1R6[ma_tmp]	<<	", st0_Tw[6] = " << st0_Tw[6] << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R7["<< ma_tmp <<"]: " << A_B1R7[ma_tmp]	<<	", st0_Tw[7] = " << st0_Tw[7] << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R8["<< ma_tmp <<"]: " << A_B1R8[ma_tmp]	<<	", st0_Tw[8] = " << st0_Tw[8] << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R9["<< ma_tmp <<"]: " << A_B1R9[ma_tmp]	<<	", st0_Tw[9] = " << st0_Tw[9] << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R10["<< ma_tmp <<"]: "<< A_B1R10[ma_tmp]	<<	", st0_Tw[10] = " << st0_Tw[10] << endl;
							if(display == 1)DIF_DATARECORD << "A_B1R11["<< ma_tmp <<"]: "<< A_B1R11[ma_tmp]	<<	", st0_Tw[11] = " << st0_Tw[11] << endl;
							if(display == 1)DIF_DATARECORD << "A_B1R12["<< ma_tmp <<"]: "<< A_B1R12[ma_tmp]	<<	", st0_Tw[12] = " << st0_Tw[12] << endl;
							if(display == 1)DIF_DATARECORD << "A_B1R13["<< ma_tmp <<"]: "<< A_B1R13[ma_tmp]	<<	", st0_Tw[13] = " << st0_Tw[13] << endl;
							if(display == 1)DIF_DATARECORD << "A_B1R14["<< ma_tmp <<"]: "<< A_B1R14[ma_tmp]	<<	", st0_Tw[14] = " << st0_Tw[14] << endl;
							if(display == 1)DIF_DATARECORD << "A_B1R15["<< ma_tmp <<"]: "<< A_B1R15[ma_tmp]	<<	", st0_Tw[15] = " << st0_Tw[15] << endl;				
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
							if(display == 1)DIF_DATARECORD << "A_B1R0["<< ma_tmp <<"]: " << A_B1R0[ma_tmp]	<<	", st1_Tw[0] = " << st1_Tw[0] << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R1["<< ma_tmp <<"]: " << A_B1R1[ma_tmp]	<<	", st1_Tw[1] = " << st1_Tw[1] << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R2["<< ma_tmp <<"]: " << A_B1R2[ma_tmp]	<<	", st1_Tw[2] = " << st1_Tw[2] << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R3["<< ma_tmp <<"]: " << A_B1R3[ma_tmp]	<<	", st1_Tw[3] = " << st1_Tw[3] << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R4["<< ma_tmp <<"]: " << A_B1R4[ma_tmp]	<<	", st1_Tw[4] = " << st1_Tw[4] << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R5["<< ma_tmp <<"]: " << A_B1R5[ma_tmp]	<<	", st1_Tw[5] = " << st1_Tw[5] << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R6["<< ma_tmp <<"]: " << A_B1R6[ma_tmp]	<<	", st1_Tw[6] = " << st1_Tw[6] << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R7["<< ma_tmp <<"]: " << A_B1R7[ma_tmp]	<<	", st1_Tw[7] = " << st1_Tw[7] << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R8["<< ma_tmp <<"]: " << A_B1R8[ma_tmp]	<<	", st1_Tw[8] = " << st1_Tw[8] << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R9["<< ma_tmp <<"]: " << A_B1R9[ma_tmp]	<<	", st1_Tw[9] = " << st1_Tw[9] << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R10["<< ma_tmp <<"]: "<< A_B1R10[ma_tmp]	<<	", st1_Tw[10] = " << st1_Tw[10] << endl;
							if(display == 1)DIF_DATARECORD << "A_B1R11["<< ma_tmp <<"]: "<< A_B1R11[ma_tmp]	<<	", st1_Tw[11] = " << st1_Tw[11] << endl;
							if(display == 1)DIF_DATARECORD << "A_B1R12["<< ma_tmp <<"]: "<< A_B1R12[ma_tmp]	<<	", st1_Tw[12] = " << st1_Tw[12] << endl;
							if(display == 1)DIF_DATARECORD << "A_B1R13["<< ma_tmp <<"]: "<< A_B1R13[ma_tmp]	<<	", st1_Tw[13] = " << st1_Tw[13] << endl;
							if(display == 1)DIF_DATARECORD << "A_B1R14["<< ma_tmp <<"]: "<< A_B1R14[ma_tmp]	<<	", st1_Tw[14] = " << st1_Tw[14] << endl;
							if(display == 1)DIF_DATARECORD << "A_B1R15["<< ma_tmp <<"]: "<< A_B1R15[ma_tmp]	<<	", st1_Tw[15] = " << st1_Tw[15] << endl;				
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
							if(display == 1)DIF_DATARECORD << "A_B1R0["<< ma_tmp <<"]: " << A_B1R0[ma_tmp]	<<	", st2_Tw[0] = " << st2_Tw[0] << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R1["<< ma_tmp <<"]: " << A_B1R1[ma_tmp]	<<	", st2_Tw[1] = " << st2_Tw[1] << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R2["<< ma_tmp <<"]: " << A_B1R2[ma_tmp]	<<	", st2_Tw[2] = " << st2_Tw[2] << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R3["<< ma_tmp <<"]: " << A_B1R3[ma_tmp]	<<	", st2_Tw[3] = " << st2_Tw[3] << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R4["<< ma_tmp <<"]: " << A_B1R4[ma_tmp]	<<	", st2_Tw[4] = " << st2_Tw[4] << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R5["<< ma_tmp <<"]: " << A_B1R5[ma_tmp]	<<	", st2_Tw[5] = " << st2_Tw[5] << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R6["<< ma_tmp <<"]: " << A_B1R6[ma_tmp]	<<	", st2_Tw[6] = " << st2_Tw[6] << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R7["<< ma_tmp <<"]: " << A_B1R7[ma_tmp]	<<	", st2_Tw[7] = " << st2_Tw[7] << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R8["<< ma_tmp <<"]: " << A_B1R8[ma_tmp]	<<	", st2_Tw[8] = " << st2_Tw[8] << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R9["<< ma_tmp <<"]: " << A_B1R9[ma_tmp]	<<	", st2_Tw[9] = " << st2_Tw[9] << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R10["<< ma_tmp <<"]: "<< A_B1R10[ma_tmp]	<<	", st2_Tw[10] = " << st2_Tw[10] << endl;
							if(display == 1)DIF_DATARECORD << "A_B1R11["<< ma_tmp <<"]: "<< A_B1R11[ma_tmp]	<<	", st2_Tw[11] = " << st2_Tw[11] << endl;
							if(display == 1)DIF_DATARECORD << "A_B1R12["<< ma_tmp <<"]: "<< A_B1R12[ma_tmp]	<<	", st2_Tw[12] = " << st2_Tw[12] << endl;
							if(display == 1)DIF_DATARECORD << "A_B1R13["<< ma_tmp <<"]: "<< A_B1R13[ma_tmp]	<<	", st2_Tw[13] = " << st2_Tw[13] << endl;
							if(display == 1)DIF_DATARECORD << "A_B1R14["<< ma_tmp <<"]: "<< A_B1R14[ma_tmp]	<<	", st2_Tw[14] = " << st2_Tw[14] << endl;
							if(display == 1)DIF_DATARECORD << "A_B1R15["<< ma_tmp <<"]: "<< A_B1R15[ma_tmp]	<<	", st2_Tw[15] = " << st2_Tw[15] << endl;				
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
							if(display == 1)DIF_DATARECORD << "A_B1R0["<< ma_tmp <<"]: " << A_B1R0[ma_tmp]	<<	", st3_Tw[0] = " << 1 << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R1["<< ma_tmp <<"]: " << A_B1R1[ma_tmp]	<<	", st3_Tw[1] = " << 1 << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R2["<< ma_tmp <<"]: " << A_B1R2[ma_tmp]	<<	", st3_Tw[2] = " << 1 << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R3["<< ma_tmp <<"]: " << A_B1R3[ma_tmp]	<<	", st3_Tw[3] = " << 1 << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R4["<< ma_tmp <<"]: " << A_B1R4[ma_tmp]	<<	", st3_Tw[4] = " << 1 << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R5["<< ma_tmp <<"]: " << A_B1R5[ma_tmp]	<<	", st3_Tw[5] = " << 1 << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R6["<< ma_tmp <<"]: " << A_B1R6[ma_tmp]	<<	", st3_Tw[6] = " << 1 << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R7["<< ma_tmp <<"]: " << A_B1R7[ma_tmp]	<<	", st3_Tw[7] = " << 1 << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R8["<< ma_tmp <<"]: " << A_B1R8[ma_tmp]	<<	", st3_Tw[8] = " << 1 << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R9["<< ma_tmp <<"]: " << A_B1R9[ma_tmp]	<<	", st3_Tw[9] = " << 1 << endl; 
							if(display == 1)DIF_DATARECORD << "A_B1R10["<< ma_tmp <<"]: "<< A_B1R10[ma_tmp]	<<	", st3_Tw[10] = " << 1 << endl;
							if(display == 1)DIF_DATARECORD << "A_B1R11["<< ma_tmp <<"]: "<< A_B1R11[ma_tmp]	<<	", st3_Tw[11] = " << 1 << endl;
							if(display == 1)DIF_DATARECORD << "A_B1R12["<< ma_tmp <<"]: "<< A_B1R12[ma_tmp]	<<	", st3_Tw[12] = " << 1 << endl;
							if(display == 1)DIF_DATARECORD << "A_B1R13["<< ma_tmp <<"]: "<< A_B1R13[ma_tmp]	<<	", st3_Tw[13] = " << 1 << endl;
							if(display == 1)DIF_DATARECORD << "A_B1R14["<< ma_tmp <<"]: "<< A_B1R14[ma_tmp]	<<	", st3_Tw[14] = " << 1 << endl;
							if(display == 1)DIF_DATARECORD << "A_B1R15["<< ma_tmp <<"]: "<< A_B1R15[ma_tmp]	<<	", st3_Tw[15] = " << 1 << endl;				
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
					
					if(display == 1)DIF_DATARECORD <<" After mult by twiddle factor!!! \n";
					if(display == 1)DIF_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]<<"\n";
					if(display == 1)DIF_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]<<"\n";
					if(display == 1)DIF_DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]<<"\n";
					if(display == 1)DIF_DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]<<"\n";
					if(display == 1)DIF_DATARECORD <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]<<"\n";
					if(display == 1)DIF_DATARECORD <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]<<"\n";
					if(display == 1)DIF_DATARECORD <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]<<"\n";
					if(display == 1)DIF_DATARECORD <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]<<"\n";
					if(display == 1)DIF_DATARECORD <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]<<"\n";
					if(display == 1)DIF_DATARECORD <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]<<"\n";
					if(display == 1)DIF_DATARECORD <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";
					if(display == 1)DIF_DATARECORD <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";
					if(display == 1)DIF_DATARECORD <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";
					if(display == 1)DIF_DATARECORD <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";
					if(display == 1)DIF_DATARECORD <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";
					if(display == 1)DIF_DATARECORD <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";		
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
	
	std::cout << "radix-16 FFT computing over!!\n";
	
	
	for(int i = 0;i < group; i++){
        for(int j=0; j < radix; j++){
        	gray_i  = Gray(i,group);
        	BC_tmp  = j * group + gray_i;
        	RR_R16_R4(BC_tmp,((4 * (Stage-1)) - 2),BC);
        	AGU_R16(BC,bn_tmp,ma_tmp);
        	if(bn_tmp == 0){
				if(display == 1)DIF_DATARECORD <<" Before Radix-4 Butterfly unit computing!!! \n";
				if(display == 1)DIF_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";					
                Radix4_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp],A_B0R2[ma_tmp],A_B0R3[ma_tmp]);
        		Radix4_BU(A_B0R4[ma_tmp],A_B0R5[ma_tmp],A_B0R6[ma_tmp],A_B0R7[ma_tmp]);
        		Radix4_BU(A_B0R8[ma_tmp],A_B0R9[ma_tmp],A_B0R10[ma_tmp],A_B0R11[ma_tmp]);
        		Radix4_BU(A_B0R12[ma_tmp],A_B0R13[ma_tmp],A_B0R14[ma_tmp],A_B0R15[ma_tmp]);
				if(display == 1)DIF_DATARECORD <<" After Radix-4 Butterfly unit computing!!! \n";
				if(display == 1)DIF_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";			
        	}else {
				if(display == 1)DIF_DATARECORD <<" Before Radix-4 Butterfly unit computing!!! \n";
				if(display == 1)DIF_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";			
        	    Radix4_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp],A_B1R2[ma_tmp],A_B1R3[ma_tmp]);
        	    Radix4_BU(A_B1R4[ma_tmp],A_B1R5[ma_tmp],A_B1R6[ma_tmp],A_B1R7[ma_tmp]);
        	    Radix4_BU(A_B1R8[ma_tmp],A_B1R9[ma_tmp],A_B1R10[ma_tmp],A_B1R11[ma_tmp]);
        	    Radix4_BU(A_B1R12[ma_tmp],A_B1R13[ma_tmp],A_B1R14[ma_tmp],A_B1R15[ma_tmp]);
				if(display == 1)DIF_DATARECORD <<" After Radix-4 Butterfly unit computing!!! \n";
				if(display == 1)DIF_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";
				if(display == 1)DIF_DATARECORD <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";			
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

void DIF_NTTSPMB::DIF_NTT_r16_r8(std::vector<ZZ> &A,std::vector<ZZ> &B0R0,std::vector<ZZ> &B0R1,std::vector<ZZ> &B0R2,std::vector<ZZ> &B0R3,
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
	std::ofstream DIF_DATARECORD("./NTT_R16_R8_SPMB.txt");
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
	int radix_r2 = 8;
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
	////----------------------------------------


	//radix-16 Stage
    Stage_double  = log2(N);
    Stage       =  (int)floor(Stage_double/4);
	BC_WIDTH    =  (int)ceil(log2(N/16));
	offset      =  (int)N /  16;
	word_size   =  (int)N / (2 * 16);
	group       =  (int)N / (256);
	tw_modulus  =  (int)N /  16;
	bit_array_tmp.resize(BC_WIDTH);

    if(display == 1)DIF_DATARECORD << "group: "    << group << "\n";
    if(display == 1)DIF_DATARECORD << "BC_WIDTH: " << BC_WIDTH << "\n";

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
	if(display == 1)DIF_DATARECORD <<"Stage: " << Stage << "\n";	
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
		if(display == 1)DIF_DATARECORD <<"**************************************\n";										
		tw_modulus_tmp  = tw_modulus >> ( 4 * s);
		for(int i = 0 ;i < group;i++){
			bn0_bc_tmp  = 0;
			bn1_bc_tmp  = 0;
			if(display == 1)DIF_DATARECORD <<"-------------------------------------\n";										
			if(display == 1)DIF_DATARECORD <<"NOW stage: " << s << "!!!!\n";							
			for(int j = 0;j < radix;j++){
				gray_i  = Gray(i,group);
			    BC_tmp  = j * group + gray_i;
				if(display == 1)DIF_DATARECORD <<"---------------------------------\n";
				if(display == 1)DIF_DATARECORD << "i = " << i << ", j = " << j << std::endl;
				if(display == 1)DIF_DATARECORD << "BC_tmp: " << BC_tmp << "\n";					
				if(s == Stage - 1) RR_R16_R8(BC_tmp,(4 * s - 1),BC);
				else RR_R16(BC_tmp,s,BC);
				if(display == 1)DIF_DATARECORD << "After RR_R16 , BC : " << BC << "\n";		
				length = BC % tw_modulus_tmp;
				PowerMod(factor_t,factor,length,p);
				AGU_R16(BC,bn_tmp,ma_tmp);
				if(display == 1)DIF_DATARECORD << "BN : " << bn_tmp << "\n";
				if(display == 1)DIF_DATARECORD << "MA : " << ma_tmp << "\n";	

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
				if(bn_tmp == 0){
					if(j < 2)bn0_bc_tmp = BC_tmp;
					switch(s){
						case 0:
                    		if(display == 1)DIF_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]	 << ", st0_Tw[0] = " << st0_Tw[0] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]	 << ", st0_Tw[1] = " << st0_Tw[1] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]	 << ", st0_Tw[2] = " << st0_Tw[2] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]	 << ", st0_Tw[3] = " << st0_Tw[3] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]	 << ", st0_Tw[4] = " << st0_Tw[4] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]	 << ", st0_Tw[5] = " << st0_Tw[5] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]	 << ", st0_Tw[6] = " << st0_Tw[6] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]	 << ", st0_Tw[7] = " << st0_Tw[7] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]	 << ", st0_Tw[8] = " << st0_Tw[8] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]	 << ", st0_Tw[9] = " << st0_Tw[9] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp] << ", st0_Tw[10] = " << st0_Tw[10] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp] << ", st0_Tw[11] = " << st0_Tw[11] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp] << ", st0_Tw[12] = " << st0_Tw[12] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp] << ", st0_Tw[13] = " << st0_Tw[13] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp] << ", st0_Tw[14] = " << st0_Tw[14] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp] << ", st0_Tw[15] = " << st0_Tw[15] << endl;			
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
							if(display == 1)DIF_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]	 << ", st1_Tw[0] = " << st1_Tw[0] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]	 << ", st1_Tw[1] = " << st1_Tw[1] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]	 << ", st1_Tw[2] = " << st1_Tw[2] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]	 << ", st1_Tw[3] = " << st1_Tw[3] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]	 << ", st1_Tw[4] = " << st1_Tw[4] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]	 << ", st1_Tw[5] = " << st1_Tw[5] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]	 << ", st1_Tw[6] = " << st1_Tw[6] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]	 << ", st1_Tw[7] = " << st1_Tw[7] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]	 << ", st1_Tw[8] = " << st1_Tw[8] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]	 << ", st1_Tw[9] = " << st1_Tw[9] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp] << ", st1_Tw[10] = " << st1_Tw[10] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp] << ", st1_Tw[11] = " << st1_Tw[11] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp] << ", st1_Tw[12] = " << st1_Tw[12] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp] << ", st1_Tw[13] = " << st1_Tw[13] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp] << ", st1_Tw[14] = " << st1_Tw[14] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp] << ", st1_Tw[15] = " << st1_Tw[15] << endl;			
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
							if(display == 1)DIF_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]	 << ", st2_Tw[0] = " << st2_Tw[0] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]	 << ", st2_Tw[1] = " << st2_Tw[1] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]	 << ", st2_Tw[2] = " << st2_Tw[2] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]	 << ", st2_Tw[3] = " << st2_Tw[3] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]	 << ", st2_Tw[4] = " << st2_Tw[4] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]	 << ", st2_Tw[5] = " << st2_Tw[5] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]	 << ", st2_Tw[6] = " << st2_Tw[6] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]	 << ", st2_Tw[7] = " << st2_Tw[7] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]	 << ", st2_Tw[8] = " << st2_Tw[8] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]	 << ", st2_Tw[9] = " << st2_Tw[9] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp] << ", st2_Tw[10] = " << st2_Tw[10] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp] << ", st2_Tw[11] = " << st2_Tw[11] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp] << ", st2_Tw[12] = " << st2_Tw[12] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp] << ", st2_Tw[13] = " << st2_Tw[13] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp] << ", st2_Tw[14] = " << st2_Tw[14] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp] << ", st2_Tw[15] = " << st2_Tw[15] << endl;			
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
							if(display == 1)DIF_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]	 << ", st3_Tw[0] = " << 1 << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]	 << ", st3_Tw[1] = " << 1 << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]	 << ", st3_Tw[2] = " << 1 << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]	 << ", st3_Tw[3] = " << 1 << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]	 << ", st3_Tw[4] = " << 1 << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]	 << ", st3_Tw[5] = " << 1 << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]	 << ", st3_Tw[6] = " << 1 << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]	 << ", st3_Tw[7] = " << 1 << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]	 << ", st3_Tw[8] = " << 1 << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]	 << ", st3_Tw[9] = " << 1 << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp] << ", st3_Tw[10] = " << 1 << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp] << ", st3_Tw[11] = " << 1 << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp] << ", st3_Tw[12] = " << 1 << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp] << ", st3_Tw[13] = " << 1 << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp] << ", st3_Tw[14] = " << 1 << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp] << ", st3_Tw[15] = " << 1 << endl;			
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
                    
                    if(display == 1)DIF_DATARECORD <<" After multiplication!!!\n";					
                    if(display == 1)DIF_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<< A_B0R0[ma_tmp]<<"\n";					
                    if(display == 1)DIF_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<< A_B0R1[ma_tmp]<<"\n";					
                    if(display == 1)DIF_DATARECORD <<"A_B0R2["<<ma_tmp<<"]: "<< A_B0R2[ma_tmp]<<"\n";					
                    if(display == 1)DIF_DATARECORD <<"A_B0R3["<<ma_tmp<<"]: "<< A_B0R3[ma_tmp]<<"\n";					
                    if(display == 1)DIF_DATARECORD <<"A_B0R4["<<ma_tmp<<"]: "<< A_B0R4[ma_tmp]<<"\n";					
                    if(display == 1)DIF_DATARECORD <<"A_B0R5["<<ma_tmp<<"]: "<< A_B0R5[ma_tmp]<<"\n";					
                    if(display == 1)DIF_DATARECORD <<"A_B0R6["<<ma_tmp<<"]: "<< A_B0R6[ma_tmp]<<"\n";					
                    if(display == 1)DIF_DATARECORD <<"A_B0R7["<<ma_tmp<<"]: "<< A_B0R7[ma_tmp]<<"\n";					
                    if(display == 1)DIF_DATARECORD <<"A_B0R8["<<ma_tmp<<"]: "<< A_B0R8[ma_tmp]<<"\n";					
                    if(display == 1)DIF_DATARECORD <<"A_B0R9["<<ma_tmp<<"]: "<< A_B0R9[ma_tmp]<<"\n";					
                    if(display == 1)DIF_DATARECORD <<"A_B0R10["<<ma_tmp<<"]: "<< A_B0R10[ma_tmp]<<"\n";					
                    if(display == 1)DIF_DATARECORD <<"A_B0R11["<<ma_tmp<<"]: "<< A_B0R11[ma_tmp]<<"\n";					
                    if(display == 1)DIF_DATARECORD <<"A_B0R12["<<ma_tmp<<"]: "<< A_B0R12[ma_tmp]<<"\n";					
                    if(display == 1)DIF_DATARECORD <<"A_B0R13["<<ma_tmp<<"]: "<< A_B0R13[ma_tmp]<<"\n";					
                    if(display == 1)DIF_DATARECORD <<"A_B0R14["<<ma_tmp<<"]: "<< A_B0R14[ma_tmp]<<"\n";					
                    if(display == 1)DIF_DATARECORD <<"A_B0R15["<<ma_tmp<<"]: "<< A_B0R15[ma_tmp]<<"\n";											
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
                    		if(display == 1)DIF_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]	 << ", st0_Tw[0] = " << st0_Tw[0] << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]	 << ", st0_Tw[1] = " << st0_Tw[1] << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]	 << ", st0_Tw[2] = " << st0_Tw[2] << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]	 << ", st0_Tw[3] = " << st0_Tw[3] << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]	 << ", st0_Tw[4] = " << st0_Tw[4] << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]	 << ", st0_Tw[5] = " << st0_Tw[5] << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]	 << ", st0_Tw[6] = " << st0_Tw[6] << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]	 << ", st0_Tw[7] = " << st0_Tw[7] << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]	 << ", st0_Tw[8] = " << st0_Tw[8] << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]	 << ", st0_Tw[9] = " << st0_Tw[9] << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp] << ", st0_Tw[10] = " << st0_Tw[10] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp] << ", st0_Tw[11] = " << st0_Tw[11] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp] << ", st0_Tw[12] = " << st0_Tw[12] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp] << ", st0_Tw[13] = " << st0_Tw[13] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp] << ", st0_Tw[14] = " << st0_Tw[14] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp] << ", st0_Tw[15] = " << st0_Tw[15] << endl;						
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
							if(display == 1)DIF_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]	 << ", st1_Tw[0] = " << st1_Tw[0] << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]	 << ", st1_Tw[1] = " << st1_Tw[1] << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]	 << ", st1_Tw[2] = " << st1_Tw[2] << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]	 << ", st1_Tw[3] = " << st1_Tw[3] << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]	 << ", st1_Tw[4] = " << st1_Tw[4] << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]	 << ", st1_Tw[5] = " << st1_Tw[5] << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]	 << ", st1_Tw[6] = " << st1_Tw[6] << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]	 << ", st1_Tw[7] = " << st1_Tw[7] << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]	 << ", st1_Tw[8] = " << st1_Tw[8] << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]	 << ", st1_Tw[9] = " << st1_Tw[9] << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp] << ", st1_Tw[10] = " << st1_Tw[10] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp] << ", st1_Tw[11] = " << st1_Tw[11] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp] << ", st1_Tw[12] = " << st1_Tw[12] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp] << ", st1_Tw[13] = " << st1_Tw[13] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp] << ", st1_Tw[14] = " << st1_Tw[14] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp] << ", st1_Tw[15] = " << st1_Tw[15] << endl;						
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
							if(display == 1)DIF_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]	 << ", st2_Tw[0] = " << st2_Tw[0] << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]	 << ", st2_Tw[1] = " << st2_Tw[1] << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]	 << ", st2_Tw[2] = " << st2_Tw[2] << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]	 << ", st2_Tw[3] = " << st2_Tw[3] << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]	 << ", st2_Tw[4] = " << st2_Tw[4] << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]	 << ", st2_Tw[5] = " << st2_Tw[5] << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]	 << ", st2_Tw[6] = " << st2_Tw[6] << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]	 << ", st2_Tw[7] = " << st2_Tw[7] << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]	 << ", st2_Tw[8] = " << st2_Tw[8] << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]	 << ", st2_Tw[9] = " << st2_Tw[9] << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp] << ", st2_Tw[10] = " << st2_Tw[10] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp] << ", st2_Tw[11] = " << st2_Tw[11] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp] << ", st2_Tw[12] = " << st2_Tw[12] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp] << ", st2_Tw[13] = " << st2_Tw[13] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp] << ", st2_Tw[14] = " << st2_Tw[14] << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp] << ", st2_Tw[15] = " << st2_Tw[15] << endl;						
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
							if(display == 1)DIF_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]	 << ", st3_Tw[0] = " << 1 << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]	 << ", st3_Tw[1] = " << 1 << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]	 << ", st3_Tw[2] = " << 1 << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]	 << ", st3_Tw[3] = " << 1 << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]	 << ", st3_Tw[4] = " << 1 << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]	 << ", st3_Tw[5] = " << 1 << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]	 << ", st3_Tw[6] = " << 1 << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]	 << ", st3_Tw[7] = " << 1 << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]	 << ", st3_Tw[8] = " << 1 << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]	 << ", st3_Tw[9] = " << 1 << endl;				
                    		if(display == 1)DIF_DATARECORD <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp] << ", st3_Tw[10] = " << 1 << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp] << ", st3_Tw[11] = " << 1 << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp] << ", st3_Tw[12] = " << 1 << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp] << ", st3_Tw[13] = " << 1 << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp] << ", st3_Tw[14] = " << 1 << endl;					
                    		if(display == 1)DIF_DATARECORD <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp] << ", st3_Tw[15] = " << 1 << endl;						
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
        
                    if(display == 1)DIF_DATARECORD <<" After multiplication!!!\n";					
                    if(display == 1)DIF_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<< A_B1R0[ma_tmp]<<"\n";					
                    if(display == 1)DIF_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<< A_B1R1[ma_tmp]<<"\n";					
                    if(display == 1)DIF_DATARECORD <<"A_B1R2["<<ma_tmp<<"]: "<< A_B1R2[ma_tmp]<<"\n";					
                    if(display == 1)DIF_DATARECORD <<"A_B1R3["<<ma_tmp<<"]: "<< A_B1R3[ma_tmp]<<"\n";					
                    if(display == 1)DIF_DATARECORD <<"A_B1R4["<<ma_tmp<<"]: "<< A_B1R4[ma_tmp]<<"\n";					
                    if(display == 1)DIF_DATARECORD <<"A_B1R5["<<ma_tmp<<"]: "<< A_B1R5[ma_tmp]<<"\n";					
                    if(display == 1)DIF_DATARECORD <<"A_B1R6["<<ma_tmp<<"]: "<< A_B1R6[ma_tmp]<<"\n";					
                    if(display == 1)DIF_DATARECORD <<"A_B1R7["<<ma_tmp<<"]: "<< A_B1R7[ma_tmp]<<"\n";					
                    if(display == 1)DIF_DATARECORD <<"A_B1R8["<<ma_tmp<<"]: "<< A_B1R8[ma_tmp]<<"\n";					
                    if(display == 1)DIF_DATARECORD <<"A_B1R9["<<ma_tmp<<"]: "<< A_B1R9[ma_tmp]<<"\n";					
                    if(display == 1)DIF_DATARECORD <<"A_B1R10["<<ma_tmp<<"]: "<< A_B1R10[ma_tmp]<<"\n";					
                    if(display == 1)DIF_DATARECORD <<"A_B1R11["<<ma_tmp<<"]: "<< A_B1R11[ma_tmp]<<"\n";					
                    if(display == 1)DIF_DATARECORD <<"A_B1R12["<<ma_tmp<<"]: "<< A_B1R12[ma_tmp]<<"\n";					
                    if(display == 1)DIF_DATARECORD <<"A_B1R13["<<ma_tmp<<"]: "<< A_B1R13[ma_tmp]<<"\n";					
                    if(display == 1)DIF_DATARECORD <<"A_B1R14["<<ma_tmp<<"]: "<< A_B1R14[ma_tmp]<<"\n";					
                    if(display == 1)DIF_DATARECORD <<"A_B1R15["<<ma_tmp<<"]: "<< A_B1R15[ma_tmp]<<"\n";						
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
        	if(bn_tmp == 0){
                Radix8_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp],A_B0R2[ma_tmp],A_B0R3[ma_tmp],A_B0R4[ma_tmp],A_B0R5[ma_tmp],A_B0R6[ma_tmp],A_B0R7[ma_tmp]);
        		Radix8_BU(A_B0R8[ma_tmp],A_B0R9[ma_tmp],A_B0R10[ma_tmp],A_B0R11[ma_tmp],A_B0R12[ma_tmp],A_B0R13[ma_tmp],A_B0R14[ma_tmp],A_B0R15[ma_tmp]);
        	}else {
        	    Radix8_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp],A_B1R2[ma_tmp],A_B1R3[ma_tmp],A_B1R4[ma_tmp],A_B1R5[ma_tmp],A_B1R6[ma_tmp],A_B1R7[ma_tmp]);
        	    Radix8_BU(A_B1R8[ma_tmp],A_B1R9[ma_tmp],A_B1R10[ma_tmp],A_B1R11[ma_tmp],A_B1R12[ma_tmp],A_B1R13[ma_tmp],A_B1R14[ma_tmp],A_B1R15[ma_tmp]);
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
	
	DIF_DATARECORD <<"***************************************************************\n";
	DIF_DATARECORD <<"***** DATA OUTPUT!!                                         ***\n";
	DIF_DATARECORD <<"***************************************************************\n";
	//data output
	// SPMB data output , output function is "NTT_REORDERINDEX_R16_R8_OUT".
	for(int i = 0;i < group;i++){
		for(int j = 0;j < radix;j++){
			//gray_i  = Gray(i,group);
			BC_tmp  = j * group + i;
			if(display == 1)DIF_DATARECORD <<"---------------------------------------------------\n";
			if(display == 1)DIF_DATARECORD <<" BC_tmp: "<< BC_tmp <<"\n";
			RR_R16_R8(BC_tmp,((4 * (Stage-1)) - 1),BC);
			if(display == 1)DIF_DATARECORD <<" BC: "<< BC <<"\n";
			index_128_flag   = BC >> 3;
			index_3_BITS_FLAG = BC % 8;
			AGU_R16(BC,bn_tmp,ma_tmp);
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