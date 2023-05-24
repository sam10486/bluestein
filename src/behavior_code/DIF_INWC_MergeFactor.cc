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

void DIF_INWC::DIF_INWC_MergeFactor_radix2(std::vector<ZZ> &A){
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

	std::ofstream INWC_DATARECORD("./NWC_PrintData/INWC_Merge_R2_SPMB.txt");
	std::ofstream INWC_radix2("./NWC_PrintData/INWC_Merge_r2.txt");
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
	//DTFAG.DTFAG_ROM_init(
    //    radix_r1, radix_r2, fft_twiddle, fft_prime, debug,
    //    ROM0, ROM1, ROM2);
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

        //-----------for INWC----------------
        DTFAG.DTFAG_INWC_MergeFactor_ROM_init(
        radix_r1, radix_r2, fft_twiddle, fft_prime, debug, s, InvPhi,
        ROM0, ROM1, ROM2);
        //-----------------------------------

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
                ZZ stage3_factor;
                if (s == 3)
                {
                    stage3_factor = PowerMod(InvPhi, InvPhi_deg, p);
                    MulMod(stage3_factor, stage3_factor, InvTwo, p);
                }
                
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
					        //if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],InvTwo,p);
							// down
							if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],st0_Tw[1],p);
							//if(!debug) MulMod(InvPhi_dot_IW, InvPhi_Order, st0_Tw[1], p);
							//if(!debug) MulMod(InvPhi_dot_IW_dot_InvTwo, InvPhi_dot_IW, InvTwo, p);
							//if(!debug) MulMod(A_B0R1[ma_tmp], A_B0R1[ma_tmp], InvPhi_dot_IW_dot_InvTwo, p);
                            break;
                        case 1:
                            INWC_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp] << ", st1_Tw[0] = " << st1_Tw[0] << endl;
					        INWC_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp] << ", st1_Tw[1] = " << st1_Tw[1] << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_Order = " << InvPhi_Order << endl;
							if(debug) INWC_DATARECORD << "InvPhi_Order = " << InvPhi_deg << endl;
					        if(!debug) INWC_Radix2_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp]);
							// upper
                            if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],st1_Tw[0],p);
							//if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],InvTwo,p);
							// down
					        if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],st1_Tw[1],p);
							//if(!debug) MulMod(InvPhi_dot_IW, InvPhi_Order, st1_Tw[1], p);
							//if(!debug) MulMod(InvPhi_dot_IW_dot_InvTwo, InvPhi_dot_IW, InvTwo, p);
							//if(!debug) MulMod(A_B0R1[ma_tmp], A_B0R1[ma_tmp], InvPhi_dot_IW_dot_InvTwo, p);
                            break;
                        case 2:
                            INWC_DATARECORD <<"A_B0R0["<<ma_tmp<<"]: "<<A_B0R0[ma_tmp] << ", st2_Tw[0] = " << st2_Tw[0] << endl;
					        INWC_DATARECORD <<"A_B0R1["<<ma_tmp<<"]: "<<A_B0R1[ma_tmp] << ", st2_Tw[1] = " << st2_Tw[1] << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_Order = " << InvPhi_Order << endl;
							if(debug) INWC_DATARECORD << "InvPhi_Order = " << InvPhi_deg << endl;
					        if(!debug) INWC_Radix2_BU(A_B0R0[ma_tmp],A_B0R1[ma_tmp]);
							// upper
                            if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],st2_Tw[0],p);
							//if(!debug) MulMod(A_B0R0[ma_tmp],A_B0R0[ma_tmp],InvTwo,p);
							// down
					        if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],st2_Tw[1],p);
							//if(!debug) MulMod(InvPhi_dot_IW, InvPhi_Order, st2_Tw[1], p);
							//if(!debug) MulMod(InvPhi_dot_IW_dot_InvTwo, InvPhi_dot_IW, InvTwo, p);
							//if(!debug) MulMod(A_B0R1[ma_tmp], A_B0R1[ma_tmp], InvPhi_dot_IW_dot_InvTwo, p);
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
					        if(!debug) MulMod(A_B0R1[ma_tmp],A_B0R1[ma_tmp],stage3_factor,p);	
							//if(!debug) MulMod(InvPhi_dot_IW, InvPhi_Order, 1, p);
							//if(!debug) MulMod(InvPhi_dot_IW_dot_InvTwo, InvPhi_dot_IW, InvTwo, p);
							//if(!debug) MulMod(A_B0R1[ma_tmp], A_B0R1[ma_tmp], InvPhi_dot_IW_dot_InvTwo, p);	
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
							//if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],InvTwo,p);
							// down
					        if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],st0_Tw[1],p);  
							//if(!debug) MulMod(InvPhi_dot_IW, InvPhi_Order, st0_Tw[1], p);
							//if(!debug) MulMod(InvPhi_dot_IW_dot_InvTwo, InvPhi_dot_IW, InvTwo, p);
							//if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],InvPhi_dot_IW_dot_InvTwo,p);
                            break;
                        case 1:
                            INWC_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<< ", st1_Tw[0] = " << st1_Tw[0] << endl;
					        INWC_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<< ", st1_Tw[1] = " << st1_Tw[1] << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_Order = " << InvPhi_Order << endl;
							if(debug) INWC_DATARECORD << "InvPhi_Order = " << InvPhi_deg << endl;
					        if(!debug) INWC_Radix2_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp]);
							// upper
                            if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],st1_Tw[0],p);
							//if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],InvTwo,p);
							// down
					        if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],st1_Tw[1],p);
							//if(!debug) MulMod(InvPhi_dot_IW, InvPhi_Order, st1_Tw[1], p);
							//if(!debug) MulMod(InvPhi_dot_IW_dot_InvTwo, InvPhi_dot_IW, InvTwo, p);
							//if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],InvPhi_dot_IW_dot_InvTwo,p);
                            break;
                        case 2:
                            INWC_DATARECORD <<"A_B1R0["<<ma_tmp<<"]: "<<A_B1R0[ma_tmp]<< ", st2_Tw[0] = " << st2_Tw[0] << endl;
					        INWC_DATARECORD <<"A_B1R1["<<ma_tmp<<"]: "<<A_B1R1[ma_tmp]<< ", st2_Tw[1] = " << st2_Tw[1] << endl;
							if(!debug) INWC_DATARECORD << "InvPhi_Order = " << InvPhi_Order << endl;
							if(debug) INWC_DATARECORD << "InvPhi_Order = " << InvPhi_deg << endl;
					        if(!debug) INWC_Radix2_BU(A_B1R0[ma_tmp],A_B1R1[ma_tmp]);
							// upper
                            if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],st2_Tw[0],p);
							//if(!debug) MulMod(A_B1R0[ma_tmp],A_B1R0[ma_tmp],InvTwo,p);
							// down
					        if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],st2_Tw[1],p);
							//if(!debug) MulMod(InvPhi_dot_IW, InvPhi_Order, st2_Tw[1], p);
							//if(!debug) MulMod(InvPhi_dot_IW_dot_InvTwo, InvPhi_dot_IW, InvTwo, p);
							//if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],InvPhi_dot_IW_dot_InvTwo,p);
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
					        if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],stage3_factor,p);
							//if(!debug) MulMod(InvPhi_dot_IW, InvPhi_Order, 1, p);
							//if(!debug) MulMod(InvPhi_dot_IW_dot_InvTwo, InvPhi_dot_IW, InvTwo, p);
							//if(!debug) MulMod(A_B1R1[ma_tmp],A_B1R1[ma_tmp],InvPhi_dot_IW_dot_InvTwo,p);
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
