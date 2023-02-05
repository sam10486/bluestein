#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <NTL/ZZ.h>

#include "SPMB.h"

using namespace NTL;
using namespace std;


void SPMB::init(unsigned long fft_p , unsigned long r,int bc_w,unsigned long CP_w,
	 ZZ cyclotomic_p,long m_th,long phi_m){
     double radix_bit;
	 double FFTP_doub;
	 double stage_ceil;
	 double stage_floor;
	 unsigned long Mixed_radix;
	 radix_bit   = log2(r);
	 FFTP_doub   = log2(fft_p);
	 stage_ceil  = ceil(FFTP_doub/radix_bit);
	 stage_floor = floor(FFTP_doub/radix_bit);
	 if(stage_ceil == stage_floor) IsMixed = 0;
	 else IsMixed = 1;
	 std::ofstream parameter_o("./fft_parameter.txt");
	 parameter_o <<"---------------------------------\n";
	 parameter_o <<"fft_point: " << fft_p <<"\n";
	 parameter_o <<"radix: " << r <<"\n";
	 parameter_o <<"butterfly counter width: " << bc_w <<"\n";
	 parameter_o <<"cyclotomic prime width: " << CP_w <<"\n";
	 parameter_o <<"cyclotomic prime : " << cyclotomic_p <<"\n";
	 parameter_o <<"m-th : " << m_th <<"\n";
	 parameter_o <<"phi_m : " << phi_m <<"\n";
	 
	 Mixed_radix = fft_p;
	 while((Mixed_radix % 16) == 0){
		 Mixed_radix = (unsigned long) Mixed_radix / 16; 
	 }
	 
	 std::cout << "--------------------------------------\n";
	 std::cout << "IsMixed: "     << IsMixed <<"\n";
	 std::cout << "Mixed_radix: " << Mixed_radix <<"\n";
	 
	 if(r==4){
	    if(IsMixed == 0)init_r4(fft_p,r,bc_w,CP_w,cyclotomic_p,m_th,phi_m);
		else init_r4_r2(fft_p,r,bc_w,CP_w,cyclotomic_p,m_th,phi_m);
	 }
	 if(r==8)init_r8(fft_p,r,bc_w,CP_w,cyclotomic_p,m_th,phi_m); 
	 if(r==16){
		 if(IsMixed == 0)init_r16(fft_p,r,bc_w,CP_w,cyclotomic_p,m_th,phi_m,0);
		 else init_r16_Mixed_radix(fft_p,r,bc_w,CP_w,cyclotomic_p,m_th,phi_m,0);
	 }
	 parameter_o.close();
}

void SPMB::init_Reconfigure(unsigned long fft_p , unsigned long r,int bc_w,unsigned long CP_w,
	 ZZ cyclotomic_p,long m_th,long phi_m){
     double radix_bit;
	 double FFTP_doub;
	 double stage_ceil;
	 double stage_floor;
	 unsigned long Mixed_radix;
	 radix_bit   = log2(r);
	 FFTP_doub   = log2(fft_p);
	 stage_ceil  = ceil(FFTP_doub/radix_bit);
	 stage_floor = floor(FFTP_doub/radix_bit);
	 if(stage_ceil == stage_floor) IsMixed = 0;
	 else IsMixed = 1;
	 std::ofstream parameter_o("./fft_parameter.txt");
	 parameter_o <<"---------------------------------\n";
	 parameter_o <<"fft_point: " << fft_p <<"\n";
	 parameter_o <<"radix: " << r <<"\n";
	 parameter_o <<"butterfly counter width: " << bc_w <<"\n";
	 parameter_o <<"cyclotomic prime width: " << CP_w <<"\n";
	 parameter_o <<"cyclotomic prime : " << cyclotomic_p <<"\n";
	 parameter_o <<"m-th : " << m_th <<"\n";
	 parameter_o <<"phi_m : " << phi_m <<"\n";
	 
	 Mixed_radix = fft_p;
	 while((Mixed_radix % 16) == 0){
		 Mixed_radix = (unsigned long) Mixed_radix / 16; 
	 }
	 
	 std::cout << "--------------------------------------\n";
	 std::cout << "IsMixed: "     << IsMixed <<"\n";
	 std::cout << "Mixed_radix: " << Mixed_radix <<"\n";
	 

     if(IsMixed == 0)init_r16(fft_p,r,bc_w,CP_w,cyclotomic_p,m_th,phi_m,1);
     else init_r16_Mixed_radix(fft_p,r,bc_w,CP_w,cyclotomic_p,m_th,phi_m,1);
	 
	 parameter_o.close();
}

void SPMB::time_o(std::vector<ZZ> time_data,std::string string_in){
   if(radix == 4 )  time_o_r4(time_data,string_in);
   if(radix == 8 )  time_o_r8(time_data,string_in);
   if(radix == 16)  time_o_r16(time_data,string_in);
   if(radix == 0 )  std::cout << "parameter need to init!\n";
}

void SPMB::H_freq_o(std::vector<ZZ> B_NTT){
     unsigned long Mixed_radix;
     
	 Mixed_radix = fft_point;
	 
	 while(Mixed_radix % 16 == 0){
		 Mixed_radix = Mixed_radix / 16; 
	 }

	 //radix-4
     if(radix == 4 )  {
	     if(IsMixed == 0){
	       H_freq_o_r4(B_NTT);
	       H_freq_o_r4_Mux16(B_NTT);	   
	     }else {
	        H_freq_o_r4_r2(B_NTT); 
            H_freq_o_r4_r2_Mux16(B_NTT);	   
	    }
    }
	//radix-8
    if(radix == 8 )H_freq_o_r8(B_NTT);
	//radix-16 and mixed radix generate
    if(radix == 16){ 
	  if(IsMixed == 0){
		 H_freq_o_r16(B_NTT);
	  }else {
		 if(Mixed_radix == 2) H_freq_o_r16_r2(B_NTT);
		 if(Mixed_radix == 4) H_freq_o_r16_r4(B_NTT);
		 if(Mixed_radix == 8) H_freq_o_r16_r8(B_NTT);
	  }
	}
    if(radix == 0 )  std:cout << "parameter need to init!\n"; 
}

void SPMB::re_order_factor(ZZ m_2_rou){
	if(radix == 4 ) re_order_factor_r4(m_2_rou);
	if(radix == 8 ) re_order_factor_r8(m_2_rou);
	if(radix == 16) re_order_factor_r16(m_2_rou);
	if(radix == 0 )  std::cout << "parameter need to init!\n";
}

//===================================================================================================================
//radix-4 data
void SPMB::init_r4(unsigned long fft_p , unsigned long r,int bc_w,unsigned long CP_w,
	 ZZ cyclotomic_p,long m_th,long phi_m){
    fft_point        = fft_p;
    radix            = r;
    bc_width         = bc_w;
	CP_width         = CP_w;
	cyclotomic_prime = cyclotomic_p;
	m                = m_th;
	phim             = phi_m;
    
    offset             =  (int)fft_point/radix;
    group              =  (int)fft_point/(radix*radix);
    fft_point_bit_size =  (int)log2(fft_point); 
    counter_iteration = fft_point/(radix*2);
    ma.resize(fft_point);
    bn.resize(fft_point);
    bit_array_tmp.resize(fft_point_bit_size);
    
    int ma_tmp;
    int bit_tmp;
    int bn_tmp;
    int tmp = 0;
    int addr_tmp = 0;
    
    for(int i=0; i < group; i++){        
        for(int ss=0 ; ss < radix; ss++){
            bn_tmp = 0;
            ma_tmp = 0;
            tmp = ss * group + i;    
            //bit calculate
            for(int j=0; j < fft_point_bit_size;j++){
                bit_tmp = tmp % 2;
                tmp = tmp >> 1;
                bit_array_tmp[j] = bit_tmp;
            } 
            for(int rs = 0; rs < fft_point_bit_size; rs++){
                if((bit_array_tmp[rs] == 1) && (rs != 0)) ma_tmp = ma_tmp + exp2((rs-1)); 
                bn_tmp = bn_tmp ^ bit_array_tmp[rs];
            }
            for(int kk=0; kk<radix; kk++){
                addr_tmp = ss*group + i + kk * (offset);
                if(addr_tmp >= fft_point)std::cout << addr_tmp <<"\n";
                ma[addr_tmp] = ma_tmp;
                bn[addr_tmp] = bn_tmp;
            }
        }    
    }
    
    ofstream index_file_bn0;
    ofstream index_file_bn1;
    index_file_bn0.open("./SPMB/SPMB_bn0.txt");
    index_file_bn1.open("./SPMB/SPMB_bn1.txt");
    for(int i=0;i < fft_point; i++){
       if(bn[i]==0) index_file_bn0 <<" index:" << i << ", MA:" << ma[i] <<"\n";
       else index_file_bn1 <<" index:" << i << ", MA:" << ma[i] << "\n";
    }
    
   
    //calculate inverse SPMB bn and ma 
    ima.resize(fft_point);
    ibn.resize(fft_point);
    ibit_br.resize(bc_width);
    ibit_array_tmp.resize(fft_point_bit_size);
    int ima_tmp;
    int ibit_tmp;
    int ibn_tmp;
    int itmp = 0;
    int iaddr_tmp=0;
    
    for(int i=0; i < group; i++){
        for(int ss = 0;ss < radix; ss++){
            ibn_tmp = 0;
            ima_tmp = 0; 
            itmp = ss * group + i; 
        //bit calculate
            for(int j=0; j < fft_point_bit_size;j++){
                ibit_tmp = itmp % 2;
                itmp = itmp >> 1;
                ibit_array_tmp[j] = ibit_tmp;
            }
            //radix-4 bit-reverse
			for(int br_index=0;br_index < (bc_width/2); br_index++){
				ibit_br[(bc_width-1) - 2*br_index     ]= ibit_array_tmp[2 * br_index + 1]; 
				ibit_br[(bc_width-1) - 2*br_index - 1 ]= ibit_array_tmp[2 * br_index + 0]; 
			}
			
            //bit-reverse over
            for(int rs = 0; rs < bc_width; rs++){
                if((ibit_br[rs] == 1) && (rs != 0)) ima_tmp = ima_tmp + exp2((rs-1)); 
                ibn_tmp = ibn_tmp ^ ibit_br[rs];
            }
            for(int kk=0; kk< radix; kk++){
                iaddr_tmp = (ss*group) + i + (kk * (offset));
                if(iaddr_tmp >= fft_point)std::cout << iaddr_tmp <<"\n";
                ima[iaddr_tmp] = ima_tmp;
                ibn[iaddr_tmp] = ibn_tmp;  
            }
        }   
    }
    
    ofstream index_file_ibn0;
    ofstream index_file_ibn1;
    index_file_ibn0.open("./SPMB/SPMB_ibn0.txt");
    index_file_ibn1.open("./SPMB/SPMB_ibn1.txt");
    for(int i=0;i < fft_point; i++){
        if(ibn[i]==0) index_file_ibn0 <<" index:" << i << ", MA:" << ima[i] <<"\n";
        else index_file_ibn1 <<" index:" << i << ", MA:" << ima[i] << "\n";
    }
    //=========================================================================================================
    //ROM DATA Generate 
	// rom word size : total worde size is divided  into 4 banks  
    std::cout << "*************************************************\n";
	std::cout << "IMPORTANT!!! ROM Data Generate command ===> ";
     
	 unsigned long ROM_word_size;
	 ROM_word_size = fft_p / ( 2  * radix);
	 if(ROM_word_size > 4096 ) std::cout << " make R4_ROM_ALL \n";
	 else std::cout << " make R4_ALL \n";
	 
	 std::cout << "**************************************************\n";
	if(ROM_word_size > 4096 ){
		r4_FFT_TW_ROM_D4();
		r4_IFFT_TW_ROM_D4();
	}else {
		r4_FFT_TW_ROM();
		r4_IFFT_TW_ROM();
	}
}
void SPMB::init_r4_r2(unsigned long fft_p , unsigned long r,int bc_w,unsigned long CP_w,
	 ZZ cyclotomic_p,long m_th,long phi_m){
    fft_point        = fft_p;
    radix            = r;
    bc_width         = bc_w;
	CP_width         = CP_w;
	cyclotomic_prime = cyclotomic_p;
	m                = m_th;
	phim             = phi_m;
    
    offset             =  (int)fft_point/radix;
    group              =  (int)fft_point/(radix*radix);
    fft_point_bit_size =  (int)log2(fft_point); 
    counter_iteration = fft_point/(radix*2);
    ma.resize(fft_point);
    bn.resize(fft_point);
    bit_array_tmp.resize(bc_width);
    
    int ma_tmp;
    int bit_tmp;
    int bn_tmp;
    int BC_tmp;
    int addr_tmp = 0;
    
	BC_tmp = 0;
	
    for(int i=0; i < group; i++){        
        for(int ss=0 ; ss < radix; ss++){
            bn_tmp = 0;
            ma_tmp = 0;
            BC_tmp = ss * group + i;    
            //bit calculate
            for(int j=0; j < bc_width;j++){
                bit_tmp = BC_tmp % 2;
                BC_tmp = BC_tmp >> 1;
                bit_array_tmp[j] = bit_tmp;
            } 
            for(int rs = 0; rs < bc_width; rs++){
                if((bit_array_tmp[rs] == 1) && (rs != 0)) ma_tmp = ma_tmp + exp2((rs-1)); 
                bn_tmp = bn_tmp ^ bit_array_tmp[rs];
            }
            for(int kk=0; kk<radix; kk++){
                addr_tmp = ss*group + i + kk * (offset);
                if(addr_tmp >= fft_point)std::cout << addr_tmp <<"\n";
                ma[addr_tmp] = ma_tmp;
                bn[addr_tmp] = bn_tmp;
            }
        }    
    }
    
    ofstream index_file_bn0;
    ofstream index_file_bn1;
    index_file_bn0.open("./SPMB/SPMB_bn0.txt");
    index_file_bn1.open("./SPMB/SPMB_bn1.txt");
    for(int i=0;i < fft_point; i++){
       if(bn[i]==0) index_file_bn0 <<" index:" << i << ", MA:" << ma[i] <<"\n";
       else index_file_bn1 <<" index:" << i << ", MA:" << ma[i] << "\n";
    }
    
   
    //calculate inverse SPMB bn and ma 
	ofstream R4_R2_INTT_Index_AGU("./SPMB/R4_R2_INTT_Index_AGU.txt");
	
    ima.resize(fft_point);
    ibn.resize(fft_point);
    ibit_br.resize(bc_width);
    ibit_array_tmp.resize(bc_width);
    int ima_tmp;
    int ibit_tmp;
    int ibn_tmp;
    int iaddr_tmp;
	
	for(int i = 0;i < fft_point;i++){
		ibn_tmp = 0;
		ima_tmp = 0;
		R4_R2_INTT_Index_AGU << "--------------------------\n";
		R4_R2_INTT_Index_AGU << "i : " << i << "\n";
		R4_R2_Pointwise_Mult_Index_AGU(i,ima_tmp,ibn_tmp);
		R4_R2_INTT_Index_AGU << "ibn_tmp : " << ibn_tmp << "\n";
		R4_R2_INTT_Index_AGU << "ima_tmp : " << ima_tmp << "\n";
		ima[i] = ima_tmp;
		ibn[i] = ibn_tmp;
	}
      
    ofstream index_file_ibn0;
    ofstream index_file_ibn1;
    index_file_ibn0.open("./SPMB/SPMB_ibn0.txt");
    index_file_ibn1.open("./SPMB/SPMB_ibn1.txt");
    for(int i=0;i < fft_point; i++){
        if(ibn[i]==0) index_file_ibn0 <<" index:" << i << ", MA:" << ima[i] <<"\n";
        else index_file_ibn1 <<" index:" << i << ", MA:" << ima[i] << "\n";
    }
    //=========================================================================================================
    //ROM DATA Generate 
	// rom word size : total worde size is divided  into 4 banks  
    std::cout << "*************************************************\n";
	std::cout << "IMPORTANT!!! ROM Data Generate command ===> ";
     
	 unsigned long ROM_word_size;
	 ROM_word_size = fft_p / ( 2  * radix);
	 if(ROM_word_size > 4096 ) std::cout << " make ROM \n";
	 else std::cout << " make ROM_ND \n";
	 
	 std::cout << "**************************************************\n";
	if(ROM_word_size > 4096 ){
		r4_FFT_TW_ROM_D4();
		r4_IFFT_TW_ROM_D4();
	}else {
		r4_FFT_TW_ROM();
		r4_IFFT_TW_ROM();
	}
}
void SPMB::R4_R2_Pointwise_Mult_Index_AGU(int index,int &MA,int &BN){
    int Isodd;
	int BC_tmp;
	int index_tmp;
	int index_bit_tmp;
	int weight_tmp;
	int bn_tmp;
	int odd_index_tmp;
	int even_index_tmp;
    int index_reorder_lsb;	
    int index_reorder_p1;	
	
	std::vector<int> Index_bit_array;
	std::vector<int> Index_bit_array_tmp;
    Index_bit_array.resize(fft_point_bit_size);
    Index_bit_array_tmp.resize(fft_point_bit_size);
	
	index_tmp = index;
    for(int i=0;i<fft_point_bit_size;i++){
       index_bit_tmp = index_tmp % 2;
       index_tmp = index_tmp >> 1; // right shift 1 bits
	   Index_bit_array[i]      = index_bit_tmp;
	   Index_bit_array_tmp[i]  = index_bit_tmp;
    }
	
     //-----------------
	 Isodd = 1;
	 index_reorder_lsb = fft_point_bit_size - 1;
	 index_reorder_p1  = index_reorder_lsb  - 2;
	 odd_index_tmp  = fft_point_bit_size - 2;
	 even_index_tmp = index_reorder_p1 - 2;

     //data index mapping to memory address 
     //for example 
     for(int j = 0; j < fft_point_bit_size;j++){
        if(j==0) Index_bit_array[j] = Index_bit_array_tmp[index_reorder_lsb];
        else if(j == 1)Index_bit_array[j] = Index_bit_array_tmp[index_reorder_p1];
		else if(j == fft_point_bit_size - 1 )Index_bit_array[j] = Index_bit_array_tmp[1];
		else if(j == fft_point_bit_size - 2 )Index_bit_array[j] = Index_bit_array_tmp[0];
		else {
			if(Isodd == 1){
			  	Index_bit_array[j] = Index_bit_array_tmp[odd_index_tmp];
				odd_index_tmp = odd_index_tmp - 2;
				Isodd = 0;
			}else {
				Index_bit_array[j] = Index_bit_array_tmp[even_index_tmp];
				even_index_tmp = even_index_tmp - 2;
				Isodd = 1;
			}
        }
     }
	 
	
	 BC_tmp = 0;
	 bn_tmp = 0;
	 weight_tmp = 0;
     for(int j = 2; j < fft_point_bit_size;j++){	
        if(Index_bit_array[j] == 1) weight_tmp = 1 << (j - 2);
		else weight_tmp = 0;
		bn_tmp = bn_tmp ^ Index_bit_array[j];
		BC_tmp = BC_tmp + weight_tmp;
	 }
	 
	//output  
	MA     = BC_tmp / 2;
    BN     = bn_tmp;
	 
}
void SPMB::time_o_r4(std::vector<ZZ> time_data,std::string string_in){
    ofstream b0radix0(string_in + "b0radix0.txt");
    ofstream b0radix1(string_in + "b0radix1.txt");
    ofstream b0radix2(string_in + "b0radix2.txt");
    ofstream b0radix3(string_in + "b0radix3.txt");
    ofstream b1radix0(string_in + "b1radix0.txt");
    ofstream b1radix1(string_in + "b1radix1.txt");
    ofstream b1radix2(string_in + "b1radix2.txt");
    ofstream b1radix3(string_in + "b1radix3.txt");
   
    std::string string_buf; 
    
    ZZ time_data_tmp; 

    int counter =0;
    int loop_index = 0;
    int loop_address = 0;
 
    //bank0
    while(counter < counter_iteration){ 
       if((ma[loop_index] == counter )&& (bn[loop_index] == 0) ){
          for(int j=0;j<radix;j++){
            loop_address = loop_index + offset * j;
            time_data_tmp = time_data[loop_address];
            string_buf  = ZZtohex_cp(time_data_tmp);
            
            if(j==0) {
              b0radix0 << string_buf;
              b0radix0 << "\n";
            }     
            if(j==1) { 
              b0radix1 << string_buf;
              b0radix1 << "\n";
            }
            if(j==2) {
              b0radix2 << string_buf;
              b0radix2 << "\n";
            } 
            if(j==3) {
              b0radix3 << string_buf;
              b0radix3 << "\n";           
            }
          }           
           loop_index = 0;
           counter = counter + 1;
       }
       else loop_index = loop_index + 1;
    }
  
    counter      = 0;
    loop_index   = 0;
    loop_address = 0;
    //bank 1
    while(counter < counter_iteration){ 
        if((ma[loop_index] == counter )&& (bn[loop_index] == 1) ){
            for(int j=0;j<radix;j++){
            loop_address = loop_index + offset*j;
            time_data_tmp = time_data[loop_address];
            string_buf  = ZZtohex_cp(time_data_tmp);
            
            if(j==0) {
                b1radix0 << string_buf;
                b1radix0 << "\n";
            }     
            if(j==1) {
                b1radix1 << string_buf;
                b1radix1 << "\n";
            }
            if(j==2) {
                b1radix2 << string_buf;
                b1radix2 << "\n";
            } 
            if(j==3) {
                b1radix3 << string_buf;
                b1radix3 << "\n";           
            }
            }           
            loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
    
    //h2_time_of.close();
    b0radix0.close();
    b0radix1.close();
    b0radix2.close();
    b0radix3.close();
    b1radix0.close();
    b1radix1.close();
    b1radix2.close();
    b1radix3.close();
}
//synthesis mux-8 max word size =4096  , bit size = 128
void SPMB::H_freq_o_r4(std::vector<ZZ> H_NTT){
    //bank 0
    std::ofstream H_freq_out("./ROM_Data/H_freq_out.txt");
    std::ofstream H_b0ROM0("./ROM_Data/H_b0ROM0.txt"); // radix0 ,radix1
    std::ofstream H_b0ROM1("./ROM_Data/H_b0ROM1.txt"); //radix2 radix3
    //bank1                     ROM_Data
    std::ofstream H_b1ROM0("./ROM_Data/H_b1ROM0.txt"); //radix0 radix1
    std::ofstream H_b1ROM1("./ROM_Data/H_b1ROM1.txt"); // radix2 radix3
	
	//bank0
	std::ofstream H_NTT_B0R0("./SPMB/H_NTT_B0R0.txt");
	std::ofstream H_NTT_B0R1("./SPMB/H_NTT_B0R1.txt");
	std::ofstream H_NTT_B0R2("./SPMB/H_NTT_B0R2.txt");
	std::ofstream H_NTT_B0R3("./SPMB/H_NTT_B0R3.txt");
    //bank1 
	std::ofstream H_NTT_B1R0("./SPMB/H_NTT_B1R0.txt");
	std::ofstream H_NTT_B1R1("./SPMB/H_NTT_B1R1.txt");
	std::ofstream H_NTT_B1R2("./SPMB/H_NTT_B1R2.txt");
	std::ofstream H_NTT_B1R3("./SPMB/H_NTT_B1R3.txt");
	//------------------------------------------------
    std::vector<int> tobit_r0;
    std::vector<int> tobit_r1;
    std::vector<int> tobit_r2;
    std::vector<int> tobit_r3;
    tobit_r0.resize(64);
    tobit_r1.resize(64);
    tobit_r2.resize(64);
    tobit_r3.resize(64);
    
    ZZ H_NTT_r0;
    ZZ H_NTT_r1;
    ZZ H_NTT_r2;
    ZZ H_NTT_r3;
    
    int counter =0;
    int loop_index = 0;
    
	std::cout <<"SPMB H_freq_o func star!\n";
    for(int i=0;i<fft_point;i++){
        H_freq_out << H_NTT[i];
        H_freq_out << "\n";
    }
    
    //bank0
    while(counter < counter_iteration){ 
	    if(loop_index > fft_point)std::cout << "loop_index: " << loop_index <<"\n";
	    if(loop_index > fft_point)std::cout << "counter: " << counter <<"\n";
        if((ima[loop_index] == counter )&& (ibn[loop_index] == 0) ){
			
            H_NTT_r0 = H_NTT[loop_index];
            H_NTT_r1 = H_NTT[loop_index + offset];
            H_NTT_r2 = H_NTT[loop_index + 2 * offset];
            H_NTT_r3 = H_NTT[loop_index + 3 * offset];
            
			H_NTT_B0R0 << H_NTT_r0 << "\n";
			H_NTT_B0R1 << H_NTT_r1 << "\n";
			H_NTT_B0R2 << H_NTT_r2 << "\n";
			H_NTT_B0R3 << H_NTT_r3 << "\n";
			
            for(int bit_index=0; bit_index < 64 ;bit_index++){
                if(H_NTT_r0 % 2 == 1) tobit_r0[bit_index] = 1;
                else tobit_r0[bit_index] = 0;
				if(H_NTT_r1 % 2 == 1) tobit_r1[bit_index] = 1;
                else tobit_r1[bit_index] = 0;
				if(H_NTT_r2 % 2 == 1) tobit_r2[bit_index] = 1;
                else tobit_r2[bit_index] = 0;
				if(H_NTT_r3 % 2 == 1) tobit_r3[bit_index] = 1;
                else tobit_r3[bit_index] = 0;
                H_NTT_r0 = H_NTT_r0 >> 1;
                H_NTT_r1 = H_NTT_r1 >> 1;
                H_NTT_r2 = H_NTT_r2 >> 1;
                H_NTT_r3 = H_NTT_r3 >> 1;
            }
            
            for(int g=0;g < 64; g++){
                H_b0ROM0 << tobit_r0[63-g];
            } 

            for(int g=0;g < 64; g++){
                H_b0ROM0 << tobit_r1[63-g];
            }
			
			H_b0ROM0 << "\n";
            
            for(int g=0;g < 64; g++){
                H_b0ROM1 << tobit_r2[63-g];
            }
			
            for(int g=0;g < 64; g++){   
                H_b0ROM1 << tobit_r3[63-g];
            }
            H_b0ROM1 << "\n";
				
            loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
    
	H_NTT_B0R0.close();
	H_NTT_B0R1.close();
	H_NTT_B0R2.close();
	H_NTT_B0R3.close();
	
    counter = 0;
    loop_index = 0;
    std::cout <<"SPMB bn0 H_freq_o func over!\n";
	//bank 1
    while(counter < counter_iteration){ 
        if((ima[loop_index] == counter )&& (ibn[loop_index] == 1) ){
            
			H_NTT_r0 = H_NTT[loop_index];
            H_NTT_r1 = H_NTT[loop_index + offset];
            H_NTT_r2 = H_NTT[loop_index + 2 * offset];
            H_NTT_r3 = H_NTT[loop_index + 3 * offset];
			
			H_NTT_B1R0 << H_NTT_r0 << "\n";
			H_NTT_B1R1 << H_NTT_r1 << "\n";
			H_NTT_B1R2 << H_NTT_r2 << "\n";
			H_NTT_B1R3 << H_NTT_r3 << "\n";
			
            for(int bit_index=0; bit_index < 64 ;bit_index++){
                if(H_NTT_r0 % 2 == 1) tobit_r0[bit_index] = 1;
                else tobit_r0[bit_index] = 0;
				if(H_NTT_r1 % 2 == 1) tobit_r1[bit_index] = 1;
                else tobit_r1[bit_index] = 0;
				if(H_NTT_r2 % 2 == 1) tobit_r2[bit_index] = 1;
                else tobit_r2[bit_index] = 0;
				if(H_NTT_r3 % 2 == 1) tobit_r3[bit_index] = 1;
                else tobit_r3[bit_index] = 0;
                H_NTT_r0 = H_NTT_r0 >> 1;
                H_NTT_r1 = H_NTT_r1 >> 1;
                H_NTT_r2 = H_NTT_r2 >> 1;
                H_NTT_r3 = H_NTT_r3 >> 1;
            }
            
            
            for(int g=0;g < 64; g++){
                H_b1ROM0 << tobit_r0[63-g];
            } 
            
            for(int g=0;g < 64; g++){
                H_b1ROM0 << tobit_r1[63-g];
            }
            
			H_b1ROM0 << "\n"; 
            
            for(int g=0;g < 64; g++){
                H_b1ROM1 << tobit_r2[63-g];
            }
                    
            for(int g=0;g < 64; g++){   
                H_b1ROM1 << tobit_r3[63-g];
            }
    
            H_b1ROM1 << "\n";           
            
			loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
	
	H_NTT_B1R0.close();
	H_NTT_B1R1.close();
	H_NTT_B1R2.close();
	H_NTT_B1R3.close();
}

void SPMB::H_freq_o_r4_r2(std::vector<ZZ> H_NTT){
    //bank 0
    std::ofstream H_freq_out("./ROM_Data/H_freq_out.txt");
    std::ofstream H_b0ROM0("./ROM_Data/H_b0ROM0.txt"); // radix0 ,radix1
    std::ofstream H_b0ROM1("./ROM_Data/H_b0ROM1.txt"); //radix2 radix3
    //bank1                     ROM_Data
    std::ofstream H_b1ROM0("./ROM_Data/H_b1ROM0.txt"); //radix0 radix1
    std::ofstream H_b1ROM1("./ROM_Data/H_b1ROM1.txt"); // radix2 radix3
	
	//bank0
	std::ofstream H_NTT_B0R0("./SPMB/H_NTT_B0R0.txt");
	std::ofstream H_NTT_B0R1("./SPMB/H_NTT_B0R1.txt");
	std::ofstream H_NTT_B0R2("./SPMB/H_NTT_B0R2.txt");
	std::ofstream H_NTT_B0R3("./SPMB/H_NTT_B0R3.txt");
    //bank1 
	std::ofstream H_NTT_B1R0("./SPMB/H_NTT_B1R0.txt");
	std::ofstream H_NTT_B1R1("./SPMB/H_NTT_B1R1.txt");
	std::ofstream H_NTT_B1R2("./SPMB/H_NTT_B1R2.txt");
	std::ofstream H_NTT_B1R3("./SPMB/H_NTT_B1R3.txt");
	//------------------------------------------------
    std::vector<int> tobit_r0;
    std::vector<int> tobit_r1;
    std::vector<int> tobit_r2;
    std::vector<int> tobit_r3;
    tobit_r0.resize(64);
    tobit_r1.resize(64);
    tobit_r2.resize(64);
    tobit_r3.resize(64);
    
    ZZ H_NTT_r0;
    ZZ H_NTT_r1;
    ZZ H_NTT_r2;
    ZZ H_NTT_r3;
    
    int counter =0;
    int loop_index = 0;
	int INTT_R1_offset; //between the R0 index and R1 index
	int INTT_R2_offset; //between the R0 index and R2 index
	int INTT_R3_offset; //between the R0 index and R2 index
    
	
	INTT_R1_offset = (int)fft_point / 2;
	INTT_R2_offset = (int)fft_point / 8;
	INTT_R3_offset = INTT_R1_offset + INTT_R2_offset;
	
	
	std::cout <<"SPMB H_freq_o func star!\n";
    for(int i=0;i<fft_point;i++){
        H_freq_out << H_NTT[i];
        H_freq_out << "\n";
    }
    
    //bank0
    while(counter < counter_iteration){ 
        if((ima[loop_index] == counter )&& (ibn[loop_index] == 0) ){
			
            H_NTT_r0 = H_NTT[loop_index];
            H_NTT_r1 = H_NTT[loop_index + INTT_R1_offset];
            H_NTT_r2 = H_NTT[loop_index + INTT_R2_offset];
            H_NTT_r3 = H_NTT[loop_index + INTT_R3_offset];
            
			H_NTT_B0R0 << H_NTT_r0 << "\n";
			H_NTT_B0R1 << H_NTT_r1 << "\n";
			H_NTT_B0R2 << H_NTT_r2 << "\n";
			H_NTT_B0R3 << H_NTT_r3 << "\n";
			
            for(int bit_index=0; bit_index < 64 ;bit_index++){
                if(H_NTT_r0 % 2 == 1) tobit_r0[bit_index] = 1;
                else tobit_r0[bit_index] = 0;
				if(H_NTT_r1 % 2 == 1) tobit_r1[bit_index] = 1;
                else tobit_r1[bit_index] = 0;
				if(H_NTT_r2 % 2 == 1) tobit_r2[bit_index] = 1;
                else tobit_r2[bit_index] = 0;
				if(H_NTT_r3 % 2 == 1) tobit_r3[bit_index] = 1;
                else tobit_r3[bit_index] = 0;
                H_NTT_r0 = H_NTT_r0 >> 1;
                H_NTT_r1 = H_NTT_r1 >> 1;
                H_NTT_r2 = H_NTT_r2 >> 1;
                H_NTT_r3 = H_NTT_r3 >> 1;
            }
            
            for(int g=0;g < 64; g++){
                H_b0ROM0 << tobit_r0[63-g];
            } 

            for(int g=0;g < 64; g++){
                H_b0ROM0 << tobit_r1[63-g];
            }
			
			H_b0ROM0 << "\n";
            
            for(int g=0;g < 64; g++){
                H_b0ROM1 << tobit_r2[63-g];
            }
			
            for(int g=0;g < 64; g++){   
                H_b0ROM1 << tobit_r3[63-g];
            }
            H_b0ROM1 << "\n";
				
            loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
    
	H_NTT_B0R0.close();
	H_NTT_B0R1.close();
	H_NTT_B0R2.close();
	H_NTT_B0R3.close();
	
    counter = 0;
    loop_index = 0;
    std::cout <<"SPMB bn0 H_freq_o func over!\n";
	//bank 1
    while(counter < counter_iteration){ 
        if((ima[loop_index] == counter )&& (ibn[loop_index] == 1) ){
            
			H_NTT_r0 = H_NTT[loop_index];
            H_NTT_r1 = H_NTT[loop_index + INTT_R1_offset];
            H_NTT_r2 = H_NTT[loop_index + INTT_R2_offset];
            H_NTT_r3 = H_NTT[loop_index + INTT_R3_offset];
			
			H_NTT_B1R0 << H_NTT_r0 << "\n";
			H_NTT_B1R1 << H_NTT_r1 << "\n";
			H_NTT_B1R2 << H_NTT_r2 << "\n";
			H_NTT_B1R3 << H_NTT_r3 << "\n";
			
            for(int bit_index=0; bit_index < 64 ;bit_index++){
                if(H_NTT_r0 % 2 == 1) tobit_r0[bit_index] = 1;
                else tobit_r0[bit_index] = 0;
				if(H_NTT_r1 % 2 == 1) tobit_r1[bit_index] = 1;
                else tobit_r1[bit_index] = 0;
				if(H_NTT_r2 % 2 == 1) tobit_r2[bit_index] = 1;
                else tobit_r2[bit_index] = 0;
				if(H_NTT_r3 % 2 == 1) tobit_r3[bit_index] = 1;
                else tobit_r3[bit_index] = 0;
                H_NTT_r0 = H_NTT_r0 >> 1;
                H_NTT_r1 = H_NTT_r1 >> 1;
                H_NTT_r2 = H_NTT_r2 >> 1;
                H_NTT_r3 = H_NTT_r3 >> 1;
            }
            
            
            for(int g=0;g < 64; g++){
                H_b1ROM0 << tobit_r0[63-g];
            } 
            
            for(int g=0;g < 64; g++){
                H_b1ROM0 << tobit_r1[63-g];
            }
            
			H_b1ROM0 << "\n"; 
            
            for(int g=0;g < 64; g++){
                H_b1ROM1 << tobit_r2[63-g];
            }
                    
            for(int g=0;g < 64; g++){   
                H_b1ROM1 << tobit_r3[63-g];
            }
    
            H_b1ROM1 << "\n";           
            
			loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
	
	H_NTT_B1R0.close();
	H_NTT_B1R1.close();
	H_NTT_B1R2.close();
	H_NTT_B1R3.close();
}

//synthesis mux-16 max word size = 8192 , bit size = 64
void SPMB::H_freq_o_r4_Mux16(std::vector<ZZ> H_NTT){
    //bank 0
    //std::ofstream H_freq_out("./ROM_Data/H_freq_out.txt");
    std::ofstream H_b0r0("./ROM_Data/H_b0r0.txt"); // radix0 
    std::ofstream H_b0r1("./ROM_Data/H_b0r1.txt"); //radix1
    std::ofstream H_b0r2("./ROM_Data/H_b0r2.txt"); //radix2 
    std::ofstream H_b0r3("./ROM_Data/H_b0r3.txt"); //radix3
    //bank1                     ROM_Data
    std::ofstream H_b1r0("./ROM_Data/H_b1r0.txt"); //radix0 
    std::ofstream H_b1r1("./ROM_Data/H_b1r1.txt"); // radix1 
    std::ofstream H_b1r2("./ROM_Data/H_b1r2.txt"); // radix2 
    std::ofstream H_b1r3("./ROM_Data/H_b1r3.txt"); // radix3 
    
    std::vector<int> tobit_r0;
    std::vector<int> tobit_r1;
    std::vector<int> tobit_r2;
    std::vector<int> tobit_r3;
    tobit_r0.resize(64);
    tobit_r1.resize(64);
    tobit_r2.resize(64);
    tobit_r3.resize(64);
    
    ZZ H_NTT_r0;
    ZZ H_NTT_r1;
    ZZ H_NTT_r2;
    ZZ H_NTT_r3;
    
    int counter =0;
    int loop_index = 0;
   
    //bank0
    while(counter < counter_iteration){ 
        if((ima[loop_index] == counter )&& (ibn[loop_index] == 0) ){
			
            H_NTT_r0 = H_NTT[loop_index];
            H_NTT_r1 = H_NTT[loop_index + offset];
            H_NTT_r2 = H_NTT[loop_index + 2 * offset];
            H_NTT_r3 = H_NTT[loop_index + 3 * offset];
            
            for(int bit_index=0; bit_index < 64 ;bit_index++){
                if(H_NTT_r0 % 2 == 1) tobit_r0[bit_index] = 1;
                else tobit_r0[bit_index] = 0;
				if(H_NTT_r1 % 2 == 1) tobit_r1[bit_index] = 1;
                else tobit_r1[bit_index] = 0;
				if(H_NTT_r2 % 2 == 1) tobit_r2[bit_index] = 1;
                else tobit_r2[bit_index] = 0;
				if(H_NTT_r3 % 2 == 1) tobit_r3[bit_index] = 1;
                else tobit_r3[bit_index] = 0;
                H_NTT_r0 = H_NTT_r0 >> 1;
                H_NTT_r1 = H_NTT_r1 >> 1;
                H_NTT_r2 = H_NTT_r2 >> 1;
                H_NTT_r3 = H_NTT_r3 >> 1;
            }
            
            for(int g=0;g < 64; g++){
                H_b0r0 << tobit_r0[63-g];
            } 
			H_b0r0 << "\n";
			
            for(int g=0;g < 64; g++){
                H_b0r1 << tobit_r1[63-g];
            }
			
			H_b0r1 << "\n";
            
            for(int g=0;g < 64; g++){
                H_b0r2 << tobit_r2[63-g];
            }
			H_b0r2 << "\n";
			
            for(int g=0;g < 64; g++){   
                H_b0r3 << tobit_r3[63-g];
            }
            H_b0r3 << "\n";
				
            loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
    
    counter = 0;
    loop_index = 0;
    
	//bank 1
    while(counter < counter_iteration){ 
        if((ima[loop_index] == counter )&& (ibn[loop_index] == 1) ){
            
			H_NTT_r0 = H_NTT[loop_index];
            H_NTT_r1 = H_NTT[loop_index + offset];
            H_NTT_r2 = H_NTT[loop_index + 2 * offset];
            H_NTT_r3 = H_NTT[loop_index + 3 * offset];
			
            for(int bit_index=0; bit_index < 64 ;bit_index++){
                if(H_NTT_r0 % 2 == 1) tobit_r0[bit_index] = 1;
                else tobit_r0[bit_index] = 0;
				if(H_NTT_r1 % 2 == 1) tobit_r1[bit_index] = 1;
                else tobit_r1[bit_index] = 0;
				if(H_NTT_r2 % 2 == 1) tobit_r2[bit_index] = 1;
                else tobit_r2[bit_index] = 0;
				if(H_NTT_r3 % 2 == 1) tobit_r3[bit_index] = 1;
                else tobit_r3[bit_index] = 0;
                H_NTT_r0 = H_NTT_r0 >> 1;
                H_NTT_r1 = H_NTT_r1 >> 1;
                H_NTT_r2 = H_NTT_r2 >> 1;
                H_NTT_r3 = H_NTT_r3 >> 1;
            }
            
            
            for(int g=0;g < 64; g++){
                H_b1r0 << tobit_r0[63-g];
            } 
            H_b1r0 << "\n";
			
            for(int g=0;g < 64; g++){
                H_b1r1 << tobit_r1[63-g];
            }
            
			H_b1r1 << "\n"; 
            
            for(int g=0;g < 64; g++){
                H_b1r2 << tobit_r2[63-g];
            }
            H_b1r2 << "\n"; 
			
            for(int g=0;g < 64; g++){   
                H_b1r3 << tobit_r3[63-g];
            }
    
            H_b1r3 << "\n";           
            
			loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
}

void SPMB::H_freq_o_r4_r2_Mux16(std::vector<ZZ> H_NTT){
    //bank 0
    //std::ofstream H_freq_out("./ROM_Data/H_freq_out.txt");
    std::ofstream H_b0r0("./ROM_Data/H_b0r0.txt"); // radix0 
    std::ofstream H_b0r1("./ROM_Data/H_b0r1.txt"); //radix1
    std::ofstream H_b0r2("./ROM_Data/H_b0r2.txt"); //radix2 
    std::ofstream H_b0r3("./ROM_Data/H_b0r3.txt"); //radix3
    //bank1                     ROM_Data
    std::ofstream H_b1r0("./ROM_Data/H_b1r0.txt"); //radix0 
    std::ofstream H_b1r1("./ROM_Data/H_b1r1.txt"); // radix1 
    std::ofstream H_b1r2("./ROM_Data/H_b1r2.txt"); // radix2 
    std::ofstream H_b1r3("./ROM_Data/H_b1r3.txt"); // radix3 
    
    std::vector<int> tobit_r0;
    std::vector<int> tobit_r1;
    std::vector<int> tobit_r2;
    std::vector<int> tobit_r3;
    tobit_r0.resize(64);
    tobit_r1.resize(64);
    tobit_r2.resize(64);
    tobit_r3.resize(64);
    
    ZZ H_NTT_r0;
    ZZ H_NTT_r1;
    ZZ H_NTT_r2;
    ZZ H_NTT_r3;
    
    int counter =0;
    int loop_index = 0;
   
    int INTT_R1_offset; //between the R0 index and R1 index
	int INTT_R2_offset; //between the R0 index and R2 index
	int INTT_R3_offset; //between the R0 index and R2 index
    
	
	INTT_R1_offset = (int)fft_point / 2;
	INTT_R2_offset = (int)fft_point / 8;
	INTT_R3_offset = INTT_R1_offset + INTT_R2_offset;
   
    //bank0
    while(counter < counter_iteration){ 
        if((ima[loop_index] == counter )&& (ibn[loop_index] == 0) ){
			
            H_NTT_r0 = H_NTT[loop_index];
            H_NTT_r1 = H_NTT[loop_index + INTT_R1_offset];
            H_NTT_r2 = H_NTT[loop_index + INTT_R2_offset];
            H_NTT_r3 = H_NTT[loop_index + INTT_R3_offset];
            
            for(int bit_index=0; bit_index < 64 ;bit_index++){
                if(H_NTT_r0 % 2 == 1) tobit_r0[bit_index] = 1;
                else tobit_r0[bit_index] = 0;
				if(H_NTT_r1 % 2 == 1) tobit_r1[bit_index] = 1;
                else tobit_r1[bit_index] = 0;
				if(H_NTT_r2 % 2 == 1) tobit_r2[bit_index] = 1;
                else tobit_r2[bit_index] = 0;
				if(H_NTT_r3 % 2 == 1) tobit_r3[bit_index] = 1;
                else tobit_r3[bit_index] = 0;
                H_NTT_r0 = H_NTT_r0 >> 1;
                H_NTT_r1 = H_NTT_r1 >> 1;
                H_NTT_r2 = H_NTT_r2 >> 1;
                H_NTT_r3 = H_NTT_r3 >> 1;
            }
            
            for(int g=0;g < 64; g++){
                H_b0r0 << tobit_r0[63-g];
            } 
			H_b0r0 << "\n";
			
            for(int g=0;g < 64; g++){
                H_b0r1 << tobit_r1[63-g];
            }
			
			H_b0r1 << "\n";
            
            for(int g=0;g < 64; g++){
                H_b0r2 << tobit_r2[63-g];
            }
			H_b0r2 << "\n";
			
            for(int g=0;g < 64; g++){   
                H_b0r3 << tobit_r3[63-g];
            }
            H_b0r3 << "\n";
				
            loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
    
    counter = 0;
    loop_index = 0;
    
	//bank 1
    while(counter < counter_iteration){ 
        if((ima[loop_index] == counter )&& (ibn[loop_index] == 1) ){
            
			H_NTT_r0 = H_NTT[loop_index];
            H_NTT_r1 = H_NTT[loop_index + INTT_R1_offset];
            H_NTT_r2 = H_NTT[loop_index + INTT_R2_offset];
            H_NTT_r3 = H_NTT[loop_index + INTT_R3_offset];
			
            for(int bit_index=0; bit_index < 64 ;bit_index++){
                if(H_NTT_r0 % 2 == 1) tobit_r0[bit_index] = 1;
                else tobit_r0[bit_index] = 0;
				if(H_NTT_r1 % 2 == 1) tobit_r1[bit_index] = 1;
                else tobit_r1[bit_index] = 0;
				if(H_NTT_r2 % 2 == 1) tobit_r2[bit_index] = 1;
                else tobit_r2[bit_index] = 0;
				if(H_NTT_r3 % 2 == 1) tobit_r3[bit_index] = 1;
                else tobit_r3[bit_index] = 0;
                H_NTT_r0 = H_NTT_r0 >> 1;
                H_NTT_r1 = H_NTT_r1 >> 1;
                H_NTT_r2 = H_NTT_r2 >> 1;
                H_NTT_r3 = H_NTT_r3 >> 1;
            }
            
            
            for(int g=0;g < 64; g++){
                H_b1r0 << tobit_r0[63-g];
            } 
            H_b1r0 << "\n";
			
            for(int g=0;g < 64; g++){
                H_b1r1 << tobit_r1[63-g];
            }
            
			H_b1r1 << "\n"; 
            
            for(int g=0;g < 64; g++){
                H_b1r2 << tobit_r2[63-g];
            }
            H_b1r2 << "\n"; 
			
            for(int g=0;g < 64; g++){   
                H_b1r3 << tobit_r3[63-g];
            }
    
            H_b1r3 << "\n";           
            
			loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
}
//re-order factor output are store as rom data
void SPMB::re_order_factor_r4(ZZ m_2_rou){

    std::vector<ZZ> re_order_factor_array;
    re_order_factor_array.resize(fft_point);
    std::ofstream  ReROM_o("./ROM_Data/ReORDER.txt");    
	
    for(unsigned long i=0;i<m;i++){
       unsigned long exp;
       ZZ w_2m_tmp_i;
       exp = pow(i,2);
       exp = exp % (2*m);
       PowerMod(w_2m_tmp_i,m_2_rou,exp,cyclotomic_prime);
       re_order_factor_array[i] = w_2m_tmp_i;
	   //re_order_factor_array[i] = 1;
	   ReROM_o << re_order_factor_array[i] << "\n";
    }
    

	//cyclotomic polynomial prime must be 22~24. 
    //ROM 0 bit size = 64  // stroing bank0 radix-0 and  radix-1
	//ROM 1 bit size = 64 
	//re_order_factor_array,while index > m then data is ZERO.
	//because using SPMB Memory addressing,then this reorder factor  
	//=> ROM0[63:42]=bank0_radix0,ROM0[41:20]=bank0_radix1,ROM0[19:0]=20'b0;
    //=> ROM1[63:42]=bank1_radix0,ROM1[41:20]=bank1_radix1,ROM1[19:0]=20'b0;
	//SPMB Memory addressing, for N = 4096  
	//bank0 MA[0]: index 0 , index 1024 ,index 2048,index 3072
	//bank1 MA[0]: index 1 , index 1025 ,index 2049,index 3073
	std::ofstream reorder_ROM0("./ROM_Data/reorder_ROM0.txt");
    std::ofstream reorder_ROM1("./ROM_Data/reorder_ROM1.txt");
   
    std::vector<int> tobit_r0;
    std::vector<int> tobit_r1;

    tobit_r0.resize(CP_width);
    tobit_r1.resize(CP_width);
    
    ZZ re_order_r0_tmp;
    ZZ re_order_r1_tmp;
    
    int counter =0;
    int loop_index = 0;
    int loop_address = 0;

    //bank0
    while(counter < counter_iteration){ 
        if((ma[loop_index] == counter )&& (bn[loop_index] == 0) ){
            re_order_r0_tmp = re_order_factor_array[loop_index];
            re_order_r1_tmp = re_order_factor_array[loop_index + offset];
            
			//calculate bit
            for(int bit_index=0; bit_index < CP_width ;bit_index++){
                //radix0
				if(re_order_r0_tmp % 2 == 1) tobit_r0[bit_index] = 1;
                else tobit_r0[bit_index] = 0;
				//radix1
				if(re_order_r1_tmp % 2 == 1) tobit_r1[bit_index] = 1;
                else tobit_r1[bit_index] = 0;
				
				re_order_r0_tmp = re_order_r0_tmp >> 1;
				re_order_r1_tmp = re_order_r1_tmp >> 1;
            }
            
           
            for(int g=0;g < CP_width; g++){
				reorder_ROM0 << tobit_r0[CP_width-1-g];
            } 
            
            for(int g=0;g < CP_width; g++){
                reorder_ROM0 << tobit_r1[CP_width-1-g];
            } 
            
			int padding_zero_rom0 ;
			
			padding_zero_rom0 = 64 - 2*CP_width;
            
            for(int g=0;g < padding_zero_rom0; g++){
                reorder_ROM0 << 0;
			}
            reorder_ROM0 << "\n";
                                
            loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
    reorder_ROM0.close();
 
    counter = 0;
    loop_index = 0;
    //bank 1
    while(counter < counter_iteration){ 
        if((ma[loop_index] == counter )&& (bn[loop_index] == 1) ){
            re_order_r0_tmp = re_order_factor_array[loop_index];
            re_order_r1_tmp = re_order_factor_array[loop_index + offset];
            
			//calculate bit
            for(int bit_index=0; bit_index < CP_width ;bit_index++){
                //radix0
				if(re_order_r0_tmp % 2 == 1) tobit_r0[bit_index] = 1;
                else tobit_r0[bit_index] = 0;
				//radix1
				if(re_order_r1_tmp % 2 == 1) tobit_r1[bit_index] = 1;
                else tobit_r1[bit_index] = 0;
				
				re_order_r0_tmp = re_order_r0_tmp >> 1;
				re_order_r1_tmp = re_order_r1_tmp >> 1;
            }
                    
            for(int g=0;g < CP_width; g++){
                reorder_ROM1 << tobit_r0[CP_width-1-g];
            } 
             
            for(int g=0;g < CP_width; g++){
                reorder_ROM1 << tobit_r1[CP_width-1-g];
            } 
            
			int padding_zero_rom1 = 64;
			
			padding_zero_rom1 = padding_zero_rom1 - 2*CP_width;
            
            for(int g=0;g < padding_zero_rom1; g++){
                reorder_ROM1 << 0;
			}
            reorder_ROM1 << "\n";
                 
            loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
    reorder_ROM1.close();

    loop_address = 0;
    
	//inverse order factor 
	std::vector<int> itobit_r0;
	std::vector<int> itobit_r1;
	std::vector<int> itobit_r2;
	std::vector<int> itobit_r3;
	itobit_r0.resize(CP_width);
	itobit_r1.resize(CP_width);
	itobit_r2.resize(CP_width);
	itobit_r3.resize(CP_width);
	ZZ re_iorder_r0_tmp;
	ZZ re_iorder_r1_tmp;
	ZZ re_iorder_r2_tmp;
	ZZ re_iorder_r3_tmp;
	/*************************************************/
	//m-th cyclotomic polynomial 
	//all polynomial at this thing,the coefficient is zero while the degree > phi(m);
    //and inverse reorder factor index > phi(m) is zero.
	//using SPMB (Final stage type) ,the data index is continuous
	//we can divide inverse re-order factor into two banks ROM.
	//ROM0 bit width= 64 ,ROM1 bit width = 32
	//cyclotomic prime = 22 bits
	//ROM0[63:42] = radix0 ,ROM0[41:20] = radix1,ROM0[19:0] = radix2[21:2]
	//ROM1[31:30] = raidx2[1:0],ROM1[29:8] = radix3,ROM1[7:0] = 8'b0;
	//ROM Word size ,original word size is N/radix => 4096/4 = 1024,
	//a_0x^(0)+........+a_phi(m)-1*x^(phi(m)-1)+0*x^phi(m)+...+0^x(m-1)
	//so word size can redution to phi(m)/radix 
    //for example m=1705 phi(m) = 1200 , then ma > 300 data is zero	
	//then we set word size is 512
	/*************************************************/
    std::ofstream ireorder_ROM0("./ROM_Data/ireorder_ROM0.txt");
    std::ofstream ireorder_ROM1("./ROM_Data/ireorder_ROM1.txt");
    
    double IReROM_WORD_BIT;
	int    ROM_WORD_SIZE;
    IReROM_WORD_BIT = ceil((double)m / 4);
    IReROM_WORD_BIT = log2(IReROM_WORD_BIT);
    IReROM_WORD_BIT = ceil(IReROM_WORD_BIT);	
    ROM_WORD_SIZE   = exp2(IReROM_WORD_BIT);
	    
	for(int i=0;i < ROM_WORD_SIZE; i++){
         re_iorder_r0_tmp = re_order_factor_array[radix * i + 0];
         re_iorder_r1_tmp = re_order_factor_array[radix * i + 1];
         re_iorder_r2_tmp = re_order_factor_array[radix * i + 2];
         re_iorder_r3_tmp = re_order_factor_array[radix * i + 3];
        
        for(int bit_index=0; bit_index < CP_width ;bit_index++){
            //radix0
			if(re_iorder_r0_tmp % 2 == 1) itobit_r0[bit_index] = 1;
            else itobit_r0[bit_index] = 0;
			//radix1
			if(re_iorder_r1_tmp % 2 == 1) itobit_r1[bit_index] = 1;
			else itobit_r1[bit_index] = 0;
			//radix2
			if(re_iorder_r2_tmp % 2 == 1) itobit_r2[bit_index] = 1;
			else itobit_r2[bit_index] = 0;
			//radix3
            if(re_iorder_r3_tmp % 2 == 1) itobit_r3[bit_index] = 1;
			else itobit_r3[bit_index] = 0;
			
			re_iorder_r0_tmp = re_iorder_r0_tmp >> 1;
            re_iorder_r1_tmp = re_iorder_r1_tmp >> 1;
            re_iorder_r2_tmp = re_iorder_r2_tmp >> 1;
            re_iorder_r3_tmp = re_iorder_r3_tmp >> 1;
        }
     
        for(int g=0;g < CP_width; g++){
            ireorder_ROM0 << itobit_r0[CP_width-1-g];
        } 
         
        for(int g=0;g < CP_width; g++){
            ireorder_ROM0 << itobit_r1[CP_width-1-g];
        } 
        
		int IROM0_Remaining_bits;
		IROM0_Remaining_bits = 64 - 2 * CP_width;
		
        for(int g=0;g < CP_width; g++){
            if(IROM0_Remaining_bits > g)ireorder_ROM0 << itobit_r2[CP_width-1-g];
			else ireorder_ROM1 << itobit_r2[CP_width-1-g];
        } 
        ireorder_ROM0 << "\n"; 
        
        for(int g=0;g < CP_width; g++){
            ireorder_ROM1 << itobit_r3[CP_width-1-g];
        } 
		
		int IROM1_Padding_zero;
		if(4 * CP_width > 96) {
		    IROM1_Padding_zero = 64 - CP_width -(CP_width -IROM0_Remaining_bits);
		}
		else {
			IROM1_Padding_zero = 32 - CP_width -(CP_width -IROM0_Remaining_bits);
		}
		
		for(int g=0;g < IROM1_Padding_zero; g++){
            ireorder_ROM1 << 0;
        } 
		
        ireorder_ROM1 << "\n"; 
         
    }
    ireorder_ROM0.close();
    ireorder_ROM1.close();

}
//twiddle factor rom is divided into 4  banks
void SPMB::r4_FFT_TW_ROM_D4(){
    //Bank0
    std::ofstream  ROMB0R0("./ROM_Data/R4_FFTROM0_D0_4Bank.txt");
    std::ofstream  ROMB0R1("./ROM_Data/R4_FFTROM0_D1_4Bank.txt");
    std::ofstream  ROMB0R2("./ROM_Data/R4_FFTROM0_D2_4Bank.txt");
    std::ofstream  ROMB0R3("./ROM_Data/R4_FFTROM0_D3_4Bank.txt");
    //Bank1
    std::ofstream  ROMB1R0("./ROM_Data/R4_FFTROM1_D0_4Bank.txt");
    std::ofstream  ROMB1R1("./ROM_Data/R4_FFTROM1_D1_4Bank.txt");
    std::ofstream  ROMB1R2("./ROM_Data/R4_FFTROM1_D2_4Bank.txt");
    std::ofstream  ROMB1R3("./ROM_Data/R4_FFTROM1_D3_4Bank.txt");
    //FFT Prime
    ZZ FFT_Prime;
    ZZ twiddle_65536;
    ZZ R1_tmp;
    ZZ R2_tmp;
    ZZ R3_tmp;
    long exp;
    long order;  // order = 65536 / fft_point
    long addr_length; // FFT_point / radix
    conv(FFT_Prime,"18446744069414584321");     //FFT_Prime = 2^64 - 2^32 +1
    conv(twiddle_65536,"14603442835287214144"); //65536-th root of unity 
    order = (long) 65536 / fft_point;
    addr_length  = (long) fft_point / radix;

    std::vector<int> tobit_r1;
    std::vector<int> tobit_r2;
    std::vector<int> tobit_r3;
    tobit_r1.resize(64);
    tobit_r2.resize(64);
    tobit_r3.resize(64);
    
    std::vector<ZZ> tw_Table;  //twiddle factor table 
    tw_Table.resize(fft_point);
    
    exp = 0;
    for(int i=0;i<fft_point;i++){
       tw_Table[i] = PowerMod(twiddle_65536,exp,FFT_Prime);
       exp = exp + order;
    }

    //radix 0
    for(int i=0;i<(addr_length/4);i++){
        R1_tmp = tw_Table[i];
        R2_tmp = tw_Table[2*i];
        R3_tmp = tw_Table[3*i];
        for(int bit_index=0; bit_index < 64 ;bit_index++){
            if(R1_tmp % 2 == 1) tobit_r1[bit_index] = 1;
            else tobit_r1[bit_index] = 0;
            
            if(R2_tmp % 2 == 1) tobit_r2[bit_index] = 1;
            else tobit_r2[bit_index] = 0;             

            if(R3_tmp % 2 == 1) tobit_r3[bit_index] = 1;
            else tobit_r3[bit_index] = 0;             
            
            R1_tmp = R1_tmp >> 1;
            R2_tmp = R2_tmp >> 1;
            R3_tmp = R3_tmp >> 1;
        }
        
        for(int g=0;g < 64; g++){
            ROMB0R0 << tobit_r1[63-g];
            ROMB1R0 << tobit_r2[63-g];
        }
        for(int g=0;g < 64; g++){
            ROMB1R0 << tobit_r3[63-g];
        }
        ROMB0R0 << "\n";
        ROMB1R0 << "\n";
    }
    ROMB0R0.close();
    ROMB1R0.close();
    //radix1
    for(int i=(addr_length/4);i<(addr_length/2);i++){
        R1_tmp = tw_Table[i];
        R2_tmp = tw_Table[2*i];
        R3_tmp = tw_Table[3*i];
        for(int bit_index=0; bit_index < 64 ;bit_index++){
            if(R1_tmp % 2 == 1) tobit_r1[bit_index] = 1;
            else tobit_r1[bit_index] = 0;
            
            if(R2_tmp % 2 == 1) tobit_r2[bit_index] = 1;
            else tobit_r2[bit_index] = 0;             

            if(R3_tmp % 2 == 1) tobit_r3[bit_index] = 1;
            else tobit_r3[bit_index] = 0;             
            
            R1_tmp = R1_tmp >> 1;
            R2_tmp = R2_tmp >> 1;
            R3_tmp = R3_tmp >> 1;
        }
        
        for(int g=0;g < 64; g++){
            ROMB0R1 << tobit_r1[63-g];
            ROMB1R1 << tobit_r2[63-g];
        }
        for(int g=0;g < 64; g++){
            ROMB1R1 << tobit_r3[63-g];
        }
        ROMB0R1 << "\n";
        ROMB1R1 << "\n";
    }
    ROMB0R1.close();
    ROMB1R1.close();
    //radix2
    for(int i=(addr_length/2);i<((3*addr_length)/4);i++){
        R1_tmp = tw_Table[i];
        R2_tmp = tw_Table[2*i];
        R3_tmp = tw_Table[3*i];
        for(int bit_index=0; bit_index < 64 ;bit_index++){
            if(R1_tmp % 2 == 1) tobit_r1[bit_index] = 1;
            else tobit_r1[bit_index] = 0;
            
            if(R2_tmp % 2 == 1) tobit_r2[bit_index] = 1;
            else tobit_r2[bit_index] = 0;             

            if(R3_tmp % 2 == 1) tobit_r3[bit_index] = 1;
            else tobit_r3[bit_index] = 0;             
            
            R1_tmp = R1_tmp >> 1;
            R2_tmp = R2_tmp >> 1;
            R3_tmp = R3_tmp >> 1;
        }
        
        for(int g=0;g < 64; g++){
            ROMB0R2 << tobit_r1[63-g];
            ROMB1R2 << tobit_r2[63-g];
        }
        for(int g=0;g < 64; g++){
            ROMB1R2 << tobit_r3[63-g];
        }
        ROMB0R2 << "\n";
        ROMB1R2 << "\n";
    }
    ROMB0R2.close();
    ROMB1R2.close();
    //radix3
    for(int i=((3*addr_length)/4);i<(addr_length);i++){
        R1_tmp = tw_Table[i];
        R2_tmp = tw_Table[2*i];
        R3_tmp = tw_Table[3*i];
        for(int bit_index=0; bit_index < 64 ;bit_index++){
            if(R1_tmp % 2 == 1) tobit_r1[bit_index] = 1;
            else tobit_r1[bit_index] = 0;
            
            if(R2_tmp % 2 == 1) tobit_r2[bit_index] = 1;
            else tobit_r2[bit_index] = 0;             

            if(R3_tmp % 2 == 1) tobit_r3[bit_index] = 1;
            else tobit_r3[bit_index] = 0;             
            
            R1_tmp = R1_tmp >> 1;
            R2_tmp = R2_tmp >> 1;
            R3_tmp = R3_tmp >> 1;
        }
        
        for(int g=0;g < 64; g++){
            ROMB0R3 << tobit_r1[63-g];
            ROMB1R3 << tobit_r2[63-g];
        }
        for(int g=0;g < 64; g++){
            ROMB1R3 << tobit_r3[63-g];
        }
        ROMB0R3 << "\n";
        ROMB1R3 << "\n";
    }    
    ROMB0R3.close();
    ROMB1R3.close();
}
void SPMB::r4_FFT_TW_ROM(){
    //---------siang print data -----------
    std::ofstream siang_ROM0("./my_print_data/siang_R4_FFTROM0.txt");
    std::ofstream siang_ROM1("./my_print_data/siang_R4_FFTROM1.txt");
    //----------
	//ROM0 bitsize = 64 , ROM1 bitsize = 128
	//word size = FFT_point / radix
	//Bank0
    std::ofstream  ROM0("./ROM_Data/R4_FFTROM0.txt");
    //Bank1
    std::ofstream  ROM1("./ROM_Data/R4_FFTROM1.txt");
    //FFT Prime
    ZZ FFT_Prime;
    ZZ twiddle_65536;
    ZZ R1_tmp;
    ZZ R2_tmp;
    ZZ R3_tmp;
    long exp;
    long order;  // order = 65536 / fft_point
    long addr_length; // FFT_point / radix
    conv(FFT_Prime,"18446744069414584321");     //FFT_Prime = 2^64 - 2^32 +1
    conv(twiddle_65536,"14603442835287214144"); //65536-th root of unity 
    order = (long) 65536 / fft_point;
    addr_length  = (long) fft_point / radix;

    std::vector<int> tobit_r1;
    std::vector<int> tobit_r2;
    std::vector<int> tobit_r3;
    tobit_r1.resize(64);
    tobit_r2.resize(64);
    tobit_r3.resize(64);
    
    std::vector<ZZ> tw_Table;  //twiddle factor table 
    tw_Table.resize(fft_point);
    
    exp = 0;
    for(int i=0;i<fft_point;i++){
       tw_Table[i] = PowerMod(twiddle_65536,exp,FFT_Prime);
       exp = exp + order;
    }
    
    //radix 0
    for(int i=0;i<(addr_length);i++){
        R1_tmp = tw_Table[i];
        R2_tmp = tw_Table[2*i];
        R3_tmp = tw_Table[3*i];

        //siang twiddle print out
        //i=128
        siang_ROM0 << "exp = "  << i*order    << ", " << R1_tmp;
        siang_ROM1 << "exp = "  << 2*i*order  << ", " << R2_tmp;

        siang_ROM1 << " ,   ";

        siang_ROM1 << "exp = "  << 3*i*order << ", " << R3_tmp;
        
        siang_ROM0 << "\n";
        siang_ROM1 << "\n";
        //------------------------  

        for(int bit_index=0; bit_index < 64 ;bit_index++){
            if(R1_tmp % 2 == 1) tobit_r1[bit_index] = 1;
            else tobit_r1[bit_index] = 0;
            
            if(R2_tmp % 2 == 1) tobit_r2[bit_index] = 1;
            else tobit_r2[bit_index] = 0;             

            if(R3_tmp % 2 == 1) tobit_r3[bit_index] = 1;
            else tobit_r3[bit_index] = 0;             
            
            R1_tmp = R1_tmp >> 1;
            R2_tmp = R2_tmp >> 1;
            R3_tmp = R3_tmp >> 1;
        }
        
        for(int g=0;g < 64; g++){
            ROM0 << tobit_r1[63-g];
            ROM1 << tobit_r2[63-g];
        }
        for(int g=0;g < 64; g++){
            ROM1 << tobit_r3[63-g];
        }
        ROM0 << "\n";
        ROM1 << "\n";
    }
    ROM0.close();
    ROM1.close();
    //---------------
    siang_ROM0.close();
	siang_ROM1.close();
    //----------------
}
//twiddle factor rom is divided into 4  banks
void SPMB::r4_IFFT_TW_ROM_D4(){
    //Bank0
    std::ofstream  ROMB0R0("./ROM_Data/R4_IFFTROM0_D0_4Bank.txt");
    std::ofstream  ROMB0R1("./ROM_Data/R4_IFFTROM0_D1_4Bank.txt");
    std::ofstream  ROMB0R2("./ROM_Data/R4_IFFTROM0_D2_4Bank.txt");
    std::ofstream  ROMB0R3("./ROM_Data/R4_IFFTROM0_D3_4Bank.txt");
    //Bank1                           
    std::ofstream  ROMB1R0("./ROM_Data/R4_IFFTROM1_D0_4Bank.txt");
    std::ofstream  ROMB1R1("./ROM_Data/R4_IFFTROM1_D1_4Bank.txt");
    std::ofstream  ROMB1R2("./ROM_Data/R4_IFFTROM1_D2_4Bank.txt");
    std::ofstream  ROMB1R3("./ROM_Data/R4_IFFTROM1_D3_4Bank.txt");
    //FFT Prime
    ZZ FFT_Prime;
    ZZ twiddle_65536;
    ZZ R1_tmp;
    ZZ R2_tmp;
    ZZ R3_tmp;
    long exp;
    long order;  // order = 65536 / fft_point
    long addr_length; // FFT_point / radix
    conv(FFT_Prime,"18446744069414584321");     //FFT_Prime = 2^64 - 2^32 +1
    conv(twiddle_65536,"14603442835287214144"); //65536-th root of unity 
    order = (long) 65536 / fft_point;
    addr_length  = (long) fft_point / radix;

    std::vector<int> tobit_r1;
    std::vector<int> tobit_r2;
    std::vector<int> tobit_r3;
    tobit_r1.resize(64);
    tobit_r2.resize(64);
    tobit_r3.resize(64);
    
    std::vector<ZZ> tw_Table;  //twiddle factor table 
    tw_Table.resize(fft_point);
    
    exp = 0;
    for(int i=0;i<fft_point;i++){
       tw_Table[i] = PowerMod(twiddle_65536,exp,FFT_Prime);
       exp = exp + order;
    }
    
    //radix 0
    for(int i=0;i<(addr_length/4);i++){
        if(i==0){    
            R1_tmp = tw_Table[0];
            R2_tmp = tw_Table[0];
            R3_tmp = tw_Table[0];
        }
        else {
            R1_tmp = tw_Table[fft_point - i];
            R2_tmp = tw_Table[fft_point - 2*i];
            R3_tmp = tw_Table[fft_point - 3*i];
        }
        for(int bit_index=0; bit_index < 64 ;bit_index++){
            if(R1_tmp % 2 == 1) tobit_r1[bit_index] = 1;
            else tobit_r1[bit_index] = 0;
            
            if(R2_tmp % 2 == 1) tobit_r2[bit_index] = 1;
            else tobit_r2[bit_index] = 0;             

            if(R3_tmp % 2 == 1) tobit_r3[bit_index] = 1;
            else tobit_r3[bit_index] = 0;             
            
            R1_tmp = R1_tmp >> 1;
            R2_tmp = R2_tmp >> 1;
            R3_tmp = R3_tmp >> 1;
        }
        
        for(int g=0;g < 64; g++){
            ROMB0R0 << tobit_r1[63-g];
            ROMB1R0 << tobit_r2[63-g];
        }
        for(int g=0;g < 64; g++){
            ROMB1R0 << tobit_r3[63-g];
        }
        ROMB0R0 << "\n";
        ROMB1R0 << "\n";
    }
    ROMB0R0.close();
    ROMB1R0.close();
    //radix1
    for(int i=(addr_length/4);i<(addr_length/2);i++){
        if(i==0){
            R1_tmp = tw_Table[0];
            R2_tmp = tw_Table[0];
            R3_tmp = tw_Table[0];
        }
        else {
            R1_tmp = tw_Table[fft_point - i];
            R2_tmp = tw_Table[fft_point - 2*i];
            R3_tmp = tw_Table[fft_point - 3*i];
        }
        for(int bit_index=0; bit_index < 64 ;bit_index++){
            if(R1_tmp % 2 == 1) tobit_r1[bit_index] = 1;
            else tobit_r1[bit_index] = 0;
            
            if(R2_tmp % 2 == 1) tobit_r2[bit_index] = 1;
            else tobit_r2[bit_index] = 0;             

            if(R3_tmp % 2 == 1) tobit_r3[bit_index] = 1;
            else tobit_r3[bit_index] = 0;             
            
            R1_tmp = R1_tmp >> 1;
            R2_tmp = R2_tmp >> 1;
            R3_tmp = R3_tmp >> 1;
        }
        
        for(int g=0;g < 64; g++){
            ROMB0R1 << tobit_r1[63-g];
            ROMB1R1 << tobit_r2[63-g];
        }
        for(int g=0;g < 64; g++){
            ROMB1R1 << tobit_r3[63-g];
        }
        ROMB0R1 << "\n";
        ROMB1R1 << "\n";
    }
    ROMB0R1.close();
    ROMB1R1.close();
    //radix2
    for(int i=(addr_length/2);i<((3*addr_length)/4);i++){
        if(i==0){
            R1_tmp = tw_Table[0];
            R2_tmp = tw_Table[0];
            R3_tmp = tw_Table[0]; 
        }
        else {
            R1_tmp = tw_Table[fft_point - i];
            R2_tmp = tw_Table[fft_point - 2*i];
            R3_tmp = tw_Table[fft_point - 3*i];   
        }
        for(int bit_index=0; bit_index < 64 ;bit_index++){
            if(R1_tmp % 2 == 1) tobit_r1[bit_index] = 1;
            else tobit_r1[bit_index] = 0;
            
            if(R2_tmp % 2 == 1) tobit_r2[bit_index] = 1;
            else tobit_r2[bit_index] = 0;             

            if(R3_tmp % 2 == 1) tobit_r3[bit_index] = 1;
            else tobit_r3[bit_index] = 0;             
            
            R1_tmp = R1_tmp >> 1;
            R2_tmp = R2_tmp >> 1;
            R3_tmp = R3_tmp >> 1;
        }
        
        for(int g=0;g < 64; g++){
            ROMB0R2 << tobit_r1[63-g];
            ROMB1R2 << tobit_r2[63-g];
        }
        for(int g=0;g < 64; g++){
            ROMB1R2 << tobit_r3[63-g];
        }
        ROMB0R2 << "\n";
        ROMB1R2 << "\n";
    }
    ROMB0R2.close();
    ROMB1R2.close();
    //radix3
    for(int i=((3*addr_length)/4);i<(addr_length);i++){
        if(i==0){
            R1_tmp = tw_Table[0];
            R2_tmp = tw_Table[0];
            R3_tmp = tw_Table[0];
        }
        else {
            R1_tmp = tw_Table[fft_point - i];
            R2_tmp = tw_Table[fft_point - 2*i];
            R3_tmp = tw_Table[fft_point - 3*i];
        }
        for(int bit_index=0; bit_index < 64 ;bit_index++){
            if(R1_tmp % 2 == 1) tobit_r1[bit_index] = 1;
            else tobit_r1[bit_index] = 0;
            
            if(R2_tmp % 2 == 1) tobit_r2[bit_index] = 1;
            else tobit_r2[bit_index] = 0;             

            if(R3_tmp % 2 == 1) tobit_r3[bit_index] = 1;
            else tobit_r3[bit_index] = 0;             
            
            R1_tmp = R1_tmp >> 1;
            R2_tmp = R2_tmp >> 1;
            R3_tmp = R3_tmp >> 1;
        }
        
        for(int g=0;g < 64; g++){
            ROMB0R3 << tobit_r1[63-g];
            ROMB1R3 << tobit_r2[63-g];
        }
        for(int g=0;g < 64; g++){
            ROMB1R3 << tobit_r3[63-g];
        }
        ROMB0R3 << "\n";
        ROMB1R3 << "\n";
    }    
    ROMB0R3.close();
    ROMB1R3.close();    
}
void SPMB::r4_IFFT_TW_ROM(){
    
	//Bank0
    std::ofstream  ROM0("./ROM_Data/R4_IFFTROM0.txt");
    //Bank1                           
    std::ofstream  ROM1("./ROM_Data/R4_IFFTROM1.txt");
    //FFT Prime
    ZZ FFT_Prime;
    ZZ twiddle_65536;
    ZZ R1_tmp;
    ZZ R2_tmp;
    ZZ R3_tmp;
    long exp;
    long order;  // order = 65536 / fft_point
    long addr_length; // FFT_point / radix
    conv(FFT_Prime,"18446744069414584321");     //FFT_Prime = 2^64 - 2^32 +1
    conv(twiddle_65536,"14603442835287214144"); //65536-th root of unity 
    order = (long) 65536 / fft_point;
    addr_length  = (long) fft_point / radix;

    std::vector<int> tobit_r1;
    std::vector<int> tobit_r2;
    std::vector<int> tobit_r3;
    tobit_r1.resize(64);
    tobit_r2.resize(64);
    tobit_r3.resize(64);
    
    std::vector<ZZ> tw_Table;  //twiddle factor table 
    tw_Table.resize(fft_point);
    
    exp = 0;
    for(int i=0;i<fft_point;i++){
       tw_Table[i] = PowerMod(twiddle_65536,exp,FFT_Prime);
       exp = exp + order;
    }
    
    for(int i=0;i<(addr_length);i++){
        if(i==0){    
            R1_tmp = tw_Table[0];
            R2_tmp = tw_Table[0];
            R3_tmp = tw_Table[0];
        }
        else {
            R1_tmp = tw_Table[fft_point - i];
            R2_tmp = tw_Table[fft_point - 2*i];
            R3_tmp = tw_Table[fft_point - 3*i];
        }
        for(int bit_index=0; bit_index < 64 ;bit_index++){
            if(R1_tmp % 2 == 1) tobit_r1[bit_index] = 1;
            else tobit_r1[bit_index] = 0;
            
            if(R2_tmp % 2 == 1) tobit_r2[bit_index] = 1;
            else tobit_r2[bit_index] = 0;             

            if(R3_tmp % 2 == 1) tobit_r3[bit_index] = 1;
            else tobit_r3[bit_index] = 0;             
            
            R1_tmp = R1_tmp >> 1;
            R2_tmp = R2_tmp >> 1;
            R3_tmp = R3_tmp >> 1;
        }
        
        for(int g=0;g < 64; g++){
            ROM0 << tobit_r1[63-g];
            ROM1 << tobit_r2[63-g];
        }
        for(int g=0;g < 64; g++){
            ROM1 << tobit_r3[63-g];
        }
        ROM0 << "\n";
        ROM1 << "\n";
    }
    ROM0.close();
    ROM1.close();
}
//===================================================================================================================
//radix-8
void SPMB::init_r8(unsigned long fft_p , unsigned long r,int bc_w,unsigned long CP_w,
	 ZZ cyclotomic_p,long m_th,long phi_m){
    fft_point        = fft_p;
    radix            = r;
    bc_width         = bc_w;
	CP_width         = CP_w;
	cyclotomic_prime = cyclotomic_p;
	m                = m_th;
	phim             = phi_m;
    
    offset             =  (int)fft_point/radix;
    group              =  (int)fft_point/(radix*radix);
    fft_point_bit_size =  (int)log2(fft_point); 
    counter_iteration = fft_point/(radix*2);
    ma.resize(fft_point);
    bn.resize(fft_point);
    bit_array_tmp.resize(fft_point_bit_size);
    
    int ma_tmp;
    int bit_tmp;
    int bn_tmp;
    int tmp = 0;
    int addr_tmp = 0;
    
    for(int i=0; i < group; i++){        
        for(int ss=0 ; ss < radix; ss++){
            bn_tmp = 0;
            ma_tmp = 0;
            tmp = ss * group + i;    
            //bit calculate
            for(int j=0; j < fft_point_bit_size;j++){
                bit_tmp = tmp % 2;
                tmp = tmp >> 1;
                bit_array_tmp[j] = bit_tmp;
            } 
            for(int rs = 0; rs < fft_point_bit_size; rs++){
                if((bit_array_tmp[rs] == 1) && (rs != 0)) ma_tmp = ma_tmp + exp2((rs-1)); 
                bn_tmp = bn_tmp ^ bit_array_tmp[rs];
            }
            for(int kk=0; kk<radix; kk++){
                addr_tmp = ss*group + i + kk * (offset);
                if(addr_tmp >= fft_point)std::cout << addr_tmp <<"\n";
                ma[addr_tmp] = ma_tmp;
                bn[addr_tmp] = bn_tmp;
            }
        }    
    }
    
    ofstream index_file_bn0;
    ofstream index_file_bn1;
    index_file_bn0.open("./SPMB/SPMB_bn0.txt");
    index_file_bn1.open("./SPMB/SPMB_bn1.txt");
    for(int i=0;i < fft_point; i++){
       if(bn[i]==0) index_file_bn0 <<" index:" << i << ", MA:" << ma[i] <<"\n";
       else index_file_bn1 <<" index:" << i << ", MA:" << ma[i] << "\n";
    }
    
   
    //calculate inverse SPMB bn and ma 
    ima.resize(fft_point);
    ibn.resize(fft_point);
    ibit_br.resize(bc_width);
    ibit_array_tmp.resize(fft_point_bit_size);
    int ima_tmp;
    int ibit_tmp;
    int ibn_tmp;
    int itmp = 0;
    int iaddr_tmp=0;
    
    for(int i=0; i < group; i++){
        for(int ss = 0;ss < radix; ss++){
            ibn_tmp = 0;
            ima_tmp = 0; 
            itmp = ss * group + i; 
        //bit calculate
            for(int j=0; j < fft_point_bit_size;j++){
                ibit_tmp = itmp % 2;
                itmp = itmp >> 1;
                ibit_array_tmp[j] = ibit_tmp;
            }
            //radix-8 bit-reverse
			for(int br_index=0;br_index < (bc_width/3); br_index++){
				ibit_br[(bc_width-1) - 3*br_index     ]= ibit_array_tmp[3 * br_index + 2]; 
				ibit_br[(bc_width-1) - 3*br_index - 1 ]= ibit_array_tmp[3 * br_index + 1]; 
				ibit_br[(bc_width-1) - 3*br_index - 2 ]= ibit_array_tmp[3 * br_index + 0]; 
			}
			
            //bit-reverse over
            for(int rs = 0; rs < bc_width; rs++){
                if((ibit_br[rs] == 1) && (rs != 0)) ima_tmp = ima_tmp + exp2((rs-1)); 
                ibn_tmp = ibn_tmp ^ ibit_br[rs];
            }
            for(int kk=0; kk< radix; kk++){
                iaddr_tmp = (ss*group) + i + (kk * (offset));
                if(iaddr_tmp >= fft_point)std::cout << iaddr_tmp <<"\n";
                ima[iaddr_tmp] = ima_tmp;
                ibn[iaddr_tmp] = ibn_tmp;  
            }
        }   
    }
    
    ofstream index_file_ibn0;
    ofstream index_file_ibn1;
    index_file_ibn0.open("./SPMB/SPMB_ibn0.txt");
    index_file_ibn1.open("./SPMB/SPMB_ibn1.txt");
    for(int i=0;i < fft_point; i++){
        if(ibn[i]==0) index_file_ibn0 <<" index:" << i << ", MA:" << ima[i] <<"\n";
        else index_file_ibn1 <<" index:" << i << ", MA:" << ima[i] << "\n";
    }
    //=========================================================================================================
    //ROM DATA Generate 
	// rom word size : total worde size is divided  into 4 banks  
    std::cout << "*************************************************\n";
	std::cout << "IMPORTANT!!! ROM Data Generate command ===> ";
    std::cout << " make R8_ALL \n";
	std::cout << "**************************************************\n";
    r8_FFT_TW_ROM();
    r8_IFFT_TW_ROM();

}
void SPMB::time_o_r8(std::vector<ZZ> time_data,std::string string_in){
    ofstream b0radix0(string_in + "b0radix0.txt");
    ofstream b0radix1(string_in + "b0radix1.txt");
    ofstream b0radix2(string_in + "b0radix2.txt");
    ofstream b0radix3(string_in + "b0radix3.txt");
    ofstream b0radix4(string_in + "b0radix4.txt");
    ofstream b0radix5(string_in + "b0radix5.txt");
    ofstream b0radix6(string_in + "b0radix6.txt");
    ofstream b0radix7(string_in + "b0radix7.txt");
    ofstream b1radix0(string_in + "b1radix0.txt");
    ofstream b1radix1(string_in + "b1radix1.txt");
    ofstream b1radix2(string_in + "b1radix2.txt");
    ofstream b1radix3(string_in + "b1radix3.txt");
    ofstream b1radix4(string_in + "b1radix4.txt");
    ofstream b1radix5(string_in + "b1radix5.txt");
    ofstream b1radix6(string_in + "b1radix6.txt");
    ofstream b1radix7(string_in + "b1radix7.txt");
   
    std::string string_buf; 
    
    ZZ time_data_tmp; 

    int counter =0;
    int loop_index = 0;
    int loop_address = 0;
 
    //bank0
    while(counter < counter_iteration){ 
       if((ma[loop_index] == counter )&& (bn[loop_index] == 0) ){
          for(int j=0;j<radix;j++){
            loop_address = loop_index + offset * j;
            time_data_tmp = time_data[loop_address];
            string_buf  = ZZtohex_cp(time_data_tmp);
            
            if(j==0) {b0radix0 << string_buf;b0radix0 << "\n";}     
            if(j==1) {b0radix1 << string_buf;b0radix1 << "\n";}
            if(j==2) {b0radix2 << string_buf;b0radix2 << "\n";} 
            if(j==3) {b0radix3 << string_buf;b0radix3 << "\n";}
            if(j==4) {b0radix4 << string_buf;b0radix4 << "\n";}     
            if(j==5) {b0radix5 << string_buf;b0radix5 << "\n";}
            if(j==6) {b0radix6 << string_buf;b0radix6 << "\n";} 
            if(j==7) {b0radix7 << string_buf;b0radix7 << "\n";}			
          }           
           loop_index = 0;
           counter = counter + 1;
       }
       else loop_index = loop_index + 1;
    }
  
    counter      = 0;
    loop_index   = 0;
    loop_address = 0;
    //bank 1
    while(counter < counter_iteration){ 
        if((ma[loop_index] == counter )&& (bn[loop_index] == 1) ){
            for(int j=0;j<radix;j++){
               loop_address = loop_index + offset*j;
               time_data_tmp = time_data[loop_address];
               string_buf  = ZZtohex_cp(time_data_tmp);
               
               if(j==0) {b1radix0 << string_buf;b1radix0 << "\n";}     
               if(j==1) {b1radix1 << string_buf;b1radix1 << "\n";}
               if(j==2) {b1radix2 << string_buf;b1radix2 << "\n";} 
               if(j==3) {b1radix3 << string_buf;b1radix3 << "\n";}
			   if(j==4) {b1radix4 << string_buf;b1radix4 << "\n";}
			   if(j==5) {b1radix5 << string_buf;b1radix5 << "\n";}
			   if(j==6) {b1radix6 << string_buf;b1radix6 << "\n";}
			   if(j==7) {b1radix7 << string_buf;b1radix7 << "\n";}
            }           
            loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
    
    //h2_time_of.close();
    b0radix0.close();
    b0radix1.close();
    b0radix2.close();
    b0radix3.close();
	b0radix4.close();
    b0radix5.close();
    b0radix6.close();
    b0radix7.close();
    b1radix0.close();
    b1radix1.close();
    b1radix2.close();
    b1radix3.close();
	b1radix4.close();
    b1radix5.close();
    b1radix6.close();
    b1radix7.close();
}
void SPMB::re_order_factor_r8(ZZ m_2_rou){

    std::vector<ZZ> re_order_factor_array;
    re_order_factor_array.resize(fft_point);
    std::ofstream  ReROM_o("./ROM_Data/ReORDER.txt");    
	
    for(unsigned long i=0;i<m;i++){
       unsigned long exp;
       ZZ w_2m_tmp_i;
       exp = pow(i,2);
       exp = exp % (2*m);
       PowerMod(w_2m_tmp_i,m_2_rou,exp,cyclotomic_prime);
       re_order_factor_array[i] = w_2m_tmp_i;
	   //testing 2020/07/04

	   ReROM_o << re_order_factor_array[i] << "\n";
    }
    

	//cyclotomic polynomial prime must be 22~24. 
    //ROM 0 bit size = 64  // stroing bank0 radix-0 and  radix-1
	//ROM 1 bit size = 64 
	//re_order_factor_array,while index > m then data is ZERO.
	//because using SPMB Memory addressing,then this reorder factor  
	//=> ROM0[63:42]=bank0_radix0,ROM0[41:20]=bank0_radix1,ROM0[19:0]=20'b0;
    //=> ROM1[63:42]=bank1_radix0,ROM1[41:20]=bank1_radix1,ROM1[19:0]=20'b0;
	//SPMB Memory addressing, for N = 4096  
	//bank0 MA[0]: index 0 , index 1024 ,index 2048,index 3072
	//bank1 MA[0]: index 1 , index 1025 ,index 2049,index 3073
	
	//bank0
	std::ofstream reorder_ROM0("./ROM_Data/reorder_ROM0.txt");
	std::ofstream reorder_ROM1("./ROM_Data/reorder_ROM1.txt");
	//bank1
    std::ofstream reorder_ROM2("./ROM_Data/reorder_ROM2.txt");
    std::ofstream reorder_ROM3("./ROM_Data/reorder_ROM3.txt");
   
    std::vector<int> tobit_r0;
    std::vector<int> tobit_r1;
    std::vector<int> tobit_r2;
    std::vector<int> tobit_r3;

    tobit_r0.resize(CP_width);
    tobit_r1.resize(CP_width);
    tobit_r2.resize(CP_width);
    tobit_r3.resize(CP_width);
    
    ZZ re_order_r0_tmp;
    ZZ re_order_r1_tmp;
    ZZ re_order_r2_tmp;
    ZZ re_order_r3_tmp;
    
    int counter =0;
    int loop_index = 0;
    int loop_address = 0;

    //bank0
    while(counter < counter_iteration){ 
        if((ma[loop_index] == counter )&& (bn[loop_index] == 0) ){
            re_order_r0_tmp = re_order_factor_array[loop_index];
            re_order_r1_tmp = re_order_factor_array[loop_index + offset];
            re_order_r2_tmp = re_order_factor_array[loop_index + 2 * offset];
            re_order_r3_tmp = re_order_factor_array[loop_index + 3 * offset];
            
			//calculate bit
            for(int bit_index=0; bit_index < CP_width ;bit_index++){
                //radix0
				if(re_order_r0_tmp % 2 == 1) tobit_r0[bit_index] = 1;
                else tobit_r0[bit_index] = 0;
				//radix1
				if(re_order_r1_tmp % 2 == 1) tobit_r1[bit_index] = 1;
                else tobit_r1[bit_index] = 0;
				//radix-2
			    if(re_order_r2_tmp % 2 == 1) tobit_r2[bit_index] = 1;
                else tobit_r2[bit_index] = 0;
				//radix-3
				if(re_order_r3_tmp % 2 == 1) tobit_r3[bit_index] = 1;
                else tobit_r3[bit_index] = 0;
				
				re_order_r0_tmp = re_order_r0_tmp >> 1;
				re_order_r1_tmp = re_order_r1_tmp >> 1;
				re_order_r2_tmp = re_order_r2_tmp >> 1;
				re_order_r3_tmp = re_order_r3_tmp >> 1;
            }
            
           
            for(int g=0;g < CP_width; g++){
				reorder_ROM0 << tobit_r0[CP_width-1-g];
            } 
            for(int g=0;g < CP_width; g++){
                reorder_ROM0 << tobit_r1[CP_width-1-g];
            } 
            
			int ROM0_Remaining_bits ;
			
			ROM0_Remaining_bits = 64 - 2*CP_width;
            for(int g=0;g < CP_width; g++){
                if(ROM0_Remaining_bits > g)reorder_ROM0 << tobit_r2[CP_width-1-g] ;
				else reorder_ROM1 << tobit_r2[CP_width-1-g];
			}
            reorder_ROM0 << "\n";

			for(int g=0;g < CP_width; g++){
                reorder_ROM1 << tobit_r3[CP_width-1-g] ;
			}                                
			
			int ROM1_padding_zero;
			
            if((4*CP_width) > 96){
				ROM1_padding_zero = 64 - (2 * CP_width) + ROM0_Remaining_bits;
			    
				for(int g=0;g < ROM1_padding_zero; g++){
			        reorder_ROM1 << 0 ;
			    }				
			}
            if((4*CP_width)<= 96){ 
				ROM1_padding_zero = 32 - (2 * CP_width) + ROM0_Remaining_bits;
			    
				for(int g=0;g < ROM1_padding_zero; g++){
			        reorder_ROM1 << 0 ;
			    }				
			}
			reorder_ROM1 << "\n";
			
            loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
    reorder_ROM0.close();
    reorder_ROM1.close();
    counter = 0;
    loop_index = 0;
    //bank 1
    while(counter < counter_iteration){ 
        if((ma[loop_index] == counter )&& (bn[loop_index] == 1) ){
            re_order_r0_tmp = re_order_factor_array[loop_index];
            re_order_r1_tmp = re_order_factor_array[loop_index + offset];
            re_order_r2_tmp = re_order_factor_array[loop_index + 2 * offset];
            re_order_r3_tmp = re_order_factor_array[loop_index + 3 * offset];
            
			//calculate bit
            for(int bit_index=0; bit_index < CP_width ;bit_index++){
                //radix0
				if(re_order_r0_tmp % 2 == 1) tobit_r0[bit_index] = 1;
                else tobit_r0[bit_index] = 0;
				//radix1
				if(re_order_r1_tmp % 2 == 1) tobit_r1[bit_index] = 1;
                else tobit_r1[bit_index] = 0;
				//radix-2
			    if(re_order_r2_tmp % 2 == 1) tobit_r2[bit_index] = 1;
                else tobit_r2[bit_index] = 0;
				//radix-3
				if(re_order_r3_tmp % 2 == 1) tobit_r3[bit_index] = 1;
                else tobit_r3[bit_index] = 0;
				
				re_order_r0_tmp = re_order_r0_tmp >> 1;
				re_order_r1_tmp = re_order_r1_tmp >> 1;
				re_order_r2_tmp = re_order_r2_tmp >> 1;
				re_order_r3_tmp = re_order_r3_tmp >> 1;
            }
            
           
            for(int g=0;g < CP_width; g++){
				reorder_ROM2 << tobit_r0[CP_width-1-g];
            } 
            for(int g=0;g < CP_width; g++){
                reorder_ROM2 << tobit_r1[CP_width-1-g];
            } 
            
			int ROM2_Remaining_bits ;
			
			ROM2_Remaining_bits = 64 - 2*CP_width;
            for(int g=0;g < CP_width; g++){
                if(ROM2_Remaining_bits > g)reorder_ROM2 << tobit_r2[CP_width-1-g] ;
				else reorder_ROM3 << tobit_r2[CP_width-1-g];
			}
            reorder_ROM2 << "\n";

			for(int g=0;g < CP_width; g++){
                reorder_ROM3 << tobit_r3[CP_width-1-g] ;
			}                                
			
			int ROM3_padding_zero;
			
            if((4*CP_width) > 96){
				ROM3_padding_zero = 64 - (2 * CP_width) + ROM2_Remaining_bits;
			    
				for(int g=0;g < ROM3_padding_zero; g++){
			        reorder_ROM3 << 0 ;
			    }				
			}
            if((4*CP_width)< 96){ 
				ROM3_padding_zero = 32 - (2 * CP_width) + ROM2_Remaining_bits;
			    
				for(int g=0;g < ROM3_padding_zero; g++){
			        reorder_ROM3 << 0 ;
			    }				
			}
			reorder_ROM3 << "\n";
			
            loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
    reorder_ROM2.close();
    reorder_ROM3.close();

    loop_address = 0;
    
	//inverse order factor 
	std::vector<int> itobit_r0;
	std::vector<int> itobit_r1;
	std::vector<int> itobit_r2;
	std::vector<int> itobit_r3;
	std::vector<int> itobit_r4;
	std::vector<int> itobit_r5;
	std::vector<int> itobit_r6;
	std::vector<int> itobit_r7;
	itobit_r0.resize(CP_width);
	itobit_r1.resize(CP_width);
	itobit_r2.resize(CP_width);
	itobit_r3.resize(CP_width);
	itobit_r4.resize(CP_width);
	itobit_r5.resize(CP_width);
	itobit_r6.resize(CP_width);
	itobit_r7.resize(CP_width);
	ZZ re_iorder_r0_tmp;
	ZZ re_iorder_r1_tmp;
	ZZ re_iorder_r2_tmp;
	ZZ re_iorder_r3_tmp;
	ZZ re_iorder_r4_tmp;
	ZZ re_iorder_r5_tmp;
	ZZ re_iorder_r6_tmp;
	ZZ re_iorder_r7_tmp;
	/*************************************************/
	//m-th cyclotomic polynomial 
	//all polynomial at this thing,the coefficient is zero while the degree > phi(m);
    //and inverse reorder factor index > phi(m) is zero.
	//using SPMB (Final stage type) ,the data index is continuous
	//we can divide inverse re-order factor into two banks ROM.
	//ROM0 bit width= 64 ,ROM1 bit width = 32
	//cyclotomic prime = 22 bits
	//ROM0[63:42] = radix0 ,ROM0[41:20] = radix1,ROM0[19:0] = radix2[21:2]
	//ROM1[31:30] = raidx2[1:0],ROM1[29:8] = radix3,ROM1[7:0] = 8'b0;
	//ROM Word size ,original word size is N/radix => 4096/4 = 1024,
	//a_0x^(0)+........+a_phi(m)-1*x^(phi(m)-1)+0*x^phi(m)+...+0^x(m-1)
	//so word size can redution to phi(m)/radix 
    //for example m=1705 phi(m) = 1200 , then ma > 300 data is zero	
	//then we set word size is 512
	/*************************************************/
    std::ofstream ireorder_ROM0("./ROM_Data/ireorder_ROM0.txt");
    std::ofstream ireorder_ROM1("./ROM_Data/ireorder_ROM1.txt");
    std::ofstream ireorder_ROM2("./ROM_Data/ireorder_ROM2.txt");
    std::ofstream ireorder_ROM3("./ROM_Data/ireorder_ROM3.txt");
	
	
    double IReROM_WORD_BIT;
	int    ROM_WORD_SIZE;
    IReROM_WORD_BIT = ceil((double)m / 8);
    IReROM_WORD_BIT = log2(IReROM_WORD_BIT);
    IReROM_WORD_BIT = ceil(IReROM_WORD_BIT);	
    ROM_WORD_SIZE   = exp2(IReROM_WORD_BIT);
	    
	for(int i=0;i < ROM_WORD_SIZE; i++){
         re_iorder_r0_tmp = re_order_factor_array[radix * i + 0];
         re_iorder_r1_tmp = re_order_factor_array[radix * i + 1];
         re_iorder_r2_tmp = re_order_factor_array[radix * i + 2];
         re_iorder_r3_tmp = re_order_factor_array[radix * i + 3];
         re_iorder_r4_tmp = re_order_factor_array[radix * i + 4];
         re_iorder_r5_tmp = re_order_factor_array[radix * i + 5];
         re_iorder_r6_tmp = re_order_factor_array[radix * i + 6];
         re_iorder_r7_tmp = re_order_factor_array[radix * i + 7];
        
        for(int bit_index=0; bit_index < CP_width ;bit_index++){
            //radix0
			if(re_iorder_r0_tmp % 2 == 1) itobit_r0[bit_index] = 1;
            else itobit_r0[bit_index] = 0;
			//radix1
			if(re_iorder_r1_tmp % 2 == 1) itobit_r1[bit_index] = 1;
			else itobit_r1[bit_index] = 0;
			//radix2
			if(re_iorder_r2_tmp % 2 == 1) itobit_r2[bit_index] = 1;
			else itobit_r2[bit_index] = 0;
			//radix3
            if(re_iorder_r3_tmp % 2 == 1) itobit_r3[bit_index] = 1;
			else itobit_r3[bit_index] = 0;
			//radix4
            if(re_iorder_r4_tmp % 2 == 1) itobit_r4[bit_index] = 1;
			else itobit_r4[bit_index] = 0;			
			//radix5
            if(re_iorder_r5_tmp % 2 == 1) itobit_r5[bit_index] = 1;
			else itobit_r5[bit_index] = 0;			
			//radix6
            if(re_iorder_r6_tmp % 2 == 1) itobit_r6[bit_index] = 1;
			else itobit_r6[bit_index] = 0;
			//radix3
            if(re_iorder_r7_tmp % 2 == 1) itobit_r7[bit_index] = 1;
			else itobit_r7[bit_index] = 0;			
			re_iorder_r0_tmp = re_iorder_r0_tmp >> 1;
            re_iorder_r1_tmp = re_iorder_r1_tmp >> 1;
            re_iorder_r2_tmp = re_iorder_r2_tmp >> 1;
            re_iorder_r3_tmp = re_iorder_r3_tmp >> 1;
            re_iorder_r4_tmp = re_iorder_r4_tmp >> 1;
            re_iorder_r5_tmp = re_iorder_r5_tmp >> 1;
            re_iorder_r6_tmp = re_iorder_r6_tmp >> 1;
            re_iorder_r7_tmp = re_iorder_r7_tmp >> 1;
        }
     
        for(int g=0;g < CP_width; g++){
            ireorder_ROM0 << itobit_r0[CP_width-1-g];
        } 
         
        for(int g=0;g < CP_width; g++){
            ireorder_ROM0 << itobit_r1[CP_width-1-g];
        } 
        
		int IROM0_Remaining_bits;
		IROM0_Remaining_bits = 64 - 2 * CP_width;
		
        for(int g=0;g < CP_width; g++){
            if(IROM0_Remaining_bits > g)ireorder_ROM0 << itobit_r2[CP_width-1-g];
			else ireorder_ROM1 << itobit_r2[CP_width-1-g];
        } 
        ireorder_ROM0 << "\n"; 
        
        for(int g=0;g < CP_width; g++){
            ireorder_ROM1 << itobit_r3[CP_width-1-g];
        } 
		
		for(int g=0;g < CP_width; g++){
            ireorder_ROM1 << itobit_r4[CP_width-1-g];
        } 
		
		int IROM1_Remaining_bits;
		IROM1_Remaining_bits = 64 - 3 * CP_width + IROM0_Remaining_bits;
		for(int g=0;g < CP_width; g++){
            if(IROM1_Remaining_bits > g)ireorder_ROM1 << itobit_r5[CP_width-1-g];
			else ireorder_ROM2 << itobit_r5[CP_width-1-g];
        }		
        ireorder_ROM1 << "\n"; 
		//ROM2[59:38] = radix6
		for(int g=0;g < CP_width; g++){
            ireorder_ROM2 << itobit_r6[CP_width-1-g];
        } 
		
		if(8*CP_width < 192){
		   for(int g=0;g < CP_width; g++){
              ireorder_ROM2 << itobit_r7[CP_width-1-g];
           } 
	       int padding_zero_rom2;
		   padding_zero_rom2 = 64 - 3 * CP_width + IROM1_Remaining_bits;
		   for(int g=0; g < padding_zero_rom2;g++){
			   ireorder_ROM2 << "0";
		   }
		}
		if(8*CP_width == 192){
		    for(int g=0;g < CP_width; g++){
              ireorder_ROM2 << itobit_r7[CP_width-1-g];
            } 
		}
		if(8 * CP_width > 192){
			for(int g=0;g < CP_width-8; g++){
              ireorder_ROM2 << itobit_r7[CP_width-1-g];
            } 
			for(int g=0;g<8;g++){
			  ireorder_ROM3 << itobit_r7[7-g];
			}
			ireorder_ROM3 << "\n";
		}
		ireorder_ROM2 << "\n";
    }
    ireorder_ROM0.close();
    ireorder_ROM1.close();
    ireorder_ROM2.close();
    ireorder_ROM3.close();
}
void SPMB::H_freq_o_r8(std::vector<ZZ> H_NTT){
    //bank 0
    std::ofstream H_freq_out("./test_input/H_freq_out.txt");
    std::ofstream H_b0ROM0("./ROM_Data/H_b0ROM0.txt"); // radix0 ,radix1
    std::ofstream H_b0ROM1("./ROM_Data/H_b0ROM1.txt"); //radix2 radix3
	std::ofstream H_b0ROM2("./ROM_Data/H_b0ROM2.txt"); // radix4 ,radix5
	std::ofstream H_b0ROM3("./ROM_Data/H_b0ROM3.txt"); //radix6 radix7
    //bank1                     ROM_Data
    std::ofstream H_b1ROM0("./ROM_Data/H_b1ROM0.txt"); //radix0 radix1
    std::ofstream H_b1ROM1("./ROM_Data/H_b1ROM1.txt"); // radix2 radix3
    std::ofstream H_b1ROM2("./ROM_Data/H_b1ROM2.txt"); // radix4 radix5
    std::ofstream H_b1ROM3("./ROM_Data/H_b1ROM3.txt"); // radix6 radix7
    
    std::vector<int> tobit_r0;
    std::vector<int> tobit_r1;
    std::vector<int> tobit_r2;
    std::vector<int> tobit_r3;
    std::vector<int> tobit_r4;
    std::vector<int> tobit_r5;
    std::vector<int> tobit_r6;
    std::vector<int> tobit_r7;
    tobit_r0.resize(64);
    tobit_r1.resize(64);
    tobit_r2.resize(64);
    tobit_r3.resize(64);
    tobit_r4.resize(64);
    tobit_r5.resize(64);
    tobit_r6.resize(64);
    tobit_r7.resize(64);
    
    ZZ H_NTT_r0;
    ZZ H_NTT_r1;
    ZZ H_NTT_r2;
    ZZ H_NTT_r3;
    ZZ H_NTT_r4;
    ZZ H_NTT_r5;
    ZZ H_NTT_r6;
    ZZ H_NTT_r7;
    
    int counter =0;
    int loop_index = 0;
    
    for(int i=0;i<fft_point;i++){
		H_freq_out << H_NTT[i];
        H_freq_out << "\n";
    }
    
    //bank0
    while(counter < counter_iteration){ 
        if((ima[loop_index] == counter )&& (ibn[loop_index] == 0) ){
			
            H_NTT_r0 = H_NTT[loop_index];
            H_NTT_r1 = H_NTT[loop_index + offset];
            H_NTT_r2 = H_NTT[loop_index + 2 * offset];
            H_NTT_r3 = H_NTT[loop_index + 3 * offset];
            H_NTT_r4 = H_NTT[loop_index + 4 * offset];
            H_NTT_r5 = H_NTT[loop_index + 5 * offset];
            H_NTT_r6 = H_NTT[loop_index + 6 * offset];
            H_NTT_r7 = H_NTT[loop_index + 7 * offset];
            
            for(int bit_index=0; bit_index < 64 ;bit_index++){
                if(H_NTT_r0 % 2 == 1) tobit_r0[bit_index] = 1;
                else tobit_r0[bit_index] = 0;
				if(H_NTT_r1 % 2 == 1) tobit_r1[bit_index] = 1;
                else tobit_r1[bit_index] = 0;
				if(H_NTT_r2 % 2 == 1) tobit_r2[bit_index] = 1;
                else tobit_r2[bit_index] = 0;
				if(H_NTT_r3 % 2 == 1) tobit_r3[bit_index] = 1;
                else tobit_r3[bit_index] = 0;
				if(H_NTT_r4 % 2 == 1) tobit_r4[bit_index] = 1;
                else tobit_r4[bit_index] = 0;
				if(H_NTT_r5 % 2 == 1) tobit_r5[bit_index] = 1;
                else tobit_r5[bit_index] = 0;
				if(H_NTT_r6 % 2 == 1) tobit_r6[bit_index] = 1;
                else tobit_r6[bit_index] = 0;
				if(H_NTT_r7 % 2 == 1) tobit_r7[bit_index] = 1;
                else tobit_r7[bit_index] = 0;				
                H_NTT_r0 = H_NTT_r0 >> 1;
                H_NTT_r1 = H_NTT_r1 >> 1;
                H_NTT_r2 = H_NTT_r2 >> 1;
                H_NTT_r3 = H_NTT_r3 >> 1;
                H_NTT_r4 = H_NTT_r4 >> 1;
                H_NTT_r5 = H_NTT_r5 >> 1;
                H_NTT_r6 = H_NTT_r6 >> 1;
                H_NTT_r7 = H_NTT_r7 >> 1;
            }
            
            for(int g=0;g < 64; g++){
                H_b0ROM0 << tobit_r0[63-g];
            } 

            for(int g=0;g < 64; g++){
                H_b0ROM0 << tobit_r1[63-g];
            }
			
			H_b0ROM0 << "\n";
            
            for(int g=0;g < 64; g++){
                H_b0ROM1 << tobit_r2[63-g];
            }
			
            for(int g=0;g < 64; g++){   
                H_b0ROM1 << tobit_r3[63-g];
            }
            H_b0ROM1 << "\n";
			
            for(int g=0;g < 64; g++){
                H_b0ROM2 << tobit_r4[63-g];
            }
            for(int g=0;g < 64; g++){   
                H_b0ROM2 << tobit_r5[63-g];
            }
            H_b0ROM2 << "\n";				
            
			
            for(int g=0;g < 64; g++){
                H_b0ROM3 << tobit_r6[63-g];
            }
            for(int g=0;g < 64; g++){   
                H_b0ROM3 << tobit_r7[63-g];
            }
            H_b0ROM3 << "\n";				
            			
			loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
    
    counter = 0;
    loop_index = 0;
    
	//bank 1
    while(counter < counter_iteration){ 
        if((ima[loop_index] == counter )&& (ibn[loop_index] == 1) ){
            
			H_NTT_r0 = H_NTT[loop_index];
            H_NTT_r1 = H_NTT[loop_index + offset];
            H_NTT_r2 = H_NTT[loop_index + 2 * offset];
            H_NTT_r3 = H_NTT[loop_index + 3 * offset];
            H_NTT_r4 = H_NTT[loop_index + 4 * offset];
            H_NTT_r5 = H_NTT[loop_index + 5 * offset];
            H_NTT_r6 = H_NTT[loop_index + 6 * offset];
            H_NTT_r7 = H_NTT[loop_index + 7 * offset];
			
            for(int bit_index=0; bit_index < 64 ;bit_index++){
                if(H_NTT_r0 % 2 == 1) tobit_r0[bit_index] = 1;
                else tobit_r0[bit_index] = 0;
				if(H_NTT_r1 % 2 == 1) tobit_r1[bit_index] = 1;
                else tobit_r1[bit_index] = 0;
				if(H_NTT_r2 % 2 == 1) tobit_r2[bit_index] = 1;
                else tobit_r2[bit_index] = 0;
				if(H_NTT_r3 % 2 == 1) tobit_r3[bit_index] = 1;
                else tobit_r3[bit_index] = 0;
				if(H_NTT_r4 % 2 == 1) tobit_r4[bit_index] = 1;
                else tobit_r4[bit_index] = 0;				
				if(H_NTT_r5 % 2 == 1) tobit_r5[bit_index] = 1;
                else tobit_r5[bit_index] = 0;
				if(H_NTT_r6 % 2 == 1) tobit_r6[bit_index] = 1;
                else tobit_r6[bit_index] = 0;
				if(H_NTT_r7 % 2 == 1) tobit_r7[bit_index] = 1;
                else tobit_r7[bit_index] = 0;				
                H_NTT_r0 = H_NTT_r0 >> 1;
                H_NTT_r1 = H_NTT_r1 >> 1;
                H_NTT_r2 = H_NTT_r2 >> 1;
                H_NTT_r3 = H_NTT_r3 >> 1;
                H_NTT_r4 = H_NTT_r4 >> 1;
                H_NTT_r5 = H_NTT_r5 >> 1;
                H_NTT_r6 = H_NTT_r6 >> 1;
                H_NTT_r7 = H_NTT_r7 >> 1;
            }
            
            
            for(int g=0;g < 64; g++){
                H_b1ROM0 << tobit_r0[63-g];
            } 
            
            for(int g=0;g < 64; g++){
                H_b1ROM0 << tobit_r1[63-g];
            }
            
			H_b1ROM0 << "\n"; 
            
            for(int g=0;g < 64; g++){
                H_b1ROM1 << tobit_r2[63-g];
            }
                    
            for(int g=0;g < 64; g++){   
                H_b1ROM1 << tobit_r3[63-g];
            }
    
            H_b1ROM1 << "\n";           
            
            for(int g=0;g < 64; g++){
                H_b1ROM2 << tobit_r4[63-g];
            }       
            for(int g=0;g < 64; g++){   
                H_b1ROM2 << tobit_r5[63-g];
            }
            H_b1ROM2 << "\n";    
            
            for(int g=0;g < 64; g++){
                H_b1ROM3 << tobit_r6[63-g];
            }       
            for(int g=0;g < 64; g++){   
                H_b1ROM3 << tobit_r7[63-g];
            }
            H_b1ROM3 << "\n";  
            
			loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
	H_b0ROM0.close();
	H_b0ROM1.close();
	H_b0ROM2.close();
	H_b0ROM3.close();
	
	H_b1ROM0.close();	
	H_b1ROM1.close();	
	H_b1ROM2.close();	
	H_b1ROM3.close();	
	
	
}
void SPMB::r8_FFT_TW_ROM(){

    //---------siang print data -----------
    std::ofstream siang_ROM0("./my_print_data/siang_R8_FFTROM0.txt");
    std::ofstream siang_ROM1("./my_print_data/siang_R8_FFTROM1.txt");
    std::ofstream siang_ROM2("./my_print_data/siang_R8_FFTROM2.txt");
    std::ofstream siang_ROM3("./my_print_data/siang_R8_FFTROM3.txt");
    //----------

	//ROM0 bitsize = 64 , ROM1 bitsize = 128
	//word size = FFT_point / radix
    std::ofstream  ROM0("./ROM_Data/R8_FFTROM0.txt");
    std::ofstream  ROM1("./ROM_Data/R8_FFTROM1.txt");
    std::ofstream  ROM2("./ROM_Data/R8_FFTROM2.txt");
    std::ofstream  ROM3("./ROM_Data/R8_FFTROM3.txt");
    //FFT Prime
    ZZ FFT_Prime;
    ZZ twiddle_65536;
    ZZ R1_tmp;
    ZZ R2_tmp;
    ZZ R3_tmp;
    ZZ R4_tmp;
    ZZ R5_tmp;
    ZZ R6_tmp;
    ZZ R7_tmp;
    long exp;
    long order;  // order = 65536 / fft_point
    long addr_length; // FFT_point / radix
    conv(FFT_Prime,"18446744069414584321");     //FFT_Prime = 2^64 - 2^32 +1
    conv(twiddle_65536,"14603442835287214144"); //65536-th root of unity 
    order = (long) 65536 / fft_point;
    addr_length  = (long) fft_point / radix;

    std::vector<int> tobit_r1;
    std::vector<int> tobit_r2;
    std::vector<int> tobit_r3;
    std::vector<int> tobit_r4;
    std::vector<int> tobit_r5;
    std::vector<int> tobit_r6;
    std::vector<int> tobit_r7;
    tobit_r1.resize(64);
    tobit_r2.resize(64);
    tobit_r3.resize(64);
    tobit_r4.resize(64);
    tobit_r5.resize(64);
    tobit_r6.resize(64);
    tobit_r7.resize(64);
    
    std::vector<ZZ> tw_Table;  //twiddle factor table 
    tw_Table.resize(fft_point);
    
    exp = 0;
    for(int i=0;i<fft_point;i++){
       tw_Table[i] = PowerMod(twiddle_65536,exp,FFT_Prime);
       exp = exp + order;
    }
    
    //radix 0
    for(int i=0;i<(addr_length);i++){
        R1_tmp = tw_Table[i];
        R2_tmp = tw_Table[2*i];
        R3_tmp = tw_Table[3*i];
        R4_tmp = tw_Table[4*i];
        R5_tmp = tw_Table[5*i];
        R6_tmp = tw_Table[6*i];
        R7_tmp = tw_Table[7*i];


        //siang twiddle print out
        //i=128
        siang_ROM0 << "exp = "  << i*order    << ", " << R1_tmp;
        siang_ROM1 << "exp = "  << 2*i*order  << ", " << R2_tmp;
        siang_ROM2 << "exp = "  << 4*i*order  << ", " << R4_tmp;
        siang_ROM3 << "exp = "  << 6*i*order  << ", " << R6_tmp;

        siang_ROM1 << " ,   ";
        siang_ROM2 << " ,   ";
        siang_ROM3 << " ,   ";

        siang_ROM1 << "exp = "  << 3*i*order << ", " << R3_tmp;
        siang_ROM2 << "exp = "  << 5*i*order << ", " << R5_tmp;
        siang_ROM3 << "exp = "  << 7*i*order << ", " << R7_tmp;
        
        siang_ROM0 << "\n";
        siang_ROM1 << "\n";
        siang_ROM2 << "\n";
        siang_ROM3 << "\n";
        //------------------------  




        for(int bit_index=0; bit_index < 64 ;bit_index++){
            if(R1_tmp % 2 == 1) tobit_r1[bit_index] = 1;
            else tobit_r1[bit_index] = 0;
            if(R2_tmp % 2 == 1) tobit_r2[bit_index] = 1;
            else tobit_r2[bit_index] = 0;             
            if(R3_tmp % 2 == 1) tobit_r3[bit_index] = 1;
            else tobit_r3[bit_index] = 0;             
            if(R4_tmp % 2 == 1) tobit_r4[bit_index] = 1;
			else tobit_r4[bit_index] = 0;             
			if(R5_tmp % 2 == 1) tobit_r5[bit_index] = 1;
			else tobit_r5[bit_index] = 0;
 			if(R6_tmp % 2 == 1) tobit_r6[bit_index] = 1;
 			else tobit_r6[bit_index] = 0;
			if(R7_tmp % 2 == 1) tobit_r7[bit_index] = 1;
			else tobit_r7[bit_index] = 0;
			
            R1_tmp = R1_tmp >> 1;
            R2_tmp = R2_tmp >> 1;
            R3_tmp = R3_tmp >> 1;
			R4_tmp = R4_tmp >> 1;
			R5_tmp = R5_tmp >> 1;
			R6_tmp = R6_tmp >> 1;
			R7_tmp = R7_tmp >> 1;
        }
        
        for(int g=0;g < 64; g++){
            ROM0 << tobit_r1[63-g];
            ROM1 << tobit_r2[63-g];
            ROM2 << tobit_r4[63-g];
            ROM3 << tobit_r6[63-g];
			
        }
        for(int g=0;g < 64; g++){
            ROM1 << tobit_r3[63-g];
            ROM2 << tobit_r5[63-g];
            ROM3 << tobit_r7[63-g];
        }
        ROM0 << "\n";
        ROM1 << "\n";
        ROM2 << "\n";
        ROM3 << "\n";
    }
    ROM0.close();
    ROM1.close();
    ROM2.close();
    ROM3.close();

    //---------------
    siang_ROM0.close();
	siang_ROM1.close();
	siang_ROM2.close();
	siang_ROM3.close();
    //----------------
}
void SPMB::r8_IFFT_TW_ROM(){
    
	//Bank0
    std::ofstream  ROM0("./ROM_Data/R8_IFFTROM0.txt");                          
    std::ofstream  ROM1("./ROM_Data/R8_IFFTROM1.txt");
    std::ofstream  ROM2("./ROM_Data/R8_IFFTROM2.txt");
    std::ofstream  ROM3("./ROM_Data/R8_IFFTROM3.txt");
    //FFT Prime
    ZZ FFT_Prime;
    ZZ twiddle_65536;
	ZZ Inverse_twiddle_65536;
	ZZ inverse_twiddle_order;
    ZZ R1_tmp;
    ZZ R2_tmp;
    ZZ R3_tmp;
    ZZ R4_tmp;
    ZZ R5_tmp;
    ZZ R6_tmp;
    ZZ R7_tmp;
    //long exp;
    long order;  // order = 65536 / fft_point
    long addr_length; // FFT_point / radix
    conv(FFT_Prime,"18446744069414584321");     //FFT_Prime = 2^64 - 2^32 +1
    conv(twiddle_65536,"14603442835287214144"); //65536-th root of unity 
	InvMod(Inverse_twiddle_65536,twiddle_65536,FFT_Prime);
    order = (long) 65536 / fft_point;
    addr_length  = (long) fft_point / radix;


    inverse_twiddle_order = PowerMod(Inverse_twiddle_65536,order,FFT_Prime);
	
    std::vector<int> tobit_r1;
    std::vector<int> tobit_r2;
    std::vector<int> tobit_r3;
    std::vector<int> tobit_r4;
    std::vector<int> tobit_r5;
    std::vector<int> tobit_r6;
    std::vector<int> tobit_r7;
    tobit_r1.resize(64);
    tobit_r2.resize(64);
    tobit_r3.resize(64);
    tobit_r4.resize(64);
    tobit_r5.resize(64);
    tobit_r6.resize(64);
    tobit_r7.resize(64);
    
    std::vector<ZZ> tw_Table;  //twiddle factor table 
    tw_Table.resize(fft_point);
    
    for(long i=0;i<fft_point;i++){
       tw_Table[i] = PowerMod(inverse_twiddle_order,i,FFT_Prime);
    }
    
    for(int i=0;i<(addr_length);i++){
        
        R1_tmp = tw_Table[i];
        R2_tmp = tw_Table[2 * i];
        R3_tmp = tw_Table[3 * i];
        R4_tmp = tw_Table[4 * i];
        R5_tmp = tw_Table[5 * i];
        R6_tmp = tw_Table[6 * i];
        R7_tmp = tw_Table[7 * i];
        
        for(int bit_index=0; bit_index < 64 ;bit_index++){
            if(R1_tmp % 2 == 1) tobit_r1[bit_index] = 1;
            else tobit_r1[bit_index] = 0;
            
            if(R2_tmp % 2 == 1) tobit_r2[bit_index] = 1;
            else tobit_r2[bit_index] = 0;             

            if(R3_tmp % 2 == 1) tobit_r3[bit_index] = 1;
            else tobit_r3[bit_index] = 0;             

            if(R4_tmp % 2 == 1) tobit_r4[bit_index] = 1;
            else tobit_r4[bit_index] = 0;  

            if(R5_tmp % 2 == 1) tobit_r5[bit_index] = 1;
            else tobit_r5[bit_index] = 0;  

            if(R6_tmp % 2 == 1) tobit_r6[bit_index] = 1;
            else tobit_r6[bit_index] = 0;

            if(R7_tmp % 2 == 1) tobit_r7[bit_index] = 1;
            else tobit_r7[bit_index] = 0;  			
            R1_tmp = R1_tmp >> 1;
            R2_tmp = R2_tmp >> 1;
            R3_tmp = R3_tmp >> 1;
            R4_tmp = R4_tmp >> 1;
            R5_tmp = R5_tmp >> 1;
            R6_tmp = R6_tmp >> 1;
            R7_tmp = R7_tmp >> 1;
        }
        
        for(int g=0;g < 64; g++){
            ROM0 << tobit_r1[63-g];
            ROM1 << tobit_r2[63-g];
            ROM2 << tobit_r4[63-g];
            ROM3 << tobit_r6[63-g];
        }
        for(int g=0;g < 64; g++){
            ROM1 << tobit_r3[63-g];
            ROM2 << tobit_r5[63-g];
            ROM3 << tobit_r7[63-g];
        }
        ROM0 << "\n";
        ROM1 << "\n";
		ROM2 << "\n";
        ROM3 << "\n";
		
    }
    ROM0.close();
    ROM1.close();
	ROM2.close();
	ROM3.close();
}
//===================================================================================================================
//radix-16
void SPMB::init_r16(unsigned long fft_p , unsigned long r,int bc_w,unsigned long CP_w,
	ZZ cyclotomic_p,long m_th,long phi_m,int IsReconfig){
    fft_point        = fft_p;
    radix            = r;
    bc_width         = bc_w;
	CP_width         = CP_w;
	cyclotomic_prime = cyclotomic_p;
	m                = m_th;
	phim             = phi_m;
    
    offset             =  (int)fft_point/radix;
    group              =  (int)fft_point/(radix*radix);
    fft_point_bit_size =  (int)log2(fft_point); 
    counter_iteration = fft_point/(radix*2);
    ma.resize(fft_point);
    bn.resize(fft_point);
    bit_array_tmp.resize(fft_point_bit_size);
    
    int ma_tmp;
    int bit_tmp;
    int bn_tmp;
    int tmp = 0;
    int addr_tmp = 0;
    
    for(int i=0; i < group; i++){        
        for(int ss=0 ; ss < radix; ss++){
            bn_tmp = 0;
            ma_tmp = 0;
            tmp = ss * group + i;    
            //bit calculate
            for(int j=0; j < fft_point_bit_size;j++){
                bit_tmp = tmp % 2;
                tmp = tmp >> 1;
                bit_array_tmp[j] = bit_tmp;
            } 
            for(int rs = 0; rs < fft_point_bit_size; rs++){
                if((bit_array_tmp[rs] == 1) && (rs != 0)) ma_tmp = ma_tmp + exp2((rs-1)); 
                bn_tmp = bn_tmp ^ bit_array_tmp[rs];
            }
            for(int kk=0; kk<radix; kk++){
                addr_tmp = ss*group + i + kk * (offset);
                if(addr_tmp >= fft_point)std::cout << addr_tmp <<"\n";
                ma[addr_tmp] = ma_tmp;
                bn[addr_tmp] = bn_tmp;
            }
        }    
    }
    
    ofstream index_file_bn0;
    ofstream index_file_bn1;
    index_file_bn0.open("./SPMB/SPMB_bn0.txt");
    index_file_bn1.open("./SPMB/SPMB_bn1.txt");
    for(int i=0;i < fft_point; i++){
       if(bn[i]==0) index_file_bn0 <<" index:" << i << ", MA:" << ma[i] <<"\n";
       else index_file_bn1 <<" index:" << i << ", MA:" << ma[i] << "\n";
    }
    
   
    //calculate inverse SPMB bn and ma 
    ima.resize(fft_point);
    ibn.resize(fft_point);
    ibit_br.resize(bc_width);
    ibit_array_tmp.resize(fft_point_bit_size);
    int ima_tmp;
    int ibit_tmp;
    int ibn_tmp;
    int itmp = 0;
    int iaddr_tmp=0;
    
    for(int i=0; i < group; i++){
        for(int ss = 0;ss < radix; ss++){
            ibn_tmp = 0;
            ima_tmp = 0; 
            itmp = ss * group + i; 
        //bit calculate
            for(int j=0; j < fft_point_bit_size;j++){
                ibit_tmp = itmp % 2;
                itmp = itmp >> 1;
                ibit_array_tmp[j] = ibit_tmp;
            }
			//radix-16 , then 4-bits 
			for(int br_index=0;br_index < (bc_width/4); br_index++){
				ibit_br[(bc_width-1) - 4*br_index     ]= ibit_array_tmp[4 * br_index + 3]; 
				ibit_br[(bc_width-1) - 4*br_index - 1 ]= ibit_array_tmp[4 * br_index + 2]; 
				ibit_br[(bc_width-1) - 4*br_index - 2 ]= ibit_array_tmp[4 * br_index + 1]; 
				ibit_br[(bc_width-1) - 4*br_index - 3 ]= ibit_array_tmp[4 * br_index + 0]; 
			}
            //bit-reverse over
            for(int rs = 0; rs < bc_width; rs++){
                if((ibit_br[rs] == 1) && (rs != 0)) ima_tmp = ima_tmp + exp2((rs-1)); 
                ibn_tmp = ibn_tmp ^ ibit_br[rs];
            }
            for(int kk=0; kk< radix; kk++){
                iaddr_tmp = (ss*group) + i + (kk * (offset));
                if(iaddr_tmp >= fft_point)std::cout << iaddr_tmp <<"\n";
                ima[iaddr_tmp] = ima_tmp;
                ibn[iaddr_tmp] = ibn_tmp;  
            }
        }   
    }
    
    ofstream index_file_ibn0;
    ofstream index_file_ibn1;
    index_file_ibn0.open("./SPMB/SPMB_ibn0.txt");
    index_file_ibn1.open("./SPMB/SPMB_ibn1.txt");
    for(int i=0;i < fft_point; i++){
        if(ibn[i]==0) index_file_ibn0 <<" index:" << i << ", MA:" << ima[i] <<"\n";
        else index_file_ibn1 <<" index:" << i << ", MA:" << ima[i] << "\n";
    }
    //=========================================================================================================
    //ROM DATA Generate 
        //=========================================================================================================
    //ROM DATA Generate 
	if( IsReconfig == 0){
         std::cout << "*************************************************************************\n";
		 std::cout << "IMPORTANT!!! SRAM Data Generate command ===> make SRAM\n";
	     std::cout << "IMPORTANT!!! ROM  Data Generate command ===> make ROM_R16\n";
	     std::cout << "IMPORTANT!!! SRAM Data ,ROM Data and simulation command ===> make R16_ALL\n";
	     std::cout << "**************************************************************************\n";
	     r16_FFT_TW_ROM();
         r16_IFFT_TW_ROM();
		 //------------------------------------------------------------------------------------------
	}else {
         std::cout << "*************************************************************************\n";
	     std::cout << "IMPORTANT!!! SRAM Data Generate command ===> make SRAM_CONFIG\n";
	     std::cout << "IMPORTANT!!! ROM  Data Generate command ===> make ROM_CONFIG_R16\n";
		 std::cout << "IMPORTANT!!! SRAM Data ,ROM Data and simulation command ===> make R16_CONFIG_ALL\n";
	     std::cout << "**************************************************************************\n";
		 r16_FFT_TW_ROM_Reconfig();
		 r16_IFFT_TW_ROM_Reconfig();
		 //------------------------------------------------------------------------------------------
	}
}
void SPMB::init_r16_Mixed_radix(unsigned long fft_p , unsigned long r,int bc_w,unsigned long CP_w,
	ZZ cyclotomic_p,long m_th,long phi_m,int IsReconfig){
    fft_point        = fft_p;
    radix            = r;
    bc_width         = bc_w;
	CP_width         = CP_w;
	cyclotomic_prime = cyclotomic_p;
	m                = m_th;
	phim             = phi_m;
    
    offset             =  (int)fft_point/radix;
    group              =  (int)fft_point/(radix*radix);
    fft_point_bit_size =  (int)log2(fft_point); 
    counter_iteration = fft_point/(radix*2);
    ma.resize(fft_point);
    bn.resize(fft_point);
    bit_array_tmp.resize(fft_point_bit_size);
    
    int ma_tmp;
    int bit_tmp;
    int bn_tmp;
    int tmp = 0;
    int addr_tmp = 0;
    std::cout << "--------------------------------------- \n";
	std::cout << "Now in init_r16_Mixed_radix function!!! \n";
	
    for(int i=0; i < group; i++){        
        for(int ss=0 ; ss < radix; ss++){
            bn_tmp = 0;
            ma_tmp = 0;
            tmp = ss * group + i;    
            //bit calculate
            for(int j=0; j < fft_point_bit_size;j++){
                bit_tmp = tmp % 2;
                tmp = tmp >> 1;
                bit_array_tmp[j] = bit_tmp;
            } 
            for(int rs = 0; rs < fft_point_bit_size; rs++){
                if((bit_array_tmp[rs] == 1) && (rs != 0)) ma_tmp = ma_tmp + exp2((rs-1)); 
                bn_tmp = bn_tmp ^ bit_array_tmp[rs];
            }
            for(int kk=0; kk<radix; kk++){
                addr_tmp = ss*group + i + kk * (offset);
                if(addr_tmp >= fft_point)std::cout << addr_tmp <<"\n";
                ma[addr_tmp] = ma_tmp;
                bn[addr_tmp] = bn_tmp;
            }
        }    
    }
    
    ofstream index_file_bn0;
    ofstream index_file_bn1;
    index_file_bn0.open("./SPMB/SPMB_bn0.txt");
    index_file_bn1.open("./SPMB/SPMB_bn1.txt");
    for(int i=0;i < fft_point; i++){
       if(bn[i]==0) index_file_bn0 <<" index:" << i << ", MA:" << ma[i] <<"\n";
       else index_file_bn1 <<" index:" << i << ", MA:" << ma[i] << "\n";
    }
    
   
    //calculate inverse SPMB bn and ma 
    ima.resize(fft_point);
    ibn.resize(fft_point);
    ibit_br.resize(bc_width);
    ibit_array_tmp.resize(fft_point_bit_size);
    int ima_tmp;
    int ibit_tmp;
    int ibn_tmp;
    int itmp = 0;
    int iaddr_tmp=0;

	for(int i = 0;i < fft_point;i++){
		ibn_tmp = 0;
		ima_tmp = 0;
		R16_Mixed_Radix_Pointwise_Mult_Index_AGU(i,ima_tmp,ibn_tmp);
		ima[i] = ima_tmp;
		ibn[i] = ibn_tmp;
	}
    

    ofstream index_file_ibn0;
    ofstream index_file_ibn1;
    index_file_ibn0.open("./SPMB/SPMB_ibn0.txt");
    index_file_ibn1.open("./SPMB/SPMB_ibn1.txt");
    for(int i=0;i < fft_point; i++){
        if(ibn[i]==0) index_file_ibn0 <<" index:" << i << ", MA:" << ima[i] <<"\n";
        else index_file_ibn1 <<" index:" << i << ", MA:" << ima[i] << "\n";
    }
    //=========================================================================================================
    //ROM DATA Generate 
        //=========================================================================================================
    //ROM DATA Generate 
    //=========================================================================================================
    //ROM DATA Generate 
        //=========================================================================================================
    //ROM DATA Generate 
	if( IsReconfig == 0){
         std::cout << "*************************************************************************\n";
		 std::cout << "IMPORTANT!!! SRAM Data Generate command ===> make SRAM\n";
	     std::cout << "IMPORTANT!!! ROM  Data Generate command ===> make ROM_R16\n";
		 std::cout << "IMPORTANT!!! SRAM Data ,ROM Data and simulation command ===> make R16_ALL\n";
	     std::cout << "**************************************************************************\n";
	     r16_FFT_TW_ROM();
		 
         r16_IFFT_TW_ROM();
		 //------------------------------------------------------------------------------------------
	}else {
         std::cout << "*************************************************************************\n";
	     std::cout << "IMPORTANT!!! SRAM Data Generate command ===> make SRAM_CONFIG\n";
	     std::cout << "IMPORTANT!!! ROM  Data Generate command ===> make ROM_CONFIG_R16\n";
		  std::cout << "IMPORTANT!!! SRAM Data ,ROM Data and simulation command ===> make R16_CONFIG_ALL\n";
	     std::cout << "**************************************************************************\n";
		 r16_FFT_TW_ROM_Reconfig();
		 std::cout << "Twiddle factor ROM generate over!!!!\n";
		 r16_IFFT_TW_ROM_Reconfig();
		 std::cout << "Inverse Twiddle factor ROM generate over!!!!\n";
		 //------------------------------------------------------------------------------------------
	}
}
void SPMB::time_o_r16(std::vector<ZZ> time_data,std::string string_in){
    ofstream b0radix0(string_in   +  "b0radix0.txt");
    ofstream b0radix1(string_in   +  "b0radix1.txt");
    ofstream b0radix2(string_in   +  "b0radix2.txt");
    ofstream b0radix3(string_in   +  "b0radix3.txt");
    ofstream b0radix4(string_in   +  "b0radix4.txt");
    ofstream b0radix5(string_in   +  "b0radix5.txt");
    ofstream b0radix6(string_in   +  "b0radix6.txt");
    ofstream b0radix7(string_in   +  "b0radix7.txt");
    ofstream b0radix8(string_in   +  "b0radix8.txt");
    ofstream b0radix9(string_in   +  "b0radix9.txt");
    ofstream b0radix10(string_in  +  "b0radix10.txt");
    ofstream b0radix11(string_in  +  "b0radix11.txt");
    ofstream b0radix12(string_in  +  "b0radix12.txt");
    ofstream b0radix13(string_in  +  "b0radix13.txt");
    ofstream b0radix14(string_in  +  "b0radix14.txt");
    ofstream b0radix15(string_in  +  "b0radix15.txt");
    ofstream b1radix0(string_in   +  "b1radix0.txt");
    ofstream b1radix1(string_in   +  "b1radix1.txt");
    ofstream b1radix2(string_in   +  "b1radix2.txt");
    ofstream b1radix3(string_in   +  "b1radix3.txt");
    ofstream b1radix4(string_in   +  "b1radix4.txt");
    ofstream b1radix5(string_in   +  "b1radix5.txt");
    ofstream b1radix6(string_in   +  "b1radix6.txt");
    ofstream b1radix7(string_in   +  "b1radix7.txt");
    ofstream b1radix8(string_in   +  "b1radix8.txt");
    ofstream b1radix9(string_in   +  "b1radix9.txt");
    ofstream b1radix10(string_in  +  "b1radix10.txt");
    ofstream b1radix11(string_in  +  "b1radix11.txt");
    ofstream b1radix12(string_in  +  "b1radix12.txt");
    ofstream b1radix13(string_in  +  "b1radix13.txt");
    ofstream b1radix14(string_in  +  "b1radix14.txt");
    ofstream b1radix15(string_in  +  "b1radix15.txt");
   
    std::string string_buf; 
    
    ZZ time_data_tmp; 

    int counter =0;
    int loop_index = 0;
    int loop_address = 0;
 
    //bank0
    while(counter < counter_iteration){ 
       if((ma[loop_index] == counter )&& (bn[loop_index] == 0) ){
          for(int j=0;j<radix;j++){
            loop_address   =  loop_index + offset*j;
            time_data_tmp  =  time_data[loop_address];
            string_buf     =  ZZtohex_cp(time_data_tmp);
            
            if(j==0 ) { b0radix0  << string_buf; b0radix0  << "\n"; }     
            if(j==1 ) { b0radix1  << string_buf; b0radix1  << "\n"; }
            if(j==2 ) { b0radix2  << string_buf; b0radix2  << "\n"; } 
            if(j==3 ) { b0radix3  << string_buf; b0radix3  << "\n"; }
            if(j==4 ) { b0radix4  << string_buf; b0radix4  << "\n"; }
            if(j==5 ) { b0radix5  << string_buf; b0radix5  << "\n"; }
            if(j==6 ) { b0radix6  << string_buf; b0radix6  << "\n"; }
            if(j==7 ) { b0radix7  << string_buf; b0radix7  << "\n"; }
            if(j==8 ) { b0radix8  << string_buf; b0radix8  << "\n"; }
            if(j==9 ) { b0radix9  << string_buf; b0radix9  << "\n"; }
            if(j==10) { b0radix10 << string_buf; b0radix10 << "\n"; }
            if(j==11) { b0radix11 << string_buf; b0radix11 << "\n"; }
            if(j==12) { b0radix12 << string_buf; b0radix12 << "\n"; }
            if(j==13) { b0radix13 << string_buf; b0radix13 << "\n"; }
            if(j==14) { b0radix14 << string_buf; b0radix14 << "\n"; }
            if(j==15) { b0radix15 << string_buf; b0radix15 << "\n"; }

          }           
           loop_index = 0;
           counter = counter + 1;
       }
       else loop_index = loop_index + 1;
    }
  
    counter      = 0;
    loop_index   = 0;
    loop_address = 0;
    //bank 1
    while(counter < counter_iteration){ 
        if((ma[loop_index] == counter )&& (bn[loop_index] == 1) ){
            for(int j=0;j<radix;j++){
            loop_address   = loop_index + offset*j;
            time_data_tmp  = time_data[loop_address];
            string_buf     = ZZtohex_cp(time_data_tmp);

            if(j==0 ) { b1radix0  << string_buf; b1radix0  << "\n"; }
            if(j==1 ) { b1radix1  << string_buf; b1radix1  << "\n"; }
            if(j==2 ) { b1radix2  << string_buf; b1radix2  << "\n"; }
            if(j==3 ) { b1radix3  << string_buf; b1radix3  << "\n"; }
            if(j==4 ) { b1radix4  << string_buf; b1radix4  << "\n"; }
            if(j==5 ) { b1radix5  << string_buf; b1radix5  << "\n"; }
            if(j==6 ) { b1radix6  << string_buf; b1radix6  << "\n"; }
            if(j==7 ) { b1radix7  << string_buf; b1radix7  << "\n"; }
            if(j==8 ) { b1radix8  << string_buf; b1radix8  << "\n"; }
            if(j==9 ) { b1radix9  << string_buf; b1radix9  << "\n"; }
            if(j==10) { b1radix10 << string_buf; b1radix10 << "\n"; }
            if(j==11) { b1radix11 << string_buf; b1radix11 << "\n"; }
            if(j==12) { b1radix12 << string_buf; b1radix12 << "\n"; }
            if(j==13) { b1radix13 << string_buf; b1radix13 << "\n"; }
            if(j==14) { b1radix14 << string_buf; b1radix14 << "\n"; }
            if(j==15) { b1radix15 << string_buf; b1radix15 << "\n"; }
            
            }           
             loop_index = 0;
             counter    = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
    
    //time_of.close();
    b0radix0.close();
	b0radix1.close();
	b0radix2.close();
	b0radix3.close();
	b0radix4.close();
	b0radix5.close();
	b0radix6.close();
	b0radix7.close();
	b0radix8.close();
	b0radix9.close();
	b0radix10.close();
	b0radix11.close();
	b0radix12.close();
	b0radix13.close();
	b0radix14.close();
	b0radix15.close();
	b1radix0.close();
	b1radix1.close();
	b1radix2.close();
	b1radix3.close();
	b1radix4.close();
	b1radix5.close();
	b1radix6.close();
	b1radix7.close();
	b1radix8.close();
	b1radix9.close();
	b1radix10.close();
	b1radix11.close();
	b1radix12.close();
	b1radix13.close();
	b1radix14.close();
	b1radix15.close();
	
}
void SPMB::H_freq_o_r16(std::vector<ZZ> H_NTT){
    //1 word  128 bits
	//bank 0
    std::ofstream H_freq_out("./ROM_Data/H_freq_out.txt");
    std::ofstream H_b0ROM0("./ROM_Data/H_b0ROM0.txt");
    std::ofstream H_b0ROM1("./ROM_Data/H_b0ROM1.txt");
    std::ofstream H_b0ROM2("./ROM_Data/H_b0ROM2.txt");
    std::ofstream H_b0ROM3("./ROM_Data/H_b0ROM3.txt");
    std::ofstream H_b0ROM4("./ROM_Data/H_b0ROM4.txt");
    std::ofstream H_b0ROM5("./ROM_Data/H_b0ROM5.txt");
    std::ofstream H_b0ROM6("./ROM_Data/H_b0ROM6.txt");
    std::ofstream H_b0ROM7("./ROM_Data/H_b0ROM7.txt");
    //bank1                     ROM_Data
    std::ofstream H_b1ROM0("./ROM_Data/H_b1ROM0.txt");
    std::ofstream H_b1ROM1("./ROM_Data/H_b1ROM1.txt");
    std::ofstream H_b1ROM2("./ROM_Data/H_b1ROM2.txt");
    std::ofstream H_b1ROM3("./ROM_Data/H_b1ROM3.txt");
    std::ofstream H_b1ROM4("./ROM_Data/H_b1ROM4.txt");
    std::ofstream H_b1ROM5("./ROM_Data/H_b1ROM5.txt");
    std::ofstream H_b1ROM6("./ROM_Data/H_b1ROM6.txt");
    std::ofstream H_b1ROM7("./ROM_Data/H_b1ROM7.txt");
    
    std::vector<int> tobit_r0;
    std::vector<int> tobit_r1;
    std::vector<int> tobit_r2;
    std::vector<int> tobit_r3;
    std::vector<int> tobit_r4;
    std::vector<int> tobit_r5;
    std::vector<int> tobit_r6;
    std::vector<int> tobit_r7;
    std::vector<int> tobit_r8;
    std::vector<int> tobit_r9;
    std::vector<int> tobit_r10;
    std::vector<int> tobit_r11;
    std::vector<int> tobit_r12;
    std::vector<int> tobit_r13;
    std::vector<int> tobit_r14;
    std::vector<int> tobit_r15;
    tobit_r0.resize(64);
    tobit_r1.resize(64);
    tobit_r2.resize(64);
    tobit_r3.resize(64);
    tobit_r4.resize(64);
    tobit_r5.resize(64);
    tobit_r6.resize(64);
    tobit_r7.resize(64);
    tobit_r8.resize(64);
    tobit_r9.resize(64);
    tobit_r10.resize(64);
    tobit_r11.resize(64);
    tobit_r12.resize(64);
    tobit_r13.resize(64);
    tobit_r14.resize(64);
    tobit_r15.resize(64);
    
    ZZ H_NTT_tmp_r0;
    ZZ H_NTT_tmp_r1;
    ZZ H_NTT_tmp_r2;
    ZZ H_NTT_tmp_r3;
    ZZ H_NTT_tmp_r4;
    ZZ H_NTT_tmp_r5;
    ZZ H_NTT_tmp_r6;
    ZZ H_NTT_tmp_r7;
    ZZ H_NTT_tmp_r8;
    ZZ H_NTT_tmp_r9;
    ZZ H_NTT_tmp_r10;
    ZZ H_NTT_tmp_r11;
    ZZ H_NTT_tmp_r12;
    ZZ H_NTT_tmp_r13;
    ZZ H_NTT_tmp_r14;
    ZZ H_NTT_tmp_r15;
    
    int counter =0;
    int loop_index = 0;
     
    for(int i=0;i<fft_point;i++){
        H_freq_out << H_NTT[i];
        H_freq_out << "\n";
    }
    
    //bank0
    while(counter < counter_iteration){ 
        if((ima[loop_index] == counter )&& (ibn[loop_index] == 0) ){
            
            H_NTT_tmp_r0  = H_NTT[loop_index              ];
            H_NTT_tmp_r1  = H_NTT[loop_index + 1  * offset ];
            H_NTT_tmp_r2  = H_NTT[loop_index + 2  * offset ];
            H_NTT_tmp_r3  = H_NTT[loop_index + 3  * offset ];
            H_NTT_tmp_r4  = H_NTT[loop_index + 4  * offset ];
            H_NTT_tmp_r5  = H_NTT[loop_index + 5  * offset ];
            H_NTT_tmp_r6  = H_NTT[loop_index + 6  * offset ];
            H_NTT_tmp_r7  = H_NTT[loop_index + 7  * offset ];
            H_NTT_tmp_r8  = H_NTT[loop_index + 8  * offset ];
            H_NTT_tmp_r9  = H_NTT[loop_index + 9  * offset ];
            H_NTT_tmp_r10 = H_NTT[loop_index + 10 * offset ];
            H_NTT_tmp_r11 = H_NTT[loop_index + 11 * offset ];
            H_NTT_tmp_r12 = H_NTT[loop_index + 12 * offset ];
            H_NTT_tmp_r13 = H_NTT[loop_index + 13 * offset ];
            H_NTT_tmp_r14 = H_NTT[loop_index + 14 * offset ];
            H_NTT_tmp_r15 = H_NTT[loop_index + 15 * offset ];
            
            for(int bit_index=0; bit_index < 64 ;bit_index++){
                if(H_NTT_tmp_r0 % 2 == 1) tobit_r0[bit_index] = 1;
                else tobit_r0[bit_index] = 0;
				if(H_NTT_tmp_r1 % 2 == 1) tobit_r1[bit_index] = 1;
				else tobit_r1[bit_index] = 0;
				if(H_NTT_tmp_r2 % 2 == 1) tobit_r2[bit_index] = 1;
				else tobit_r2[bit_index] = 0;
				if(H_NTT_tmp_r3 % 2 == 1) tobit_r3[bit_index] = 1;
				else tobit_r3[bit_index] = 0;
				if(H_NTT_tmp_r4 % 2 == 1) tobit_r4[bit_index] = 1;
				else tobit_r4[bit_index] = 0;
				if(H_NTT_tmp_r5 % 2 == 1) tobit_r5[bit_index] = 1;
				else tobit_r5[bit_index] = 0;
				if(H_NTT_tmp_r6 % 2 == 1) tobit_r6[bit_index] = 1;
				else tobit_r6[bit_index] = 0;
				if(H_NTT_tmp_r7 % 2 == 1) tobit_r7[bit_index] = 1;
				else tobit_r7[bit_index] = 0;
				if(H_NTT_tmp_r8 % 2 == 1) tobit_r8[bit_index] = 1;
				else tobit_r8[bit_index] = 0;
				if(H_NTT_tmp_r9 % 2 == 1) tobit_r9[bit_index] = 1;
				else tobit_r9[bit_index] = 0;
				if(H_NTT_tmp_r10 % 2 == 1) tobit_r10[bit_index] = 1;
				else tobit_r10[bit_index] = 0;
				if(H_NTT_tmp_r11 % 2 == 1) tobit_r11[bit_index] = 1;
				else tobit_r11[bit_index] = 0;
				if(H_NTT_tmp_r12 % 2 == 1) tobit_r12[bit_index] = 1;
				else tobit_r12[bit_index] = 0;
				if(H_NTT_tmp_r13 % 2 == 1) tobit_r13[bit_index] = 1;
				else tobit_r13[bit_index] = 0;
				if(H_NTT_tmp_r14 % 2 == 1) tobit_r14[bit_index] = 1;
				else tobit_r14[bit_index] = 0;
				if(H_NTT_tmp_r15 % 2 == 1) tobit_r15[bit_index] = 1;
				else tobit_r15[bit_index] = 0;
				
                H_NTT_tmp_r0  = H_NTT_tmp_r0 >> 1;
                H_NTT_tmp_r1  = H_NTT_tmp_r1 >> 1;
                H_NTT_tmp_r2  = H_NTT_tmp_r2 >> 1;
                H_NTT_tmp_r3  = H_NTT_tmp_r3 >> 1;
                H_NTT_tmp_r4  = H_NTT_tmp_r4 >> 1;
                H_NTT_tmp_r5  = H_NTT_tmp_r5 >> 1;
                H_NTT_tmp_r6  = H_NTT_tmp_r6 >> 1;
                H_NTT_tmp_r7  = H_NTT_tmp_r7 >> 1;
                H_NTT_tmp_r8  = H_NTT_tmp_r8 >> 1;
                H_NTT_tmp_r9  = H_NTT_tmp_r9 >> 1;
                H_NTT_tmp_r10 = H_NTT_tmp_r10 >> 1;
                H_NTT_tmp_r11 = H_NTT_tmp_r11 >> 1;
                H_NTT_tmp_r12 = H_NTT_tmp_r12 >> 1;
                H_NTT_tmp_r13 = H_NTT_tmp_r13 >> 1;
                H_NTT_tmp_r14 = H_NTT_tmp_r14 >> 1;
                H_NTT_tmp_r15 = H_NTT_tmp_r15 >> 1;
            }
            
            for(int g=0;g < 64; g++){
				H_b0ROM0 << tobit_r0[63-g];
			}
            for(int g=0;g < 64; g++){ 
				H_b0ROM0 << tobit_r1[63-g];
			} 
			H_b0ROM0  << "\n";
            
			for(int g=0;g < 64; g++){
				H_b0ROM1 << tobit_r2[63-g];
			} 
            for(int g=0;g < 64; g++){
            	H_b0ROM1 << tobit_r3[63-g];
            } 
            H_b0ROM1  << "\n";
			
			for(int g=0;g < 64; g++){
            	H_b0ROM2 << tobit_r4[63-g];
            } 
            for(int g=0;g < 64; g++){
            	H_b0ROM2 << tobit_r5[63-g];
            } 
			H_b0ROM2 << "\n";
			
            
			for(int g=0;g < 64; g++){
            	H_b0ROM3 << tobit_r6[63-g];
            } 
            for(int g=0;g < 64; g++){
				H_b0ROM3 << tobit_r7[63-g];
            } 
			H_b0ROM3 << "\n";

			for(int g=0;g < 64; g++){
            	H_b0ROM4 << tobit_r8[63-g];
            } 
			for(int g=0;g < 64; g++){
            	H_b0ROM4 << tobit_r9[63-g];
            } 
			H_b0ROM4 << "\n";
			
			for(int g=0;g < 64; g++){
				H_b0ROM5 << tobit_r10[63-g];
			} 
			for(int g=0;g < 64; g++){
				H_b0ROM5 << tobit_r11[63-g];
			} 
			H_b0ROM5 << "\n";
			
			for(int g=0;g < 64; g++){
				H_b0ROM6 << tobit_r12[63-g];
			} 
			for(int g=0;g < 64; g++){
				H_b0ROM6 << tobit_r13[63-g];
			} 
			H_b0ROM6 << "\n";
			
			for(int g=0;g < 64; g++){
				H_b0ROM7 << tobit_r14[63-g];
			} 
			for(int g=0;g < 64; g++){
				H_b0ROM7 << tobit_r15[63-g];
			} 
			H_b0ROM7 << "\n";
			
			
			
            loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
    
    counter = 0;
    loop_index = 0;
    //bank 1
    while(counter < counter_iteration){ 
        if((ima[loop_index] == counter )&& (ibn[loop_index] == 1) ){
            H_NTT_tmp_r0  = H_NTT[loop_index              ];
            H_NTT_tmp_r1  = H_NTT[loop_index + 1  * offset ];
            H_NTT_tmp_r2  = H_NTT[loop_index + 2  * offset ];
            H_NTT_tmp_r3  = H_NTT[loop_index + 3  * offset ];
            H_NTT_tmp_r4  = H_NTT[loop_index + 4  * offset ];
            H_NTT_tmp_r5  = H_NTT[loop_index + 5  * offset ];
            H_NTT_tmp_r6  = H_NTT[loop_index + 6  * offset ];
            H_NTT_tmp_r7  = H_NTT[loop_index + 7  * offset ];
            H_NTT_tmp_r8  = H_NTT[loop_index + 8  * offset ];
            H_NTT_tmp_r9  = H_NTT[loop_index + 9  * offset ];
            H_NTT_tmp_r10 = H_NTT[loop_index + 10 * offset ];
            H_NTT_tmp_r11 = H_NTT[loop_index + 11 * offset ];
            H_NTT_tmp_r12 = H_NTT[loop_index + 12 * offset ];
            H_NTT_tmp_r13 = H_NTT[loop_index + 13 * offset ];
            H_NTT_tmp_r14 = H_NTT[loop_index + 14 * offset ];
            H_NTT_tmp_r15 = H_NTT[loop_index + 15 * offset ];
            
            for(int bit_index=0; bit_index < 64 ;bit_index++){
                if(H_NTT_tmp_r0 % 2 == 1) tobit_r0[bit_index] = 1;
                else tobit_r0[bit_index] = 0;
				if(H_NTT_tmp_r1 % 2 == 1) tobit_r1[bit_index] = 1;
				else tobit_r1[bit_index] = 0;
				if(H_NTT_tmp_r2 % 2 == 1) tobit_r2[bit_index] = 1;
				else tobit_r2[bit_index] = 0;
				if(H_NTT_tmp_r3 % 2 == 1) tobit_r3[bit_index] = 1;
				else tobit_r3[bit_index] = 0;
				if(H_NTT_tmp_r4 % 2 == 1) tobit_r4[bit_index] = 1;
				else tobit_r4[bit_index] = 0;
				if(H_NTT_tmp_r5 % 2 == 1) tobit_r5[bit_index] = 1;
				else tobit_r5[bit_index] = 0;
				if(H_NTT_tmp_r6 % 2 == 1) tobit_r6[bit_index] = 1;
				else tobit_r6[bit_index] = 0;
				if(H_NTT_tmp_r7 % 2 == 1) tobit_r7[bit_index] = 1;
				else tobit_r7[bit_index] = 0;
				if(H_NTT_tmp_r8 % 2 == 1) tobit_r8[bit_index] = 1;
				else tobit_r8[bit_index] = 0;
				if(H_NTT_tmp_r9 % 2 == 1) tobit_r9[bit_index] = 1;
				else tobit_r9[bit_index] = 0;
				if(H_NTT_tmp_r10 % 2 == 1) tobit_r10[bit_index] = 1;
				else tobit_r10[bit_index] = 0;
				if(H_NTT_tmp_r11 % 2 == 1) tobit_r11[bit_index] = 1;
				else tobit_r11[bit_index] = 0;
				if(H_NTT_tmp_r12 % 2 == 1) tobit_r12[bit_index] = 1;
				else tobit_r12[bit_index] = 0;
				if(H_NTT_tmp_r13 % 2 == 1) tobit_r13[bit_index] = 1;
				else tobit_r13[bit_index] = 0;
				if(H_NTT_tmp_r14 % 2 == 1) tobit_r14[bit_index] = 1;
				else tobit_r14[bit_index] = 0;
				if(H_NTT_tmp_r15 % 2 == 1) tobit_r15[bit_index] = 1;
				else tobit_r15[bit_index] = 0;
				
                H_NTT_tmp_r0  = H_NTT_tmp_r0 >> 1;
                H_NTT_tmp_r1  = H_NTT_tmp_r1 >> 1;
                H_NTT_tmp_r2  = H_NTT_tmp_r2 >> 1;
                H_NTT_tmp_r3  = H_NTT_tmp_r3 >> 1;
                H_NTT_tmp_r4  = H_NTT_tmp_r4 >> 1;
                H_NTT_tmp_r5  = H_NTT_tmp_r5 >> 1;
                H_NTT_tmp_r6  = H_NTT_tmp_r6 >> 1;
                H_NTT_tmp_r7  = H_NTT_tmp_r7 >> 1;
                H_NTT_tmp_r8  = H_NTT_tmp_r8 >> 1;
                H_NTT_tmp_r9  = H_NTT_tmp_r9 >> 1;
                H_NTT_tmp_r10 = H_NTT_tmp_r10 >> 1;
                H_NTT_tmp_r11 = H_NTT_tmp_r11 >> 1;
                H_NTT_tmp_r12 = H_NTT_tmp_r12 >> 1;
                H_NTT_tmp_r13 = H_NTT_tmp_r13 >> 1;
                H_NTT_tmp_r14 = H_NTT_tmp_r14 >> 1;
                H_NTT_tmp_r15 = H_NTT_tmp_r15 >> 1;
            }
            
            for(int g=0;g < 64; g++){
				H_b1ROM0 << tobit_r0[63-g];
			}
            for(int g=0;g < 64; g++){ 
				H_b1ROM0 << tobit_r1[63-g];
			} 
			H_b1ROM0  << "\n";
            
			for(int g=0;g < 64; g++){
				H_b1ROM1 << tobit_r2[63-g];
			} 
            for(int g=0;g < 64; g++){
            	H_b1ROM1 << tobit_r3[63-g];
            } 
            H_b1ROM1  << "\n";
			
			for(int g=0;g < 64; g++){
            	H_b1ROM2 << tobit_r4[63-g];
            } 
            for(int g=0;g < 64; g++){
            	H_b1ROM2 << tobit_r5[63-g];
            } 
			H_b1ROM2 << "\n";
			
            
			for(int g=0;g < 64; g++){
            	H_b1ROM3 << tobit_r6[63-g];
            } 
            for(int g=0;g < 64; g++){
				H_b1ROM3 << tobit_r7[63-g];
            } 
			H_b1ROM3 << "\n";

			for(int g=0;g < 64; g++){
            	H_b1ROM4 << tobit_r8[63-g];
            } 
			for(int g=0;g < 64; g++){
            	H_b1ROM4 << tobit_r9[63-g];
            } 
			H_b1ROM4 << "\n";
			
			for(int g=0;g < 64; g++){
				H_b1ROM5 << tobit_r10[63-g];
			} 
			for(int g=0;g < 64; g++){
				H_b1ROM5 << tobit_r11[63-g];
			} 
			H_b1ROM5 << "\n";
			
			for(int g=0;g < 64; g++){
				H_b1ROM6 << tobit_r12[63-g];
			} 
			for(int g=0;g < 64; g++){
				H_b1ROM6 << tobit_r13[63-g];
			} 
			H_b1ROM6 << "\n";
			
			for(int g=0;g < 64; g++){
				H_b1ROM7 << tobit_r14[63-g];
			} 
			for(int g=0;g < 64; g++){
				H_b1ROM7 << tobit_r15[63-g];
			} 
			H_b1ROM7 << "\n";
            
            loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
}
void SPMB::H_freq_o_r16_r2(std::vector<ZZ> H_NTT){
    //1 word  128 bits
	//bank 0
    std::ofstream H_freq_out("./ROM_Data/H_freq_out.txt");
    std::ofstream H_b0ROM0("./ROM_Data/H_b0ROM0.txt");
    std::ofstream H_b0ROM1("./ROM_Data/H_b0ROM1.txt");
    std::ofstream H_b0ROM2("./ROM_Data/H_b0ROM2.txt");
    std::ofstream H_b0ROM3("./ROM_Data/H_b0ROM3.txt");
    std::ofstream H_b0ROM4("./ROM_Data/H_b0ROM4.txt");
    std::ofstream H_b0ROM5("./ROM_Data/H_b0ROM5.txt");
    std::ofstream H_b0ROM6("./ROM_Data/H_b0ROM6.txt");
    std::ofstream H_b0ROM7("./ROM_Data/H_b0ROM7.txt");
    //bank1                     ROM_Data
    std::ofstream H_b1ROM0("./ROM_Data/H_b1ROM0.txt");
    std::ofstream H_b1ROM1("./ROM_Data/H_b1ROM1.txt");
    std::ofstream H_b1ROM2("./ROM_Data/H_b1ROM2.txt");
    std::ofstream H_b1ROM3("./ROM_Data/H_b1ROM3.txt");
    std::ofstream H_b1ROM4("./ROM_Data/H_b1ROM4.txt");
    std::ofstream H_b1ROM5("./ROM_Data/H_b1ROM5.txt");
    std::ofstream H_b1ROM6("./ROM_Data/H_b1ROM6.txt");
    std::ofstream H_b1ROM7("./ROM_Data/H_b1ROM7.txt");
	
	//*************************************
	//H_freq data output
	for(int i=0;i<fft_point;i++){
        H_freq_out << H_NTT[i];
        H_freq_out << "\n";
    }
	//*************************************
    
    std::vector<int> tobit_r0;
    std::vector<int> tobit_r1;
    std::vector<int> tobit_r2;
    std::vector<int> tobit_r3;
    std::vector<int> tobit_r4;
    std::vector<int> tobit_r5;
    std::vector<int> tobit_r6;
    std::vector<int> tobit_r7;
    std::vector<int> tobit_r8;
    std::vector<int> tobit_r9;
    std::vector<int> tobit_r10;
    std::vector<int> tobit_r11;
    std::vector<int> tobit_r12;
    std::vector<int> tobit_r13;
    std::vector<int> tobit_r14;
    std::vector<int> tobit_r15;
    tobit_r0.resize(64);
    tobit_r1.resize(64);
    tobit_r2.resize(64);
    tobit_r3.resize(64);
    tobit_r4.resize(64);
    tobit_r5.resize(64);
    tobit_r6.resize(64);
    tobit_r7.resize(64);
    tobit_r8.resize(64);
    tobit_r9.resize(64);
    tobit_r10.resize(64);
    tobit_r11.resize(64);
    tobit_r12.resize(64);
    tobit_r13.resize(64);
    tobit_r14.resize(64);
    tobit_r15.resize(64);
	
    ZZ H_NTT_tmp_r0;
    ZZ H_NTT_tmp_r1;
    ZZ H_NTT_tmp_r2;
    ZZ H_NTT_tmp_r3;
    ZZ H_NTT_tmp_r4;
    ZZ H_NTT_tmp_r5;
    ZZ H_NTT_tmp_r6;
    ZZ H_NTT_tmp_r7;
    ZZ H_NTT_tmp_r8;
    ZZ H_NTT_tmp_r9;
    ZZ H_NTT_tmp_r10;
    ZZ H_NTT_tmp_r11;
    ZZ H_NTT_tmp_r12;
    ZZ H_NTT_tmp_r13;
    ZZ H_NTT_tmp_r14;
    ZZ H_NTT_tmp_r15;
    
    int counter =0;
    int loop_index = 0;
	int INTT_R1_offset;
	int INTT_R2_offset;
	int INTT_R3_offset;
	int INTT_R4_offset;
	int INTT_R5_offset;
	int INTT_R6_offset;
	int INTT_R7_offset;
	int INTT_R8_offset;
	int INTT_R9_offset;
	int INTT_R10_offset;
	int INTT_R11_offset;
	int INTT_R12_offset;
	int INTT_R13_offset;
	int INTT_R14_offset;
	int INTT_R15_offset;
     
	 
	//calculate INTT offset index 
	INTT_R1_offset  = (int) fft_point / 2;
	INTT_R2_offset  = (int) fft_point / 16;
	INTT_R3_offset  = INTT_R1_offset + INTT_R2_offset;
	INTT_R4_offset  = (int)fft_point / 8;
	INTT_R5_offset  = INTT_R4_offset + INTT_R1_offset;
	INTT_R6_offset  = INTT_R4_offset + INTT_R2_offset;
	INTT_R7_offset  = INTT_R4_offset + INTT_R3_offset;
	INTT_R8_offset  = (int)fft_point / 4;
	INTT_R9_offset  = INTT_R8_offset + INTT_R1_offset ;
	INTT_R10_offset = INTT_R8_offset + INTT_R2_offset ;
	INTT_R11_offset = INTT_R8_offset + INTT_R3_offset ;
	INTT_R12_offset = INTT_R8_offset + INTT_R4_offset ;
	INTT_R13_offset = INTT_R8_offset + INTT_R5_offset ;
	INTT_R14_offset = INTT_R8_offset + INTT_R6_offset ;
	INTT_R15_offset = INTT_R8_offset + INTT_R7_offset ;
    
    //bank0
    while(counter < counter_iteration){ 
        if((ima[loop_index] == counter )&& (ibn[loop_index] == 0) ){
			
            H_NTT_tmp_r0  = H_NTT[loop_index ];
            H_NTT_tmp_r1  = H_NTT[loop_index + INTT_R1_offset ];
            H_NTT_tmp_r2  = H_NTT[loop_index + INTT_R2_offset ];
            H_NTT_tmp_r3  = H_NTT[loop_index + INTT_R3_offset ];
            H_NTT_tmp_r4  = H_NTT[loop_index + INTT_R4_offset ];
            H_NTT_tmp_r5  = H_NTT[loop_index + INTT_R5_offset ];
            H_NTT_tmp_r6  = H_NTT[loop_index + INTT_R6_offset ];
            H_NTT_tmp_r7  = H_NTT[loop_index + INTT_R7_offset ];
            H_NTT_tmp_r8  = H_NTT[loop_index + INTT_R8_offset ];
            H_NTT_tmp_r9  = H_NTT[loop_index + INTT_R9_offset ];
            H_NTT_tmp_r10 = H_NTT[loop_index + INTT_R10_offset ];
            H_NTT_tmp_r11 = H_NTT[loop_index + INTT_R11_offset ];
            H_NTT_tmp_r12 = H_NTT[loop_index + INTT_R12_offset ];
            H_NTT_tmp_r13 = H_NTT[loop_index + INTT_R13_offset ];
            H_NTT_tmp_r14 = H_NTT[loop_index + INTT_R14_offset ];
            H_NTT_tmp_r15 = H_NTT[loop_index + INTT_R15_offset ];
            
            for(int bit_index=0; bit_index < 64 ;bit_index++){
                if(H_NTT_tmp_r0 % 2 == 1) tobit_r0[bit_index] = 1;
                else tobit_r0[bit_index] = 0;
				if(H_NTT_tmp_r1 % 2 == 1) tobit_r1[bit_index] = 1;
				else tobit_r1[bit_index] = 0;
				if(H_NTT_tmp_r2 % 2 == 1) tobit_r2[bit_index] = 1;
				else tobit_r2[bit_index] = 0;
				if(H_NTT_tmp_r3 % 2 == 1) tobit_r3[bit_index] = 1;
				else tobit_r3[bit_index] = 0;
				if(H_NTT_tmp_r4 % 2 == 1) tobit_r4[bit_index] = 1;
				else tobit_r4[bit_index] = 0;
				if(H_NTT_tmp_r5 % 2 == 1) tobit_r5[bit_index] = 1;
				else tobit_r5[bit_index] = 0;
				if(H_NTT_tmp_r6 % 2 == 1) tobit_r6[bit_index] = 1;
				else tobit_r6[bit_index] = 0;
				if(H_NTT_tmp_r7 % 2 == 1) tobit_r7[bit_index] = 1;
				else tobit_r7[bit_index] = 0;
				if(H_NTT_tmp_r8 % 2 == 1) tobit_r8[bit_index] = 1;
				else tobit_r8[bit_index] = 0;
				if(H_NTT_tmp_r9 % 2 == 1) tobit_r9[bit_index] = 1;
				else tobit_r9[bit_index] = 0;
				if(H_NTT_tmp_r10 % 2 == 1) tobit_r10[bit_index] = 1;
				else tobit_r10[bit_index] = 0;
				if(H_NTT_tmp_r11 % 2 == 1) tobit_r11[bit_index] = 1;
				else tobit_r11[bit_index] = 0;
				if(H_NTT_tmp_r12 % 2 == 1) tobit_r12[bit_index] = 1;
				else tobit_r12[bit_index] = 0;
				if(H_NTT_tmp_r13 % 2 == 1) tobit_r13[bit_index] = 1;
				else tobit_r13[bit_index] = 0;
				if(H_NTT_tmp_r14 % 2 == 1) tobit_r14[bit_index] = 1;
				else tobit_r14[bit_index] = 0;
				if(H_NTT_tmp_r15 % 2 == 1) tobit_r15[bit_index] = 1;
				else tobit_r15[bit_index] = 0;
				
                H_NTT_tmp_r0  = H_NTT_tmp_r0 >> 1;
                H_NTT_tmp_r1  = H_NTT_tmp_r1 >> 1;
                H_NTT_tmp_r2  = H_NTT_tmp_r2 >> 1;
                H_NTT_tmp_r3  = H_NTT_tmp_r3 >> 1;
                H_NTT_tmp_r4  = H_NTT_tmp_r4 >> 1;
                H_NTT_tmp_r5  = H_NTT_tmp_r5 >> 1;
                H_NTT_tmp_r6  = H_NTT_tmp_r6 >> 1;
                H_NTT_tmp_r7  = H_NTT_tmp_r7 >> 1;
                H_NTT_tmp_r8  = H_NTT_tmp_r8 >> 1;
                H_NTT_tmp_r9  = H_NTT_tmp_r9 >> 1;
                H_NTT_tmp_r10 = H_NTT_tmp_r10 >> 1;
                H_NTT_tmp_r11 = H_NTT_tmp_r11 >> 1;
                H_NTT_tmp_r12 = H_NTT_tmp_r12 >> 1;
                H_NTT_tmp_r13 = H_NTT_tmp_r13 >> 1;
                H_NTT_tmp_r14 = H_NTT_tmp_r14 >> 1;
                H_NTT_tmp_r15 = H_NTT_tmp_r15 >> 1;
            }
            
            for(int g=0;g < 64; g++){
				H_b0ROM0 << tobit_r0[63-g];
			}
            for(int g=0;g < 64; g++){ 
				H_b0ROM0 << tobit_r1[63-g];
			} 
			H_b0ROM0  << "\n";
            
			for(int g=0;g < 64; g++){
				H_b0ROM1 << tobit_r2[63-g];
			} 
            for(int g=0;g < 64; g++){
            	H_b0ROM1 << tobit_r3[63-g];
            } 
            H_b0ROM1  << "\n";
			
			for(int g=0;g < 64; g++){
            	H_b0ROM2 << tobit_r4[63-g];
            } 
            for(int g=0;g < 64; g++){
            	H_b0ROM2 << tobit_r5[63-g];
            } 
			H_b0ROM2 << "\n";
			
            
			for(int g=0;g < 64; g++){
            	H_b0ROM3 << tobit_r6[63-g];
            } 
            for(int g=0;g < 64; g++){
				H_b0ROM3 << tobit_r7[63-g];
            } 
			H_b0ROM3 << "\n";

			for(int g=0;g < 64; g++){
            	H_b0ROM4 << tobit_r8[63-g];
            } 
			for(int g=0;g < 64; g++){
            	H_b0ROM4 << tobit_r9[63-g];
            } 
			H_b0ROM4 << "\n";
			
			for(int g=0;g < 64; g++){
				H_b0ROM5 << tobit_r10[63-g];
			} 
			for(int g=0;g < 64; g++){
				H_b0ROM5 << tobit_r11[63-g];
			} 
			H_b0ROM5 << "\n";
			
			for(int g=0;g < 64; g++){
				H_b0ROM6 << tobit_r12[63-g];
			} 
			for(int g=0;g < 64; g++){
				H_b0ROM6 << tobit_r13[63-g];
			} 
			H_b0ROM6 << "\n";
			
			for(int g=0;g < 64; g++){
				H_b0ROM7 << tobit_r14[63-g];
			} 
			for(int g=0;g < 64; g++){
				H_b0ROM7 << tobit_r15[63-g];
			} 
			H_b0ROM7 << "\n";
			
            loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
    
    counter = 0;
    loop_index = 0;
    //bank 1
    while(counter < counter_iteration){ 
        if((ima[loop_index] == counter )&& (ibn[loop_index] == 1) ){
            H_NTT_tmp_r0  = H_NTT[loop_index              ];
            H_NTT_tmp_r1  = H_NTT[loop_index + INTT_R1_offset  ];
            H_NTT_tmp_r2  = H_NTT[loop_index + INTT_R2_offset  ];
            H_NTT_tmp_r3  = H_NTT[loop_index + INTT_R3_offset  ];
            H_NTT_tmp_r4  = H_NTT[loop_index + INTT_R4_offset  ];
            H_NTT_tmp_r5  = H_NTT[loop_index + INTT_R5_offset  ];
            H_NTT_tmp_r6  = H_NTT[loop_index + INTT_R6_offset  ];
            H_NTT_tmp_r7  = H_NTT[loop_index + INTT_R7_offset  ];
            H_NTT_tmp_r8  = H_NTT[loop_index + INTT_R8_offset  ];
            H_NTT_tmp_r9  = H_NTT[loop_index + INTT_R9_offset  ];
            H_NTT_tmp_r10 = H_NTT[loop_index + INTT_R10_offset ];
            H_NTT_tmp_r11 = H_NTT[loop_index + INTT_R11_offset ];
            H_NTT_tmp_r12 = H_NTT[loop_index + INTT_R12_offset ];
            H_NTT_tmp_r13 = H_NTT[loop_index + INTT_R13_offset ];
            H_NTT_tmp_r14 = H_NTT[loop_index + INTT_R14_offset ];
            H_NTT_tmp_r15 = H_NTT[loop_index + INTT_R15_offset ];
            
            for(int bit_index=0; bit_index < 64 ;bit_index++){
                if(H_NTT_tmp_r0 % 2 == 1) tobit_r0[bit_index] = 1;
                else tobit_r0[bit_index] = 0;
				if(H_NTT_tmp_r1 % 2 == 1) tobit_r1[bit_index] = 1;
				else tobit_r1[bit_index] = 0;
				if(H_NTT_tmp_r2 % 2 == 1) tobit_r2[bit_index] = 1;
				else tobit_r2[bit_index] = 0;
				if(H_NTT_tmp_r3 % 2 == 1) tobit_r3[bit_index] = 1;
				else tobit_r3[bit_index] = 0;
				if(H_NTT_tmp_r4 % 2 == 1) tobit_r4[bit_index] = 1;
				else tobit_r4[bit_index] = 0;
				if(H_NTT_tmp_r5 % 2 == 1) tobit_r5[bit_index] = 1;
				else tobit_r5[bit_index] = 0;
				if(H_NTT_tmp_r6 % 2 == 1) tobit_r6[bit_index] = 1;
				else tobit_r6[bit_index] = 0;
				if(H_NTT_tmp_r7 % 2 == 1) tobit_r7[bit_index] = 1;
				else tobit_r7[bit_index] = 0;
				if(H_NTT_tmp_r8 % 2 == 1) tobit_r8[bit_index] = 1;
				else tobit_r8[bit_index] = 0;
				if(H_NTT_tmp_r9 % 2 == 1) tobit_r9[bit_index] = 1;
				else tobit_r9[bit_index] = 0;
				if(H_NTT_tmp_r10 % 2 == 1) tobit_r10[bit_index] = 1;
				else tobit_r10[bit_index] = 0;
				if(H_NTT_tmp_r11 % 2 == 1) tobit_r11[bit_index] = 1;
				else tobit_r11[bit_index] = 0;
				if(H_NTT_tmp_r12 % 2 == 1) tobit_r12[bit_index] = 1;
				else tobit_r12[bit_index] = 0;
				if(H_NTT_tmp_r13 % 2 == 1) tobit_r13[bit_index] = 1;
				else tobit_r13[bit_index] = 0;
				if(H_NTT_tmp_r14 % 2 == 1) tobit_r14[bit_index] = 1;
				else tobit_r14[bit_index] = 0;
				if(H_NTT_tmp_r15 % 2 == 1) tobit_r15[bit_index] = 1;
				else tobit_r15[bit_index] = 0;
				
                H_NTT_tmp_r0  = H_NTT_tmp_r0 >> 1;
                H_NTT_tmp_r1  = H_NTT_tmp_r1 >> 1;
                H_NTT_tmp_r2  = H_NTT_tmp_r2 >> 1;
                H_NTT_tmp_r3  = H_NTT_tmp_r3 >> 1;
                H_NTT_tmp_r4  = H_NTT_tmp_r4 >> 1;
                H_NTT_tmp_r5  = H_NTT_tmp_r5 >> 1;
                H_NTT_tmp_r6  = H_NTT_tmp_r6 >> 1;
                H_NTT_tmp_r7  = H_NTT_tmp_r7 >> 1;
                H_NTT_tmp_r8  = H_NTT_tmp_r8 >> 1;
                H_NTT_tmp_r9  = H_NTT_tmp_r9 >> 1;
                H_NTT_tmp_r10 = H_NTT_tmp_r10 >> 1;
                H_NTT_tmp_r11 = H_NTT_tmp_r11 >> 1;
                H_NTT_tmp_r12 = H_NTT_tmp_r12 >> 1;
                H_NTT_tmp_r13 = H_NTT_tmp_r13 >> 1;
                H_NTT_tmp_r14 = H_NTT_tmp_r14 >> 1;
                H_NTT_tmp_r15 = H_NTT_tmp_r15 >> 1;
            }
            
            for(int g=0;g < 64; g++){
				H_b1ROM0 << tobit_r0[63-g];
			}
            for(int g=0;g < 64; g++){ 
				H_b1ROM0 << tobit_r1[63-g];
			} 
			H_b1ROM0  << "\n";
            
			for(int g=0;g < 64; g++){
				H_b1ROM1 << tobit_r2[63-g];
			} 
            for(int g=0;g < 64; g++){
            	H_b1ROM1 << tobit_r3[63-g];
            } 
            H_b1ROM1  << "\n";
			
			for(int g=0;g < 64; g++){
            	H_b1ROM2 << tobit_r4[63-g];
            } 
            for(int g=0;g < 64; g++){
            	H_b1ROM2 << tobit_r5[63-g];
            } 
			H_b1ROM2 << "\n";
			
            
			for(int g=0;g < 64; g++){
            	H_b1ROM3 << tobit_r6[63-g];
            } 
            for(int g=0;g < 64; g++){
				H_b1ROM3 << tobit_r7[63-g];
            } 
			H_b1ROM3 << "\n";

			for(int g=0;g < 64; g++){
            	H_b1ROM4 << tobit_r8[63-g];
            } 
			for(int g=0;g < 64; g++){
            	H_b1ROM4 << tobit_r9[63-g];
            } 
			H_b1ROM4 << "\n";
			
			for(int g=0;g < 64; g++){
				H_b1ROM5 << tobit_r10[63-g];
			} 
			for(int g=0;g < 64; g++){
				H_b1ROM5 << tobit_r11[63-g];
			} 
			H_b1ROM5 << "\n";
			
			for(int g=0;g < 64; g++){
				H_b1ROM6 << tobit_r12[63-g];
			} 
			for(int g=0;g < 64; g++){
				H_b1ROM6 << tobit_r13[63-g];
			} 
			H_b1ROM6 << "\n";
			
			for(int g=0;g < 64; g++){
				H_b1ROM7 << tobit_r14[63-g];
			} 
			for(int g=0;g < 64; g++){
				H_b1ROM7 << tobit_r15[63-g];
			} 
			H_b1ROM7 << "\n";
            
            loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
}
void SPMB::H_freq_o_r16_r4(std::vector<ZZ> H_NTT){
    //1 word  128 bits
	//bank 0
    std::ofstream H_freq_out("./ROM_Data/H_freq_out.txt");
    std::ofstream H_b0ROM0("./ROM_Data/H_b0ROM0.txt");
    std::ofstream H_b0ROM1("./ROM_Data/H_b0ROM1.txt");
    std::ofstream H_b0ROM2("./ROM_Data/H_b0ROM2.txt");
    std::ofstream H_b0ROM3("./ROM_Data/H_b0ROM3.txt");
    std::ofstream H_b0ROM4("./ROM_Data/H_b0ROM4.txt");
    std::ofstream H_b0ROM5("./ROM_Data/H_b0ROM5.txt");
    std::ofstream H_b0ROM6("./ROM_Data/H_b0ROM6.txt");
    std::ofstream H_b0ROM7("./ROM_Data/H_b0ROM7.txt");
    //bank1                     ROM_Data
    std::ofstream H_b1ROM0("./ROM_Data/H_b1ROM0.txt");
    std::ofstream H_b1ROM1("./ROM_Data/H_b1ROM1.txt");
    std::ofstream H_b1ROM2("./ROM_Data/H_b1ROM2.txt");
    std::ofstream H_b1ROM3("./ROM_Data/H_b1ROM3.txt");
    std::ofstream H_b1ROM4("./ROM_Data/H_b1ROM4.txt");
    std::ofstream H_b1ROM5("./ROM_Data/H_b1ROM5.txt");
    std::ofstream H_b1ROM6("./ROM_Data/H_b1ROM6.txt");
    std::ofstream H_b1ROM7("./ROM_Data/H_b1ROM7.txt");
    
    for(int i=0;i<fft_point;i++){
        H_freq_out << H_NTT[i];
        H_freq_out << "\n";
    }	
	
    std::vector<int> tobit_r0;
    std::vector<int> tobit_r1;
    std::vector<int> tobit_r2;
    std::vector<int> tobit_r3;
    std::vector<int> tobit_r4;
    std::vector<int> tobit_r5;
    std::vector<int> tobit_r6;
    std::vector<int> tobit_r7;
    std::vector<int> tobit_r8;
    std::vector<int> tobit_r9;
    std::vector<int> tobit_r10;
    std::vector<int> tobit_r11;
    std::vector<int> tobit_r12;
    std::vector<int> tobit_r13;
    std::vector<int> tobit_r14;
    std::vector<int> tobit_r15;
    tobit_r0.resize(64);
    tobit_r1.resize(64);
    tobit_r2.resize(64);
    tobit_r3.resize(64);
    tobit_r4.resize(64);
    tobit_r5.resize(64);
    tobit_r6.resize(64);
    tobit_r7.resize(64);
    tobit_r8.resize(64);
    tobit_r9.resize(64);
    tobit_r10.resize(64);
    tobit_r11.resize(64);
    tobit_r12.resize(64);
    tobit_r13.resize(64);
    tobit_r14.resize(64);
    tobit_r15.resize(64);
    
    ZZ H_NTT_tmp_r0;
    ZZ H_NTT_tmp_r1;
    ZZ H_NTT_tmp_r2;
    ZZ H_NTT_tmp_r3;
    ZZ H_NTT_tmp_r4;
    ZZ H_NTT_tmp_r5;
    ZZ H_NTT_tmp_r6;
    ZZ H_NTT_tmp_r7;
    ZZ H_NTT_tmp_r8;
    ZZ H_NTT_tmp_r9;
    ZZ H_NTT_tmp_r10;
    ZZ H_NTT_tmp_r11;
    ZZ H_NTT_tmp_r12;
    ZZ H_NTT_tmp_r13;
    ZZ H_NTT_tmp_r14;
    ZZ H_NTT_tmp_r15;
    
    int counter =0;
    int loop_index = 0;
	int INTT_R1_offset;
	int INTT_R2_offset;
	int INTT_R3_offset;
	int INTT_R4_offset;
	int INTT_R5_offset;
	int INTT_R6_offset;
	int INTT_R7_offset;
	int INTT_R8_offset;
	int INTT_R9_offset;
	int INTT_R10_offset;
	int INTT_R11_offset;
	int INTT_R12_offset;
	int INTT_R13_offset;
	int INTT_R14_offset;
	int INTT_R15_offset;     

	//calculate INTT offset index 
	INTT_R1_offset  = (int) fft_point / 4;
	INTT_R2_offset  = (int) fft_point / 2;
	INTT_R3_offset  = INTT_R1_offset + INTT_R2_offset;
	INTT_R4_offset  = (int)fft_point / 16;
	INTT_R5_offset  = INTT_R4_offset + INTT_R1_offset;
	INTT_R6_offset  = INTT_R4_offset + INTT_R2_offset;
	INTT_R7_offset  = INTT_R4_offset + INTT_R3_offset;
	INTT_R8_offset  = (int)fft_point / 8;
	INTT_R9_offset  = INTT_R8_offset + INTT_R1_offset ;
	INTT_R10_offset = INTT_R8_offset + INTT_R2_offset ;
	INTT_R11_offset = INTT_R8_offset + INTT_R3_offset ;
	INTT_R12_offset = INTT_R8_offset + INTT_R4_offset ;
	INTT_R13_offset = INTT_R8_offset + INTT_R5_offset ;
	INTT_R14_offset = INTT_R8_offset + INTT_R6_offset ;
	INTT_R15_offset = INTT_R8_offset + INTT_R7_offset ;

    
    //bank0
    while(counter < counter_iteration){ 
        if((ima[loop_index] == counter )&& (ibn[loop_index] == 0) ){
            
            H_NTT_tmp_r0  = H_NTT[loop_index              ];
            H_NTT_tmp_r1  = H_NTT[loop_index + INTT_R1_offset ];
            H_NTT_tmp_r2  = H_NTT[loop_index + INTT_R2_offset ];
            H_NTT_tmp_r3  = H_NTT[loop_index + INTT_R3_offset ];
            H_NTT_tmp_r4  = H_NTT[loop_index + INTT_R4_offset ];
            H_NTT_tmp_r5  = H_NTT[loop_index + INTT_R5_offset ];
            H_NTT_tmp_r6  = H_NTT[loop_index + INTT_R6_offset ];
            H_NTT_tmp_r7  = H_NTT[loop_index + INTT_R7_offset ];
            H_NTT_tmp_r8  = H_NTT[loop_index + INTT_R8_offset ];
            H_NTT_tmp_r9  = H_NTT[loop_index + INTT_R9_offset ];
            H_NTT_tmp_r10 = H_NTT[loop_index + INTT_R10_offset];
            H_NTT_tmp_r11 = H_NTT[loop_index + INTT_R11_offset];
            H_NTT_tmp_r12 = H_NTT[loop_index + INTT_R12_offset];
            H_NTT_tmp_r13 = H_NTT[loop_index + INTT_R13_offset];
            H_NTT_tmp_r14 = H_NTT[loop_index + INTT_R14_offset];
            H_NTT_tmp_r15 = H_NTT[loop_index + INTT_R15_offset];
            
            for(int bit_index=0; bit_index < 64 ;bit_index++){
                if(H_NTT_tmp_r0 % 2 == 1) tobit_r0[bit_index] = 1;
                else tobit_r0[bit_index] = 0;
				if(H_NTT_tmp_r1 % 2 == 1) tobit_r1[bit_index] = 1;
				else tobit_r1[bit_index] = 0;
				if(H_NTT_tmp_r2 % 2 == 1) tobit_r2[bit_index] = 1;
				else tobit_r2[bit_index] = 0;
				if(H_NTT_tmp_r3 % 2 == 1) tobit_r3[bit_index] = 1;
				else tobit_r3[bit_index] = 0;
				if(H_NTT_tmp_r4 % 2 == 1) tobit_r4[bit_index] = 1;
				else tobit_r4[bit_index] = 0;
				if(H_NTT_tmp_r5 % 2 == 1) tobit_r5[bit_index] = 1;
				else tobit_r5[bit_index] = 0;
				if(H_NTT_tmp_r6 % 2 == 1) tobit_r6[bit_index] = 1;
				else tobit_r6[bit_index] = 0;
				if(H_NTT_tmp_r7 % 2 == 1) tobit_r7[bit_index] = 1;
				else tobit_r7[bit_index] = 0;
				if(H_NTT_tmp_r8 % 2 == 1) tobit_r8[bit_index] = 1;
				else tobit_r8[bit_index] = 0;
				if(H_NTT_tmp_r9 % 2 == 1) tobit_r9[bit_index] = 1;
				else tobit_r9[bit_index] = 0;
				if(H_NTT_tmp_r10 % 2 == 1) tobit_r10[bit_index] = 1;
				else tobit_r10[bit_index] = 0;
				if(H_NTT_tmp_r11 % 2 == 1) tobit_r11[bit_index] = 1;
				else tobit_r11[bit_index] = 0;
				if(H_NTT_tmp_r12 % 2 == 1) tobit_r12[bit_index] = 1;
				else tobit_r12[bit_index] = 0;
				if(H_NTT_tmp_r13 % 2 == 1) tobit_r13[bit_index] = 1;
				else tobit_r13[bit_index] = 0;
				if(H_NTT_tmp_r14 % 2 == 1) tobit_r14[bit_index] = 1;
				else tobit_r14[bit_index] = 0;
				if(H_NTT_tmp_r15 % 2 == 1) tobit_r15[bit_index] = 1;
				else tobit_r15[bit_index] = 0;
				
                H_NTT_tmp_r0  = H_NTT_tmp_r0 >> 1;
                H_NTT_tmp_r1  = H_NTT_tmp_r1 >> 1;
                H_NTT_tmp_r2  = H_NTT_tmp_r2 >> 1;
                H_NTT_tmp_r3  = H_NTT_tmp_r3 >> 1;
                H_NTT_tmp_r4  = H_NTT_tmp_r4 >> 1;
                H_NTT_tmp_r5  = H_NTT_tmp_r5 >> 1;
                H_NTT_tmp_r6  = H_NTT_tmp_r6 >> 1;
                H_NTT_tmp_r7  = H_NTT_tmp_r7 >> 1;
                H_NTT_tmp_r8  = H_NTT_tmp_r8 >> 1;
                H_NTT_tmp_r9  = H_NTT_tmp_r9 >> 1;
                H_NTT_tmp_r10 = H_NTT_tmp_r10 >> 1;
                H_NTT_tmp_r11 = H_NTT_tmp_r11 >> 1;
                H_NTT_tmp_r12 = H_NTT_tmp_r12 >> 1;
                H_NTT_tmp_r13 = H_NTT_tmp_r13 >> 1;
                H_NTT_tmp_r14 = H_NTT_tmp_r14 >> 1;
                H_NTT_tmp_r15 = H_NTT_tmp_r15 >> 1;
            }
            
            for(int g=0;g < 64; g++){
				H_b0ROM0 << tobit_r0[63-g];
			}
            for(int g=0;g < 64; g++){ 
				H_b0ROM0 << tobit_r1[63-g];
			} 
			H_b0ROM0  << "\n";
            
			for(int g=0;g < 64; g++){
				H_b0ROM1 << tobit_r2[63-g];
			} 
            for(int g=0;g < 64; g++){
            	H_b0ROM1 << tobit_r3[63-g];
            } 
            H_b0ROM1  << "\n";
			
			for(int g=0;g < 64; g++){
            	H_b0ROM2 << tobit_r4[63-g];
            } 
            for(int g=0;g < 64; g++){
            	H_b0ROM2 << tobit_r5[63-g];
            } 
			H_b0ROM2 << "\n";
			
            
			for(int g=0;g < 64; g++){
            	H_b0ROM3 << tobit_r6[63-g];
            } 
            for(int g=0;g < 64; g++){
				H_b0ROM3 << tobit_r7[63-g];
            } 
			H_b0ROM3 << "\n";

			for(int g=0;g < 64; g++){
            	H_b0ROM4 << tobit_r8[63-g];
            } 
			for(int g=0;g < 64; g++){
            	H_b0ROM4 << tobit_r9[63-g];
            } 
			H_b0ROM4 << "\n";
			
			for(int g=0;g < 64; g++){
				H_b0ROM5 << tobit_r10[63-g];
			} 
			for(int g=0;g < 64; g++){
				H_b0ROM5 << tobit_r11[63-g];
			} 
			H_b0ROM5 << "\n";
			
			for(int g=0;g < 64; g++){
				H_b0ROM6 << tobit_r12[63-g];
			} 
			for(int g=0;g < 64; g++){
				H_b0ROM6 << tobit_r13[63-g];
			} 
			H_b0ROM6 << "\n";
			
			for(int g=0;g < 64; g++){
				H_b0ROM7 << tobit_r14[63-g];
			} 
			for(int g=0;g < 64; g++){
				H_b0ROM7 << tobit_r15[63-g];
			} 
			H_b0ROM7 << "\n";
			
			
			
            loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
    
    counter = 0;
    loop_index = 0;
    //bank 1
    while(counter < counter_iteration){ 
        if((ima[loop_index] == counter )&& (ibn[loop_index] == 1) ){
            H_NTT_tmp_r0  = H_NTT[loop_index              ];
            H_NTT_tmp_r1  = H_NTT[loop_index + INTT_R1_offset ];
            H_NTT_tmp_r2  = H_NTT[loop_index + INTT_R2_offset ];
            H_NTT_tmp_r3  = H_NTT[loop_index + INTT_R3_offset ];
            H_NTT_tmp_r4  = H_NTT[loop_index + INTT_R4_offset ];
            H_NTT_tmp_r5  = H_NTT[loop_index + INTT_R5_offset ];
            H_NTT_tmp_r6  = H_NTT[loop_index + INTT_R6_offset ];
            H_NTT_tmp_r7  = H_NTT[loop_index + INTT_R7_offset ];
            H_NTT_tmp_r8  = H_NTT[loop_index + INTT_R8_offset ];
            H_NTT_tmp_r9  = H_NTT[loop_index + INTT_R9_offset ];
            H_NTT_tmp_r10 = H_NTT[loop_index + INTT_R10_offset];
            H_NTT_tmp_r11 = H_NTT[loop_index + INTT_R11_offset];
            H_NTT_tmp_r12 = H_NTT[loop_index + INTT_R12_offset];
            H_NTT_tmp_r13 = H_NTT[loop_index + INTT_R13_offset];
            H_NTT_tmp_r14 = H_NTT[loop_index + INTT_R14_offset];
            H_NTT_tmp_r15 = H_NTT[loop_index + INTT_R15_offset];
            
            for(int bit_index=0; bit_index < 64 ;bit_index++){
                if(H_NTT_tmp_r0 % 2 == 1) tobit_r0[bit_index] = 1;
                else tobit_r0[bit_index] = 0;
				if(H_NTT_tmp_r1 % 2 == 1) tobit_r1[bit_index] = 1;
				else tobit_r1[bit_index] = 0;
				if(H_NTT_tmp_r2 % 2 == 1) tobit_r2[bit_index] = 1;
				else tobit_r2[bit_index] = 0;
				if(H_NTT_tmp_r3 % 2 == 1) tobit_r3[bit_index] = 1;
				else tobit_r3[bit_index] = 0;
				if(H_NTT_tmp_r4 % 2 == 1) tobit_r4[bit_index] = 1;
				else tobit_r4[bit_index] = 0;
				if(H_NTT_tmp_r5 % 2 == 1) tobit_r5[bit_index] = 1;
				else tobit_r5[bit_index] = 0;
				if(H_NTT_tmp_r6 % 2 == 1) tobit_r6[bit_index] = 1;
				else tobit_r6[bit_index] = 0;
				if(H_NTT_tmp_r7 % 2 == 1) tobit_r7[bit_index] = 1;
				else tobit_r7[bit_index] = 0;
				if(H_NTT_tmp_r8 % 2 == 1) tobit_r8[bit_index] = 1;
				else tobit_r8[bit_index] = 0;
				if(H_NTT_tmp_r9 % 2 == 1) tobit_r9[bit_index] = 1;
				else tobit_r9[bit_index] = 0;
				if(H_NTT_tmp_r10 % 2 == 1) tobit_r10[bit_index] = 1;
				else tobit_r10[bit_index] = 0;
				if(H_NTT_tmp_r11 % 2 == 1) tobit_r11[bit_index] = 1;
				else tobit_r11[bit_index] = 0;
				if(H_NTT_tmp_r12 % 2 == 1) tobit_r12[bit_index] = 1;
				else tobit_r12[bit_index] = 0;
				if(H_NTT_tmp_r13 % 2 == 1) tobit_r13[bit_index] = 1;
				else tobit_r13[bit_index] = 0;
				if(H_NTT_tmp_r14 % 2 == 1) tobit_r14[bit_index] = 1;
				else tobit_r14[bit_index] = 0;
				if(H_NTT_tmp_r15 % 2 == 1) tobit_r15[bit_index] = 1;
				else tobit_r15[bit_index] = 0;
				
                H_NTT_tmp_r0  = H_NTT_tmp_r0 >> 1;
                H_NTT_tmp_r1  = H_NTT_tmp_r1 >> 1;
                H_NTT_tmp_r2  = H_NTT_tmp_r2 >> 1;
                H_NTT_tmp_r3  = H_NTT_tmp_r3 >> 1;
                H_NTT_tmp_r4  = H_NTT_tmp_r4 >> 1;
                H_NTT_tmp_r5  = H_NTT_tmp_r5 >> 1;
                H_NTT_tmp_r6  = H_NTT_tmp_r6 >> 1;
                H_NTT_tmp_r7  = H_NTT_tmp_r7 >> 1;
                H_NTT_tmp_r8  = H_NTT_tmp_r8 >> 1;
                H_NTT_tmp_r9  = H_NTT_tmp_r9 >> 1;
                H_NTT_tmp_r10 = H_NTT_tmp_r10 >> 1;
                H_NTT_tmp_r11 = H_NTT_tmp_r11 >> 1;
                H_NTT_tmp_r12 = H_NTT_tmp_r12 >> 1;
                H_NTT_tmp_r13 = H_NTT_tmp_r13 >> 1;
                H_NTT_tmp_r14 = H_NTT_tmp_r14 >> 1;
                H_NTT_tmp_r15 = H_NTT_tmp_r15 >> 1;
            }
            
            for(int g=0;g < 64; g++){
				H_b1ROM0 << tobit_r0[63-g];
			}
            for(int g=0;g < 64; g++){ 
				H_b1ROM0 << tobit_r1[63-g];
			} 
			H_b1ROM0  << "\n";
            
			for(int g=0;g < 64; g++){
				H_b1ROM1 << tobit_r2[63-g];
			} 
            for(int g=0;g < 64; g++){
            	H_b1ROM1 << tobit_r3[63-g];
            } 
            H_b1ROM1  << "\n";
			
			for(int g=0;g < 64; g++){
            	H_b1ROM2 << tobit_r4[63-g];
            } 
            for(int g=0;g < 64; g++){
            	H_b1ROM2 << tobit_r5[63-g];
            } 
			H_b1ROM2 << "\n";
			
            
			for(int g=0;g < 64; g++){
            	H_b1ROM3 << tobit_r6[63-g];
            } 
            for(int g=0;g < 64; g++){
				H_b1ROM3 << tobit_r7[63-g];
            } 
			H_b1ROM3 << "\n";

			for(int g=0;g < 64; g++){
            	H_b1ROM4 << tobit_r8[63-g];
            } 
			for(int g=0;g < 64; g++){
            	H_b1ROM4 << tobit_r9[63-g];
            } 
			H_b1ROM4 << "\n";
			
			for(int g=0;g < 64; g++){
				H_b1ROM5 << tobit_r10[63-g];
			} 
			for(int g=0;g < 64; g++){
				H_b1ROM5 << tobit_r11[63-g];
			} 
			H_b1ROM5 << "\n";
			
			for(int g=0;g < 64; g++){
				H_b1ROM6 << tobit_r12[63-g];
			} 
			for(int g=0;g < 64; g++){
				H_b1ROM6 << tobit_r13[63-g];
			} 
			H_b1ROM6 << "\n";
			
			for(int g=0;g < 64; g++){
				H_b1ROM7 << tobit_r14[63-g];
			} 
			for(int g=0;g < 64; g++){
				H_b1ROM7 << tobit_r15[63-g];
			} 
			H_b1ROM7 << "\n";
            
            loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
}
void SPMB::H_freq_o_r16_r8(std::vector<ZZ> H_NTT){
    //1 word  128 bits
	//bank 0
    std::ofstream H_freq_out("./ROM_Data/H_freq_out.txt");
    std::ofstream H_b0ROM0("./ROM_Data/H_b0ROM0.txt");
    std::ofstream H_b0ROM1("./ROM_Data/H_b0ROM1.txt");
    std::ofstream H_b0ROM2("./ROM_Data/H_b0ROM2.txt");
    std::ofstream H_b0ROM3("./ROM_Data/H_b0ROM3.txt");
    std::ofstream H_b0ROM4("./ROM_Data/H_b0ROM4.txt");
    std::ofstream H_b0ROM5("./ROM_Data/H_b0ROM5.txt");
    std::ofstream H_b0ROM6("./ROM_Data/H_b0ROM6.txt");
    std::ofstream H_b0ROM7("./ROM_Data/H_b0ROM7.txt");
    //bank1                     ROM_Data
    std::ofstream H_b1ROM0("./ROM_Data/H_b1ROM0.txt");
    std::ofstream H_b1ROM1("./ROM_Data/H_b1ROM1.txt");
    std::ofstream H_b1ROM2("./ROM_Data/H_b1ROM2.txt");
    std::ofstream H_b1ROM3("./ROM_Data/H_b1ROM3.txt");
    std::ofstream H_b1ROM4("./ROM_Data/H_b1ROM4.txt");
    std::ofstream H_b1ROM5("./ROM_Data/H_b1ROM5.txt");
    std::ofstream H_b1ROM6("./ROM_Data/H_b1ROM6.txt");
    std::ofstream H_b1ROM7("./ROM_Data/H_b1ROM7.txt");

    for(int i=0;i<fft_point;i++){
        H_freq_out << H_NTT[i];
        H_freq_out << "\n";
    }
    
    std::vector<int> tobit_r0;
    std::vector<int> tobit_r1;
    std::vector<int> tobit_r2;
    std::vector<int> tobit_r3;
    std::vector<int> tobit_r4;
    std::vector<int> tobit_r5;
    std::vector<int> tobit_r6;
    std::vector<int> tobit_r7;
    std::vector<int> tobit_r8;
    std::vector<int> tobit_r9;
    std::vector<int> tobit_r10;
    std::vector<int> tobit_r11;
    std::vector<int> tobit_r12;
    std::vector<int> tobit_r13;
    std::vector<int> tobit_r14;
    std::vector<int> tobit_r15;
    tobit_r0.resize(64);
    tobit_r1.resize(64);
    tobit_r2.resize(64);
    tobit_r3.resize(64);
    tobit_r4.resize(64);
    tobit_r5.resize(64);
    tobit_r6.resize(64);
    tobit_r7.resize(64);
    tobit_r8.resize(64);
    tobit_r9.resize(64);
    tobit_r10.resize(64);
    tobit_r11.resize(64);
    tobit_r12.resize(64);
    tobit_r13.resize(64);
    tobit_r14.resize(64);
    tobit_r15.resize(64);
    
    ZZ H_NTT_tmp_r0;
    ZZ H_NTT_tmp_r1;
    ZZ H_NTT_tmp_r2;
    ZZ H_NTT_tmp_r3;
    ZZ H_NTT_tmp_r4;
    ZZ H_NTT_tmp_r5;
    ZZ H_NTT_tmp_r6;
    ZZ H_NTT_tmp_r7;
    ZZ H_NTT_tmp_r8;
    ZZ H_NTT_tmp_r9;
    ZZ H_NTT_tmp_r10;
    ZZ H_NTT_tmp_r11;
    ZZ H_NTT_tmp_r12;
    ZZ H_NTT_tmp_r13;
    ZZ H_NTT_tmp_r14;
    ZZ H_NTT_tmp_r15;
    
    int counter =0;
    int loop_index = 0;
	int INTT_R1_offset;
	int INTT_R2_offset;
	int INTT_R3_offset;
	int INTT_R4_offset;
	int INTT_R5_offset;
	int INTT_R6_offset;
	int INTT_R7_offset;
	int INTT_R8_offset;
	int INTT_R9_offset;
	int INTT_R10_offset;
	int INTT_R11_offset;
	int INTT_R12_offset;
	int INTT_R13_offset;
	int INTT_R14_offset;
	int INTT_R15_offset;     

	//calculate INTT offset index 
	INTT_R1_offset  = (int) fft_point / 8;
	INTT_R2_offset  = (int) fft_point / 4;
	INTT_R3_offset  = INTT_R1_offset + INTT_R2_offset;
	INTT_R4_offset  = (int)fft_point / 2;
	INTT_R5_offset  = INTT_R4_offset + INTT_R1_offset;
	INTT_R6_offset  = INTT_R4_offset + INTT_R2_offset;
	INTT_R7_offset  = INTT_R4_offset + INTT_R3_offset;
	INTT_R8_offset  = (int)fft_point / 16;
	INTT_R9_offset  = INTT_R8_offset + INTT_R1_offset ;
	INTT_R10_offset = INTT_R8_offset + INTT_R2_offset ;
	INTT_R11_offset = INTT_R8_offset + INTT_R3_offset ;
	INTT_R12_offset = INTT_R8_offset + INTT_R4_offset ;
	INTT_R13_offset = INTT_R8_offset + INTT_R5_offset ;
	INTT_R14_offset = INTT_R8_offset + INTT_R6_offset ;
	INTT_R15_offset = INTT_R8_offset + INTT_R7_offset ;     

    
    //bank0
    while(counter < counter_iteration){ 
        if((ima[loop_index] == counter )&& (ibn[loop_index] == 0) ){
            
            H_NTT_tmp_r0  = H_NTT[loop_index              ];
            H_NTT_tmp_r1  = H_NTT[loop_index + INTT_R1_offset  ];
            H_NTT_tmp_r2  = H_NTT[loop_index + INTT_R2_offset  ];
            H_NTT_tmp_r3  = H_NTT[loop_index + INTT_R3_offset  ];
            H_NTT_tmp_r4  = H_NTT[loop_index + INTT_R4_offset  ];
            H_NTT_tmp_r5  = H_NTT[loop_index + INTT_R5_offset  ];
            H_NTT_tmp_r6  = H_NTT[loop_index + INTT_R6_offset  ];
            H_NTT_tmp_r7  = H_NTT[loop_index + INTT_R7_offset  ];
            H_NTT_tmp_r8  = H_NTT[loop_index + INTT_R8_offset  ];
            H_NTT_tmp_r9  = H_NTT[loop_index + INTT_R9_offset  ];
            H_NTT_tmp_r10 = H_NTT[loop_index + INTT_R10_offset ];
            H_NTT_tmp_r11 = H_NTT[loop_index + INTT_R11_offset ];
            H_NTT_tmp_r12 = H_NTT[loop_index + INTT_R12_offset ];
            H_NTT_tmp_r13 = H_NTT[loop_index + INTT_R13_offset ];
            H_NTT_tmp_r14 = H_NTT[loop_index + INTT_R14_offset ];
            H_NTT_tmp_r15 = H_NTT[loop_index + INTT_R15_offset ];
            
            for(int bit_index=0; bit_index < 64 ;bit_index++){
                if(H_NTT_tmp_r0 % 2 == 1) tobit_r0[bit_index] = 1;
                else tobit_r0[bit_index] = 0;
				if(H_NTT_tmp_r1 % 2 == 1) tobit_r1[bit_index] = 1;
				else tobit_r1[bit_index] = 0;
				if(H_NTT_tmp_r2 % 2 == 1) tobit_r2[bit_index] = 1;
				else tobit_r2[bit_index] = 0;
				if(H_NTT_tmp_r3 % 2 == 1) tobit_r3[bit_index] = 1;
				else tobit_r3[bit_index] = 0;
				if(H_NTT_tmp_r4 % 2 == 1) tobit_r4[bit_index] = 1;
				else tobit_r4[bit_index] = 0;
				if(H_NTT_tmp_r5 % 2 == 1) tobit_r5[bit_index] = 1;
				else tobit_r5[bit_index] = 0;
				if(H_NTT_tmp_r6 % 2 == 1) tobit_r6[bit_index] = 1;
				else tobit_r6[bit_index] = 0;
				if(H_NTT_tmp_r7 % 2 == 1) tobit_r7[bit_index] = 1;
				else tobit_r7[bit_index] = 0;
				if(H_NTT_tmp_r8 % 2 == 1) tobit_r8[bit_index] = 1;
				else tobit_r8[bit_index] = 0;
				if(H_NTT_tmp_r9 % 2 == 1) tobit_r9[bit_index] = 1;
				else tobit_r9[bit_index] = 0;
				if(H_NTT_tmp_r10 % 2 == 1) tobit_r10[bit_index] = 1;
				else tobit_r10[bit_index] = 0;
				if(H_NTT_tmp_r11 % 2 == 1) tobit_r11[bit_index] = 1;
				else tobit_r11[bit_index] = 0;
				if(H_NTT_tmp_r12 % 2 == 1) tobit_r12[bit_index] = 1;
				else tobit_r12[bit_index] = 0;
				if(H_NTT_tmp_r13 % 2 == 1) tobit_r13[bit_index] = 1;
				else tobit_r13[bit_index] = 0;
				if(H_NTT_tmp_r14 % 2 == 1) tobit_r14[bit_index] = 1;
				else tobit_r14[bit_index] = 0;
				if(H_NTT_tmp_r15 % 2 == 1) tobit_r15[bit_index] = 1;
				else tobit_r15[bit_index] = 0;
				
                H_NTT_tmp_r0  = H_NTT_tmp_r0 >> 1;
                H_NTT_tmp_r1  = H_NTT_tmp_r1 >> 1;
                H_NTT_tmp_r2  = H_NTT_tmp_r2 >> 1;
                H_NTT_tmp_r3  = H_NTT_tmp_r3 >> 1;
                H_NTT_tmp_r4  = H_NTT_tmp_r4 >> 1;
                H_NTT_tmp_r5  = H_NTT_tmp_r5 >> 1;
                H_NTT_tmp_r6  = H_NTT_tmp_r6 >> 1;
                H_NTT_tmp_r7  = H_NTT_tmp_r7 >> 1;
                H_NTT_tmp_r8  = H_NTT_tmp_r8 >> 1;
                H_NTT_tmp_r9  = H_NTT_tmp_r9 >> 1;
                H_NTT_tmp_r10 = H_NTT_tmp_r10 >> 1;
                H_NTT_tmp_r11 = H_NTT_tmp_r11 >> 1;
                H_NTT_tmp_r12 = H_NTT_tmp_r12 >> 1;
                H_NTT_tmp_r13 = H_NTT_tmp_r13 >> 1;
                H_NTT_tmp_r14 = H_NTT_tmp_r14 >> 1;
                H_NTT_tmp_r15 = H_NTT_tmp_r15 >> 1;
            }
            
            for(int g=0;g < 64; g++){
				H_b0ROM0 << tobit_r0[63-g];
			}
            for(int g=0;g < 64; g++){ 
				H_b0ROM0 << tobit_r1[63-g];
			} 
			H_b0ROM0  << "\n";
            
			for(int g=0;g < 64; g++){
				H_b0ROM1 << tobit_r2[63-g];
			} 
            for(int g=0;g < 64; g++){
            	H_b0ROM1 << tobit_r3[63-g];
            } 
            H_b0ROM1  << "\n";
			
			for(int g=0;g < 64; g++){
            	H_b0ROM2 << tobit_r4[63-g];
            } 
            for(int g=0;g < 64; g++){
            	H_b0ROM2 << tobit_r5[63-g];
            } 
			H_b0ROM2 << "\n";
			
            
			for(int g=0;g < 64; g++){
            	H_b0ROM3 << tobit_r6[63-g];
            } 
            for(int g=0;g < 64; g++){
				H_b0ROM3 << tobit_r7[63-g];
            } 
			H_b0ROM3 << "\n";

			for(int g=0;g < 64; g++){
            	H_b0ROM4 << tobit_r8[63-g];
            } 
			for(int g=0;g < 64; g++){
            	H_b0ROM4 << tobit_r9[63-g];
            } 
			H_b0ROM4 << "\n";
			
			for(int g=0;g < 64; g++){
				H_b0ROM5 << tobit_r10[63-g];
			} 
			for(int g=0;g < 64; g++){
				H_b0ROM5 << tobit_r11[63-g];
			} 
			H_b0ROM5 << "\n";
			
			for(int g=0;g < 64; g++){
				H_b0ROM6 << tobit_r12[63-g];
			} 
			for(int g=0;g < 64; g++){
				H_b0ROM6 << tobit_r13[63-g];
			} 
			H_b0ROM6 << "\n";
			
			for(int g=0;g < 64; g++){
				H_b0ROM7 << tobit_r14[63-g];
			} 
			for(int g=0;g < 64; g++){
				H_b0ROM7 << tobit_r15[63-g];
			} 
			H_b0ROM7 << "\n";
			
			
			
            loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
    
    counter = 0;
    loop_index = 0;
    //bank 1
    while(counter < counter_iteration){ 
        if((ima[loop_index] == counter )&& (ibn[loop_index] == 1) ){
            H_NTT_tmp_r0  = H_NTT[loop_index              ];
            H_NTT_tmp_r1  = H_NTT[loop_index + INTT_R1_offset ];
            H_NTT_tmp_r2  = H_NTT[loop_index + INTT_R2_offset ];
            H_NTT_tmp_r3  = H_NTT[loop_index + INTT_R3_offset ];
            H_NTT_tmp_r4  = H_NTT[loop_index + INTT_R4_offset ];
            H_NTT_tmp_r5  = H_NTT[loop_index + INTT_R5_offset ];
            H_NTT_tmp_r6  = H_NTT[loop_index + INTT_R6_offset ];
            H_NTT_tmp_r7  = H_NTT[loop_index + INTT_R7_offset ];
            H_NTT_tmp_r8  = H_NTT[loop_index + INTT_R8_offset ];
            H_NTT_tmp_r9  = H_NTT[loop_index + INTT_R9_offset ];
            H_NTT_tmp_r10 = H_NTT[loop_index + INTT_R10_offset];
            H_NTT_tmp_r11 = H_NTT[loop_index + INTT_R11_offset];
            H_NTT_tmp_r12 = H_NTT[loop_index + INTT_R12_offset];
            H_NTT_tmp_r13 = H_NTT[loop_index + INTT_R13_offset];
            H_NTT_tmp_r14 = H_NTT[loop_index + INTT_R14_offset];
            H_NTT_tmp_r15 = H_NTT[loop_index + INTT_R15_offset];
            
            for(int bit_index=0; bit_index < 64 ;bit_index++){
                if(H_NTT_tmp_r0 % 2 == 1) tobit_r0[bit_index] = 1;
                else tobit_r0[bit_index] = 0;
				if(H_NTT_tmp_r1 % 2 == 1) tobit_r1[bit_index] = 1;
				else tobit_r1[bit_index] = 0;
				if(H_NTT_tmp_r2 % 2 == 1) tobit_r2[bit_index] = 1;
				else tobit_r2[bit_index] = 0;
				if(H_NTT_tmp_r3 % 2 == 1) tobit_r3[bit_index] = 1;
				else tobit_r3[bit_index] = 0;
				if(H_NTT_tmp_r4 % 2 == 1) tobit_r4[bit_index] = 1;
				else tobit_r4[bit_index] = 0;
				if(H_NTT_tmp_r5 % 2 == 1) tobit_r5[bit_index] = 1;
				else tobit_r5[bit_index] = 0;
				if(H_NTT_tmp_r6 % 2 == 1) tobit_r6[bit_index] = 1;
				else tobit_r6[bit_index] = 0;
				if(H_NTT_tmp_r7 % 2 == 1) tobit_r7[bit_index] = 1;
				else tobit_r7[bit_index] = 0;
				if(H_NTT_tmp_r8 % 2 == 1) tobit_r8[bit_index] = 1;
				else tobit_r8[bit_index] = 0;
				if(H_NTT_tmp_r9 % 2 == 1) tobit_r9[bit_index] = 1;
				else tobit_r9[bit_index] = 0;
				if(H_NTT_tmp_r10 % 2 == 1) tobit_r10[bit_index] = 1;
				else tobit_r10[bit_index] = 0;
				if(H_NTT_tmp_r11 % 2 == 1) tobit_r11[bit_index] = 1;
				else tobit_r11[bit_index] = 0;
				if(H_NTT_tmp_r12 % 2 == 1) tobit_r12[bit_index] = 1;
				else tobit_r12[bit_index] = 0;
				if(H_NTT_tmp_r13 % 2 == 1) tobit_r13[bit_index] = 1;
				else tobit_r13[bit_index] = 0;
				if(H_NTT_tmp_r14 % 2 == 1) tobit_r14[bit_index] = 1;
				else tobit_r14[bit_index] = 0;
				if(H_NTT_tmp_r15 % 2 == 1) tobit_r15[bit_index] = 1;
				else tobit_r15[bit_index] = 0;
				
                H_NTT_tmp_r0  = H_NTT_tmp_r0 >> 1;
                H_NTT_tmp_r1  = H_NTT_tmp_r1 >> 1;
                H_NTT_tmp_r2  = H_NTT_tmp_r2 >> 1;
                H_NTT_tmp_r3  = H_NTT_tmp_r3 >> 1;
                H_NTT_tmp_r4  = H_NTT_tmp_r4 >> 1;
                H_NTT_tmp_r5  = H_NTT_tmp_r5 >> 1;
                H_NTT_tmp_r6  = H_NTT_tmp_r6 >> 1;
                H_NTT_tmp_r7  = H_NTT_tmp_r7 >> 1;
                H_NTT_tmp_r8  = H_NTT_tmp_r8 >> 1;
                H_NTT_tmp_r9  = H_NTT_tmp_r9 >> 1;
                H_NTT_tmp_r10 = H_NTT_tmp_r10 >> 1;
                H_NTT_tmp_r11 = H_NTT_tmp_r11 >> 1;
                H_NTT_tmp_r12 = H_NTT_tmp_r12 >> 1;
                H_NTT_tmp_r13 = H_NTT_tmp_r13 >> 1;
                H_NTT_tmp_r14 = H_NTT_tmp_r14 >> 1;
                H_NTT_tmp_r15 = H_NTT_tmp_r15 >> 1;
            }
            
            for(int g=0;g < 64; g++){
				H_b1ROM0 << tobit_r0[63-g];
			}
            for(int g=0;g < 64; g++){ 
				H_b1ROM0 << tobit_r1[63-g];
			} 
			H_b1ROM0  << "\n";
            
			for(int g=0;g < 64; g++){
				H_b1ROM1 << tobit_r2[63-g];
			} 
            for(int g=0;g < 64; g++){
            	H_b1ROM1 << tobit_r3[63-g];
            } 
            H_b1ROM1  << "\n";
			
			for(int g=0;g < 64; g++){
            	H_b1ROM2 << tobit_r4[63-g];
            } 
            for(int g=0;g < 64; g++){
            	H_b1ROM2 << tobit_r5[63-g];
            } 
			H_b1ROM2 << "\n";
			
            
			for(int g=0;g < 64; g++){
            	H_b1ROM3 << tobit_r6[63-g];
            } 
            for(int g=0;g < 64; g++){
				H_b1ROM3 << tobit_r7[63-g];
            } 
			H_b1ROM3 << "\n";

			for(int g=0;g < 64; g++){
            	H_b1ROM4 << tobit_r8[63-g];
            } 
			for(int g=0;g < 64; g++){
            	H_b1ROM4 << tobit_r9[63-g];
            } 
			H_b1ROM4 << "\n";
			
			for(int g=0;g < 64; g++){
				H_b1ROM5 << tobit_r10[63-g];
			} 
			for(int g=0;g < 64; g++){
				H_b1ROM5 << tobit_r11[63-g];
			} 
			H_b1ROM5 << "\n";
			
			for(int g=0;g < 64; g++){
				H_b1ROM6 << tobit_r12[63-g];
			} 
			for(int g=0;g < 64; g++){
				H_b1ROM6 << tobit_r13[63-g];
			} 
			H_b1ROM6 << "\n";
			
			for(int g=0;g < 64; g++){
				H_b1ROM7 << tobit_r14[63-g];
			} 
			for(int g=0;g < 64; g++){
				H_b1ROM7 << tobit_r15[63-g];
			} 
			H_b1ROM7 << "\n";
            
            loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
}
//Mixed radix function
//those function calculate the index is be located at No. address and BN
//(BFFT)At Point-wise multplication phase 
void SPMB::R16_Mixed_Radix_Pointwise_Mult_Index_AGU(int index,int &MA,int &BN){
	 unsigned long Mixed_radix;
	 
	 Mixed_radix = fft_point;
	 
	 while( (Mixed_radix % 16) == 0){
		 Mixed_radix = Mixed_radix / 16;
	 }
	
	 //std::cout << "R16_Mixed_Radix_Pointwise_Mult_Index_AGU!!\n";
	 if(Mixed_radix == 2) R16_R2_Index_AGU(index,MA,BN);
	 else if(Mixed_radix == 4) R16_R4_Index_AGU(index,MA,BN);
	 else if(Mixed_radix == 8) R16_R8_Index_AGU(index,MA,BN);
	 else std::cout << "error!!! \n";
}

void SPMB::R16_R2_Index_AGU(int index,int &MA,int &BN){
    int     fft_bit_length;
	int     BC_bit_length;
	int     BC_R16_R2_tmp;
	int     BN_R16_R2_tmp;
    int     index_tmp;
	int     bit_tmp;
	int     weight_tmp;
	int     Number_of_digit;	
	
	//std::cout << "R16_R2_Index_AGU!\n" ;
	//-----------------------------------------------------
	std::vector<int> index_BA;
	std::vector<int> BC_BR_array;
	std::vector<int> BC_BR_array_tmp;
	//-----------------------------------------------------	
    fft_bit_length = (int)ceil(log2(fft_point));
	BC_bit_length  = fft_bit_length - 4;
	index_BA.resize(fft_bit_length);
	BC_BR_array.resize(BC_bit_length);
	BC_BR_array_tmp.resize(BC_bit_length);
    //Number of digit , four bits as one digit.
	Number_of_digit = (BC_bit_length - 1) / 4;
	//-----------------------------------------------------
	index_tmp = index;
	//bit array calculate
	for(int i = 0; i < fft_bit_length; i++){
        bit_tmp    = index_tmp % 2;
        index_tmp  = index_tmp >> 1;
        index_BA[i]     = bit_tmp;
	}
	
	for(int i = 0; i < (BC_bit_length); i++){
		BC_BR_array[i]  =  index_BA[i];
		BC_BR_array_tmp[i] = index_BA[i];
	}
	
	//-----------------------------------------------------
	//index relocation
    BC_BR_array[0] = BC_BR_array_tmp[BC_bit_length-1];
	for(int i = 0;i < Number_of_digit; i++){
		BC_BR_array[(4 * i) + 4] = BC_BR_array_tmp[BC_bit_length - (4 * i) - 2];
		BC_BR_array[(4 * i) + 3] = BC_BR_array_tmp[BC_bit_length - (4 * i) - 3];
		BC_BR_array[(4 * i) + 2] = BC_BR_array_tmp[BC_bit_length - (4 * i) - 4];
		BC_BR_array[(4 * i) + 1] = BC_BR_array_tmp[BC_bit_length - (4 * i) - 5];
	}
    
	//calculate BC
    BC_R16_R2_tmp = 0;
	BN_R16_R2_tmp = 0;
    for(int j = 0; j < BC_bit_length ; j++){
		if(BC_BR_array[j] == 1) weight_tmp = 1 << j;
		else weight_tmp = 0;
		BN_R16_R2_tmp = BN_R16_R2_tmp ^ BC_BR_array[j];
		BC_R16_R2_tmp = BC_R16_R2_tmp + weight_tmp;
    }
	
	BN = BN_R16_R2_tmp;
	MA = BC_R16_R2_tmp / 2;
}
void SPMB::R16_R4_Index_AGU(int index,int &MA,int &BN){
    int     fft_bit_length;
	int     BC_bit_length;
	int     BC_R16_R4_tmp;
	int     BN_R16_R4_tmp;
    int     index_tmp;
	int     bit_tmp;
	int     weight_tmp;
	int     Number_of_digit;	
	
	//-----------------------------------------------------
	std::vector<int> index_BA;
	std::vector<int> BC_BR_array;
	std::vector<int> BC_BR_array_tmp;	

	//-----------------------------------------------------	
    fft_bit_length = (int)ceil(log2(fft_point));
	BC_bit_length  = fft_bit_length - 4;
	index_BA.resize(fft_bit_length);
	BC_BR_array.resize(BC_bit_length);
	BC_BR_array_tmp.resize(BC_bit_length);

    //Number of digit , four bits as one digit.
	Number_of_digit = (BC_bit_length - 2) / 4;
	//-----------------------------------------------------
	index_tmp = index;
	//bit array calculate
	for(int i = 0; i < fft_bit_length; i++){
        bit_tmp    = index_tmp % 2;
        index_tmp  = index_tmp >> 1;
        index_BA[i]     = bit_tmp;
	}	
	
	for(int i = 0; i < BC_bit_length; i++){
		BC_BR_array[i] = index_BA[i];
		BC_BR_array_tmp[i] = index_BA[i];
	}
	//-----------------------------------------------------
	//index relocation
    BC_BR_array[1] = BC_BR_array_tmp[BC_bit_length-1];
    BC_BR_array[0] = BC_BR_array_tmp[BC_bit_length-2];
	for(int i = 0;i < Number_of_digit; i++){
		BC_BR_array[(4 * i) + 5] = BC_BR_array_tmp[BC_bit_length - (4 * i) - 3];
		BC_BR_array[(4 * i) + 4] = BC_BR_array_tmp[BC_bit_length - (4 * i) - 4];
		BC_BR_array[(4 * i) + 3] = BC_BR_array_tmp[BC_bit_length - (4 * i) - 5];
		BC_BR_array[(4 * i) + 2] = BC_BR_array_tmp[BC_bit_length - (4 * i) - 6];
	}
    
	//calculate BC
    BC_R16_R4_tmp = 0;
	BN_R16_R4_tmp = 0;
    for(int j = 0; j < BC_bit_length ; j++){
		if(BC_BR_array[j] == 1) weight_tmp = 1 << j ;
		else weight_tmp = 0;
		BN_R16_R4_tmp = BN_R16_R4_tmp ^ BC_BR_array[j];
		BC_R16_R4_tmp = BC_R16_R4_tmp + weight_tmp;
    }
	
	BN = BN_R16_R4_tmp;
	MA = BC_R16_R4_tmp / 2;	
}
void SPMB::R16_R8_Index_AGU(int index,int &MA,int &BN){
    int     fft_bit_length;
	int     BC_bit_length;
	int     BC_R16_R8_tmp;
	int     BN_R16_R8_tmp;
    int     index_tmp;
	int     bit_tmp;
	int     weight_tmp;
	int     Number_of_digit;	
	
	//-----------------------------------------------------
	std::vector<int> index_BA;
	std::vector<int> BC_BR_array;
	std::vector<int> BC_BR_array_tmp;	
	//-----------------------------------------------------	
    fft_bit_length = (int)ceil(log2(fft_point));
	BC_bit_length  = fft_bit_length - 4;
	index_BA.resize(fft_bit_length);
	BC_BR_array.resize(BC_bit_length);
	BC_BR_array_tmp.resize(BC_bit_length);

    //Number of digit , four bits as one digit.
	Number_of_digit = (BC_bit_length - 3) / 4;
	//-----------------------------------------------------
	index_tmp = index;
	//bit array calculate
	for(int i = 0; i < fft_bit_length; i++){
        bit_tmp    = index_tmp % 2;
        index_tmp  = index_tmp >> 1;
        index_BA[i]     = bit_tmp;
	}	
	
	for(int i = 0; i < BC_bit_length; i++){
		BC_BR_array[i]     = index_BA[i];
		BC_BR_array_tmp[i] = index_BA[i];
	}	
	//-----------------------------------------------------
	//index relocation
	BC_BR_array[2]  = BC_BR_array_tmp[BC_bit_length - 1];
	BC_BR_array[1]  = BC_BR_array_tmp[BC_bit_length - 2];
	BC_BR_array[0]  = BC_BR_array_tmp[BC_bit_length - 3];
	for(int i = 0;i < Number_of_digit; i++){
		BC_BR_array[(4 * i) + 6] = BC_BR_array_tmp[BC_bit_length - (4 * i) - 4];
		BC_BR_array[(4 * i) + 5] = BC_BR_array_tmp[BC_bit_length - (4 * i) - 5];
		BC_BR_array[(4 * i) + 4] = BC_BR_array_tmp[BC_bit_length - (4 * i) - 6];
		BC_BR_array[(4 * i) + 3] = BC_BR_array_tmp[BC_bit_length - (4 * i) - 7];
	}
    
    
	//calculate BC
    BC_R16_R8_tmp = 0;
	BN_R16_R8_tmp = 0;
    for(int j = 0; j < BC_bit_length; j++){
		if(BC_BR_array[j] == 1) weight_tmp = 1 << j;
		else weight_tmp = 0;
		BN_R16_R8_tmp = BN_R16_R8_tmp ^ BC_BR_array[j];
		BC_R16_R8_tmp = BC_R16_R8_tmp + weight_tmp;
    }
	
	BN = BN_R16_R8_tmp;
	MA = BC_R16_R8_tmp / 2;		
	
}
//re-order factor output are store as rom data
//CP_width must be 22~25
//********************************************
// ROM bit size : 128 bits
//********************************************
void SPMB::re_order_factor_r16(ZZ m_2_rou){

     std::vector<ZZ> re_order_factor_array;
     re_order_factor_array.resize(fft_point);
     
     for(unsigned long i=0;i<m;i++){
        unsigned long exp;
        ZZ w_2m_tmp_i;
        exp = pow(i,2);
        exp = exp % (2*m);
        PowerMod(w_2m_tmp_i,m_2_rou,exp,cyclotomic_prime);
		//conv(w_2m_tmp_i,"1");
        //std::cout << "siang w_2m_tmp_i = " << w_2m_tmp_i << std::endl; //---------siang
        re_order_factor_array[i] = w_2m_tmp_i;
     }
    
	//cyclotomic polynomial prime must be 22~25.   
	//re_order_factor_array,while index > m then data is ZERO.
	//because using SPMB Memory addressing,then this reorder factor all index 8 ~ index 15 are zero for radix-16 FFT
    //so we just store reorder factor index 0 ~ index 7 	
    //for example prime bit = 22 ,then total 22 * 8 = 172 bits,so there are 2 ROMs ,word size = 128 bits 
    //ROM0[127:106] = radix0,ROM0[105:84] = radix1, ROM0[83:62] = radix2 , ROM0[61:40] = radix3 , ROM0[39:18] = radix4
    //ROM0[17:0],ROM1[127:124] = radix5 , ROM1[123:102] = radix6 , ROM1[101,80] = radix7

	//bank0
	std::ofstream reorder_ROM0("./ROM_Data/reorder_ROM0.txt");
    std::ofstream reorder_ROM1("./ROM_Data/reorder_ROM1.txt");
    //bank1
    std::ofstream reorder_ROM2("./ROM_Data/reorder_ROM2.txt");
    std::ofstream reorder_ROM3("./ROM_Data/reorder_ROM3.txt");
   
    std::vector<int> tobit_r0;
    std::vector<int> tobit_r1;
    std::vector<int> tobit_r2;
    std::vector<int> tobit_r3;
    std::vector<int> tobit_r4;
    std::vector<int> tobit_r5;
    std::vector<int> tobit_r6;
    std::vector<int> tobit_r7;

    tobit_r0.resize(CP_width);
    tobit_r1.resize(CP_width);
    tobit_r2.resize(CP_width);
    tobit_r3.resize(CP_width);
    tobit_r4.resize(CP_width);
    tobit_r5.resize(CP_width);
    tobit_r6.resize(CP_width);
    tobit_r7.resize(CP_width);
    
    ZZ re_order_r0_tmp;
    ZZ re_order_r1_tmp;
    ZZ re_order_r2_tmp;
    ZZ re_order_r3_tmp;
    ZZ re_order_r4_tmp;
    ZZ re_order_r5_tmp;
    ZZ re_order_r6_tmp;
    ZZ re_order_r7_tmp;
    
    int counter =0;
    int loop_index = 0;
    int loop_address = 0;

    //bank0
    while(counter < counter_iteration){ 
        if((ma[loop_index] == counter )&& (bn[loop_index] == 0) ){
            re_order_r0_tmp = re_order_factor_array[loop_index + 0];
            re_order_r1_tmp = re_order_factor_array[loop_index + offset];
            re_order_r2_tmp = re_order_factor_array[loop_index + 2 * offset];
			re_order_r3_tmp = re_order_factor_array[loop_index + 3 * offset];
			re_order_r4_tmp = re_order_factor_array[loop_index + 4 * offset];
			re_order_r5_tmp = re_order_factor_array[loop_index + 5 * offset];
			re_order_r6_tmp = re_order_factor_array[loop_index + 6 * offset];
			re_order_r7_tmp = re_order_factor_array[loop_index + 7 * offset];
			
			//calculate bit
            for(int bit_index=0; bit_index < CP_width ;bit_index++){
                //radix0
				if(re_order_r0_tmp % 2 == 1) tobit_r0[bit_index] = 1;
                else tobit_r0[bit_index] = 0;
				//radix1
				if(re_order_r1_tmp % 2 == 1) tobit_r1[bit_index] = 1;
                else tobit_r1[bit_index] = 0;
				
				if(re_order_r2_tmp % 2 == 1) tobit_r2[bit_index] = 1;
                else tobit_r2[bit_index] = 0;
				
				if(re_order_r3_tmp % 2 == 1) tobit_r3[bit_index] = 1;
                else tobit_r3[bit_index] = 0;
				
				if(re_order_r4_tmp % 2 == 1) tobit_r4[bit_index] = 1;
                else tobit_r4[bit_index] = 0;
				
				if(re_order_r5_tmp % 2 == 1) tobit_r5[bit_index] = 1;
                else tobit_r5[bit_index] = 0;
				
				if(re_order_r6_tmp % 2 == 1) tobit_r6[bit_index] = 1;
                else tobit_r6[bit_index] = 0;
				
				if(re_order_r7_tmp % 2 == 1) tobit_r7[bit_index] = 1;
                else tobit_r7[bit_index] = 0;
				
				re_order_r0_tmp = re_order_r0_tmp >> 1;
				re_order_r1_tmp = re_order_r1_tmp >> 1;
				re_order_r2_tmp = re_order_r2_tmp >> 1;
				re_order_r3_tmp = re_order_r3_tmp >> 1;
				re_order_r4_tmp = re_order_r4_tmp >> 1;
				re_order_r5_tmp = re_order_r5_tmp >> 1;
				re_order_r6_tmp = re_order_r6_tmp >> 1;
				re_order_r7_tmp = re_order_r7_tmp >> 1;
            }
            
           
            for(int g=0;g < CP_width; g++){
				reorder_ROM0 << tobit_r0[CP_width-1-g];
            } 
            for(int g=0;g < CP_width; g++){
                reorder_ROM0 << tobit_r1[CP_width-1-g];
            } 
            for(int g=0;g < CP_width; g++){
                reorder_ROM0 << tobit_r2[CP_width-1-g];
			}
			for(int g=0;g < CP_width; g++){
                reorder_ROM0 << tobit_r3[CP_width-1-g];
			}
			for(int g=0;g < CP_width; g++){
                reorder_ROM0 << tobit_r4[CP_width-1-g] ;
			}
			int ROM0_Remaining_bits;
		    ROM0_Remaining_bits = 128 - 5 * CP_width;
			
			for(int g=0;g < CP_width; g++){
                if(ROM0_Remaining_bits > g)reorder_ROM0 << tobit_r5[CP_width-1-g] ;
                else reorder_ROM1 << tobit_r5[CP_width-1-g] ;
			}
			reorder_ROM0 << "\n";
			
			for(int g=0;g < CP_width; g++){
			    reorder_ROM1 << tobit_r6[CP_width-1-g] ;
			}
			
			for(int g=0;g < CP_width; g++){
			    reorder_ROM1 << tobit_r7[CP_width-1-g] ;
			}
			
			int ROM1_padding_zero;
		    ROM1_padding_zero = 128 - 3 * CP_width + ROM0_Remaining_bits ;
			
			for(int g=0;g < ROM1_padding_zero; g++){
			    reorder_ROM1 << 0 ;
			}
			reorder_ROM1 << "\n";
			
            loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
    reorder_ROM0.close();
    reorder_ROM1.close();
	
    counter = 0;
    loop_index = 0;
    //bank 1
    while(counter < counter_iteration){ 
        if((ma[loop_index] == counter )&& (bn[loop_index] == 1) ){
            re_order_r0_tmp = re_order_factor_array[loop_index + 0];
            re_order_r1_tmp = re_order_factor_array[loop_index + offset];
            re_order_r2_tmp = re_order_factor_array[loop_index + 2 * offset];
			re_order_r3_tmp = re_order_factor_array[loop_index + 3 * offset];
			re_order_r4_tmp = re_order_factor_array[loop_index + 4 * offset];
			re_order_r5_tmp = re_order_factor_array[loop_index + 5 * offset];
			re_order_r6_tmp = re_order_factor_array[loop_index + 6 * offset];
			re_order_r7_tmp = re_order_factor_array[loop_index + 7 * offset];
			
			//calculate bit
            for(int bit_index=0; bit_index < CP_width ;bit_index++){
                //radix0
				if(re_order_r0_tmp % 2 == 1) tobit_r0[bit_index] = 1;
                else tobit_r0[bit_index] = 0;
				//radix1
				if(re_order_r1_tmp % 2 == 1) tobit_r1[bit_index] = 1;
                else tobit_r1[bit_index] = 0;
				
				if(re_order_r2_tmp % 2 == 1) tobit_r2[bit_index] = 1;
                else tobit_r2[bit_index] = 0;
				
				if(re_order_r3_tmp % 2 == 1) tobit_r3[bit_index] = 1;
                else tobit_r3[bit_index] = 0;
				
				if(re_order_r4_tmp % 2 == 1) tobit_r4[bit_index] = 1;
                else tobit_r4[bit_index] = 0;
				
				if(re_order_r5_tmp % 2 == 1) tobit_r5[bit_index] = 1;
                else tobit_r5[bit_index] = 0;
				
				if(re_order_r6_tmp % 2 == 1) tobit_r6[bit_index] = 1;
                else tobit_r6[bit_index] = 0;
				
				if(re_order_r7_tmp % 2 == 1) tobit_r7[bit_index] = 1;
                else tobit_r7[bit_index] = 0;
				
				re_order_r0_tmp = re_order_r0_tmp >> 1;
				re_order_r1_tmp = re_order_r1_tmp >> 1;
				re_order_r2_tmp = re_order_r2_tmp >> 1;
				re_order_r3_tmp = re_order_r3_tmp >> 1;
				re_order_r4_tmp = re_order_r4_tmp >> 1;
				re_order_r5_tmp = re_order_r5_tmp >> 1;
				re_order_r6_tmp = re_order_r6_tmp >> 1;
				re_order_r7_tmp = re_order_r7_tmp >> 1;
            }
            
           
            for(int g=0;g < CP_width; g++){
				reorder_ROM2 << tobit_r0[CP_width-1-g];
            } 
            for(int g=0;g < CP_width; g++){
                reorder_ROM2 << tobit_r1[CP_width-1-g];
            } 
            for(int g=0;g < CP_width; g++){
                reorder_ROM2 << tobit_r2[CP_width-1-g];
			}			
			for(int g=0;g < CP_width; g++){
                reorder_ROM2 << tobit_r3[CP_width-1-g] ;
			}
			
			for(int g=0;g < CP_width; g++){
                reorder_ROM2 << tobit_r4[CP_width-1-g] ;
			}
			int ROM2_Remaining_bits;
		    ROM2_Remaining_bits = 128 - 5 * CP_width;
			
			for(int g=0;g < CP_width; g++){
                if(ROM2_Remaining_bits > g)reorder_ROM2 << tobit_r5[CP_width-1-g] ;
                else reorder_ROM3 << tobit_r5[CP_width-1-g] ;
			}
			reorder_ROM2 << "\n";
			
			for(int g=0;g < CP_width; g++){
			    reorder_ROM3 << tobit_r6[CP_width-1-g] ;
			}
			
			for(int g=0;g < CP_width; g++){
			    reorder_ROM3 << tobit_r7[CP_width-1-g] ;
			}
			
			int ROM3_padding_zero;
		    ROM3_padding_zero = 128 - 3 * CP_width + ROM2_Remaining_bits ;
			
			for(int g=0;g < ROM3_padding_zero; g++){
			    reorder_ROM3 << 0 ;
			}
			reorder_ROM3 << "\n"; 
                 
            loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
    reorder_ROM2.close();
    reorder_ROM3.close();
	
    loop_address = 0;
    
	//inverse order factor 
	std::vector<int> itobit_r0;
	std::vector<int> itobit_r1;
	std::vector<int> itobit_r2;
	std::vector<int> itobit_r3;
	std::vector<int> itobit_r4;
	std::vector<int> itobit_r5;
	std::vector<int> itobit_r6;
	std::vector<int> itobit_r7;
	std::vector<int> itobit_r8;
	std::vector<int> itobit_r9;
	std::vector<int> itobit_r10;
	std::vector<int> itobit_r11;
	std::vector<int> itobit_r12;
	std::vector<int> itobit_r13;
	std::vector<int> itobit_r14;
	std::vector<int> itobit_r15;
	itobit_r0.resize(CP_width);
	itobit_r1.resize(CP_width);
	itobit_r2.resize(CP_width);
	itobit_r3.resize(CP_width);
	itobit_r4.resize(CP_width);
	itobit_r5.resize(CP_width);
	itobit_r6.resize(CP_width);
	itobit_r7.resize(CP_width);
	itobit_r8.resize(CP_width);
	itobit_r9.resize(CP_width);
	itobit_r10.resize(CP_width);
	itobit_r11.resize(CP_width);
	itobit_r12.resize(CP_width);
	itobit_r13.resize(CP_width);
	itobit_r14.resize(CP_width);
	itobit_r15.resize(CP_width);
	ZZ re_iorder_r0_tmp;
	ZZ re_iorder_r1_tmp;
	ZZ re_iorder_r2_tmp;
	ZZ re_iorder_r3_tmp;
	ZZ re_iorder_r4_tmp;
	ZZ re_iorder_r5_tmp;
	ZZ re_iorder_r6_tmp;
	ZZ re_iorder_r7_tmp;
	ZZ re_iorder_r8_tmp;
	ZZ re_iorder_r9_tmp;
	ZZ re_iorder_r10_tmp;
	ZZ re_iorder_r11_tmp;
	ZZ re_iorder_r12_tmp;
	ZZ re_iorder_r13_tmp;
	ZZ re_iorder_r14_tmp;
	ZZ re_iorder_r15_tmp;
	/*************************************************/
	//m-th cyclotomic polynomial 
	//all polynomial at this thing,the coefficient is zero while the degree > phi(m);
    //and inverse reorder factor index > phi(m) is zero.
	//using SPMB (Final stage type) ,the data index is continuous
	//we can divide inverse re-order factor into two banks ROM.
	//ROM0 bit width= 64 ,ROM1 bit width = 32
	//cyclotomic prime = 22 bits
	//ROM0[63:42] = radix0 ,ROM0[41:20] = radix1,ROM0[19:0] = radix2[21:2]
	//ROM1[63:62] = raidx2[1:0],ROM1[61:40] = radix3 , ROM1[39:18] = radix4 ,ROM1[17:0] = radix5[21:4]
	//ROM2[63:60] = radix5[3:0],ROM2[59:38] = radix6 , ROM2[37:16] = radix7 ,ROM2[15:0] = radix8[21:6]
	//ROM3[63:58] = radix8[5:0],ROM3[57:36] = radix9 , ROM3[35:14] = radix10,ROM3[13:0] = radix11[21:8]
	//ROM4[63:56] = radix11[7:0],ROM4[55:34]= radix12, ROM4[33:12] = radix13,ROM4[11:0] = radix14[21:10]
	//ROM5[31:22] = radix14[9:0],ROM5 = radix15
	//ROM Word size ,original word size is N/radix => 4096/4 = 1024,
	//a_0x^(0)+........+a_phi(m)-1*x^(phi(m)-1)+0*x^phi(m)+...+0^x(m-1)
	//so word size can redution to phi(m)/radix 
    //for example m=1705 phi(m) = 1200 , then ma > 300 data is zero	
	//then we set word size is 512
	//**************************************************
	// modify at 2021 / 02 / 18 
	// change rom bit size from 64 to 128
	//
	//  ROM0~ ROM2 : 128 bits , ROM3 : 32 bits
	//  ROM0 : radix0  ~ radix4
	//  ROM1 : radix5  ~ radix9
	//  ROM2 : radix10 ~ radix14
	//  ROM3 : radix15
	/*************************************************/
	std::ofstream ireorder_out("./ROM_Data/ireorder_outfile.txt");
	//128-bits
    std::ofstream ireorder_ROM0("./ROM_Data/ireorder_ROM0.txt");
    std::ofstream ireorder_ROM1("./ROM_Data/ireorder_ROM1.txt");
    std::ofstream ireorder_ROM2("./ROM_Data/ireorder_ROM2.txt");
    //64bits
	std::ofstream ireorder_ROM3("./ROM_Data/ireorder_ROM3.txt");
    
    double IReROM_WORD_BIT;
	int    ROM_WORD_SIZE;
    IReROM_WORD_BIT = ceil((double)m / 16);
    IReROM_WORD_BIT = log2(IReROM_WORD_BIT);
    IReROM_WORD_BIT = ceil(IReROM_WORD_BIT);	
    ROM_WORD_SIZE   = exp2(IReROM_WORD_BIT);
	
	for(int i=0;i<fft_point;i++){
		ireorder_out << re_order_factor_array[i] << "\n";
	}
	ireorder_out.close();
	
   
	for(int i=0;i < ROM_WORD_SIZE; i++){
         re_iorder_r0_tmp = re_order_factor_array[radix * i + 0];
         re_iorder_r1_tmp = re_order_factor_array[radix * i + 1];
         re_iorder_r2_tmp = re_order_factor_array[radix * i + 2];
         re_iorder_r3_tmp = re_order_factor_array[radix * i + 3];
         re_iorder_r4_tmp = re_order_factor_array[radix * i + 4];
         re_iorder_r5_tmp = re_order_factor_array[radix * i + 5];
         re_iorder_r6_tmp = re_order_factor_array[radix * i + 6];
         re_iorder_r7_tmp = re_order_factor_array[radix * i + 7];
         re_iorder_r8_tmp = re_order_factor_array[radix * i + 8];
         re_iorder_r9_tmp = re_order_factor_array[radix * i + 9];
         re_iorder_r10_tmp = re_order_factor_array[radix * i + 10];
         re_iorder_r11_tmp = re_order_factor_array[radix * i + 11];
         re_iorder_r12_tmp = re_order_factor_array[radix * i + 12];
         re_iorder_r13_tmp = re_order_factor_array[radix * i + 13];
         re_iorder_r14_tmp = re_order_factor_array[radix * i + 14];
         re_iorder_r15_tmp = re_order_factor_array[radix * i + 15];
        
        for(int bit_index=0; bit_index < CP_width ;bit_index++){
            //radix0
			if(re_iorder_r0_tmp % 2 == 1) itobit_r0[bit_index] = 1;
            else itobit_r0[bit_index] = 0;
			//radix1
			if(re_iorder_r1_tmp % 2 == 1) itobit_r1[bit_index] = 1;
			else itobit_r1[bit_index] = 0;
			//radix2
			if(re_iorder_r2_tmp % 2 == 1) itobit_r2[bit_index] = 1;
			else itobit_r2[bit_index] = 0;
			//radix3
            if(re_iorder_r3_tmp % 2 == 1) itobit_r3[bit_index] = 1;
			else itobit_r3[bit_index] = 0;
			//radix4
			if(re_iorder_r4_tmp % 2 == 1) itobit_r4[bit_index] = 1;
			else itobit_r4[bit_index] = 0;
			//radix5
			if(re_iorder_r5_tmp % 2 == 1) itobit_r5[bit_index] = 1;
			else itobit_r5[bit_index] = 0;
			//radix6
			if(re_iorder_r6_tmp % 2 == 1) itobit_r6[bit_index] = 1;
			else itobit_r6[bit_index] = 0;
			//radix7
			if(re_iorder_r7_tmp % 2 == 1) itobit_r7[bit_index] = 1;
			else itobit_r7[bit_index] = 0;
			//radix8
			if(re_iorder_r8_tmp % 2 == 1) itobit_r8[bit_index] = 1;
			else itobit_r8[bit_index] = 0;
			//radix9
			if(re_iorder_r9_tmp % 2 == 1) itobit_r9[bit_index] = 1;
			else itobit_r9[bit_index] = 0;
			//radix10
			if(re_iorder_r10_tmp % 2 == 1) itobit_r10[bit_index] = 1;
			else itobit_r10[bit_index] = 0;
			//radix11
			if(re_iorder_r11_tmp % 2 == 1) itobit_r11[bit_index] = 1;
			else itobit_r11[bit_index] = 0;
			//radix12
			if(re_iorder_r12_tmp % 2 == 1) itobit_r12[bit_index] = 1;
			else itobit_r12[bit_index] = 0;
			//radix13
			if(re_iorder_r13_tmp % 2 == 1) itobit_r13[bit_index] = 1;
			else itobit_r13[bit_index] = 0;
			//radix14
			if(re_iorder_r14_tmp % 2 == 1) itobit_r14[bit_index] = 1;
			else itobit_r14[bit_index] = 0;
			//radix15
			if(re_iorder_r15_tmp % 2 == 1) itobit_r15[bit_index] = 1;
			else itobit_r15[bit_index] = 0;
			
			re_iorder_r0_tmp  = re_iorder_r0_tmp >> 1;
            re_iorder_r1_tmp  = re_iorder_r1_tmp >> 1;
            re_iorder_r2_tmp  = re_iorder_r2_tmp >> 1;
            re_iorder_r3_tmp  = re_iorder_r3_tmp >> 1;
            re_iorder_r4_tmp  = re_iorder_r4_tmp >> 1;
            re_iorder_r5_tmp  = re_iorder_r5_tmp >> 1;
            re_iorder_r6_tmp  = re_iorder_r6_tmp >> 1;
            re_iorder_r7_tmp  = re_iorder_r7_tmp >> 1;
            re_iorder_r8_tmp  = re_iorder_r8_tmp >> 1;
            re_iorder_r9_tmp  = re_iorder_r9_tmp >> 1;
            re_iorder_r10_tmp = re_iorder_r10_tmp >> 1;
            re_iorder_r11_tmp = re_iorder_r11_tmp >> 1;
            re_iorder_r12_tmp = re_iorder_r12_tmp >> 1;
            re_iorder_r13_tmp = re_iorder_r13_tmp >> 1;
            re_iorder_r14_tmp = re_iorder_r14_tmp >> 1;
            re_iorder_r15_tmp = re_iorder_r15_tmp >> 1;
        }
        
		//--------------------------------------------------------
		//radix0 ~ radix4
		//ROM0[127:106] = radix0
        for(int g=0;g < CP_width; g++){
            ireorder_ROM0 << itobit_r0[CP_width-1-g];
        } 
        //ROM0[105:84] = radix1
		for(int g=0;g < CP_width; g++){
            ireorder_ROM0 << itobit_r1[CP_width-1-g];
        } 
        //ROM0[83:62] = radix2 
        for(int g=0;g < CP_width; g++){
            ireorder_ROM0 << itobit_r2[CP_width-1-g];
        } 
        //ROM0[61:40] = radix3 
        for(int g=0;g < CP_width; g++){
            ireorder_ROM0 << itobit_r3[CP_width-1-g];
        } 
		//ROM0[39:18] = radix4
		for(int g=0;g < CP_width; g++){
            ireorder_ROM0 << itobit_r4[CP_width-1-g];
        } 
		int IROM0_padding_zero;
		IROM0_padding_zero = 128 - 5 * CP_width;

		for(int g=0;g < IROM0_padding_zero; g++){
		    ireorder_ROM0 << 0 ;
		}		
		ireorder_ROM0 << "\n";
		//----------------------------------------------------
		//radix5 ~ radix9
		//ROM1[127:106] = radix5
        for(int g=0;g < CP_width; g++){
            ireorder_ROM1 << itobit_r5[CP_width-1-g];
        } 
        //ROM1[105:84] = radix6
		for(int g=0;g < CP_width; g++){
            ireorder_ROM1 << itobit_r6[CP_width-1-g];
        } 
        //ROM1[83:62] = radix7
        for(int g=0;g < CP_width; g++){
            ireorder_ROM1 << itobit_r7[CP_width-1-g];
        } 
        //ROM1[61:40] = radix8 
        for(int g=0;g < CP_width; g++){
            ireorder_ROM1 << itobit_r8[CP_width-1-g];
        } 
		//ROM1[39:18] = radix9
		for(int g=0;g < CP_width; g++){
            ireorder_ROM1 << itobit_r9[CP_width-1-g];
        } 
		int IROM1_padding_zero;
		IROM1_padding_zero = 128 - 5 * CP_width;

		for(int g=0;g < IROM1_padding_zero; g++){
		    ireorder_ROM1 << 0 ;
		}		
		ireorder_ROM1 << "\n";
		//----------------------------------------------------
		//radix10 ~ radix14
		//ROM1[127:106] = radix10
        for(int g=0;g < CP_width; g++){
            ireorder_ROM2 << itobit_r10[CP_width-1-g];
        } 
        //ROM0[105:84] = radix1
		for(int g=0;g < CP_width; g++){
            ireorder_ROM2 << itobit_r11[CP_width-1-g];
        } 
        //ROM0[83:62] = radix2 
        for(int g=0;g < CP_width; g++){
            ireorder_ROM2 << itobit_r12[CP_width-1-g];
        } 
        //ROM0[61:40] = radix3 
        for(int g=0;g < CP_width; g++){
            ireorder_ROM2 << itobit_r13[CP_width-1-g];
        } 
		//ROM0[39:18] = radix4
		for(int g=0;g < CP_width; g++){
            ireorder_ROM2 << itobit_r14[CP_width-1-g];
        } 
		int IROM2_padding_zero;
		IROM2_padding_zero = 128 - 5 * CP_width;

		for(int g=0;g < IROM2_padding_zero; g++){
		    ireorder_ROM2 << 0 ;
		}		
		ireorder_ROM2 << "\n";
        //------------------------------------------------------
        //radix15		
        for(int g=0;g < CP_width; g++){
            ireorder_ROM3 << itobit_r15[CP_width-1-g];
        } 		
		int IROM3_padding_zero;
		IROM3_padding_zero = 64 - CP_width;

		for(int g=0;g < IROM3_padding_zero; g++){
		    ireorder_ROM3 << 0 ;
		}		
		ireorder_ROM3 << "\n";
    }
    ireorder_ROM0.close();
    ireorder_ROM1.close();
    ireorder_ROM2.close();
    ireorder_ROM3.close();
}
//--------------------
//Reconfigure data generate
void SPMB::radix16_Reconfigure_DATA(std::vector<ZZ> B_NTT,ZZ m_2_rou){
     unsigned long Mixed_radix;
     
	 Mixed_radix = fft_point;
	 
	 while(Mixed_radix % 16 == 0){
		 Mixed_radix = Mixed_radix / 16; 
	 }
	 
	 //std::cout <<"H_Freq Output!!!\n";
	 if(IsMixed == 0){
		 H_freq_Reconfigure_r16(B_NTT);
	 }else {
		 if(Mixed_radix == 2) H_freq_Reconfigure_r16_r2(B_NTT);
		 if(Mixed_radix == 4) H_freq_Reconfigure_r16_r4(B_NTT);
		 if(Mixed_radix == 8) H_freq_Reconfigure_r16_r8(B_NTT);
	 }
     //std::cout <<"Reorder factor data Output!!!\n";
	 re_order_factor_Reconfigure_r16(m_2_rou);
	 //std::cout <<"Reorder factor data Output OVER!!!\n";
}
void SPMB::H_freq_Reconfigure_r16(std::vector<ZZ> H_NTT){
    //1 word  128 bits
	//bank 0
    std::ofstream H_freq_out("./Reconfigure/H_freq_out.txt");
    std::ofstream H_b0SRAM0("./Reconfigure/H_b0SRAM0.txt");
    std::ofstream H_b0SRAM1("./Reconfigure/H_b0SRAM1.txt");
    std::ofstream H_b0SRAM2("./Reconfigure/H_b0SRAM2.txt");
    std::ofstream H_b0SRAM3("./Reconfigure/H_b0SRAM3.txt");
    std::ofstream H_b0SRAM4("./Reconfigure/H_b0SRAM4.txt");
    std::ofstream H_b0SRAM5("./Reconfigure/H_b0SRAM5.txt");
    std::ofstream H_b0SRAM6("./Reconfigure/H_b0SRAM6.txt");
    std::ofstream H_b0SRAM7("./Reconfigure/H_b0SRAM7.txt");
    //bank1                     ROM_Data
    std::ofstream H_b1SRAM0("./Reconfigure/H_b1SRAM0.txt");
    std::ofstream H_b1SRAM1("./Reconfigure/H_b1SRAM1.txt");
    std::ofstream H_b1SRAM2("./Reconfigure/H_b1SRAM2.txt");
    std::ofstream H_b1SRAM3("./Reconfigure/H_b1SRAM3.txt");
    std::ofstream H_b1SRAM4("./Reconfigure/H_b1SRAM4.txt");
    std::ofstream H_b1SRAM5("./Reconfigure/H_b1SRAM5.txt");
    std::ofstream H_b1SRAM6("./Reconfigure/H_b1SRAM6.txt");
    std::ofstream H_b1SRAM7("./Reconfigure/H_b1SRAM7.txt");
    
    std::string H_NTT_HEX_r0;
    std::string H_NTT_HEX_r1;
    std::string H_NTT_HEX_r2;
    std::string H_NTT_HEX_r3;
    std::string H_NTT_HEX_r4;
    std::string H_NTT_HEX_r5;
    std::string H_NTT_HEX_r6;
    std::string H_NTT_HEX_r7;
    std::string H_NTT_HEX_r8;
    std::string H_NTT_HEX_r9;
    std::string H_NTT_HEX_r10;
    std::string H_NTT_HEX_r11;
    std::string H_NTT_HEX_r12;
    std::string H_NTT_HEX_r13;
    std::string H_NTT_HEX_r14;
    std::string H_NTT_HEX_r15;
    
    ZZ H_NTT_tmp_r0;
    ZZ H_NTT_tmp_r1;
    ZZ H_NTT_tmp_r2;
    ZZ H_NTT_tmp_r3;
    ZZ H_NTT_tmp_r4;
    ZZ H_NTT_tmp_r5;
    ZZ H_NTT_tmp_r6;
    ZZ H_NTT_tmp_r7;
    ZZ H_NTT_tmp_r8;
    ZZ H_NTT_tmp_r9;
    ZZ H_NTT_tmp_r10;
    ZZ H_NTT_tmp_r11;
    ZZ H_NTT_tmp_r12;
    ZZ H_NTT_tmp_r13;
    ZZ H_NTT_tmp_r14;
    ZZ H_NTT_tmp_r15;
    
    int counter =0;
    int loop_index = 0;
     
    for(int i=0;i<fft_point;i++){
        H_freq_out << H_NTT[i];
        H_freq_out << "\n";
    }
    
    //bank0
    while(counter < counter_iteration){ 
        if((ima[loop_index] == counter )&& (ibn[loop_index] == 0) ){
            
            H_NTT_tmp_r0  = H_NTT[loop_index              ];
            H_NTT_tmp_r1  = H_NTT[loop_index + 1  * offset ];
            H_NTT_tmp_r2  = H_NTT[loop_index + 2  * offset ];
            H_NTT_tmp_r3  = H_NTT[loop_index + 3  * offset ];
            H_NTT_tmp_r4  = H_NTT[loop_index + 4  * offset ];
            H_NTT_tmp_r5  = H_NTT[loop_index + 5  * offset ];
            H_NTT_tmp_r6  = H_NTT[loop_index + 6  * offset ];
            H_NTT_tmp_r7  = H_NTT[loop_index + 7  * offset ];
            H_NTT_tmp_r8  = H_NTT[loop_index + 8  * offset ];
            H_NTT_tmp_r9  = H_NTT[loop_index + 9  * offset ];
            H_NTT_tmp_r10 = H_NTT[loop_index + 10 * offset ];
            H_NTT_tmp_r11 = H_NTT[loop_index + 11 * offset ];
            H_NTT_tmp_r12 = H_NTT[loop_index + 12 * offset ];
            H_NTT_tmp_r13 = H_NTT[loop_index + 13 * offset ];
            H_NTT_tmp_r14 = H_NTT[loop_index + 14 * offset ];
            H_NTT_tmp_r15 = H_NTT[loop_index + 15 * offset ];
            
            //conver to hex
            H_NTT_HEX_r0  = ZZtohex(H_NTT_tmp_r0);            
            H_NTT_HEX_r1  = ZZtohex(H_NTT_tmp_r1);            
            H_NTT_HEX_r2  = ZZtohex(H_NTT_tmp_r2);            
            H_NTT_HEX_r3  = ZZtohex(H_NTT_tmp_r3);            
            H_NTT_HEX_r4  = ZZtohex(H_NTT_tmp_r4);            
            H_NTT_HEX_r5  = ZZtohex(H_NTT_tmp_r5);            
            H_NTT_HEX_r6  = ZZtohex(H_NTT_tmp_r6);            
            H_NTT_HEX_r7  = ZZtohex(H_NTT_tmp_r7);            
            H_NTT_HEX_r8  = ZZtohex(H_NTT_tmp_r8);            
            H_NTT_HEX_r9  = ZZtohex(H_NTT_tmp_r9);            
            H_NTT_HEX_r10 = ZZtohex(H_NTT_tmp_r10);            
            H_NTT_HEX_r11 = ZZtohex(H_NTT_tmp_r11);            
            H_NTT_HEX_r12 = ZZtohex(H_NTT_tmp_r12);            
            H_NTT_HEX_r13 = ZZtohex(H_NTT_tmp_r13);            
            H_NTT_HEX_r14 = ZZtohex(H_NTT_tmp_r14);            
            H_NTT_HEX_r15 = ZZtohex(H_NTT_tmp_r15);            
 
			H_b0SRAM0 << H_NTT_HEX_r0;
			H_b0SRAM0 << H_NTT_HEX_r1;
			H_b0SRAM0 << "\n";
            
			H_b0SRAM1 << H_NTT_HEX_r2;
			H_b0SRAM1 << H_NTT_HEX_r3;
			H_b0SRAM1 << "\n";
			
			H_b0SRAM2 << H_NTT_HEX_r4;
			H_b0SRAM2 << H_NTT_HEX_r5;
			H_b0SRAM2 << "\n";
			
            
			H_b0SRAM3 << H_NTT_HEX_r6;
			H_b0SRAM3 << H_NTT_HEX_r7;
			H_b0SRAM3 << "\n";

			H_b0SRAM4 << H_NTT_HEX_r8;
			H_b0SRAM4 << H_NTT_HEX_r9;
			H_b0SRAM4 << "\n";
			
			H_b0SRAM5 << H_NTT_HEX_r10;
			H_b0SRAM5 << H_NTT_HEX_r11;
			H_b0SRAM5 << "\n";
			
			
			H_b0SRAM6 << H_NTT_HEX_r12;
			H_b0SRAM6 << H_NTT_HEX_r13;
			H_b0SRAM6 << "\n";
			
			H_b0SRAM7 << H_NTT_HEX_r14;
			H_b0SRAM7 << H_NTT_HEX_r15;
			H_b0SRAM7 << "\n";
			
            loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
    
    counter = 0;
    loop_index = 0;
    //bank 1
    while(counter < counter_iteration){ 
        if((ima[loop_index] == counter )&& (ibn[loop_index] == 1) ){
            H_NTT_tmp_r0  = H_NTT[loop_index              ];
            H_NTT_tmp_r1  = H_NTT[loop_index + 1  * offset ];
            H_NTT_tmp_r2  = H_NTT[loop_index + 2  * offset ];
            H_NTT_tmp_r3  = H_NTT[loop_index + 3  * offset ];
            H_NTT_tmp_r4  = H_NTT[loop_index + 4  * offset ];
            H_NTT_tmp_r5  = H_NTT[loop_index + 5  * offset ];
            H_NTT_tmp_r6  = H_NTT[loop_index + 6  * offset ];
            H_NTT_tmp_r7  = H_NTT[loop_index + 7  * offset ];
            H_NTT_tmp_r8  = H_NTT[loop_index + 8  * offset ];
            H_NTT_tmp_r9  = H_NTT[loop_index + 9  * offset ];
            H_NTT_tmp_r10 = H_NTT[loop_index + 10 * offset ];
            H_NTT_tmp_r11 = H_NTT[loop_index + 11 * offset ];
            H_NTT_tmp_r12 = H_NTT[loop_index + 12 * offset ];
            H_NTT_tmp_r13 = H_NTT[loop_index + 13 * offset ];
            H_NTT_tmp_r14 = H_NTT[loop_index + 14 * offset ];
            H_NTT_tmp_r15 = H_NTT[loop_index + 15 * offset ];
            
            //conver to hex
            H_NTT_HEX_r0  = ZZtohex(H_NTT_tmp_r0);            
            H_NTT_HEX_r1  = ZZtohex(H_NTT_tmp_r1);            
            H_NTT_HEX_r2  = ZZtohex(H_NTT_tmp_r2);            
            H_NTT_HEX_r3  = ZZtohex(H_NTT_tmp_r3);            
            H_NTT_HEX_r4  = ZZtohex(H_NTT_tmp_r4);            
            H_NTT_HEX_r5  = ZZtohex(H_NTT_tmp_r5);            
            H_NTT_HEX_r6  = ZZtohex(H_NTT_tmp_r6);            
            H_NTT_HEX_r7  = ZZtohex(H_NTT_tmp_r7);            
            H_NTT_HEX_r8  = ZZtohex(H_NTT_tmp_r8);            
            H_NTT_HEX_r9  = ZZtohex(H_NTT_tmp_r9);            
            H_NTT_HEX_r10 = ZZtohex(H_NTT_tmp_r10);            
            H_NTT_HEX_r11 = ZZtohex(H_NTT_tmp_r11);            
            H_NTT_HEX_r12 = ZZtohex(H_NTT_tmp_r12);            
            H_NTT_HEX_r13 = ZZtohex(H_NTT_tmp_r13);            
            H_NTT_HEX_r14 = ZZtohex(H_NTT_tmp_r14);            
            H_NTT_HEX_r15 = ZZtohex(H_NTT_tmp_r15);            
 
			H_b1SRAM0 << H_NTT_HEX_r0;
			H_b1SRAM0 << H_NTT_HEX_r1;
			H_b1SRAM0 << "\n";
            
			H_b1SRAM1 << H_NTT_HEX_r2;
			H_b1SRAM1 << H_NTT_HEX_r3;
			H_b1SRAM1 << "\n";
			
			H_b1SRAM2 << H_NTT_HEX_r4;
			H_b1SRAM2 << H_NTT_HEX_r5;
			H_b1SRAM2 << "\n";
			
            
			H_b1SRAM3 << H_NTT_HEX_r6;
			H_b1SRAM3 << H_NTT_HEX_r7;
			H_b1SRAM3 << "\n";

			H_b1SRAM4 << H_NTT_HEX_r8;
			H_b1SRAM4 << H_NTT_HEX_r9;
			H_b1SRAM4 << "\n";
			
			H_b1SRAM5 << H_NTT_HEX_r10;
			H_b1SRAM5 << H_NTT_HEX_r11;
			H_b1SRAM5 << "\n";
			
			
			H_b1SRAM6 << H_NTT_HEX_r12;
			H_b1SRAM6 << H_NTT_HEX_r13;
			H_b1SRAM6 << "\n";
			
			H_b1SRAM7 << H_NTT_HEX_r14;
			H_b1SRAM7 << H_NTT_HEX_r15;
			H_b1SRAM7 << "\n";
            
            loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
}
void SPMB::H_freq_Reconfigure_r16_r2(std::vector<ZZ> H_NTT){
    //1 word  128 bits
	//bank 0
    std::ofstream H_freq_out("./Reconfigure/H_freq_out.txt");
    std::ofstream H_b0SRAM0("./Reconfigure/H_b0SRAM0.txt");
    std::ofstream H_b0SRAM1("./Reconfigure/H_b0SRAM1.txt");
    std::ofstream H_b0SRAM2("./Reconfigure/H_b0SRAM2.txt");
    std::ofstream H_b0SRAM3("./Reconfigure/H_b0SRAM3.txt");
    std::ofstream H_b0SRAM4("./Reconfigure/H_b0SRAM4.txt");
    std::ofstream H_b0SRAM5("./Reconfigure/H_b0SRAM5.txt");
    std::ofstream H_b0SRAM6("./Reconfigure/H_b0SRAM6.txt");
    std::ofstream H_b0SRAM7("./Reconfigure/H_b0SRAM7.txt");
    //bank1                     ROM_Data
    std::ofstream H_b1SRAM0("./Reconfigure/H_b1SRAM0.txt");
    std::ofstream H_b1SRAM1("./Reconfigure/H_b1SRAM1.txt");
    std::ofstream H_b1SRAM2("./Reconfigure/H_b1SRAM2.txt");
    std::ofstream H_b1SRAM3("./Reconfigure/H_b1SRAM3.txt");
    std::ofstream H_b1SRAM4("./Reconfigure/H_b1SRAM4.txt");
    std::ofstream H_b1SRAM5("./Reconfigure/H_b1SRAM5.txt");
    std::ofstream H_b1SRAM6("./Reconfigure/H_b1SRAM6.txt");
    std::ofstream H_b1SRAM7("./Reconfigure/H_b1SRAM7.txt");
	
	//*************************************
	//H_freq data output
	for(int i=0;i<fft_point;i++){
        H_freq_out << H_NTT[i];
        H_freq_out << "\n";
    }
	//*************************************

    std::string H_NTT_HEX_r0;
    std::string H_NTT_HEX_r1;
    std::string H_NTT_HEX_r2;
    std::string H_NTT_HEX_r3;
    std::string H_NTT_HEX_r4;
    std::string H_NTT_HEX_r5;
    std::string H_NTT_HEX_r6;
    std::string H_NTT_HEX_r7;
    std::string H_NTT_HEX_r8;
    std::string H_NTT_HEX_r9;
    std::string H_NTT_HEX_r10;
    std::string H_NTT_HEX_r11;
    std::string H_NTT_HEX_r12;
    std::string H_NTT_HEX_r13;
    std::string H_NTT_HEX_r14;
    std::string H_NTT_HEX_r15;    
	
    ZZ H_NTT_tmp_r0;
    ZZ H_NTT_tmp_r1;
    ZZ H_NTT_tmp_r2;
    ZZ H_NTT_tmp_r3;
    ZZ H_NTT_tmp_r4;
    ZZ H_NTT_tmp_r5;
    ZZ H_NTT_tmp_r6;
    ZZ H_NTT_tmp_r7;
    ZZ H_NTT_tmp_r8;
    ZZ H_NTT_tmp_r9;
    ZZ H_NTT_tmp_r10;
    ZZ H_NTT_tmp_r11;
    ZZ H_NTT_tmp_r12;
    ZZ H_NTT_tmp_r13;
    ZZ H_NTT_tmp_r14;
    ZZ H_NTT_tmp_r15;

    int counter =0;
    int loop_index = 0;
	int INTT_R1_offset;
	int INTT_R2_offset;
	int INTT_R3_offset;
	int INTT_R4_offset;
	int INTT_R5_offset;
	int INTT_R6_offset;
	int INTT_R7_offset;
	int INTT_R8_offset;
	int INTT_R9_offset;
	int INTT_R10_offset;
	int INTT_R11_offset;
	int INTT_R12_offset;
	int INTT_R13_offset;
	int INTT_R14_offset;
	int INTT_R15_offset;
     
	 
	//calculate INTT offset index 
	INTT_R1_offset  = (int) fft_point / 2;
	INTT_R2_offset  = (int) fft_point / 16;
	INTT_R3_offset  = INTT_R1_offset + INTT_R2_offset;
	INTT_R4_offset  = (int)fft_point / 8;
	INTT_R5_offset  = INTT_R4_offset + INTT_R1_offset;
	INTT_R6_offset  = INTT_R4_offset + INTT_R2_offset;
	INTT_R7_offset  = INTT_R4_offset + INTT_R3_offset;
	INTT_R8_offset  = (int)fft_point / 4;
	INTT_R9_offset  = INTT_R8_offset + INTT_R1_offset ;
	INTT_R10_offset = INTT_R8_offset + INTT_R2_offset ;
	INTT_R11_offset = INTT_R8_offset + INTT_R3_offset ;
	INTT_R12_offset = INTT_R8_offset + INTT_R4_offset ;
	INTT_R13_offset = INTT_R8_offset + INTT_R5_offset ;
	INTT_R14_offset = INTT_R8_offset + INTT_R6_offset ;
	INTT_R15_offset = INTT_R8_offset + INTT_R7_offset ;
    
    //bank0
    while(counter < counter_iteration){ 
        if((ima[loop_index] == counter )&& (ibn[loop_index] == 0) ){
			
            H_NTT_tmp_r0  = H_NTT[loop_index ];
            H_NTT_tmp_r1  = H_NTT[loop_index + INTT_R1_offset ];
            H_NTT_tmp_r2  = H_NTT[loop_index + INTT_R2_offset ];
            H_NTT_tmp_r3  = H_NTT[loop_index + INTT_R3_offset ];
            H_NTT_tmp_r4  = H_NTT[loop_index + INTT_R4_offset ];
            H_NTT_tmp_r5  = H_NTT[loop_index + INTT_R5_offset ];
            H_NTT_tmp_r6  = H_NTT[loop_index + INTT_R6_offset ];
            H_NTT_tmp_r7  = H_NTT[loop_index + INTT_R7_offset ];
            H_NTT_tmp_r8  = H_NTT[loop_index + INTT_R8_offset ];
            H_NTT_tmp_r9  = H_NTT[loop_index + INTT_R9_offset ];
            H_NTT_tmp_r10 = H_NTT[loop_index + INTT_R10_offset ];
            H_NTT_tmp_r11 = H_NTT[loop_index + INTT_R11_offset ];
            H_NTT_tmp_r12 = H_NTT[loop_index + INTT_R12_offset ];
            H_NTT_tmp_r13 = H_NTT[loop_index + INTT_R13_offset ];
            H_NTT_tmp_r14 = H_NTT[loop_index + INTT_R14_offset ];
            H_NTT_tmp_r15 = H_NTT[loop_index + INTT_R15_offset ];

            //conver to hex
            H_NTT_HEX_r0  = ZZtohex(H_NTT_tmp_r0);            
            H_NTT_HEX_r1  = ZZtohex(H_NTT_tmp_r1);            
            H_NTT_HEX_r2  = ZZtohex(H_NTT_tmp_r2);            
            H_NTT_HEX_r3  = ZZtohex(H_NTT_tmp_r3);            
            H_NTT_HEX_r4  = ZZtohex(H_NTT_tmp_r4);            
            H_NTT_HEX_r5  = ZZtohex(H_NTT_tmp_r5);            
            H_NTT_HEX_r6  = ZZtohex(H_NTT_tmp_r6);            
            H_NTT_HEX_r7  = ZZtohex(H_NTT_tmp_r7);            
            H_NTT_HEX_r8  = ZZtohex(H_NTT_tmp_r8);            
            H_NTT_HEX_r9  = ZZtohex(H_NTT_tmp_r9);            
            H_NTT_HEX_r10 = ZZtohex(H_NTT_tmp_r10);            
            H_NTT_HEX_r11 = ZZtohex(H_NTT_tmp_r11);            
            H_NTT_HEX_r12 = ZZtohex(H_NTT_tmp_r12);            
            H_NTT_HEX_r13 = ZZtohex(H_NTT_tmp_r13);            
            H_NTT_HEX_r14 = ZZtohex(H_NTT_tmp_r14);            
            H_NTT_HEX_r15 = ZZtohex(H_NTT_tmp_r15);            
 
			H_b0SRAM0 << H_NTT_HEX_r0;
			H_b0SRAM0 << H_NTT_HEX_r1;
			H_b0SRAM0 << "\n";
            
			H_b0SRAM1 << H_NTT_HEX_r2;
			H_b0SRAM1 << H_NTT_HEX_r3;
			H_b0SRAM1 << "\n";
			
			H_b0SRAM2 << H_NTT_HEX_r4;
			H_b0SRAM2 << H_NTT_HEX_r5;
			H_b0SRAM2 << "\n";
			
            
			H_b0SRAM3 << H_NTT_HEX_r6;
			H_b0SRAM3 << H_NTT_HEX_r7;
			H_b0SRAM3 << "\n";

			H_b0SRAM4 << H_NTT_HEX_r8;
			H_b0SRAM4 << H_NTT_HEX_r9;
			H_b0SRAM4 << "\n";
			
			H_b0SRAM5 << H_NTT_HEX_r10;
			H_b0SRAM5 << H_NTT_HEX_r11;
			H_b0SRAM5 << "\n";
			
			
			H_b0SRAM6 << H_NTT_HEX_r12;
			H_b0SRAM6 << H_NTT_HEX_r13;
			H_b0SRAM6 << "\n";
			
			H_b0SRAM7 << H_NTT_HEX_r14;
			H_b0SRAM7 << H_NTT_HEX_r15;
			H_b0SRAM7 << "\n";
			
            loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
    
    counter = 0;
    loop_index = 0;
    //bank 1
    while(counter < counter_iteration){ 
        if((ima[loop_index] == counter )&& (ibn[loop_index] == 1) ){
            H_NTT_tmp_r0  = H_NTT[loop_index              ];
            H_NTT_tmp_r1  = H_NTT[loop_index + INTT_R1_offset  ];
            H_NTT_tmp_r2  = H_NTT[loop_index + INTT_R2_offset  ];
            H_NTT_tmp_r3  = H_NTT[loop_index + INTT_R3_offset  ];
            H_NTT_tmp_r4  = H_NTT[loop_index + INTT_R4_offset  ];
            H_NTT_tmp_r5  = H_NTT[loop_index + INTT_R5_offset  ];
            H_NTT_tmp_r6  = H_NTT[loop_index + INTT_R6_offset  ];
            H_NTT_tmp_r7  = H_NTT[loop_index + INTT_R7_offset  ];
            H_NTT_tmp_r8  = H_NTT[loop_index + INTT_R8_offset  ];
            H_NTT_tmp_r9  = H_NTT[loop_index + INTT_R9_offset  ];
            H_NTT_tmp_r10 = H_NTT[loop_index + INTT_R10_offset ];
            H_NTT_tmp_r11 = H_NTT[loop_index + INTT_R11_offset ];
            H_NTT_tmp_r12 = H_NTT[loop_index + INTT_R12_offset ];
            H_NTT_tmp_r13 = H_NTT[loop_index + INTT_R13_offset ];
            H_NTT_tmp_r14 = H_NTT[loop_index + INTT_R14_offset ];
            H_NTT_tmp_r15 = H_NTT[loop_index + INTT_R15_offset ];
            
            H_NTT_HEX_r0  = ZZtohex(H_NTT_tmp_r0);            
            H_NTT_HEX_r1  = ZZtohex(H_NTT_tmp_r1);            
            H_NTT_HEX_r2  = ZZtohex(H_NTT_tmp_r2);            
            H_NTT_HEX_r3  = ZZtohex(H_NTT_tmp_r3);            
            H_NTT_HEX_r4  = ZZtohex(H_NTT_tmp_r4);            
            H_NTT_HEX_r5  = ZZtohex(H_NTT_tmp_r5);            
            H_NTT_HEX_r6  = ZZtohex(H_NTT_tmp_r6);            
            H_NTT_HEX_r7  = ZZtohex(H_NTT_tmp_r7);            
            H_NTT_HEX_r8  = ZZtohex(H_NTT_tmp_r8);            
            H_NTT_HEX_r9  = ZZtohex(H_NTT_tmp_r9);            
            H_NTT_HEX_r10 = ZZtohex(H_NTT_tmp_r10);            
            H_NTT_HEX_r11 = ZZtohex(H_NTT_tmp_r11);            
            H_NTT_HEX_r12 = ZZtohex(H_NTT_tmp_r12);            
            H_NTT_HEX_r13 = ZZtohex(H_NTT_tmp_r13);            
            H_NTT_HEX_r14 = ZZtohex(H_NTT_tmp_r14);            
            H_NTT_HEX_r15 = ZZtohex(H_NTT_tmp_r15);  
				
			H_b1SRAM0 << H_NTT_HEX_r0;
			H_b1SRAM0 << H_NTT_HEX_r1;
			H_b1SRAM0 << "\n";
            
			H_b1SRAM1 << H_NTT_HEX_r2;
			H_b1SRAM1 << H_NTT_HEX_r3;
			H_b1SRAM1 << "\n";
			
			H_b1SRAM2 << H_NTT_HEX_r4;
			H_b1SRAM2 << H_NTT_HEX_r5;
			H_b1SRAM2 << "\n";
			
            
			H_b1SRAM3 << H_NTT_HEX_r6;
			H_b1SRAM3 << H_NTT_HEX_r7;
			H_b1SRAM3 << "\n";

			H_b1SRAM4 << H_NTT_HEX_r8;
			H_b1SRAM4 << H_NTT_HEX_r9;
			H_b1SRAM4 << "\n";
			
			H_b1SRAM5 << H_NTT_HEX_r10;
			H_b1SRAM5 << H_NTT_HEX_r11;
			H_b1SRAM5 << "\n";
						
			H_b1SRAM6 << H_NTT_HEX_r12;
			H_b1SRAM6 << H_NTT_HEX_r13;
			H_b1SRAM6 << "\n";
			
			H_b1SRAM7 << H_NTT_HEX_r14;
			H_b1SRAM7 << H_NTT_HEX_r15;
			H_b1SRAM7 << "\n";            

            loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
}
void SPMB::H_freq_Reconfigure_r16_r4(std::vector<ZZ> H_NTT){
    //1 word  128 bits
	//bank 0
    std::ofstream H_freq_out("./Reconfigure/H_freq_out.txt");
    std::ofstream H_b0SRAM0("./Reconfigure/H_b0SRAM0.txt");
    std::ofstream H_b0SRAM1("./Reconfigure/H_b0SRAM1.txt");
    std::ofstream H_b0SRAM2("./Reconfigure/H_b0SRAM2.txt");
    std::ofstream H_b0SRAM3("./Reconfigure/H_b0SRAM3.txt");
    std::ofstream H_b0SRAM4("./Reconfigure/H_b0SRAM4.txt");
    std::ofstream H_b0SRAM5("./Reconfigure/H_b0SRAM5.txt");
    std::ofstream H_b0SRAM6("./Reconfigure/H_b0SRAM6.txt");
    std::ofstream H_b0SRAM7("./Reconfigure/H_b0SRAM7.txt");
    //bank1                     ROM_Data
    std::ofstream H_b1SRAM0("./Reconfigure/H_b1SRAM0.txt");
    std::ofstream H_b1SRAM1("./Reconfigure/H_b1SRAM1.txt");
    std::ofstream H_b1SRAM2("./Reconfigure/H_b1SRAM2.txt");
    std::ofstream H_b1SRAM3("./Reconfigure/H_b1SRAM3.txt");
    std::ofstream H_b1SRAM4("./Reconfigure/H_b1SRAM4.txt");
    std::ofstream H_b1SRAM5("./Reconfigure/H_b1SRAM5.txt");
    std::ofstream H_b1SRAM6("./Reconfigure/H_b1SRAM6.txt");
    std::ofstream H_b1SRAM7("./Reconfigure/H_b1SRAM7.txt");
    
    for(int i=0;i<fft_point;i++){
        H_freq_out << H_NTT[i];
        H_freq_out << "\n";
    }	
	
    std::string H_NTT_HEX_r0;
    std::string H_NTT_HEX_r1;
    std::string H_NTT_HEX_r2;
    std::string H_NTT_HEX_r3;
    std::string H_NTT_HEX_r4;
    std::string H_NTT_HEX_r5;
    std::string H_NTT_HEX_r6;
    std::string H_NTT_HEX_r7;
    std::string H_NTT_HEX_r8;
    std::string H_NTT_HEX_r9;
    std::string H_NTT_HEX_r10;
    std::string H_NTT_HEX_r11;
    std::string H_NTT_HEX_r12;
    std::string H_NTT_HEX_r13;
    std::string H_NTT_HEX_r14;
    std::string H_NTT_HEX_r15;    
    
    ZZ H_NTT_tmp_r0;
    ZZ H_NTT_tmp_r1;
    ZZ H_NTT_tmp_r2;
    ZZ H_NTT_tmp_r3;
    ZZ H_NTT_tmp_r4;
    ZZ H_NTT_tmp_r5;
    ZZ H_NTT_tmp_r6;
    ZZ H_NTT_tmp_r7;
    ZZ H_NTT_tmp_r8;
    ZZ H_NTT_tmp_r9;
    ZZ H_NTT_tmp_r10;
    ZZ H_NTT_tmp_r11;
    ZZ H_NTT_tmp_r12;
    ZZ H_NTT_tmp_r13;
    ZZ H_NTT_tmp_r14;
    ZZ H_NTT_tmp_r15;
    
    int counter =0;
    int loop_index = 0;
	int INTT_R1_offset;
	int INTT_R2_offset;
	int INTT_R3_offset;
	int INTT_R4_offset;
	int INTT_R5_offset;
	int INTT_R6_offset;
	int INTT_R7_offset;
	int INTT_R8_offset;
	int INTT_R9_offset;
	int INTT_R10_offset;
	int INTT_R11_offset;
	int INTT_R12_offset;
	int INTT_R13_offset;
	int INTT_R14_offset;
	int INTT_R15_offset;     

	//calculate INTT offset index 
	INTT_R1_offset  = (int) fft_point / 4;
	INTT_R2_offset  = (int) fft_point / 2;
	INTT_R3_offset  = INTT_R1_offset + INTT_R2_offset;
	INTT_R4_offset  = (int)fft_point / 16;
	INTT_R5_offset  = INTT_R4_offset + INTT_R1_offset;
	INTT_R6_offset  = INTT_R4_offset + INTT_R2_offset;
	INTT_R7_offset  = INTT_R4_offset + INTT_R3_offset;
	INTT_R8_offset  = (int)fft_point / 8;
	INTT_R9_offset  = INTT_R8_offset + INTT_R1_offset ;
	INTT_R10_offset = INTT_R8_offset + INTT_R2_offset ;
	INTT_R11_offset = INTT_R8_offset + INTT_R3_offset ;
	INTT_R12_offset = INTT_R8_offset + INTT_R4_offset ;
	INTT_R13_offset = INTT_R8_offset + INTT_R5_offset ;
	INTT_R14_offset = INTT_R8_offset + INTT_R6_offset ;
	INTT_R15_offset = INTT_R8_offset + INTT_R7_offset ;

    
    //bank0
    while(counter < counter_iteration){ 
        if((ima[loop_index] == counter )&& (ibn[loop_index] == 0) ){
            
            H_NTT_tmp_r0  = H_NTT[loop_index              ];
            H_NTT_tmp_r1  = H_NTT[loop_index + INTT_R1_offset ];
            H_NTT_tmp_r2  = H_NTT[loop_index + INTT_R2_offset ];
            H_NTT_tmp_r3  = H_NTT[loop_index + INTT_R3_offset ];
            H_NTT_tmp_r4  = H_NTT[loop_index + INTT_R4_offset ];
            H_NTT_tmp_r5  = H_NTT[loop_index + INTT_R5_offset ];
            H_NTT_tmp_r6  = H_NTT[loop_index + INTT_R6_offset ];
            H_NTT_tmp_r7  = H_NTT[loop_index + INTT_R7_offset ];
            H_NTT_tmp_r8  = H_NTT[loop_index + INTT_R8_offset ];
            H_NTT_tmp_r9  = H_NTT[loop_index + INTT_R9_offset ];
            H_NTT_tmp_r10 = H_NTT[loop_index + INTT_R10_offset];
            H_NTT_tmp_r11 = H_NTT[loop_index + INTT_R11_offset];
            H_NTT_tmp_r12 = H_NTT[loop_index + INTT_R12_offset];
            H_NTT_tmp_r13 = H_NTT[loop_index + INTT_R13_offset];
            H_NTT_tmp_r14 = H_NTT[loop_index + INTT_R14_offset];
            H_NTT_tmp_r15 = H_NTT[loop_index + INTT_R15_offset];
            
            //conver to hex
            H_NTT_HEX_r0  = ZZtohex(H_NTT_tmp_r0);            
            H_NTT_HEX_r1  = ZZtohex(H_NTT_tmp_r1);            
            H_NTT_HEX_r2  = ZZtohex(H_NTT_tmp_r2);            
            H_NTT_HEX_r3  = ZZtohex(H_NTT_tmp_r3);            
            H_NTT_HEX_r4  = ZZtohex(H_NTT_tmp_r4);            
            H_NTT_HEX_r5  = ZZtohex(H_NTT_tmp_r5);            
            H_NTT_HEX_r6  = ZZtohex(H_NTT_tmp_r6);            
            H_NTT_HEX_r7  = ZZtohex(H_NTT_tmp_r7);            
            H_NTT_HEX_r8  = ZZtohex(H_NTT_tmp_r8);            
            H_NTT_HEX_r9  = ZZtohex(H_NTT_tmp_r9);            
            H_NTT_HEX_r10 = ZZtohex(H_NTT_tmp_r10);            
            H_NTT_HEX_r11 = ZZtohex(H_NTT_tmp_r11);            
            H_NTT_HEX_r12 = ZZtohex(H_NTT_tmp_r12);            
            H_NTT_HEX_r13 = ZZtohex(H_NTT_tmp_r13);            
            H_NTT_HEX_r14 = ZZtohex(H_NTT_tmp_r14);            
            H_NTT_HEX_r15 = ZZtohex(H_NTT_tmp_r15);            
 
			H_b0SRAM0 << H_NTT_HEX_r0;
			H_b0SRAM0 << H_NTT_HEX_r1;
			H_b0SRAM0 << "\n";
            
			H_b0SRAM1 << H_NTT_HEX_r2;
			H_b0SRAM1 << H_NTT_HEX_r3;
			H_b0SRAM1 << "\n";
			
			H_b0SRAM2 << H_NTT_HEX_r4;
			H_b0SRAM2 << H_NTT_HEX_r5;
			H_b0SRAM2 << "\n";
			
            
			H_b0SRAM3 << H_NTT_HEX_r6;
			H_b0SRAM3 << H_NTT_HEX_r7;
			H_b0SRAM3 << "\n";

			H_b0SRAM4 << H_NTT_HEX_r8;
			H_b0SRAM4 << H_NTT_HEX_r9;
			H_b0SRAM4 << "\n";
			
			H_b0SRAM5 << H_NTT_HEX_r10;
			H_b0SRAM5 << H_NTT_HEX_r11;
			H_b0SRAM5 << "\n";
			
			
			H_b0SRAM6 << H_NTT_HEX_r12;
			H_b0SRAM6 << H_NTT_HEX_r13;
			H_b0SRAM6 << "\n";
			
			H_b0SRAM7 << H_NTT_HEX_r14;
			H_b0SRAM7 << H_NTT_HEX_r15;
			H_b0SRAM7 << "\n";
			
			
			
            loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
    
    counter = 0;
    loop_index = 0;
    //bank 1
    while(counter < counter_iteration){ 
        if((ima[loop_index] == counter )&& (ibn[loop_index] == 1) ){
            H_NTT_tmp_r0  = H_NTT[loop_index              ];
            H_NTT_tmp_r1  = H_NTT[loop_index + INTT_R1_offset ];
            H_NTT_tmp_r2  = H_NTT[loop_index + INTT_R2_offset ];
            H_NTT_tmp_r3  = H_NTT[loop_index + INTT_R3_offset ];
            H_NTT_tmp_r4  = H_NTT[loop_index + INTT_R4_offset ];
            H_NTT_tmp_r5  = H_NTT[loop_index + INTT_R5_offset ];
            H_NTT_tmp_r6  = H_NTT[loop_index + INTT_R6_offset ];
            H_NTT_tmp_r7  = H_NTT[loop_index + INTT_R7_offset ];
            H_NTT_tmp_r8  = H_NTT[loop_index + INTT_R8_offset ];
            H_NTT_tmp_r9  = H_NTT[loop_index + INTT_R9_offset ];
            H_NTT_tmp_r10 = H_NTT[loop_index + INTT_R10_offset];
            H_NTT_tmp_r11 = H_NTT[loop_index + INTT_R11_offset];
            H_NTT_tmp_r12 = H_NTT[loop_index + INTT_R12_offset];
            H_NTT_tmp_r13 = H_NTT[loop_index + INTT_R13_offset];
            H_NTT_tmp_r14 = H_NTT[loop_index + INTT_R14_offset];
            H_NTT_tmp_r15 = H_NTT[loop_index + INTT_R15_offset];
            
            //conver to hex
            H_NTT_HEX_r0  = ZZtohex(H_NTT_tmp_r0);            
            H_NTT_HEX_r1  = ZZtohex(H_NTT_tmp_r1);            
            H_NTT_HEX_r2  = ZZtohex(H_NTT_tmp_r2);            
            H_NTT_HEX_r3  = ZZtohex(H_NTT_tmp_r3);            
            H_NTT_HEX_r4  = ZZtohex(H_NTT_tmp_r4);            
            H_NTT_HEX_r5  = ZZtohex(H_NTT_tmp_r5);            
            H_NTT_HEX_r6  = ZZtohex(H_NTT_tmp_r6);            
            H_NTT_HEX_r7  = ZZtohex(H_NTT_tmp_r7);            
            H_NTT_HEX_r8  = ZZtohex(H_NTT_tmp_r8);            
            H_NTT_HEX_r9  = ZZtohex(H_NTT_tmp_r9);            
            H_NTT_HEX_r10 = ZZtohex(H_NTT_tmp_r10);            
            H_NTT_HEX_r11 = ZZtohex(H_NTT_tmp_r11);            
            H_NTT_HEX_r12 = ZZtohex(H_NTT_tmp_r12);            
            H_NTT_HEX_r13 = ZZtohex(H_NTT_tmp_r13);            
            H_NTT_HEX_r14 = ZZtohex(H_NTT_tmp_r14);            
            H_NTT_HEX_r15 = ZZtohex(H_NTT_tmp_r15);            
 
			H_b1SRAM0 << H_NTT_HEX_r0;
			H_b1SRAM0 << H_NTT_HEX_r1;
			H_b1SRAM0 << "\n";
            
			H_b1SRAM1 << H_NTT_HEX_r2;
			H_b1SRAM1 << H_NTT_HEX_r3;
			H_b1SRAM1 << "\n";
			
			H_b1SRAM2 << H_NTT_HEX_r4;
			H_b1SRAM2 << H_NTT_HEX_r5;
			H_b1SRAM2 << "\n";
			
            
			H_b1SRAM3 << H_NTT_HEX_r6;
			H_b1SRAM3 << H_NTT_HEX_r7;
			H_b1SRAM3 << "\n";

			H_b1SRAM4 << H_NTT_HEX_r8;
			H_b1SRAM4 << H_NTT_HEX_r9;
			H_b1SRAM4 << "\n";
			
			H_b1SRAM5 << H_NTT_HEX_r10;
			H_b1SRAM5 << H_NTT_HEX_r11;
			H_b1SRAM5 << "\n";
			
			
			H_b1SRAM6 << H_NTT_HEX_r12;
			H_b1SRAM6 << H_NTT_HEX_r13;
			H_b1SRAM6 << "\n";
			
			H_b1SRAM7 << H_NTT_HEX_r14;
			H_b1SRAM7 << H_NTT_HEX_r15;
			H_b1SRAM7 << "\n";
            
            loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
}
void SPMB::H_freq_Reconfigure_r16_r8(std::vector<ZZ> H_NTT){
    //1 word  128 bits
	//bank 0
    std::ofstream H_freq_out("./Reconfigure/H_freq_out.txt");
    std::ofstream H_b0SRAM0("./Reconfigure/H_b0SRAM0.txt");
    std::ofstream H_b0SRAM1("./Reconfigure/H_b0SRAM1.txt");
    std::ofstream H_b0SRAM2("./Reconfigure/H_b0SRAM2.txt");
    std::ofstream H_b0SRAM3("./Reconfigure/H_b0SRAM3.txt");
    std::ofstream H_b0SRAM4("./Reconfigure/H_b0SRAM4.txt");
    std::ofstream H_b0SRAM5("./Reconfigure/H_b0SRAM5.txt");
    std::ofstream H_b0SRAM6("./Reconfigure/H_b0SRAM6.txt");
    std::ofstream H_b0SRAM7("./Reconfigure/H_b0SRAM7.txt");
    //bank1                     ROM_Data
    std::ofstream H_b1SRAM0("./Reconfigure/H_b1SRAM0.txt");
    std::ofstream H_b1SRAM1("./Reconfigure/H_b1SRAM1.txt");
    std::ofstream H_b1SRAM2("./Reconfigure/H_b1SRAM2.txt");
    std::ofstream H_b1SRAM3("./Reconfigure/H_b1SRAM3.txt");
    std::ofstream H_b1SRAM4("./Reconfigure/H_b1SRAM4.txt");
    std::ofstream H_b1SRAM5("./Reconfigure/H_b1SRAM5.txt");
    std::ofstream H_b1SRAM6("./Reconfigure/H_b1SRAM6.txt");
    std::ofstream H_b1SRAM7("./Reconfigure/H_b1SRAM7.txt");

    for(int i=0;i<fft_point;i++){
        H_freq_out << H_NTT[i];
        H_freq_out << "\n";
    }
    
    std::string H_NTT_HEX_r0;
    std::string H_NTT_HEX_r1;
    std::string H_NTT_HEX_r2;
    std::string H_NTT_HEX_r3;
    std::string H_NTT_HEX_r4;
    std::string H_NTT_HEX_r5;
    std::string H_NTT_HEX_r6;
    std::string H_NTT_HEX_r7;
    std::string H_NTT_HEX_r8;
    std::string H_NTT_HEX_r9;
    std::string H_NTT_HEX_r10;
    std::string H_NTT_HEX_r11;
    std::string H_NTT_HEX_r12;
    std::string H_NTT_HEX_r13;
    std::string H_NTT_HEX_r14;
    std::string H_NTT_HEX_r15;    
    
    ZZ H_NTT_tmp_r0;
    ZZ H_NTT_tmp_r1;
    ZZ H_NTT_tmp_r2;
    ZZ H_NTT_tmp_r3;
    ZZ H_NTT_tmp_r4;
    ZZ H_NTT_tmp_r5;
    ZZ H_NTT_tmp_r6;
    ZZ H_NTT_tmp_r7;
    ZZ H_NTT_tmp_r8;
    ZZ H_NTT_tmp_r9;
    ZZ H_NTT_tmp_r10;
    ZZ H_NTT_tmp_r11;
    ZZ H_NTT_tmp_r12;
    ZZ H_NTT_tmp_r13;
    ZZ H_NTT_tmp_r14;
    ZZ H_NTT_tmp_r15;
    
    int counter =0;
    int loop_index = 0;
	int INTT_R1_offset;
	int INTT_R2_offset;
	int INTT_R3_offset;
	int INTT_R4_offset;
	int INTT_R5_offset;
	int INTT_R6_offset;
	int INTT_R7_offset;
	int INTT_R8_offset;
	int INTT_R9_offset;
	int INTT_R10_offset;
	int INTT_R11_offset;
	int INTT_R12_offset;
	int INTT_R13_offset;
	int INTT_R14_offset;
	int INTT_R15_offset;     

	//calculate INTT offset index 
	INTT_R1_offset  = (int) fft_point / 8;
	INTT_R2_offset  = (int) fft_point / 4;
	INTT_R3_offset  = INTT_R1_offset + INTT_R2_offset;
	INTT_R4_offset  = (int)fft_point / 2;
	INTT_R5_offset  = INTT_R4_offset + INTT_R1_offset;
	INTT_R6_offset  = INTT_R4_offset + INTT_R2_offset;
	INTT_R7_offset  = INTT_R4_offset + INTT_R3_offset;
	INTT_R8_offset  = (int)fft_point / 16;
	INTT_R9_offset  = INTT_R8_offset + INTT_R1_offset ;
	INTT_R10_offset = INTT_R8_offset + INTT_R2_offset ;
	INTT_R11_offset = INTT_R8_offset + INTT_R3_offset ;
	INTT_R12_offset = INTT_R8_offset + INTT_R4_offset ;
	INTT_R13_offset = INTT_R8_offset + INTT_R5_offset ;
	INTT_R14_offset = INTT_R8_offset + INTT_R6_offset ;
	INTT_R15_offset = INTT_R8_offset + INTT_R7_offset ;     

    
    //bank0
    while(counter < counter_iteration){ 
        if((ima[loop_index] == counter )&& (ibn[loop_index] == 0) ){
            
            H_NTT_tmp_r0  = H_NTT[loop_index              ];
            H_NTT_tmp_r1  = H_NTT[loop_index + INTT_R1_offset  ];
            H_NTT_tmp_r2  = H_NTT[loop_index + INTT_R2_offset  ];
            H_NTT_tmp_r3  = H_NTT[loop_index + INTT_R3_offset  ];
            H_NTT_tmp_r4  = H_NTT[loop_index + INTT_R4_offset  ];
            H_NTT_tmp_r5  = H_NTT[loop_index + INTT_R5_offset  ];
            H_NTT_tmp_r6  = H_NTT[loop_index + INTT_R6_offset  ];
            H_NTT_tmp_r7  = H_NTT[loop_index + INTT_R7_offset  ];
            H_NTT_tmp_r8  = H_NTT[loop_index + INTT_R8_offset  ];
            H_NTT_tmp_r9  = H_NTT[loop_index + INTT_R9_offset  ];
            H_NTT_tmp_r10 = H_NTT[loop_index + INTT_R10_offset ];
            H_NTT_tmp_r11 = H_NTT[loop_index + INTT_R11_offset ];
            H_NTT_tmp_r12 = H_NTT[loop_index + INTT_R12_offset ];
            H_NTT_tmp_r13 = H_NTT[loop_index + INTT_R13_offset ];
            H_NTT_tmp_r14 = H_NTT[loop_index + INTT_R14_offset ];
            H_NTT_tmp_r15 = H_NTT[loop_index + INTT_R15_offset ];
            
            //conver to hex
            H_NTT_HEX_r0  = ZZtohex(H_NTT_tmp_r0);            
            H_NTT_HEX_r1  = ZZtohex(H_NTT_tmp_r1);            
            H_NTT_HEX_r2  = ZZtohex(H_NTT_tmp_r2);            
            H_NTT_HEX_r3  = ZZtohex(H_NTT_tmp_r3);            
            H_NTT_HEX_r4  = ZZtohex(H_NTT_tmp_r4);            
            H_NTT_HEX_r5  = ZZtohex(H_NTT_tmp_r5);            
            H_NTT_HEX_r6  = ZZtohex(H_NTT_tmp_r6);            
            H_NTT_HEX_r7  = ZZtohex(H_NTT_tmp_r7);            
            H_NTT_HEX_r8  = ZZtohex(H_NTT_tmp_r8);            
            H_NTT_HEX_r9  = ZZtohex(H_NTT_tmp_r9);            
            H_NTT_HEX_r10 = ZZtohex(H_NTT_tmp_r10);            
            H_NTT_HEX_r11 = ZZtohex(H_NTT_tmp_r11);            
            H_NTT_HEX_r12 = ZZtohex(H_NTT_tmp_r12);            
            H_NTT_HEX_r13 = ZZtohex(H_NTT_tmp_r13);            
            H_NTT_HEX_r14 = ZZtohex(H_NTT_tmp_r14);            
            H_NTT_HEX_r15 = ZZtohex(H_NTT_tmp_r15);            
 
			H_b0SRAM0 << H_NTT_HEX_r0;
			H_b0SRAM0 << H_NTT_HEX_r1;
			H_b0SRAM0 << "\n";
            
			H_b0SRAM1 << H_NTT_HEX_r2;
			H_b0SRAM1 << H_NTT_HEX_r3;
			H_b0SRAM1 << "\n";
			
			H_b0SRAM2 << H_NTT_HEX_r4;
			H_b0SRAM2 << H_NTT_HEX_r5;
			H_b0SRAM2 << "\n";
			
            
			H_b0SRAM3 << H_NTT_HEX_r6;
			H_b0SRAM3 << H_NTT_HEX_r7;
			H_b0SRAM3 << "\n";

			H_b0SRAM4 << H_NTT_HEX_r8;
			H_b0SRAM4 << H_NTT_HEX_r9;
			H_b0SRAM4 << "\n";
			
			H_b0SRAM5 << H_NTT_HEX_r10;
			H_b0SRAM5 << H_NTT_HEX_r11;
			H_b0SRAM5 << "\n";
			
			
			H_b0SRAM6 << H_NTT_HEX_r12;
			H_b0SRAM6 << H_NTT_HEX_r13;
			H_b0SRAM6 << "\n";
			
			H_b0SRAM7 << H_NTT_HEX_r14;
			H_b0SRAM7 << H_NTT_HEX_r15;
			H_b0SRAM7 << "\n";
			
			
			
            loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
    
    counter = 0;
    loop_index = 0;
    //bank 1
    while(counter < counter_iteration){ 
        if((ima[loop_index] == counter )&& (ibn[loop_index] == 1) ){
            H_NTT_tmp_r0  = H_NTT[loop_index              ];
            H_NTT_tmp_r1  = H_NTT[loop_index + INTT_R1_offset ];
            H_NTT_tmp_r2  = H_NTT[loop_index + INTT_R2_offset ];
            H_NTT_tmp_r3  = H_NTT[loop_index + INTT_R3_offset ];
            H_NTT_tmp_r4  = H_NTT[loop_index + INTT_R4_offset ];
            H_NTT_tmp_r5  = H_NTT[loop_index + INTT_R5_offset ];
            H_NTT_tmp_r6  = H_NTT[loop_index + INTT_R6_offset ];
            H_NTT_tmp_r7  = H_NTT[loop_index + INTT_R7_offset ];
            H_NTT_tmp_r8  = H_NTT[loop_index + INTT_R8_offset ];
            H_NTT_tmp_r9  = H_NTT[loop_index + INTT_R9_offset ];
            H_NTT_tmp_r10 = H_NTT[loop_index + INTT_R10_offset];
            H_NTT_tmp_r11 = H_NTT[loop_index + INTT_R11_offset];
            H_NTT_tmp_r12 = H_NTT[loop_index + INTT_R12_offset];
            H_NTT_tmp_r13 = H_NTT[loop_index + INTT_R13_offset];
            H_NTT_tmp_r14 = H_NTT[loop_index + INTT_R14_offset];
            H_NTT_tmp_r15 = H_NTT[loop_index + INTT_R15_offset];
            
            //conver to hex
            H_NTT_HEX_r0  = ZZtohex(H_NTT_tmp_r0);            
            H_NTT_HEX_r1  = ZZtohex(H_NTT_tmp_r1);            
            H_NTT_HEX_r2  = ZZtohex(H_NTT_tmp_r2);            
            H_NTT_HEX_r3  = ZZtohex(H_NTT_tmp_r3);            
            H_NTT_HEX_r4  = ZZtohex(H_NTT_tmp_r4);            
            H_NTT_HEX_r5  = ZZtohex(H_NTT_tmp_r5);            
            H_NTT_HEX_r6  = ZZtohex(H_NTT_tmp_r6);            
            H_NTT_HEX_r7  = ZZtohex(H_NTT_tmp_r7);            
            H_NTT_HEX_r8  = ZZtohex(H_NTT_tmp_r8);            
            H_NTT_HEX_r9  = ZZtohex(H_NTT_tmp_r9);            
            H_NTT_HEX_r10 = ZZtohex(H_NTT_tmp_r10);            
            H_NTT_HEX_r11 = ZZtohex(H_NTT_tmp_r11);            
            H_NTT_HEX_r12 = ZZtohex(H_NTT_tmp_r12);            
            H_NTT_HEX_r13 = ZZtohex(H_NTT_tmp_r13);            
            H_NTT_HEX_r14 = ZZtohex(H_NTT_tmp_r14);            
            H_NTT_HEX_r15 = ZZtohex(H_NTT_tmp_r15);            
 
			H_b1SRAM0 << H_NTT_HEX_r0;
			H_b1SRAM0 << H_NTT_HEX_r1;
			H_b1SRAM0 << "\n";
            
			H_b1SRAM1 << H_NTT_HEX_r2;
			H_b1SRAM1 << H_NTT_HEX_r3;
			H_b1SRAM1 << "\n";
			
			H_b1SRAM2 << H_NTT_HEX_r4;
			H_b1SRAM2 << H_NTT_HEX_r5;
			H_b1SRAM2 << "\n";
			
            
			H_b1SRAM3 << H_NTT_HEX_r6;
			H_b1SRAM3 << H_NTT_HEX_r7;
			H_b1SRAM3 << "\n";

			H_b1SRAM4 << H_NTT_HEX_r8;
			H_b1SRAM4 << H_NTT_HEX_r9;
			H_b1SRAM4 << "\n";
			
			H_b1SRAM5 << H_NTT_HEX_r10;
			H_b1SRAM5 << H_NTT_HEX_r11;
			H_b1SRAM5 << "\n";
			
			
			H_b1SRAM6 << H_NTT_HEX_r12;
			H_b1SRAM6 << H_NTT_HEX_r13;
			H_b1SRAM6 << "\n";
			
			H_b1SRAM7 << H_NTT_HEX_r14;
			H_b1SRAM7 << H_NTT_HEX_r15;
			H_b1SRAM7 << "\n";
            
            loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
}

//re-order factor output are store as rom data
//CP_width must be 22~24
//=======================================================
//need to add cp_width 25
//========================================================
void SPMB::re_order_factor_Reconfigure_r16(ZZ m_2_rou){

     std::vector<ZZ> re_order_factor_array;
     re_order_factor_array.resize(fft_point);
     
     for(unsigned long i=0;i<m;i++){
        unsigned long exp;
        ZZ w_2m_tmp_i;
        exp = pow(i,2);
        exp = exp % (2*m);
        PowerMod(w_2m_tmp_i,m_2_rou,exp,cyclotomic_prime);
		//conv(w_2m_tmp_i,"1");
        re_order_factor_array[i] = w_2m_tmp_i;
     }
	//cyclotomic polynomial prime must be 22~24.   
	//re_order_factor_array,while index > m then data is ZERO.
	//because using SPMB Memory addressing,then this reorder factor all index 8 ~ index 15 are zero for radix-16 FFT
    //so we just store reorder factor index 0 ~ index 7 	
    //for example prime bit = 22 ,then total 22 * 8 = 172 bits,so there are 3 ROMs ,word size = 64 bits 
    //ROM0[63:42] = radix0,ROM0[41:20] = radix1, ROM0[19:0] = radix2[21:2]
	//ROM1[63:62] = radix2[1:0], ROM1[61:40] = radix3 ROM1[39:18] = radix4,ROM1[17:0]=radix5[21:4]
	//ROM2[63:60] = radix5[3:0], ROM2[59:38] = radix6,ROM2[37:16] = radix7
	
	
	//bank0
	//128 bits
	std::ofstream reorder_SRAM0("./Reconfigure/Reorder_SRAM0.txt");
	//128 bits
    std::ofstream reorder_SRAM1("./Reconfigure/Reorder_SRAM1.txt");
	//bank1
	//128 bits
    std::ofstream reorder_SRAM2("./Reconfigure/Reorder_SRAM2.txt");
	//128 bits
    std::ofstream reorder_SRAM3("./Reconfigure/Reorder_SRAM3.txt");
   
    std::vector<int> tobit_r0;
    std::vector<int> tobit_r1;
    std::vector<int> tobit_r2;
    std::vector<int> tobit_r3;
    std::vector<int> tobit_r4;
    std::vector<int> tobit_r5;
    std::vector<int> tobit_r6;
    std::vector<int> tobit_r7;
	//128 bits
    std::vector<int> tobit_SRAM0;
	//128 bits
    std::vector<int> tobit_SRAM1;
	//128 bits
    std::vector<int> tobit_SRAM2;
	//128 bits
    std::vector<int> tobit_SRAM3;

    std::string SRAM_128_HEX_string;
    std::string SRAM_64_HEX_string;
    
    tobit_r0.resize(CP_width);
    tobit_r1.resize(CP_width);
    tobit_r2.resize(CP_width);
    tobit_r3.resize(CP_width);
    tobit_r4.resize(CP_width);
    tobit_r5.resize(CP_width);
    tobit_r6.resize(CP_width);
    tobit_r7.resize(CP_width);
    //
	tobit_SRAM0.resize(128);
	tobit_SRAM1.resize(128);
	tobit_SRAM2.resize(128);
	tobit_SRAM3.resize(128);
	
    ZZ re_order_r0_tmp;
    ZZ re_order_r1_tmp;
    ZZ re_order_r2_tmp;
    ZZ re_order_r3_tmp;
    ZZ re_order_r4_tmp;
    ZZ re_order_r5_tmp;
    ZZ re_order_r6_tmp;
    ZZ re_order_r7_tmp;
    
    int counter =0;
    int loop_index = 0;
    int loop_address = 0;

    std::cout << " Reorder_SRAM_Word_Size : " << counter_iteration << "\n";
    //bank0
    while(counter < counter_iteration){ 
        if((ma[loop_index] == counter )&& (bn[loop_index] == 0) ){
            re_order_r0_tmp = re_order_factor_array[loop_index + 0];
            re_order_r1_tmp = re_order_factor_array[loop_index + offset];
            re_order_r2_tmp = re_order_factor_array[loop_index + 2 * offset];
			re_order_r3_tmp = re_order_factor_array[loop_index + 3 * offset];
			re_order_r4_tmp = re_order_factor_array[loop_index + 4 * offset];
			re_order_r5_tmp = re_order_factor_array[loop_index + 5 * offset];
			re_order_r6_tmp = re_order_factor_array[loop_index + 6 * offset];
			re_order_r7_tmp = re_order_factor_array[loop_index + 7 * offset];
			
			//calculate bit
            for(int bit_index=0; bit_index < CP_width ;bit_index++){
                //radix0
				if(re_order_r0_tmp % 2 == 1) tobit_r0[bit_index] = 1;
                else tobit_r0[bit_index] = 0;
				//radix1
				if(re_order_r1_tmp % 2 == 1) tobit_r1[bit_index] = 1;
                else tobit_r1[bit_index] = 0;
				
				if(re_order_r2_tmp % 2 == 1) tobit_r2[bit_index] = 1;
                else tobit_r2[bit_index] = 0;
				
				if(re_order_r3_tmp % 2 == 1) tobit_r3[bit_index] = 1;
                else tobit_r3[bit_index] = 0;
				
				if(re_order_r4_tmp % 2 == 1) tobit_r4[bit_index] = 1;
                else tobit_r4[bit_index] = 0;
				
				if(re_order_r5_tmp % 2 == 1) tobit_r5[bit_index] = 1;
                else tobit_r5[bit_index] = 0;
				
				if(re_order_r6_tmp % 2 == 1) tobit_r6[bit_index] = 1;
                else tobit_r6[bit_index] = 0;
				
				if(re_order_r7_tmp % 2 == 1) tobit_r7[bit_index] = 1;
                else tobit_r7[bit_index] = 0;
				
				re_order_r0_tmp = re_order_r0_tmp >> 1;
				re_order_r1_tmp = re_order_r1_tmp >> 1;
				re_order_r2_tmp = re_order_r2_tmp >> 1;
				re_order_r3_tmp = re_order_r3_tmp >> 1;
				re_order_r4_tmp = re_order_r4_tmp >> 1;
				re_order_r5_tmp = re_order_r5_tmp >> 1;
				re_order_r6_tmp = re_order_r6_tmp >> 1;
				re_order_r7_tmp = re_order_r7_tmp >> 1;
            }
            
            
            for(int g=0;g < CP_width; g++){
				tobit_SRAM0[127-g] = tobit_r0[CP_width-1-g];
            } 
            
            for(int g=0;g < CP_width; g++){
				tobit_SRAM0[127-CP_width-g] = tobit_r1[CP_width-1-g];
            } 
                        
            for(int g=0;g < CP_width; g++){
                tobit_SRAM0[127 - (2 *CP_width)-g] = tobit_r2[CP_width-1-g];
			}
            
          
			for(int g=0;g < CP_width; g++){
                tobit_SRAM0[127 - (3 *CP_width)-g] = tobit_r3[CP_width-1-g] ;
			}
			
			for(int g=0;g < CP_width; g++){
                tobit_SRAM0[127 - (4 *CP_width)-g] = tobit_r4[CP_width-1-g] ;
			}
			
			int SRAM0_Remaining_bits;
		    SRAM0_Remaining_bits = 128 - 5 * CP_width;
			
			for(int g=0;g < CP_width; g++){
                if(SRAM0_Remaining_bits > g)tobit_SRAM0[127 - (5 *CP_width)-g] = tobit_r5[CP_width-1-g] ;
                else tobit_SRAM1[127 - g + SRAM0_Remaining_bits] = tobit_r5[CP_width-1-g] ;
			}
		    
			//BITtohex
			//SRAM 128bit per word
			SRAM_128_HEX_string = BITtohex(tobit_SRAM0,128);
			reorder_SRAM0 << SRAM_128_HEX_string;
			reorder_SRAM0 << "\n";
			
			
			for(int g=0;g < CP_width; g++){
			    tobit_SRAM1[127 - CP_width + SRAM0_Remaining_bits - g] = tobit_r6[CP_width-1-g] ;
			}
			
			for(int g=0;g < CP_width; g++){
			    tobit_SRAM1[127 - (2 * CP_width) + SRAM0_Remaining_bits - g] = tobit_r7[CP_width-1-g] ;
			}
			
		    //BITtohex
			SRAM_128_HEX_string = BITtohex(tobit_SRAM1,128);
			reorder_SRAM1 << SRAM_128_HEX_string;
			reorder_SRAM1 << "\n";
			
			
            loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
    reorder_SRAM0.close();
    reorder_SRAM1.close();

	
    counter = 0;
    loop_index = 0;
    //bank 1
    while(counter < counter_iteration){ 
        if((ma[loop_index] == counter )&& (bn[loop_index] == 1) ){
            re_order_r0_tmp = re_order_factor_array[loop_index + 0];
            re_order_r1_tmp = re_order_factor_array[loop_index + offset];
            re_order_r2_tmp = re_order_factor_array[loop_index + 2 * offset];
			re_order_r3_tmp = re_order_factor_array[loop_index + 3 * offset];
			re_order_r4_tmp = re_order_factor_array[loop_index + 4 * offset];
			re_order_r5_tmp = re_order_factor_array[loop_index + 5 * offset];
			re_order_r6_tmp = re_order_factor_array[loop_index + 6 * offset];
			re_order_r7_tmp = re_order_factor_array[loop_index + 7 * offset];
			
			//calculate bit
            for(int bit_index=0; bit_index < CP_width ;bit_index++){
                //radix0
				if(re_order_r0_tmp % 2 == 1) tobit_r0[bit_index] = 1;
                else tobit_r0[bit_index] = 0;
				//radix1
				if(re_order_r1_tmp % 2 == 1) tobit_r1[bit_index] = 1;
                else tobit_r1[bit_index] = 0;
				
				if(re_order_r2_tmp % 2 == 1) tobit_r2[bit_index] = 1;
                else tobit_r2[bit_index] = 0;
				
				if(re_order_r3_tmp % 2 == 1) tobit_r3[bit_index] = 1;
                else tobit_r3[bit_index] = 0;
				
				if(re_order_r4_tmp % 2 == 1) tobit_r4[bit_index] = 1;
                else tobit_r4[bit_index] = 0;
				
				if(re_order_r5_tmp % 2 == 1) tobit_r5[bit_index] = 1;
                else tobit_r5[bit_index] = 0;
				
				if(re_order_r6_tmp % 2 == 1) tobit_r6[bit_index] = 1;
                else tobit_r6[bit_index] = 0;
				
				if(re_order_r7_tmp % 2 == 1) tobit_r7[bit_index] = 1;
                else tobit_r7[bit_index] = 0;
				
				re_order_r0_tmp = re_order_r0_tmp >> 1;
				re_order_r1_tmp = re_order_r1_tmp >> 1;
				re_order_r2_tmp = re_order_r2_tmp >> 1;
				re_order_r3_tmp = re_order_r3_tmp >> 1;
				re_order_r4_tmp = re_order_r4_tmp >> 1;
				re_order_r5_tmp = re_order_r5_tmp >> 1;
				re_order_r6_tmp = re_order_r6_tmp >> 1;
				re_order_r7_tmp = re_order_r7_tmp >> 1;
            }
            
           
            for(int g=0;g < CP_width; g++){
				tobit_SRAM2[127 - g] = tobit_r0[CP_width-1-g];
            } 
            
            for(int g=0;g < CP_width; g++){
				tobit_SRAM2[127- CP_width - g] = tobit_r1[CP_width-1-g];
            } 
            
            for(int g=0;g < CP_width; g++){
                tobit_SRAM2[127- (2 * CP_width) - g] = tobit_r2[CP_width-1-g];
			}

			for(int g=0;g < CP_width; g++){
                tobit_SRAM2[127- (3 * CP_width) - g] = tobit_r3[CP_width-1-g];
			}
			
			for(int g=0;g < CP_width; g++){
                tobit_SRAM2[127- (4 * CP_width) - g] = tobit_r4[CP_width-1-g];
			}
			
			int SRAM2_Remaining_bits;
		    SRAM2_Remaining_bits = 128 - 5 * CP_width;
			
			for(int g=0;g < CP_width; g++){
                if(SRAM2_Remaining_bits > g)tobit_SRAM2[127- (5 * CP_width) - g] = tobit_r5[CP_width-1-g] ;
                else tobit_SRAM3[127 - g + SRAM2_Remaining_bits] = tobit_r5[CP_width-1-g] ;
			}

			//BITtohex
			//SRAM 128bit per word
			SRAM_128_HEX_string = BITtohex(tobit_SRAM2,128);
			reorder_SRAM2 << SRAM_128_HEX_string;
			reorder_SRAM2 << "\n";

			for(int g=0;g < CP_width; g++){
			    tobit_SRAM3[127 - CP_width + SRAM2_Remaining_bits - g] = tobit_r6[CP_width-1-g] ;
			}
			
			for(int g=0;g < CP_width; g++){
			    tobit_SRAM3[127 - (2 * CP_width ) + SRAM2_Remaining_bits - g] = tobit_r7[CP_width-1-g] ;
			}

		    //BITtohex
			SRAM_128_HEX_string = BITtohex(tobit_SRAM3,128);
			reorder_SRAM3 << SRAM_128_HEX_string;
			reorder_SRAM3 << "\n";			
                 
            loop_index = 0;
            counter = counter + 1;
        }
        else loop_index = loop_index + 1;
    }
    reorder_SRAM2.close();
    reorder_SRAM3.close();
    
    loop_address = 0;
	    
	//inverse order factor 
	std::vector<int> itobit_r0;
	std::vector<int> itobit_r1;
	std::vector<int> itobit_r2;
	std::vector<int> itobit_r3;
	std::vector<int> itobit_r4;
	std::vector<int> itobit_r5;
	std::vector<int> itobit_r6;
	std::vector<int> itobit_r7;
	std::vector<int> itobit_r8;
	std::vector<int> itobit_r9;
	std::vector<int> itobit_r10;
	std::vector<int> itobit_r11;
	std::vector<int> itobit_r12;
	std::vector<int> itobit_r13;
	std::vector<int> itobit_r14;
	std::vector<int> itobit_r15;
	itobit_r0.resize(CP_width);
	itobit_r1.resize(CP_width);
	itobit_r2.resize(CP_width);
	itobit_r3.resize(CP_width);
	itobit_r4.resize(CP_width);
	itobit_r5.resize(CP_width);
	itobit_r6.resize(CP_width);
	itobit_r7.resize(CP_width);
	itobit_r8.resize(CP_width);
	itobit_r9.resize(CP_width);
	itobit_r10.resize(CP_width);
	itobit_r11.resize(CP_width);
	itobit_r12.resize(CP_width);
	itobit_r13.resize(CP_width);
	itobit_r14.resize(CP_width);
	itobit_r15.resize(CP_width);
	//--------------------------
    std::vector<int> tobit_I_SRAM0;
    std::vector<int> tobit_I_SRAM1;
    std::vector<int> tobit_I_SRAM2;
    std::vector<int> tobit_I_SRAM3;
	
	tobit_I_SRAM0.resize(128);
	tobit_I_SRAM1.resize(128);
	tobit_I_SRAM2.resize(128);
	tobit_I_SRAM3.resize(64);
	
	//-----------
	ZZ re_iorder_r0_tmp;
	ZZ re_iorder_r1_tmp;
	ZZ re_iorder_r2_tmp;
	ZZ re_iorder_r3_tmp;
	ZZ re_iorder_r4_tmp;
	ZZ re_iorder_r5_tmp;
	ZZ re_iorder_r6_tmp;
	ZZ re_iorder_r7_tmp;
	ZZ re_iorder_r8_tmp;
	ZZ re_iorder_r9_tmp;
	ZZ re_iorder_r10_tmp;
	ZZ re_iorder_r11_tmp;
	ZZ re_iorder_r12_tmp;
	ZZ re_iorder_r13_tmp;
	ZZ re_iorder_r14_tmp;
	ZZ re_iorder_r15_tmp;
	/*************************************************/
	//m-th cyclotomic polynomial 
	//all polynomial at this thing,the coefficient is zero while the degree > phi(m);
    //and inverse reorder factor index > phi(m) is zero.
	//using SPMB (Final stage type) ,the data index is continuous
	//we can divide inverse re-order factor into two banks ROM.
	//ROM0 bit width= 64 ,ROM1 bit width = 32
	//cyclotomic prime = 22 bits
	//ROM0[63:42] = radix0 ,ROM0[41:20] = radix1,ROM0[19:0] = radix2[21:2]
	//ROM1[63:62] = raidx2[1:0],ROM1[61:40] = radix3 , ROM1[39:18] = radix4 ,ROM1[17:0] = radix5[21:4]
	//ROM2[63:60] = radix5[3:0],ROM2[59:38] = radix6 , ROM2[37:16] = radix7 ,ROM2[15:0] = radix8[21:6]
	//ROM3[63:58] = radix8[5:0],ROM3[57:36] = radix9 , ROM3[35:14] = radix10,ROM3[13:0] = radix11[21:8]
	//ROM4[63:56] = radix11[7:0],ROM4[55:34]= radix12, ROM4[33:12] = radix13,ROM4[11:0] = radix14[21:10]
	//ROM5[31:22] = radix14[9:0],ROM5 = radix15
	//ROM Word size ,original word size is N/radix => 4096/4 = 1024,
	//a_0x^(0)+........+a_phi(m)-1*x^(phi(m)-1)+0*x^phi(m)+...+0^x(m-1)
	//so word size can redution to phi(m)/radix 
    //for example m=1705 phi(m) = 1200 , then ma > 300 data is zero	
	//then we set word size is 512
	/*************************************************/
	
	std::ofstream ireorder_out("./ROM_Data/ireorder_outfile.txt");
    //128 bits
	std::ofstream ireorder_SRAM0("./Reconfigure/IReorder_SRAM0.txt");
	std::ofstream ireorder_SRAM1("./Reconfigure/IReorder_SRAM1.txt");
	std::ofstream ireorder_SRAM2("./Reconfigure/IReorder_SRAM2.txt");
	//64 bits
	std::ofstream ireorder_SRAM3("./Reconfigure/IReorder_SRAM3.txt");
	
    double IReROM_WORD_BIT;
	int    ROM_WORD_SIZE;
    IReROM_WORD_BIT = ceil((double)m / 16);
    IReROM_WORD_BIT = log2(IReROM_WORD_BIT);
    IReROM_WORD_BIT = ceil(IReROM_WORD_BIT);	
    ROM_WORD_SIZE   = exp2(IReROM_WORD_BIT);
	
	for(int i=0;i<fft_point;i++){
		ireorder_out << re_order_factor_array[i] << "\n";
	}
	ireorder_out.close();
	
    std::cout << "Inverse ROM_WORD_SIZE : " << ROM_WORD_SIZE << "\n";
	for(int i = 0;i < ROM_WORD_SIZE; i++){
		 //std::cout << " index: " << i << "\n";
         re_iorder_r0_tmp = re_order_factor_array[radix * i + 0];
         re_iorder_r1_tmp = re_order_factor_array[radix * i + 1];
         re_iorder_r2_tmp = re_order_factor_array[radix * i + 2];
         re_iorder_r3_tmp = re_order_factor_array[radix * i + 3];
         re_iorder_r4_tmp = re_order_factor_array[radix * i + 4];
         re_iorder_r5_tmp = re_order_factor_array[radix * i + 5];
         re_iorder_r6_tmp = re_order_factor_array[radix * i + 6];
         re_iorder_r7_tmp = re_order_factor_array[radix * i + 7];
         re_iorder_r8_tmp = re_order_factor_array[radix * i + 8];
         re_iorder_r9_tmp = re_order_factor_array[radix * i + 9];
         re_iorder_r10_tmp = re_order_factor_array[radix * i + 10];
         re_iorder_r11_tmp = re_order_factor_array[radix * i + 11];
         re_iorder_r12_tmp = re_order_factor_array[radix * i + 12];
         re_iorder_r13_tmp = re_order_factor_array[radix * i + 13];
         re_iorder_r14_tmp = re_order_factor_array[radix * i + 14];
         re_iorder_r15_tmp = re_order_factor_array[radix * i + 15];

        
        for(int bit_index=0; bit_index < CP_width ;bit_index++){
            //radix0
			if(re_iorder_r0_tmp % 2 == 1) itobit_r0[bit_index] = 1;
            else itobit_r0[bit_index] = 0;
			//radix1
			if(re_iorder_r1_tmp % 2 == 1) itobit_r1[bit_index] = 1;
			else itobit_r1[bit_index] = 0;
			//radix2
			if(re_iorder_r2_tmp % 2 == 1) itobit_r2[bit_index] = 1;
			else itobit_r2[bit_index] = 0;
			//radix3
            if(re_iorder_r3_tmp % 2 == 1) itobit_r3[bit_index] = 1;
			else itobit_r3[bit_index] = 0;
			//radix4
			if(re_iorder_r4_tmp % 2 == 1) itobit_r4[bit_index] = 1;
			else itobit_r4[bit_index] = 0;
			//radix5
			if(re_iorder_r5_tmp % 2 == 1) itobit_r5[bit_index] = 1;
			else itobit_r5[bit_index] = 0;
			//radix6
			if(re_iorder_r6_tmp % 2 == 1) itobit_r6[bit_index] = 1;
			else itobit_r6[bit_index] = 0;
			//radix7
			if(re_iorder_r7_tmp % 2 == 1) itobit_r7[bit_index] = 1;
			else itobit_r7[bit_index] = 0;
			//radix8
			if(re_iorder_r8_tmp % 2 == 1) itobit_r8[bit_index] = 1;
			else itobit_r8[bit_index] = 0;
			//radix9
			if(re_iorder_r9_tmp % 2 == 1) itobit_r9[bit_index] = 1;
			else itobit_r9[bit_index] = 0;
			//radix10
			if(re_iorder_r10_tmp % 2 == 1) itobit_r10[bit_index] = 1;
			else itobit_r10[bit_index] = 0;
			//radix11
			if(re_iorder_r11_tmp % 2 == 1) itobit_r11[bit_index] = 1;
			else itobit_r11[bit_index] = 0;
			//radix12
			if(re_iorder_r12_tmp % 2 == 1) itobit_r12[bit_index] = 1;
			else itobit_r12[bit_index] = 0;
			//radix13
			if(re_iorder_r13_tmp % 2 == 1) itobit_r13[bit_index] = 1;
			else itobit_r13[bit_index] = 0;
			//radix14
			if(re_iorder_r14_tmp % 2 == 1) itobit_r14[bit_index] = 1;
			else itobit_r14[bit_index] = 0;
			//radix15
			if(re_iorder_r15_tmp % 2 == 1) itobit_r15[bit_index] = 1;
			else itobit_r15[bit_index] = 0;
			
			re_iorder_r0_tmp  = re_iorder_r0_tmp >> 1;
            re_iorder_r1_tmp  = re_iorder_r1_tmp >> 1;
            re_iorder_r2_tmp  = re_iorder_r2_tmp >> 1;
            re_iorder_r3_tmp  = re_iorder_r3_tmp >> 1;
            re_iorder_r4_tmp  = re_iorder_r4_tmp >> 1;
            re_iorder_r5_tmp  = re_iorder_r5_tmp >> 1;
            re_iorder_r6_tmp  = re_iorder_r6_tmp >> 1;
            re_iorder_r7_tmp  = re_iorder_r7_tmp >> 1;
            re_iorder_r8_tmp  = re_iorder_r8_tmp >> 1;
            re_iorder_r9_tmp  = re_iorder_r9_tmp >> 1;
            re_iorder_r10_tmp = re_iorder_r10_tmp >> 1;
            re_iorder_r11_tmp = re_iorder_r11_tmp >> 1;
            re_iorder_r12_tmp = re_iorder_r12_tmp >> 1;
            re_iorder_r13_tmp = re_iorder_r13_tmp >> 1;
            re_iorder_r14_tmp = re_iorder_r14_tmp >> 1;
            re_iorder_r15_tmp = re_iorder_r15_tmp >> 1;
        }
        //***************************************
		// r0 ~ r4
        for(int g=0;g < CP_width; g++){
			tobit_I_SRAM0[127 - g] = itobit_r0[CP_width-1-g];
        } 
        // 22 * 2 = 44
		for(int g=0;g < CP_width; g++){
            tobit_I_SRAM0[127 -CP_width - g] = itobit_r1[CP_width-1-g];
        } 
        // 22 * 3 = 66
        for(int g=0;g < CP_width; g++){
            tobit_I_SRAM0[127 - (2 * CP_width) - g] = itobit_r2[CP_width-1-g];
        } 
		 // 22 * 4 = 88
        for(int g=0;g < CP_width; g++){
            tobit_I_SRAM0[127 - (3 * CP_width) - g] = itobit_r3[CP_width-1-g];
        } 
		//22 * 5 = 110
		for(int g=0;g < CP_width; g++){
             tobit_I_SRAM0[127 - (4 * CP_width) - g] = itobit_r4[CP_width-1-g];
        } 
		
	    //BITtohex
	    //SRAM 128bit per word
	    SRAM_128_HEX_string = BITtohex(tobit_I_SRAM0,128);
	    ireorder_SRAM0 << SRAM_128_HEX_string;
	    ireorder_SRAM0 << "\n";		
		//**************************************************
		// r5 ~ r9
        for(int g=0;g < CP_width; g++){
			tobit_I_SRAM1[127 - g] = itobit_r5[CP_width-1-g];
        } 
        // 22 * 2 = 44
		for(int g=0;g < CP_width; g++){
            tobit_I_SRAM1[127 -CP_width - g] = itobit_r6[CP_width-1-g];
        } 
        // 22 * 3 = 66
        for(int g=0;g < CP_width; g++){
            tobit_I_SRAM1[127 - (2 * CP_width) - g] = itobit_r7[CP_width-1-g];
        } 
		 // 22 * 4 = 88
        for(int g=0;g < CP_width; g++){
            tobit_I_SRAM1[127 - (3 * CP_width) - g] = itobit_r8[CP_width-1-g];
        } 
		//22 * 5 = 110
		for(int g=0;g < CP_width; g++){
             tobit_I_SRAM1[127 - (4 * CP_width) - g] = itobit_r9[CP_width-1-g];
        } 
		
	    //BITtohex
	    //SRAM 128bit per word
	    SRAM_128_HEX_string = BITtohex(tobit_I_SRAM1,128);
	    ireorder_SRAM1 << SRAM_128_HEX_string;
	    ireorder_SRAM1 << "\n";			
		//**************************************************
		// r10 ~ r14
        for(int g=0;g < CP_width; g++){
			tobit_I_SRAM2[127 - g] = itobit_r10[CP_width-1-g];
        } 
        // 22 * 2 = 44
		for(int g=0;g < CP_width; g++){
            tobit_I_SRAM2[127 -CP_width - g] = itobit_r11[CP_width-1-g];
        } 
        // 22 * 3 = 66
        for(int g=0;g < CP_width; g++){
            tobit_I_SRAM2[127 - (2 * CP_width) - g] = itobit_r12[CP_width-1-g];
        } 
		 // 22 * 4 = 88
        for(int g=0;g < CP_width; g++){
            tobit_I_SRAM2[127 - (3 * CP_width) - g] = itobit_r13[CP_width-1-g];
        } 
		//22 * 5 = 110
		for(int g=0;g < CP_width; g++){
             tobit_I_SRAM2[127 - (4 * CP_width) - g] = itobit_r14[CP_width-1-g];
        } 
	    //BITtohex
	    //SRAM 128bit per word
	    SRAM_128_HEX_string = BITtohex(tobit_I_SRAM2,128);
	    ireorder_SRAM2 << SRAM_128_HEX_string;
	    ireorder_SRAM2 << "\n";
        //*****************************************************
        // r15
        for(int g=0;g < CP_width; g++){
			tobit_I_SRAM3[63 - g] = itobit_r15[CP_width-1-g];
        } 
	    SRAM_64_HEX_string = BITtohex(tobit_I_SRAM3,64);
	    ireorder_SRAM3 << SRAM_64_HEX_string;
	    ireorder_SRAM3 << "\n";
    }
    ireorder_SRAM0.close();
    ireorder_SRAM1.close();
    ireorder_SRAM2.close();
    ireorder_SRAM3.close();
	
}

//Twiddle factor
void SPMB::r16_FFT_TW_ROM(){

    // siang print data -----------
    std::ofstream siang_ROM0("./my_print_data/siang_R16_FFTROM0.txt");
    std::ofstream siang_ROM1("./my_print_data/siang_R16_FFTROM1.txt");
    std::ofstream siang_ROM2("./my_print_data/siang_R16_FFTROM2.txt");
    std::ofstream siang_ROM3("./my_print_data/siang_R16_FFTROM3.txt");
    std::ofstream siang_ROM4("./my_print_data/siang_R16_FFTROM4.txt");
    std::ofstream siang_ROM5("./my_print_data/siang_R16_FFTROM5.txt");
    std::ofstream siang_ROM6("./my_print_data/siang_R16_FFTROM6.txt");
    std::ofstream siang_ROM7("./my_print_data/siang_R16_FFTROM7.txt");
 
    //-----------------------------

    std::ofstream ROM0("./ROM_Data/R16_FFTROM0.txt");  //64bits
    std::ofstream ROM1("./ROM_Data/R16_FFTROM1.txt");  //128bits
    std::ofstream ROM2("./ROM_Data/R16_FFTROM2.txt");  //128bits
    std::ofstream ROM3("./ROM_Data/R16_FFTROM3.txt");  //128bits
    std::ofstream ROM4("./ROM_Data/R16_FFTROM4.txt");  //128bits
    std::ofstream ROM5("./ROM_Data/R16_FFTROM5.txt");  //128bits
    std::ofstream ROM6("./ROM_Data/R16_FFTROM6.txt");  //128bits
    std::ofstream ROM7("./ROM_Data/R16_FFTROM7.txt");  //128bits
    //FFT Prime
    ZZ FFT_Prime;
    ZZ twiddle_65536;
	//radix-0 is zero  
    ZZ R1_tmp;
    ZZ R2_tmp;
	ZZ R3_tmp;
	ZZ R4_tmp;
	ZZ R5_tmp;
	ZZ R6_tmp;
	ZZ R7_tmp;
	ZZ R8_tmp;
	ZZ R9_tmp;
	ZZ R10_tmp;
	ZZ R11_tmp;
	ZZ R12_tmp;
	ZZ R13_tmp;
	ZZ R14_tmp;
	ZZ R15_tmp;
	
    long exp;
    long order;  // order = 65536 / fft_point
    long addr_length; // FFT_point / radix
    conv(FFT_Prime,"18446744069414584321");     //FFT_Prime = 2^64 - 2^32 +1
    conv(twiddle_65536,"14603442835287214144"); //65536-th root of unity 
    order = (long) 65536 / fft_point;
    addr_length  = (long) fft_point / radix;

    std::vector<int> tobit_r1;
    std::vector<int> tobit_r2;
    std::vector<int> tobit_r3;
    std::vector<int> tobit_r4;
    std::vector<int> tobit_r5;
    std::vector<int> tobit_r6;
    std::vector<int> tobit_r7;
    std::vector<int> tobit_r8;
    std::vector<int> tobit_r9;
    std::vector<int> tobit_r10;
    std::vector<int> tobit_r11;
    std::vector<int> tobit_r12;
    std::vector<int> tobit_r13;
    std::vector<int> tobit_r14;
    std::vector<int> tobit_r15;
    tobit_r1.resize(64);
    tobit_r2.resize(64);
    tobit_r3.resize(64);
    tobit_r4.resize(64);
    tobit_r5.resize(64);
    tobit_r6.resize(64);
    tobit_r7.resize(64);
    tobit_r8.resize(64);
    tobit_r9.resize(64);
    tobit_r10.resize(64);
    tobit_r11.resize(64);
    tobit_r12.resize(64);
    tobit_r13.resize(64);
    tobit_r14.resize(64);
    tobit_r15.resize(64);
    
    std::vector<ZZ> tw_Table;  //twiddle factor table 
    tw_Table.resize(fft_point);
    
    exp = 0;
    std::cout << "siang twiddle_65536 = " << twiddle_65536 << ", FFT_Prime = " << FFT_Prime << ", order = " << order << endl;
    for(int i=0;i<fft_point;i++){
       tw_Table[i] = PowerMod(twiddle_65536,exp,FFT_Prime);
       //std::cout << "i = " << i << ", exp = " << exp << ", tw_Table[" << i << "] = " << tw_Table[i] << endl;
       exp = exp + order;
    }
    
    //radix 0
    for(int i=0;i < (fft_point/radix) ;i++){
        R1_tmp  = tw_Table[i];
        R2_tmp  = tw_Table[2*i];
        R3_tmp  = tw_Table[3*i];
        R4_tmp  = tw_Table[4*i];
        R5_tmp  = tw_Table[5*i];
        R6_tmp  = tw_Table[6*i];
        R7_tmp  = tw_Table[7*i];
        R8_tmp  = tw_Table[8*i];
        R9_tmp  = tw_Table[9*i];
        R10_tmp = tw_Table[10*i];
        R11_tmp = tw_Table[11*i];
        R12_tmp = tw_Table[12*i];
        R13_tmp = tw_Table[13*i];
        R14_tmp = tw_Table[14*i];
        R15_tmp = tw_Table[15*i];
        
        //siang twiddle print out
        //i=128
        siang_ROM0 << "exp = "  << i*order    << ", " << R1_tmp;
        siang_ROM1 << "exp = "  << 2*i*order  << ", " << R2_tmp;
        siang_ROM2 << "exp = "  << 4*i*order  << ", " << R4_tmp;
        siang_ROM3 << "exp = "  << 6*i*order  << ", " << R6_tmp;
        siang_ROM4 << "exp = "  << 8*i*order  << ", " << R8_tmp;
        siang_ROM5 << "exp = "  << 10*i*order << ", " << R10_tmp;
        siang_ROM6 << "exp = "  << 12*i*order << ", " << R12_tmp;
        siang_ROM7 << "exp = "  << 14*i*order << ", " << R14_tmp;

        siang_ROM1 << " ,   ";
        siang_ROM2 << " ,   ";
        siang_ROM3 << " ,   ";
        siang_ROM4 << " ,   ";
        siang_ROM5 << " ,   ";
        siang_ROM6 << " ,   ";
        siang_ROM7 << " ,   ";

        siang_ROM1 << "exp = "  << 3*i*order << ", " << R3_tmp;
        siang_ROM2 << "exp = "  << 5*i*order << ", " << R5_tmp;
        siang_ROM3 << "exp = "  << 7*i*order << ", " << R7_tmp;
        siang_ROM4 << "exp = "  << 9*i*order << ", " << R9_tmp;
        siang_ROM5 << "exp  = " << 11*i*order << ", " << R11_tmp;
        siang_ROM6 << "exp  = " << 13*i*order << ", " << R13_tmp;
        siang_ROM7 << "exp  = " << 15*i*order << ", " << R15_tmp;
        
        siang_ROM0 << "\n";
        siang_ROM1 << "\n";
        siang_ROM2 << "\n";
        siang_ROM3 << "\n";
        siang_ROM4 << "\n";
        siang_ROM5 << "\n";
        siang_ROM6 << "\n";
        siang_ROM7 << "\n";
        //------------------------  
		
        for(int bit_index=0; bit_index < 64 ;bit_index++){
            if(R1_tmp % 2 == 1) tobit_r1[bit_index] = 1;
            else tobit_r1[bit_index] = 0;
            
            if(R2_tmp % 2 == 1) tobit_r2[bit_index] = 1;
            else tobit_r2[bit_index] = 0;             

            if(R3_tmp % 2 == 1) tobit_r3[bit_index] = 1;
            else tobit_r3[bit_index] = 0;             
            
			if(R4_tmp % 2 == 1) tobit_r4[bit_index] = 1;
			else tobit_r4[bit_index] = 0;

		    if(R5_tmp % 2 == 1) tobit_r5[bit_index] = 1;
			else tobit_r5[bit_index] = 0;             
            
			if(R6_tmp % 2 == 1) tobit_r6[bit_index] = 1;
			else tobit_r6[bit_index] = 0;             
			
			if(R7_tmp % 2 == 1) tobit_r7[bit_index] = 1;
			else tobit_r7[bit_index] = 0;             
			
			if(R8_tmp % 2 == 1) tobit_r8[bit_index] = 1;
			else tobit_r8[bit_index] = 0;             
			
			if(R9_tmp % 2 == 1) tobit_r9[bit_index] = 1;
			else tobit_r9[bit_index] = 0;             
			
			if(R10_tmp % 2 == 1) tobit_r10[bit_index] = 1;
			else tobit_r10[bit_index] = 0;             
			
			if(R11_tmp % 2 == 1) tobit_r11[bit_index] = 1;
			else tobit_r11[bit_index] = 0;

            if(R12_tmp % 2 == 1) tobit_r12[bit_index] = 1;
			else tobit_r12[bit_index] = 0;             
			
			if(R13_tmp % 2 == 1) tobit_r13[bit_index] = 1;
			else tobit_r13[bit_index] = 0;

            if(R14_tmp % 2 == 1) tobit_r14[bit_index] = 1;			
			else tobit_r14[bit_index] = 0;             
			
			if(R15_tmp % 2 == 1) tobit_r15[bit_index] = 1;
			else tobit_r15[bit_index] = 0;           

			R1_tmp = R1_tmp >> 1;
            R2_tmp = R2_tmp >> 1;
            R3_tmp = R3_tmp >> 1;
            R4_tmp = R4_tmp >> 1;
            R5_tmp = R5_tmp >> 1;
            R6_tmp = R6_tmp >> 1;
            R7_tmp = R7_tmp >> 1;
            R8_tmp = R8_tmp >> 1;
            R9_tmp = R9_tmp >> 1;
            R10_tmp = R10_tmp >> 1;
            R11_tmp = R11_tmp >> 1;
            R12_tmp = R12_tmp >> 1;
            R13_tmp = R13_tmp >> 1;
            R14_tmp = R14_tmp >> 1;
            R15_tmp = R15_tmp >> 1;
        }
        
		for(int g=0;g < 64; g++){
            ROM0 << tobit_r1[63-g];
            ROM1 << tobit_r2[63-g];
            ROM2 << tobit_r4[63-g];
            ROM3 << tobit_r6[63-g];
            ROM4 << tobit_r8[63-g];
            ROM5 << tobit_r10[63-g];
            ROM6 << tobit_r12[63-g];
            ROM7 << tobit_r14[63-g];
        }
        ROM1 << ", ";
        ROM2 << ", ";
        ROM3 << ", ";
        ROM4 << ", ";
        ROM5 << ", ";
        ROM6 << ", ";
        ROM7 << ", ";
		for(int g=0;g < 64; g++){
            ROM1 << tobit_r3[63-g];
            ROM2 << tobit_r5[63-g];
            ROM3 << tobit_r7[63-g];
            ROM4 << tobit_r9[63-g];
            ROM5 << tobit_r11[63-g];
            ROM6 << tobit_r13[63-g];
            ROM7 << tobit_r15[63-g];
        }
		ROM0 << "\n";
        ROM1 << "\n";
        ROM2 << "\n";
        ROM3 << "\n";
        ROM4 << "\n";
        ROM5 << "\n";
        ROM6 << "\n";
        ROM7 << "\n";
	}
	
	ROM0.close();
	ROM1.close();
	ROM2.close();
	ROM3.close();
	ROM4.close();
	ROM5.close();
	ROM6.close();
	ROM7.close();

    //---------------
    siang_ROM0.close();
	siang_ROM1.close();
	siang_ROM2.close();
	siang_ROM3.close();
	siang_ROM4.close();
	siang_ROM5.close();
	siang_ROM6.close();
	siang_ROM7.close();

    //----------------
}
void SPMB::r16_IFFT_TW_ROM(){
    
    std::ofstream  ROM0("./ROM_Data/R16_IFFTROM0.txt");
    std::ofstream  ROM1("./ROM_Data/R16_IFFTROM1.txt");
    std::ofstream  ROM2("./ROM_Data/R16_IFFTROM2.txt");
    std::ofstream  ROM3("./ROM_Data/R16_IFFTROM3.txt");                       
    std::ofstream  ROM4("./ROM_Data/R16_IFFTROM4.txt");
    std::ofstream  ROM5("./ROM_Data/R16_IFFTROM5.txt");
    std::ofstream  ROM6("./ROM_Data/R16_IFFTROM6.txt");
    std::ofstream  ROM7("./ROM_Data/R16_IFFTROM7.txt");
    //FFT Prime
    ZZ FFT_Prime;
    ZZ twiddle_65536;
    ZZ R1_tmp;
    ZZ R2_tmp;
    ZZ R3_tmp;
    ZZ R4_tmp;
    ZZ R5_tmp;
    ZZ R6_tmp;
    ZZ R7_tmp;
    ZZ R8_tmp;
    ZZ R9_tmp;
    ZZ R10_tmp;
    ZZ R11_tmp;
    ZZ R12_tmp;
    ZZ R13_tmp;
    ZZ R14_tmp;
    ZZ R15_tmp;
    long exp;
    long order;  // order = 65536 / fft_point
    long addr_length; // FFT_point / radix
    conv(FFT_Prime,"18446744069414584321");     //FFT_Prime = 2^64 - 2^32 +1
    conv(twiddle_65536,"14603442835287214144"); //65536-th root of unity 
    order = (long) 65536 / fft_point;
    addr_length  = (long) fft_point / radix;

    std::vector<int> tobit_r1;
    std::vector<int> tobit_r2;
    std::vector<int> tobit_r3;
    std::vector<int> tobit_r4;
    std::vector<int> tobit_r5;
    std::vector<int> tobit_r6;
    std::vector<int> tobit_r7;
    std::vector<int> tobit_r8;
    std::vector<int> tobit_r9;
    std::vector<int> tobit_r10;
    std::vector<int> tobit_r11;
    std::vector<int> tobit_r12;
    std::vector<int> tobit_r13;
    std::vector<int> tobit_r14;
    std::vector<int> tobit_r15;
    tobit_r1.resize(64);
    tobit_r2.resize(64);
    tobit_r3.resize(64);
    tobit_r4.resize(64);
    tobit_r5.resize(64);
    tobit_r6.resize(64);
    tobit_r7.resize(64);
    tobit_r8.resize(64);
    tobit_r9.resize(64);
    tobit_r10.resize(64);
    tobit_r11.resize(64);
    tobit_r12.resize(64);
    tobit_r13.resize(64);
    tobit_r14.resize(64);
    tobit_r15.resize(64);
    
    std::vector<ZZ> tw_Table;  //twiddle factor table 
    tw_Table.resize(fft_point);
    
	std::ofstream IFFTROM("./ROM_Data/ITW_ROM.txt");
	
    exp = 0;
    for(int i=0;i<fft_point;i++){
       tw_Table[i] = PowerMod(twiddle_65536,exp,FFT_Prime);
       exp = exp + order;
	   IFFTROM << tw_Table[i] << "\n";
    }
	IFFTROM.close();
	
	
    for(int i=0;i<(fft_point/radix);i++){
        if(i==0){    
            R1_tmp  = tw_Table[0];
            R2_tmp  = tw_Table[0];
            R3_tmp  = tw_Table[0];
            R4_tmp  = tw_Table[0];
            R5_tmp  = tw_Table[0];
            R6_tmp  = tw_Table[0];
            R7_tmp  = tw_Table[0];
            R8_tmp  = tw_Table[0];
            R9_tmp  = tw_Table[0];
            R10_tmp = tw_Table[0];
            R11_tmp = tw_Table[0];
            R12_tmp = tw_Table[0];
            R13_tmp = tw_Table[0];
            R14_tmp = tw_Table[0];
            R15_tmp = tw_Table[0];
        }
        else {
            R1_tmp  = tw_Table[fft_point - i];
            R2_tmp  = tw_Table[fft_point - 2*i];
            R3_tmp  = tw_Table[fft_point - 3*i];
            R4_tmp  = tw_Table[fft_point - 4*i];
            R5_tmp  = tw_Table[fft_point - 5*i];
            R6_tmp  = tw_Table[fft_point - 6*i];
            R7_tmp  = tw_Table[fft_point - 7*i];
            R8_tmp  = tw_Table[fft_point - 8*i];
            R9_tmp  = tw_Table[fft_point - 9*i];
            R10_tmp = tw_Table[fft_point - 10*i];
            R11_tmp = tw_Table[fft_point - 11*i];
            R12_tmp = tw_Table[fft_point - 12*i];
            R13_tmp = tw_Table[fft_point - 13*i];
            R14_tmp = tw_Table[fft_point - 14*i];
            R15_tmp = tw_Table[fft_point - 15*i];
        }
        for(int bit_index=0; bit_index < 64 ;bit_index++){
            if(R1_tmp % 2 == 1) tobit_r1[bit_index] = 1;
            else tobit_r1[bit_index] = 0;
            
            if(R2_tmp % 2 == 1) tobit_r2[bit_index] = 1;
            else tobit_r2[bit_index] = 0;             

            if(R3_tmp % 2 == 1) tobit_r3[bit_index] = 1;
            else tobit_r3[bit_index] = 0;             
            
			if(R4_tmp % 2 == 1) tobit_r4[bit_index] = 1;
            else tobit_r4[bit_index] = 0; 
			
			if(R5_tmp % 2 == 1) tobit_r5[bit_index] = 1;
            else tobit_r5[bit_index] = 0; 
			
		    if(R6_tmp % 2 == 1) tobit_r6[bit_index] = 1;
            else tobit_r6[bit_index] = 0; 
			
			if(R7_tmp % 2 == 1) tobit_r7[bit_index] = 1;
			else tobit_r7[bit_index] = 0; 
            
			if(R8_tmp % 2 == 1) tobit_r8[bit_index] = 1;
			else tobit_r8[bit_index] = 0; 
			
			if(R9_tmp % 2 == 1) tobit_r9[bit_index] = 1;
			else tobit_r9[bit_index] = 0; 
			
			if(R10_tmp % 2 == 1) tobit_r10[bit_index] = 1;
			else tobit_r10[bit_index] = 0; 
			
			if(R11_tmp % 2 == 1) tobit_r11[bit_index] = 1;
			else tobit_r11[bit_index] = 0; 
			
			if(R12_tmp % 2 == 1) tobit_r12[bit_index] = 1;
			else tobit_r12[bit_index] = 0; 
			
			if(R13_tmp % 2 == 1) tobit_r13[bit_index] = 1;
			else tobit_r13[bit_index] = 0; 
			
			if(R14_tmp % 2 == 1) tobit_r14[bit_index] = 1;
			else tobit_r14[bit_index] = 0; 
			
			if(R15_tmp % 2 == 1) tobit_r15[bit_index] = 1;
			else tobit_r15[bit_index] = 0; 
			
			R1_tmp  = R1_tmp >> 1;
            R2_tmp  = R2_tmp >> 1;
            R3_tmp  = R3_tmp >> 1;
            R4_tmp  = R4_tmp >> 1;
            R5_tmp  = R5_tmp >> 1;
            R6_tmp  = R6_tmp >> 1;
            R7_tmp  = R7_tmp >> 1;
            R8_tmp  = R8_tmp >> 1;
            R9_tmp  = R9_tmp >> 1;
            R10_tmp = R10_tmp >> 1;
            R11_tmp = R11_tmp >> 1;
            R12_tmp = R12_tmp >> 1;
            R13_tmp = R13_tmp >> 1;
            R14_tmp = R14_tmp >> 1;
            R15_tmp = R15_tmp >> 1;
        }
        
        for(int g=0;g < 64; g++){
            ROM0 << tobit_r1[63-g];
            ROM1 << tobit_r2[63-g];
            ROM2 << tobit_r4[63-g];
            ROM3 << tobit_r6[63-g];
            ROM4 << tobit_r8[63-g];
            ROM5 << tobit_r10[63-g];
            ROM6 << tobit_r12[63-g];
            ROM7 << tobit_r14[63-g];
        }
        for(int g=0;g < 64; g++){
           ROM1 << tobit_r3[63-g];
		   ROM2 << tobit_r5[63-g];
		   ROM3 << tobit_r7[63-g];
		   ROM4 << tobit_r9[63-g];
		   ROM5 << tobit_r11[63-g];
		   ROM6 << tobit_r13[63-g];
		   ROM7 << tobit_r15[63-g];
        }
		ROM0 << "\n";
        ROM1 << "\n";
        ROM2 << "\n";
        ROM3 << "\n";
        ROM4 << "\n";
        ROM5 << "\n";
        ROM6 << "\n";
        ROM7 << "\n";
    }
	ROM0.close();
	ROM1.close();
	ROM2.close();
	ROM3.close();
	ROM4.close();
	ROM5.close();
	ROM6.close();
	ROM7.close();
}
//Twiddle factor
void SPMB::r16_FFT_TW_ROM_Reconfig(){
    std::ofstream ROM0("./ROM_Data/R16_FFTROM0.txt");  //64bits
    std::ofstream ROM1("./ROM_Data/R16_FFTROM1.txt");  //128bits
    std::ofstream ROM2("./ROM_Data/R16_FFTROM2.txt");  //128bits
    std::ofstream ROM3("./ROM_Data/R16_FFTROM3.txt");  //128bits
    std::ofstream ROM4("./ROM_Data/R16_FFTROM4.txt");  //128bits
    std::ofstream ROM5("./ROM_Data/R16_FFTROM5.txt");  //128bits
    std::ofstream ROM6("./ROM_Data/R16_FFTROM6.txt");  //128bits
    std::ofstream ROM7("./ROM_Data/R16_FFTROM7.txt");  //128bits
    //FFT Prime
    ZZ FFT_Prime;
    ZZ twiddle_65536;
	//radix-0 is zero  
    ZZ R1_tmp;
    ZZ R2_tmp;
	ZZ R3_tmp;
	ZZ R4_tmp;
	ZZ R5_tmp;
	ZZ R6_tmp;
	ZZ R7_tmp;
	ZZ R8_tmp;
	ZZ R9_tmp;
	ZZ R10_tmp;
	ZZ R11_tmp;
	ZZ R12_tmp;
	ZZ R13_tmp;
	ZZ R14_tmp;
	ZZ R15_tmp;
	
    long exp;
    conv(FFT_Prime,"18446744069414584321");     //FFT_Prime = 2^64 - 2^32 +1
    conv(twiddle_65536,"14603442835287214144"); //65536-th root of unity 
   
    std::vector<int> tobit_r1;
    std::vector<int> tobit_r2;
    std::vector<int> tobit_r3;
    std::vector<int> tobit_r4;
    std::vector<int> tobit_r5;
    std::vector<int> tobit_r6;
    std::vector<int> tobit_r7;
    std::vector<int> tobit_r8;
    std::vector<int> tobit_r9;
    std::vector<int> tobit_r10;
    std::vector<int> tobit_r11;
    std::vector<int> tobit_r12;
    std::vector<int> tobit_r13;
    std::vector<int> tobit_r14;
    std::vector<int> tobit_r15;
    tobit_r1.resize(64);
    tobit_r2.resize(64);
    tobit_r3.resize(64);
    tobit_r4.resize(64);
    tobit_r5.resize(64);
    tobit_r6.resize(64);
    tobit_r7.resize(64);
    tobit_r8.resize(64);
    tobit_r9.resize(64);
    tobit_r10.resize(64);
    tobit_r11.resize(64);
    tobit_r12.resize(64);
    tobit_r13.resize(64);
    tobit_r14.resize(64);
    tobit_r15.resize(64);
    
    std::vector<ZZ> tw_Table;  //twiddle factor table 
    tw_Table.resize(65536);
    
    exp = 0;
    for(int i=0;i<65536;i++){
       tw_Table[i] = PowerMod(twiddle_65536,exp,FFT_Prime);
       exp = exp + 1;
    }
    
    //radix 0
    for(int i=0;i < 4096 ;i++){
        R1_tmp  = tw_Table[i];
        R2_tmp  = tw_Table[2*i];
        R3_tmp  = tw_Table[3*i];
        R4_tmp  = tw_Table[4*i];
        R5_tmp  = tw_Table[5*i];
        R6_tmp  = tw_Table[6*i];
        R7_tmp  = tw_Table[7*i];
        R8_tmp  = tw_Table[8*i];
        R9_tmp  = tw_Table[9*i];
        R10_tmp = tw_Table[10*i];
        R11_tmp = tw_Table[11*i];
        R12_tmp = tw_Table[12*i];
        R13_tmp = tw_Table[13*i];
        R14_tmp = tw_Table[14*i];
        R15_tmp = tw_Table[15*i];
		
        for(int bit_index=0; bit_index < 64 ;bit_index++){
            if(R1_tmp % 2 == 1) tobit_r1[bit_index] = 1;
            else tobit_r1[bit_index] = 0;
            
            if(R2_tmp % 2 == 1) tobit_r2[bit_index] = 1;
            else tobit_r2[bit_index] = 0;             

            if(R3_tmp % 2 == 1) tobit_r3[bit_index] = 1;
            else tobit_r3[bit_index] = 0;             
            
			if(R4_tmp % 2 == 1) tobit_r4[bit_index] = 1;
			else tobit_r4[bit_index] = 0;

		    if(R5_tmp % 2 == 1) tobit_r5[bit_index] = 1;
			else tobit_r5[bit_index] = 0;             
            
			if(R6_tmp % 2 == 1) tobit_r6[bit_index] = 1;
			else tobit_r6[bit_index] = 0;             
			
			if(R7_tmp % 2 == 1) tobit_r7[bit_index] = 1;
			else tobit_r7[bit_index] = 0;             
			
			if(R8_tmp % 2 == 1) tobit_r8[bit_index] = 1;
			else tobit_r8[bit_index] = 0;             
			
			if(R9_tmp % 2 == 1) tobit_r9[bit_index] = 1;
			else tobit_r9[bit_index] = 0;             
			
			if(R10_tmp % 2 == 1) tobit_r10[bit_index] = 1;
			else tobit_r10[bit_index] = 0;             
			
			if(R11_tmp % 2 == 1) tobit_r11[bit_index] = 1;
			else tobit_r11[bit_index] = 0;

            if(R12_tmp % 2 == 1) tobit_r12[bit_index] = 1;
			else tobit_r12[bit_index] = 0;             
			
			if(R13_tmp % 2 == 1) tobit_r13[bit_index] = 1;
			else tobit_r13[bit_index] = 0;

            if(R14_tmp % 2 == 1) tobit_r14[bit_index] = 1;			
			else tobit_r14[bit_index] = 0;             
			
			if(R15_tmp % 2 == 1) tobit_r15[bit_index] = 1;
			else tobit_r15[bit_index] = 0;             
			
			R1_tmp = R1_tmp >> 1;
            R2_tmp = R2_tmp >> 1;
            R3_tmp = R3_tmp >> 1;
            R4_tmp = R4_tmp >> 1;
            R5_tmp = R5_tmp >> 1;
            R6_tmp = R6_tmp >> 1;
            R7_tmp = R7_tmp >> 1;
            R8_tmp = R8_tmp >> 1;
            R9_tmp = R9_tmp >> 1;
            R10_tmp = R10_tmp >> 1;
            R11_tmp = R11_tmp >> 1;
            R12_tmp = R12_tmp >> 1;
            R13_tmp = R13_tmp >> 1;
            R14_tmp = R14_tmp >> 1;
            R15_tmp = R15_tmp >> 1;
        }
        
		for(int g=0;g < 64; g++){
            ROM0 << tobit_r1[63-g];
            ROM1 << tobit_r2[63-g];
            ROM2 << tobit_r4[63-g];
            ROM3 << tobit_r6[63-g];
            ROM4 << tobit_r8[63-g];
            ROM5 << tobit_r10[63-g];
            ROM6 << tobit_r12[63-g];
            ROM7 << tobit_r14[63-g];
        }
		for(int g=0;g < 64; g++){
            ROM1 << tobit_r3[63-g];
            ROM2 << tobit_r5[63-g];
            ROM3 << tobit_r7[63-g];
            ROM4 << tobit_r9[63-g];
            ROM5 << tobit_r11[63-g];
            ROM6 << tobit_r13[63-g];
            ROM7 << tobit_r15[63-g];
        }
		ROM0 << "\n";
        ROM1 << "\n";
        ROM2 << "\n";
        ROM3 << "\n";
        ROM4 << "\n";
        ROM5 << "\n";
        ROM6 << "\n";
        ROM7 << "\n";
	}
	
	ROM0.close();
	ROM1.close();
	ROM2.close();
	ROM3.close();
	ROM4.close();
	ROM5.close();
	ROM6.close();
	ROM7.close();
}
void SPMB::r16_IFFT_TW_ROM_Reconfig(){
    
    std::ofstream  ROM0("./ROM_Data/R16_IFFTROM0.txt");
    std::ofstream  ROM1("./ROM_Data/R16_IFFTROM1.txt");
    std::ofstream  ROM2("./ROM_Data/R16_IFFTROM2.txt");
    std::ofstream  ROM3("./ROM_Data/R16_IFFTROM3.txt");                       
    std::ofstream  ROM4("./ROM_Data/R16_IFFTROM4.txt");
    std::ofstream  ROM5("./ROM_Data/R16_IFFTROM5.txt");
    std::ofstream  ROM6("./ROM_Data/R16_IFFTROM6.txt");
    std::ofstream  ROM7("./ROM_Data/R16_IFFTROM7.txt");
    //FFT Prime
    ZZ FFT_Prime;
    ZZ twiddle_65536;
    ZZ R1_tmp;
    ZZ R2_tmp;
    ZZ R3_tmp;
    ZZ R4_tmp;
    ZZ R5_tmp;
    ZZ R6_tmp;
    ZZ R7_tmp;
    ZZ R8_tmp;
    ZZ R9_tmp;
    ZZ R10_tmp;
    ZZ R11_tmp;
    ZZ R12_tmp;
    ZZ R13_tmp;
    ZZ R14_tmp;
    ZZ R15_tmp;
    long exp;
    conv(FFT_Prime,"18446744069414584321");     //FFT_Prime = 2^64 - 2^32 +1
    conv(twiddle_65536,"14603442835287214144"); //65536-th root of unity 

    std::vector<int> tobit_r1;
    std::vector<int> tobit_r2;
    std::vector<int> tobit_r3;
    std::vector<int> tobit_r4;
    std::vector<int> tobit_r5;
    std::vector<int> tobit_r6;
    std::vector<int> tobit_r7;
    std::vector<int> tobit_r8;
    std::vector<int> tobit_r9;
    std::vector<int> tobit_r10;
    std::vector<int> tobit_r11;
    std::vector<int> tobit_r12;
    std::vector<int> tobit_r13;
    std::vector<int> tobit_r14;
    std::vector<int> tobit_r15;
    tobit_r1.resize(64);
    tobit_r2.resize(64);
    tobit_r3.resize(64);
    tobit_r4.resize(64);
    tobit_r5.resize(64);
    tobit_r6.resize(64);
    tobit_r7.resize(64);
    tobit_r8.resize(64);
    tobit_r9.resize(64);
    tobit_r10.resize(64);
    tobit_r11.resize(64);
    tobit_r12.resize(64);
    tobit_r13.resize(64);
    tobit_r14.resize(64);
    tobit_r15.resize(64);
    
    std::vector<ZZ> tw_Table;  //twiddle factor table 
    tw_Table.resize(65536);
    
	std::ofstream IFFTROM("./ROM_Data/ITW_ROM.txt");
	
    exp = 0;
    for(int i=0;i<65536;i++){
       tw_Table[i] = PowerMod(twiddle_65536,exp,FFT_Prime);
       exp = exp + 1;
	   IFFTROM << tw_Table[i] << "\n";
    }
	IFFTROM.close();
	
	
    for(int i=0;i < (4096) ;i++){
        if(i==0){    
            R1_tmp  = tw_Table[0];
            R2_tmp  = tw_Table[0];
            R3_tmp  = tw_Table[0];
            R4_tmp  = tw_Table[0];
            R5_tmp  = tw_Table[0];
            R6_tmp  = tw_Table[0];
            R7_tmp  = tw_Table[0];
            R8_tmp  = tw_Table[0];
            R9_tmp  = tw_Table[0];
            R10_tmp = tw_Table[0];
            R11_tmp = tw_Table[0];
            R12_tmp = tw_Table[0];
            R13_tmp = tw_Table[0];
            R14_tmp = tw_Table[0];
            R15_tmp = tw_Table[0];
        }
        else {
            R1_tmp  = tw_Table[65536 - i];
            R2_tmp  = tw_Table[65536 - 2*i];
            R3_tmp  = tw_Table[65536 - 3*i];
            R4_tmp  = tw_Table[65536 - 4*i];
            R5_tmp  = tw_Table[65536 - 5*i];
            R6_tmp  = tw_Table[65536 - 6*i];
            R7_tmp  = tw_Table[65536 - 7*i];
            R8_tmp  = tw_Table[65536 - 8*i];
            R9_tmp  = tw_Table[65536 - 9*i];
            R10_tmp = tw_Table[65536 - 10*i];
            R11_tmp = tw_Table[65536 - 11*i];
            R12_tmp = tw_Table[65536 - 12*i];
            R13_tmp = tw_Table[65536 - 13*i];
            R14_tmp = tw_Table[65536 - 14*i];
            R15_tmp = tw_Table[65536 - 15*i];
        }
        for(int bit_index=0; bit_index < 64 ;bit_index++){
            if(R1_tmp % 2 == 1) tobit_r1[bit_index] = 1;
            else tobit_r1[bit_index] = 0;
            
            if(R2_tmp % 2 == 1) tobit_r2[bit_index] = 1;
            else tobit_r2[bit_index] = 0;             

            if(R3_tmp % 2 == 1) tobit_r3[bit_index] = 1;
            else tobit_r3[bit_index] = 0;             
            
			if(R4_tmp % 2 == 1) tobit_r4[bit_index] = 1;
            else tobit_r4[bit_index] = 0; 
			
			if(R5_tmp % 2 == 1) tobit_r5[bit_index] = 1;
            else tobit_r5[bit_index] = 0; 
			
		    if(R6_tmp % 2 == 1) tobit_r6[bit_index] = 1;
            else tobit_r6[bit_index] = 0; 
			
			if(R7_tmp % 2 == 1) tobit_r7[bit_index] = 1;
			else tobit_r7[bit_index] = 0; 
            
			if(R8_tmp % 2 == 1) tobit_r8[bit_index] = 1;
			else tobit_r8[bit_index] = 0; 
			
			if(R9_tmp % 2 == 1) tobit_r9[bit_index] = 1;
			else tobit_r9[bit_index] = 0; 
			
			if(R10_tmp % 2 == 1) tobit_r10[bit_index] = 1;
			else tobit_r10[bit_index] = 0; 
			
			if(R11_tmp % 2 == 1) tobit_r11[bit_index] = 1;
			else tobit_r11[bit_index] = 0; 
			
			if(R12_tmp % 2 == 1) tobit_r12[bit_index] = 1;
			else tobit_r12[bit_index] = 0; 
			
			if(R13_tmp % 2 == 1) tobit_r13[bit_index] = 1;
			else tobit_r13[bit_index] = 0; 
			
			if(R14_tmp % 2 == 1) tobit_r14[bit_index] = 1;
			else tobit_r14[bit_index] = 0; 
			
			if(R15_tmp % 2 == 1) tobit_r15[bit_index] = 1;
			else tobit_r15[bit_index] = 0; 
			
			R1_tmp  = R1_tmp >> 1;
            R2_tmp  = R2_tmp >> 1;
            R3_tmp  = R3_tmp >> 1;
            R4_tmp  = R4_tmp >> 1;
            R5_tmp  = R5_tmp >> 1;
            R6_tmp  = R6_tmp >> 1;
            R7_tmp  = R7_tmp >> 1;
            R8_tmp  = R8_tmp >> 1;
            R9_tmp  = R9_tmp >> 1;
            R10_tmp = R10_tmp >> 1;
            R11_tmp = R11_tmp >> 1;
            R12_tmp = R12_tmp >> 1;
            R13_tmp = R13_tmp >> 1;
            R14_tmp = R14_tmp >> 1;
            R15_tmp = R15_tmp >> 1;
        }
        
        for(int g=0;g < 64; g++){
            ROM0 << tobit_r1[63-g];
            ROM1 << tobit_r2[63-g];
            ROM2 << tobit_r4[63-g];
            ROM3 << tobit_r6[63-g];
            ROM4 << tobit_r8[63-g];
            ROM5 << tobit_r10[63-g];
            ROM6 << tobit_r12[63-g];
            ROM7 << tobit_r14[63-g];
        }
        for(int g=0;g < 64; g++){
           ROM1 << tobit_r3[63-g];
		   ROM2 << tobit_r5[63-g];
		   ROM3 << tobit_r7[63-g];
		   ROM4 << tobit_r9[63-g];
		   ROM5 << tobit_r11[63-g];
		   ROM6 << tobit_r13[63-g];
		   ROM7 << tobit_r15[63-g];
        }
		ROM0 << "\n";
        ROM1 << "\n";
        ROM2 << "\n";
        ROM3 << "\n";
        ROM4 << "\n";
        ROM5 << "\n";
        ROM6 << "\n";
        ROM7 << "\n";
    }
	ROM0.close();
	ROM1.close();
	ROM2.close();
	ROM3.close();
	ROM4.close();
	ROM5.close();
	ROM6.close();
	ROM7.close();
}
std::string SPMB::ZZtohex(ZZ zz_tmp){
    std::string string_tmp;
    std::vector<char> tmp_hex; 
    std::stringstream ss;
    
    long tmp;
    int length;
    length = 16;
    tmp_hex.resize(16);
    
    for(int i =0; i < length; i++){
        tmp  =  to_long(zz_tmp % 16);
        
        if(tmp == 0) tmp_hex[i] = '0'; 
        if(tmp == 1) tmp_hex[i] = '1'; 
        if(tmp == 2) tmp_hex[i] = '2'; 
        if(tmp == 3) tmp_hex[i] = '3'; 
        if(tmp == 4) tmp_hex[i] = '4'; 
        if(tmp == 5) tmp_hex[i] = '5'; 
        if(tmp == 6) tmp_hex[i] = '6'; 
        if(tmp == 7) tmp_hex[i] = '7'; 
        if(tmp == 8) tmp_hex[i] = '8'; 
        if(tmp == 9) tmp_hex[i] = '9'; 
        if(tmp == 10) tmp_hex[i] = 'a'; 
        if(tmp == 11) tmp_hex[i] = 'b'; 
        if(tmp == 12) tmp_hex[i] = 'c'; 
        if(tmp == 13) tmp_hex[i] = 'd'; 
        if(tmp == 14) tmp_hex[i] = 'e'; 
        if(tmp == 15) tmp_hex[i] = 'f'; 
        
        zz_tmp = zz_tmp >> 4;
    }
    for(int i = 1;i <= length; i++){
        ss << tmp_hex[length-i];
    }
    
    string_tmp = ss.str();
    
    return string_tmp;
}
std::string SPMB::ZZtohex_cp(ZZ zz_tmp){
    std::string string_tmp;
    std::vector<char> tmp_hex; 
    std::stringstream ss;
    
    long tmp;
    int length;
	double cp_w;
	cp_w = (double) CP_width;
	cp_w = cp_w / 4; // number digit in hex representment
    cp_w = ceil(cp_w);
	length = (long)cp_w;
    tmp_hex.resize(16);
    
    for(int i =0; i < length; i++){
        tmp  =  to_long(zz_tmp % 16);
        
        if(tmp == 0) tmp_hex[i] = '0'; 
        if(tmp == 1) tmp_hex[i] = '1'; 
        if(tmp == 2) tmp_hex[i] = '2'; 
        if(tmp == 3) tmp_hex[i] = '3'; 
        if(tmp == 4) tmp_hex[i] = '4'; 
        if(tmp == 5) tmp_hex[i] = '5'; 
        if(tmp == 6) tmp_hex[i] = '6'; 
        if(tmp == 7) tmp_hex[i] = '7'; 
        if(tmp == 8) tmp_hex[i] = '8'; 
        if(tmp == 9) tmp_hex[i] = '9'; 
        if(tmp == 10) tmp_hex[i] = 'a'; 
        if(tmp == 11) tmp_hex[i] = 'b'; 
        if(tmp == 12) tmp_hex[i] = 'c'; 
        if(tmp == 13) tmp_hex[i] = 'd'; 
        if(tmp == 14) tmp_hex[i] = 'e'; 
        if(tmp == 15) tmp_hex[i] = 'f'; 
        
        zz_tmp = zz_tmp >> 4;
    }
    for(int i = 1;i <= length; i++){
        ss << tmp_hex[length-i];
    }
    
    string_tmp = ss.str();
    
    return string_tmp;
}
std::string SPMB::BITtohex(std::vector<int> A_in,int length){
    std::string string_tmp;
    std::vector<char> tmp_hex; 
    std::stringstream ss;
    
    int tmp;
	int weight_tmp;
    int Number_digit;
    Number_digit = length / 4;
    tmp_hex.resize(Number_digit);
	
    for(int i =0; i < Number_digit; i++){        
		for(int j = 0; j < 4; j++){
			weight_tmp = A_in[4*i + j];
			weight_tmp = weight_tmp << j;
			tmp = tmp + weight_tmp;
		}
        if(tmp == 0) tmp_hex[i] = '0'; 
        if(tmp == 1) tmp_hex[i] = '1'; 
        if(tmp == 2) tmp_hex[i] = '2'; 
        if(tmp == 3) tmp_hex[i] = '3'; 
        if(tmp == 4) tmp_hex[i] = '4'; 
        if(tmp == 5) tmp_hex[i] = '5'; 
        if(tmp == 6) tmp_hex[i] = '6'; 
        if(tmp == 7) tmp_hex[i] = '7'; 
        if(tmp == 8) tmp_hex[i] = '8'; 
        if(tmp == 9) tmp_hex[i] = '9'; 
        if(tmp == 10) tmp_hex[i] = 'a'; 
        if(tmp == 11) tmp_hex[i] = 'b'; 
        if(tmp == 12) tmp_hex[i] = 'c'; 
        if(tmp == 13) tmp_hex[i] = 'd'; 
        if(tmp == 14) tmp_hex[i] = 'e'; 
        if(tmp == 15) tmp_hex[i] = 'f'; 
        
		tmp = 0;
		weight_tmp = 0;
    }
    for(int i = (Number_digit-1);i >= 0 ; i--){
        ss << tmp_hex[i];
    }
    
    string_tmp = ss.str();
    
    return string_tmp;
}











