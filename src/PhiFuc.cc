#include <vector>
#include <iostream>
#include <NTL/ZZ.h>

using namespace NTL;

#include "PhiFuc.h"

void PhiFuc::calculate_phim(long m_th){
     m = m_th;
     long number = 0;
     long tmp;
     for(long i=0;i < m_th; i++){
         tmp = i + 1;
         tmp = GCD(tmp,m);
         if(tmp == 1) number ++;
     }
     phim = number;
}

void PhiFuc::calculate_index(){
     index.resize(phim);
     
     long tmp = 0;
     long number = 0;
     for(long i=0;i<m;i++){
         tmp = i+1;
         tmp = GCD(tmp,m);
         if(tmp == 1){
             index[number] = (i+1);
             number ++;
         }
     }
    
}
