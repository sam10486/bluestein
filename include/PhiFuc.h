#ifndef _PHIFUC_H_
#define _PHIFUC_H_

#include <vector>
#include <iostream>
#include <NTL/ZZ.h>

using namespace NTL;

/*
calculate Euler’s Totient Function
input m,then output phi(m);
*/
class PhiFuc{
public:
  long m;
  long phim;
  std::vector<long> index; // record index i which gcd(i,m) = 1
  
  void calculate_phim(long m_th); 
  void calculate_index();
};

#endif
