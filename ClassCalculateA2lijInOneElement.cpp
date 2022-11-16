/*
 *
 *  CalculateA2lijInOneElement.cpp
 *
 */


#include "ClassCalculateA2lijInOneElement.hpp"
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>


void CalculateA2lij::CalculateA2lijInOneElement(const double elements[26784], const double
  nodes[9206], double gamma1, double gamma2, double l)
{
  double a2xlij[9];
  double a2ylij[9];
  double Di;
  double Si;
  
  double d_xi_idx_0_tmp;
  double eta_idx_0;
  double eta_idx_0_tmp;
  double eta_idx_1;
  double eta_idx_1_tmp;
  double eta_idx_2;
  
  double gx; //gx
  
  double xi_idx_0;
  double xi_idx_1;
  double xi_idx_2;
  
  int xi_idx_0_tmp;
  int xi_idx_1_tmp;
    
  int b_xi_idx_0_tmp;
  int c_xi_idx_0_tmp;
  
  
  
  xi_idx_0_tmp = 3 * (static_cast<int>(l) - 1); //typecast double l to int l
  
  b_xi_idx_0_tmp = (static_cast<int>(elements[xi_idx_0_tmp + 2]) - 1) << 1;
  Di = nodes[b_xi_idx_0_tmp];
  
  c_xi_idx_0_tmp = (static_cast<int>(elements[xi_idx_0_tmp + 1]) - 1) << 1;
  d_xi_idx_0_tmp = nodes[c_xi_idx_0_tmp];
  xi_idx_0 = d_xi_idx_0_tmp - Di;
  xi_idx_1_tmp = (static_cast<int>(elements[xi_idx_0_tmp]) - 1) << 1;
  
  Si = nodes[xi_idx_1_tmp];
  
  xi_idx_1 = Di - Si;
  xi_idx_2 = Si - d_xi_idx_0_tmp;
  gx = nodes[b_xi_idx_0_tmp + 1];
  
  eta_idx_0_tmp = nodes[c_xi_idx_0_tmp + 1];
  eta_idx_0 = eta_idx_0_tmp - gx;
  eta_idx_1_tmp = nodes[xi_idx_1_tmp + 1];
  eta_idx_1 = gx - eta_idx_1_tmp;
  eta_idx_2 = eta_idx_1_tmp - eta_idx_0_tmp;
  Di = ((d_xi_idx_0_tmp * gx - Di * eta_idx_0_tmp) + (Di * eta_idx_1_tmp - Si *
         gx)) + (Si * eta_idx_0_tmp - d_xi_idx_0_tmp * eta_idx_1_tmp);
  
  Si = std::abs(Di) / 2.0;

 
  // Wall of Gamma /////////////
  gx = gamma1 + 2.0 * gamma2;

  // e2 parallel with H0
  // /////////////////////////////
  // //////////////////////////////
  // 
  d_xi_idx_0_tmp = 1.0 / (2.0 * (Di * Di));
  Di = d_xi_idx_0_tmp * eta_idx_0;
  a2xlij[0] = Di * eta_idx_0 * Si * gx;
  a2xlij[3] = Di * eta_idx_1 * Si * gx;
  a2xlij[6] = Di * eta_idx_2 * Si * gx;
  Di = d_xi_idx_0_tmp * xi_idx_0;
  
  a2ylij[0] = Di * xi_idx_0 * Si * gamma1;
  a2ylij[3] = Di * xi_idx_1 * Si * gamma1;
  a2ylij[6] = Di * xi_idx_2 * Si * gamma1;
  Di = d_xi_idx_0_tmp * eta_idx_1;
  
  a2xlij[1] = Di * eta_idx_0 * Si * gx;
  a2xlij[4] = Di * eta_idx_1 * Si * gx;
  a2xlij[7] = Di * eta_idx_2 * Si * gx;
  Di = d_xi_idx_0_tmp * xi_idx_1;
  
  a2ylij[1] = Di * xi_idx_0 * Si * gamma1;
  a2ylij[4] = Di * xi_idx_1 * Si * gamma1;
  a2ylij[7] = Di * xi_idx_2 * Si * gamma1;
  Di = d_xi_idx_0_tmp * eta_idx_2;
  
  a2xlij[2] = Di * eta_idx_0 * Si * gx;
  a2xlij[5] = Di * eta_idx_1 * Si * gx;
  a2xlij[8] = Di * eta_idx_2 * Si * gx;
  Di = d_xi_idx_0_tmp * xi_idx_2;
  
  a2ylij[2] = Di * xi_idx_0 * Si * gamma1;
  a2ylij[5] = Di * xi_idx_1 * Si * gamma1;
  a2ylij[8] = Di * xi_idx_2 * Si * gamma1;

  //  matrix elements symmetrized 
  for (xi_idx_0_tmp = 0; xi_idx_0_tmp < 9; xi_idx_0_tmp++) {
      
    a2xlij[xi_idx_0_tmp] += a2ylij[xi_idx_0_tmp];
    
  }

  for (xi_idx_0_tmp = 0; xi_idx_0_tmp < 3; xi_idx_0_tmp++) {
    A2lij[3 * xi_idx_0_tmp] = 0.5 * (a2xlij[3 * xi_idx_0_tmp] +
      a2xlij[xi_idx_0_tmp]);
    
    xi_idx_1_tmp = 3 * xi_idx_0_tmp + 1;
    A2lij[xi_idx_1_tmp] = 0.5 * (a2xlij[xi_idx_1_tmp] + a2xlij[xi_idx_0_tmp + 3]);
    
    xi_idx_1_tmp = 3 * xi_idx_0_tmp + 2;
    A2lij[xi_idx_1_tmp] = 0.5 * (a2xlij[xi_idx_1_tmp] + a2xlij[xi_idx_0_tmp + 6]);
  }
  
  arr=A2lij; //transfer the pointer of array A2lij to object data member pointer arr
}

const double * const CalculateA2lij::returnResult(size_t length) const
{
   
    return this->arr; // this pointer of function object return address of result
}


