/*
 * 
 * CalculateA2lijInOneElement.h
 *
 *
 *
 */

#ifndef CLASSCALCULATEA2LIJINONEELEMENT_HPP
#define CLASSCALCULATEA2LIJINONEELEMENT_HPP


#include "ClassCalculateA2lijInOneElement.hpp"
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>

class CalculateA2lij {
     
      private:
         double *arr;
         double A2lij[9];
         size_t length;    
        
      public:
          // double xi;
          // double xi1;
          
      // reveive the patition data as constant array    
      void CalculateA2lijInOneElement(const double elements[26784], const double nodes[9206], 
                                      double gamma1, double gamma2, double l);
      
      // return the A2lij[9] by pointer arr
      const double * const returnResult(size_t length) const;
};      

#endif


