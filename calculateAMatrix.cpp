/*
 * this code is the member function of c++ class CalculateAMtrix, which 
 * gets elemts and nodes array to calculate the "A matrix"
 *
 */

#include "calculateAMatrix.hpp"
#include <cstddef>
#include <cmath>

void CalculateAMatrix::getAllArray(int ** array1ptr, double ** array2ptr, const size_t rowDim1, const size_t colDim1, const size_t rowDim2, const size_t colDim2)
{

  // make a object copy of Element array
  if (this->lengthRowElement == 0 && this->lengthColElement == 0 && this->elementArrayPtr == NULL)
        {
	  // delete[] this->arr;
          // this->len_row = 0;
	  // this->len_col = 0;
	  
          // allocate the row pointers array, and allocate column array for every
	  // pointer elments of row pointer array.
	  // If you forget this step, you will get "Segmentation fault (core dumped)"
	  // error
             this->elementArrayPtr = new int * [rowDim1];
             for(size_t idx_row = 0; idx_row < rowDim1; ++idx_row)
	     {
	         elementArrayPtr[idx_row] = new int[colDim1];
	     }

	  
        }
  // copy array element into the Object-local Element array 
  for (size_t idx_row = 0; idx_row < rowDim1; ++idx_row)
	  {
	    for (size_t idx_col = 0; idx_col < colDim1; ++idx_col)
                this->elementArrayPtr[idx_row][idx_col] = array1ptr[idx_row][idx_col];
	  }  

  this->lengthRowElement = rowDim1;
  this->lengthColElement = colDim1;

  
  // make a object copy of Node array
  if (this->lengthRowNode == 0 && this->lengthColNode == 0 && this->nodeArrayPtr == NULL)
        {
	  // delete[] this->arr;
          // this->len_row = 0;
	  // this->len_col = 0;
	  
          // allocate the row pointers array, and allocate column array for every
	  // pointer elments of row pointer array.
	  // If you forget this step, you will get "Segmentation fault (core dumped)"
	  // error
             this->nodeArrayPtr = new double * [rowDim2];
             for(size_t idx_row = 0; idx_row < rowDim2; ++idx_row)
	     {
	         nodeArrayPtr[idx_row] = new double[colDim2];
	     }

	  
        }

  // copy array element into the Object-local Node array 
  for (size_t idx_row = 0; idx_row < rowDim2; ++idx_row)
	  {
	    for (size_t idx_col = 0; idx_col < colDim2; ++idx_col)
                this->nodeArrayPtr[idx_row][idx_col] = array2ptr[idx_row][idx_col];
	  }  

  this->lengthRowNode = rowDim2;
  this->lengthColNode = colDim2;

 
}

// void Mat::setMat(double * src, const size_t rowDim, const size_t colDim)
// void Mat::setMat(double * src, const size_t rowDim, const size_t colDim)
/*void Mat::CheckArray(double ** src, const size_t rowDim, const size_t colDim)   
{
  //  std::cout << *src << std::endl;
  std::cout << src[0] << std::endl;
  std::cout << src[1] << std::endl;
  //
  std::cout << " first column pointer : " << src[0] << std::endl;
  std::cout << " second column pointer : " << src[1] << std::endl;

  
  //  using the column array pinter (src[rowidx]) to fetch the column elemets 
   
  for(size_t rowidx = 0; rowidx < rowDim; ++rowidx)
    {
      for(size_t colidx = 0; colidx < colDim; ++colidx)
	std::cout << src[rowidx][colidx] << "  ";
      std::cout << std::endl;
    }

}*/

void CalculateAMatrix::showOriginalArr(int ** array1ptr, double ** array2ptr, const size_t rowDim1,const size_t colDim1, const size_t rowDim2, const size_t colDim2)
{

  std::cout << " first pointer to pointer of orignal Element array: " << array1ptr << std::endl;
  std::cout << " second pointer to pointer of original Node array: " << array2ptr << std::endl;

  std::cout << " Show the first array passed in " << std::endl;
  for(size_t rowidx = 0; rowidx < rowDim1; ++rowidx)
    {
      for(size_t colidx = 0; colidx < colDim1; ++colidx)
	std::cout << array1ptr[rowidx][colidx] << "  ";
      std::cout << std::endl;
    }

  std::cout << " Show the second array passed in " << std::endl;
  for(size_t rowidx = 0; rowidx < rowDim2; ++rowidx)
    {
      for(size_t colidx = 0; colidx < colDim2; ++colidx)
	std::cout << array2ptr[rowidx][colidx] << "  ";
      std::cout << std::endl;
    }

}

void CalculateAMatrix::showObjectLocalArray()   
{
  std::cout << " pointer to pointer of Object-local Element array: " << elementArrayPtr << std::endl; 
  for (size_t idx_row = 0; idx_row < lengthRowElement; ++idx_row)
	  {
	    for (size_t idx_col = 0; idx_col < lengthColElement; ++idx_col)
	      std::cout <<  this->elementArrayPtr[idx_row][idx_col] << " ";
	    std::cout << std::endl;
	  }  

  std::cout << " pointer to pointer of Object-local Node array: " << nodeArrayPtr << std::endl; 
  for (size_t idx_row = 0; idx_row < lengthRowNode; ++idx_row)
	  {
	    for (size_t idx_col = 0; idx_col < lengthColNode; ++idx_col)
	      std::cout <<  this->nodeArrayPtr[idx_row][idx_col] << " ";
	    std::cout << std::endl;
	  }  
  
}

// memeber function AMatrixCalculatorRun() is  the main functionality of this class
// i.e., calculate the A Matrix with object-local arraies
void CalculateAMatrix::AMatrixCalculatorRun()
{
  std::cout << " calculator starts " << std::endl;
  
  // dynamically allocate the array of A-Matrix, using pointer array for store the Row
  AMatrixArrayPtr = new double * [lengthRowNode];

  std::cout << " ptr array allocate successfully " << std::endl;
  
  for(size_t idx_row = 0; idx_row < lengthRowNode; ++idx_row)
    {
      AMatrixArrayPtr[idx_row] = new double[lengthRowNode];
    }

  std::cout << " A-Mtrix columns allocate successfully " << std::endl;
  
  // put 0.0 into the allocated A-Matrix, innitilization
  for(size_t idx_row = 0; idx_row < lengthRowNode; ++idx_row)
    {
      std::cout << " assign loop 1-level, idx_row is  " << idx_row << std::endl;
      for(size_t idx_col = 0; idx_col < lengthRowNode; ++idx_col)
	{
          std::cout << " assign loop 2-level, idx_col is  " << idx_col << std::endl;
	  AMatrixArrayPtr[idx_row][idx_col] = 0.0;
	}  
    }

  std::cout << " dynamical allocated successfull " << std::endl;

  // calculate the Matrix-elements and save the results into the initilized A-Matrix array
  size_t idx_element = 0;
  while( idx_element < lengthRowElement )
    {

     std::cout << " loop start: idx_element is " << idx_element << std::endl;
      
     // data prerparations  
     int nodeIndex0Elem = elementArrayPtr[idx_element][0];
     int nodeIndex1Elem = elementArrayPtr[idx_element][1];
     int nodeIndex2Elem = elementArrayPtr[idx_element][2];
     int nodeIndex[] = {nodeIndex0Elem, nodeIndex1Elem, nodeIndex2Elem};

     std::cout << " node indcies with idx_element " << idx_element << std::endl;
     std::cout << " nodeIndex in element " << idx_element << " are " << nodeIndex[0] << " " << nodeIndex[1] << " " << nodeIndex[2] << std::endl;

     double x0 = nodeArrayPtr[nodeIndex0Elem][0];
     double y0 = nodeArrayPtr[nodeIndex0Elem][1];

     std::cout << " node indcies 0th coordinates with idx_element " << idx_element << std::endl;
     
     double x1 = nodeArrayPtr[nodeIndex1Elem][0];
     double y1 = nodeArrayPtr[nodeIndex1Elem][1];

     std::cout << " node indcies 1th coordinates with idx_element " << idx_element << std::endl;

     double x2 = nodeArrayPtr[nodeIndex2Elem][0];
     double y2 = nodeArrayPtr[nodeIndex2Elem][1];

     std::cout << " node indcies and coordinates with idx_element " << idx_element << std::endl;

     // xi, eta, omega
     double xi[] = {x1-x2, x2-x0, x0-x1};
     double eta[] = {y1-y2, y2-y0, y0-y1};
     double omega[] = {x1*y2-x2*y1, x2*y0-x0*y2, x0*y1-x1*y0};
     
     double Di = 0; // initilization Di
       for(size_t i = 0; i < 3; ++i) { Di += omega[i]; }
     double Si = (std::abs (Di))/2;

     std::cout << " primary data been prepared with " << idx_element << std::endl;

     // Wall effect
     // double gx = gamma1 + 2*gamma2;
     // double gy = gamma1;
        double gx = 0.1;
	double gy = 0.05;

	/*
           all data have been prepared, next is build 
           the element matrix ax2ij, ay2ij, a2ij
         */

     //	define a 3 by 3 array for ax2ij
	double ax2ij[3][3] = {
	  {0.0, 0.0, 0.0},
	  {0.0, 0.0, 0.0},
	  {0.0, 0.0, 0.0}
	};

     //	define a 3 by 3 array for ay2ij
	double ay2ij[3][3] = {
	  {0.0, 0.0, 0.0},
	  {0.0, 0.0, 0.0},
	  {0.0, 0.0, 0.0}
	};
	
     //	define a 3 by 3 array for a2ij
	double a2ij[3][3] = {
	  {0.0, 0.0, 0.0},
	  {0.0, 0.0, 0.0},
	  {0.0, 0.0, 0.0}
	};

     // calculate ax2ij and ay2ij
	for(size_t o = 0; o < 3; ++o)
	  {
	    ax2ij[o][0] = (1.0/(2.0*(Di*Di)))*eta[o]*eta[0]*Si*gx;
	    ax2ij[o][1] = (1.0/(2.0*(Di*Di)))*eta[o]*eta[1]*Si*gx;
	    ax2ij[o][2] = (1.0/(2.0*(Di*Di)))*eta[o]*eta[2]*Si*gx;

	    ay2ij[o][0] = (1.0/(2.0*(Di*Di)))*xi[o]*xi[0]*Si*gy;
	    ay2ij[o][1] = (1.0/(2.0*(Di*Di)))*xi[o]*xi[1]*Si*gy;
	    ay2ij[o][2] = (1.0/(2.0*(Di*Di)))*xi[o]*xi[2]*Si*gy;
	  }
     //	a2ij = ax2ij + ay2ij
	for(size_t i = 0; i < 3; ++i)
	  {
	    for(size_t j = 0; j < 3; ++j)
	      {
		a2ij[i][j] = ax2ij[i][j] + ay2ij[i][j];
	      }
         
	  }

	std::cout << " a2ij matrix now is " << std::endl;
        for(size_t i = 0; i < 3; ++i)
	  {
	    for(size_t j = 0; j < 3; ++j)
	      std::cout << a2ij[i][j] << " ";
	    std::cout << std::endl;
	  }

    // traspose the a2ij
	double a2ji[3][3] = {
	  {0.0, 0.0, 0.0},
	  {0.0, 0.0, 0.0},
	  {0.0, 0.0, 0.0}
	};

       for(size_t i = 0; i < 3; ++i)
	 {
	   for(size_t j = 0; j < 3; ++j)
	     {
	       a2ji[j][i] =a2ij[i][j]; // matrix transpose 
	     }
	 }

       std::cout << " a2ji matrix now is " << std::endl;
        for(size_t i = 0; i < 3; ++i)
	  {
	    for(size_t j = 0; j < 3; ++j)
	      std::cout << a2ji[i][j] << " ";
	    std::cout << std::endl;
	  }

    // symmetrize the element array
	double A2ij[3][3] = {
	  {0.0, 0.0, 0.0},
	  {0.0, 0.0, 0.0},
	  {0.0, 0.0, 0.0}
	};

        for(size_t i = 0; i < 3; ++i)
	  {
	    for(size_t j = 0; j < 3; ++j)
	      {
	        A2ij[i][j] = (1.0/2.0)*(a2ij[i][j] + a2ji[i][j]);
		std::cout << A2ij[i][j] << " and " << (1/2)*(a2ij[i][j]+a2ji[i][j]) << std::endl;
	      }	
	    
	  }

	std::cout << " A2ij matrix now is " << std::endl;
        for(size_t i = 0; i < 3; ++i)
	  {
	    for(size_t j = 0; j < 3; ++j)
	      std::cout << A2ij[i][j] << " ";
	    std::cout << std::endl;
	  }
        
	
	std::cout << " idx_element is " << idx_element << " till now " << std::endl; 

	
    // add the element reslults onto the A-Matrix array
	for(size_t row = 0; row < 3; ++row)
	  {

           std::cout << " row is " << row << std::endl;
	    
	   AMatrixArrayPtr[nodeIndex[row]][nodeIndex[0]] = A2ij[row][0] + AMatrixArrayPtr[nodeIndex[row]][nodeIndex[0]];
	   AMatrixArrayPtr[nodeIndex[row]][nodeIndex[1]] = A2ij[row][1] + AMatrixArrayPtr[nodeIndex[row]][nodeIndex[1]];
	   AMatrixArrayPtr[nodeIndex[row]][nodeIndex[2]] = A2ij[row][2] + AMatrixArrayPtr[nodeIndex[row]][nodeIndex[2]];
	    
	  }

	std::cout << "idx_element is " << idx_element << " after A-Matrix assignment "  << std::endl;
     
        ++idx_element;
    }
}

void  CalculateAMatrix::showObjectLocalAMatrix()
{
  std::cout << " pointer to pointer of Object-local A-Matrix array: " << AMatrixArrayPtr << std::endl; 
  for (size_t idx_row = 0; idx_row < lengthRowNode; ++idx_row)
	  {
	    for (size_t idx_col = 0; idx_col < lengthRowNode; ++idx_col)
	      std::cout <<  this->AMatrixArrayPtr[idx_row][idx_col] << " ";
	    std::cout << std::endl;
	  }  
}

// Here, you can define the return function as "const double **", but
// you can difine it as double ** const, i.e., you can return const pointer,
// but you can not require pointer to point const object
// double ** const Mat::returnPointerOfMatrix() const
// {
//  return this->arr;
// }

double ** const CalculateAMatrix::returnAMatrixPtr() const
{
  return this->AMatrixArrayPtr;
}


