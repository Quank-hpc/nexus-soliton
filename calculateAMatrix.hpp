/*
 * this code is the header file of c++ class CalculateAMtrix, which 
 * gets elemts and nodes array to calculate the "A matrix"
 *
 * the inputs are the pointer arraies of element array (2D) and node array (2D)
 */


#ifndef calculateAMatrix_hpp
#define calculateAMatrix_hpp

#include <cstddef>
#include <iostream>

class CalculateAMatrix {
    int ** elementArrayPtr;
    double ** nodeArrayPtr;
    double ** AMatrixArrayPtr;

      
public:

  size_t lengthRowElement;
  size_t lengthColElement;

  size_t lengthRowNode;
  size_t lengthColNode;

  CalculateAMatrix(): lengthRowElement(0), lengthColElement(0), elementArrayPtr(NULL), lengthRowNode(0), lengthColNode(0), nodeArrayPtr(NULL) { }; // default constructor to initilaze data member

    // pass in the elements array and node array
  void getAllArray(int ** array1ptr, double ** array2ptr, const size_t rowDim1,const size_t colDim1, const size_t rowDim2, const size_t colDim2);
  
    // Show what array have been passed in this "A-Matrix" calculator
  void showOriginalArr(int ** array1ptr, double ** array2ptr, const size_t rowDim1,const size_t colDim1, const size_t rowDim2, const size_t colDim2);
    // void setMat(double * src, const size_t rowDim, const size_t colDim);

  void showObjectLocalArray();

  void AMatrixCalculatorRun();

  void showObjectLocalAMatrix();

  double ** const returnAMatrixPtr() const;
  
  
    // make an object copy of PrimaryArray by allocate new memmory
  // void calculateA();
  

  // double ** const returnPointerOfTheAMatrix() const;
  
    // // get arr pointer
    // const double * const getMat(size_t len) const;

    //  // get arr length
    // size_t getLength() const;

    // // copy arr to dest
    // void copyMat(double * dest, size_t len);
};


#endif

