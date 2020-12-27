/*
BSD 2-Clause License

Copyright (c) 2020, Andrey Kudryavtsev (andrewkoudr@hotmail.com)
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#pragma once

#include <iostream>
#include <assert.h>
#include "Types.h"
#include "Allocator.h"


                              // printout progress
extern bool progressprint;

/**
  Banded matrix (maybe non-symmetric). Essentially this is mapping from [i,j],
  i - row, counted from top to bottom, j - column.
  Matrix is stored by row parts within the band. Band itself is rounded to
  4 float values or 2 double values to facilitate SIMD operations. It also
  stores both asymmetric and symmetric parts. That is,
  XXXXDXXX is one-row representaton for band width 1 (pure diagonal) 2,3 or 4
  symmetric band widths.
  Maximum order of matrix is defined by shared memory available, e.g. this
  is 2GB for a 2GB computer with 1GB memory free. It means a banded matrix 
  50000x5000 of doubles or 70000x7000 of floats is possible.
*/

template <class T> class BandedMatrix 
{
public:
                              // memory to hold big matrix
  Allocator *mem = nullptr;
                              // memory allocated?
  bool OK;
                              // matrix order
  size_t height;
                              // original half band width, including diagonal
  size_t bandwidth;
                              // half band width, rounded up to 4 for
                              // floats and to 2 for doubles (including
                              // diagonal)
  size_t ebandwidth;
                              // full width (ebandwidth * 2)
  size_t width;
                              // constructors
  BandedMatrix() : OK(false),bandwidth(0),height(0),
    ebandwidth(0),width(0)
  {
  }
                              // n - matrix order, bandwidth -
                              // original bandwidth (including diagonal)
  BandedMatrix(Allocator* pmem, size_t n, size_t pbandwidth);
                              // destructor
  ~BandedMatrix()
  {
  }
                              // SET/GET
															// get element address; i,j are positions in GLOBAL
                              // square matrix, each being in the range 
                              // 0..(height - 1)
	T *getElement(size_t i, size_t j);
                              // set element (no checks on i,j, be careful)
	void setElement(size_t i, size_t j, T value);
                              // add to element value (no checks on i,j, be careful)
	void addToElement(size_t i, size_t j, T value);
                              // k = i * N + j (i,j are for square matrix)
  T *operator[](size_t k);
                              // change equation to all zeroes and 1.0 at diagonal
  void degenerateEquation(size_t index);
};

template <class T> class FEMVirtMatrix : public BandedMatrix<T> 
{
public:
                              // CONSTRUCTION/DESTRUCTION
															// constructor; n is number of unknowns,
                              // bandwidth is half band width plus 1
	FEMVirtMatrix() : BandedMatrix()
  {
  }
	FEMVirtMatrix(Allocator* pmem, size_t n, size_t bandwidth) : BandedMatrix<T>(pmem,n,bandwidth)
  {
  }
															// destructor
	~FEMVirtMatrix()
  {
  }
                              // SOLVER
															// solve system by Gauss elimination;
															// system matrix IS SPOILT; use StoreMatrix()/RestoreMatrix()
															// to avoid the problem; returns false if close to zero
                              // pivot is encountered
	bool solveSystem(T *B, T zero, int numthreads);
                              // simple slow straightforward solution; use for
                              // testing only
  bool solveSystemSimple(T *B, T zero);

                              // MATRIX MULTIPLIER
                              // use it to e.g. check residuals after solution
															// multiply matrix by vector V with the result in R;
															// vectors V and R MUST BE 4 FLOATS (2 DOUBLES) 
                              // LONGER than the order of the matrix, with the 
                              // rest filled by zeroes
	bool multiply(T *V, T *R, int numthreads);

                              // AUX : STORE/RESTORE MATRIX
															// store matrix
	bool storeMatrix();
															// restore matrix
	bool restoreMatrix(bool freestored);
};


template <class T> bool FEMVirtMatrix<T>::storeMatrix()
{
  return BandedMatrix<T>::mem->storeCopy("FEMVirtMatrix.bin");
}

template <class T> bool FEMVirtMatrix<T>::restoreMatrix(bool freestored)
{
  return BandedMatrix<T>::mem->restoreCopy(freestored);
}

template <class T> 
bool FEMVirtMatrix<T>::solveSystemSimple(T *B, T zero)
{
  LINT K,K1,J,I,Ie,Je;
  T AKK;

  LINT N = static_cast<LINT>(BandedMatrix<T>::height);
  LINT MB = static_cast<LINT>(BandedMatrix<T>::bandwidth) - 1;

  for (K = 0; K < N; K++)
  {
    if (progressprint && K % 100 == 0)
    {
      std::cout << K << " of " << N << "\r";
    }

    K1 = K + 1;
    AKK = *BandedMatrix<T>::getElement(K,K);
    if (fabs(AKK) < zero)
		{	
			return false;
    }

    B[K] /= AKK;
    if (K == (N-1)) 
      break;

		Je = K1 + MB; if (Je > (N - 1)) Je = N - 1;
    for (J = K1; J <= Je; J++)
    {
      *BandedMatrix<T>::getElement(K,J) /= AKK; // *this[K,J] /= AKK

			Ie = K1 + MB; if (Ie > (N - 1)) Ie = N - 1;
      for (I = K1; I <= Ie; I++) 
      {
        *BandedMatrix<T>::getElement(I,J) -= (*BandedMatrix<T>::getElement(I,K) * (*BandedMatrix<T>::getElement(K,J))); // *this[I,J] -= *this[I,K] * *this[K,J]
      };
      B[J] -= (*BandedMatrix<T>::getElement(J,K) * B[K]);
    }
  }

  if (progressprint)
    std::cout << std::endl;

  for (;;)
  {
    K1 = K;
    K -= 1;
    if (K < 0) break;
    LINT K2 = K1 + MB;
    if (K2 > (N - 1)) K2 = N - 1;
    for (J = K1; J <= K2; J++) 
    {
      B[K] -= (*BandedMatrix<T>::getElement(K,J) * B[J]);
    }
  }

  return true;
}


/** Fast SIMD subtraction of floats arbitrarily aligned,
  v0 = v0 - v1, size must be % 4 = 0 for floats
  and % 2 = 0 for doubles */
template <class T> void vectorsSubtract(T* v0, T* v1, size_t length);
