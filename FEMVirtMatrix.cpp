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

                              // headers
#include "Types.h"
#include "FEMVirtMatrix.h"
#include "VirtMatrixThreads.h"
                              // intrinsics
#include <xmmintrin.h>
#include <emmintrin.h>
#include <assert.h>
#include <vector>
#include <thread>


bool progressprint = false;

// Generates ranges array of numranges + 1 elements
template <class T>
void getRanges(T max, int numranges, std::vector<T>& ranges)
{
  T step = max / numranges;
  LIMIT_MIN(step, 1);

  ranges.resize(numranges + 1);
  T value = T(0);

  for (size_t i = 0; i <= numranges; i++)
  {
    if (value >= max)
      value = max;

    ranges[i] = value;

    value += step;
  }

  ranges[numranges] = max;
}

BandedMatrix<float>::BandedMatrix(Allocator* pmem, size_t n, size_t pbandwidth) 
{
  mem = pmem;
                              // compute real matrix band width
  height = n;
  bandwidth = pbandwidth;

  ebandwidth = bandwidth / 4;
  if (bandwidth % 4) 
    ebandwidth++;
  ebandwidth *= 4;
  width = ebandwidth * 2;
                              // allocate floats
  size_t size = width * height * sizeof(float);
                              // get matrix size in bytes
                              //!!! this will create anonymous shared memory
  OK = mem->init("",size,true);

  assert(OK);
}

BandedMatrix<double>::BandedMatrix(Allocator* pmem, size_t n, size_t pbandwidth) 
{
  mem = pmem;
                              // compute real matrix band width
  height = n;
  bandwidth = pbandwidth;

  ebandwidth = bandwidth / 2;
  if (bandwidth % 2) 
    ebandwidth++;
  ebandwidth *= 2;
  width = ebandwidth * 2;
                              // allocate doubles
  size_t size = width * height * sizeof(double);
                              // get matrix size in bytes
                              //!!! this will create anonymous shared memory
  OK = mem->init("",size,true);

  assert(OK);
}

float *BandedMatrix<float>::getElement(size_t i, size_t j)
{
  assert(OK);
  assert(i >= 0 && i < height);
  assert(j >= 0 && j < height);

	return reinterpret_cast<float *>(mem->buffer) + i * width + j + ebandwidth - i;
}

double *BandedMatrix<double>::getElement(size_t i, size_t j)
{
  assert(OK);
  assert(i >= 0 && i < height);
  assert(j >= 0 && j < height);

	return reinterpret_cast<double *>(mem->buffer) + i * width + j + ebandwidth - i;
}

void BandedMatrix<float>::setElement(size_t i, size_t j, float value)
{
  *getElement(i,j) = value;
}

void BandedMatrix<double>::setElement(size_t i, size_t j, double value)
{
  *getElement(i,j) = value;
}

void BandedMatrix<float>::degenerateEquation(size_t index)
{
  unsigned char *e = mem->buffer + width * index * sizeof(float);
  memset(e,0,width * sizeof(float));
  setElement(index,index,float(1.0));
}

void BandedMatrix<double>::degenerateEquation(size_t index)
{
  unsigned char *e = mem->buffer + width * index * sizeof(double);
  memset(e,0,width * sizeof(double));
  setElement(index,index,double(1.0));
}

void BandedMatrix<float>::addToElement(size_t i, size_t j, float value)
{
  *getElement(i,j) += value;
}

void BandedMatrix<double>::addToElement(size_t i, size_t j, double value)
{
  *getElement(i,j) += value;
}

float *BandedMatrix<float>::operator[](size_t k)
{
  size_t i = k / height;
  size_t j = k % height;
  return getElement(i,j);
}

double *BandedMatrix<double>::operator[](size_t k)
{
  size_t i = k / height;
  size_t j = k % height;
  return getElement(i,j);
}

bool FEMVirtMatrix<float>::solveSystem(float *B, float zero, int numthreads)
{
  if (!OK)
    return false;

  assert(numthreads > 0);
                              // data to pass to threads
  std::vector<ThreadData> data(numthreads);
  data[0].init(numthreads);
                              // create threads
  if (numthreads > 1)
  {
                              // start threads
    for (int i = 0; i < numthreads; i++)
    {
                              // arguments to threads
      data[i].init(numthreads);
    }
  }
															// get first RHS vector
	float *b = B;
															// pointer to diagonal element of leading row
	float *upperrow = getElement(0,0);
															// go down
	size_t height1 = height - 1;
                              // offset to the next row in floats
	size_t shiftdown = width - 1;

	for (size_t k = 0; k < height; k++)
	{
		if (std::abs(*upperrow) < zero)
		{
      if (numthreads > 1)
      {
        for (int i = 0; i < numthreads; i++)
        {
          data[i].upperrow = nullptr;
        }
      }

			return false;
		}
															// make triangle
		size_t stepsright = ebandwidth >> 2;
		size_t numstepsdown = bandwidth - 1;
		if ((k + numstepsdown) > height1)
		{
			numstepsdown = height1 - k;
		}
                              // divide leading row by pivot
                              // inverted diagonal element for leading row
	  float diag = float(1.0) / *upperrow;
                                // divide leading row by its pivot
    register __m128 xdiag = _mm_set1_ps(diag);
    float *x = upperrow;
    for (size_t j = 0; j < stepsright; j++)
    {
      _mm_storeu_ps(x,_mm_mul_ps(_mm_loadu_ps(x),xdiag));
      x += 4; 
    }
                              // divide RHS by pivot
    *b *= diag;
                              // here we have :
                              // upperrow     - pointer to lead row
                              // b            - pointer to RHS for leading row
                              // stepsright   - active row length in XMM registers
                              // shiftdown    - offset to next row
                              // numstepsdown - number of rows below lead row
    if (numthreads <= 1)
    {
      data[0].upperrow = upperrow;
      data[0].b = b;
      data[0].stepsright = stepsright;
      data[0].shiftdown = shiftdown;
      data[0].numstepsdown = numstepsdown;
      data[0].rowoffset = 1;

      ThreadProc4(&data[0]);
    } else
    {
                              // divide numstepsdown rows between threads
      std::vector<int> ranges;
      getRanges<int>(static_cast<int>(numstepsdown),numthreads,ranges);
                              // prepare data for each thread and signal to start
      for (int i = 0; i < numthreads; i++)
      {
        data[i].upperrow = upperrow;
        data[i].b = b;
        data[i].stepsright = stepsright;
        data[i].shiftdown = shiftdown;
        data[i].numstepsdown = ranges[i + 1] - ranges[i];
        data[i].rowoffset = ranges[i] + 1;
      }

      // threads
      std::vector<std::thread> threads(numthreads);

                              // reset all events
      for (int i = 0; i < numthreads; i++)
      {
        // create threads
        threads[i] = std::thread(&ThreadProcFloat, &data[i]);
      }

      for (auto& t : threads)
      {
        t.join();
      }
    }
															// get next diagonal element
		upperrow += width;
															// get next RHS vector
		b++;

    if (progressprint && k % 100 == 0)
      printf("%d\r",static_cast<int>(k));
	}

  if (progressprint)
    printf("\n");
                              // stop all threads
  if (numthreads > 0)
  {
                              // set data to stop all threads
    for (int i = 0; i < numthreads; i++)
    {
      data[i].upperrow = nullptr;
    }
  }
															// back substitution, not yet parallel
	b--;
	float *bj;
	float *aj;
	for (LINT k = (height - 2); k >= 0; k--)
	{
		b--;
		aj = getElement(k,k + 1);
		bj = b + 1;
		size_t numsteps = height1 - k;
		LIMIT_MAX(numsteps,bandwidth - 1);
		for (size_t j = 0; j < numsteps; j++) 
		{
			*b -= *bj * (*aj);
			aj++;
			bj++;
		}
	}

	return true;
};

bool FEMVirtMatrix<double>::solveSystem(double *B, double zero, int numthreads)
{
  if (!OK)
    return false;

  assert(numthreads > 0);
                              // data to pass to threads
  std::vector<ThreadData> data(numthreads);
  data[0].init(numthreads);
                              // create semaphore and threads
  if (numthreads > 1)
  {
                              // start threads
    for (int i = 0; i < numthreads; i++)
    {
                              // arguments to threads
      data[i].init(numthreads);
    }
  }
															// get first RHS vector
	double *b = B;
															// pointer to diagonal element of leading row
	double *upperrow = getElement(0,0);
															// go down
	size_t height1 = height - 1;
                              // offset to the next row in floats
	size_t shiftdown = width - 1;

	for (size_t k = 0; k < height; k++)
	{
		if (std::abs(*upperrow) < zero)
		{
      if (numthreads > 1)
      {
        for (int i = 0; i < numthreads; i++)
        {
          data[i].upperrow = nullptr;
        }
      }

			return false;
		}
															// make triangle
		size_t stepsright = ebandwidth >> 1;
		size_t numstepsdown = bandwidth - 1;
		if ((k + numstepsdown) > height1)
		{
			numstepsdown = height1 - k;
		}
                              // divide leading row by pivot
                              // inverted diagonal element for leading row
	  double diag = double(1.0) / *upperrow;
                                // divide leading row by its pivot
    register __m128d xdiag = _mm_set1_pd(diag);
    double *x = upperrow;
    for (size_t j = 0; j < stepsright; j++)
    {
      _mm_storeu_pd(x,_mm_mul_pd(_mm_loadu_pd(x),xdiag));
      x += 2; 
    }
                              // divide RHS by pivot
    *b *= diag;
                              // here we have :
                              // upperrow     - pointer to lead row
                              // b            - pointer to RHS for leading row
                              // stepsright   - active row length in XMM registers
                              // shiftdown    - offset to next row
                              // numstepsdown - number of rows below lead row
    if (numthreads <= 1)
    {
      data[0].upperrow = upperrow;
      data[0].b = b;
      data[0].stepsright = stepsright;
      data[0].shiftdown = shiftdown;
      data[0].numstepsdown = numstepsdown;
      data[0].rowoffset = 1;

      ThreadProc8(&data[0]);
    } else
    {
                              // divide numstepsdown rows between threads
      std::vector<int> ranges;
      getRanges<int>(static_cast<int>(numstepsdown),numthreads,ranges);
                              // prepare data for each thread and signal to start
      for (int i = 0; i < numthreads; i++)
      {
        data[i].upperrow = upperrow;
        data[i].b = b;
        data[i].stepsright = stepsright;
        data[i].shiftdown = shiftdown;
        data[i].numstepsdown = ranges[i + 1] - ranges[i];
        data[i].rowoffset = ranges[i] + 1;
      }

      // threads
      std::vector<std::thread> threads(numthreads);

                              // reset all events
      for (int i = 0; i < numthreads; i++)
      {
        // create threads
        threads[i] = std::thread(&ThreadProcDouble, &data[i]);
      }

      for (auto& t : threads)
      {
        t.join();
      }
    }
															// get next diagonal element
		upperrow += width;
															// get next RHS vector
		b++;

    if (progressprint && k % 100 == 0)
      printf("%zd\r",k);
	}

  if (progressprint)
    printf("\n");
                              // stop all threads
  if (numthreads > 0)
  {
                              // set data to stop all threads
    for (int i = 0; i < numthreads; i++)
    {
      data[i].upperrow = nullptr;
    }
  }
															// back substitution, not yet parallel
	b--;
	double *bj;
	double *aj;
	for (LINT k = (height - 2); k >= 0; k--)
	{
		b--;
		aj = getElement(k,k + 1);
		bj = b + 1;
		size_t numsteps = height1 - k;
		LIMIT_MAX(numsteps,bandwidth - 1);
		for (size_t j = 0; j < numsteps; j++) 
		{
			*b -= *bj * (*aj);
			aj++;
			bj++;
		}
	}


	return true;
};


                              // data for thread operation over matrix and vector
typedef struct ThreadVectorData {
  size_t numthreads;          // total number of threads
  void *matrix;               // pointer to VirtSquareMatrix<float> 
  void *vector;               // vector to multiply by (float *)
  void *result;               // result vector (float *)
  LINT row;                   // starting row number
  LINT numrows;               // num rows to multiply

  void init(size_t numThreads)
  {
    numthreads = numThreads;
    matrix = vector = result = NULL;
    row = numrows = 0;
  }
} ThreadVectorData;

															// multiply numrows of matrix by vector starting
                              // from row
void ThreadMultFloat(ThreadVectorData* data)
{

  if (data->numrows == 0)
    return;
                              // pointer to matrix
  FEMVirtMatrix<float> *matrix = (FEMVirtMatrix<float> *) data->matrix;
  float *vector = (float *) data->vector;
  float *result = (float *) data->result;

	LINT j1,j2;
	LINT height1 = static_cast<LINT>(matrix->height) - 1;
                              // starting result pointer
  float *rs = result + data->row;
                              // loop along rows
	for (LINT n = 0; n < data->numrows; n++)
	{
                              // actual row number
    LINT i = data->row + n;

		j1 = i - matrix->ebandwidth;
		LIMIT_MIN(j1,0);
		j2 = i + matrix->ebandwidth - 1;
		LIMIT_MAX(j2,height1);
		LINT numsteps = (j2 - j1 + 1) / 4;
		if ((j2 - j1 + 1) % 4) numsteps++;

		float *as = matrix->getElement(i,j1);
		float *vs = vector + j1;
		*rs = 0;

    register __m128 sum = _mm_set1_ps(0.0);
    float *a = as;
    float *v = vs;
    for (int j = 0; j < numsteps; j++)
    {
      register __m128 m = _mm_mul_ps(_mm_loadu_ps(a),_mm_loadu_ps(v));
                              // dot product
      register __m128 m1 = m;
      register __m128 m2 = m;
      register __m128 m3 = m;
      m1 = _mm_shuffle_ps(m1,m1,1);
      m2 = _mm_shuffle_ps(m2,m2,2);
      m3 = _mm_shuffle_ps(m3,m3,3);
                              // increment sum
      sum = _mm_add_ss(_mm_add_ss(_mm_add_ss(_mm_add_ss(sum,m),m1),m2),m3);

      a += 4; 
      v += 4; 
    }
                              // store
    _mm_store_ss(rs,sum);
                              // next XMM register
		rs++;
	}

  return;
}
															// multiply numrows of matrix by vector starting
                              // from row
void ThreadMultDouble(ThreadVectorData* data)
{
  if (data->numrows == 0)
    return;
                              // pointer to matrix
  FEMVirtMatrix<double> *matrix = (FEMVirtMatrix<double> *) data->matrix;
  double *vector = (double *) data->vector;
  double *result = (double *) data->result;

	LINT j1,j2;
	LINT height1 = static_cast<LINT>(matrix->height) - 1;
                              // starting result pointer
  double *rs = result + data->row;
                              // loop along rows
	for (LINT n = 0; n < data->numrows; n++)
	{
                              // actual row number
    LINT i = data->row + n;

		j1 = i - matrix->ebandwidth;
		LIMIT_MIN(j1,0);
		j2 = i + matrix->ebandwidth - 1;
		LIMIT_MAX(j2,height1);
		LINT numsteps = (j2 - j1 + 1) / 2;
		if ((j2 - j1 + 1) % 2) numsteps++;

		double *as = matrix->getElement(i,j1);
		double *vs = vector + j1;
		*rs = 0;

    register __m128d sum = _mm_set1_pd(0.0);
    double *a = as;
    double *v = vs;
    for (int j = 0; j < numsteps; j++)
    {
      register __m128d m = _mm_mul_pd(_mm_loadu_pd(a),_mm_loadu_pd(v));
                              // dot product
      register __m128d m1 = m;
      m1 = _mm_shuffle_pd(m1,m1,1);
                              // increment sum
      sum = _mm_add_sd(_mm_add_sd(sum,m),m1);

      a += 2;
      v += 2;
    }
                              // store
    _mm_store_sd(rs,sum);
                              // next XMM register
		rs++;
	}

  return;
}
															// multiply by vector V with the result in R;
															// vectors V and R MUST BE 4 FLOATS LONGER than the order
															// of the matrix, with the rest filled by zeroes
bool FEMVirtMatrix<float>::multiply(float *V, float *R, int numthreads)
{
  if (!OK)
    return false;

  if (numthreads <= 1)
  {
    ThreadVectorData data;
    data.init(numthreads);
    data.matrix = this;
    data.vector = V;
    data.result = R;
    data.row = 0;
    data.numrows = height;

    ThreadMultFloat(&data);
  } else
  {
                              // threads
    std::vector<thread> threads(numthreads);
                              // data to pass to threads
    std::vector<ThreadVectorData> data(numthreads);
                              // divide numstepsdown rows between threads
    std::vector<int> ranges;
    getRanges<int>(static_cast<int>(height),numthreads,ranges);
                              // prepare data for each thread and signal to start
    for (int i = 0; i < numthreads; i++)
    {
      data[i].matrix = this;
      data[i].vector = V;
      data[i].result = R;
      data[i].row = ranges[i];
      data[i].numrows = ranges[i + 1] - ranges[i];
      threads[i] = std::thread(&ThreadMultFloat, &data[i]);
    }
    // start thread
    for (auto& t : threads)
    {
      t.join();
    }
  }

  return true;
}
															// multiply by vector V with the result in R;
															// vectors V and R MUST BE 2 DOUBLES LONGER than the order
															// of the matrix, with the rest filled by zeroes
bool FEMVirtMatrix<double>::multiply(double *V, double *R, int numthreads)
{
  if (!OK)
    return false;

  if (numthreads <= 1)
  {
    ThreadVectorData data;
    data.init(numthreads);
    data.matrix = this;
    data.vector = V;
    data.result = R;
    data.row = 0;
    data.numrows = height;

    ThreadMultDouble(&data);
  } else
  {
                              // threads
    std::vector<thread> threads(numthreads);
                              // data to pass to threads
    std::vector<ThreadVectorData> data(numthreads);
                              // divide numstepsdown rows between threads
    std::vector<int> ranges;
    getRanges<int>(static_cast<int>(height),numthreads,ranges);
                              // prepare data for each thread and signal to start
    for (int i = 0; i < numthreads; i++)
    {
      data[i].matrix = this;
      data[i].vector = V;
      data[i].result = R;
      data[i].row = ranges[i];
      data[i].numrows = ranges[i + 1] - ranges[i];
                              // start worker threads
      threads[i] = std::thread(&ThreadMultDouble, &data[i]);
    }
                              // start threads
    for (auto& t : threads)
    {
      t.join();
    }
  }

  return true;
}

template<> void vectorsSubtract(float* v0, float* v1, size_t length)
{
  assert(length % 4 == 0);

  size_t numsteps = length / 4;
  for (size_t i = 0; i < numsteps; i++)
  {
    _mm_storeu_ps(v0, _mm_sub_ps(_mm_loadu_ps(v0), _mm_loadu_ps(v1)));
    v0 += 4;
    v1 += 4;
  }
}

template<> void vectorsSubtract(double* v0, double* v1, size_t length)
{
  assert(length % 2 == 0);

  size_t numsteps = length / 2;
  for (size_t i = 0; i < numsteps; i++)
  {
    _mm_storeu_pd(v0, _mm_sub_pd(_mm_loadu_pd(v0), _mm_loadu_pd(v1)));
    v0 += 2;
    v1 += 2;
  }
}
