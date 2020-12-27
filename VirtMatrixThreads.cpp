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
                              // class header
#include "VirtMatrixThreads.h"
                              // intrinsics
#include <xmmintrin.h>
#include <emmintrin.h>
#include <assert.h>

using namespace std;

                              // thread function
void ThreadProc4(ThreadData *data)
{
  float *upperrow = (float *) data->upperrow;
  assert(upperrow != nullptr);
                              // this element should be 1.0
  assert(std::abs(*upperrow - 1.0) < 0.000001);
                              // nothing to do
  if (data->numstepsdown == 0)
    return;

  float *b = (float *) data->b;
                              // go down below lead row excluding one unknown
  float *row = upperrow + data->shiftdown * data->rowoffset;
  float *b1 = b + data->rowoffset;
  for (size_t i = 0; i < data->numstepsdown; i++)
  {
    float *x = upperrow;
    float *r = row;
                              // RHS
    *b1 -= (*b) * (*row);
                              // first row element
    register __m128 xfirst = _mm_set1_ps(*row);
    for (size_t j = 0; j < data->stepsright; j++)
    {     
      _mm_storeu_ps(r,_mm_sub_ps(_mm_loadu_ps(r),_mm_mul_ps(_mm_loadu_ps(x),xfirst)));
      x += 4; 
      r += 4; 
    }
                              // this element should be zero
//    assert(std::abs(*row) < 0.000001);
                              // to the start of next row
    row += data->shiftdown;
                              // ..and next RHS
    b1++;
  }

  return;
}
                              // thread function
void ThreadProc8(ThreadData *data)
{
  double *upperrow = (double *) data->upperrow;
  assert(upperrow != nullptr);
                              // this element should be 1.0
  assert(std::abs(*upperrow - 1.0) < 0.000001);
                              // nothing to do
  if (data->numstepsdown == 0)
    return;

  double *b = (double *) data->b;
                              // go down below lead row excluding one unknown
  double *row = upperrow + data->shiftdown * data->rowoffset; //!!! rowoffset should be in doubles
  double *b1 = b + data->rowoffset;
  for (size_t i = 0; i < data->numstepsdown; i++) 
  {
    double *x = upperrow;
    double *r = row;
                              // RHS
    *b1 -= (*b) * (*row);
                              // first row element
    register __m128d xfirst = _mm_set1_pd(*row);
    for (size_t j = 0; j < data->stepsright; j++)
    {     
      _mm_storeu_pd(r,_mm_sub_pd(_mm_loadu_pd(r),_mm_mul_pd(_mm_loadu_pd(x),xfirst)));
      x += 2; 
      r += 2; 
    }
                              // this element should be zero
//    assert(std::abs(*row) < 0.000001);
                              // to the start of next row
    row += data->shiftdown;   //!!! shiftdown should be in doubles
                              // ..and next RHS
    b1++;
  }

  return;
}
                              // thread function
void ThreadProcFloat(ThreadData *data)
{
  if (data->upperrow != nullptr)
  {
                            // do the job
    ThreadProc4(data);
  }

  return;
}
                              // thread function
void ThreadProcDouble(ThreadData *data)
{
  if (data->upperrow != nullptr)
  {
    // do the job
    ThreadProc8(data);
  }

  return;
}

