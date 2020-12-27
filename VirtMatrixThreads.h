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

#include <assert.h>
#include "Allocator.h"

using namespace std;


/**
  Thread code for VirtSquareMatrix
*/

                                  // data for thread
typedef struct ThreadData {
  size_t numthreads;              // number of threads
  void *upperrow;                 // pointer to lead row
  void *b;                        // pointer to RHS for leading row
  size_t stepsright;              // active row length in XMM registers
  size_t shiftdown;               // offset to next row
  size_t numstepsdown;            // number of rows below lead row
  size_t rowoffset;               // offset down in rows from leading row

  void init(const size_t numThreads)
  {
    numthreads = numThreads;
    upperrow = b = NULL;
    stepsright = shiftdown = numstepsdown = 0;
    rowoffset = 1;
  }
} ThreadData;

                              // thread function
void ThreadProc4(ThreadData *data);
                              // thread function
void ThreadProc8(ThreadData *data);
                              // thread function
void ThreadProcFloat(ThreadData *data);
                              // thread function
void ThreadProcDouble(ThreadData *data);

