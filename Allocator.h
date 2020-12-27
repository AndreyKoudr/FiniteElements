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
                              // server is defined for char * only
#undef UNICODE

#include <string>

/**
    Basic allocator. This class allocates a big 16-byte aligned buffer for a matrix. 
  In case of failure, buffer remains nullptr.
    The buffer contents can also be stored in a temporary file by storeCopy() and restored
  by restoreCopy().
*/

class Allocator {
protected:

	/** File name for data copy stored in file */
  std::string storedFileName_;

	/** Free stored copy (delete file) */
	void freeStoredCopy();

public:

  /** Allocated buffer, no access via method like buffer() for speed */
  unsigned char* buffer = nullptr;

  /** Buffer size in bytes, no access via method like size() for speed */
  size_t size = 0;


  /** Default constructor */
  Allocator() = default;

  /** Initialise */
  virtual bool init(const char* pname, const size_t psize, const bool clearMemory);

  /** Move constructor is only allowed */
  Allocator(Allocator&& other);

  /** Move assignment is only allowed */
  Allocator& operator=(Allocator&& other);

  /** No copy constructors */
  Allocator(const Allocator&) = delete;

  /** No assignments */
  Allocator& operator=(const Allocator&) = delete;

  /** Destructor, must NOT be virtual */
  ~Allocator();

  /** Store copy of data in temp file. The file directory must be writable and 
    only one copy can be stored. */
  bool storeCopy(const std::string &tempFileName);

  /** Restore copy from file; delete file after restore if freeStored is specified */
  bool restoreCopy(const bool freeStored);
};

