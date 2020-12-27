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

#include <assert.h>
                              // class header
#include "Allocator.h"

using namespace std;

                              // report progress?
extern bool progressprint;

bool Allocator::init(const char* pname, const size_t psize, const bool clearMemory)
{
  assert(buffer == nullptr);
                              // on x64, buffer is 16-byte aligned by default
  buffer = (unsigned char*) malloc(psize);
                              // just in case check
  assert((reinterpret_cast<size_t>(buffer) & 0x000000000000000F) == 0);

  if (buffer != nullptr)
  {
    size = psize;

    if (clearMemory)
      memset(buffer, 0, size);

    return true;
  } else
  {
    return false;
  }
}

Allocator::Allocator(Allocator&& other)
{
  buffer = other.buffer;
  size = other.size;
  storedFileName_ = other.storedFileName_;

  other.buffer = nullptr;
  other.size = 0;
  other.storedFileName_.clear();
}

Allocator& Allocator::operator=(Allocator&& other)
{
  if (this != &other)
  {
    buffer = other.buffer;
    size = other.size;
    storedFileName_ = other.storedFileName_;

    other.buffer = nullptr;
    other.size = 0;
    other.storedFileName_.clear();
  }

  return *this;
}

Allocator::~Allocator()
{
  freeStoredCopy();

  if (buffer != nullptr)
  {
    free(buffer);
    buffer = nullptr;
    size = 0;
  }
}

/** File exists? */
static bool FileExists(const string& path)
{
  FILE* file = nullptr;

  if (fopen_s(&file,path.c_str(), "rb") == 0)
  {
    fclose(file);
    return true;
  }
  else
  {
    return false;
  }
}

void Allocator::freeStoredCopy()
{
  if (storedFileName_.length() && FileExists(storedFileName_))
  {
    if (progressprint)
      printf("Allocator stored file deleted : %s\n",storedFileName_.c_str());

    remove(storedFileName_.c_str());
  }
}

bool Allocator::storeCopy(const std::string &tempFileName)
{
  assert(buffer != NULL);
  if (buffer == NULL)
    return false;
														  // free stored matrix, only one copy can be stored				
  freeStoredCopy();
														  // allocate
  storedFileName_ = tempFileName;

  if (progressprint)
    printf("Allocator stored file : %s, size %zd\n",storedFileName_.c_str(), size);
														  // save to file
  FILE *fp = nullptr;
                              // save by portions of 16M
  if (fopen_s(&fp,storedFileName_.c_str(),"wb") == 0)
  {
    size_t portion = 16777216;
    size_t left = size;

    while (left > 0)
    {
      size_t writesize = (left > portion) ? portion : left;
      size_t written = fwrite(buffer + size - left,1,writesize,fp);
      if (written != writesize)
      {
        fclose(fp);
        return false;
      }
      left -= written;
    }

    fclose(fp);
  }

  return true;
}

/** Get file size */
static size_t fileSize(string &path)
{
  struct stat st;
  stat(path.c_str(), &st);
  size_t size = st.st_size;
  return size;
}

bool Allocator::restoreCopy(bool freestored)
{
  if (storedFileName_.length() && FileExists(storedFileName_))
  {
                            // get stored file size
    size_t filesize = fileSize(storedFileName_);

    if (filesize == size)
    {
                            // load from file
      FILE *fp = nullptr;

      if (fopen_s(&fp,storedFileName_.c_str(),"rb") == 0)
      {
                            // read by portions of 16M
        size_t portion = 16777216;
        size_t left = size;

        while (left > 0)
        {
          size_t readsize = (left > portion) ? portion : left;
          size_t read = fread(buffer + size - left,1,readsize,fp);
          if (read != readsize)
          {
            fclose(fp);
            return false;
          }   
          left -= read;
        }

        fclose(fp);
      }
													  // free stored matrix				
	    if (freestored) 
        freeStoredCopy();

      return true;
    } else
    {
      return false;
    }
  } else
  {
	  return false;
  }
}

