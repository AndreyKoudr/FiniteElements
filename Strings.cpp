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
#include "Strings.h"
#include "Types.h"

#include <algorithm>
#include <assert.h>

std::string to_string(double d, int len, int len_after_dot)
{
  assert(len > len_after_dot);

  char buf[1024];
  std::string format = std::string("%") + to_string(len) + '.' + to_string(len_after_dot) + 'f';
  sprintf_s(buf,sizeof(buf),format.c_str(),d);
  return std::string(buf);
}

std::string to_string(int i, int numdigits)
{
  std::string s = to_string(i);
  std::string res = padFromLeft(s,numdigits,'0');
  return res;
}

std::string padFromLeft(std::string s, unsigned int len, char ch)
{
  int n = len - static_cast<int>(s.length());
  if (n > 0) s = s.insert(0,n,ch);

  return s;
}

std::string upCase(std::string str)
{
  std::transform(str.begin(),str.end(),str.begin(),::toupper);
  return str;
}

std::string lowerCase(std::string str)
{
  std::transform(str.begin(),str.end(),str.begin(),::tolower);
  return str;
}

std::string trimLeft(std::string str, std::string chars)
{
  str.erase(0,str.find_first_not_of(chars));
  return str;
}

std::string trimRight(std::string str, std::string chars)
{
  str.erase(str.find_last_not_of(chars)+1);
  return str;
}

std::string trim(std::string str, std::string chars)
{
  return trimRight(trimLeft(str,chars),chars);
}

std::string replace(std::string str, char from, char to)
{
  std::string s = str;
  std::replace(s.begin(),s.end(),from,to);
  return s;
}

int find(std::string str, char ch)
{
  size_t pos = str.find(ch);
  if (pos == std::string::npos)
  {
    return NOT_FOUND;
  } else
  {
    return static_cast<int>(pos);
  }
}

int find(std::string str, std::string substr)
{
  size_t pos = str.find(substr);
  if (pos == std::string::npos)
  {
    return NOT_FOUND;
  } else
  {
    return static_cast<int>(pos);
  }
}

void deleteChars(std::string &str, int pos0, int numchars)
{
  str.erase(pos0,numchars);
}

std::string getSubString(std::string str, int pos0, int pos1)
{
  std::string res = str.substr(pos0,(pos1 - pos0 + 1));
  return res;
}

int parseWords(std::string str, char divider, int pos1[], int pos2[], int maxcount)
{
  int numwords = 0;
  int len = static_cast<int>(str.length());
  if (len > 0)
  {
		int curpos1 = NOT_DEFINED;
		if (str[0] != divider) curpos1 = 0;
		for (int i = 0; i < len; i++)
		{
			if (str[i] == divider)
      {
				if (curpos1 != NOT_DEFINED)
				{
          if (numwords < maxcount)
          {
            pos1[numwords] = curpos1;
            pos2[numwords] = i - 1;
					  curpos1 = NOT_DEFINED;
					  numwords++;
          } else
          {
            break;
          }
				}
			} else
			{
				if (curpos1 == NOT_DEFINED)
				{
          curpos1 = i;
				}
			}
		}

		if (curpos1 != NOT_DEFINED)
		{
      if (numwords < maxcount)
      {
        pos1[numwords] = curpos1;
        pos2[numwords] = len;
			  numwords++;
      }
		}
  }

  return numwords;
}

std::string forceExtension(std::string s, std::string ext)
{
  std::size_t pos = s.rfind('.');
	if (pos != std::string::npos)
	{
		if (ext.length())
		{
			return s.substr(0,pos + 1) + ext;
		} else
		{
			return s.substr(0,pos);
		}
	} else
	{
		if (ext.length())
		{
			return s + "." + ext;
		} else
		{
			return s;
		}
	}
}

std::string getExtension(std::string s)
{
  std::string ext;
  std::size_t pos = s.rfind('.');
	if (pos != std::string::npos)
	{
    ext = s.substr(pos + 1);
	}
  return ext;
}

std::string addBackslash(std::string directory)
{
  if (directory.length() == 0)
  {
    return std::string("");
  } else 
  {
#ifdef _WIN32
    if (directory.back() != '\\')
    {
      return directory + "\\";
#else
    if (directory.back() != '/')
    {
      return directory + "/";
#endif
    } else
    {
      return directory;
    }
  }
}

std::string justFileName(const std::string &path)
{
  std::string name = path;
  name.erase(0,path.find_last_of("\\/") + 1);

	return name;
}

std::string justDirectory(const std::string &path)
{
  std::string name = path;
  name.erase(path.find_last_of("\\/") + 1);

	return name;
}

bool containsOnlyChars(const std::string &s, char from, char to)
{
  for (size_t i = 0; i < s.length(); i++)
  {
    if (s[i] < from || s[i] > to) return false;
  }

  return true;
}

