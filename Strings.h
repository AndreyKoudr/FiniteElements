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
                              
#include <string>

using namespace std;

/**
  Some basic string operations needed for this project.

  Why not ...(const string &s... ?

  Example :

  (1)
  std::string upCase(std::string str)
  {
    std::transform(str.begin(),str.end(),str.begin(),::toupper);
    return str;
  }

  (2)
  std::string upCase(const std::string& str)
  {
    string s1 = str;
    std::transform(s1.begin(),s1.end(),s1.begin(),::toupper);
    return s1;
  }

  (1) : 2 constructors (parameter by value + on return); operations on stack variable str are fast
  (2) : 2 constructors (assignment + on return); operations on stack variable s1 are fast

*/

                              // double to string
string to_string(double d, int len, int len_after_dot);
                              // get string of specified length
                              // with leading zeroes
string to_string(int i, int numdigits);

                              // pad string from left by chars
string padFromLeft(string s, unsigned int len, char ch = ' ');

                              // convert to upper case
string upCase(string str);
                              // convert to lower case
string lowerCase(string str);



                              // delete characters (" \n\r\t" in the worst case)
                              // from string head
string trimLeft(string str, string chars = " ");
                              // delete characters (" \n\r\t" in the worst case)
                              // from string tail
string trimRight(string str, string chars = " ");
                              // delete characters (" \n\r\t" in the worst case)
                              // from both string tail and head
string trim(string str, string chars = " ");
                              // replace characters
string replace(string str, char from, char to);
                              // delete numchars characters from pos0
void deleteChars(std::string &str, int pos0, int numchars);

                              // find first in string, return position,
                              // -1 if not found
int find(string str, char ch);
int find(string str, string substr);
                              // get substring form pos0 to pos1
string getSubString(string str, int pos0, int pos1);
                              // parse string by divider; returns number of extracted words;
                              // pos1 and pos2 on exit contain starting and ending positions
                              // of each word; maxcount is dimension of pos1,pos2 arrays;
                              // positions are 0-based (?)
int parseWords(string str, char divider, int pos1[], int pos2[], int maxcount);

															// force extension
string forceExtension(string s, string ext);
                              // get extension if any, otherwise returns empty string
string getExtension(string s);
                              // add backslash to directory name
string addBackslash(string directory);
                              // just file name
std::string justFileName(const std::string &path);
                              // just directory
std::string justDirectory(const std::string &path);

                              // string contains only chars from to
bool containsOnlyChars(const string &s, char from, char to);

