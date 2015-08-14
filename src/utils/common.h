/*  common.h -- misc utilities

    Copyright (c) 2015, The Griffith Lab

    Author: Avinash Ramu <aramu@genome.wustl.edu>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#ifndef COMMON_H_
#define COMMON_H_

#include <iostream>
#include <sstream>
#include <map>

using namespace std;

//Convert a number to a string
template <typename T>
string num_to_str(T num) {
    stringstream ss;
    ss << num;
    return ss.str();
}

//Convert a string to a number
template <typename T>
string num_to_str(string str) {
    stringstream ss(str);
    T num;
    ss << num;
    if(!ss.eof())
        throw "Unable to convert string to number!";
    return ss.str();
}

//Reverse complement short DNA seqs
inline string rev_comp(string s1) {
    string rc;
    for(int i = s1.length() - 1; i >= 0; i--) {
        char rc_char;
        switch(s1[i]) {
            case 'A':
                rc_char = 'T';
                break;
            case 'C':
                rc_char = 'G';
                break;
            case 'G':
                rc_char = 'C';
                break;
            case 'T':
                rc_char = 'A';
                break;
            default:
                rc_char = 'N';
                break;
        }
        rc.insert(rc.end(), rc_char);
    }
    return rc;
}

//Remove quotes from strings
inline void unquote(string & s1) {
    if(s1.empty())
        return;
    if(s1[0] == '"' && s1[s1.length() - 1] == '"') {
        s1.erase(s1.begin());
        s1.erase(s1.end() - 1);
    }
}

#endif

