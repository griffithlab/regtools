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
#include <map>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>
#include "stdint.h"
#include "bedFile.h"
#include "htslib/hts.h"
#include "htslib/tbx.h"

using namespace std;

namespace common {
    //Convert a number to a string
    template <typename T>
        string num_to_str(T num) {
            stringstream ss;
            ss << num;
            return ss.str();
    }

    //Convert a number to a string
    inline uint32_t str_to_num(string num) {
            stringstream ss;
            uint32_t num_uint;
            ss << num;
            ss >> num_uint;
            return num_uint;
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

    //Define cmdline_help_exception - Thanks tabbott!
    class cmdline_help_exception : public std::runtime_error {
        public:
            cmdline_help_exception(std::string const& msg)
                : std::runtime_error(msg) {
                }
    };

    //Check if file exists
    inline bool file_exists(const std::string& file) {
        struct stat buf1;
        return (stat(file.c_str(), &buf1) == 0);
    }

    //Difference in CHRPOS coordinates
    inline uint32_t coordinate_diff(CHRPOS pos1, CHRPOS pos2) {
        if(pos1 > pos2)
            return pos1 - pos2;
        else
            return pos2 - pos1;
    }

    //Copy one stream object into another
    inline void copy_stream(const ostream &source,
                         ostream &dest) {
            dest.copyfmt(source);
            dest.basic_ios<char>::rdbuf(source.rdbuf());
            dest.clear(source.rdstate());
    }

    //Create a region string using chr, start, end
    //this is of the form chr:start-end
    inline std::string create_region_string(const char* chr,
                                        int start, int end) {
        stringstream ss1;
        ss1 << chr << ":" << start << "-" << end;
        return ss1.str();
    }

    //Check if tabix index exists
    //Throws runtime_error if index does not exist
    inline bool check_tabix_index(string file) {
        htsFile *fp = hts_open(file.c_str(), "rb");
        if(!fp) {
            std::cerr << "Unable to open " << file;
            throw runtime_error("Unable to open file.");
        }
        tbx_t *idx = tbx_index_load(fp->fn);
        if(!idx) {
            stringstream ss;
            ss << "Unable to open tabix index for " << file << endl;
            throw runtime_error(ss.str());
        }
        hts_close(fp);
        tbx_destroy(idx);
        return true;
    }
}

#endif
