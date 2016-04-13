/*  gtf_utils.cc utility functions related to gtf files

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

#include <vector>
#include "bedFile.h"

using namespace std;

//Return True if variant within a certain window from the transcript
bool is_variant_within_transcript_window(const vector<BED> &exons, uint32_t pos,
                                                string transcript_strand,
                                                uint32_t window_size) {
    int n_exons = exons.size();
    if(transcript_strand == "+") {
        //variant inside transcript
        if(pos >= exons[0].start && pos <= exons[n_exons - 1].end) {
            return true;
        }
        //variant outside transcript
        if(exons[0].start - pos <= window_size &&
           exons[n_exons -1].start > pos) {
            return true;
        }
        //variant outside transcript
        if(pos - exons[n_exons - 1].end <= window_size &&
           exons[0].end < pos) {
            return true;
        }
    } else if(transcript_strand == "-") {
        //variant inside transcript
        if(pos >= exons[n_exons - 1].start && pos <= exons[0].end) {
            return true;
        }
        //variant outside transcript
        if(pos - exons[0].end <= window_size &&
           exons[n_exons -1].end < pos) {
            return true;
        }
        //variant outside transcript
        if(exons[n_exons - 1].start - pos <= window_size &&
           exons[0].start > pos) {
            return true;
        }
    } else {
        throw runtime_error("Unknown transcript strand.");
    }
    return false;
}
