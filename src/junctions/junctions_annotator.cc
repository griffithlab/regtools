/*  junctions_annotator.cc -- JunctionsAnnotator class

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

#include <getopt.h>
#include "junctions_annotator.h"
#include "faidx.h"

int JunctionsAnnotator::usage() {
    cout << "\nUsage:\t\t" << "regtools junctions annotate [options] -r ref.fa -j junctions.bed -g annotations.gtf";
    cout << "\nCommand:\t" << "annotate\tAnnotate the junctions.";
    cout << "\n";
    return 0;
}

int JunctionsAnnotator::parse_options(int argc, char *argv[]) {
    static struct option long_options[] = {
        {"ref", required_argument, 0, 'r'},
        {"junctions", required_argument, 0, 'j'},
        {"gtf", required_argument, 0, 'g'},
    };
    int option_index = 0;
    int c = getopt_long(argc, argv, "r:g:j:",
                    long_options, &option_index);
    while(c != -1) {
        switch(c) {
            case 'r':
                ref_ = optarg;
                break;
            case 'g':
                gtf_ = optarg;
                break;
            case 'j':
                junctions.bedFile = optarg;
                break;
            default:
                usage();
                exit(-1);
        }
        c = getopt_long(argc, argv, "r:g:j:",
                    long_options, &option_index);
    }
    if(optind < argc || ref_.empty() || junctions.bedFile.empty() || gtf_.empty()) {
        usage();
        exit(-1);
    }
    cout << "\nReference: " << ref_;
    cout << "\nGTF: " << gtf_;
    cout << "\nJunctions: " << junctions.bedFile;
}

string JunctionsAnnotator::get_reference_sequence(string position) {
    int len;
    faidx_t *fai = fai_load(ref_.c_str());
    char *s = fai_fetch(fai, position.c_str(), &len);
    std::string seq(s);
    free(s);
    fai_destroy(fai);
    return seq;
}

