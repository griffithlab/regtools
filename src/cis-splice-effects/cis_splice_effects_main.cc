/*  cis_splice_effects_main.cc -- handle the 'cis-splice-effects' commands

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

#include <iostream>
#include <getopt.h>
#include <stdexcept>

using namespace std;

//Main command for identify
int cis_splice_effects_identify() {
    cerr << "in identify";
    return 0;
}

//Usage for cis-splice-effects subcommands
int cis_splice_effects_usage(ostream &out = cout) {
    out << "\nUsage:\t\t" << "regtools cis-splice-effects <command> [options]";
    out << "\n";
    return 0;
}

//Main command for cis-splice-effects
int cis_splice_effects_main(int argc, char* argv[]) {
    cerr << endl << "in splice effects model";
    if(string(argv[0]) == "identify") {
        cis_splice_effects_identify();
    }
    return 0;
}
