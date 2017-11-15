/*  cis_ase_main.cc -- handle the 'cis-ase' commands

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

//Deal with the limits in htslib
#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
#endif

#include <iostream>
#include <getopt.h>
#include <stdexcept>
#include "cis_ase_identifier.h"
#include "common.h"

using namespace std;

//Main command for identify
int cis_ase_identify(int argc, char* argv[]) {
    CisAseIdentifier cai1;
    try {
        cai1.parse_options(argc, argv);
        cai1.run();
    } catch(const common::cmdline_help_exception& e) {
        cerr << e.what();
        return 0;
    } catch (const std::runtime_error &e) {
        cerr << e.what() << endl;
        return 1;
    } catch (const std::logic_error &e) {
        cerr << e.what() << endl;
        return 1;
    }
    return 0;
}

//Usage for cis-ase subcommands
int cis_ase_usage(ostream &out = cout) {
    out << "Usage:\t\t" << "regtools cis-ase <command> [options]" << endl;
    out << "Command:\t" << "identify\t\tIdentify cis ase." << endl;
    out << endl;
    return 0;
}

//Main command for cis-ase
int cis_ase_main(int argc, char* argv[]) {
    if(argc < 2) {
        return cis_ase_usage(std::cout);
    }
    if(string(argv[1]) == "identify") {
        return cis_ase_identify(argc - 1, argv + 1);
    }
    return cis_ase_usage(std::cout);
    return 0;
}
