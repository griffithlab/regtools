/*  variants_main.cc -- handle the 'variants' commands

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
#include <stdexcept>
#include "common.h"
#include "variants_annotator.h"

using namespace std;

//Usage for variants subcommands
int variants_usage(ostream &out = cout) {
    out << "\nUsage:\t\t" << "regtools variants <command> [options]";
    out << "\nCommand:\t" << "annotate\t\tAnnotate variants with splicing information.";
    out << "\n";
    return 0;
}

//Run 'variants annotate' subcommand
int variants_annotate(int argc, char *argv[]) {
    VariantsAnnotator va;
    try {
        va.parse_options(argc, argv);
        va.annotate_vcf();
    } catch(const common::cmdline_help_exception& e) {
        cerr << e.what();
        return 0;
    } catch (runtime_error e) {
        cerr << e.what();
        return 1;
    }
    return 0;
}

//Parse out subcommands under variants
int variants_main(int argc, char *argv[]) {
    if(argc > 1) {
        string subcmd(argv[1]);
        if(subcmd == "annotate") {
            return variants_annotate(argc - 1, argv + 1);
        }
    }
    return variants_usage();
}
