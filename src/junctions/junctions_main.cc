/*  junctions_main.cc -- handle the 'junction' commands

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
#include "common.h"
#include "gtf_parser.h"
#include "junctions_annotator.h"
#include "junctions_extractor.h"

using namespace std;

//Usage for junctions subcommands
int junctions_usage(ostream &out = cout) {
    out << "\nUsage:\t\t" << "regtools junctions <command> [options]";
    out << "\nCommand:\t" << "extract\t\tIdentify exon-exon junctions from alignments.";
    out << "\n\t\tannotate\tAnnotate the junctions.";
    out << "\n";
    return 0;
}


//Run 'junctions extract'
int junctions_extract(int argc, char *argv[]) {
    JunctionsExtractor extract;
    try {
        extract.parse_options(argc, argv);
        extract.identify_junctions_from_BAM();
        extract.print_all_junctions();
    } catch(const cmdline_help_exception& e) {
        cerr << e.what();
        return 0;
    } catch(const runtime_error& error) {
        cerr << error.what();
        extract.usage();
        return 1;
    }
    return 0;
}

//Run 'junctions annotate' subcommand
int junctions_annotate(int argc, char *argv[]) {
    JunctionsAnnotator anno;
    AnnotatedJunction line;
    line.reset();
    int linec = 0;
    ofstream out;
    try {
        anno.parse_options(argc, argv);
        anno.read_gtf();
        anno.open_junctions();
        anno.set_ofstream_object(out);
        line.print_header(out);
        while(anno.get_single_junction(line)) {
            anno.get_splice_site(line);
            anno.annotate_junction_with_gtf(line);
            line.print(out);
            line.reset();
            linec++;
        }
        anno.close_ofstream();
        cerr << endl << "Annotated " << linec << " lines.";
        anno.close_junctions();
    } catch(const cmdline_help_exception& e) {
        cerr << e.what();
        return 0;
    } catch(const runtime_error& e) {
        cerr << endl << e.what() << endl;
        return 1;
    }
    return 0;
}

//Parse out subcommands under junctions
int junctions_main(int argc, char *argv[]) {
    if(argc > 1) {
        string subcmd(argv[1]);
        if(subcmd == "extract") {
            return junctions_extract(argc - 1, argv + 1);
        }
        if(subcmd == "annotate") {
            return junctions_annotate(argc - 1, argv + 1);
        }
    }
    return junctions_usage();
}

