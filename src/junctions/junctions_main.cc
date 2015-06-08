#include <iostream>

using namespace std;

int junctions_usage() {
    cout << "\nUsage:\t\t" << "regtools junctions <command> [options]";
    cout << "\nCommand:\t" << "annotate\tAnnotate the junctions.";
    cout << "\n";
    return 0;
}

int junctions_main(int argc, char* argv[]) {
    cout << "in junctions";
    if(argc > 1) {
        string subcmd(argv[1]);
        if(subcmd == "annotate") {
            return junctions_main(argc - 1, argv + 1);
        }
    }
    return junctions_usage();
}
