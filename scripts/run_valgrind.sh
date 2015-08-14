set -o xtrace
rm -f regtools_valgrind.op; valgrind --tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=20 --track-fds=yes ./regtools junctions create ../tests/test.rnaseq.bam 2>regtools_valgrind.op
