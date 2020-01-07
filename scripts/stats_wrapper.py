import glob
import subprocess
import os
import argparse
import shutil

input_parser = argparse.ArgumentParser(
    description="Run RegTools stats script",
)
input_parser.add_argument(
    'tag',
    help="Variant tag parameter used to run RegTools.",
)

args = input_parser.parse_args()

tag = args.tag
cwd = os.getcwd()

lines_per_file = 25000
smallfile = None
with open(f'all_splicing_variants_{tag}.bed', 'r') as bigfile:
    header = bigfile.readline()
    for lineno, line in enumerate(bigfile):
        if lineno % lines_per_file == 0:
            if smallfile:
                smallfile.close()
            small_filename = 'small_file_{}.txt'.format(lineno + lines_per_file)
            smallfile = open(small_filename, "w")
            smallfile.write(header)
        smallfile.write(line)
    if smallfile:
        smallfile.close()
#get chunks
files = glob.glob('small_file_*')
files.sort()
number_of_in_files = len(files)
for file in files:
    subprocess.run(f'Rscript --vanilla /home/ec2-user/workspace/regtools/scripts/compare_junctions_hist_v2.R {tag} {file}', shell=True, check=True)
output_files = glob.glob("*_out.tsv")
output_files.sort()# glob lacks reliable ordering, so impose your own if output order matters
number_of_out_files = len(output_files)
if number_of_in_files == number_of_out_files:
    with open(f'compare_junctions/hist/junction_pvalues_{tag}.tsv', 'wb') as outfile:
        for i, fname in enumerate(output_files):
            with open(fname, 'rb') as infile:
                if i != 0:
                    infile.readline()  # Throw away header on all but first file
                # Block copy rest of file from input to output without parsing
                shutil.copyfileobj(infile, outfile)
                print(fname + " has been imported.")
else:
    print("Number of output files doesn't match the number of input files that should have been processed")
files = glob.glob('small_file_*')
for file in files:
     os.remove(file)

