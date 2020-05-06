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

target_lines_per_file = 25000
lines_per_file = 0
input_file = f'all_splicing_variants_{tag}.bed'
lines = open(input_file).readlines()
count = len(lines)
if count <= lines_per_file:
    subprocess.run(f'Rscript --vanilla /home/ec2-user/workspace/regtools/scripts/compare_junctions_hist_v2.R {tag} {input_file}')
else:
    header = lines[0]
    lines.pop(0)
    lines.sort()
    filenum = 1
    small_filename = f'small_file_{filenum}.txt'
    smallfile = open(small_filename, "w")
    smallfile.write(header)
    lines_per_file += target_lines_per_file
    for lineno, line in enumerate(lines):
        smallfile.write(line)
        if lineno >= lines_per_file:
            fields1 = line.split('\t')
            variant1 = f'{fields1[0]}_{fields1[1]}_{fields1[2]}'
            fields2 = lines[lineno+1].split('\t')
            variant2 = f'{fields2[0]}_{fields2[1]}_{fields2[2]}'
            if variant1 != variant2:
                smallfile.close()
                filenum += 1
                small_filename = f'small_file_{filenum}.txt'
                smallfile = open(small_filename, "w")
                smallfile.write(header)
                lines_per_file += target_lines_per_file
# get chunks
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
