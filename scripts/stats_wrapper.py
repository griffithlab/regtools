import glob
import subprocess
import os

tags = ['default', 'i50e5', 'E', 'I']
cwd = os.getcwd()


for tag in tags:
    lines_per_file = 50000
    smallfile = None
    with open(f'all_splicing_variants_{tag}.tsv', 'r') as bigfile:
        for lineno, line in enumerate(bigfile):
            if lineno % lines_per_file == 0:
                if smallfile:
                    smallfile.close()
                small_filename = 'small_file_{}.txt'.format(lineno + lines_per_file)
                smallfile = open(small_filename, "w")
            smallfile.write(line)
        if smallfile:
            smallfile.close()
    #get chunks
    files = glob.glob('small_file_*')
    for file in files:
        subprocess.run(f'Rscript --vanilla /home/ec2-user/workspace/data/compare_junctions_hist_v2.R {tag} {file}')
    subprocess.run(f"awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; } 1 {print} ' small_file*.txt > junction_pvalues_{tag}.tsv")
    os.remove('small_file*')
