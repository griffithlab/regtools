import subprocess
import os
import glob
import shutil
import sys

def run(cmd: str) -> None:
    subprocess.run(cmd, shell=True, check=True, stdout=sys.stdout)

cohorts = ['SKCM', 'GBM', 'READ', 'ESCA', 'PAAD', 'SARC',
          'OV', 'KIRP', 'CESC', 'KIRC', 'LIHC', 'STAD', 'BLCA',
          'COAD', 'LUSC', 'HNSC', 'LGG', 'LUAD', 'UCEC', 'BRCA']]

for cohort in cohorts:
    os.makedirs(f'{cohort}/samples', exist_ok=True)
    os.makedirs(f'{cohort}/compare_junctions/hist', exist_ok=True)
    bed_files = glob.glob(f'{cohort}*modified.bed')
    for bed_file in bed_files:
        shutil.copy(bed_file, f'{cohort}/{bed_file}')
    os.chdir(f'{cohort}/samples')
    run(f'aws s3 cp s3://regtools-results-unstranded/{cohort}/ . --recursive')
    tar_files = glob.glob('*.tar.gz')
    for tar_file in tar_files:
        run(f'tar xzf {tar_file}')
        os.remove(tar_file)
        run('rm -rf all*; rm -rf compare_junctions*')
    os.chdir('..')
    run('ls samples/ > dir_names.tsv')
    bed_files = glob.glob(f'*modified.bed')
    for bed_file in bed_files:
        tag = bed_file.split('_')[1]
        os.rename(bed_file, f'all_splicing_variants_{tag}.bed')
        run(f'python3 /home/ec2-user/workspace/regtools/scripts/stats_wrapper.py {tag}')
        run(f'Rscript --vanilla /home/ec2-user/workspace/regtools/scripts/filter_and_BH.R {tag}')
	run(f'aws s3 cp compare_junctions/ s3://regtools-results-unstranded/{cohort}/compare_junctions3/ --recursive)')
	os.chdir('..')
    shutil.rmtree(cohort)
