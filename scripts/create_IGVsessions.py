import yaml
import subprocess
import glob
import csv
import sys
import os
import requests
import shutil
import json
from pathlib import Path
import argparse

uid = os.getuid()
gid = os.getgid()
cwd = os.getcwd()


input_parser = argparse.ArgumentParser(
    description="Create IGV sessions for TCGA data",
)
input_parser.add_argument(
    'yaml file',
    help="Yaml file with inputs",
)

def run(cmd):
    subprocess.run(cmd, shell=True, check=True, stdout=sys.stdout)

def read_file(filename):
    with open(filename) as f:
        return f.read()

def run_docker_image_as_current_user(image_name_and_all_args: str, stdout=None):
    docker_args = image_name_and_all_args.split()
    print(f'Running docker container command "{docker_args[0]} {docker_args[1]}"')
    run(f'docker run --user {uid}:{gid} -v {cwd}:/manifests {image_name_and_all_args}')

def get_bam(sample, token_file, chrom, bam_window_start, bam_window_end, variant_junction):
    print('Obtaining bam manifest file')
    bam_url = 'https://api.gdc.cancer.gov/files?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_format%22%2C%22value%22%3A%5B%22BAM%22%5D%7D%7D%2C%7B%22op%22%3A%22AND%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.experimental_strategy%22%2C%22value%22%3A%5B%22RNA-Seq%22%5D%7D%7D%2C%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22IN%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.program.name%22%2C%22value%22%3A%5B%22TCGA%22%5D%7D%7D%2C%7B%22op%22%3A%22IN%22%2C%22content%22%3A%7B%22field%22%3A%22cases.samples.submitter_id%22%2C%22value%22%3A%5B%22{0}%22%5D%7D%7D%5D%7D%5D%7D%5D%7D%5D%7D&query=files.data_format%20in%20%5B%22BAM%22%5D%20and%20files.experimental_strategy%20in%20%5B%22RNA-Seq%22%5D%20AND%20cases.project.program.name%20IN%20%5BTCGA%5D%20and%20cases.samples.submitter_id%20IN%20%5B%22P{0}&return_type=manifest'.format(sample)
    response = requests.get(bam_url)
    if response.status_code != 200:
        raise Exception(sample, f'Failed to download bam manifest file: HTTP Status Code: {response.status_code}')
    manifest_contents = response.content.decode('utf-8').strip()
    if not manifest_contents.startswith('id\t'):
        raise Exception(sample, f'bam manifest did not have expected format')

    bam_manifest_path = Path(f'{sample}_bam_manifest.txt')
    bam_manifest_path.write_text(manifest_contents)
    # run_docker_image_as_current_user(
    #     f'mgibio/gdc-client gdc-client download -m /manifests/{bam_manifest_path} -t /manifests/gdc-user-token.txt -d /manifests/igv_session')
    run(f'aws s3 cp {token_file} .')
    with open(bam_manifest_path, 'r') as manifest_file:
        reader = csv.reader(manifest_file, delimiter='\t')
        next(manifest_file)
        for line in reader:
            bam_id = line[0]
            bam_file = f'{sample}_{variant_junction}.bam'
            token_file = "gdc-user-token.txt"

            file_id = f'{bam_id}'

            data_endpt = "https://api.gdc.cancer.gov/slicing/view/{}".format(file_id)

            with open(token_file, "r") as token:
                token_string = str(token.read().strip())

            params = {"region": [f'{chrom}:{bam_window_start}-{bam_window_end}']}

            response = requests.post(data_endpt,
                                     data=json.dumps(params),
                                     headers={
                                         "Content-Type": "application/json",
                                         "X-Auth-Token": token_string
                                     })


            with open(bam_file, "wb") as output_file:
                output_file.write(response.content)
            run(f'samtools index {bam_file}')


with open('test.yaml') as f:
    data = yaml.load(f, Loader=yaml.SafeLoader)

file_location = data['file_location']
cohort = data['cohort']
token_file = data['token_file']
results_files = data['results_files']
igv_session_location = data['igv_session_location']
igv_file_location = data['igv_file_location']
cluster_location = data['cluster_location']
bam_location = data['bam_location']

run(f'aws s3 cp {results_files}/{cohort}/compare_junctions/hist/junction_pvalues_significant_0.05_filtered_BH_default.tsv .')
run(f'aws s3 cp {results_files}/{cohort}/compare_junctions/hist/junction_pvalues_significant_0.05_filtered_BH_i50e5.tsv .')
run(f'aws s3 cp {results_files}/{cohort}/compare_junctions/hist/junction_pvalues_significant_0.05_filtered_BH_E.tsv .')
run(f'aws s3 cp {results_files}/{cohort}/compare_junctions/hist/junction_pvalues_significant_0.05_filtered_BH_I.tsv .')
run(f'aws s3 cp {token_file} .')



files = glob.glob('junction_pvalues_significant_0.05_filtered_BH*.tsv')
if os.path.exists('igv_session'):
    shutil.rmtree('igv_session')

for i in files:
    tag = i.split('_')[-1].split('.')[0]
    with open(i, 'r') as result_file:
        reader = csv.DictReader(result_file, delimiter='\t')
        for line in reader:
            os.mkdir('igv_session')
            os.chdir('igv_session')
            samples_field = line['variant_junction_samples']
            samples = samples_field.split(',')
            vcf_bed_sample = samples[0]
            chrom = line['chrom']
            start = line['start']
            end = line['end']
            window_start = int(start) - 500
            window_end = int(end) + 500
            bam_window_start = int(start) - 5000
            bam_window_end = int(end) + 5000
            variant_junction = line['variant_junction_info']
            variant_junction_wo_colon = variant_junction.replace(':', '_')
            run(f'aws s3 cp {results_files}/{cohort}/{vcf_bed_sample}.tar.gz .')
            run(f'tar xzf {vcf_bed_sample}.tar.gz')
            vcf = f'{vcf_bed_sample}/output/cse_identify_filtered_compare_{tag}.vcf'
            run(f'aws s3 cp {vcf} {igv_file_location}/{vcf_bed_sample}/')
            bed = f'{vcf_bed_sample}/output/cse_identify_filtered_compare_{tag}.bed'
            run(f'aws s3 cp {bed} {igv_file_location}/{vcf_bed_sample}/')
            for sample in samples:
                get_bam(sample, token_file, chrom, bam_window_start, bam_window_end, variant_junction_wo_colon)
                bam_file = glob.glob('*.bam')
                bai_file = glob.glob('*.bai')
                # run(f'ssh kcotto@virtual-workstation4.gsc.wustl.edu "mkdir -p {cluster_location}/{sample}"')
                run(f'aws s3 cp {bam_file[0]} {bam_location}/{sample}/{bam_file[0]}')
                run(f'aws s3 cp {bai_file[0]} {bam_location}/{sample}/{bai_file[0]}')
            samples_string = '_'.join(samples)
            xml_file = f'{samples_string}_{variant_junction_wo_colon}.xml'
            with open(xml_file, 'w') as outfile:
                outfile.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n')
                outfile.write(f'<Session genome="hg38" hasGeneTrack="true" hasSequenceTrack="true" locus="{chrom}:{window_start}-{window_end}" nextAutoscaleGroup="2" version="8">\n')
                outfile.write('<Resources>\n')
                for sample in samples:
                    outfile.write(f'<Resource path="https://igv.gsc.wustl.edu{cluster_location}/{sample}/{bam_file[0]}"/>\n')
                outfile.write(f'<Resource path="https://regtools-igv-files.s3.amazonaws.com/{vcf_bed_sample}/cse_identify_filtered_compare_{tag}.vcf"/>\n')
                outfile.write(f'<Resource path="https://regtools-igv-files.s3.amazonaws.com/{vcf_bed_sample}/cse_identify_filtered_compare_{tag}.bed"/>\n')
                outfile.write('</Resources>\n')
                outfile.write('<Regions>\n')
                outfile.write(f'<Region chromosome="{chrom}" description="" end="{end}" start="{start}"/>\n')
                outfile.write('</Regions>\n')
                outfile.write('<HiddenAttributes>\n')
                outfile.write('<Attribute name="DATA FILE"/>\n')
                outfile.write('<Attribute name="DATA TYPE"/>\n')
                outfile.write('<Attribute name="NAME"/>\n')
                outfile.write('</HiddenAttributes>\n')
                outfile.write('</Session>\n')
            run(f'aws s3 cp {xml_file} {igv_session_location}/{cohort}/{tag}/')
            os.chdir('..')
            if os.path.exists('igv_session'):
                shutil.rmtree('igv_session')
