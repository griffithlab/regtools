import csv
import fnmatch

regtools_file = '/Users/kcotto/Downloads/junction_pvalues_significant_0.05_filtered_BH_default.tsv'
spliceAI_file = '/Users/kcotto/Downloads/all_variants_sorted_spliceAI.vcf'
gtex_file = '/Users/kcotto/Projects/regtools/TCGA_paper/GTEx_matrices_avg/Liver.tsv'
base_file = regtools_file.split('.')[0]
tmp_file = f'{base_file}_tmp.tsv'
final_output_file = f'{base_file}_gtex_spliceai.tsv'


def annotate_GTEx(regtools_input_file, gtex_file):
    with open(regtools_input_file, 'r') as regtools, open(gtex_file, 'r') as gtex:
        gtex_values = {}
        gtex_reader = csv.reader(gtex, delimiter='\t')
        regtools_reader = csv.reader(regtools, delimiter='\t')
        next(gtex_reader)
        header = next(regtools_reader)
        outputfile = tmp_file
        match_count = 0
        with open(outputfile, 'w') as outfile:
            reformatted_header = '\t'.join(header)
            outfile.write(f'{reformatted_header}\tGTEx_mean\tGTEx_sd\n')
            for line in gtex_reader:
                mean = line[2]
                stdev = line[3]
                gtex_key = line[0]
                gtex_values[gtex_key] = (mean, stdev)
            for line in regtools_reader:
                chrom = line[1]
                junction_start = str(int(line[2]) + 1)
                junction_end = str(int(line[3]) -1)
                regtools_key = f'{chrom}_{junction_start}_{junction_end}'
                output_line = '\t'.join(line)
                if regtools_key in gtex_values:
                    match_count += 1
                    print(f'match for {output_line}')
                    outfile.write(f'{output_line}\t{gtex_values[regtools_key][0]}\t{gtex_values[regtools_key][1]}\n')
                else:
                    outfile.write(f'{output_line}\tNA\tNA\n')
    print(match_count)

def annotate_spliceAI(regtools_input_file, spliceAI_file):
    with open(regtools_input_file, 'r') as regtools, open(spliceAI_file, 'r') as spliceAI:
        spliceAI_reader = csv.reader(spliceAI, delimiter='\t')
        regtools_reader = csv.reader(regtools, delimiter='\t')
        spliceAI_values = {}
        for line in spliceAI_reader:
            li = str(line[0]).strip()
            if not li.startswith('#'):
                info_field = line[7]
                if str(info_field).__contains__('SpliceAI'):
                        pattern = 'SpliceAI*'
                        extract_SF = line[7].split(';')
                        matching = fnmatch.filter(extract_SF, pattern)
                        spliceAI_variant_key = f'{line[0]}:{line[1]}'
                        spliceAI_string = matching[0]
                        if spliceAI_string.__contains__(','):
                            split_values = spliceAI_string.split(',')
                            spliceAI_values[spliceAI_variant_key] = split_values[0]
                        else:
                            spliceAI_values[spliceAI_variant_key] = spliceAI_string
        header = next(regtools_reader)
        outputfile = final_output_file
        with open (outputfile, 'w') as outfile:
            reformatted_header = '\t'.join(header)
            outfile.write(f'{reformatted_header}\tSpliceAI_raw\tSpliceAI_match\n')
            for line in regtools_reader:
                chrom = line[6].split(':')[0]
                variant = line[6].split('-')[-1]
                regtools_variant_key = f'{chrom}:{variant}'
                output_line = '\t'.join(line)
                if regtools_variant_key in spliceAI_values:
                    junction_start_match = False
                    junction_end_match = False
                    junction_start = int(line[2])
                    junction_end = int(line[3])
                    spliceAI_info = spliceAI_values[regtools_variant_key]
                    spliceAI_split = spliceAI_info.split('|')
                    gene = spliceAI_split[1]
                    positions = spliceAI_split[-4:]
                    for position in positions:
                        genomic_location = int(variant) + int(position)
                        if junction_start == genomic_location:
                            junction_start_match = True
                        if junction_end == genomic_location:
                            junction_end_match = True
                    if junction_start_match and junction_end_match:
                        outfile.write(f'{output_line}\t{spliceAI_info}\tjunction start and end match\n')
                    elif junction_start_match:
                        outfile.write(f'{output_line}\t{spliceAI_info}\tjunction start match\n')
                    elif junction_end_match:
                        outfile.write(f'{output_line}\t{spliceAI_info}\tjunction end match\n')
                    else:
                        outfile.write(f'{output_line}\t{spliceAI_info}\tNA\n')
                else:
                    outfile.write(f'{output_line}\tNA\tNA\n')


annotate_GTEx(regtools_file, gtex_file)
annotate_spliceAI(tmp_file, spliceAI_file)
