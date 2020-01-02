import csv

default = 'junction_pvalues_significant_0.05_filtered_BH_default.tsv'
i50e5 = 'junction_pvalues_significant_0.05_filtered_BH_i50e5.tsv'
E = 'junction_pvalues_significant_0.05_filtered_BH_E.tsv'
I = 'junction_pvalues_significant_0.05_filtered_BH_I.tsv'

default_set = set()
i50e5_set = set()

new_default_file = 'mutually_exclusive_junction_pvalues_significant_0.05_filtered_BH_default.tsv'
new_i50e5_file = 'mutually_exclusive_junction_pvalues_significant_0.05_filtered_BH_i50e5.tsv'
new_E_file = 'mutually_exclusive_junction_pvalues_significant_0.05_filtered_BH_E.tsv'
new_I_file = 'mutually_exclusive_junction_pvalues_significant_0.05_filtered_BH_I.tsv'

with open(default, 'r') as d, open(new_default_file, 'w') as outfile:
    line_count = 0
    reader = csv.reader(d, delimiter='\t')
    for line in reader:
        line_count += 1
        if line_count == 1:
            outfile.write('\t'.join(line) + '\n')
        if line_count != 1:
            variant_junction = line[7]
            default_set.add(variant_junction)
            outfile.write('\t'.join(line) + '\n')

with open(i50e5, 'r') as i50e5_file, open(new_i50e5_file, 'w') as outfile:
    line_count = 0
    reader = csv.reader(i50e5_file, delimiter='\t')
    for line in reader:
        line_count += 1
        if line_count == 1:
            outfile.write('\t'.join(line) + '\n')
        if line_count != 1:
            variant_junction = line[7]
            if variant_junction not in default_set:
                outfile.write('\t'.join(line) + '\n')
                i50e5_set.add(variant_junction)
            else:
                print(f'This record ({variant_junction}) exists in the default results file.')

combined_set = default_set.union(i50e5_set)

with open(E, 'r') as E_file, open(new_E_file, 'w') as outfile:
    reader = csv.reader(E_file, delimiter='\t')
    for line in reader:
        line_count += 1
        if line_count == 1:
            outfile.write('\t'.join(line) + '\n')
        if line_count != 1:
            variant_junction = line[7]
            if variant_junction not in combined_set:
                outfile.write('\t'.join(line) + '\n')
            else:
                print(f'This record ({variant_junction}) exists in either the default or i50e5 results file.')

with open(I, 'r') as I_file, open(new_I_file, 'w') as outfile:
    reader = csv.reader(I_file, delimiter='\t')
    for line in reader:
        line_count += 1
        if line_count == 1:
            outfile.write('\t'.join(line) + '\n')
        if line_count != 1:
            variant_junction = line[7]
            if variant_junction not in combined_set:
                outfile.write('\t'.join(line) + '\n')
            else:
                print(f'This record ({variant_junction}) exists in either the default or i50e5 results file.')


print(len(default_set))
print(len(i50e5_set))
print(len(combined_set))