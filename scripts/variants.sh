#!/bin/bash
# usage: 
# 	bash variants.sh <input.tsv> <output.bed>
awk 'NR>1 {print $17}' $1  > $2
sed -ie 's/[,]/\n/g' $2
sed -ie 's/[:\-]/\t/g' $2
sort $2 | uniq > ${2}.tmp 
cat ${2}.tmp > $2
rm -f ${2}.tmp
rm ${2}e
