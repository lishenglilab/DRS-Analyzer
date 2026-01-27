#!/bin/bash
####bin/referenceGenerate/ORF_generate.sh
for i in `seq 1 40`
do
{
python /home/yvzeng/Project/mouseDRS/data/nextflow/DRS-Analyzer/bin/referenceGenerate/1.ORF_generate.py -i /home/yvzeng/Project/mouseDRS/data/nextflow/testdata/result/work/db/b5b559c44e30fe47b2f15d6945ca59/split/trancript/filtered_isoforms_${i}.fa -o /home/yvzeng/Project/mouseDRS/data/nextflow/testdata/result/work/db/b5b559c44e30fe47b2f15d6945ca59/split/ORF/ -p ORF_${i}
}&
done
wait