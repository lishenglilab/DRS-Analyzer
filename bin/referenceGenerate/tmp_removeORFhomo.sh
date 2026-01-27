#!/bin/bash
####bin/referenceGenerate/removeORFhomo.sh
for i in $(seq 1 10)
do
{
    blastp -query '/home/yvzeng/Project/mouseDRS/data/nextflow/testdata/result/work/86/26ef6d1e49e6984e37145d6ceed17a/split/ORF/ORF_no_redundancy_'${i}'.fa' \
           -out '/home/yvzeng/Project/mouseDRS/data/nextflow/testdata/result/work/86/26ef6d1e49e6984e37145d6ceed17a/split/blast/ORF_no_redundancy_'${i}'.blast' \
           -db /home/yvzeng/Project/mouseDRS/data/nextflow/testdata/result/work/86/26ef6d1e49e6984e37145d6ceed17a/db/ref.db \
           -outfmt 6 -evalue 1e-5 -num_threads 1
}&
done
wait