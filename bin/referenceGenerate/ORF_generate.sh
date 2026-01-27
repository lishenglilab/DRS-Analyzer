#!/bin/bash
####bin/referenceGenerate/ORF_generate.sh
for i in `seq 1 number`
do
{
python 1.ORF_generate.py -i transcript_${i}.fa -o path -p ORF_${i}
}&
done
wait