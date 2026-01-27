#!/bin/bash
####bin/referenceGenerate/removeORFhomo.sh
for i in $(seq 1 {n})
do
{{
    blastp -query '{orf}_'${{i}}'.fa' \
           -out '{path}_'${{i}}'.blast' \
           -db {db} \
           -outfmt 6 -evalue 1e-5 -num_threads 1
}}&
done
wait