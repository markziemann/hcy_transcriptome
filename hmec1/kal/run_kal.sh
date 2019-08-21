#!/bin/bash
for FQ1 in *pair1.fastq ; do
  FQ2=$( echo $FQ1 | sed 's/pair1.fastq/pair2.fastq/' )
  ../sw/kallisto/kallisto quant \
  -i ../ref/gencode.v31.transcripts.fa.gz.idx \
  -o ${FQ1}_kal -t 16 $FQ1 $FQ2
done


for TSV in */*abundance.tsv ; do
  NAME=$(echo $TSV | cut -d '_' -f1) ; cut -f1,4 $TSV | sed 1d | sed "s/^/${NAME}\t/"
done > 3col.tsv
