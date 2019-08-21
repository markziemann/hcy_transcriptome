#!/bin/bash

for FQZ1 in *R1_001.fastq.gz ; do
  FQZ2=$(echo $FQZ1 | sed 's/_R1_/_R2_/' )
  skewer -q 10 -t 8 $FQZ1 $FQZ2
done
