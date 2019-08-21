#!/bin/bash

fastqc -t 8 *fastq.gz

multiqc .

