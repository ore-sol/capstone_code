#!/usr/bin/bash

mkdir fastqc

eval "$(~/miniforge3/bin/conda shell.bash hook)"

conda activate fastqc

fastqc --noextract -o fastqc *.fastq.gz

cd fastqc 

multiqc .