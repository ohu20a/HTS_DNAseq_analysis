#!/usr/bin/env/python

## Check if multiple lanes of bam files present for given sample
### If no, rename bam file to remove lane number
### If yes, merge them and output a single bam file

import subprocess

sample = snakemake.wildcards.sample
files = snakemake.input
output = snakemake.output
if len(files) == 1:
    subprocess.run(["mv", files[0], output])
else:
    subprocess.run(["samtools", output].extend(files))
