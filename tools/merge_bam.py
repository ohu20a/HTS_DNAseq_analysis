#!/usr/bin/env/python

## Check if multiple lanes of bam files present for given sample
### If no, rename bam file to remove lane number
### If yes, merge them and output a single bam file

import subprocess

sample = snakemake.wildcards.sample
files = snakemake.input
output = snakemake.output
threads = snakemake.threads
print(output)
if len(files) == 1:
    subprocess.run(["mv", files[0], output[0]])
else:
    subprocess.run(["sambamba", "-t", threads, output[0]].extend(files))
