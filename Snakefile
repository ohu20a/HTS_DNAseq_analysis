localrules: all, link_ref, collect_stats, unzip_ref, raw_multiqc, trimmed_multiqc
import re
import os
from pprint import pprint
from collections import defaultdict
# File name should follow the format: sample_lane_read(R1 or R2).fastq.gz
# For example: NAD215_1_R1.fastq.gz
# Lane number should always exist even if a sample is only sequenced on one lane.

configfile: "config.yaml"

# Get file names
RAW_FILES, EXT,  = glob_wildcards("raw_data/{raw_file,[\w\d]+_[\w\d]+_R(1|2)}.{ext,fastq(\.gz)?}")
SAMPLES = []
SAMPLES_DICT = defaultdict(lambda: defaultdict(lambda: defaultdict(str)))
GENOME = os.path.basename(config['genome'])

for file in RAW_FILES:
    match = re.match(r"([\w\d]+)_([\w\d]+)_(R[12])", file)
    if match:
        sample, lane, read = match[1], match[2], match[3]
        SAMPLES.append(sample)
        SAMPLES_DICT[sample][lane][read] = file
# given a sample name, return all lanes of that sample
def sample_list(sample):
    list = []
    for lane in SAMPLES_DICT[sample]:
        file = SAMPLES_DICT[sample][lane]['R1']
        match = re.match(r"([\w\d]+?_[\w\d]+?)_R[12]", file)
        list.append(match[1])
#    print(list)
    return list

rule all:
    input: 
        "Results/raw_fastqc/multiqc_report.html",
        "Results/trimmed_fastqc/multiqc_report.html",
#        expand("Results/mapping/{sample}.merged.sorted.bam.bai", sample = SAMPLES)
        "Results/metrics/mapping_stats.csv",
        "Results/variants/variants.vcf"

rule link_ref:
    input: expand("{genome}.{ext}", genome = config['genome'], ext = config['genome_extension'])
    output: expand("ref/{genome}.{ext}", genome = GENOME, ext = config['genome_extension'])
    shell:
     """
     ln -s {input} {output}
     """

rule unzip_ref:
    input: expand("{{file}}.{ext}", ext = config['genome_extension'])
    output: "{file}.fa"
    shell:
     """
     gunzip -k -c {input} > {output}
     """

rule fa_idx:
    input: "ref/{file}.fa"
    output: "ref/{file}.fa.fai"
    params: time = "60"
    shell:
     """
     samtools faidx {input}
     """

rule raw_fastqc:
    input: "raw_data/{sample}.fastq.gz"
    output: "Results/raw_fastqc/{sample}.html"
    params: time="60"
    conda: "envs/fastqc.yaml"
    shell:
     """
     fastqc -o Results/raw_fastqc/ {input}
     """

rule raw_multiqc:
    input: expand("Results/raw_fastqc/{sample}.html", sample = RAW_FILES)
    output: "Results/raw_fastqc/multiqc_report.html"
    conda: "envs/fastqc.yaml"
    shell:
     """
     multiqc -o Results/raw_fastqc/ Results/raw_fastqc/
     """

rule trim:
    input: 
        r1 = "raw_data/{sample}_R1.fastq.gz",
        r2 = "raw_data/{sample}_R2.fastq.gz"
    output: 
        r1 = "Results/trimmed_data/{sample}_R1.qc.fq.gz",
        r2 = "Results/trimmed_data/{sample}_R2.qc.fq.gz"
    conda: "envs/fastp.yaml"
    params: time = "240"
    threads: 4
    shell:
     """
     fastp -i {input.r1} -o {output.r1} -I {input.r2} -O {output.r2} --detect_adapter_for_pe --thread {threads}
     """

rule trim_fastqc:
    input: "Results/trimmed_data/{sample}.qc.fq.gz"
    output: "Results/trimmed_fastqc/{sample}.html"
    conda: "envs/fastqc.yaml"
    params: time = "60"
    shell:
     """
     fastqc -o Results/trimmed_fastqc {input}
     """

rule trimmed_multiqc:
    input: expand("Results/trimmed_fastqc/{sample}.html", sample = RAW_FILES)
    output: "Results/trimmed_fastqc/multiqc_report.html"
    conda: "envs/fastqc.yaml"
    shell:
     """
     multiqc -o Results/trimmed_fastqc/ Results/trimmed_fastqc/
     """

rule bwa_index:
    input: expand("ref/{genome}.fa", genome = GENOME)
    output: expand("ref/{genome}.{ext}", genome = GENOME, ext = ["amb", "ann", "bwt", "pac", "sa"])
    params: time = "120"
    conda: "envs/bwa.yaml"
    shell:
     """
     bwa index {input}
     """

rule bwa:
    input:
        rules.bwa_index.output, 
        r1 = "Results/trimmed_data/{sample}_{lane}_R1.qc.fq.gz",
        r2 = "Results/trimmed_data/{sample}_{lane}_R2.qc.fq.gz"
    output: temp("Results/mapping/{sample}_{lane}.aligned.bam")
    conda: "envs/bwa.yaml"
    params: time="2-0"
    threads: 4
    log: "logs/mapping/bwa_{sample}_{lane}.log"
    benchmark: "benchmarks/mapping/{sample}_{lane}.tsv"
    shell:
     """
     bwa mem -t {threads} -R "@RG\\tID:{wildcards.sample}_{wildcards.lane}\\tSM:{wildcards.sample}\\tPL:illumina" -o {output} {config[genome]} | samtools -bS - > {output}
     """

rule combine_bam:
    input: lambda wildcards: expand("Results/mapping/{file}.aligned.sorted.bam", file = sample_list(wildcards.sample))
    output: temp("Results/mapping/{sample}.merged.bam")
    threads: 2
    params: time = "300"
    conda: "envs/sambamba.yaml"
    script:
     """
     tools/merge_bam.py
     """

rule sort_bam:
    input: "Results/mapping/{sample}.{stage}.bam"
    output: temp("Results/mapping/{sample}.{stage}.sorted.bam")
    threads: 4
    params: time = "240"
    conda: "envs/samtools.yaml"
    shell:
     """
     samtools sort -o {output} -@ {threads} -m 1G {input}
     """

rule index_bam:
    input: "Results/mapping/{sample}.bam"
    output: "Results/mapping/{sample}.bam.bai"
    conda: "envs/samtools.yaml"
    params: time = "240"
    shell:
     """
     samtools index {input}
     """

rule bam_stat:
    input: 
        bai = "Results/mapping/{sample}.merged.sorted.bam.bai",
        bam = "Results/mapping/{sample}.merged.sorted.bam"
    output: "Results/mapping/stats/{sample}.stats.txt"
    conda: "envs/sambamba.yaml"
    params: time = "60"
    threads: 2
    log: "Results/logs/mapping/stat_bam_{sample}.log"
    benchmark: "Results/benchmarks/mapping/stat_bam_{sample}.tsv"
    shell:
     """
     sambamba flagstat -t {threads} {input.bam} 2>{log} 1>{output}
     """

rule collec_stats:
    input: expand("Results/mapping/stats/{samples}.stats.txt", samples = SAMPLES)
    output: "Results/metrics/mapping_stats.csv"
    params: dir = "Results/mapping/stats/"
    conda: "envs/stat_curator.yaml"
    script: "tools/flagstat_curator.py"

rule remove_duplicates:
    input: "Results/mapping/{sample}.merged.sorted.bam"
    output: "Results/mapping/{sample}.dedup.bam"
    conda: "envs/sambamba.yaml"
    threads: 4
    params: time = "240"
    log: "Results/logs/mapping/dedup_{sample}.log"
    benchmark: "Results/benchmarks/mapping/dedup_{sample}.tsv"
    shell:
     """
     sambamba markdup -r -t {threads} {input} {output} &> {log}    
     """

rule merge_bam:
# merge all bams. This helps reduce RAM footprint for freebayes
    input: expand("Results/mapping/{sample}.dedup.bam", sample = SAMPLES)
    output: "Results/mapping/merged.bam"
    conda: "envs/sambamba.yaml"
    threads: 5
    params: time = "300"
    log: "Results/logs/mapping/merge_bam.log"
    benchmark: "Results/benchmarks/mapping/merge_bam.tsv"
    shell:
     """
     sambamba merge -t {threads} -l 9 {output} {input} 
     """

rule freebayes:
    input:
        bam = expand("Results/mapping/{sample}.dedup.bam", sample = SAMPLES),
        bai = expand("Results/mapping/{sample}.dedup.bam.bai", sample = SAMPLES),
        fai = expand("ref/{genome}.fa.fai", genome = GENOME),
        fa = expand("ref/{genome}.fa", genome = GENOME)
    output: "Results/variants/variants.vcf"
    conda: "envs/freebayes.yaml"
    params: time="3-0"
    threads: 36
    log: "Results/logs/variant_calling/freebayes.log"
    benchmark: "Results/benchmarks/variant_calling/freebayes.tsv"
    shell:
     """
     freebayes-parallel <(fasta_generate_regions.py {input.fai} 100000) {threads} -f {input.fa} {input.bam} > {output}
     """

