

from pathlib import Path
from collections import OrderedDict, defaultdict
import json


# configfile: "./config_seCLIP_YJL.yaml"
assert len(config) > 0, "--configfile is required"

outdir = Path(config["outdir"])

if "bin" in config:
    umi_tools = config["bin"].get("umi_tools", "umi_tools")

SAMPLE_DATA = OrderedDict()
path = config["data"]["path"]
PAIRED_SAMPLES = dict()
for sn, data in config["data"]["samples"].items():
    if sn.startswith("demo"):
        continue
    data["R1"] = os.path.join(path, data["R1"])
    data["R2"] = os.path.join(path, data["R2"])
    SAMPLE_DATA[sn] = data
    cell = data["cell"]
    rbp = data["rbp"]
    rep = data["rep"]
    treatment = data["treatment"]
    group = data["group"]
    assert group in ["IP", "input"]
    key = '_'.join([cell, rbp, treatment])
    if key not in PAIRED_SAMPLES:
        PAIRED_SAMPLES[key] = {}
    if rep not in PAIRED_SAMPLES[key]:
        PAIRED_SAMPLES[key][rep] = {}
    assert group not in PAIRED_SAMPLES[key][rep], "{}".format((PAIRED_SAMPLES[key][rep], sn))
    PAIRED_SAMPLES[key][rep][group] = sn
IP_SAMPLES = [sn for sn, data in SAMPLE_DATA.items() if data["group"] == "IP"]

print(SAMPLE_DATA, flush=True)
print(json.dumps(PAIRED_SAMPLES, indent=4), flush=True)

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
TOOL_DIR = os.path.join(FILE_DIR, "tools")
CLIPPER_DATA = '/scratch/midway3/kenchen/Documents/seCLIP-snakemake/reference/clipper_data'
BIN = {
    "clipper_img": config.get("bin", {}).get("clipper_img", os.path.join(TOOL_DIR, "clipper_5d865bb.sif")),
}

config["clipper_gff"] = os.path.expanduser(config["clipper_gff"])


rule all:
    input:
        expand(outdir / "alignment_star/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam", sample=SAMPLE_DATA.keys()),
        #bed = outdir / "peak_cluster/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDupSo.peakClusters.bed"
        [
            outdir / f"peak_cluster/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDupSo.peakClusters.bed"
            for sample in IP_SAMPLES
        ]

rule link_fastq:
    input:
        r1 = lambda wildcards: SAMPLE_DATA[wildcards.sample]["R1"],
    output:
        r1 = outdir / "reads_raw/{sample}/{sample}.r1.fq.gz"
    shell:
        "ln -s {input.r1} {output.r1}"

rule step84_umitools_extract:
    input:
        r1 = outdir / "reads_raw/{sample}/{sample}.r1.fq.gz"
    output:
        # r1 = "rep1.IP.umi.r1.fq.gz"
        r1 = outdir / "reads_clean/{sample}/{sample}.umi.r1.fq.gz"
    log:
        outdir / "reads_clean/{sample}/{sample}.umitools.metrics"
    params:
        bc_pattern = "NNNNNNNNNN"
    threads:
        1
    resources:
        mem_mb = 8 * 1024
    shell:
        """
        umi_tools extract \
            --random-seed 1 \
            --bc-pattern {params.bc_pattern} \
            --stdin {input.r1} \
            --stdout {output.r1} \
            --log {log}
        """

rule step85_cutadapt_trim:
    input:
        r1 = outdir / "reads_clean/{sample}/{sample}.umi.r1.fq.gz"
    output:
        r1 = outdir / "reads_clean/{sample}/{sample}.umi.r1.fqTr.fq.gz",
        metrics = outdir / "reads_clean/{sample}/{sample}.umi.r1.fqTr.metrics"
    resources:
        mem_mb = 8 * 1024
    threads:
        8
    shell:
        """
        cutadapt -O 1 \
            --match-read-wildcards \
            --times 1 \
            -e 0.1 \
            --quality-cutoff 6 \
            -m 18 \
            -a AGATCGGAAGAGCAC \
            -a GATCGGAAGAGCACA \
            -a ATCGGAAGAGCACAC \
            -a TCGGAAGAGCACACG \
            -a CGGAAGAGCACACGT \
            -a GGAAGAGCACACGTC \
            -a GAAGAGCACACGTCT \
            -a AAGAGCACACGTCTG \
            -a AGAGCACACGTCTGA \
            -a GAGCACACGTCTGAA \
            -a AGCACACGTCTGAAC \
            -a GCACACGTCTGAACT \
            -a CACACGTCTGAACTC \
            -a ACACGTCTGAACTCC \
            -a CACGTCTGAACTCCA \
            -a ACGTCTGAACTCCAG \
            -a CGTCTGAACTCCAGT \
            -a GTCTGAACTCCAGTC \
            -a TCTGAACTCCAGTCA \
            -a CTGAACTCCAGTCAC \
            -o {output.r1} \
            -j {threads} \
            {input.r1} > {output.metrics}
        """



rule step86_cutadapt_trim_again:
    input:
        r1 = outdir / "reads_clean/{sample}/{sample}.umi.r1.fqTr.fq.gz"
    output:
        r1 = outdir / "reads_clean/{sample}/{sample}.umi.r1.fqTrTr.fq",
        metrics = outdir / "reads_clean/{sample}/{sample}.umi.r1.fqTrTr.metrics"
    resources:
        mem_mb = 8 * 1024
    threads:
        8
    shell:
        """
        cutadapt \
            -O 5 \
            --match-read-wildcards \
            --times 1 \
            -e 0.1 \
            --quality-cutoff 6 \
            -m 18 \
            -a AGATCGGAAGAGCAC \
            -a GATCGGAAGAGCACA \
            -a ATCGGAAGAGCACAC \
            -a TCGGAAGAGCACACG \
            -a CGGAAGAGCACACGT \
            -a GGAAGAGCACACGTC \
            -a GAAGAGCACACGTCT \
            -a AAGAGCACACGTCTG \
            -a AGAGCACACGTCTGA \
            -a GAGCACACGTCTGAA \
            -a AGCACACGTCTGAAC \
            -a GCACACGTCTGAACT \
            -a CACACGTCTGAACTC \
            -a ACACGTCTGAACTCC \
            -a CACGTCTGAACTCCA \
            -a ACGTCTGAACTCCAG \
            -a CGTCTGAACTCCAGT \
            -a GTCTGAACTCCAGTC \
            -a TCTGAACTCCAGTCA \
            -a CTGAACTCCAGTCAC \
            -o {output.r1} \
            -j {threads} \
            {input.r1} > {output.metrics}
        """

rule step87_fastq_sort:
    input:
        r1 = outdir / "reads_clean/{sample}/{sample}.umi.r1.fqTrTr.fq"
    output:
        r1 = outdir / "reads_clean/{sample}/{sample}.umi.r1.fqTrTr.sorted.fq.gz"
    resources:
        mem_mb = 4 * 1024
    shell:
        """
        fastq-sort --id {input.r1} | gzip > {output.r1}
        """

rule step88_star_to_repeat:
    input:
        r1 = outdir / "reads_clean/{sample}/{sample}.umi.r1.fqTrTr.sorted.fq.gz"
    output:
        bam = outdir / "alignment_star/{sample}/{sample}.umi.r1.fqTrTr.sorted.STAR.Aligned.out.bam",
        unmap = outdir / "alignment_star/{sample}/{sample}.umi.r1.fqTrTr.sorted.STAR.Unmapped.out.mate1"
    params:
        prefix = lambda wildcards: outdir / f"alignment_star/{wildcards.sample}/{wildcards.sample}.umi.r1.fqTrTr.sorted.STAR.",
        tmpdir = lambda wildcards: outdir / f"alignment_star/{wildcards.sample}/_tmpdir",
        index = "/project/mengjiechen/kenchen/db/humrep_repbase"
    resources:
        mem_mb = 36 * 1024
    threads:
        16
    shell:
        """
        test -d {params.tmpdir} && rm -rf {params.tmpdir}
        STAR \
            --alignEndsType EndToEnd \
            --genomeDir {params.index} \
            --genomeLoad NoSharedMemory \
            --outBAMcompression 10 \
            --outFilterMultimapNmax 30 \
            --outFilterMultimapScoreRange 1 \
            --outFilterScoreMin 10 \
            --outFilterType BySJout \
            --outReadsUnmapped Fastx \
            --outSAMattrRGline ID:{wildcards.sample} \
            --outSAMattributes All \
            --outSAMmode Full \
            --outSAMtype BAM Unsorted \
            --outSAMunmapped Within \
            --outStd Log \
            --runMode alignReads \
            --readFilesCommand zcat \
            --runThreadN {threads} \
            --readFilesIn {input.r1} \
            --outTmpDir {params.tmpdir} \
            --outFileNamePrefix {params.prefix} &> {output.bam}.log
        """

rule step89_star_to_genome:
    input:
        r1 = outdir / "alignment_star/{sample}/{sample}.umi.r1.fqTrTr.sorted.STAR.Unmapped.out.mate1"
    output:
        bam = outdir / "alignment_star/{sample}/{sample}.umi.r1.fq.genome-mapped.Aligned.out.bam"
    params:
        prefix = lambda wildcards: outdir / f"alignment_star/{wildcards.sample}/{wildcards.sample}.umi.r1.fq.genome-mapped.",
        tmpdir = lambda wildcards: outdir / f"alignment_star/{wildcards.sample}/_tmpdir",
        index = "/project/mengjiechen/kenchen/db/GRCh38-gencode.STAR_overhang100"
    resources:
        mem_mb = 36 * 1024
    threads:
        16
    shell:
        """
        test -d {params.tmpdir} && rm -rf {params.tmpdir}
        STAR \
            --alignEndsType EndToEnd \
            --genomeDir {params.index} \
            --genomeLoad NoSharedMemory \
            --outBAMcompression 10 \
            --outFilterMultimapNmax 1 \
            --outFilterMultimapScoreRange 1 \
            --outFilterScoreMin 10 \
            --outFilterType BySJout \
            --outReadsUnmapped Fastx \
            --outSAMattrRGline ID:{wildcards.sample}\
            --outSAMattributes All \
            --outSAMmode Full \
            --outSAMtype BAM Unsorted \
            --outSAMunmapped Within \
            --outStd Log \
            --runMode alignReads \
            --runThreadN {threads} \
            --outFileNamePrefix {params.prefix} \
            --outTmpDir {params.tmpdir} \
            --readFilesIn {input.r1} &> {output.bam}.log
        """

rule step90_sort_uniquely_mapped_reads:
    input:
        bam = outdir / "alignment_star/{sample}/{sample}.umi.r1.fq.genome-mapped.Aligned.out.bam"
    output:
        bam = temp(outdir / "alignment_star/{sample}/{sample}.umi.r1.fq.genome-mappedSo.bam"),
        bam_sorted = outdir / "alignment_star/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.bam",
        bam_sorted_bai = outdir / "alignment_star/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.bam.bai"
    resources:
        mem_mb = 1024 * 24
    threads:
        8
    shell:
        """
        samtools sort -@ {threads} -n {input.bam} > {output.bam} && \
        samtools sort -@ {threads} {output.bam} > {output.bam_sorted} && \
        samtools index -@ {threads} {output.bam_sorted}
        """

rule step91_umitools_dedup:
    input:
        bam = outdir / "alignment_star/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.bam"
    output:
        bam = outdir / "alignment_star/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDup.bam",
        bam_sorted = outdir / "alignment_star/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam",
        bam_sorted_bai = outdir / "alignment_star/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam.bai",
    params:
        prefix = lambda wildcards: outdir / f"alignment_star/{wildcards.sample}/{wildcards.sample}.umi.r1.fq.genome-mappedSoSo"
    resources:
        mem_mb = 1024 * 24
    threads:
        1
    shell:
        """
        umi_tools dedup \
            --random-seed 1 \
            --method unique \
            --output-stats {params.prefix} \
            -I {input.bam} \
            -S {output.bam} && \
        samtools sort \
            -@ {threads} \
            -o {output.bam_sorted} \
            {output.bam} && \
        samtools index -@ {threads} {output.bam_sorted}
        """



# 

# CHROMOSOMES = ["chr{}".format(c) for c in range(1, 23)] + ["chrX", "chrY", "chrM"]
CHROMOSOMES = ["chr1"]
CLIPPER_GFF_SUFFIX = ".AS.STRUCTURE.COMPILED.gff"
#def split_gtf(annotation: str, chromosomes: list):
#    # assert "GRCh38_v29e.AS.STRUCTURE.COMPILED.gff"
#    assert annotation.endswith(CLIPPER_GFF_SUFFIX), "Annotation file must be named as '[species].AS.STRUCTURE.COMPILED.gff'"
#    prefix = annotation.replace(CLIPPER_GFF_SUFFIX, "")
#    outputs = {c: open(prefix + "_" + c + CLIPPER_GFF_SUFFIX, 'wt') for c in chromosomes}
#    found = set()
#    with open(annotation) as f:
#        for line in f:
#            if line.startswith("#"):
#                for c in chromosomes:
#                    print(line, file=outputs[c], end="")
#            else:
#                chrom = line.split("\t")[0]
#                found.add(chrom)
#                if chrom in chromosomes:
#                    print(line, file=outputs[chrom], end="")
#    for c in chromosomes:
#        outputs[c].close()

def split_gtf(annotation: str, N=100):
    # assert "GRCh38_v29e.AS.STRUCTURE.COMPILED.gff"
    assert annotation.endswith(CLIPPER_GFF_SUFFIX), "Annotation file must be named as '[species].AS.STRUCTURE.COMPILED.gff'"
    prefix = annotation.replace(CLIPPER_GFF_SUFFIX, "")
    outputs = {i: open(prefix + "_line-r{}".format(i) + CLIPPER_GFF_SUFFIX, 'wt') for i in range(N)}
    found = set()
    with open(annotation) as f:
        for nr, line in enumerate(f):
            if line.startswith("#"):
                for i in range(N):
                    print(line, file=outputs[i], end="")
            else:
                r = nr % N
                print(line, file=outputs[r], end="")
    for c in range(N):
        outputs[c].close()


rule custom_step_split_annotation:
    input:
        gff = config["clipper_gff"],
    output:
        expand(
            config["clipper_gff"].replace(CLIPPER_GFF_SUFFIX, "_line-r{nr}" + CLIPPER_GFF_SUFFIX),
            nr=range(100)
        )
    threads:
        1
    resources:
            mem_mb = 1024
    run:
        split_gtf(input.gff, 100)

rule step92_clipper_peakcluster:
    input:
        bam = outdir / "alignment_star/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam",
        gff = config["clipper_gff"].replace(CLIPPER_GFF_SUFFIX, "_line-r{nr}" + CLIPPER_GFF_SUFFIX)
    params:
        species = os.path.basename(config["clipper_gff"]).replace(CLIPPER_GFF_SUFFIX, ""),
        clipper_img = BIN["clipper_img"],
        ref_path = CLIPPER_DATA
    output:
        chrom_bed = temp(outdir / "peak_cluster/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDupSo.peakClusters.line-r{nr}.bed")
    log:
        temp(outdir / "peak_cluster/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDupSo.peakClusters.line-r{nr}.log")
    threads:
        28
    resources:
        mem_mb = 1024 * 48
    shell:
        """
        singularity run \
            --bind {params.ref_path}:/opt/conda/envs/clipper3/lib/python3.7/site-packages/clipper-2.0.0-py3.7-linux-x86_64.egg/clipper/data \
            {params.clipper_img} clipper \
            --species {params.species} \
            -g ENSG00000227232 \
            --bam {input.bam} \
            --outfile {output.chrom_bed} \
            --processors {threads} &> {log}
        """
rule step92_merge_peak_clusters:
    input:
        bed = [
            outdir / f"peak_cluster/{{sample}}/{{sample}}.umi.r1.fq.genome-mappedSoSo.rmDupSo.peakClusters.line-r{nr}.bed"
            for nr in range(100)
        ],
        log = [
            outdir / f"peak_cluster/{{sample}}/{{sample}}.umi.r1.fq.genome-mappedSoSo.rmDupSo.peakClusters.line-r{nr}.log"
            for nr in range(100)
        ]
    output:
        bed = outdir / "peak_cluster/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDupSo.peakClusters.bed"
    log:
        outdir / "peak_cluster/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDupSo.peakClusters.log"
    resources:
        mem_mb = 1024 // 2
    threads:
        1
    shell:
        """
        cat {input.bed} > {output} && \
        cat {input.log} > {log}
        """
