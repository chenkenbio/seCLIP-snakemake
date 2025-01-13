

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
IP2IN = dict() 
for key in PAIRED_SAMPLES:
    for rep in PAIRED_SAMPLES[key]:
        ip_sn = PAIRED_SAMPLES[key][rep]["IP"]
        in_sn = PAIRED_SAMPLES[key][rep]["input"]
        IP2IN[ip_sn] = in_sn


FILE_DIR = os.path.dirname(os.path.abspath(__file__))
TOOL_DIR = os.path.join(FILE_DIR, "tools")
CLIPPER_DATA = '/scratch/midway3/kenchen/Documents/seCLIP-snakemake/reference/clipper_data'
BIN = {
    "clipper_img": config.get("bin", {}).get("clipper_img", os.path.join(TOOL_DIR, "clipper_5d865bb.sif")),
    "makebigwigfiles_img": config.get("bin", {}).get("makebigwigfiles_img", os.path.join(TOOL_DIR, "makebigwigfiles_0.0.3.sif")),
    "repeatitive_element_mapping_img": config.get("bin", {}).get("repeatitive_element_mapping_img", os.path.join(TOOL_DIR, "repeatitive_element_mapping_1.0.0.sif")),
    "eclip_perl_img": "/home/kenchen/Documents/seCLIP-snakemake/tools/eclip_0.7.0_perl.sif"
}
for k, v in BIN.items():
    assert os.path.exists(v), f"File not found: {v}"

config["clipper_gff"] = os.path.expanduser(config["clipper_gff"])


rule all:
    input:
        expand(outdir / "alignment_star/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam", sample=SAMPLE_DATA.keys()),
        ## bed = outdir / "peak_cluster/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDupSo.peakClusters.bed"
        [
            outdir / f"peak_cluster/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDupSo.peakClusters.bed"
            for sample in IP_SAMPLES
        ],
        [
            outdir / f"alignment_star/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDupSo.norm.{strand}.bw"
            for sample in SAMPLE_DATA
            for strand in ["pos", "neg"]
        ],
        # sam = outdir / "alignment_repeat/{sample}/{sample}.umi.r1.fqTrTr.fq.Rep.sam"
        [
            outdir / f"alignment_repeat/{sample}/{sample}.umi.r1.fqTrTr.fq.Rep.sam"
            for sample in SAMPLE_DATA
        ],
        [
            outdir / f"alignment_repeat/{sample}/{n1}{n2}.{sample}.umi.r1.fqTrTr.fq.Rep.sam.tmp"
            for n1 in "ACGTN"
            for n2 in "ACGTN"
            for sample in SAMPLE_DATA
        ],
        [
            outdir / f"alignment_star/{sample}/{n1}{n2}.{sample}.umi.r1.fq.genome-mappedSoSo.bam.tmp"
            for n1 in "ACGTN"
            for n2 in "ACGTN"
            for sample in SAMPLE_DATA
        ],
        # normed_peaks_full = outdir / "peak_cluster/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDupSo.normed_peaks.full",
        [
            outdir / f"peaks/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDupSo.normed_peaks.bed"
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
        tmpdir = temp(directory(outdir / "reads_clean/{sample}_tmp")),
        r1 = outdir / "reads_clean/{sample}/{sample}.umi.r1.fqTrTr.sorted.fq.gz"
    resources:
        mem_mb = 4 * 1024
    shell:
        """
        rm -rf {output.tmpdir} && mkdir -p {output.tmpdir} && \
        fastq-sort --id {input.r1} --temporary-directory {output.tmpdir} | gzip > {output.r1}
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
        samtools sort -m 2G -@ {threads} -n {input.bam} > {output.bam} && \
        samtools sort -m 2G -@ {threads} {output.bam} > {output.bam_sorted} && \
        samtools index -@ {threads} {output.bam_sorted}
        """

rule step91_umitools_dedup:
    input:
        bam = outdir / "alignment_star/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.bam"
    output:
        bam = outdir / "alignment_star/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDup.bam",
        bam_sorted = protected(outdir / "alignment_star/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam"),
        bam_sorted_bai = protected(outdir / "alignment_star/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam.bai"),
    params:
        prefix = lambda wildcards: outdir / f"alignment_star/{wildcards.sample}/{wildcards.sample}.umi.r1.fq.genome-mappedSoSo"
    resources:
        mem_mb = 1024 * 48
    threads:
        16
    shell:
        """
        umi_tools dedup \
            --random-seed 1 \
            --method unique \
            --output-stats {params.prefix} \
            -I {input.bam} \
            -S {output.bam} && \
        samtools sort \
            -m 2G \
            -@ {threads} \
            -o {output.bam_sorted} \
            {output.bam} && \
        samtools index -@ {threads} {output.bam_sorted}
        """



# 

# CHROMOSOMES = ["chr{}".format(c) for c in range(1, 23)] + ["chrX", "chrY", "chrM"]
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

CLIPPER_GFF_SUFFIX = ".AS.STRUCTURE.COMPILED.gff"
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
    priority: -10
    input:
        bam = outdir / "alignment_star/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam",
        gff = config["clipper_gff"].replace(CLIPPER_GFF_SUFFIX, "_line-r{nr}" + CLIPPER_GFF_SUFFIX)
    params:
        species = lambda wildcards: os.path.basename(config["clipper_gff"]).replace(CLIPPER_GFF_SUFFIX, "") + f"_line-r{wildcards.nr}",
        clipper_img = BIN["clipper_img"],
        ref_path = CLIPPER_DATA
    output:
        partial_bed = protected(outdir / "peak_cluster/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDupSo.peakClusters.line-r{nr}.bed")
    log:
        protected(outdir / "peak_cluster/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDupSo.peakClusters.line-r{nr}.log")
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
            --bam {input.bam} \
            --outfile {output.partial_bed} \
            --processors {threads} &> {log} && touch {output.partial_bed}
        """
rule step92_merge_peak_clusters:
    priority: -10
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
        bed = protected(outdir / "peak_cluster/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDupSo.peakClusters.bed")
    log:
        protected(outdir / "peak_cluster/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDupSo.peakClusters.log")
    resources:
        mem_mb = 1024 * 16
    threads:
        16
    shell:
        """
        cat {input.bed} | sort -k1,1 -k2,2n -k3,3n --parallel {threads} -S 12G > {output} && \
        cat {input.log} > {log}
        """


rule step93_overlap_peakfi_with_bam:
    priority: -10
    wildcard_constraints:
        ip_sn = "|".join(IP_SAMPLES),
    input:
        ip_bam = outdir / "alignment_star/{ip_sn}/{ip_sn}.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam",
        in_bam = lambda wildcards: outdir / f"alignment_star/{IP2IN[wildcards.ip_sn]}/{IP2IN[wildcards.ip_sn]}.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam",
    output:
        normed_bed = outdir / "peak_cluster/{ip_sn}/{ip_sn}.umi.r1.fq.genome-mappedSoSo.rmDupSo.peakClusters.normed.bed",
        compressed_bed = outdir / "peak_cluster/{ip_sn}/{ip_sn}.umi.r1.fq.genome-mappedSoSo.rmDupSo.peakClusters.normed.compressed.bed"
    resources:
        mem_mb = 1024 * 8
    threads:
        4
    log:
        outdir / "peak_cluster/{ip_sn}/{ip_sn}.umi.r1.fq.genome-mappedSoSo.rmDupSo.peakClusters.normed.log"
    shell:
        """
        samtools view -@ {threads} -cF 4 {input.ip_bam} > {input.ip_bam}.readnum.txt && \
        samtools view -@ {threads} -cF 4 {input.in_bam} > {input.in_bam}.readnum.txt && \
        perl overlap_peakfi_with_bam.pl \
            {input.ip_bam} \
            {input.in_bam} \
            {input.peak_bed} \
            {input.ip_bam}.readnum.txt \
            {input.in_bam}.readnum.txt \
            {output.normed_bed}
        perl compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat.pl \
            {output.normed_bed} \
            {output.compressed_bed}
        """

rule step94_make_bigwig:
    input:
        bam = outdir / "alignment_star/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam",
    output:
        pos_bw = outdir / "alignment_star/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDupSo.norm.pos.bw",
        neg_bw = outdir / "alignment_star/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDupSo.norm.neg.bw"
    params:
        makebigwigfiles_img = BIN["makebigwigfiles_img"],
        genome = config["chrom_sizes"] # chrNameLength.txt
    resources:
        mem_mb = 1024 * 8
    threads:
        1
    shell:
        """
        singularity run \
            --bind /project2/mengjiechen/kenchen/db:/mnt/db \
            {params.makebigwigfiles_img} \
            makebigwigfiles \
            --bw_pos {output.pos_bw} \
            --bw_neg {output.neg_bw} \
            --bam {input.bam} \
            --genome {params.genome} \
            --direction f
        """

rule step95_align_to_repeat_families:
    input:
        fq = outdir / "reads_clean/{sample}/{sample}.umi.r1.fqTrTr.sorted.fq.gz",
    output:
        sam = outdir / "alignment_repeat/{sample}/{sample}.umi.r1.fqTrTr.fq.Rep.sam"
    params:
        repeatitive_element_mapping_img = BIN["repeatitive_element_mapping_img"],
        bowtie2_index = "/home/kenchen/Documents/seCLIP-snakemake/tools/repetitive-element-mapping/repeat-family-mapping-grch38/bowtie2_index/MASTER_FILELIST.20201203.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat",
        wrepbaseandtRNA = "/home/kenchen/Documents/seCLIP-snakemake/tools/repetitive-element-mapping/repeat-family-mapping-grch38/MASTER_FILELIST.20201203.wrepbaseandtRNA.enst2id.fixed.UpdatedSimpleRepeat.wmiRs.tsv"
    resources:
        mem_mb = 1024 * 24
    threads:
        16
    shell:
        """
       /home/kenchen/Documents/seCLIP-snakemake/tools/repetitive-element-mapping/bin/perl/parse_bowtie2_output_realtime_includemultifamily_SE.pl \
            {input.fq} \
            {params.bowtie2_index} \
            {output.sam} \
            {params.wrepbaseandtRNA}
        # singularity run {params.repeatitive_element_mapping_img} \
            #parse_bowtie2_output_realtime_includemultifamily_SE.pl \
        """

rule step96_split_repeat_sam:
    input:
        sam = outdir / "alignment_repeat/{sample}/{sample}.umi.r1.fqTrTr.fq.Rep.sam"
    output:
        sam = [
            outdir / f"alignment_repeat/{{sample}}/{n1}{n2}.{{sample}}.umi.r1.fqTrTr.fq.Rep.sam.tmp"
            for n1 in "ACGTN"
            for n2 in "ACGTN"
        ]
    params:
        outdir = lambda wildcards: outdir / f"alignment_repeat/{wildcards.sample}",
        repeatitive_element_mapping_img = BIN["repeatitive_element_mapping_img"]
    resources:
        mem_mb = 1024 * 8
    threads:
        1
    shell:
        """
        singularity run {params.repeatitive_element_mapping_img} \
            split_bam_to_subfiles_SEorPE.pl \
            {input.sam} \
            SE && \
        mv -f [ACGTN][ACGTN].{wildcards.sample}.umi.r1.fqTrTr.fq.Rep.sam.tmp {params.outdir}
        """

rule step97_split_genome_bam:
    input:
        bam = outdir / "alignment_star/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.bam"
    output:
        bam = [
            outdir / f"alignment_star/{{sample}}/{n1}{n2}.{{sample}}.umi.r1.fq.genome-mappedSoSo.bam.tmp"
            for n1 in "ACGTN"
            for n2 in "ACGTN"
        ]
    params:
        outdir = lambda wildcards: outdir / f"alignment_star/{wildcards.sample}",
        repeatitive_element_mapping_img = BIN["repeatitive_element_mapping_img"]
    resources:
        mem_mb = 1024 * 8
    threads:
        1
    shell:
        """
        singularity run {params.repeatitive_element_mapping_img} \
            split_bam_to_subfiles_SEorPE.pl \
            {input.bam} \
            SE && \
        mv -f [ACGTN][ACGTN].{wildcards.sample}.umi.r1.fq.genome-mappedSoSo.bam.tmp {params.outdir}
        """

rule step98_merge_repeat_genome_sams:
    input:
        rep_sam = outdir / "alignment_repeat/{sample}/{prefix}.{sample}.umi.r1.fqTrTr.fq.Rep.sam.tmp",
        gen_bam = outdir / "alignment_star/{sample}/{prefix}.{sample}.umi.r1.fq.genome-mappedSoSo.bam.tmp"
    params:
        gtf = "/home/kenchen/Documents/seCLIP-snakemake/tools/repetitive-element-mapping/repeat-family-mapping-grch38/gencode.v33.chr_patch_hapl_scaff.annotation.gtf",
        ucsc_table = "/home/kenchen/Documents/seCLIP-snakemake/tools/repetitive-element-mapping/repeat-family-mapping-grch38/gencode.v33.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat.tsv",
        unique_genomic_elements = "/home/kenchen/Documents/seCLIP-snakemake/tools/repetitive-element-mapping/repeat-family-mapping-grch38/UniqueGenomicElements.hg38.bed",
        repbase_trna_table = "/home/kenchen/Documents/seCLIP-snakemake/tools/repetitive-element-mapping/repeat-family-mapping-grch38/MASTER_FILELIST.20201203.wrepbaseandtRNA.enst2id.fixed.UpdatedSimpleRepeat.wmiRs.tsv"
    output:
        sam_rmdup = outdir / "alignment_repeat/{sample}/{prefix}.{sample}.umi.r1.fqTrTr.fq.Rep.sam.tmp.combined_w_uniquemap.rmDup.sam",
        sam_prerm = outdir / "alignment_repeat/{sample}/{prefix}.{sample}.umi.r1.fqTrTr.fq.Rep.sam.tmp.combined_w_uniquemap.prermDup.sam",
        parsed = outdir / "alignment_repeat/{sample}/{prefix}.{sample}.umi.r1.fqTrTr.fq.Rep.sam.tmp.combined_w_uniquemap.rmDup.sam.parsed",
    resources:
        mem_mb = 1024 * 8
    threads:
        1
    shell:
        """
        /home/kenchen/Documents/seCLIP-snakemake/tools/repetitive-element-mapping/bin/perl/duplicate_removal.pl \
            {input.rep_sam} \
            {input.gen_bam} \
            SE \
            {params.gtf} \
            {params.ucsc_table} \
            {params.unique_genomic_elements} \
            {params.repbase_trna_table}
        # singularity run {params.repeatitive_element_mapping_img} \
            # duplicate_removal.pl \
        """


rule step99_merge_multiple_parsed_files:
    input:
        parsed = [
            outdir / "alignment_repeat/{sample}/{prefix}.{sample}.umi.r1.fqTrTr.fq.Rep.sam.tmp.combined_w_uniquemap.rmDup.sam.parsed"
            for prefix in [f"{n1}{n2}" for n1 in "ACGTN" for n2 in "ACGTN"]
        ]
    output:
        parsed = outdir / "alignment_repeat/{sample}/{sample}.umi.r1.fqTrTr.fq.Rep.sam.tmp.combined_w_uniquemap.rmDup.sam.parsed"
    resources:
        mem_mb = 1024 * 8
    threads:
        1
    shell:
        """
        merge_multiple_parsed_files.simplified_20191022.pl \
            {output.parsed} \
            {input.parsed}
        """

rule step100_overlap_peakfi_with_bam:
    wildcard_constraints:
        sample = "|".join(IP_SAMPLES)
    input:
        clip_bam = outdir / "alignment_star/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDupSo.bam",
        input_bam = outdir / "alignment_star/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDup.bam",
        clip_bed = outdir / "peak_cluster/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDupSo.peakClusters.bed"
    output:
        normedbed = outdir / "peaks/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDupSo.normed_peaks.bed",
        # normedbed_full = outdir / "peak_cluster/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDupSo.normed.full.bed"
    params:
        eclip_img = BIN.get("eclip_perl_img")
    resources:
        mem_mb = 1024 * 8
    threads:
        4
    shell:
        """
        samtools view -@ {threads} -c -F 4 {input.clip_bam} > {input.clip_bam}.readnum && \
        samtools view -@ {threads} -c -F 4 {input.input_bam} > {input.input_bam}.readnum && \
        singularity run {params.eclip_img} overlap_peakfi_with_bam.pl \
            {input.clip_bam} \
            {input.input_bam} \
            {input.clip_bed} \
            {input.clip_bam}.readnum \
            {input.input_bam}.readnum \
            {output.normedbed}
        """

rule step101_compress_peaks:
    wildcard_constraints:
        sample = "|".join(IP_SAMPLES)
    input:
        bed = outdir / "peaks/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDupSo.normed_peaks.bed"
    output:
        compressed_bed = outdir / "peaks/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDupSo.normed_peaks.compressed.bed",
        compressed_bed_full = outdir / "peaks/{sample}/{sample}.umi.r1.fq.genome-mappedSoSo.rmDupSo.normed_peaks.compressed.bed.full"
    resources:
        mem_mb = 1024 * 8
    threads:
        1
    shell:
        """
        compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat_outputfull.pl \
            {input.bed}.full \
            {output.compressed_bed} \
            {output.compressed_bed_full}
        """
