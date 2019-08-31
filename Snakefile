import pandas as pd

configfile: "cfg/config.yaml"
configfile: "cfg/cluster.yaml"

sample_names = pd.read_table("sample_names.tsv").set_index("samples", drop=False)


rule all:
    input:
        'reference/{}.fasta.amb'.format(config["genome"]),
        'reference/{}.fasta.ann'.format(config["genome"]),
        'reference/{}.fasta.bwt'.format(config["genome"]),
        'reference/{}.fasta.pac'.format(config["genome"]),
        'reference/{}.fasta.sa'.format(config["genome"]),
        #expand("sorted_reads/{sample}.bam.bai", sample=sample_names.index),
        #expand("vardict_calls/{sample}.vcf", sample=sample_names.index)
        expand("vep_annotate/{sample}.vcf", sample=sample_names.index)

       
rule bwa_index:
    input:
        reference = "reference/{genome}.fasta"
    output:
        "reference/{genome}.fasta.amb",
        "reference/{genome}.fasta.ann",
        "reference/{genome}.fasta.bwt",
        "reference/{genome}.fasta.pac",
        "reference/{genome}.fasta.sa"
    log:
        "logs/bwa_index/{genome}.log"
    params:
        algorithm = "bwtsw"
    shell:
        "(bwa index -a {params.algorithm} {input.reference}) 2> {log}"


rule bwa_align:
    input:
        reference = "reference/{}.fasta".format(config["genome"]),
        fastqs = ["fastqs/{sample}_L001_R1_001.fastq.gz", "fastqs/{sample}_L001_R2_001.fastq.gz"]
    output:
        temp("mapped_reads/{sample}.bam")
    params:
        rg = r"@RG\tID:{sample}\tSM:{sample}",
        cores = config["bwa_align"]["cores"],
    log:
        "logs/bwa_mem/{sample}.log"
    shell:
        "(bwa mem -R '{params.rg}' -t {params.cores} {input.reference} {input.fastqs} | "
        "samtools view -Sb - > {output}) 2> {log}"


rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"


rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"


rule vardict:
    input:
        bam = "sorted_reads/{sample}.bam",
        reference = "reference/{}.fasta".format(config["genome"]),
    output:
        "vardict_calls/{sample}.tsv"
    params:
        min_vaf = config["min_vaf"],
        amplicons_bed = config["amplicons"]
    log:
        "logs/vardict/{sample}.log"
    shell:
        "(vardict-java -G {input.reference} -f {params.min_vaf} -N {wildcards.sample} -b {input.bam} TAMSEQ_TP53_amplicons.bed > {output}) 2> {log}"


rule vardict_strand_bias:
    input:
        tsv = "vardict_calls/{sample}.tsv",
    output:
        "vardict_calls/{sample}.strandbias.tsv"
    log:
        "logs/vardict/{sample}.strandbias.log"
    shell:
        "(cat {input.tsv} | teststrandbias.R > {output}) 2> {log}"


rule vardict_vcf:
    input:
        tsv = "vardict_calls/{sample}.strandbias.tsv",
    output:
        "vardict_calls/{sample}.vcf"
    params:
        min_vaf = config["min_vaf"] 
    log:
        "logs/vardict/{sample}.vcf.log"
    shell:
        "(cat {input.tsv} | var2vcf_valid.pl -a -A -N {wildcards.sample} -E -f {params.min_vaf} > {output}) 2> {log}"


rule vep_annotate:
    input:
        vcf = "vardict_calls/{sample}.vcf",
        reference = "reference/{}.fasta".format(config["genome"]),
    output:
        "vep_annotate/{sample}.vcf"
    params:
    log:
        "logs/vep/{sample}.vcf.log"
    shell:
        "(vep --cache --offline --dir_cache ./vep_cache/ --fasta {input.reference} -i {input.vcf} --sift b --polyphen b --symbol --numbers --biotype --total_length --hgvs --exclude_predicted --af_gnomad --format vcf --force_overwrite --vcf --fields Consequence,IMPACT,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,HGVSc,HGVSp,cDNA_position,CDS_position,gnomAD_AF,gnomAD_AFR_AF,gnomAD_AMR_AF,gnomAD_ASJ_AF,gnomAD_EAS_AF,gnomAD_FIN_AF,gnomAD_NFE_AF,gnomAD_OTH_AF,gnomAD_SAS_AF,MaxEntScan_alt,MaxEntScan_diff,MaxEntScan_ref,PICK --flag_pick -o {wildcards.sample}.vcf) 2> {log}"
