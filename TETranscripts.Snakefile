from os.path import join
import re,os
from os import listdir

report: "Reports/workflow.rst"

configfile: "run.json"
workpath = config['project']['workpath']
exp_name= config['project']['experiment_name']

# create output directories

fastq_dir='FASTQ'
trim_dir='trim'
bam_dir='bam'
tecounts_dir='counts'
qc_dir='QC'

# define samples inputs 
samples, = glob_wildcards(join(fastq_dir, '{sample}.R1.fastq.gz'))

for d in [trim_dir,bam_dir,tecounts_dir,qc_dir]:
    if not os.path.exists(join(workpath,d)):
        os.mkdir(join(workpath,d))

rule all:
    input:
        expand(join(workpath,"multiqc_report.html")),     
        expand(join(workpath,tecounts_dir,"{name}.cntTable"),name=samples),

rule trim:
    input:
        file1=join(workpath,fastq_dir,"{name}.R1.fastq.gz"),
        file2=join(workpath,fastq_dir,"{name}.R2.fastq.gz"),
    output:
        outfq1=join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),
        outfq2=join(workpath,trim_dir,"{name}.R2.trim.fastq.gz"),
    params:
        rname="TETranscripts:trim",
        cutadaptver=config['bin']['CUTADAPTVER'],
        adapters=config['references']['FASTAWITHADAPTERSETD']
    threads: 16
    shell:"""
    module load {params.cutadaptver};
    cutadapt --pair-filter=any \
    --nextseq-trim=2 --trim-n -n 5 -O 5 \
    -q 10,10 -m 35:35 -j {threads} \
    -b {params.adapters} -B {params.adapters} \
    -o {output.outfq1} -p {output.outfq2} \
    {input.file1} {input.file2}
    """

rule fastqc:
    input:
        expand(join(workpath,trim_dir,"{name}.R{rn}.trim.fastq.gz"), name=samples,rn=[1,2]),
    output:
        expand(join(workpath,qc_dir,"{name}.R{rn}.trim_fastqc.html"), name=samples,rn=[1,2]),
    params:
        fastqcver=config['bin']['FASTQCVER'],
        rname="TETranscripts:fastqc",
    threads: 12
    shell:"""
    module load {params.fastqcver};
    fastqc {input} -t {threads} -o {qc_dir}
    """

rule multiqc:
    input:
        expand(join(workpath,qc_dir,"{name}.R{rn}.trim_fastqc.html"), name=samples,rn=[1,2]),
        expand(join(workpath,bam_dir,"{name}.Aligned.out.sorted.bam"), name=samples),
        expand(join(workpath,qc_dir,"{name}.infer_exp.txt"),name=samples),
        expand(join(workpath,qc_dir,"qualimap.{name}","qualimap.html"),name=samples),
    params:
        rname="multiqc",
        multiqc=config['bin']['MULTIQCVER'],
    output:
        join(workpath,"multiqc_report.html")
    shell:"""
    module load {params.multiqc};
    multiqc -f {trim_dir} {bam_dir} {qc_dir}
    """

rule align:
    input:
        file1=join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),
        file2=join(workpath,trim_dir,"{name}.R2.trim.fastq.gz"),
    output:
        out1=join(workpath,bam_dir,"{name}.STAR.Aligned.out.bam"),
        out2=join(workpath,bam_dir,"{name}.Aligned.out.sorted.bam"),
        out3=join(workpath,bam_dir,"{name}.Aligned.out.sorted.bam.flagstat"),
        out4=join(workpath,bam_dir,"{name}.Aligned.out.sorted.bam.idxstats")
    params:
        rname='TETranscripts:align',
        index=config['references']['STARINDEX'],
        starver=config['bin']['STARVER'],
        samtoolsver=config['bin']['SAMTOOLSVER'],
        genesgtf=config['references']['GENESGTF'],
    threads: 32
    shell:"""
    module load {params.starver};
    STAR --runThreadN 32 --genomeDir {params.index} \
    --sjdbGTFfile {params.genesgtf} \
    --sjdbOverhang 74 \
    --readFilesIn {input.file1} {input.file2} \
    --outMultimapperOrder Random \
    --runRNGseed 777 \
    --readFilesCommand zcat \
    --outSAMtype BAM Unsorted \
    --winAnchorMultimapNmax 100 \
    --outFilterMultimapNmax 100 \
    --outFileNamePrefix {workpath}/{bam_dir}/{wildcards.name}.STAR.
    module load {params.samtoolsver}
    samtools sort -@{threads} {output.out1} -o {output.out2}
    samtools index {output.out2}
    samtools flagstat {output.out2} >{output.out3}
    samtools idxstats {output.out2} >{output.out4}
    """

rule qualimap:
    input:
        bam=join(workpath,bam_dir,"{name}.Aligned.out.sorted.bam"),
    output:
        join(workpath,qc_dir,"qualimap.{name}","qualimap.html"),
        outdir=join(workpath,qc_dir,"qualimap.{name}")
    params:
        qualimapver=config['bin']['QUALIMAPVER'],
        genesgtf=config['references']['GENESGTF'],
        rname='TETranscripts:qualimap',
    shell:"""
    unset DISPLAY
    module load {params.qualimapver};
    mkdir -p {output.outdir};
    qualimap rnaseq \
    -outdir {output.outdir} \
    -a proportional \
    -bam {input.bam} \
    -gtf {params.genesgtf} \
    -pe \
    --java-mem-size=12G \
    """

rule inferstrand:
    input:
        bam=join(workpath,bam_dir,"{name}.STAR.Aligned.out.bam"),
    output:
        out1=join(workpath,qc_dir,"{name}.infer_exp.txt"),
        out2=temp(join(workpath,qc_dir,"{name}.exptype.txt"))
    params:
        rnaseqcver=config['bin']['RSEQCVER'],
        genesbed=config['references']['GENESBED'],
        scriptsdir=config['bin']['SCRIPTSDIR'],
        pythonver=config['bin']['PYTHONVER'],
        rname='TETranscripts:inferstrand'
    shell:"""
    module load {params.rnaseqcver}
    infer_experiment.py -i {input.bam} -s 400000 -r {params.genesbed} >{output.out1}
    module load {params.pythonver}
    python {params.scriptsdir}/get_strandedness_for_tecount.py -i {output.out1} >{output.out2}
    """

rule tecount:
    input:
        bam=join(workpath,bam_dir,"{name}.Aligned.out.sorted.bam"),
        exptype=join(workpath,qc_dir,"{name}.exptype.txt")
    output:
        join(workpath,tecounts_dir,"{name}.cntTable"),
    params:
        tetookitver=config['bin']['TETOOLKITVER'],
        genesgtf=config['references']['GENESGTF'],
        repeatsgtf=config['references']['TEGTF'],
        rname='TETranscripts:tecount'
    shell:"""
    exptype=$(head -n 1 {input.exptype})
    module load {params.tetookitver}
    TEcount \
    -b {input.bam} \
    --GTF {params.genesgtf} \
    --TE {params.repeatsgtf} \
    --format BAM \
    --mode multi \
    --project {workpath}/{tecounts_dir}/{wildcards.name} \
    --stranded $exptype
    """