# TE pipeline

Snakemake pipeline for running TEcount on [NIH's Biowulf cluster](http://hpc.nih.gov).

Performs adapter trimming, qc, experiment inference, alignment (STAR) and quantifcation of genes/TEs using TEcount from the TEtoolkit.


## Dependencies
	
	snakemake/7.7.0
	cutadapt/1.18
  	fastqc/0.11.9
  	multiqc/1.9
  	python/3.9
  	qualimap/2.2.1
  	rseqc/4.0.0
  	samtools/1.13
  	STAR/2.7.6a
  	tetoolkit/2.1.4
   
## Setup

Dependencies and paths to reference annotations and alignment indices are set in the file `Templates/template_TEtranscripts.json` and will need to be edited to suit your purposes.


## Setting up a run on Biowulf

Fastqs must be in directory FASTQ/

read1 and read2 must end with .R?.gz

Construct a meta.tab file, tab-delimited file containing sample IDs and a treatment index. Samples with same treatment index will be combined in downstream analyses (e.g. differential analysis, not yet implemented in this workflow).

	#samplename<tab>treatment_index
	sample1<tab>1
	sample2<tab>2

Then run:

`bash setup_run.sh <experiment_name>`

This generates the "run.json" file, which contains all of the parameters for the run.


To do a dry run: 

	module load snakemake/7.7.0
	snakemake -n -s TEtranscripts.Snakefile

To submit the workflow on biowulf (NIH HPC):

`sbatch run_analysis.sh`

