#!/bin/bash

module load snakemake/7.7.0

snakemake all --delete-all-output -s TETranscripts.Snakefile --cores all
