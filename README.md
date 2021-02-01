# cutntag
**Nextflow pipeline for Cut&amp;Tag analysis**

## Introduction
This Nextflow pipeline for analysis of Cut&Tag data is based on the Cut&Tag analysis workflow presented in [CUTTag_tutorial](https://github.com/yezhengSTAT/CUTTag_tutorial) with some modifications.
The peak calling is performed by [SEACR](https://github.com/FredHutch/SEACR), intended to call peaks from sparse data.

## Pipeline summary

Download the cutntag_env.yml file and create a conda environemnt that contain all packages neccisary to run the pieline.
```bash
conda create -f cutntag_env.yml
```

Dowload the pipleine and run
```bash
nextflow run main.nf --samples samplex.txt --outdir outdir --bowtie2_index /path/to/bowtie_index/prefix --chromsize /path/to/chrom.sizes --design design.txt
```

## Run the pipeline

### Input

### Options

### Output

### Configuraiton
