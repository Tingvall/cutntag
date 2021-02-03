# cutntag
**Nextflow pipeline for Cut&amp;Tag analysis**

## Introduction
This Nextflow pipeline for analysis of Cut&Tag data is based on the Cut&Tag analysis workflow presented in [CUTTag_tutorial](https://github.com/yezhengSTAT/CUTTag_tutorial) with some modifications.
The peak calling is performed by [SEACR](https://github.com/FredHutch/SEACR), intended to call peaks from sparse data.

## Pipeline summary

The pipeline consist of the following processes:

1. Unzip and Merge fastq
2. Reads QC - Untrimmed ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
3. Read trimming ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)) & FastQC of trimmed reads
4. Alignment ([`Bowtie2`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
5. Duplicate removal ([`picard`](https://broadinstitute.github.io/picard/))
6. Sample stats & plots
7. Filtering
8. Sample correlation
9. Peak calling (`SEACR](https://github.com/FredHutch/SEACR))
10. Peak calling stats & plots
11. BigWig generation

## Run the pipeline

Download the cutntag_env.yml file and create a conda environemnt that contain all packages neccisary to run the pieline.
```bash
conda create -f cutntag_env.yaml
```

Dowload the pipleine and run
```bash
nextflow run main.nf --samples samplex.txt --outdir outdir --bowtie2_index /path/to/bowtie_index/prefix --chromsize /path/to/chrom.sizes --design design.txt
```

### Input

### Options

### Output

### Configuraiton
