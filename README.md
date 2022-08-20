# novoBreak-rna

## Description

novoBreak-rna is a tool used for novel splice junctions and fusion transcripts detection from RNA-seq data. novoBreak-rna is based on the local assembly model, which offers a tradeoff between the alignment-based and the de novo whole transcriptome assembly (WTA) approach, namely, being more sensitive in assembling novel junctions that cannot be directly aligned, and more efficient due to the strategy that focuses on junctions rather than full-length transcripts. 

## System requirements and dependency

novoBreak-rna runs on a x86_64 Linux system with a ~40GB physical memory. It depends on picard for duplicates removal, novoBreak for novel k-mer detection, SSAKE for local assembly, bwa-mem for contig mapping, samtools (v1.3 and above) to extract reads, and annovar for gene annotation. All of those are in app directory.

## Installation

```
git clone https://github.com/KChen-lab/novoBreak-rna.git
```

## Usage

### Preprocessing:

Input files: 

- Bam file is required. If alignment is required, please use STAR with "--outSAMunmapped Within" option to keep all the unmapped sequences, e.g.,
```
STAR --runMode alignReads --genomeDir <hg38.fa> --outSAMunmapped Within --twopassMode Basic --outSAMtype BAM SortedByCoordinate --readFilesIn <fastq1> <fastq2> --runThreadN <n> --outFileNamePrefix <name>
```

- Duplicates should be removed from the bam files using picard, e.g.,
```
java -jar Xmg30g picard.jar MarkDuplicates -I bamfile -O output -M markdup.metrics -AS true -REMOVE_DUPLICATES true
```

### Runrun_novobreak_rna.sh

```
bash run_novobreak_rna.sh -m <mode> -i <tumor.bam> -c <normal.bam> -r <reference.transcript> -g <reference.genome> -n <CUPs> -d <novoBreak-rna.directory> -o <output.directory> [options]
Options:
  -m <string>    "Fusion"/"Splice"
  -i <string>    Tumor bam file
  -c <string>    Normal bam file [optional; default is None]
  -r <string>    Transcript reference file in fasta format (reads larger than 31bp) [optional]
  -g <string>    Genome reference file in fasta format (BWA) [optional]
  -n <int>       The number of CUPs
  -d <string>    Path of novoBreak-rna directory
  -o <string>    Path of output directory
```
