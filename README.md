# novoRNABreak

## Description

novoRNABreak is a tool used for novel splice junctions and fusion transcripts detection from RNA-seq data. novoRNABreak is based on the local assembly model, which offers a tradeoff between the alignment-based and the de novo whole transcriptome assembly (WTA) approach, namely, being more sensitive in assembling novel junctions that cannot be directly aligned, and more efficient due to the strategy that focuses on junctions rather than full-length transcripts. 

## System requirements and dependency

novoRNABreak runs on a x86_64 Linux system with a ~40GB physical memory. It depends on picard for duplicates removal, novoBreak for novel k-mer detection, SSAKE for local assembly, bwa-mem for contig mapping, samtools (v1.3 and above) to extract reads, annovar for gene annotation, and regtools for junction annotation. All of those are in app directory.

The human hg38 genome index (bwa) and transcriptome references (k-mer) are stored at https://doi.org/10.5281/zenodo.7735433 Please download and decompress to the novoRNABreak directory before running the pipeline. If you want to run the "splice" mode, please create the genome index using STAR and add the path to "-g" option. 

```
tar -xvzf database.tar.gz -C <novoRNABreak.directory>
STAR --runMode genomeGenerate --genomeDir ref/ --genomeFastaFiles Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa --sjdbGTFfile Homo_sapiens.GRCh38.105.gtf --runThreadN 16
```
The files for --genomeFastaFiles and --sjdbGTFfile can be downloaded at https://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna/ and https://ftp.ensembl.org/pub/release-105/gtf/homo_sapiens/

## Installation

```
git clone https://github.com/KChen-lab/novoRNABreak.git
```

## Usage

### Preprocessing:

Input files: 

- Bam file is required. If alignment is needed, please use STAR with "--outSAMunmapped Within" option to keep all the unmapped sequences, e.g.,
```
STAR --runMode alignReads --genomeDir <hg38.fa> --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate --readFilesIn <fastq1> <fastq2> --runThreadN <n> --outFileNamePrefix <name>
```

- Duplicates should be removed from the bam files using picard, e.g.,
```
java -jar Xmg30g picard.jar MarkDuplicates -I bamfile -O output -M markdup.metrics -AS true -REMOVE_DUPLICATES true
```

### Run run_novobreak_rna.sh

```
bash run_novobreak_rna.sh -m <mode> -i <tumor.bam> -c <normal.bam> -r <reference.transcript> -g <reference.genome> -a <annotation.genome> -n <CUPs> -d <novoRNABreak.directory> -o <output.directory> -k <chromcheck> [options]
Options:
  -m <string>    "fusion"/"splice"
  -i <string>    Tumor bam file
  -c <string>    Normal bam file [optional; default is None]
  -r <string>    Transcript reference file in fasta format (reads larger than 31bp) [optional; included in database]
  -g <string>    Genome reference index file (BWA index for fusion mode; STAR index for splice mode) [optional; included in database]
  -a <string>    Genome annotation gtf file [optional; included in database]
  -n <int>       The number of CUPs
  -d <string>    Path of novoRNABreak directory
  -o <string>    Path of output directory
  -k <Boolean>   "True/False", whether the chromosomes are well formated, e.g., chr1, chr2, ...(True), or 1, 2, ...(False) [optional; default is True]  
```
