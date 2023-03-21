#! /bin/bash
#
#Author: Yukun Tan
#

while getopts g:r:i:c:n:o:d:m:k:a: flag
do
    case "${flag}" in
        r) ref=${OPTARG};;
        i) tumor=${OPTARG};;
        c) normal=${OPTARG};;
        n) n_cpus=${OPTARG};;
        o) output=${OPTARG};;
	d) directory=${OPTARG};;
	g) genome_ref=${OPTARG};;
	a) genome_annotation=${OPTARG};;
	m) mode=${OPTARG};;
	k) check=${OPTARG}
    esac
done

if [[ -z $normal ]]
then
    normal=$directory/source/null.bam
fi

if [[ -z $ref ]]
then
    ref=${directory}/database/mrna.fa
fi

if [[ -z $genome_ref && $mode == "fusion" ]]
then
    genome_ref=${directory}/database/hg38_bwa/genome
elif [[ -z $genome_ref && $mode == "splice" ]]
then
    genome_ref=${directory}/database/hg38_star
fi

if [[ -z $genome_annotation ]]
then
    genome_annotation=$directory/database/hg38_bwa/genome.gtf
fi

if [[ -z $check ]]
then
    check="True"
fi

novobreak=${directory}/app/novoBreak
bwa=$directory/app/bwa
STAR=$directory/app/STAR
STARlong=$directory/app/STARlong
samtools=$directory/app/samtools

mkdir -p $output
cd $output

$novobreak -i $tumor -c $normal -r ${ref} -o kmer.stat
$samtools view -H somaticreads.bam | $samtools reheader - somaticreads.bam | $samtools sort -n -@ $n_cpus - -o somaticreads.srtnm.bam
$samtools bam2fq -1 read1.fq -2 read2.fq  somaticreads.srtnm.bam

#clustering
perl $directory/src/group_bp_reads_dev.pl kmer.stat read1.fq read2.fq > bp_reads.txt

#split clusters and do de novo assembly. For parallel purpose
mkdir split
cd split
awk '{print > int($1/10000)".txt" }' ../bp_reads.txt
for file in *.txt; do perl $directory/src/run_ssake2.pl $file 2\>\&1 \> /dev/null ; done
cd ..

# collect contigs
cat split/*.asm.out | awk 'length($1)>1' > ssake.ctg.fa


# align contigs to reference
if [[ $mode == "fusion" ]]
then
	$bwa mem -t12 -M -T0 $genome_ref ssake.ctg.fa > ssake.ctg_hg38.sam
elif [[ $mode == "splice" ]]
then
	python $directory/src/split_fasta.py $output
	$STAR --runMode alignReads --genomeDir $genome_ref --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate --readFilesIn ssake.ctg.short.fa --runThreadN $n_cpus --outFileNamePrefix all --sjdbGTFfile $genome_annotation --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --clip3pAdapterSeq AAAAAAAAAAAAAAA --alignSJoverhangMin 8 --alignSJDBoverhangMin 3
	awk '{ if ($1 ~ /^([1-9][0-9]?|X|Y)$/) { print } }' allSJ.out.tab | awk '$6 == 0' > allSJ_filtered.out.tab
	$STAR --runMode alignReads --genomeDir $genome_ref --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate --readFilesIn ssake.ctg.short.fa --runThreadN $n_cpus --outFileNamePrefix pass2all --sjdbGTFfile $genome_annotation --sjdbFileChrStartEnd allSJ_filtered.out.tab --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --clip3pAdapterSeq AAAAAAAAAAAAAAA --alignSJoverhangMin 8 --alignSJDBoverhangMin 3
	samtools view -f 4 pass2allAligned.sortedByCoord.out.bam | awk '{print $1}' > unmapped.txt
	python $directory/src/extract_fasta.py $output
	cat ssake.ctg.unmap.fa ssake.ctg.long.fa > ssake.ctg.rest.fa
	$STARlong --runMode alignReads --genomeDir $genome_ref --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate --readFilesIn ssake.ctg.rest.fa --runThreadN $n_cpus --outFileNamePrefix unmap --outFilterMultimapScoreRange 20 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0.66 --outFilterMismatchNmax 1000 --winAnchorMultimapNmax 200 --seedSearchStartLmax 12 --seedPerReadNmax 100000 --seedPerWindowNmax 100 --alignTranscriptsPerReadNmax 100000 --alignTranscriptsPerWindowNmax 10000 --sjdbGTFfile $genome_annotation
	awk '{ if ($1 ~ /^([1-9][0-9]?|X|Y)$/) { print } }' unmapSJ.out.tab | awk '$6 == 0' > unmapSJ_filtered.out.tab
	$STARlong --runMode alignReads --genomeDir $genome_ref --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate --readFilesIn ssake.ctg.rest.fa --runThreadN $n_cpus --outFileNamePrefix pass2unmap --outFilterMultimapScoreRange 20 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0.66 --outFilterMismatchNmax 1000 --winAnchorMultimapNmax 200 --seedSearchStartLmax 12 --seedPerReadNmax 100000 --seedPerWindowNmax 100 --alignTranscriptsPerReadNmax 100000 --alignTranscriptsPerWindowNmax 10000 --sjdbGTFfile $genome_annotation --sjdbFileChrStartEnd unmapSJ_filtered.out.tab
	cat pass2unmapSJ.out.tab pass2allSJ.out.tab | sort -u -k1,1n -k2,2n -k3,3n > final_SJ.out.tab
	awk '{ if ($1 ~ /^([1-9][0-9]?|X|Y)$/) { print } }' final_SJ.out.tab > final_filtered.out.tab
fi


if [[ $mode == "fusion" ]]
then
	### infer fusions
	perl $directory/src/infer_genefusion.pl ssake.ctg_hg38.sam | grep -v KI | grep -v GL | grep -v MT > ssake.ctg_hg38.vcf
	### annotate and generate raw fusions
	python $directory/src/annovar_annotation.py ssake.ctg_hg38.vcf $directory $mode
	grep -v '^#' ssake.ctg_hg38.all.annovar | awk  -v OFS="\t" '{gsub(/\|/, "\t", $6); print}' | sed 's/read//' | perl -ne '$chr2=$1 if /CHR2=(\S+?);/; $pos2=$1 if /END=(\d+);/; print $chr2,"\t",$pos2,"\t",$_' | awk '{if(!x[$1$2$3$4]){y[$1$2$3$4]=$6;x[$1$2$3$4]=$0}else{if($6>y[$1$2$3$4]){y[$1$2$3$4]=$6; x[$1$2$3$4]=$0}}}END{for(i in x){print x[i]}}' | cut -f3- | sort -k1,1 -k2,2n  | perl -ne 'if(/TRA/){print}elsif(/SVLEN=(\d+)/){if($1>50000){print $_}}elsif(/SVLEN=-(\d+)/){if($1>50000){print}}' > ssake.pass.vcf
	python $directory/src/mark_multi_hit.py ssake.pass.vcf $directory
	grep '^#' ssake.ctg_hg38.all.annovar > header.txt	
	num=`wc -l mark.vcf | cut -f1 -d' '`
	rec=`echo $num/$n_cpus | bc`
	rec=$((rec+1))
	mkdir filter_split
	cd filter_split
	split -l $rec ../mark.vcf # set proper split parameters when needed
	for file in x??
	do
		echo $file
		perl $directory/src/infer_bp_fusion.pl $file $tumor $check > $file.sp.vcf &
	done
	wait
	cd ..

	perl $directory/src/filter_fusion.pl filter_split/*.sp.vcf | cat header.txt - > fusion_novoRNABreak_pass.tsv
	python $directory/src/fusion_genes.py $output

fi
cd ..

