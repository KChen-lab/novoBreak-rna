#! /bin/bash
#
#Author: Yukun Tan
#

while getopts g:r:i:c:n:o:d:m: flag
do
    case "${flag}" in
        r) ref=${OPTARG};;
        i) tumor=${OPTARG};;
        c) normal=${OPTARG};;
        n) n_cpus=${OPTARG};;
        o) output=${OPTARG};;
	d) directory=${OPTARG};;
	g) genome_ref=${OPTARG};;
	m) mode=${OPTARG}	
    esac
done

if [[ -z $normal ]]
then
    normal=$directory/source/null.bam
fi

if [[ -z $ref ]]
then
    ref=${directory}/source/mrna.fa
fi

if [[ -z $genome_ref ]]
then
    genome_ref=/rsrch3/scratch/bcb/ytan1/prj/reference/BWA/hg38/genome
fi

novobreak=${directory}/app/novoBreak
bwa=$directory/app/bwa
samtools=$directory/app/samtools

mkdir -p $output
cd $output

$novobreak -i $tumor -c $normal -r ${ref} -o kmer.stat
$samtools view -H somaticreads.bam | $samtools reheader - somaticreads.bam | $samtools sort -n -@ 4 - -o somaticreads.srtnm.bam
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
$bwa mem -t12 -M -T0 $genome_ref ssake.ctg.fa > ssake.ctg_hg38.sam


if [[ $mode == "fusion" ]]
then
	## infer fusions
	perl $directory/src/infer_genefusion.pl ssake.ctg_hg38.sam | grep -v KI | grep -v GL | grep -v MT > ssake.ctg_hg38.vcf
	## annotate and generate raw fusions
	python $directory/src/annovar_annotation.py ssake.ctg_hg38.vcf $directory
	grep -v '^#' ssake.ctg_hg38.all.annovar | awk  -v OFS="\t" '{gsub(/\|/, "\t", $6); print}' | sed 's/read//' | perl -ne '$chr2=$1 if /CHR2=(\S+?);/; $pos2=$1 if /END=(\d+);/; print $chr2,"\t",$pos2,"\t",$_' | awk '{if(!x[$1$2$3$4]){y[$1$2$3$4]=$6;x[$1$2$3$4]=$0}else{if($6>y[$1$2$3$4]){y[$1$2$3$4]=$6; x[$1$2$3$4]=$0}}}END{for(i in x){print x[i]}}' | cut -f3- | sort -k1,1 -k2,2n  | perl -ne 'if(/TRA/){print}elsif(/SVLEN=(\d+)/){if($1>50000){print $_}}elsif(/SVLEN=-(\d+)/){if($1>50000){print}}' > ssake.pass.vcf
	python $directory/src/mark_multi_hit.py ssake.pass.vcf $directory
elif [[ $mode == "splice" ]]
then
	## infer splice junctions
	perl $directory/src/infer_splicejunction.pl ssake.ctg_hg38.sam | grep -v KI | grep -v GL | grep -v MT > ssake.ctg_hg38.vcf
	
	## annotate and generate raw splice junctions
	python $directory/src/annovar_annotation.py ssake.ctg_hg38.vcf $directory
	grep -v '^#' ssake.ctg_hg38.all.annovar | awk  -v OFS="\t" '{gsub(/\|/, "\t", $6); print}' | sed 's/read//' | perl -ne '$chr2=$1 if /CHR2=(\S+?);/; $pos2=$1 if /END=(\d+);/; print $chr2,"\t",$pos2,"\t",$_' | awk '{if(!x[$1$2$3$4]){y[$1$2$3$4]=$6;x[$1$2$3$4]=$0}else{if($6>y[$1$2$3$4]){y[$1$2$3$4]=$6; x[$1$2$3$4]=$0}}}END{for(i in x){print x[i]}}' | cut -f3- | sort -k1,1 -k2,2n  | perl -ne 'if(/SVLEN=(\d+)/){if($1>50 and $1<400000){print $_}}elsif(/SVLEN=-(\d+)/){if($1>50 and $1<400000){print}}' > ssake.pass.vcf
else
	echo "not a valid mode options (support 'fusion' and 'splice' only)"
fi

num=`wc -l mark.vcf | cut -f1 -d' '`
rec=`echo $num/$n_cpus | bc`
rec=$((rec+1))
mkdir filter_split
cd filter_split
split -l $rec ../mark.vcf # set proper split parameters when needed
for file in x??
do
	echo $file
	perl $directory/src/infer_bp_v4.pl $file $tumor > $file.sp.vcf &
done
wait
cd ..

#below is a naive filter, pay attention to it
grep '^#' ssake.ctg_hg38.all.annovar > header.txt	

if [[ $mode == "fusion" ]]
then
	perl $directory/src/filter_sv_icgc.pl filter_split/*.sp.vcf | cat header.txt - > novoBreak.rna.pass.vcf 
	perl $directory/src/full_filter_sv_icgc.pl filter_split/*.sp.vcf | cat header.txt - > full_novoBreak.rna.pass.vcf 
	python $directory/src/fusion_genes.py $output
elif [[ $mode == "splice" ]]
then
	perl $directory/src/filter_splice_junction.pl filter_split/*.sp.vcf | cat header.txt - > novoBreak.rna.pass.vcf
fi
cd ..
