#! /bin/bash
#
#Author: Yukun Tan
#

while getopts g:r:i:c:n:o:d:m:k: flag
do
    case "${flag}" in
        r) ref=${OPTARG};;
        i) tumor=${OPTARG};;
        c) normal=${OPTARG};;
        n) n_cpus=${OPTARG};;
        o) output=${OPTARG};;
	d) directory=${OPTARG};;
	g) genome_ref=${OPTARG};;
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

if [[ -z $genome_ref ]]
then
    genome_ref=${directory}/database/hg38/genome
fi

if [[ -z $check ]]
then
    check="True"
fi

novobreak=${directory}/app/novoBreak
bwa=$directory/app/bwa
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
$bwa mem -t12 -M -T0 $genome_ref ssake.ctg.fa > ssake.ctg_hg38.sam


if [[ $mode == "fusion" ]]
then
	### infer fusions
	perl $directory/src/infer_genefusion.pl ssake.ctg_hg38.sam | grep -v KI | grep -v GL | grep -v MT > ssake.ctg_hg38.vcf
	### annotate and generate raw fusions
	python $directory/src/annovar_annotation.py ssake.ctg_hg38.vcf $directory $mode
	grep -v '^#' ssake.ctg_hg38.all.annovar | awk  -v OFS="\t" '{gsub(/\|/, "\t", $6); print}' | sed 's/read//' | perl -ne '$chr2=$1 if /CHR2=(\S+?);/; $pos2=$1 if /END=(\d+);/; print $chr2,"\t",$pos2,"\t",$_' | awk '{if(!x[$1$2$3$4]){y[$1$2$3$4]=$6;x[$1$2$3$4]=$0}else{if($6>y[$1$2$3$4]){y[$1$2$3$4]=$6; x[$1$2$3$4]=$0}}}END{for(i in x){print x[i]}}' | cut -f3- | sort -k1,1 -k2,2n  | perl -ne 'if(/TRA/){print}elsif(/SVLEN=(\d+)/){if($1>50000){print $_}}elsif(/SVLEN=-(\d+)/){if($1>50000){print}}' > ssake.pass.vcf
	python $directory/src/mark_multi_hit.py ssake.pass.vcf $directory
elif [[ $mode == "splice" ]]
then
	## infer splice junctions
	perl $directory/src/infer_splicejunction.pl ssake.ctg_hg38.sam | grep -v KI | grep -v GL | grep -v MT > ssake.ctg_hg38.vcf
	
	## annotate and generate raw splice junctions
	python $directory/src/annovar_annotation.py ssake.ctg_hg38.vcf $directory $mode
	grep -v '^#' ssake.ctg_hg38.all.annovar | awk  -v OFS="\t" '{gsub(/\|/, "\t", $6); print}' | sed 's/read//' | perl -ne '$chr2=$1 if /CHR2=(\S+?);/; $pos2=$1 if /END=(\d+);/; print $chr2,"\t",$pos2,"\t",$_' | awk '{if(!x[$1$2$3$4$5]){y[$1$2$3$4$5]=$6;x[$1$2$3$4$5]=$0}else{if($6>y[$1$2$3$4$5]){y[$1$2$3$4$5]=$6; x[$1$2$3$4$5]=$0}}}END{for(i in x){print x[i]}}' | cut -f3- | sort -k1,1 -k2,2n  | perl -ne 'if(/SVLEN=(\d+)/){if($1>50 and $1<400000){print $_}}elsif(/SVLEN=-(\d+)/){if($1>50 and $1<400000){print}}' > mark.vcf
else
	echo "not a valid mode options (support 'fusion' and 'splice' only)"
fi

grep '^#' ssake.ctg_hg38.all.annovar > header.txt	

if [[ $mode == "fusion" ]]
then
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

elif [[ $mode == "splice" ]]
then
	cat mark.vcf | perl -ne '$chr2=$1 if /CHR2=(\S+?);/; $pos2=$1 if /END=(\d+);/; print $chr2,"\t",$pos2,"\t",$_' | cut -f1,2,3,4,5,6,11 > splice_junction.txt

	for i in 3 2 1 0 -1 -2 -3
	do
		for j in 3 2 1 0 -1 -2 -3
		do
			awk -v i=$i -v j=$j '{$4+=i; $2+=j; if ($2 > $4) print $1"\t"$4"\t"$2-1"\t"".""\t"".""\t"$5"\t"$4"\t"$2-1"\t""255,0,0""\t""2""\t"".,.""\t"".,."; else print $1"\t"$2"\t"$4-1"\t"".""\t"".""\t"$5"\t"$2"\t"$4-1"\t""255,0,0""\t""2""\t"".,.""\t"".,.";}' splice_junction.txt > sj_${i}_${j}.bed
			$directory/app/regtools junctions annotate -S -o sj_${i}_${j}_annotated_novobreak.tsv sj_${i}_${j}.bed $directory/database/hg38/genome.fa $directory/database/hg38/genome.gtf
			sed -i '1d' sj_${i}_${j}_annotated_novobreak.tsv
			paste sj_${i}_${j}_annotated_novobreak.tsv <(cut -f6,7 splice_junction.txt) > sj_new_${i}_${j}_annotated_novobreak.tsv
			flag=$(( ($i<0?-$i:$i)+($j<0?-$j:$j) ))
			sed -i "s/$/\t${flag}/" sj_new_${i}_${j}_annotated_novobreak.tsv
		done
	done
	python $directory/src/merge.py $output
	cat merged_total.csv | sort -u -k1,1 -k2,2n -k3,3n > merged_splice.csv

	num=`wc -l merged_splice.csv | cut -f1 -d' '`
	rec=`echo $num/$n_cpus | bc`
	rec=$((rec+1))
	mkdir filter_split_splice
	cd filter_split_splice
	split -l $rec ../merged_splice.csv # set proper split parameters when needed
	for file in x??
	do
		echo $file
		perl $directory/src/infer_bp_splice.pl $file $tumor $check > $file.sp.vcf &
	done
	wait
	cd ..

	perl $directory/src/filter_splice_junction.pl filter_split_splice/*.sp.vcf | cat header.txt - > temp_splice_pass.tsv
	python $directory/src/merge_junction.py $output
	rm sj_*
fi
cd ..

