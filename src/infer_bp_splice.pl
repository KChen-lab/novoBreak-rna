#!/usr/bin/perl -w
#
#Author: Yukun Tan
#
use strict;
use warnings;
use English;

my $vcf = shift or die "Usage: $0 <vcf> <tumor.bam> <bool>\n";
my $tbam = shift or die "Usage: $0 <vcf> <tumor.bam> <bool>\n";
my $bool = shift or die  "Usage: $0 <vcf> <tumor.bam> <bool>\n";

open IN, $vcf or die $!;
while (<IN>) {
	print if /^#/;
	next if /^#/;
	chomp;
	my $line = $_;
	my @e = split /\s+/, $_;
	next if $e[14] eq 'mark';
	my $fh;
	my $beg = $e[1]-500;
	my $fin = $e[1]+500;
	my $chr1 = $e[0];
	my $beg2 = $e[2]-500;
	my $fin2 = $e[2]+500;
	$chr1=~ tr/chr//d;
	if ($bool eq 'True'){
		system("samtools", "view", $tbam, "chr".$chr1.":".$beg."-".$fin, "-o", "$PID.regionsam2.tmp");
	}else{
		system("samtools", "view", $tbam, $chr1.":".$beg."-".$fin, "-o", "$PID.regionsam2.tmp");
	}
	open $fh , "$PID.regionsam2.tmp" or die $!;
	my @pos1_tum = count_sp($fh, $e[1]);
	seek $fh, 0, 0;
	close $fh;
	if ($bool eq 'True'){
		system("samtools", "view", $tbam, "chr".$chr1.":".$beg2."-".$fin2, "-o", "$PID.regionsam2.tmp");
	}else{
		system("samtools", "view", $tbam, $chr1.":".$beg2."-".$fin2, "-o", "$PID.regionsam2.tmp");
	}
	open $fh , "$PID.regionsam2.tmp" or die $!;
	my @pos2_tum = count_sp($fh, $e[2]);
	seek $fh, 0, 0;
	close $fh;
	print join("\t", ($line, @pos1_tum, @pos2_tum)), "\n";
}
close IN;

1;

sub count_sp {
	my ($fd, $target) = @_;
	my $cnt = 0;
	my $qual = 0;
	my $cnt2 = 0;
	my $qual2 = 0;
        my $cnt3 = 0;
        my $qual3 = 0;
	my $total = 0;
	while (<$fd>) {
		chomp;
		my @e = split /\s+/, $_;
		#$total++ if $e[3] < $target and $e[3]+100>$target;
		$total++ if $e[3] <= $target and $e[3]+abs($e[8])>=$target;
		next if $e[1] == 4 or $e[5] !~ /[SHN]/;
		#next if $e[4] < 10;
		my $bp;
		my ($m1, $s1) = (0, 0);
		if ($e[5] =~ /N/) {
			my $sum = 0;
			my @f=split(/M/,$e[5]);
			for(my $i=0; $i<@f; $i++){
				$f[$i]="$f[$i]M";
				if($i==0){
                                       if ($f[$i] =~ /(\d+)[SHD](\d+)[M]/g) {
					     $bp = $e[3];
                                             if ($bp >= $target-1 and $bp <= $target+1) {
						if ($2 > 10){
                                             		$cnt ++;
				             		$qual += $e[4];
				             		if ($e[4] >= 29) {
				             			$cnt2 ++;
				             			$qual2 += $e[4];
				             		}
					     		if ($e[4] == 255){
					     			$cnt3 ++;
					     			$qual3 += $e[4];
                                                	}                     
                                              		last;
						}
                                              } else {
                                              	$bp = $e[3] + $2 - 1;
						if ($2 > 10){
                                             		if ($bp >= $target-1 and $bp <= $target+1) {
                                             			$cnt ++;
				             			$qual += $e[4];
				             			if ($e[4] >= 29) {
				             				$cnt2 ++;
				             				$qual2 += $e[4];
				             			}
					     			if ($e[4] == 255){
					     				$cnt3 ++;
					     				$qual3 += $e[4];
                                             		   	}                     
                                             		 	last;
                                             		 }
						}
                                              }
				       } elsif ($f[$i] =~ /(\d+)[M]/g) {
                                             $bp = $e[3] + $1 - 1;
					     if ($1 > 10){
				             	if ($bp >= $target-1 and $bp <= $target+1) {
				             			$cnt ++;
				             			$qual += $e[4];
				             			if ($e[4] >= 29) {
				             				$cnt2 ++;
				             				$qual2 += $e[4];
				             			}
								if ($e[4] == 255){
									$cnt3 ++;
									$qual3 += $e[4];
								}
				             	last; 
				             	}
                                             } 
				       }
                                }else{
					$sum = 0;
					while ($f[$i] =~ /(\d+)[SHDN]/g) {
                                        	$sum = $sum + $1;
					}
					$bp = $bp + $sum;
					if ($f[$i] =~ /(\d+)[M]/g) {
                                                if ($1 > 10){
							if ($bp >= $target-1 and $bp <= $target+1) {
								$cnt ++;
								$qual += $e[4];
								if ($e[4] >= 29) {
									$cnt2 ++;
									$qual2 += $e[4];
										
								}
								if ($e[4] == 255){
									$cnt3 ++;
									$qual3 += $e[4];
								}
				        		last; 
							}
							if ($i == scalar @f - 1) {
								last;
							}
							else{
								$bp = $bp + $1;
								if ($bp >= $target-1 and $bp <= $target+1) {
									$cnt ++;
									$qual += $e[4];
									if ($e[4] >= 29) {
										$cnt2 ++;
										$qual2 += $e[4];
											
									}
									if ($e[4] == 255){
										$cnt3 ++;
										$qual3 += $e[4];
									}
				        			last; 
                        					}
							}
						} else {
							$bp = $bp + $1;
						}
					}
                                }
			}	
		} else {
			while ($e[5] =~ /(\d+)[SH]/g) {
					$s1 = $1 if $1 > $s1;
			}
			while ($e[5] =~ /(\d+)M/g) {
					$m1 = $1 if $1 > $m1;
			}
			if ($e[5] =~ /$m1[M].*?$s1[SH]/) {
				$bp = $e[3]+$m1-1;
			} else {
				$bp = $e[3];
			}
                        if ($m1 > 10){
				if ($bp >= $target-1 and $bp <= $target+1) {
					$cnt ++;
					$qual += $e[4];
					if ($e[4] >= 29) {
						$cnt2 ++;
						$qual2 += $e[4];
							
					}
					if ($e[4] == 255){
						$cnt3 ++;
						$qual3 += $e[4];
					}
                        	}
                         }
		}	

	}
	$qual = $qual/$cnt if $cnt > 0;
	$qual2 = $qual2/$cnt2 if $cnt2 > 0;
	$qual3 = $qual3/$cnt3 if $cnt3 > 0;

	return ($total, $cnt, $qual, $cnt2, $qual2, $cnt3, $qual3);
}
