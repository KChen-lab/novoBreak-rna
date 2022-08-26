#!/usr/bin/perl -w
#
#Author: Yukun Tan
#
use strict;
use warnings;

my @slines = ();
my @dlines = ();

my $id = 0;
print qq(#CHROM	POS	ALT	QUAL	INFO	ENSEMBL
);
my $pre = "";
my @results = ();
while (<>) {
	next if /^@/;
	my @e = split /\s+/, $_;
	next if $e[1] == 4 or $e[5]!~/[SH]/;
	my $newcigar = "";
	my $match = 0;
	while ($e[5] =~ /(\d+)([SMIDH])/g) {
		my ($l, $t) = ($1, $2);
		if ($t eq 'S' or $t eq 'H') {
			if (!$match) {
				$newcigar .= $l.$t;
			} else {
				$newcigar .= $match."M";
				$newcigar .= $l.$t;
				$match = 0;
			}
		}
		$match += $l if $t eq 'M';
		if ($t eq 'I' or $t eq 'D') {
			if ($l > 20) { # TODO: magic number
				$newcigar .= $match."M";
				$newcigar .= $l.$t;
				$match = 0;
			} else {
				$match += $l if $t eq 'I';
			}
		}
	}
	$newcigar .= $match."M" if $match;
	$e[5] = $newcigar;
	#if ($e[5] =~ /(\d+)[SH].+M(\d+)[SH]/) {
	#	next if ($1>5 and $2>5);
	#}
	if ($e[5] =~ /(\d+)M.+?(\d+)M/) {
		next if ($1>10 and $2>10);
	}
	$_ = join("\t", @e);
	if ($e[0] ne $pre) {
		if ($pre ne "") {
			if (@dlines >= 2) {
				my $seq = &pick_consensus(@dlines);
				for (my $i = 0; $i < scalar @dlines - 1; $i++) {
					for (my $j = $i+1; $j < scalar @dlines; $j++) {
						push @results, &parse_bp1($dlines[$i], $dlines[$j], $seq);
					}
				}
			}
			@dlines = ();
		} 
		$pre = $e[0];
		if ($e[5] =~ /(\d+)[SH].+M(\d+)[SH]/) {
			my $tmp = $e[5];
			if ($1 > 20) {
				$e[5] =~ s/$2[SH]//;
				$_ = join ("\t", @e);
				push @dlines, $_;
				$e[5] = $tmp;
			}
			if ($2 > 20) {
				$e[5] =~ s/$1[SH]//;
				$_ = join ("\t", @e);
				push @dlines, $_;
			}
			$e[5] = $tmp;
		} else {
			push @dlines, $_;
		}
	} else {
		if ($e[5] =~ /(\d+)[SH].+M(\d+)[SH]/) {
			my $tmp = $e[5];
			if ($1 > 20) {
				$e[5] =~ s/$2[SH]//;
				$_ = join ("\t", @e);
				push @dlines, $_;
				$e[5] = $tmp;
			}
			if ($2 > 20) {
				$e[5] =~ s/$1[SH]//;
				$_ = join ("\t", @e);
				push @dlines, $_;
			}
			$e[5] = $tmp;
		} else {
			push @dlines, $_;
		}
	}
}

if (@dlines >= 2) {
	my $seq = &pick_consensus(@dlines);
	for (my $i = 0; $i < scalar @dlines - 1; $i++) {
		for (my $j = $i+1; $j < scalar @dlines; $j++) {
			push @results, &parse_bp1($dlines[$i], $dlines[$j], $seq);
		}
	}
}

@results = sort {$a->[0] cmp $b->[0] or $a->[1] <=> $b->[1]} @results;
my $prepos = 0;
for (my $i = 0; $i < @results; $i++) {
	print join("\t", @{$results[$i]}), "\n";
}

1;

sub pick_consensus {
	my @lines = @_;
	foreach my $l (@lines) {
		my @e = split /\s+/, $l;
		next if $e[1] & 256;
		return $e[9];
	}
}

sub parse_bp1 { 
	my $line1 = shift;
	my $line2 = shift;
	my $seq = shift;
	my @ret = ();
	my ($m1, $s1, $m2, $s2) = (0, 0, 0, 0);
	my @e1 = split /\s+/, $line1;
	my @e2 = split /\s+/, $line2;
	if ($e1[5] =~ /[SMH]/ and $e2[5] =~ /[SMH]/) {
			while ($e1[5] =~ /(\d+)[SH]/g) {
					$s1 = $1 if $1 > $s1;
			}
			while ($e1[5] =~ /(\d+)M/g) {
					$m1 = $1 if $1 > $m1;
			}

			while ($e2[5] =~ /(\d+)[SH]/g) {
					$s2 = $1 if $1 > $s2;
			}
			while ($e2[5] =~ /(\d+)M/g) {
					$m2 = $1 if $1 > $m2;
			}
	}
	my ($pos1, $pos2) = (0, 0);
	if ($e1[5] =~ /$s1[SH].*?$m1[M]/) {
		$pos1 = $e1[3];
	} else {
		$pos1 = $e1[3]+$m1-1;
	}
	if ($e2[5] =~ /$s2[SH].*?$m2[M]/) {
		$pos2 = $e2[3];
	} else {
		$pos2 = $e2[3]+$m2-1;
	}
	if (abs($m1-$s2)<=25 or abs($m2-$s1)<=25) {
		if ($e1[2] ne $e2[2]) { # trans
			if ((($e1[1] & 0x10) ^ ($e2[1] & 0x10)) == 0) { # +,+ or -,-
				if (($pos1 == $e1[3] and $pos2 != $e2[3]) or ($pos2 == $e2[3] and $pos1 != $e1[3])) {
					if ($pos1 < $pos2) {
						push @ret, [ ($e1[2], $pos1, "<TRA>", ($e1[4]+$e2[4])/2, "CONSENSUS=$seq;SVTYPE=TRA;CHR2=$e2[2];END=$pos2;SVLEN=0", $e1[0]) ];
					} else {
						push @ret, [ ($e2[2], $pos2, "<TRA>", ($e1[4]+$e2[4])/2, "CONSENSUS=$seq;SVTYPE=TRA;CHR2=$e1[2];END=$pos1;SVLEN=0", $e2[0]) ];
					}	
				}
			} else { # +,- or -,+
				if (($pos1 == $e1[3] and $pos2 == $e2[3]) or ($pos2 != $e2[3] and $pos1 != $e1[3])) {
					if ($pos1 < $pos2) {
						push @ret, [ ($e1[2], $pos1, "<TRA>", ($e1[4]+$e2[4])/2, "CONSENSUS=$seq;SVTYPE=TRA;CHR2=$e2[2];END=$pos2;SVLEN=0", $e1[0]) ];
					} else {
						push @ret, [ ($e2[2], $pos2, "<TRA>", ($e1[4]+$e2[4])/2, "CONSENSUS=$seq;SVTYPE=TRA;CHR2=$e1[2];END=$pos1;SVLEN=0", $e2[0]) ];
					}
				}
			}
		} else {
			if ((($e1[1] & 16) and !($e2[1] & 16)) or ((!($e1[1] & 16)) and ($e2[1] & 16))) { # inv
				if(abs($pos2 - $pos1) < 10) { # TODO magic number controling inv size
				} 
				else {
					if (($pos1 == $e1[3] and $pos2 == $e2[3]) or ($pos2 != $e2[3] and $pos1 != $e1[3])) {
						if ($pos1 < $pos2) {
							push @ret, [ ($e1[2], $pos1, "<INV>", ($e1[4]+$e2[4])/2, "CONSENSUS=$seq;SVTYPE=INV;CHR2=$e2[2];END=$pos2;SVLEN=".abs($pos2-$pos1+1), $e2[0]) ] ;
						} else {
							push @ret, [ ($e2[2], $pos2, "<INV>", ($e1[4]+$e2[4])/2, "CONSENSUS=$seq;SVTYPE=INV;CHR2=$e2[2];END=$pos1;SVLEN=".abs($pos2-$pos1+1), $e2[0]) ] ;
						}
					}
				}
			} else {
				if ($pos1 < $pos2) {
					if ((!($e1[1]&16) and (($pos1 == $e1[3] and $pos2 != $e2[3]) or ($pos1 != $e1[3] and $pos2 == $e2[3]))) or ($e1[1]&16 and (($pos1 != $e1[3] and $pos2 == $e2[3]) or ($pos1 == $e1[3] and $pos2 != $e2[3])))) { 
						if(abs($pos2 - $pos1) > 10) { # TODO magic number controling del size
							push @ret, [ ($e1[2], $pos1, "<DEL/DUP>", ($e1[4]+$e2[4])/2, "CONSENSUS=$seq;SVTYPE=DEL;CHR2=$e2[2];END=$pos2;SVLEN=".($pos1-$pos2), $e2[0])] ;
					}} #else { # dup
				} else {
					if ((!($e2[1]&16) and (($pos1 == $e1[3] and $pos2 != $e2[3]) or ($pos1 != $e1[3] and $pos2 == $e2[3]))) or ($e2[1]&16 and (($pos1 != $e1[3] and $pos2 == $e2[3]) or ($pos1 == $e1[3] and $pos2 != $e2[3])))) {
						if(abs($pos2 - $pos1) > 10) { # TODO magic number controling del size
							push @ret, [ ($e2[2], $pos2, "<DEL/DUP>", ($e1[4]+$e2[4])/2, "CONSENSUS=$seq;SVTYPE=DEL;CHR2=$e1[2];END=$pos1;SVLEN=".($pos2-$pos1), $e2[0]) ];
						}
					} 
				}
			}
		}
	}
	return @ret;
}
