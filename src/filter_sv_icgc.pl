#!/usr/bin/perl -w
#
use strict;
use warnings;

while (<>) {
	my @e = split;
	next if $e[14] eq 'mark';
        if ($e[15] >= 2 and $e[22] >= 2 and $e[3] >= 30 and $e[14]>=3 and $e[15]/$e[14]>0.1 and ($e[21]>=3 and $e[22]/$e[21]>0.1)) {
		print;
	} elsif ($e[19] >= 1 and $e[26] >= 1 and $e[3] >= 30 and $e[8] >= 3 and $e[14]>=3 and $e[15]/$e[14]>0.003 and ($e[21]>=3 and $e[22]/$e[21]>0.003)) {
		print;
	#} elsif ($e[14] == 2 and $e[3] >= 60) {
	#	print;
	} elsif ($e[8] >= 5 and $e[3] >= 60){
		print;
	}
}
