#!/usr/bin/perl -w
#
use strict;
use warnings;

while (<>) {
	my @e = split;
	next if $e[14] eq 'mark';
        if ($e[16] >= 2 and $e[23] >= 2 and $e[3] >= 30 and $e[15]>=3 and $e[16]/$e[15]>0.1 and ($e[22]>=3 and $e[23]/$e[22]>0.1)) {
		#print if $e[19]>0 and $e[20]/$e[19]>0.1 and ($e[26]>0 and $e[27]/$e[26]>0.1);
		print;
	} elsif ($e[20] >= 1 and $e[27] >= 1 and $e[3] >= 30 and $e[8] >= 3 and $e[15]>=3 and $e[16]/$e[15]>0.01 and ($e[22]>=3 and $e[23]/$e[22]>0.01)) {
		print;
	} elsif ($e[8] >= 5 and $e[3] >= 60){
		print;
	}
}
