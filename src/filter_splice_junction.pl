#!/usr/bin/perl -w
#
use strict;
use warnings;

while (<>) {
	my @e = split;
	#next if $e[6] eq 8;
        if ($e[21] >= 2 and $e[28] >= 2 and $e[20]>= 3 and $e[21]/$e[20]>0.1 and ($e[27]>= 5 and $e[28]/$e[27]>0.1)) {
		#print if $e[19]>0 and $e[20]/$e[19]>0.1 and ($e[26]>0 and $e[27]/$e[26]>0.1);
		print;
	} elsif ($e[25] >= 1 and $e[32] >= 1 and $e[20]>= 3 and $e[21]/$e[20]>0.001 and ($e[27]>= 3 and $e[28]/$e[27]>0.001)) {
		print;
	} elsif ($e[18] >= 5 and $e[17] >= 60){
		print;
	}
}
