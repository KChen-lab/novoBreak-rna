#!/usr/bin/perl -w
#
use strict;
use warnings;

while (<>) {
	my @e = split;
	print;
	#my $len = abs $1 if /SVLEN=(\S+)/;
	#if ($len < 1000 and $len > 100 and $e[22] + $e[32] > 4 and $e[25] <= 1 and $e[35] <= 1 and $e[13] > 3 and $e[5] > 10) { # small svs
	#	print if $e[19]>0 and $e[20]/$e[19]>0.1 or ($e[29]>0 and $e[30]/$e[29]>0.1);
	#	#} elsif ($e[35]>=2 and $e[36]/$e[35] < 0.1 and $e[37]>=2 and $e[38]/$e[37] < 0.1 and $e[36]+$e[38]<1  and $e[15] < 100 and $e[25] < 100 or ($e[35]<=1 and $e[37]<=1 and $e[13]>4 and $e[21]<1 and $e[31]<1 and $e[18]>=2 and $e[28]>=2)) { # this may be a better filter for real data
	#} elsif ($e[39]>=3 and $e[40]/$e[39] < 0.1 and $e[41]>=3 and $e[42]/$e[41] < 0.1 and $e[40] <=2 and $e[42]<=2 or ($e[39]<=1 and $e[41]<=1 and $e[13]>4 and $e[25]<1 and $e[35]<1 and $e[22]>=2 and $e[32]>=2)) { # this may be a better filter for real data
#	#elsif ($e[36] <= 1 and $e[38] <= 1 and $e[23] <= 1 and $e[33] <= 1) {
	#	if (/INV/) {
	#		#print if $len < 20000000;
	#		print;
	#	} elsif (/TRA/) {
	#		print if $e[39] > 3 or $e[41] > 3;
	#	} else {
	#		print;
	#	}
	#} elsif (/TRA/) {
	#	print if ($e[5] > 30 and $e[22]+$e[32] > 10 and $e[25]<=1 and $e[35]<=1);
	#}
}
