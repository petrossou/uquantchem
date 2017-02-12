#!/usr/bin/perl -w
#================================================================#
# This script takes absorpion spectra: ABSSPECTRUM_X.dat,  	 #
# ABSSPECTRUM_Y.dat and ABSSPECTRUM_Z.dat, produced from 3 dirac #
#    pulse TDDFT/THF calculations and merges them into one	 #
#       absorption spectra result is saved in the file:		 #
#                 MERGED_ABS_SPECTRA.dat			 #
#================================================================#
$n = 0;
open(DATAX,"ABSSPECTRUM_X.dat");
while(<DATAX>){
	$n++;
	@inlast = split(" ",$_);
	$OMEGA[$n] = $inlast[0]; $SX[$n] = $inlast[1]; $SSX[$n] = $inlast[2];
}
close(DATAX);

$n = 0;
open(DATAY,"ABSSPECTRUM_Y.dat");
while(<DATAY>){
	$n++;
	@inlast = split(" ",$_);
	$OMEGA[$n] = $inlast[0]; $SY[$n] = $inlast[1]; $SSY[$n] = $inlast[2];
}
close(DATAY);

$n = 0;
open(DATAZ,"ABSSPECTRUM_Z.dat");
while(<DATAZ>){
	$n++;
	@inlast = split(" ",$_);
	$OMEGA[$n] = $inlast[0]; $SZ[$n] = $inlast[1]; $SSZ[$n] = $inlast[2];
}
close(DATAZ);

open(DATA,">MERGED_ABS_SPECTRA.dat");
for ( $i=1;$i <= $n; $i++ ) {
	$OMEGA[$i] = 27.2113838*$OMEGA[$i];
	$S1 = $SX[$i] + $SY[$i] + $SZ[$i];
	$S2 = $SSX[$i] + $SSY[$i] + $SSZ[$i];
	print DATA "$OMEGA[$i] $S1 $S2\n";
}
close(DATA);
