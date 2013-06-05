#!/usr/bin/perl -w
#===============================================================#
# This script takes the uquantchem movie-file: MOLDYNMOVIE.xsf 	#
# and calculates the velocity aouto correlation function. The	#
# result is saved in the file: VELAUTO.dat.			#
#===============================================================#
# User specified input:
$NATOM = 3;				# Total number of atoms in molecule.
$DT = 5*0.02418884326505/1000; 		# Time-step in ps.
$START = 100;				# Defining the t=0, starting time.
$MAXIMUMFREQ = 150;			# Maximum frequency in THz
$FREQMESH = 10000;			# The frequency mesh.

open(DATA1,"MOLDYNMOVIE.xsf");
open(DATA2,">VELAUTO.dat");
$nlines = 0;
$norm = 0;
$count = 0;
while(<DATA1>){
	@inlast = split(" ",$_);
	if ( $_ =~ /ATOMS/ ) {$nlines = 0; $time = $inlast[1];$VCORR =0;}
	$length = @inlast;
	if ( $length == 7 ) {
		$nlines++;
		if ( $time == $START ) { $V0[1][$nlines] = $inlast[4]; $V0[2][$nlines] = $inlast[5]; $V0[3][$nlines] = $inlast[6]; }
		if ( $time >= $START ) {           $V[1] = $inlast[4];           $V[2] = $inlast[5];           $V[3] = $inlast[6]; }
		if ( $time >= $START ) { 
			$VCORR = $VCORR + $V0[1][$nlines]*$V[1] + $V0[2][$nlines]*$V[2]  + $V0[3][$nlines]*$V[3];
			if ( $time == $START ) { $norm = $norm + $V0[1][$nlines]*$V0[1][$nlines] + $V0[2][$nlines]*$V0[2][$nlines]  + $V0[3][$nlines]*$V0[3][$nlines];}
			$t = ($time-$START)*$DT;
			if ( $nlines == $NATOM  ) {$VCORR = $VCORR/$norm; print DATA2 "$t $VCORR\n"; $count++; $tid[$count] = $t; $VC[$count] = $VCORR;}
		}
	}
}
$minf = (1.0/$t);
$minfcm = $minf*33.35640946151377;

$maxf = (1.0/$DT)/10;
$maxfcm = $maxf*33.35640946151377;

print "\n";
print " ===============================================================================================\n";
print "                  Minimum frequency to be resolved by this calculation\n";
print "                  f_min = $minf [THz] = $minfcm [cm-1]\n";
print "\n";
print "  Maximum frequency to be resolved by this calculation (Asuming that at least 10 time-steps\n";
print "  are needed to resolve each period), f_max = $maxf [THz] = $maxfcm [cm-1]\n";
print " ===============================================================================================\n";
print "\n";

$pi = 3.14159265358979323846264338327950288;

#############################################################################################
# Here we calculate the temporal Fourier transform for the Velocity Auto Correlation Function
#############################################################################################
open(DATA3,">SPECTRUM_THz.dat");
open(DATA4,">SPECTRUM_cm.dat");
for ( $i= 0; $i <= $FREQMESH; $i++ ) {
	$freq = $minf + ($MAXIMUMFREQ-$minf)*$i/$FREQMESH;
	$FTRANS = 0;
	for ( $j= 1; $j <= $count-1; $j++ ) {
		$FTRANS = $FTRANS + 0.5*$DT*($VC[$j]*cos(2.0*$pi*$freq*$tid[$j])+$VC[$j+1]*cos(2.0*$pi*$freq*$tid[$j+1]))
	}
	$FTRANS = ($FTRANS/$tid[$count])**2;
	$frequencyTHz = $freq;
	$frequencycm  = $frequencyTHz*33.35640946151377;
	print DATA3 "$frequencyTHz $FTRANS\n";
	print DATA4 "$frequencycm $FTRANS\n";
}
close(DATA1);
close(DATA2);
close(DATA3);
close(DATA4);
