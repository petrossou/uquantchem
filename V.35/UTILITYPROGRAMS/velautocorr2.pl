#!/usr/bin/perl -w
#===============================================================#
# This script takes the uquantchem movie-file: MOLDYNMOVIE.xsf 	#
# and calculates the velocity aouto correlation function. The	#
# result is saved in the file: VELAUTO.dat.			#
#===============================================================#
# User specified input:
$NATOM = 4;				# Total number of atoms in molecule.
$DT = 5*0.02418884326505/1000; 		# Time-step in ps.
$MAXIMUMFREQ = 150;			# Maximum frequency in THz
$FREQMESH = 10000;			# The frequency mesh.
$NSAMPLES = 2000;			# Number of starting times for sampling of the velocity auto correlation function
$BEGIN = 1000;				# From this timestep and onward the  velocity auto correlation function will be sampled.
$NTIMESTEPS = 20000;			# Number of timesteps used to calculate the velocity auto correlation function.


$nlines = 0;
$count = 0;
$available = 0;

# counting the maximum number of time-steps:
open(DATA1,"MOLDYNMOVIE.xsf");
while(<DATA1>){
        @inlast = split(" ",$_);
	if ( $_ =~ /ATOMS/ ) {$available++;}
	if ( $_ =~ /ATOMS/ ) {$nlines = 0; $time = $inlast[1];}
	$length = @inlast;
	if ( $length == 7 ) {
		$nlines++;
		$VELOCITY_X[$time][$nlines] = $inlast[4];
		$VELOCITY_Y[$time][$nlines] = $inlast[5];
		$VELOCITY_Z[$time][$nlines] = $inlast[6];
	}
}
close(DATA1);

# Checking that everything is kosher with the input:
if ( $NTIMESTEPS + $BEGIN + $NSAMPLES > $available ) {
	print "\n";
	print "=====================================================================\n";
	print "     WARNING THE SPECIFIED NUMBER OF TIMESTEPS = $NTIMESTEPS, IS BIGGER\n";
	print "            THAN THE AVAILABLE NUMBER OF TIMESTEPS = $available\n";
	$NTIMESTEPS = $available - $BEGIN - $NSAMPLES -1;
	print "                    DECREASING NTIMESTEPS TO $NTIMESTEPS\n";
	print "=====================================================================\n";
	print "\n";
}


for ( $START = $BEGIN; $START <= $BEGIN + $NSAMPLES; $START++ ) {
	$NO = $START+1 - $BEGIN;
	print "t0-SAMPLE No: $NO\n";
	$count = 0;
	for ( $j=1; $j <= $NATOM; $j++ ) { $V0[$j][1] = $VELOCITY_X[$START][$j]; $V0[$j][2] = $VELOCITY_Y[$START][$j]; $V0[$j][3] = $VELOCITY_Z[$START][$j];}
	for ( $j=1; $j <= $NATOM; $j++ ) { $norm[$j] = $V0[$j][1]**2 + $V0[$j][2]**2 + $V0[$j][3]**2;}
	for ( $i =$START; $i <= $START+$NTIMESTEPS; $i++ ) {
		$count++;
		$tid[$count] = $count*$DT; 
		if ( $START == $BEGIN ){ $VC[$count] = 0;}
		for ( $j=1; $j <= $NATOM; $j ++ ) { 
			$VC[$count] = $VC[$count] + ($V0[$j][1]*$VELOCITY_X[$i][$j] + $V0[$j][2]*$VELOCITY_Y[$i][$j] + $V0[$j][3]*$VELOCITY_Z[$i][$j])/$norm[$j];
		}
		$VC[$count] = $VC[$count]/$NATOM;
	}
}

open(DATA2,">VELAUTO.dat");
for ( $i=1; $i <= $count; $i++){print DATA2 "$tid[$i]  $VC[$i]\n";}

$minf = (1.0/$tid[$count]);
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
close(DATA2);
close(DATA3);
close(DATA4);
