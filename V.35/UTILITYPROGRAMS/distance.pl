#!/usr/bin/perl -w
# This script takes the uquantchem movie-file: MOLDYNMOVIE.xsf 
# and calculates the interatomic distance between atom A and B
# as a function of time.
#
# User specified input:
$NATOM = 3;	# Total number of atoms in molecule.
$A = 1;   	# Index of atom A (not atomic number)
$B = 3;		# Index of atom B (not atomic number)

open(DATA1,"MOLDYNMOVIE.xsf");
open(DATA2,">RAB.dat");
$nlines = 0;
while(<DATA1>){
	@inlast = split(" ",$_);
	if ( $_ =~ /ATOMS/ ) {$nlines = 0; $time = $inlast[1];}
	$length = @inlast;
	if ( $length == 7 ) {
		$nlines++;
		if ( $nlines == $A ) { $RA[1] = $inlast[1]; $RA[2] = $inlast[2]; $RA[3] = $inlast[3]; }
		if ( $nlines == $B ) { $RB[1] = $inlast[1]; $RB[2] = $inlast[2]; $RB[3] = $inlast[3]; }
		if ( $nlines == $NATOM ) { 
			$RAB = sqrt( ($RA[1]-$RB[1])**2 + ($RA[2]-$RB[2])**2 + ($RA[3]-$RB[3])**2 );
			print DATA2 "$time $RAB\n";
		}
	}
}
