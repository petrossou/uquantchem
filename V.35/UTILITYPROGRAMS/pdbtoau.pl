#!/usr/bin/perl -w
# For converting the atomic coordinates in pdb-files
# to a coordinate file suitible for UQUANTCHEM

open(DATA1,"COORD.pdb");
open(DATA2,">COORD_IN_Au.dat");
$nlines = 0;
while(<DATA1>){
	if ( $_ =~ /ATOM/ ) {
		$nlines ++;
		@inlast = split(" ",$_);
		$r[$nlines][0] = $inlast[5]; $r[$nlines][1] = $inlast[6]; $r[$nlines][2] = $inlast[7];
		$ATOM1[$nlines] = substr($inlast[2], 0, 1 );
                $leng[$nlines] = length($inlast[2]);
		$ATOM2[$nlines] = "HH";
		if ( $leng[$nlines] >= 2 ) { $ATOM2[$nlines] = substr($inlast[2], 0, 2 );}else{$ATOM2[$nlines] = $ATOM1[$nlines];}
		$Z[$nlines] = 0;
		if ( $ATOM1[$nlines] =~ /H/ ) { $Z[$nlines] = 1; }
		if ( $ATOM1[$nlines] =~ /B/ ) { $Z[$nlines] = 5; }
		if ( $ATOM1[$nlines] =~ /C/ ) { $Z[$nlines] = 6; }
		if ( $ATOM1[$nlines] =~ /N/ ) { $Z[$nlines] = 7; }
		if ( $ATOM1[$nlines] =~ /O/ ) { $Z[$nlines] = 8; }
		if ( $ATOM1[$nlines] =~ /F/ ) { $Z[$nlines] = 9; }
		if ( $ATOM1[$nlines] =~ /P/ ) { $Z[$nlines] = 15; }
		if ( $ATOM1[$nlines] =~ /S/ ) { $Z[$nlines] = 16; }
		if ( $ATOM1[$nlines] =~ /K/ ) { $Z[$nlines] = 19; }
		if ( $ATOM1[$nlines] =~ /V/ ) { $Z[$nlines] = 23; }
		if ( $ATOM2[$nlines] =~ /He/ ) { $Z[$nlines] = 2; }
		if ( $ATOM2[$nlines] =~ /Li/ ) { $Z[$nlines] = 3; }
		if ( $ATOM2[$nlines] =~ /Be/ ) { $Z[$nlines] = 4; }
		if ( $ATOM2[$nlines] =~ /Ne/ ) { $Z[$nlines] = 10; }
		if ( $ATOM2[$nlines] =~ /Na/ ) { $Z[$nlines] = 11; }
		if ( $ATOM2[$nlines] =~ /Mg/ ) { $Z[$nlines] = 12; }
		if ( $ATOM2[$nlines] =~ /Al/ ) { $Z[$nlines] = 13; }
		if ( $ATOM2[$nlines] =~ /Si/ ) { $Z[$nlines] = 14; }
		if ( $ATOM2[$nlines] =~ /Cl/ ) { $Z[$nlines] = 17; }
		if ( $ATOM2[$nlines] =~ /Ar/ ) { $Z[$nlines] = 18; }
		if ( $ATOM2[$nlines] =~ /Ca/ ) { $Z[$nlines] = 20; }
		if ( $ATOM2[$nlines] =~ /Sc/ ) { $Z[$nlines] = 21; }
		if ( $ATOM2[$nlines] =~ /Ti/ ) { $Z[$nlines] = 22; }
		if ( $ATOM2[$nlines] =~ /Cr/ ) { $Z[$nlines] = 24; }
		if ( $ATOM2[$nlines] =~ /Mn/ ) { $Z[$nlines] = 25; }
		if ( $ATOM2[$nlines] =~ /Fe/ ) { $Z[$nlines] = 26; }
		if ( $ATOM2[$nlines] =~ /Co/ ) { $Z[$nlines] = 27; }
		if ( $ATOM2[$nlines] =~ /Ni/ ) { $Z[$nlines] = 28; }
		if ( $ATOM2[$nlines] =~ /Cu/ ) { $Z[$nlines] = 29; }
		if ( $ATOM2[$nlines] =~ /Zn/ ) { $Z[$nlines] = 30; }
		if ( $ATOM2[$nlines] =~ /Ga/ ) { $Z[$nlines] = 31; }
		if ( $ATOM2[$nlines] =~ /Ge/ ) { $Z[$nlines] = 32; }
		if ( $ATOM2[$nlines] =~ /As/ ) { $Z[$nlines] = 33; }
		if ( $ATOM2[$nlines] =~ /Se/ ) { $Z[$nlines] = 34; }
		if ( $ATOM2[$nlines] =~ /Br/ ) { $Z[$nlines] = 35; }
		if ( $ATOM2[$nlines] =~ /Kr/ ) { $Z[$nlines] = 36; }
	}
}
$nelect = 0;
$au = 0.52917720859;
$DR[0] = $r[1][0]; $DR[1] = $r[1][1]; $DR[2] = $r[1][2];
print "THE ATOMS ARE:\n";
for ( $i=1; $i <= $nlines; $i++ ) {
	 print "$ATOM2[$i] $Z[$i]\n";
	$rnew[0] = sprintf("%16f",($r[$i][0]-$DR[0])/$au);$rnew[1] = sprintf("%16f",($r[$i][1]-$DR[1])/$au);$rnew[2] = sprintf("%16f",($r[$i][2]-$DR[2])/$au);
	print DATA2 "ATOM $Z[$i] $rnew[0] $rnew[1] $rnew[2]\n";
	$nelect = $nelect + $Z[$i];
}
print "Total number of electrons in neutral molecule = $nelect\n";
close(DATA1);
close(DATA2);
