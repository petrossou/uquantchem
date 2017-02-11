#!/usr/bin/perl -w
# For converting the atomic coordinates in pdb-files
# to a coordinate file suitible for UQUANTCHEM

open(DATA1,"COORD.dat");
open(DATA2,">COORD_IN_Ang.dat");
while(<DATA1>){
		@koord = split(" ",$_);
		$koord[2] = $koord[2]*0.52917720859;
		$koord[3] = $koord[3]*0.52917720859;
		$koord[4] = $koord[4]*0.52917720859;
		print DATA2 "$koord[1] $koord[2] $koord[3] $koord[4]\n";
	}
close(DATA2);
close(DATA1);
