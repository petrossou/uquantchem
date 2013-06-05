#!/usr/bin/perl -w
use POSIX;
$header = "MG3S";

$FILENAMEF = $header.".dat";
$FILENAMERAW = $header.".raw";

open(FORTRANFILE,">$FILENAMEF");
open(RAWFILE,"$FILENAMERAW");

print FORTRANFILE "$header\n";
print FORTRANFILE "BASIS:\n";

$STARTREAD = 0;
$nstar = 0;
$Z = 0;

while(<RAWFILE>){
	if ( $_ =~ /\$basis/ ) { $STARTREAD = 1;}
	if ( $STARTREAD  == 1 ) {
		@read = split(" ",$_);
		$size = @read;
		if ( $_ =~ /\*/ && $size == 1 ) {$nstar++;}
		$_=~ s/\n//g;
		if ( $size == 2 ) {
			if ( $read[0] =~ /^\D/ && $read[0] !~ /\*/) { $Z++; $NAME[$Z] = ucfirst($read[0]); $nbas = 0;}
			if ( $read[0] =~/^\d/ && $read[1] =~ /\D/ && $read[1] !~ /\./) { $nbas++; $n = 0; $NPRIM[$nbas] = $read[0]; $L = $read[1]}
			if ( $read[0] =~/\./ && $read[1] =~ /\./ ){ $n++; $basfnk[$nbas][$n] = $_; $ANGULARMOM[$nbas][$n] = $L;};
		}
		$MOD = int($nstar-3*floor($nstar/3));
		if ( $nstar > 2 && $MOD == 0 ) {
			$nstar = 1;
			#print "$MOD\n";
			print FORTRANFILE "$NAME[$Z] $Z $nbas";
			for ( $i = 1; $i < $nbas; $i++ ) { print FORTRANFILE " $NPRIM[$i]";}
			print FORTRANFILE " $NPRIM[$nbas]\n";
			print FORTRANFILE "\n";
			for ( $i = 1; $i <= $nbas; $i++ ) {
				for ( $j = 1; $j <= $NPRIM[$i]; $j++ ) { print FORTRANFILE "$ANGULARMOM[$i][$j] $basfnk[$i][$j]\n";}
				print FORTRANFILE "\n";
			}
		}

	}
}
close(FORTRANFILE);
close(RAWFILE);

