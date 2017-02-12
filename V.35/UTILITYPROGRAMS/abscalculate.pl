#!/usr/bin/perl -w
$DT = 0.01;
$EMAX = 1.0/100000000;
$OMEGAMAX = 100;

$nlines = -1;
$MEAN = 0;
open(DATA1,"TDFTOUT.dat");
open(DATA2,">ABSSPECTRUMM.dat");
while(<DATA1>){
	@inlast = split(" ",$_);
	if ( $_ !~ /E-field/ ) {$nlines++; $time[$nlines] = $inlast[1]; $MU[$nlines] = $inlast[2];$MEAN = $MEAN + $inlast[2];}
}
$MEAN = $MEAN/$nlines;
$pi = 3.14159265358979323846264338327950288;
$Domega = 2*$pi/($DT*$nlines);
close(DATA1);

$NOMEGA = sprintf("%.0f", $nlines/2);
for ( $i= 0; $i <= $NOMEGA; $i++ ) {
	$F[$i] = 0;
	$OMEGA[$i] = $Domega*$i;
	for ( $j= 0; $j<= $nlines-1; $j++ ) {
		$F[$i] = $F[$i] + 2*$DT*sin($OMEGA[$i]*$time[$j])*($MU[$j]-$MEAN)/$EMAX;
	}
	$F[$i] = abs($F[$i]);
	print DATA2 "$OMEGA[$i] $F[$i]\n";
	if ( $OMEGA[$i] > $OMEGAMAX ) { $i = $NOMEGA+1;}
}
close(DATA2);
