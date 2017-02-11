#!/usr/bin/perl -w
######################################################################
# This script only takes every m:th dqmc energy, EL, and calculates  #
# the the mean <EL> and standard deviation. It also calculates the   #
# correlation between two consecutive walker populations             #
######################################################################

#################################################################
# 			USER SPECIFIED INPUT		        #
#################################################################

$m = 5; # We use only every m:th Local energy to avoid 
	# correlation between walker populations.

$ns = 600000; # Start calculating the standard deviation 
	      # from this time step

#################################################################
#    	       	END OF USER SPECIFIED INPUT  	                #
#################################################################
print "\n";
print "##########################################################\n";
print "# Performing statistical analysis on the Local energy of #\n";
print "#  the Diffision Quantum Monte Carlo (DQMC) Energy EL.   #\n";
print "##########################################################\n";

open(DATA1,"ENERGYDQMC.dat");
open(DATA2,">ENERGYDQMC_CHECKED.dat");
$n = 1;
$i = 0;
$EM[0]  = 0;
$EM2 = 0;

while(<DATA1>){
	$t = $n/$m;
	$nm = sprintf("%.0f",$t );
	if ( $n/$m - $nm == 0 ) {
		@E = split(" ",$_);
		$i++;
		$EM[$i] = ($EM[$i-1]*($i-1) + $E[3])/$i;
		$EM2 = ($EM2*($i-1) + $E[3]*$E[3])/$i;
		$EL[$i] = $E[3];
		print DATA2 "$i $EM[$i] $E[2] $E[3]\n";
	}
	$n++;
}
close(DATA1);
close(DATA2);

################################################################
#	CALCULATION OF MEAN LOCAL ENERGY AND THE CORRELATION   #
################################################################

$EE = 0;
$sigma = 0;

for ( $k=1; $k <= $i; $k++ ) {
	if ( $k < $i ) {$EE = $EE + $EL[$k]*$EL[$k+1];}
	if ( $k > $ns ){$sigma = $sigma + ($EM[$k]-$EM[$i])**2;}
}

$sigma = ($sigma/($i-$ns))**0.5;
$EE  = $EE/($i-1);

$EM[$i] = sprintf("%.8f",$EM[$i] );
$sigma  = sprintf("%.8f",$sigma );
print "\n";
print "------------------------------------------------------------\n";
print "| When sampling the local energy of every $m:th population  |\n";
print "|     the following statistical result is obtained:        |\n";
print "------------------------------------------------------------\n";
print "\n";
print "  -----------------------------------------------------\n";
print "  |  DQMC Energy: E = $EM[$i] +/- $sigma au  |\n";
print "  -----------------------------------------------------\n";
print "\n";
$CORR = ($EE - $EM[$i]**2)/($EM2-$EM[$i]**2);
$CORR  = sprintf("%.4f",$CORR );
print "  -----------------------------------------------------\n";
print "  |  Correlation between walker populations = $CORR  |\n";
print "  -----------------------------------------------------\n";
print "\n";
