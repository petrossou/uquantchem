#!/usr/bin/perl -w
# This script when executed runs all the test calculations
# in the directories TEST/RUN1, TEST/RUN2, ..., TEST/RUN7
$x = 0;
if (-e "SERIALVERSION/uquantchem.s" ) {
        print "\n";
	print "***********************************************\n";
	print "Performing test on serial version of UQUANTCHEM\n";
	print "***********************************************\n";
	print "\n";
 	$dir = `pwd`;
	chomp $dir;
	for ( $i =1; $i <= 9; $i++ ) {
		chdir("$dir/TESTS/RUN$i/");
		if ( $i == 1 ) { print "Testing the Hartree-Fock implementation"}
		if ( $i == 2 ) { print "Testing the DFT-PBE implementation"}
		if ( $i == 3 ) { print "Testing the DFT-B3LYP implementation"}
		if ( $i == 4 ) { print "Testing the MP2 implementation"}
		if ( $i == 5 ) { print "Testing the CISD implementation"}
		if ( $i == 6 ) { print "Testing the QDMC implementation"}
		if ( $i == 7 ) { print "Testing the relaxation implementation"}
		if ( $i == 8 ) { print "Testing the tddft implementation"}
		if ( $i == 9 ) { print "Testing the molecular dynamics implementation"}
     
		if ( -e "OUT.dat" ) {system("rm OUT.dat");}
		system("$dir/SERIALVERSION/uquantchem.s > OUT.dat");

	       	open(DATA,"OUT.dat");
		$koll = 0;
		if ( $i == 1 ) { while(<DATA>){ if ( $_ =~ /Hartree-Fock energy:/ ) { $koll = 1; } } if ( $koll == 1 ) {print " ok! \n";}else{ print "  FAILED !\n";} }
		if ( $i == 2 ) { while(<DATA>){ if ( $_ =~ /PBE energy:/ ) { $koll = 1; } } if ( $koll == 1 ) {print " ok! \n";}else{ print "  FAILED !\n";} }
		if ( $i == 3 ) { while(<DATA>){ if ( $_ =~ /CALCULATED INTERATOMC FORCES:/ ) { $koll = 1; } } if ( $koll == 1 ) {print " ok! \n";}else{ print "  FAILED !\n";} }
		if ( $i == 4 ) { while(<DATA>){ if ( $_ =~ /Results from the MP2 calculation:/ ) { $koll = 1; } } if ( $koll == 1 ) {print " ok! \n";}else{ print "  FAILED !\n";} }
		if ( $i == 5 ) { while(<DATA>){ if ( $_ =~ /RUNTIME=/ ) { $koll = 1; } } if ( $koll == 1 ) {print " ok! \n";}else{ print "  FAILED !\n";} }
		if ( $i == 6 ) { while(<DATA>){ if ( $_ =~ /RUNTIME=/ ) { $koll = 1; } } if ( $koll == 1 ) {print " ok! \n";}else{ print "  FAILED !\n";} }
		if ( $i == 7 ) { while(<DATA>){ if ( $_ =~ /The relaxed positions of the nuclea are the following:/ ) { $koll = 1; } } if ( $koll == 1 ) {print " ok! \n";}else{ print "  FAILED !\n";} }
		if ( $i == 8 ) { while(<DATA>){ if ( $_ =~ /RUNTIME=/ ) { $koll = 1; } } if ( $koll == 1 ) {print " ok! \n";}else{ print "  FAILED !\n";} }
		if ( $i == 9 ) { while(<DATA>){ if ( $_ =~ /RUNTIME=/ ) { $koll = 1; } } if ( $koll == 1 ) {print " ok! \n";}else{ print "  FAILED !\n";} }
		close(DATA);
	}
	chdir("$dir");
}

if (-e "OPENMPVERSION/uquantchem.omp" ) {
        print "\n";
	print "********************************************\n";
	print "Performing test on omp version of UQUANTCHEM\n";
	print "********************************************\n";
	print "\n";
 	$dir = `pwd`;
	chomp $dir;
	for ( $i =1; $i <= 8; $i++ ) {
		chdir("$dir/TESTS/RUN$i/");
		if ( $i == 1 ) { print "Testing the Hartree-Fock implementation"}
		if ( $i == 2 ) { print "Testing the DFT-PBE implementation"}
		if ( $i == 3 ) { print "Testing the DFT-B3LYP implementation"}
		if ( $i == 4 ) { print "Testing the MP2 implementation"}
		if ( $i == 5 ) { print "Testing the CISD implementation"}
		if ( $i == 6 ) { print "Testing the QDMC implementation"}
		if ( $i == 7 ) { print "Testing the relaxation implementation"}
		if ( $i == 8 ) { print "Testing the TDHF implementation"}
		
		if ( -e "OUT.dat" ) {system("rm OUT.dat");}
		system("$dir/OPENMPVERSION/uquantchem.omp > OUT.dat");

	       	open(DATA,"OUT.dat");
		$koll = 0;
		if ( $i == 1 ) { while(<DATA>){ if ( $_ =~ /Hartree-Fock energy:/ ) { $koll = 1; } } if ( $koll == 1 ) {print " ok! \n";}else{ print "  FAILED !\n";} }
		if ( $i == 2 ) { while(<DATA>){ if ( $_ =~ /PBE energy:/ ) { $koll = 1; } } if ( $koll == 1 ) {print " ok! \n";}else{ print "  FAILED !\n";} }
		if ( $i == 3 ) { while(<DATA>){ if ( $_ =~ /CALCULATED INTERATOMC FORCES:/ ) { $koll = 1; } } if ( $koll == 1 ) {print " ok! \n";}else{ print "  FAILED !\n";} }
		if ( $i == 4 ) { while(<DATA>){ if ( $_ =~ /Results from the MP2 calculation:/ ) { $koll = 1; } } if ( $koll == 1 ) {print " ok! \n";}else{ print "  FAILED !\n";} }
		if ( $i == 5 ) { while(<DATA>){ if ( $_ =~ /RUNTIME=/ ) { $koll = 1; } } if ( $koll == 1 ) {print " ok! \n";}else{ print "  FAILED !\n";} }
		if ( $i == 6 ) { while(<DATA>){ if ( $_ =~ /RUNTIME=/ ) { $koll = 1; } } if ( $koll == 1 ) {print " ok! \n";}else{ print "  FAILED !\n";} }
		if ( $i == 7 ) { while(<DATA>){ if ( $_ =~ /The relaxed positions of the nuclea are the following:/ ) { $koll = 1; } } if ( $koll == 1 ) {print " ok! \n";}else{ print "  FAILED !\n";} }
		if ( $i == 8 ) { while(<DATA>){ if ( $_ =~ /RUNTIME=/ ) { $koll = 1; } } if ( $koll == 1 ) {print " ok! \n";}else{ print "  FAILED !\n";} }
		close(DATA);
	}
}

if (-e "MPI_VERSION/uquantchem.mpi" ) {
        print "\n";
	print "***********************************************\n";
	print "  Performing test on MPI version of UQUANTCHEM \n";
	print "***********************************************\n";
	print "\n";
 	$dir = `pwd`;
	chomp $dir;
	for ( $i =1; $i <= 9; $i++ ) {
		chdir("$dir/TESTS/RUN$i/");
		if ( $i == 1 ) { print "Testing the Hartree-Fock implementation"}
		if ( $i == 2 ) { print "Testing the DFT-PBE implementation"}
		if ( $i == 3 ) { print "Testing the DFT-B3LYP implementation"}
		if ( $i == 4 ) { print "Testing the MP2 implementation"}
		if ( $i == 5 ) { print "Testing the CISD implementation"}
		if ( $i == 6 ) { print "Testing the QDMC implementation"}
		if ( $i == 7 ) { print "Testing the relaxation implementation"}
		if ( $i == 8 ) { print "Testing the tddft implementation"}
		if ( $i == 9 ) { print "Testing the molecular dynamics implementation"}
     
		if ( -e "OUT.dat" ) {system("rm OUT.dat");}
		system("mpirun $dir/MPI_VERSION/uquantchem.mpi > OUT.dat");

	       	open(DATA,"OUT.dat");
		$koll = 0;
		if ( $i == 1 ) { while(<DATA>){ if ( $_ =~ /Hartree-Fock energy:/ ) { $koll = 1; } } if ( $koll == 1 ) {print " ok! \n";}else{ print "  FAILED !\n";} }
		if ( $i == 2 ) { while(<DATA>){ if ( $_ =~ /PBE energy:/ ) { $koll = 1; } } if ( $koll == 1 ) {print " ok! \n";}else{ print "  FAILED !\n";} }
		if ( $i == 3 ) { while(<DATA>){ if ( $_ =~ /CALCULATED INTERATOMC FORCES:/ ) { $koll = 1; } } if ( $koll == 1 ) {print " ok! \n";}else{ print "  FAILED !\n";} }
		if ( $i == 4 ) { while(<DATA>){ if ( $_ =~ /Results from the MP2 calculation:/ ) { $koll = 1; } } if ( $koll == 1 ) {print " ok! \n";}else{ print "  FAILED !\n";} }
		if ( $i == 5 ) { while(<DATA>){ if ( $_ =~ /RUNTIME=/ ) { $koll = 1; } } if ( $koll == 1 ) {print " ok! \n";}else{ print "  FAILED !\n";} }
		if ( $i == 6 ) { while(<DATA>){ if ( $_ =~ /RUNTIME=/ ) { $koll = 1; } } if ( $koll == 1 ) {print " ok! \n";}else{ print "  FAILED !\n";} }
		if ( $i == 7 ) { while(<DATA>){ if ( $_ =~ /The relaxed positions of the nuclea are the following:/ ) { $koll = 1; } } if ( $koll == 1 ) {print " ok! \n";}else{ print "  FAILED !\n";} }
		if ( $i == 8 ) { while(<DATA>){ if ( $_ =~ /RUNTIME=/ ) { $koll = 1; } } if ( $koll == 1 ) {print " ok! \n";}else{ print "  FAILED !\n";} }
		if ( $i == 9 ) { while(<DATA>){ if ( $_ =~ /RUNTIME=/ ) { $koll = 1; } } if ( $koll == 1 ) {print " ok! \n";}else{ print "  FAILED !\n";} }
		close(DATA);
	}
	chdir("$dir");
}

