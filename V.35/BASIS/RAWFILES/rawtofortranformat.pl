#!/usr/bin/perl -w
use POSIX;
$header = "DGA1-DFT-Exchange-Fitting";

$FILENAMEF = $header.".dat";
$FILENAMERAW = $header.".raw";

open(FORTRANFILE,">$FILENAMEF");
open(RAWFILE,"$FILENAMERAW");

print FORTRANFILE "$header\n";
print FORTRANFILE "BASIS:\n";

$STARTREAD = 0;
$nstar = 0;
$Z = 0;
$ZZ[0] = 0;
while(<RAWFILE>){
        if ( $_ =~ /\$basis/ ) { $STARTREAD = 1;}
        if ( $STARTREAD  == 1 ) {
                @read = split(" ",$_);
                $size = @read;
                if ( $_ =~ /\*/ && $size == 1 ) {$nstar++;}
                $_=~ s/\n//g;
                if ( $size == 2 ) {
                        if ( $read[0] =~ /^\D/ && $read[0] !~ /\*/) { 
				$Z++; $NAME[$Z] = ucfirst($read[0]); $nbas = 0;
                               	if ( $NAME[$Z] =~ /H/ && length($NAME[$Z]) == 1 ) { $ZZ[$Z] = 1; }
                                if ( $NAME[$Z] =~ /He/ ) { $ZZ[$Z] = 2; }
                                if ( $NAME[$Z] =~ /Li/ ) { $ZZ[$Z] = 3; }
                                if ( $NAME[$Z] =~ /Be/ ) { $ZZ[$Z] = 4; }
                                if ( $NAME[$Z] =~ /B/ && length($NAME[$Z]) == 1) { $ZZ[$Z] = 5; }
                                if ( $NAME[$Z] =~ /C/ && length($NAME[$Z]) == 1) { $ZZ[$Z] = 6; }
                                if ( $NAME[$Z] =~ /N/ && length($NAME[$Z]) == 1) { $ZZ[$Z] = 7; }
                                if ( $NAME[$Z] =~ /O/ && length($NAME[$Z]) == 1) { $ZZ[$Z] = 8; }
                                if ( $NAME[$Z] =~ /F/ && length($NAME[$Z]) == 1) { $ZZ[$Z] = 9; }
                                if ( $NAME[$Z] =~ /Ne/ ) { $ZZ[$Z] = 10; }
                                if ( $NAME[$Z] =~ /Na/ ) { $ZZ[$Z] = 11; }
                                if ( $NAME[$Z] =~ /Mg/ ) { $ZZ[$Z] = 12; }
                                if ( $NAME[$Z] =~ /Al/ ) { $ZZ[$Z] = 13; }
                                if ( $NAME[$Z] =~ /Si/ ) { $ZZ[$Z] = 14; }
                                if ( $NAME[$Z] =~ /P/ && length($NAME[$Z]) == 1 ) { $ZZ[$Z] = 15; }
                                if ( $NAME[$Z] =~ /S/ && length($NAME[$Z]) == 1 ) { $ZZ[$Z] = 16; }
                                if ( $NAME[$Z] =~ /Cl/ ) { $ZZ[$Z] = 17; }
                                if ( $NAME[$Z] =~ /Ar/ ) { $ZZ[$Z] = 18; }
                                if ( $NAME[$Z] =~ /K/ && length($NAME[$Z]) == 1 ) { $ZZ[$Z] = 19; }
                                if ( $NAME[$Z] =~ /Ca/ ) { $ZZ[$Z] = 20; }
                                if ( $NAME[$Z] =~ /Sc/ ) { $ZZ[$Z] = 21; }
                                if ( $NAME[$Z] =~ /Ti/ ) { $ZZ[$Z] = 22; }
                                if ( $NAME[$Z] =~ /V/ && length($NAME[$Z]) == 1 ) { $ZZ[$Z] = 23; }
                                if ( $NAME[$Z] =~ /Cr/ ) { $ZZ[$Z] = 24; }
                                if ( $NAME[$Z] =~ /Mn/ ) { $ZZ[$Z] = 25; }
                                if ( $NAME[$Z] =~ /Fe/ ) { $ZZ[$Z] = 26; }
                                if ( $NAME[$Z] =~ /Co/ ) { $ZZ[$Z] = 27; }
                                if ( $NAME[$Z] =~ /Ni/ ) { $ZZ[$Z] = 28; }
                                if ( $NAME[$Z] =~ /Cu/ ) { $ZZ[$Z] = 29; }
                                if ( $NAME[$Z] =~ /Zn/ ) { $ZZ[$Z] = 30; }
                                if ( $NAME[$Z] =~ /Ga/ ) { $ZZ[$Z] = 31; }
                                if ( $NAME[$Z] =~ /Ge/ ) { $ZZ[$Z] = 32; }
                                if ( $NAME[$Z] =~ /As/ ) { $ZZ[$Z] = 33; }
                                if ( $NAME[$Z] =~ /Se/ ) { $ZZ[$Z] = 34; }
                                if ( $NAME[$Z] =~ /Br/ ) { $ZZ[$Z] = 35; }
                                if ( $NAME[$Z] =~ /Kr/ ) { $ZZ[$Z] = 36; }
                                if ( $NAME[$Z] =~ /Rb/ ) { $ZZ[$Z] = 37; }
                                if ( $NAME[$Z] =~ /Sr/ ) { $ZZ[$Z] = 38; }
                                if ( $NAME[$Z] =~ /Y/ && length($NAME[$Z]) == 1 ) { $ZZ[$Z] = 39; }
                                if ( $NAME[$Z] =~ /Zr/ ) { $ZZ[$Z] = 40; }
                                if ( $NAME[$Z] =~ /Nb/ ) { $ZZ[$Z] = 41; }
                                if ( $NAME[$Z] =~ /Mo/ ) { $ZZ[$Z] = 42; }
                                if ( $NAME[$Z] =~ /Tc/ ) { $ZZ[$Z] = 43; }
                                if ( $NAME[$Z] =~ /Ru/ ) { $ZZ[$Z] = 44; }
                                if ( $NAME[$Z] =~ /Rh/ ) { $ZZ[$Z] = 45; }
                                if ( $NAME[$Z] =~ /Pd/ ) { $ZZ[$Z] = 46; }
                                if ( $NAME[$Z] =~ /Ag/ ) { $ZZ[$Z] = 47; }
                                if ( $NAME[$Z] =~ /Cd/ ) { $ZZ[$Z] = 48; }
                                if ( $NAME[$Z] =~ /In/ ) { $ZZ[$Z] = 49; }
                                if ( $NAME[$Z] =~ /Sn/ ) { $ZZ[$Z] = 50; }
                                if ( $NAME[$Z] =~ /Sb/ ) { $ZZ[$Z] = 51; }
                                if ( $NAME[$Z] =~ /Te/ ) { $ZZ[$Z] = 52; }
                                if ( $NAME[$Z] =~ /I/ && length($NAME[$Z]) == 1 ) { $ZZ[$Z] = 53; }
                                if ( $NAME[$Z] =~ /Xe/ ) { $ZZ[$Z] = 54; }
                                if ( $NAME[$Z] =~ /Cs/ ) { $ZZ[$Z] = 55; }
                                if ( $NAME[$Z] =~ /Ba/ ) { $ZZ[$Z] = 56; }
                                if ( $NAME[$Z] =~ /La/ ) { $ZZ[$Z] = 57; }
                                if ( $NAME[$Z] =~ /Ce/ ) { $ZZ[$Z] = 58; }
                                if ( $NAME[$Z] =~ /Pr/ ) { $ZZ[$Z] = 59; }
                                if ( $NAME[$Z] =~ /Nd/ ) { $ZZ[$Z] = 60; }
                                if ( $NAME[$Z] =~ /Pm/ ) { $ZZ[$Z] = 61; }
                                if ( $NAME[$Z] =~ /Sm/ ) { $ZZ[$Z] = 62; }
                                if ( $NAME[$Z] =~ /Eu/ ) { $ZZ[$Z] = 63; }
                                if ( $NAME[$Z] =~ /Gd/ ) { $ZZ[$Z] = 64; }
                                if ( $NAME[$Z] =~ /Tb/ ) { $ZZ[$Z] = 65; }
                                if ( $NAME[$Z] =~ /Dy/ ) { $ZZ[$Z] = 66; }
                                if ( $NAME[$Z] =~ /Ho/ ) { $ZZ[$Z] = 67; }
                                if ( $NAME[$Z] =~ /Er/ ) { $ZZ[$Z] = 68; }
                                if ( $NAME[$Z] =~ /Tm/ ) { $ZZ[$Z] = 69; }
                                if ( $NAME[$Z] =~ /Yb/ ) { $ZZ[$Z] = 70; }
                                if ( $NAME[$Z] =~ /Lu/ ) { $ZZ[$Z] = 71; }
                                if ( $NAME[$Z] =~ /Hf/ ) { $ZZ[$Z] = 72; }
                                if ( $NAME[$Z] =~ /Ta/ ) { $ZZ[$Z] = 73; }
                                if ( $NAME[$Z] =~ /W/ && length($NAME[$Z]) == 1 ) { $ZZ[$Z] = 74; }
                                if ( $NAME[$Z] =~ /Re/ ) { $ZZ[$Z] = 75; }
                                if ( $NAME[$Z] =~ /Os/ ) { $ZZ[$Z] = 76; }
                                if ( $NAME[$Z] =~ /Ir/ ) { $ZZ[$Z] = 77; }
                                if ( $NAME[$Z] =~ /Pt/ ) { $ZZ[$Z] = 78; }
                                if ( $NAME[$Z] =~ /Au/ ) { $ZZ[$Z] = 79; }
                                if ( $NAME[$Z] =~ /Hg/ ) { $ZZ[$Z] = 80; }
                                if ( $NAME[$Z] =~ /Tl/ ) { $ZZ[$Z] = 81; }
                                if ( $NAME[$Z] =~ /Pb/ ) { $ZZ[$Z] = 82; }
                                if ( $NAME[$Z] =~ /Bi/ ) { $ZZ[$Z] = 83; }
                                if ( $NAME[$Z] =~ /Po/ ) { $ZZ[$Z] = 84; }
                                if ( $NAME[$Z] =~ /At/ ) { $ZZ[$Z] = 85; }
                                if ( $NAME[$Z] =~ /Rn/ ) { $ZZ[$Z] = 86; }
                                if ( $NAME[$Z] =~ /Fr/ ) { $ZZ[$Z] = 87; }
                                if ( $NAME[$Z] =~ /Ra/ ) { $ZZ[$Z] = 88; }
                                if ( $NAME[$Z] =~ /Ac/ ) { $ZZ[$Z] = 89; }
                                if ( $NAME[$Z] =~ /Th/ ) { $ZZ[$Z] = 90; }
                                if ( $NAME[$Z] =~ /Pa/ ) { $ZZ[$Z] = 91; }
                                if ( $NAME[$Z] =~ /U/ && length($NAME[$Z]) == 1 ) { $ZZ[$Z] = 92; }
                                if ( $NAME[$Z] =~ /Np/ ) { $ZZ[$Z] = 93; }
                                if ( $NAME[$Z] =~ /Pu/ ) { $ZZ[$Z] = 94; }
                                if ( $NAME[$Z] =~ /Am/ ) { $ZZ[$Z] = 95; }
                                if ( $NAME[$Z] =~ /Cm/ ) { $ZZ[$Z] = 96; }
                                if ( $NAME[$Z] =~ /Bk/ ) { $ZZ[$Z] = 97; }
                                if ( $NAME[$Z] =~ /Cf/ ) { $ZZ[$Z] = 98; }
                                if ( $NAME[$Z] =~ /Es/ ) { $ZZ[$Z] = 99; }
                                if ( $NAME[$Z] =~ /Fm/ ) { $ZZ[$Z] = 100; }
                                if ( $NAME[$Z] =~ /Md/ ) { $ZZ[$Z] = 101; }
                                if ( $NAME[$Z] =~ /No/ ) { $ZZ[$Z] = 102; }
                                if ( $NAME[$Z] =~ /Lw/ ) { $ZZ[$Z] = 103; }
                                if ( $NAME[$Z] =~ /Rf/ ) { $ZZ[$Z] = 104; }
                                if ( $NAME[$Z] =~ /Db/ ) { $ZZ[$Z] = 105; }
                                if ( $NAME[$Z] =~ /Sg/ ) { $ZZ[$Z] = 106; }
                                if ( $NAME[$Z] =~ /Bh/ ) { $ZZ[$Z] = 107; }
                                if ( $NAME[$Z] =~ /Hs/ ) { $ZZ[$Z] = 108; }
                                if ( $NAME[$Z] =~ /Mt/ ) { $ZZ[$Z] = 109; }
                                if ( $NAME[$Z] =~ /Ds/ ) { $ZZ[$Z] = 110; }
                                if ( $NAME[$Z] =~ /Rg/ ) { $ZZ[$Z] = 111; }
                                if ( $NAME[$Z] =~ /Uub/ ) { $ZZ[$Z] = 112; }
                                if ( $NAME[$Z] =~ /Uut/ ) { $ZZ[$Z] = 113; }
                                if ( $NAME[$Z] =~ /Uuq/ ) { $ZZ[$Z] = 114; }
                                if ( $NAME[$Z] =~ /Uup/ ) { $ZZ[$Z] = 115; }
                                if ( $NAME[$Z] =~ /Uuh/ ) { $ZZ[$Z] = 116; }
                                if ( $NAME[$Z] =~ /Uus/ ) { $ZZ[$Z] = 117; }
                                if ( $NAME[$Z] =~ /Uuo/ ) { $ZZ[$Z] = 118; }
				if ( $Z == 1 ){print "The following atoms are included in the $header basis:\n";} 
				print "$Z $ZZ[$Z] $NAME[$Z]\n";
			}
                        if ( $read[0] =~/^\d/ && $read[1] =~ /\D/ && $read[1] !~ /\./) { $nbas++; $n = 0; $NPRIM[$nbas] = $read[0]; $L = $read[1]}
                        if ( $read[0] =~/\./ && $read[1] =~ /\./ ){ $n++; $basfnk[$nbas][$n] = $_; $ANGULARMOM[$nbas][$n] = $L;};
                }
                $MOD = int($nstar-3*floor($nstar/3));
                if ( $nstar > 2 && $MOD == 0 ) {
                        $nstar = 1;
                        if ( $ZZ[$Z]-$ZZ[$Z-1] > 1 ){ for ( $ii=$ZZ[$Z-1]+1; $ii<= $ZZ[$Z]-1;$ii++){ print FORTRANFILE "XX $ii 0\n";print FORTRANFILE "\n"; }}
                        print FORTRANFILE "$NAME[$Z] $ZZ[$Z] $nbas";
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
