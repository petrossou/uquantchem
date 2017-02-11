RECURSIVE FUNCTION atomicradii(Z)
      ! This function takes the atomic number as input 
      ! and returns the empirical radii of Bragg and Slater,
      ! J. C. Slater, J. Chem. Phys. 41, 3199 (1964)
      ! The numbers in parathesis are the corresponding 
      ! theoretical estimates of the atomic radii.
      ! If there is no number in paranthesis then 
      ! the theoretical estimate has been used.
      IMPLICIT NONE
      DOUBLE PRECISION :: atomicradii
      INTEGER :: Z
      DOUBLE PRECISION :: r(100)

      ! Radii are here given in Angstrom.
      r(1) = 0.350d0 ! Here we use the data of Becke J. Chem. Phys. 88, 2547 (1988)
                     ! Instead of Slaters r = 0.25 (53)
      r(2) = 1.310d0
      r(3) = 1.450d0 ! (1.67)
      r(4) = 1.050d0 ! (1.12)
      r(5) = 0.850d0 ! (0.87)
      r(6) = 0.700d0 ! (0.67)
      r(7) = 0.650d0 ! (0.56)
      r(8) = 0.600d0 ! (0.48)
      r(9) = 0.500d0 ! (0.42)
      r(10) = 0.380d0
      r(11) = 1.800d0 ! (1.90)
      r(12) = 1.500d0 ! (1.45)
      r(13) = 1.250d0 ! (1.18)
      r(14) = 1.100d0 ! (1.11)
      r(15) = 1.000d0 ! (0.98)
      r(16) = 1.000d0 ! (0.88)
      r(17) = 1.000d0 ! (0.79)
      r(18) = 0.710d0 ! (0.71)
      r(19) = 2.200d0 ! (2.43)
      r(20) = 1.800d0 ! (1.94)
      r(21) = 1.600d0 ! (1.84)
      r(22) = 1.400d0 ! (1.76)
      r(23) = 1.350d0 ! (1.71)
      r(24) = 1.400d0 ! (1.66)
      r(25) = 1.400d0 ! (1.61)
      r(26) = 1.400d0 ! (1.56)
      r(27) = 1.350d0 ! (1.52)
      r(28) = 1.350d0 ! (1.49)
      r(29) = 1.350d0 ! (1.45)
      r(30) = 1.350d0 ! (1.42)
      r(31) = 1.300d0 ! (1.36)
      r(32) = 1.250d0 ! (1.25)
      r(33) = 1.150d0 ! (1.14)
      r(34) = 1.150d0 ! (1.03)
      r(35) = 1.150d0 ! (0.94)
      r(36) = 0.880d0 
      r(37) = 2.350d0 ! (2.65)
      r(38) = 2.000d0 ! (2.19)
      r(39) = 1.800d0 ! (2.12)
      r(40) = 1.550d0 ! (2.06)
      r(41) = 1.540d0 ! (1.98)
      r(42) = 1.450d0 ! (1.90)
      r(43) = 1.350d0 ! (1.83)
      r(44) = 1.300d0 ! (1.78)
      r(45) = 1.350d0 ! (1.73)
      r(46) = 1.400d0 ! (1.69)
      r(47) = 1.600d0 ! (1.65)
      r(48) = 1.550d0 ! (1.61)
      r(49) = 1.550d0 ! (1.56)
      r(50) = 1.450d0 ! (1.45)
      r(51) = 1.450d0 ! (1.33)
      r(52) = 1.400d0 ! (1.23)
      r(53) = 1.400d0 ! (1.15)
      r(54) = 1.080d0
      r(55) = 2.600d0 ! (2.98)
      r(56) = 2.150d0 ! (2.53)
      r(57) = 1.950d0 ! (No data)
      r(58) = 1.850d0 ! (No data)
      r(59) = 1.850d0 ! (2.47)
      r(60) = 1.850d0 ! (2.06)
      r(61) = 1.850d0 ! (2.05)
      r(62) = 1.850d0 ! (2.38)
      r(63) = 1.850d0 ! (2.31)
      r(64) = 1.800d0 ! (2.33)
      r(65) = 1.750d0 ! (2.25)
      r(66) = 1.750d0 ! (2.28)
      r(67) = 1.750d0 ! (No data)
      r(68) = 1.750d0 ! (2.26)
      r(69) = 1.750d0 ! (2.22)
      r(70) = 1.750d0 ! (2.22)
      r(71) = 1.750d0 ! (2.17)
      r(72) = 1.550d0 ! (2.08)
      r(73) = 1.450d0 ! (2.00)
      r(74) = 1.350d0 ! (1.93)
      r(75) = 1.350d0 ! (1.88)
      r(76) = 1.300d0 ! (1.85)
      r(77) = 1.350d0 ! (1.80)
      r(78) = 1.350d0 ! (1.77)
      r(79) = 1.350d0 ! (1.74)
      r(80) = 1.500d0 ! (1.71)
      r(81) = 1.900d0 ! (1.56)
      r(82) = 1.800d0 ! (1.54)
      r(83) = 1.600d0 ! (1.43)
      r(84) = 1.900d0 ! (1.35)
      r(85) = 1.500d0 ! Petros data
      r(86) = 1.200d0 
      r(87) = 1.500d0 ! Petros data
      r(88) = 2.150d0 ! (No data)
      r(89) = 1.950d0 ! (No data)
      r(90) = 1.800d0 ! (No data)
      r(91) = 1.800d0 ! (No data)
      r(92) = 1.750d0 ! (No data)
      r(93) = 1.750d0 ! (No data)
      r(94) = 1.750d0 ! (No data)
      r(95) = 1.750d0 ! (No data)
      r(96) = 1.500d0 ! Petros data
      r(97) = 1.500d0 ! Petros data
      r(98) = 1.500d0 ! Petros data
      r(99) = 1.500d0 ! Petros data
      r(100) = 1.500d0 ! Petros data

      IF ( Z .LE. 100 ) THEN
             atomicradii = r(Z)/0.529177208590d0
      ELSE
             atomicradii = 1.50d0/0.529177208590d0
      ENDIF
END FUNCTION atomicradii
