FUNCTION primoverlap(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2)
        ! This function calculates the overlap integral
        ! S = S_x*S_y*S_z as given by the expressions 
        ! on the pages 234-238 in the Cook book
        IMPLICIT NONE
        INTEGER :: L1,M1,N1,L2,M2,N2
        DOUBLE PRECISION :: primoverlap
        DOUBLE PRECISION :: alpha1,alpha2,A(3),B(3),PA(3),PB(3),P(3),C
        DOUBLE PRECISION :: Sx,Sy,Sz,gama
        INTEGER :: I,J,K
        DOUBLE PRECISION, EXTERNAL :: fj,dfac
        DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932384626433832795028841970d0

        Sx = 0.0d0
        Sy = 0.0d0
        Sz = 0.0d0
        
        gama = alpha1 + alpha2

        P = (A*alpha1 + B*alpha2)/gama

        C = (alpha1*alpha2/gama)*DOT_PRODUCT(A-B,A-B)

        PA = P-A
        PB = P-B

        DO J=0,(L1+L2-MOD(L1+L2,2))/2
               Sx = Sx+fj(2*J,L1,L2,PA(1),PB(1))*dfac(2*J-1)/((2.0d0*gama)**J)
        ENDDO
        Sx = sqrt(pi/gama)*Sx

        DO J=0,(M1+M2-MOD(M1+M2,2))/2
               Sy = Sy+fj(2*J,M1,M2,PA(2),PB(2))*dfac(2*J-1)/((2.0d0*gama)**J)
        ENDDO
        Sy = sqrt(pi/gama)*Sy
        
        DO J=0,(N1+N2-MOD(N1+N2,2))/2
               Sz = Sz+fj(2*J,N1,N2,PA(3),PB(3))*dfac(2*J-1)/((2.0d0*gama)**J)
        ENDDO
        Sz = sqrt(pi/gama)*Sz
        
        ! Temporary solution. Apparently in pyquante the primitive gaussian
        ! orbitals of corresponding to the one and the same basis function are
        ! normalized. This is a quick fix in order to check that the overlap and
        ! kinetic energy matrices are the same as in Pyquante. This only works
        ! for S-orbitals for L .NE. 0 see row 117 in PGBF.py. To do this
        ! properly the contraction coefficients should instead be redefined as
        ! Contr --> primnorm*Contr

        primoverlap = exp(-C)*Sx*Sy*Sz

 END FUNCTION primoverlap
