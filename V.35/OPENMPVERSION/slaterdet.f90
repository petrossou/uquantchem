FUNCTION slaterdet(N,N3,NS,BAS,C,r,UP)
     ! This fuction calculates the slater determiant from:
     ! (a) N = Total Number of electrons, (b) BAS = Hartree Fock basis
     ! (c) C = Fockian eigenvectors, (d) r = positions the electrons 3N-long array 
     ! (e) UP = .TRUE. Spin up slater determinant is calculated, otherwise
     ! spin-down slater determinant is to be calculated.
     ! (f) NS = slater determinant dimension
     USE datatypemodule
     IMPLICIT NONE
     DOUBLE PRECISION :: slaterdet
     INTEGER :: N,N3,NS
     TYPE(BASIS) :: BAS
     DOUBLE PRECISION :: C(BAS%NBAS,BAS%NBAS)
     DOUBLE PRECISION :: r(N3)
     LOGICAL :: UP
     DOUBLE PRECISION, EXTERNAL :: hforbitalval, det
     INTEGER :: I,J,OFFS 
     DOUBLE PRECISION :: MATRIX(NS,NS),VECT(3)

     MATRIX(:,:) = 0.0d0

     IF ( UP ) OFFS = 0
     IF ( .not. UP ) OFFS = ( N - MOD(N,2) )/2
     
     DO I=1,NS
        DO J=1,NS
            VECT(1) = r(3*(J+OFFS-1) + 1)
            VECT(2) = r(3*(J+OFFS-1) + 2)
            VECT(3) = r(3*(J+OFFS-1) + 3)
            MATRIX(I,J) = hforbitalval(BAS,C,I,VECT)
        ENDDO
     ENDDO

     slaterdet = det(MATRIX,NS)

END FUNCTION slaterdet
