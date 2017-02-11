FUNCTION laplaceslaterdet(I,N,N3,NS,BAS,C,r,UP)
     ! This fuction calculates the laplacian of a slater determiant, with 
     ! respect to the electron coordinate I, from
     ! (a) N = Total Number of electrons, (b) BAS = Hartree Fock basis
     ! (c) C = Fockian eigenvectors, (d) r = positions the electrons 3N-long array 
     ! (e) NS = Slaterdeterminant dimension.
     ! (f) UP = .TRUE. if I = spin upp orbital, else UP = .FALSE.
     USE datatypemodule
     IMPLICIT NONE
     DOUBLE PRECISION :: laplaceslaterdet
     INTEGER :: I,N,N3,NS
     TYPE(BASIS) :: BAS
     DOUBLE PRECISION :: C(BAS%NBAS,BAS%NBAS)
     DOUBLE PRECISION :: r(N3)
     LOGICAL :: UP
     DOUBLE PRECISION, EXTERNAL :: hforbitalval, laplacehforbitalval, det
     INTEGER :: K,J,OFFS,II,NONZERO
     DOUBLE PRECISION :: MATRIX(NS,NS),VECT(3)

     MATRIX(:,:) = 0.0d0
     
     NONZERO = 0

     IF ( UP ) OFFS = 0

     IF ( .not. UP ) OFFS = ( N - MOD(N,2) )/2
    
     II = I - OFFS
     
     DO K=1,NS
        DO J=1,NS
            VECT(1) = r(3*(J+OFFS-1) + 1)
            VECT(2) = r(3*(J+OFFS-1) + 2)
            VECT(3) = r(3*(J+OFFS-1) + 3)
            IF ( J .NE. II ) THEN
                MATRIX(K,J) = hforbitalval(BAS,C,K,VECT)
            ELSE
                MATRIX(K,J) = laplacehforbitalval(BAS,C,K,VECT)
                ! To make the Function return a gradient 
                ! only if the gradient coordinate equals 
                ! some of the electron coordinates. This 
                ! makes the code more robust.
                NONZERO = 1
            ENDIF
        ENDDO
     ENDDO

     laplaceslaterdet = NONZERO*det(MATRIX,NS)

END FUNCTION laplaceslaterdet
