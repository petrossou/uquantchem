SUBROUTINE gradslaterdet(I,N,N3,NS,BAS,C,r,UP,grad)
     ! This subroutine calculates the gradient of a slater determiant, with 
     ! respect to the electron coordinate I, from
     ! (a) N = Total Number of electrons, (b) BAS = Hartree Fock basis
     ! (c) C = Fockian eigenvectors, (d) r = positions the electrons 3N-long array 
     ! (e) UP = .TRUE. if I = spin upp orbital, else UP = .FALSE.
     ! (f) NS = Slaterdeterminant dimension
     USE datatypemodule
     IMPLICIT NONE
     DOUBLE PRECISION, INTENT(OUT) :: grad(3)
     INTEGER, INTENT(IN) :: I,N,N3,NS
     TYPE(BASIS), INTENT(IN) :: BAS
     DOUBLE PRECISION, INTENT(IN) :: C(BAS%NBAS,BAS%NBAS)
     DOUBLE PRECISION, INTENT(IN) :: r(N3)
     LOGICAL, INTENT(IN) :: UP
     DOUBLE PRECISION, EXTERNAL :: hforbitalval, det
     EXTERNAL :: gradhforbitalval
     DOUBLE PRECISION :: gradhf(3)
     INTEGER :: K,J,OFFS,II,NONZERO
     DOUBLE PRECISION :: MATRIX1(NS,NS),MATRIX2(NS,NS),MATRIX3(NS,NS),VECT(3)

     MATRIX1(:,:) = 0.0d0
     MATRIX2(:,:) = 0.0d0
     MATRIX3(:,:) = 0.0d0
      
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
                MATRIX1(K,J) = hforbitalval(BAS,C,K,VECT)
                MATRIX2(K,J) = hforbitalval(BAS,C,K,VECT)
                MATRIX3(K,J) = hforbitalval(BAS,C,K,VECT)
            ELSE
                CALL gradhforbitalval(BAS,C,K,VECT,gradhf)
                MATRIX1(K,J) = gradhf(1)
                MATRIX2(K,J) = gradhf(2)
                MATRIX3(K,J) = gradhf(3)
                ! To make the Function return a gradient 
                ! only if the gradient coordinate equals 
                ! some of the electron coordinates. This 
                ! makes the code more robust.
                NONZERO = 1
            ENDIF
        ENDDO
     ENDDO
     
     grad(1) = NONZERO*det(MATRIX1,NS)
     grad(2) = NONZERO*det(MATRIX2,NS)
     grad(3) = NONZERO*det(MATRIX3,NS)

END SUBROUTINE gradslaterdet
