SUBROUTINE Fnurec(NUMAX,FNUVEC,X)
      ! This subroutine calculates the values F_nu(X) (see Cook Book
      ! p.280 ) for nu=0,1,2,...,NUMAX, Through the recursion formula
      ! given on the bottom of p. 280 of the Cook Book
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NUMAX
      DOUBLE PRECISION, INTENT(IN) :: X
      DOUBLE PRECISION, INTENT(OUT) :: FNUVEC(NUMAX+1)
      INTEGER :: I
      DOUBLE PRECISION, EXTERNAL :: Fnu

      FNUVEC(NUMAX+1) = Fnu(NUMAX,X)
      DO I=1,NUMAX
        FNUVEC(NUMAX+1-I) = (EXP(-X) + 2.0d0*X*FNUVEC(NUMAX+2-I))/(2.0d0*(NUMAX+1-I)- 1.0d0)
      ENDDO
END SUBROUTINE Fnurec
       
