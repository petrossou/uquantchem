SUBROUTINE getK(P,NB,Ints,Kout)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NB
      DOUBLE PRECISION, INTENT(IN) :: P(NB,NB),Ints(NB,NB,NB,NB)
      DOUBLE PRECISION, INTENT(OUT) :: Kout(NB,NB)
      DOUBLE PRECISION :: MAT(NB,NB)
      INTEGER :: I,J,N,M

      Kout(:,:) = 0.0d0

      DO I=1,NB
        DO J=I,NB
          MAT = Ints(I,:,:,J)
          Kout(I,J) = Kout(I,J) + SUM(MAT*TRANSPOSE(P)) !TRACE(NB,MATMUL(MAT,TRANSPOSE(P)))
          IF ( I .NE. J ) Kout(J,I) = Kout(I,J)
        ENDDO
      ENDDO
END SUBROUTINE getK
           


