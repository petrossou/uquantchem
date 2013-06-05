SUBROUTINE getJ(P,NB,Ints,Jout)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NB
      DOUBLE PRECISION, INTENT(IN) :: P(NB,NB),Ints(NB,NB,NB,NB)
      DOUBLE PRECISION, INTENT(OUT) :: Jout(NB,NB)
      DOUBLE PRECISION :: MAT(NB,NB)
      INTEGER :: I,J,N,M

      Jout(:,:) = 0.0d0

      DO I=1,NB
        DO J=I,NB
          MAT = Ints(I,J,:,:)
          Jout(I,J) = Jout(I,J) + SUM(MAT*P) !TRACE(NB,MATMUL(MAT,P))
          IF ( I .NE. J ) Jout(J,I) = Jout(I,J)
        ENDDO
      ENDDO
END SUBROUTINE getJ
           


